from analysis.pipeline import PipelineComponent
from analysis.models import Analysis
from tempfile import mkdtemp, NamedTemporaryFile
from shutil import rmtree
from Bio import SeqIO, Phylo
import os
from django.conf import settings
import subprocess
from io import StringIO
import csv
from datetime import datetime
import copy


class StartPipelineComponent(PipelineComponent):
    name = "start_pipeline"

    def analysis(self, report):
        analysis_entry = Analysis.objects.get(id=report['analysis'])
        analysis_entry.start_time = datetime.now()
        analysis_entry.save()


class EndPipelineComponent(PipelineComponent):
    name = "end_pipeline"

    def analysis(self, report):
        analysis_entry = Analysis.objects.get(id=report['analysis'])
        analysis_entry.complete_time = datetime.now()
        analysis_entry.save()


class SetupGbkPipelineComponent(PipelineComponent):
    name = "setup_gbk"
    result_types = ["gbk_paths"]

    def analysis(self, report):
        analysis_entry = Analysis.objects.get(id=report['analysis'])
        genomes = analysis_entry.genomes

        report['gbk_paths'] = dict()
        for genome in genomes.all():
            report['gbk_paths'][str(genome.id)] = genome.gbk.path


class GbkMetadataComponent(PipelineComponent):
    name = "gbk_metadata"
    dependencies = ["gbk_paths"]
    result_types = ["gbk_metadata"]

    @staticmethod
    def get_genome_size(input_path):
        for record in SeqIO.parse(open(input_path),"genbank"):
            return len(record.seq)

    def analysis(self, report):
        output = dict()
        for genome_id in report["gbk_paths"].keys():
            output[str(genome_id)] = dict()
            output[str(genome_id)]["size"] = self.get_genome_size(report["gbk_paths"][genome_id])
        report["gbk_metadata"] = output


class ParsnpPipelineComponent(PipelineComponent):
    name = "parsnp"
    dependencies = ["gbk_paths"]
    result_types = ["newick"]
    temp_dir_path = None
    output_dir = settings.BIO_APP_TEMP_DIR + "parsnp/"
    temp_results_dir = None
    PARSNP_PATH = settings.PARSNP_PATH

    @staticmethod
    def convert_gbk_to_fna(input_path, output_path):
        input_handle = open(input_path, "r")
        output_handle = open(output_path, "w")

        for seq_record in SeqIO.parse(input_handle, "genbank"):
            output_handle.write(">%s %s\n%s\n" % (
                seq_record.id,
                seq_record.description,
                seq_record.seq))

        output_handle.close()
        input_handle.close()

    @staticmethod
    def parse_newick(input_path):
        tree = Phylo.read(input_path, 'newick')
        terminals = tree.get_terminals()

        for leaf in terminals:
            leaf.name = leaf.name.split('.')[0]

        processed_tree = StringIO()
        Phylo.write(tree, processed_tree, "newick")
        output = processed_tree.getvalue()
        processed_tree.close()

        return output

    def setup(self, report):
        self.temp_dir_path = mkdtemp()
        for gbk_id in report["gbk_paths"]:
            gbk_path = report["gbk_paths"][gbk_id]
            self.convert_gbk_to_fna(gbk_path,
                                    self.temp_dir_path+"/"+str(gbk_id)+".fna")

    def analysis(self, report):
        self.temp_results_dir = self.output_dir + str(report["analysis"])
        script_file = NamedTemporaryFile(delete=True)
        with open(script_file.name, 'w') as script:
            script.write("#!/bin/bash\n")
            script.write(self.PARSNP_PATH + " -r ! -d " + self.temp_dir_path +
                         " -o " + self.temp_results_dir + " -c YES")
            script.close()

        os.chmod(script_file.name, 0o0777)
        script_file.file.close()
        os.mkdir(self.temp_results_dir, 0o777)
        with open(self.temp_results_dir+"/logs", 'w') as logs:
            subprocess.check_call(script_file.name, stdout=logs)
        script_file.close()

        report["newick"] = self.parse_newick(self.temp_results_dir + "/parsnp.tree")

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)
        if self.temp_results_dir is not None and os.path.exists(self.temp_results_dir):
            rmtree(self.temp_results_dir)


class MauvePipelineComponent(PipelineComponent):
    name = "mauve"
    dependencies = ["newick", "gbk_paths"]
    result_types = ["alignment"]
    MAUVE_PATH = settings.MAUVE_PATH
    output_dir = settings.BIO_APP_TEMP_DIR + "mauve/"
    temp_dir_path = None
    backbone_file_name = "pairwise.backbone"
    minimum_homologous_region_size = 50

    def parse_mauve_pairwise_backbone(self, backbone_file):
        homologous_regions = []

        with open(backbone_file) as backbone:
            tsv_reader = csv.reader(backbone, delimiter='\t')
            next(tsv_reader)
            for row in tsv_reader:
                top_start = int(row[0])
                top_end = int(row[1])
                bottom_start = int(row[2])
                bottom_end = int(row[3])

                if abs(top_start - top_end) >= self.minimum_homologous_region_size \
                        and abs(bottom_start - bottom_end) >= self.minimum_homologous_region_size:
                    homologous_regions.append([top_start, top_end, bottom_start, bottom_end])

        return homologous_regions

    def mauve_subprocess(self, working_dir, gbk_list):
        absolute_working_dir = os.path.abspath(working_dir)
        scratch_path1 = absolute_working_dir + "/temp1"
        scratch_path2 = absolute_working_dir + "/temp2"
        os.mkdir(scratch_path1)
        os.mkdir(scratch_path2)

        script_file = NamedTemporaryFile(delete=True)
        mauve_absolute_path = os.path.abspath(self.MAUVE_PATH)

        with open(script_file.name, 'w') as script:
            script.write("#!/bin/bash\n")
            script.write(mauve_absolute_path + " --output=/dev/null  --scratch-path-1=" + scratch_path1
                         + " --scratch-path-2=" + scratch_path2 + " --backbone-output="
                         + absolute_working_dir + "/" + self.backbone_file_name + " ")

            for sequence in gbk_list:
                script.write(sequence+" ")
            script.close()

        os.chmod(script_file.name, 0o0755)
        script_file.file.close()

        with open(working_dir + "/logs", 'w') as logs:
            subprocess.check_call(script_file.name, cwd=working_dir, stdout=logs)
        script_file.close()

    def retrieve_mauve_results(self, backbone_list):
        merged_mauve_results = []

        for backbone in backbone_list:
            merged_mauve_results.append(self.parse_mauve_pairwise_backbone(backbone))

        return merged_mauve_results

    def setup(self, report):
        self.temp_dir_path = self.output_dir + str(report["analysis"])
        os.mkdir(self.temp_dir_path, 0o777)

    def analysis(self, report):
        tree = Phylo.read(StringIO(report["newick"]), 'newick')
        ordered_genome_ids = [int(clade.name) for clade in tree.get_terminals(order='preorder')]
        results_list = []

        for sequence_counter in range(len(ordered_genome_ids) - 1):
            first_genome_id = ordered_genome_ids[sequence_counter]
            second_genome_id = ordered_genome_ids[sequence_counter + 1]

            first_gbk_path = report["gbk_paths"][str(first_genome_id)]
            second_gbk_path = report["gbk_paths"][str(second_genome_id)]

            sub_results_dir = self.temp_dir_path + "/" + str(sequence_counter) + "/"
            os.mkdir(sub_results_dir, 0o777)
            results_list.append(sub_results_dir + self.backbone_file_name)

            self.mauve_subprocess(sub_results_dir,
                                  [first_gbk_path, second_gbk_path])

        report["alignment"] = self.retrieve_mauve_results(results_list)

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)


class SigiHMMPipelineComponent(PipelineComponent):
    name = "sigi"
    dependencies = ["gbk_paths"]
    result_types = ["sigi_gis"]
    output_dir = settings.BIO_APP_TEMP_DIR + "sigi/"
    temp_dir_path = None
    embl_files = {}
    SIGIHMM_PATH = settings.SIGIHMM_PATH
    SIGIHMM_EXE = settings.SIGIHMM_EXE

    @staticmethod
    def parse_sigi_gff(gff_file):
        gi_list = []

        with open(gff_file, 'r') as gff:
            start = 0
            end = 0
            island_flag = False

            for line in gff:
                # skip lines with parameters
                if line[0] == '#':
                    continue
                else:
                    cleaned_line = ' '.join(line.split())
                    gene_dict = cleaned_line.split(' ')
                    # at the start of a genomic island, set start and end of possible genomic island
                    if gene_dict[2] == 'PUTAL' and not island_flag:
                        start = gene_dict[3]
                        end = gene_dict[4]
                        island_flag = True
                    # continuation of current genomic island, change end and continue
                    elif gene_dict[2] == 'PUTAL' and island_flag:
                        end = gene_dict[4]
                    # end of genomic island, append current start and end to list
                    elif island_flag:
                        gi_list.append([int(start), int(end)])
                        island_flag = False
                    # not currently in a genomic island, continue to next line
                    elif not island_flag:
                        continue
                    # condition not included in above reached, throw an error
                    else:
                        raise Exception("Error occurred in sigi, unexpected condition reached")

        return gi_list

    def setup(self, report):
        self.temp_dir_path = self.output_dir + str(report["analysis"])
        os.mkdir(self.temp_dir_path, 0o777)

        for gbk_id in report["gbk_paths"]:
            basename = os.path.splitext(os.path.basename(report["gbk_paths"][gbk_id]))[0]
            with open(self.output_dir + str(report["analysis"]) + "/" + basename + ".embl", 'w') as embl_output:
                SeqIO.convert(report["gbk_paths"][gbk_id], "genbank", embl_output, "embl")
                self.embl_files[gbk_id] = os.path.abspath(embl_output.name)

    def analysis(self, report):
        script_file = NamedTemporaryFile(delete=True)
        report["sigi_gis"] = {}

        for embl_id in self.embl_files:
            sigi_gff_output = os.path.abspath(self.temp_dir_path + "/" + str(embl_id) + ".gff")
            sigi_output = os.path.abspath(self.temp_dir_path +"/" + str(embl_id) + ".embl")
            with open(script_file.name, 'w') as script:
                script.write("#!/bin/bash\n")
                script.write("/usr/bin/java " + self.SIGIHMM_EXE + " input=" +
                             self.embl_files[embl_id] + " output=" + sigi_output + " gff=" + sigi_gff_output)
                script.close()

            os.chmod(script_file.name, 0o755)
            script_file.file.close()

            with open(self.temp_dir_path + "/logs", 'w') as logs:
                subprocess.check_call(script_file.name, stdout=logs, cwd=self.SIGIHMM_PATH)
            script_file.close()

            report["sigi_gis"][str(embl_id)] = self.parse_sigi_gff(sigi_gff_output)

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)


class IslandPathPipelineComponent(PipelineComponent):
    name = "islandpath"
    dependencies = ["gbk_paths"]
    result_types = ["islandpath_gis"]
    output_dir = settings.BIO_APP_TEMP_DIR + "islandpath/"
    ISLANDPATH_PATH = settings.ISLANDPATH_PATH
    log_path = None
    temp_dir_path = None

    @staticmethod
    def parse_islandpath(islandpath_file):
        gi_list = []

        with open(islandpath_file, 'r') as islandpath:
            tsv_reader = csv.reader(islandpath, delimiter='\t')
            for row in tsv_reader:
                gi_list.append([row[1], row[2]])

        return gi_list

    def setup(self, report):
        self.temp_dir_path = self.output_dir + str(report["analysis"])
        os.mkdir(self.temp_dir_path, 0o777)

    def analysis(self, report):
        report["islandpath_gis"] = {}

        for gbk_id in report["gbk_paths"]:
            script_file = NamedTemporaryFile(delete=True)
            temp_path = self.temp_dir_path + "/" + str(gbk_id)

            with open(script_file.name, 'w') as script:
                script.write("#!/bin/bash\n")
                script.write(os.path.abspath(self.ISLANDPATH_PATH)
                             + " " + os.path.abspath(report["gbk_paths"][gbk_id])
                             + " " + os.path.abspath(temp_path))
                script.close()

            os.chmod(script_file.name, 0o755)
            script_file.file.close()

            self.log_path = self.temp_dir_path + "/logs_" + str(report["analysis"])
            with open(self.log_path, 'w') as logs:
                subprocess.check_call(script_file.name, stdout=logs, cwd=self.temp_dir_path)
            script_file.close()

            report["islandpath_gis"][gbk_id] = self.parse_islandpath(temp_path)

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)


class MergeIslandsPipelineComponent(PipelineComponent):
    name = "merge_gis"
    dependencies = ["islandpath_gis", "sigi_gis"]
    result_types = ["merge_gis"]
    output_dir = settings.BIO_APP_TEMP_DIR + "mash/"
    MASH_PATH = settings.MASH_PATH
    log_path = None
    temp_dir_path = None
    threshold = 500

    def merge_gi_list(self, first_list, second_list):
        merged_gis = first_list + second_list
        merged_gis.sort(key=lambda x: int(x[0]))
        output_gis = []

        current_gi = None
        for gi in merged_gis:
            if current_gi is None:
                current_gi = copy.deepcopy(gi)
            elif int(gi[0]) < (int(current_gi[1]) + self.threshold) and int(gi[1]) > int(current_gi[1]):
                current_gi[1] = gi[1]
            else:
                output_gis.append(current_gi)
                current_gi = copy.deepcopy(gi)

        if current_gi is not None:
            output_gis.append(current_gi)

        return output_gis

    def set_threshold(self, threshold):
        self.threshold = threshold

    def setup(self, report):
        self.temp_dir_path = self.output_dir + str(report["analysis"])
        os.mkdir(self.temp_dir_path, 0o777)

    def analysis(self, report):
        merged_gi_dict = dict()

        for genome_id in report["islandpath_gis"]:
            merged_gi_dict[genome_id] = self.merge_gi_list(report["islandpath_gis"][str(genome_id)],
                                                           report["sigi_gis"][str(genome_id)])

        report["merge_gis"] = merged_gi_dict

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)
