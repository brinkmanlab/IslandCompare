from analysis.pipeline import PipelineComponent
from analysis.models import Analysis, Genome
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
from Bio.SeqRecord import SeqRecord
import numpy as np
from analysis.lib.mcl_clustering import mcl
import json

class StartPipelineComponent(PipelineComponent):
    """
    Pipeline component that sets the start time of the analysis in the database.
    """
    name = "start_pipeline"

    def analysis(self, report):
        analysis_entry = Analysis.objects.get(id=report['analysis'])
        analysis_entry.start_time = datetime.now()
        analysis_entry.save()


class EndPipelineComponent(PipelineComponent):
    """
    Pipeline component that sets the end time of the analysis in the database.
    """
    name = "end_pipeline"

    def analysis(self, report):
        analysis_entry = Analysis.objects.get(id=report['analysis'])
        analysis_entry.complete_time = datetime.now()
        analysis_entry.save()


class SetupGbkPipelineComponent(PipelineComponent):
    """
    Pipeline component that adds the gbk paths to the report.
    """
    name = "setup_gbk"
    result_types = ["gbk_paths"]

    def analysis(self, report):
        analysis_entry = Analysis.objects.get(id=report['analysis'])
        genomes = analysis_entry.genomes

        report['gbk_paths'] = dict()
        for genome in genomes.all():
            report['gbk_paths'][str(genome.id)] = genome.gbk.path


class GbkMetadataComponent(PipelineComponent):
    """
    Pipeline component that adds the genome length to the report.
    """
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

class RGIPipelineComponent(PipelineComponent):
    name = "rgi"
    dependencies = ["gbk_paths"]
    result_types = ["amr_genes"]
    output_dir = settings.BIO_APP_TEMP_DIR + "rgi/"
    temp_dir_path = None
    temp_results_dir = None
    fna_files = {}
    RGI_PATH = settings.RGI_PATH

    @staticmethod
    def parse_rgi_json(json_file):
        json_dict = json.load(json_file)
        amr_genes = []
        for key in json_dict:
            for entry_key in json_dict[key]:
                entry = json_dict[key][entry_key]
                if type(entry) is dict and "orf_start" in entry:
                    amr_genes.append({k : entry.get(k) for k in ("orf_start", "orf_end", "orf_strand", "type_match")})
        # Keep only unique entries
        amr_genes = [dict(y) for y in set(tuple(x.items()) for x in amr_genes)]
        amr_sorted = sorted(amr_genes, key=lambda gene: gene['orf_start'])
        return amr_sorted

    def setup(self, report):
        # Create FASTA files from GenBank Files for use by RGI
        self.fna_files = {}
        self.temp_dir_path = mkdtemp()
        for gbk_id in report["gbk_paths"]:
            gbk_path = report["gbk_paths"][gbk_id]
            self.fna_files[gbk_id] = os.path.abspath(self.temp_dir_path+"/"+str(gbk_id)+".fna")
            ParsnpPipelineComponent.convert_gbk_to_fna(gbk_path, self.fna_files[gbk_id])

    def analysis(self, report):
        script_file = NamedTemporaryFile(delete=True)
        self.temp_results_dir = self.output_dir + str(report["analysis"])
        os.mkdir(self.temp_results_dir, 0o777)
        report["amr_genes"] = {}

        for fna_id in self.fna_files:
            rgi_output = self.temp_results_dir + "/" + str(fna_id) # RGI adds extensions
            with open(script_file.name, 'w') as script:
                script.write("#!/bin/bash\n")
                script.write("python " + self.RGI_PATH + " -i " + self.fna_files[fna_id] + " -o " + rgi_output)
                script.close()
            os.chmod(script_file.name, 0o755)
            script_file.file.close()
            # Run RGI
            with open(self.temp_results_dir + "/logs", 'w') as logs:
                subprocess.check_call(script_file.name, stdout=logs)
            script_file.close()
            # Parse RGI results and add to report
            with open(rgi_output + ".json", "r") as output:
                report["amr_genes"][str(fna_id)] = self.parse_rgi_json(output)

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)
        if self.temp_results_dir is not None and os.path.exists(self.temp_results_dir):
            rmtree(self.temp_results_dir)


class ParsnpPipelineComponent(PipelineComponent):
    """
    Pipeline component that runs and adds parsnp data to the report.
    """
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


class UserNewickPipelineComponent(PipelineComponent):
    """
    Pipeline component that adds user newick data to the report.
    The user newick file should be set before running the analysis.
    """
    name = "user_newick"
    dependencies = ["gbk_paths"]
    result_types = ["newick"]

    def set_newick(self, user_newick):
        self.logger.info("Set User Newick as:\n{}".format(user_newick))
        self.param['user_file_contents'] = user_newick

    def analysis(self, report):
        selected_genomes = Genome.objects.filter(id__in=[int(_) for _ in report["gbk_paths"].keys()])

        tree = Phylo.read(StringIO(self.param['user_file_contents']), 'newick')
        terminals = tree.get_terminals()

        for leaf in terminals:
            selected_genome = selected_genomes.filter(name__exact=leaf.name)

            leaf.name = str(selected_genome.get().id)

        processed_tree = StringIO()
        Phylo.write(tree, processed_tree, "newick")
        processed_tree.seek(0)
        output = processed_tree.getvalue()
        processed_tree.close()
        report["newick"] = output


class UserGIPipelineComponent(PipelineComponent):
    """
    Pipeline component that adds user gis to the report.
    The user gi file should be set before running the analysis.
    """
    name = "user_gi"
    dependencies = ["gbk_paths"]
    result_types = ["user_gis"]

    def set_gi(self, user_gi):
        self.logger.info("Set User GI as:\n{}".format(user_gi))
        self.param['user_file_contents'] = user_gi

    @staticmethod
    # Create a dict with contains keys as filename and value as list of genomic islands {start,end}
    def parse_gi_file(gifile):
        genomeDict = dict()
        gireader = csv.reader(gifile.read().splitlines(), dialect=csv.excel_tab)
        for row in gireader:
            genomeName = row[0]
            giStart = row[1]
            giEnd = row[2]
            # if genome name is not in genomeDict then add it
            if genomeName not in genomeDict:
                genomeDict[genomeName] = list()
            # add start and end of current genome list in genomedict
            if len(row)>3:
                giColor = row[3]
                genomeDict[genomeName].append({'start': giStart, 'end': giEnd, 'color': giColor})
            else:
                genomeDict[genomeName].append({'start': giStart, 'end': giEnd})
        return genomeDict

    def analysis(self, report):
        selected_genomes = Genome.objects.filter(id__in=[int(_) for _ in report["gbk_paths"].keys()])
        output_dict = dict()

        gi_dict = self.parse_gi_file(StringIO(self.param['user_file_contents']))

        for key in gi_dict.keys():
            selected_genome = selected_genomes.filter(name__exact=key)
            key_id = selected_genome.get().id
            output_dict[str(key_id)] = gi_dict[key]

        report["user_gis"] = output_dict


class MauvePipelineComponent(PipelineComponent):
    """
    Pipeline component that runs and adds mauve data to the report.
    """
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
    """
    Pipeline component that runs and adds sigi-hmm data for each genome in the report to the report.
    """
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
    """
    Pipeline component that runs and adds islandpath data to the report.
    """
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
    """
    Pipeline component that merges islands and adds the merged islands to the report.
    """
    name = "merge_gis"
    dependencies = ["islandpath_gis", "sigi_gis"]
    result_types = ["merge_gis"]
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

    def analysis(self, report):
        merged_gi_dict = dict()

        for genome_id in report["islandpath_gis"]:
            merged_gi_dict[genome_id] = self.merge_gi_list(report["islandpath_gis"][str(genome_id)],
                                                           report["sigi_gis"][str(genome_id)])

        report["merge_gis"] = merged_gi_dict


class MashMclClusterPipelineComponent(PipelineComponent):
    """
    Pipeline component that runs mash and clusters gis.
    """
    name = "mash_mcl"
    dependencies = ["gbk_paths", "merge_gis"]
    result_types = ["cluster_gis"]
    output_dir = settings.BIO_APP_TEMP_DIR + "mash/"
    MASH_PATH = settings.MASH_PATH
    temp_dir_path = None
    fna_dir_path = None

    def create_compound_sketch(self, fastaFileList, outputFileName):
        scriptFile = NamedTemporaryFile(delete=False)

        with open(scriptFile.name, 'w') as script:
            script.write("#!/bin/bash\n")
            script.write(self.MASH_PATH + " sketch -o " + outputFileName + " ")
            for fastaFile in fastaFileList:
                script.write(fastaFile+" ")
            script.close()

        os.chmod(scriptFile.name, 0o0755)
        scriptFile.file.close()

        self.logger.info("Running MASH script: {}".format(scriptFile.name))
        subprocess.check_call(scriptFile.name)
        scriptFile.close()

    def calculate_mash_distance(self, referenceFile, queryFastaFile, outputFile):
        scriptFile = NamedTemporaryFile(delete=True)

        with open(scriptFile.name, 'w') as script:
            script.write("#!/bin/bash\n")
            script.write(self.MASH_PATH + " dist " + referenceFile + " " + queryFastaFile + " > " + outputFile)
            script.close()

        os.chmod(scriptFile.name, 0o0755)
        scriptFile.file.close()
        subprocess.check_call(scriptFile.name)
        scriptFile.close()

    def parse_mash_output(self, outputFile):
        outputList = []
        with open(outputFile, 'r') as output:
            reader = csv.reader(output, delimiter='\t')
            for row in reader:
                outputList.append({'referenceId': row[0],
                                   'queryId': row[1],
                                   'mashDistance': row[2],
                                   'pValue': row[3],
                                   'matchingHashes': row[4]})
        return outputList

    def getSubsequence(self, genbankFile, startPosition, endPosition, islandNumber, description=None):
        record_dict = SeqIO.index(genbankFile, "genbank")
        sequenceName = list(record_dict.keys())[0]
        if description is not None:
            return SeqRecord(record_dict[sequenceName].seq[int(startPosition):int(endPosition)], id=sequenceName + "-" + str(islandNumber), description=description)
        else:
            return SeqRecord(record_dict[sequenceName].seq[int(startPosition):int(endPosition)], id=sequenceName + "-" + str(islandNumber))

    def writeFastaFile(self, outputFileName, seqRecordList):
        with open(outputFileName, 'w') as outputFileHandle:
            SeqIO.write(seqRecordList, outputFileHandle, "fasta")

    def setup(self, report):
        self.temp_dir_path = self.output_dir + str(report["analysis"])
        os.mkdir(self.temp_dir_path, 0o777)

        self.fna_dir_path = self.temp_dir_path + "/fna"
        os.mkdir(self.fna_dir_path)

    def create_gi_fasta_files(self, report):
        genome_list = report["gbk_paths"]
        island_path_list = []

        for genome_id in genome_list.keys():
            genome_fna_path = self.fna_dir_path + "/" + str(genome_id)
            os.mkdir(genome_fna_path)

            gi_counter = 0
            for gi in report['merge_gis'][genome_id]:
                self.logger.info("Adding GI: " + str(genome_id) + "-" + str(gi_counter))
                entrySequence = self.getSubsequence(report['gbk_paths'][genome_id], gi[0], gi[1], gi_counter)
                self.writeFastaFile(self.fna_dir_path + "/" + str(genome_id) + "/" + str(gi_counter), entrySequence)
                island_path_list.append(self.fna_dir_path + "/" + str(genome_id) + "/" + str(gi_counter))
                gi_counter += 1

        return island_path_list

    def analysis(self, report):
        island_path_list = self.create_gi_fasta_files(report)
        self.create_compound_sketch(island_path_list, self.temp_dir_path + "/compoundScratch")

        distance_matrix = []

        for island in island_path_list:
            self.logger.info("Processing island: " + island)
            if os.path.isfile(self.temp_dir_path + "/output"):
                os.remove(self.temp_dir_path + "/output")
            self.calculate_mash_distance(self.temp_dir_path + "/compoundScratch.msh", island, self.temp_dir_path + "/output")
            distance_matrix.append([float(i['mashDistance']) for i in self.parse_mash_output(self.temp_dir_path + "/output")])

        baseMatrix = []
        for i in range(len(island_path_list)):
            nextRow = []
            for j in range(len(island_path_list)):
                nextRow.append(1)
            baseMatrix.append(nextRow)
        numpyBaseMatrix = np.array(baseMatrix)
        numpyDistanceMatrix = np.array(distance_matrix)

        mclAdjacencyMatrix = np.subtract(numpyBaseMatrix, numpyDistanceMatrix)
        np.set_printoptions(threshold='nan')

        M, clusters = mcl(mclAdjacencyMatrix)
        outputList = {}

        for sequenceId in report["gbk_paths"].keys():
            outputList[str(sequenceId)] = {}
        islandIdList = list([i for i in range(len(island_path_list))])

        numberClusters = 0

        while len(islandIdList) > 0:
            numberClusters += 1
            currentCluster = clusters[islandIdList[0]]
            for i in currentCluster:
                self.logger.info("Assigning clusters for cluster {}".format(numberClusters))
                island = island_path_list[i]
                splitIsland = island.split('/')[-2:]
                currentSequenceId = splitIsland[0]
                islandId = splitIsland[1]
                self.logger.info("Current Sequence Id: " + str(currentSequenceId))
                self.logger.info("Current Island Id: " + str(islandId))
                outputList[str(currentSequenceId)][str(islandId)] = numberClusters - 1
            remainingIslands = filter(lambda x: x not in currentCluster, islandIdList)
            islandIdList = list(remainingIslands)

        outputList['numberClusters'] = numberClusters - 1
        report["cluster_gis"] = outputList

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)
