from analysis.pipeline import PipelineComponent
from analysis.models import Analysis, Genome
from genomes.models import Gene, GenomicIsland, UserGenomicIsland
from tempfile import mkdtemp, NamedTemporaryFile
from shutil import rmtree
from Bio import SeqIO, Phylo
import os
from django.conf import settings
import subprocess
from io import StringIO
import csv
from datetime import datetime
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
    Pipeline component that adds the gbk paths to the report and genes to the database.
    """
    name = "setup_gbk"
    result_types = ["gbk_paths"]

    @staticmethod
    def create_genes(genome):
        with open(genome.gbk.path) as genbank_handle:
            record = SeqIO.read(genbank_handle, "genbank")
            last_gene = None
            for feature in record.features:
                if feature.type in ["gene", "CDS", "tRNA", "rRNA"]:
                    start = feature.location.start
                    end = feature.location.end
                    # Check if the current gene is in the same position as the last gene.
                    if last_gene is not None and start == last_gene.start and end == last_gene.end:
                        if feature.type is not "gene":
                            # we don't want to replace the more specific types with "gene".
                            last_gene.type = feature.type
                        gene = last_gene
                    else:
                        gene = Gene(
                            type=feature.type,
                            gene="",
                            locus_tag="",
                            product="",
                            start=start,
                            end=end,
                            strand=feature.location.strand,
                            genome=genome
                        )
                        last_gene = gene
                    if gene.gene == "":
                        try:
                            gene.gene = feature.qualifiers['gene'][0]
                        except KeyError:
                            pass
                    if gene.locus_tag == "":
                        try:
                            gene.locus_tag = feature.qualifiers['locus_tag'][0]
                        except KeyError:
                            pass
                    if gene.product == "":
                        try:
                            gene.product = feature.qualifiers['product'][0]
                        except KeyError:
                            pass
                    gene.save()

    def analysis(self, report):
        analysis_entry = Analysis.objects.get(id=report['analysis'])
        genomes = analysis_entry.genomes

        report['gbk_paths'] = dict()
        for genome in genomes.all():
            report['gbk_paths'][str(genome.id)] = genome.gbk.path
            if genome.gene_set.exists():
                self.logger.info("{} gene set found".format(genome.name))
            else:
                self.logger.info("Creating gene set for genome {}".format(genome.name))
                self.create_genes(genome)


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
    """
    Pipeline component that runs and adds RGI data to the report
    """
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
                    amr_genes.append({k : entry.get(k) for k in ("orf_start", "orf_end", "orf_strand")})
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
    def parse_newick(input_path, num_genomes):
        tree = Phylo.read(input_path, 'newick')
        terminals = tree.get_terminals()

        # Parsnp may fail to include all genomes in the newick in specific cases
        assert (len(terminals) == num_genomes), "Newick produced by Parsnp does not contain all genomes"

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

        try:
            # parse_newick also validates that newick contains all genomes
            report["newick"] = self.parse_newick(self.temp_results_dir + "/parsnp.tree", len(report["gbk_paths"]))
        except AssertionError:
            # temp dir with incorrect newick file is deleted so parsnp can remake on retry
            if self.temp_results_dir is not None and os.path.exists(self.temp_results_dir):
                rmtree(self.temp_results_dir)
            raise

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
        selected_genomes = Genome.objects.filter(id__in=report["gbk_paths"].keys())

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
        analysis = Analysis.objects.get(id=report["analysis"])
        selected_genomes = Genome.objects.filter(id__in=report["gbk_paths"].keys())

        gi_dict = self.parse_gi_file(StringIO(self.param['user_file_contents']))

        for key in gi_dict.keys():
            selected_genome = selected_genomes.get(name=key)
            for island in gi_dict[key]:
                if 'color' in island:
                    color = island['color']
                else:
                    color = ""
                UserGenomicIsland(method="user",
                                  start=island['start'],
                                  end=island['end'],
                                  genome=selected_genome,
                                  analysis=analysis,
                                  color=color).save()


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
    Pipeline component that runs and adds sigi-hmm data for each genome in the report to the database.
    """
    name = "sigi"
    dependencies = ["gbk_paths"]
    output_dir = settings.BIO_APP_TEMP_DIR + "sigi/"
    temp_dir_path = None
    embl_files = {}
    SIGIHMM_PATH = settings.SIGIHMM_PATH
    SIGIHMM_EXE = settings.SIGIHMM_EXE
    exception = None

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
                    # Some results contain lines where start is 1 / end is 0. Skip to avoid erroneous GIs
                    if gene_dict[3] == '1' or gene_dict[4] == '0':
                        continue
                    if gene_dict[2] == 'PUTAL':
                        # At the start of a genomic island, set start and end of possible genomic island
                        if not island_flag:
                            start = gene_dict[3]
                            end = gene_dict[4]
                            island_flag = True
                        # Continuation of current genomic island, change end and continue
                        else:
                            end = gene_dict[4]
                    # End of genomic island, append current start and end to list
                    elif island_flag:
                        gi_list.append([int(start), int(end)])
                        island_flag = False
            # For cases where the last line is part of a genomic island
            if island_flag:
                gi_list.append([int(start), int(end)])
        return gi_list

    def setup(self, report):
        self.temp_dir_path = self.output_dir + str(report["analysis"])
        os.mkdir(self.temp_dir_path, 0o777)

        for gbk_id in report["gbk_paths"]:
            genome = Genome.objects.get(id=gbk_id)
            if genome.genomicisland_set.filter(method="sigi").exists():
                self.logger.info("{} Sigi-HMM genomic islands found".format(genome.name))
            else:
                basename = os.path.splitext(os.path.basename(report["gbk_paths"][gbk_id]))[0]
                with open(self.output_dir + str(report["analysis"]) + "/" + basename + ".embl", 'w') as embl_output:
                    SeqIO.convert(report["gbk_paths"][gbk_id], "genbank", embl_output, "embl")
                    self.embl_files[gbk_id] = os.path.abspath(embl_output.name)

    def analysis(self, report):
        script_file = NamedTemporaryFile(delete=True)

        for embl_id in self.embl_files:
            genome = Genome.objects.get(id=embl_id)
            sigi_gff_output = os.path.abspath(self.temp_dir_path + "/" + str(embl_id) + ".gff")
            sigi_output = os.path.abspath(self.temp_dir_path +"/" + str(embl_id) + ".embl")
            with open(script_file.name, 'w') as script:
                script.write("#!/bin/bash\n")
                script.write("/usr/bin/java " + self.SIGIHMM_EXE + " input=" +
                             self.embl_files[embl_id] + " output=" + sigi_output + " gff=" + sigi_gff_output)
                script.close()

            os.chmod(script_file.name, 0o755)
            script_file.file.close()

            try:
                with open(self.temp_dir_path + "/logs", 'w') as logs:
                    subprocess.check_call(script_file.name, stdout=logs, cwd=self.SIGIHMM_PATH)
                script_file.close()

                for sigi_gi in self.parse_sigi_gff(sigi_gff_output):
                    GenomicIsland(method="sigi",
                                  start=sigi_gi[0],
                                  end=sigi_gi[1],
                                  genome=genome
                                  ).save()
            except subprocess.CalledProcessError as err:
                self.logger.info("Sigi-HMM Failed! analysis {}, genome {}".format(report["analysis"], embl_id))
                self.exception = err
                if "sigi" not in report["failed_components"]:
                    report["failed_components"]["sigi"] = [genome.name]
                else:
                    report["failed_components"]["sigi"].append(genome.name)

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)
        if self.exception:
            raise self.exception


class IslandPathPipelineComponent(PipelineComponent):
    """
    Pipeline component that runs and adds islandpath data to the database.
    """
    name = "islandpath"
    dependencies = ["gbk_paths"]
    output_dir = settings.BIO_APP_TEMP_DIR + "islandpath/"
    ISLANDPATH_PATH = settings.ISLANDPATH_PATH
    log_path = None
    temp_dir_path = None
    exception = None

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

        for gbk_id in report["gbk_paths"]:
            genome = Genome.objects.get(id=gbk_id)
            if genome.genomicisland_set.filter(method="islandpath").exists():
                self.logger.info("{} IslandPath genomic islands found".format(genome.name))
            else:
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
                try:
                    with open(self.log_path, 'w') as logs:
                        subprocess.check_call(script_file.name, stdout=logs, cwd=self.temp_dir_path)
                    script_file.close()

                    for islandpath_gi in self.parse_islandpath(temp_path):
                        GenomicIsland(method="islandpath",
                                      start=islandpath_gi[0],
                                      end=islandpath_gi[1],
                                      genome=genome
                                      ).save()
                except subprocess.CalledProcessError as err:
                    self.logger.info("IslandPath Failed! analysis {}, genome {}".format(report["analysis"], gbk_id))
                    self.exception = err
                    if "islandpath" not in report["failed_components"]:
                        report["failed_components"]["islandpath"] = [genome.name]
                    else:
                        report["failed_components"]["islandpath"].append(genome.name)

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)
        if self.exception:
            raise self.exception


class MergeIslandsPipelineComponent(PipelineComponent):
    """
    Pipeline component that merges islands and adds the merged islands to the database.
    """
    name = "merge_gis"
    threshold = 500

    def merge_gi_queryset(self, gi_queryset):
        ordered_gis = gi_queryset.order_by('start')

        current_gi = ordered_gis[0]
        current_gi.pk = None # Creates a copy of the gi object
        for i in range(1, len(ordered_gis)):
            next_gi = ordered_gis[i]
            if current_gi.end + self.threshold >= next_gi.start:
                if current_gi.end < next_gi.end:
                    current_gi.end = next_gi.end
            else:
                current_gi.method = "merge"
                current_gi.save() # Will create a new entry in the database because pk was set to None
                current_gi = next_gi
                current_gi.pk = None
        current_gi.method = "merge"
        current_gi.save()

    def set_threshold(self, threshold):
        self.threshold = threshold

    def analysis(self, report):
        genomes = Genome.objects.filter(id__in=report['gbk_paths'].keys())
        for genome in genomes:
            gis = genome.genomicisland_set
            if gis.filter(method__in=["sigi", "islandpath"]).exists() and not gis.filter(method="merge").exists():
                self.merge_gi_queryset(gis.all())


class MashMclClusterPipelineComponent(PipelineComponent):
    """
    Pipeline component that runs mash and clusters gis.
    """
    name = "mash_mcl"
    dependencies = ["gbk_paths"]
    result_types = ["numberClusters"]
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

    def parse_mash_output(self, output_file):
        with open(output_file, 'r') as output:
            reader = csv.reader(output, delimiter='\t')
            # The third column contains the mash distance
            return [row[2] for row in reader]

    def get_subsequence(self, record, startPosition, endPosition, islandNumber, description=None):
        if description is not None:
            return SeqRecord(record.seq[int(startPosition):int(endPosition)], id=record.id + "-" + str(islandNumber), description=description)
        else:
            return SeqRecord(record.seq[int(startPosition):int(endPosition)], id=record.id + "-" + str(islandNumber))

    def writeFastaFile(self, outputFileName, seqRecordList):
        with open(outputFileName, 'w') as outputFileHandle:
            SeqIO.write(seqRecordList, outputFileHandle, "fasta")

    def create_gi_fasta_files(self, report):
        island_path_list = []
        if "user_gi" in report["pipeline_components"]:
            method = "user"
        else:
            method = "merge"

        for genome_id in report["gbk_paths"]:
            record = SeqIO.read(report['gbk_paths'][genome_id], "genbank")
            genome_fna_path = self.fna_dir_path + "/" + str(genome_id)
            os.mkdir(genome_fna_path)

            gi_counter = 0
            gis = Genome.objects.get(id=genome_id).genomicisland_set.filter(method=method)
            if method == "user":
                gis = gis.filter(usergenomicisland__analysis__id=report['analysis'])
            for gi in gis.all():
                self.logger.info("Adding GI: " + str(genome_id) + "-" + str(gi_counter))
                entrySequence = self.get_subsequence(record, gi.start, gi.end, gi_counter)
                self.writeFastaFile(self.fna_dir_path + "/" + str(genome_id) + "/" + str(gi_counter), entrySequence)
                island_path_list.append(self.fna_dir_path + "/" + str(genome_id) + "/" + str(gi_counter))
                gi_counter += 1

        return island_path_list

    def apply_cutoff(self, matrix, cutoff=0.96):
        n = len(matrix)
        for i in range(n - 1):
            for j in range(i + 1, n):
                if matrix[i][j] < cutoff:
                    matrix[i][j] = matrix[j][i] = 0
        return matrix

    def setup(self, report):
        self.temp_dir_path = self.output_dir + str(report["analysis"])
        os.mkdir(self.temp_dir_path, 0o777)

        self.fna_dir_path = self.temp_dir_path + "/fna"
        os.mkdir(self.fna_dir_path)

    def analysis(self, report):
        island_path_list = self.create_gi_fasta_files(report)
        if len(island_path_list) > 0:
            self.create_compound_sketch(island_path_list, self.temp_dir_path + "/compoundScratch")

            # Calculates the distance [0,1] from each genomic island to each genomic island
            distance_matrix = []
            for island in island_path_list:
                self.logger.info("Processing island: " + island)
                if os.path.isfile(self.temp_dir_path + "/output"):
                    os.remove(self.temp_dir_path + "/output")
                self.calculate_mash_distance(self.temp_dir_path + "/compoundScratch.msh", island, self.temp_dir_path + "/output")
                distance_matrix.append([float(i) for i in self.parse_mash_output(self.temp_dir_path + "/output")])

            numpy_distance_matrix = np.array(distance_matrix)
            # Convert to adjacency matrix by setting each value as 1 - value
            mcl_adjacency_matrix = np.vectorize(lambda i: 1 - i)(numpy_distance_matrix)
            np.set_printoptions(threshold='nan')

            # The cutoff ensures only very closely matches islands will be placed in a cluster together
            mcl_adjacency_matrix = self.apply_cutoff(mcl_adjacency_matrix)

            # mcl computes genomic island clusters
            M, raw_clusters = mcl(mcl_adjacency_matrix)

            processed_clusters = {str(genome_id): {} for genome_id in report["gbk_paths"]}
            cluster_number = 0
            for island in raw_clusters:
                advance_flag = 0 # Only advance the cluster_number if it is used
                for cluster_mate in raw_clusters[island]:
                    # island_path_list to get the genome ID and genome-specific GI number from the unspecific GI number
                    genome_id, island_num = island_path_list[cluster_mate].split("/")[-2:]
                    # Only add GIs to the same cluster if they both list each other in their cluster
                    if island in raw_clusters[cluster_mate] and island_num not in processed_clusters[genome_id]:
                        processed_clusters[genome_id][island_num] = cluster_number
                        advance_flag = 1
                cluster_number += advance_flag

            report['numberClusters'] = cluster_number
            analysis = Analysis.objects.get(id=report['analysis'])
            analysis.clusters = str(processed_clusters)
            analysis.save()

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)
