from analysis.pipeline import PipelineComponent
from analysis.models import Analysis
from tempfile import mkdtemp, NamedTemporaryFile
from shutil import rmtree
from Bio import SeqIO, Phylo
import os
from django.conf import settings
import subprocess


class SetupGbkPipelineComponent(PipelineComponent):
    name = "setup_gbk"
    result_types = ["gbk_paths"]

    def analysis(self, report):
        analysis_entry = Analysis.objects.get(id=report['analysis'])
        genomes = analysis_entry.genomes

        report['gbk_paths'] = dict()
        for genome in genomes.all():
            report['gbk_paths'][genome.id] = genome.gbk.path


class ParsnpPipelineComponent(PipelineComponent):
    name = "parsnp"
    dependencies = ["gbk_paths"]
    result_types = ["newick"]
    temp_dir_path = None
    output_dir = "temp/parsnp/"
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

    def read_newick(self, input_path):
        tree = Phylo.read(input_path, 'newick')
        return self.parse_tree(tree.root)

    def parse_tree(self, node):
        current_node = {}
        name_split = str(node.name).split(".")
        fna_counter = 0

        for subname in name_split:
            if "fna" not in subname:
                fna_counter += 1
            else:
                break

        current_node['genome_id'] = ".".join(name_split[0:fna_counter])
        current_node['length'] = node.branch_length
        if len(node.clades) > 0:
            current_node['children'] = []
            for childNode in node.clades:
                current_node['children'].append(self.parse_tree(childNode))
        return current_node

    def setup(self, report):
        self.temp_dir_path = mkdtemp()
        for gbk_id in report["gbk_paths"]:
            gbk_path = report["gbk_paths"][gbk_id]
            self.convert_gbk_to_fna(gbk_path,
                                    self.temp_dir_path+"/"+str(gbk_id)+".fna")

    def analysis(self, report):
        self.setup(report)

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

        report["newick"] = self.read_newick(self.temp_results_dir + "/parsnp.tree")

        self.cleanup()

    def cleanup(self):
        if self.temp_dir_path is not None and os.path.exists(self.temp_dir_path):
            rmtree(self.temp_dir_path)
        if self.temp_results_dir is not None and os.path.exists(self.temp_results_dir):
            rmtree(self.temp_results_dir)
