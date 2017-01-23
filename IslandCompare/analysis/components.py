from analysis.pipeline import PipelineComponent
from analysis.models import Analysis
from tempfile import mkdtemp
from shutil import rmtree
from Bio import SeqIO
import os


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

    def setup(self, report):
        self.temp_dir_path = mkdtemp()
        for gbk_path in report["gbk_paths"]:
            self.convert_gbk_to_fna(gbk_path,
                                    self.temp_dir_path+"/"+(os.path.splitext(gbk_path)[0]).split("/")[-1]+".fna")

    def analysis(self, report):
        self.setup(report)
        # TODO: Run PARSNP and put result in report
        self.cleanup()

    def cleanup(self):
        rmtree(self.temp_dir_path)
