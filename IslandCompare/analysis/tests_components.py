from django.test import TestCase
from django.core.files.uploadedfile import SimpleUploadedFile
from django.contrib.auth.models import User
from genomes.models import Genome
from analysis.components import SetupGbkPipelineComponent
from analysis.pipeline import Pipeline, PipelineSerializer


class GbkComponentTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_genome_1 = None
    test_genome_1_name = "genome_1"
    test_genome_1_gbk_name = "test.gbk"
    test_genome_1_gbk_contents = bytes("test", 'utf-8')
    test_genome_1_gbk = SimpleUploadedFile(test_genome_1_gbk_name, test_genome_1_gbk_contents)

    test_genome_2 = None
    test_genome_2_name = "genome_2"
    test_genome_2_gbk_name = "test_2.gbk"
    test_genome_2_gbk_contents = bytes("test2", 'utf-8')
    test_genome_2_gbk = SimpleUploadedFile(test_genome_1_gbk_name, test_genome_1_gbk_contents)

    serialized_pipeline = None
    component = SetupGbkPipelineComponent()

    def setUp(self):
        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_2_gbk)

        pipeline_name = "pipeline"
        genomes = Genome.objects.filter(owner=self.test_user)

        pipeline = Pipeline()
        pipeline.append_component(self.component)
        pipeline.create_database_entry(pipeline_name, genomes, self.test_user)
        self.serialized_pipeline = PipelineSerializer.serialize(pipeline)

    def test_run_gbk_component(self):
        result = self.component.run(self.serialized_pipeline)

        self.assertEqual(self.test_genome_1.gbk.path, result['gbk_paths'][self.test_genome_1.id])
        self.assertEqual(self.test_genome_2.gbk.path, result['gbk_paths'][self.test_genome_2.id])
        self.assertTrue('gbk_paths' in result['available_dependencies'])
        self.assertTrue('setup_gbk' in result['pipeline_components'])

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()
