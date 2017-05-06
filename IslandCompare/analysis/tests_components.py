from django.test import TestCase, mock
from django.core.files.uploadedfile import SimpleUploadedFile
from django.contrib.auth.models import User
from genomes.models import Genome
from analysis.components import SetupGbkPipelineComponent, ParsnpPipelineComponent, MauvePipelineComponent, \
    SigiHMMPipelineComponent, GbkMetadataComponent
from analysis.pipeline import Pipeline, PipelineSerializer
from django.core.files import File
import filecmp
import os


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


class ParsnpComponentTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_genome_1 = None
    test_genome_1_name = "genome_1"
    test_genome_1_gbk = File(open("../TestFiles/AE009952.gbk"))

    test_genome_2 = None
    test_genome_2_name = "genome_2"
    test_genome_2_gbk = File(open("../TestFiles/BX936398.gbk"))

    def setUp(self):
        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_2_gbk)

    def test_setup_cleanup_parsnp_component(self):
        report = {
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                self.test_genome_1.id: self.test_genome_1.gbk.path,
                self.test_genome_2.id: self.test_genome_2.gbk.path,
            },
        }
        component = ParsnpPipelineComponent()
        component.setup(report)

        expected_fna_1_path = (component.temp_dir_path +
                               "/" +
                               str(self.test_genome_1.id) +
                               ".fna")
        expected_fna_2_path = (component.temp_dir_path +
                               "/" +
                               str(self.test_genome_2.id) +
                               ".fna")

        self.assertIsNotNone(component.temp_dir_path)
        self.assertTrue(os.path.isdir(component.temp_dir_path))
        self.assertTrue(os.path.isfile(expected_fna_1_path))
        self.assertTrue(os.path.isfile(expected_fna_2_path))
        self.assertTrue(filecmp.cmp("../TestFiles/AE009952.fna", expected_fna_1_path))
        self.assertTrue(filecmp.cmp("../TestFiles/BX936398.fna", expected_fna_2_path))

        component.cleanup()

        self.assertFalse(os.path.isdir(component.temp_dir_path))

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class MauveComponentTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_genome_1 = None
    test_genome_1_name = "genome_1"
    test_genome_1_gbk = File(open("../TestFiles/AE009952.gbk"))

    test_genome_2 = None
    test_genome_2_name = "genome_2"
    test_genome_2_gbk = File(open("../TestFiles/BX936398.gbk"))

    def setUp(self):
        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_2_gbk)

    def test_mauve_setup_sub_processes(self):
        report = {
            "analysis": 1,
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                self.test_genome_1.id: self.test_genome_1.gbk.path,
                self.test_genome_2.id: self.test_genome_2.gbk.path,
            },
            "newick": "({}:200.10871,{}:200.10871):0.00000;\n".format(self.test_genome_1.id,
                                                                      self.test_genome_2.id),
        }

        component = MauvePipelineComponent()
        mauve_subprocess_mock = mock.MagicMock()
        component.mauve_subprocess = mauve_subprocess_mock
        retrieve_mauve_mock = mock.MagicMock()
        component.retrieve_mauve_results = retrieve_mauve_mock

        component.run(report)

        mauve_subprocess_mock.assert_called_once()
        retrieve_mauve_mock.assert_called_once()
        arguments, kwargs = mauve_subprocess_mock.call_args

        self.assertEqual(component.temp_dir_path + "/" + str(0) + "/", arguments[0])
        self.assertEqual({self.test_genome_1.gbk.path, self.test_genome_2.gbk.path},
                         set(arguments[1]))

        component.cleanup()

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class SigiHMMComponentTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_genome_1 = None
    test_genome_1_name = "genome_1"
    test_genome_1_gbk = File(open("../TestFiles/AE009952.gbk"))

    test_genome_2 = None
    test_genome_2_name = "genome_2"
    test_genome_2_gbk = File(open("../TestFiles/BX936398.gbk"))

    def setUp(self):
        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_2_gbk)

    def test_sigihmm_setup(self):
        report = {
            "analysis": 1,
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                self.test_genome_1.id: self.test_genome_1.gbk.path,
                self.test_genome_2.id: self.test_genome_2.gbk.path,
            },
        }

        component = SigiHMMPipelineComponent()

        component.setup(report)

        embl_files = component.embl_files
        self.assertEqual(2, len(embl_files))
        for genome_id in embl_files:
            self.assertTrue(os.path.exists(embl_files[genome_id]))

        component.cleanup()

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class GbkMetadataTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_genome_1 = None
    test_genome_1_name = "genome_1"
    test_genome_1_gbk = File(open("../TestFiles/AE009952.gbk"))
    test_genome_1_size = 4600755

    report = None

    def setUp(self):
        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.report = {
            "analysis": 1,
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                self.test_genome_1.id: self.test_genome_1.gbk.path,
            },
        }

    def test_get_genome_size(self):
        component = GbkMetadataComponent()

        component.run(self.report)

        self.assertEqual(self.test_genome_1_size, self.report["gbk_metadata"][self.test_genome_1.id]["size"])

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()
