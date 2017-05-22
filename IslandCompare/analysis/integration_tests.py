from django.test import TestCase, override_settings
from genomes.models import Genome
from django.contrib.auth.models import User
from django.core.files import File
from analysis.components import ParsnpPipelineComponent, MauvePipelineComponent, \
    SigiHMMPipelineComponent, IslandPathPipelineComponent, MashMclClusterPipelineComponent
from Bio import Phylo
from io import StringIO
from rest_framework.reverse import reverse
from rest_framework.test import force_authenticate, APITestCase, APIRequestFactory
from analysis.views import AnalysisRunView
from unittest import mock
import os


class ParsnpComponentIntegrationTestCase(TestCase):
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

    def test_run_parsnp_component(self):
        report = {
            "analysis": 1,
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                self.test_genome_1.id: self.test_genome_1.gbk.path,
                self.test_genome_2.id: self.test_genome_2.gbk.path,
            },
        }
        component = ParsnpPipelineComponent()
        component.setup(report)
        component.analysis(report)
        component.cleanup()

        self.assertTrue("newick" in report)

        tree = Phylo.read(StringIO(report["newick"]), "newick")
        terminal_names = [clade.name for clade in tree.get_terminals()]

        self.assertTrue(str(self.test_genome_1.id) in terminal_names)
        self.assertTrue(str(self.test_genome_2.id) in terminal_names)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class MauveComponentIntegrationTestCase(TestCase):
    test_username = "username"
    test_user = None

    test_genome_1 = None
    test_genome_1_name = "genome_1"
    test_genome_1_gbk = File(open("../TestFiles/AE009952.gbk"))

    test_genome_2 = None
    test_genome_2_name = "genome_2"
    test_genome_2_gbk = File(open("../TestFiles/BX936398.gbk"))

    test_genome_2 = None
    test_genome_2_name = "genome_3"
    test_genome_2_gbk = File(open("../TestFiles/CP000305.gbk"))

    def setUp(self):
        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_2_gbk)

    def test_run_single_pair_mauve_component(self):
        report = {
            "analysis": 1,
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                str(self.test_genome_1.id): self.test_genome_1.gbk.path,
                str(self.test_genome_2.id): self.test_genome_2.gbk.path,
            },
            "newick": "({}:200.10871,{}:200.10871):0.00000;\n".format(self.test_genome_1.id,
                                                                      self.test_genome_2.id),
        }

        component = MauvePipelineComponent()
        component.setup(report)
        component.analysis(report)
        component.cleanup()

        self.assertTrue("alignment" in report)

    def test_run_multiple_pair_mauve_component(self):
        report = {
            "analysis": 1,
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                str(self.test_genome_1.id): self.test_genome_1.gbk.path,
                str(self.test_genome_2.id): self.test_genome_2.gbk.path,
            },
            "newick": "({}:200.10871,{}:200.10871):0.00000;\n".format(self.test_genome_1.id,
                                                                      self.test_genome_2.id),
        }

        component = MauvePipelineComponent()
        component.setup(report)
        component.analysis(report)
        component.cleanup()

        self.assertTrue("alignment" in report)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class SigiHMMComponentIntegrationTestCase(TestCase):
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

    def test_sigihmm_analysis(self):
        report = {
            "analysis": 1,
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                str(self.test_genome_1.id): self.test_genome_1.gbk.path,
                str(self.test_genome_2.id): self.test_genome_2.gbk.path,
            },
        }

        component = SigiHMMPipelineComponent()
        component.setup(report)
        component.analysis(report)
        component.cleanup()

        self.assertTrue("sigi_gis" in report)
        self.assertEqual(2, len(report["sigi_gis"]))

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class IslandPathComponentIntegrationTestCase(TestCase):
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

    def test_islandpath_analysis(self):
        report = {
            "analysis": 1,
            "available_dependencies": "gbk_paths",
            "gbk_paths": {
                self.test_genome_1.id: self.test_genome_1.gbk.path,
                self.test_genome_2.id: self.test_genome_2.gbk.path,
            },
        }

        component = IslandPathPipelineComponent()
        component.setup(report)
        component.analysis(report)
        component.cleanup()

        self.assertTrue("islandpath_gis" in report)
        self.assertEqual(2, len(report["islandpath_gis"]))
        self.assertTrue(self.test_genome_1.id in report["islandpath_gis"])
        self.assertTrue(self.test_genome_2.id in report["islandpath_gis"])

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class AnalysisRunViewIntegrationTestCase(APITestCase):
    test_username = "username"
    test_user = None

    factory = APIRequestFactory()
    test_analysis_name = "test_analysis"

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

    @override_settings(CELERY_TASK_ALWAYS_EAGER=True)
    @mock.patch('analysis.serializers.AnalysisComponentSerializer.get_status')
    def test_valid_run_view_request(self, mock_get):
        mock_get.return_value = 'PENDING'

        url = reverse('analysis_run')

        request = self.factory.post(url,
                                    {'name': self.test_analysis_name,
                                     'genomes': [self.test_genome_1.id, self.test_genome_2.id]})
        force_authenticate(request, self.test_user)

        response = AnalysisRunView.as_view()(request)

        self.assertEqual(200, response.status_code)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class MashMCLTestCase(TestCase):
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

    def test_generate_gi_fna(self):
        component = MashMclClusterPipelineComponent()

        report = {
            "analysis": 1,
            "available_dependencies": ["merge_gis", "gbk_paths"],
            "gbk_paths": {self.test_genome_1.id: self.test_genome_1.gbk.path,
                          self.test_genome_2.id: self.test_genome_2.gbk.path},
            "merge_gis": {self.test_genome_1.id: [[0, 1000]],
                          self.test_genome_2.id: [[0, 1000]]}
        }

        component.setup(report)
        component.create_gi_fasta_files(report)

        self.assertTrue(os.path.exists(component.temp_dir_path + "/fna/1/0"))
        self.assertTrue(os.path.exists(component.temp_dir_path + "/fna/2/0"))

        component.cleanup()

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()