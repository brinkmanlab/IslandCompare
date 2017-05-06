from django.test import TestCase, override_settings
from genomes.models import Genome
from django.contrib.auth.models import User
from django.core.files import File
from analysis.components import ParsnpPipelineComponent, MauvePipelineComponent, \
    SigiHMMPipelineComponent, IslandPathPipelineComponent
from Bio import Phylo
from io import StringIO
from rest_framework.reverse import reverse
from rest_framework.test import force_authenticate, APITestCase, APIRequestFactory
from analysis.views import AnalysisRunView


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
        component.run(report)

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
        component.run(report)

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
        component.run(report)

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
        component.run(report)

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
    def test_valid_run_view_request(self):
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
