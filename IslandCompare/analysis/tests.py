from rest_framework.test import APIRequestFactory, force_authenticate, APITestCase
from django.contrib.auth.models import User
from analysis.models import Analysis, AnalysisComponent, AnalysisType
from genomes.models import Genome
from rest_framework.reverse import reverse
from analysis.views import AnalysisListView, AnalysisRunView
from django.utils import timezone
from django.core.files.uploadedfile import SimpleUploadedFile
from unittest import mock, skip
from analysis.serializers import ReportCsvSerializer, ReportVisualizationOverviewSerializer
import csv
import os


# Create your tests here.


class ListAnalysisTestCase(APITestCase):
    test_username = "test_user"
    test_user = None

    test_celery_task_id = "1"
    test_name = "test_analysis"
    test_genomes = None
    test_submit_time = timezone.now()
    test_analysis = None

    test_component_type_name = "test"
    test_component_type = None

    test_component_celery_task_id = "2"
    test_component = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_analysis = Analysis(celery_task_id=self.test_celery_task_id,
                                      name=self.test_name,
                                      submit_time=self.test_submit_time,
                                      owner=self.test_user)
        self.test_analysis.save()

        self.test_component_type = AnalysisType(name=self.test_component_type_name)
        self.test_component_type.save()

        self.test_component = AnalysisComponent(celery_task_id=self.test_component_celery_task_id,
                                                type=self.test_component_type,
                                                analysis=self.test_analysis)
        self.test_component.save()

        self.second_analysis = Analysis.objects.create(celery_task_id="2",
                                                       name="test_analysis_2",
                                                       submit_time=self.test_submit_time,
                                                       owner=self.test_user)

    def test_authenticated_list_analysis(self):
        url = reverse('analysis')

        request = self.factory.get(url)
        force_authenticate(request, self.test_user)
        response = AnalysisListView.as_view()(request)

        self.assertEqual(200, response.status_code)
        self.assertEqual(self.test_analysis.id, response.data[1]['id'])
        self.assertEqual(self.test_name, response.data[1]['name'])
        self.assertEqual(self.test_submit_time.strftime('%Y-%m-%dT%H:%M:%S.%fZ'), response.data[1]['submit_time'])
        self.assertEqual(Analysis.objects.all().count(), len(response.data))

    def test_authenticated_sublist_analysis(self):
        number_expected_arguments = 1
        url = reverse('analysis') + "?num=" + str(number_expected_arguments)

        request = self.factory.get(url)
        force_authenticate(request, self.test_user)
        response = AnalysisListView.as_view()(request)

        self.assertEqual(200, response.status_code)
        self.assertEqual(2, response.data[0]['id'])
        self.assertEqual(self.test_name + "_2", response.data[0]['name'])
        self.assertEqual(number_expected_arguments, len(response.data))

    def test_authenticated_list_analysis_components(self):
        url = reverse('analysis')

        request = self.factory.get(url)
        force_authenticate(request, self.test_user)
        response = AnalysisListView.as_view()(request)

        self.assertEqual(200, response.status_code)
        self.assertEqual(1, len(response.data[1]['analysiscomponent_set']))

    def test_unauthenticated_list_analysis(self):
        url = reverse('analysis')

        request = self.factory.get(url)
        response = AnalysisListView.as_view()(request)

        self.assertEqual(403, response.status_code)


class RetrieveAnalysisTestCase(APITestCase):
    test_username = "test_user"
    test_user = None

    test_celery_task_id = "1"
    test_name = "test_analysis"
    test_genomes = None
    test_submit_time = timezone.now()
    test_analysis = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_analysis = Analysis(celery_task_id=self.test_celery_task_id,
                                      name=self.test_name,
                                      submit_time=self.test_submit_time,
                                      owner=self.test_user)
        self.test_analysis.save()

    def test_authenticated_retrieve_analysis(self):
        url = reverse('analysis_details', kwargs={'pk': self.test_analysis.id})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.get(url)

        self.assertEqual(200, response.status_code)
        self.assertEqual(self.test_name, response.data['name'])
        self.assertEqual(self.test_submit_time.strftime('%Y-%m-%dT%H:%M:%S.%fZ'), response.data['submit_time'])

    def test_unauthenticated_retrieve_analysis(self):
        url = reverse('analysis_details', kwargs={'pk': self.test_analysis.id})

        response = self.client.get(url)

        self.assertEqual(403, response.status_code)

    def test_incorrect_user_retrieve_analysis(self):
        url = reverse('analysis_details', kwargs={'pk': self.test_analysis.id})

        incorrect_user = User(username="incorrect_user")
        incorrect_user.save()

        self.client.force_authenticate(user=incorrect_user)
        response = self.client.get(url)

        self.assertEqual(404, response.status_code)

    def test_retrieve_analysis_does_not_exist(self):
        url = reverse('analysis_details', kwargs={'pk': self.test_analysis.id+1})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.get(url)

        self.assertEqual(404, response.status_code)


class AnalysisResultsViewTestCase(APITestCase):
    test_username = "test_user"
    test_user = None

    test_analysis_name = "test_analysis"

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

    test_celery_task_id = "1"
    test_name = "test_analysis"
    test_submit_time = timezone.now()
    test_analysis = None

    report = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_2_gbk)

        self.test_analysis = Analysis.objects.create(celery_task_id=self.test_celery_task_id,
                                                     name=self.test_name,
                                                     submit_time=self.test_submit_time,
                                                     owner=self.test_user)
        self.test_analysis.genomes.add(self.test_genome_1)
        self.test_analysis.genomes.add(self.test_genome_2)
        self.test_analysis.save()

        self.report = {
            "analysis": 1,
            "available_dependencies": ["gbk_paths", "islandpath_gis"],
            "gbk_paths": {
                self.test_genome_1.id: self.test_genome_1.gbk,
                self.test_genome_2.id: self.test_genome_2.gbk,
            },
            "islandpath_gis": {
                self.test_genome_1.id: [
                    [0, 100], [110, 120]
                ],
                self.test_genome_2.id: [
                    [125, 135], [145, 155]
                ]
            }
        }

    @mock.patch('celery.result.AsyncResult.status', new_callable=mock.PropertyMock)
    @mock.patch('celery.result.AsyncResult.get')
    @skip("Figure out how to return task result in wrapper instead of task id of subtask, or rewrite test")
    def test_authenticated_completed_result_retrieval(self, mock_get, mock_status):
        mock_status.return_value = 'SUCCESS'
        mock_get.return_value = self.report
        url = reverse('analysis_results', kwargs={'pk': self.test_analysis.id})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.get(url)

        self.assertEqual(200, response.status_code)

        self.assertEqual(self.test_genome_1_name, response.data["genomes"][self.test_genome_1.id]["name"])
        self.assertEqual(self.test_genome_2_name, response.data["genomes"][self.test_genome_2.id]["name"])

    @mock.patch('celery.result.AsyncResult.status', new_callable=mock.PropertyMock)
    def test_authenticated_incomplete_analysis_retrieval(self, mock_status):
        mock_status.return_value = 'PENDING'
        url = reverse('analysis_results', kwargs={'pk': self.test_analysis.id})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.get(url)

        self.assertEqual(202, response.status_code)

    @mock.patch('celery.result.AsyncResult.status', new_callable=mock.PropertyMock)
    def test_authenticated_failed_analysis_retrieval(self, mock_status):
        mock_status.return_value = 'FAILURE'
        url = reverse('analysis_results', kwargs={'pk': self.test_analysis.id})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.get(url)

        self.assertEqual(500, response.status_code)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class ExportAnalysisResultTestCase(APITestCase):
    test_username = "test_user"
    test_user = None

    test_celery_task_id = "1"
    test_name = "test_analysis"
    test_genomes = None
    test_submit_time = timezone.now()
    test_analysis = None

    genome_1_id = 1
    genome_1_path = "mock_path_1"
    genome_2_id = 2
    genome_2_path = "mock_path_2"

    report = {
        "analysis": 1,
        "available_dependencies": ["gbk_paths", "islandpath_gis"],
        "gbk_paths": {
            genome_1_id: genome_1_path,
            genome_2_id: genome_2_path,
        },
        "islandpath_gis": {
            genome_1_id: [
                [0, 100], [110, 120]
            ],
            genome_2_id: [
                [125, 135], [145, 155]
            ]
        }
    }

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_analysis = Analysis(celery_task_id=self.test_celery_task_id,
                                      name=self.test_name,
                                      submit_time=self.test_submit_time,
                                      owner=self.test_user)
        self.test_analysis.save()

    @mock.patch('celery.result.AsyncResult.status', new_callable=mock.PropertyMock)
    @mock.patch('celery.result.AsyncResult.get')
    @skip("Figure out how to return task result in wrapper instead of task id of subtask, or rewrite test")
    def test_authenticated_completed_export_analysis(self, mock_get, mock_status):
        mock_status.return_value = 'SUCCESS'
        mock_get.return_value = self.report
        url = reverse('analysis_export', kwargs={'pk': self.test_analysis.id})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.get(url)

        self.assertEqual(200, response.status_code)

        reader = csv.reader(response.data.split("\n"))

        self.assertEqual(["IslandPath GIs"], next(reader))

        self.assertEqual([], next(reader))
        self.assertEqual([os.path.basename(self.genome_1_path)], next(reader))
        for i in self.report["islandpath_gis"][self.genome_1_id]:
            self.assertEqual(i, [int(_) for _ in next(reader)])

        self.assertEqual([], next(reader))
        self.assertEqual([os.path.basename(self.genome_2_path)], next(reader))
        for j in self.report["islandpath_gis"][self.genome_2_id]:
            self.assertEqual(j, [int(_) for _ in next(reader)])

    @mock.patch('celery.result.AsyncResult.status', new_callable=mock.PropertyMock)
    def test_authenticated_incomplete_export_analysis(self, mock_status):
        mock_status.return_value = 'PENDING'
        url = reverse('analysis_export', kwargs={'pk': self.test_analysis.id})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.get(url)

        self.assertEqual(202, response.status_code)

    @mock.patch('celery.result.AsyncResult.status', new_callable=mock.PropertyMock)
    def test_authenticated_failed_export_analysis(self, mock_status):
        mock_status.return_value = 'FAILURE'
        url = reverse('analysis_export', kwargs={'pk': self.test_analysis.id})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.get(url)

        self.assertEqual(500, response.status_code)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class UpdateAnalysisTestCase(APITestCase):
    test_username = "test_user"
    test_user = None

    test_celery_task_id = "1"
    test_name = "test_analysis"
    test_genomes = None
    test_submit_time = timezone.now()
    test_analysis = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_analysis = Analysis(celery_task_id=self.test_celery_task_id,
                                      name=self.test_name,
                                      submit_time=self.test_submit_time,
                                      owner=self.test_user)
        self.test_analysis.save()

    def test_authenticated_update_analysis(self):
        updated_name = "updated_name"

        url = reverse('analysis_details', kwargs={'pk': self.test_analysis.id})

        self.client.force_authenticate(user=self.test_user)
        response = self.client.put(url,
                                   {'name': updated_name})

        updated_analysis = Analysis.objects.get(id=self.test_analysis.id)

        self.assertEqual(200, response.status_code)
        self.assertEqual(updated_name, updated_analysis.name)

    def test_unauthenticated_update_analysis(self):
        updated_name = "updated_name"

        url = reverse('analysis_details', kwargs={'pk': self.test_analysis.id})

        response = self.client.put(url,
                                   {'name': updated_name})

        updated_analysis = Analysis.objects.get(id=self.test_analysis.id)

        self.assertEqual(403, response.status_code)
        self.assertEqual(self.test_name, updated_analysis.name)


class RunAnalysisTestCase(APITestCase):
    test_username = "test_user"
    test_user = None

    test_analysis_name = "test_analysis"

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

    expected_components = ["start_pipeline",
                           "setup_gbk",
                           "gbk_metadata",
                           "parsnp",
                           "mauve",
                           "sigi",
                           "islandpath",
                           "merge_gis",
                           "end_pipeline"]

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_2_gbk)

    def test_authenticated_run_analysis(self):
        url = reverse('analysis_run')

        request = self.factory.post(url,
                                    {'name': self.test_analysis_name,
                                     'genomes': [self.test_genome_1.id, self.test_genome_2.id]})
        force_authenticate(request, self.test_user)

        run_pipeline_callback = mock.MagicMock()
        response = AnalysisRunView.as_view(run_pipeline_callback=run_pipeline_callback)(request)

        self.assertEqual(200, response.status_code)
        self.assertTrue(Analysis.objects.filter(id=response.data['id']).exists())

        analysis = Analysis.objects.get(id=response.data['id'])
        self.assertEqual(self.test_analysis_name, analysis.name)
        self.assertTrue(analysis.genomes.filter(id=self.test_genome_1.id).exists())
        self.assertTrue(analysis.genomes.filter(id=self.test_genome_2.id).exists())

        analysis_components = analysis.analysiscomponent_set
        self.assertTrue(set(self.expected_components) == set([_.type.name for _ in analysis_components.all()]))

        run_pipeline_callback.assert_called_once()

    def test_incorrect_owner_run_analysis(self):
        url = reverse('analysis_run')

        incorrect_user = User(username="incorrect_user")
        incorrect_user.save()

        request = self.factory.post(url,
                                    {'genomes': [self.test_genome_1.id, self.test_genome_2.id]})
        force_authenticate(request, incorrect_user)

        response = AnalysisRunView.as_view()(request)

        self.assertEqual(400, response.status_code)

    def test_only_1_genome_run_analysis(self):
        url = reverse('analysis_run')

        request = self.factory.post(url,
                                    {'genomes': [self.test_genome_1.id]})
        force_authenticate(request, self.test_user)

        response = AnalysisRunView.as_view()(request)

        self.assertEqual(400, response.status_code)

    def test_unauthenticated_run_analysis(self):
        url = reverse('analysis_run')

        request = self.factory.post(url,
                                    {'genomes': [self.test_genome_1.id, self.test_genome_2.id]})

        response = AnalysisRunView.as_view()(request)

        self.assertEqual(403, response.status_code)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()


class ReportToCsvSerializerTestCase(APITestCase):
    genome_1_id = 1
    genome_1_path = "mock_path_1"
    genome_2_id = 2
    genome_2_path = "mock_path_2"

    report = {
        "analysis": 1,
        "available_dependencies": ["gbk_paths", "islandpath_gis"],
        "gbk_paths": {
            genome_1_id: genome_1_path,
            genome_2_id: genome_2_path,
        },
        "islandpath_gis": {
            genome_1_id: [
                [0, 100], [110, 120]
            ],
            genome_2_id: [
                [125, 135], [145, 155]
            ]
        }
    }

    def test_report_csv_serializer(self):
        serializer = ReportCsvSerializer(self.report)
        output = serializer.data

        reader = csv.reader(output.split("\n"))

        self.assertEqual(["IslandPath GIs"], next(reader))

        self.assertEqual([], next(reader))
        self.assertEqual([os.path.basename(self.genome_1_path)], next(reader))
        for i in self.report["islandpath_gis"][self.genome_1_id]:
            self.assertEqual(i, [int(_) for _ in next(reader)])

        self.assertEqual([], next(reader))
        self.assertEqual([os.path.basename(self.genome_2_path)], next(reader))
        for j in self.report["islandpath_gis"][self.genome_2_id]:
            self.assertEqual(j, [int(_) for _ in next(reader)])


class ReportVisualizationOverviewTestCase(APITestCase):
    test_username = "test_user"
    test_user = None

    test_analysis_name = "test_analysis"

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

    test_celery_task_id = "1"
    test_name = "test_analysis"
    test_submit_time = timezone.now()
    test_analysis = None

    report = None

    def setUp(self):
        self.factory = APIRequestFactory()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_1_gbk)

        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user,
                                                   gbk=self.test_genome_2_gbk)

        self.test_analysis = Analysis.objects.create(celery_task_id=self.test_celery_task_id,
                                                     name=self.test_name,
                                                     submit_time=self.test_submit_time,
                                                     owner=self.test_user)
        self.test_analysis.genomes.add(self.test_genome_1)
        self.test_analysis.genomes.add(self.test_genome_2)
        self.test_analysis.save()

        self.report = {
            "analysis": self.test_analysis.id,
            "available_dependencies": ["gbk_paths", "islandpath_gis"],
            "gbk_paths": {
                self.test_genome_1.id: self.test_genome_1.gbk,
                self.test_genome_2.id: self.test_genome_2.gbk,
            },
            "islandpath_gis": {
                self.test_genome_1.id: [
                    [0, 100], [110, 120]
                ],
                self.test_genome_2.id: [
                    [125, 135], [145, 155]
                ]
            }
        }

    @skip("Fix test for new implementation")
    def test_report_visualization_overview_serializer(self):
        serializer = ReportVisualizationOverviewSerializer(self.report)
        output = serializer.data

        self.assertEqual(self.test_genome_1_name, output["genomes"][self.test_genome_1.id]["name"])
        self.assertEqual(self.test_genome_2_name, output["genomes"][self.test_genome_2.id]["name"])

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()
