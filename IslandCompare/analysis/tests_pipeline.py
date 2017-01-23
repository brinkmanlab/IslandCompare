from django.test import TestCase, mock, override_settings
from analysis.pipeline import Pipeline, PipelineComponent, PipelineSerializer
from django.contrib.auth.models import User
from analysis.models import Analysis, AnalysisComponent, AnalysisType
from genomes.models import Genome
from django.core.files.uploadedfile import SimpleUploadedFile
from analysis.components import SetupGbkPipelineComponent
from analysis.tasks import run_pipeline_wrapper


class PipelineComponentStub(PipelineComponent):
    def analysis(self, report):
        return


class PipelineTestCase(TestCase):
    pipeline = None

    test_username = "test_username"
    test_user = None

    test_genome_1 = None
    test_genome_1_name = "genome_1"

    test_genome_2 = None
    test_genome_2_name = "genome_2"

    def setUp(self):
        self.pipeline = Pipeline()

        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_genome_1 = Genome.objects.create(name=self.test_genome_1_name,
                                                   owner=self.test_user)
        self.test_genome_2 = Genome.objects.create(name=self.test_genome_2_name,
                                                   owner=self.test_user)

    def test_pipeline_constructor_values_empty(self):
        self.assertEqual([], self.pipeline.pipeline_components)
        self.assertEqual([], self.pipeline.available_dependencies)
        self.assertEqual(None, self.pipeline.analysis)

    def test_pipeline_multiple_objects(self):
        pipeline_component_1 = mock.MagicMock()

        pipeline_1 = Pipeline()
        pipeline_1.pipeline_components = [pipeline_component_1]

        pipeline_2 = Pipeline()

        self.assertEqual([pipeline_component_1], pipeline_1.pipeline_components)
        self.assertEqual([], pipeline_2.pipeline_components)
        self.assertNotEqual(pipeline_1.pipeline_components, pipeline_2.pipeline_components)

    def test_pipeline_create_database_entry(self):
        name = "test_pipeline"
        genomes = Genome.objects.filter(owner=self.test_user)
        user = self.test_user
        component_mock = mock.MagicMock()
        component_mock.create_database_entry = mock.MagicMock()
        self.pipeline.pipeline_components = [component_mock]

        self.pipeline.create_database_entry(name, genomes, user)

        self.assertTrue(Analysis.objects.filter(name=name,
                                                genomes=genomes,
                                                owner=user).exists())
        component_mock.create_database_entry.assert_called_once()

    def test_pipeline_append_component_available_dependencies(self):
        new_dependency = "new_dependency"

        component_mock = mock.MagicMock()
        component_mock.result_types = [new_dependency]

        self.pipeline.append_component(component_mock)
        self.assertTrue(set(component_mock.result_types).issubset(self.pipeline.available_dependencies))

    def tearDown(self):
        for analysis in Analysis.objects.all():
            analysis.delete()

        for genome in Genome.objects.all():
            genome.delete()

        self.test_user.delete()


class PipelineComponentTestCase(TestCase):
    test_username = "test_username"
    test_user = None

    test_dependency = "test_dependency"

    test_analysis_name = "test_analysis"
    test_analysis = None

    def setUp(self):
        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_analysis = Analysis.objects.create(name=self.test_analysis_name,
                                                     owner=self.test_user)

    def test_pipeline_component_validate_dependencies(self):
        pipeline_component = PipelineComponentStub()
        pipeline_component.dependencies = [self.test_dependency]
        self.assertTrue(pipeline_component.validate_dependencies([self.test_dependency]))

    def test_pipeline_component_missing_dependency(self):
        pipeline_component = PipelineComponentStub()
        pipeline_component.dependencies = [self.test_dependency]
        self.assertFalse(pipeline_component.validate_dependencies([]))

    def test_pipeline_component_missing_dependency_raise_exception(self):
        pipeline_component = PipelineComponentStub()
        pipeline_component.dependencies = [self.test_dependency]
        with self.assertRaises(RuntimeError):
            pipeline_component.validate_dependencies([], raise_exception=True)

    def test_pipeline_component_create_database_entry(self):
        analysis_component_name = "test_analysis"

        pipeline_component = PipelineComponentStub()
        pipeline_component.name = analysis_component_name
        pipeline_component.create_database_entry(self.test_analysis)

        self.assertTrue(AnalysisType.objects.filter(name=analysis_component_name).exists())
        self.assertTrue(AnalysisComponent.objects.filter(type__name=analysis_component_name,
                                                         analysis=self.test_analysis).exists())

    def test_pipeline_component_run_all_steps(self):
        setup = mock.MagicMock()
        analysis = mock.MagicMock()
        cleanup = mock.MagicMock()

        pipeline_component = PipelineComponentStub()
        pipeline_component.setup = setup
        pipeline_component.analysis = analysis
        pipeline_component.cleanup = cleanup

        pipeline_component.run(dict())

        setup.assert_called_once()
        analysis.assert_called_once()
        cleanup.assert_called_once()

    def tearDown(self):
        self.test_user.delete()


class PipelineSerializerTestCase(TestCase):
    pipeline = None

    test_available_dependencies = ['dep_1', 'dep_2']
    test_component = PipelineComponentStub()
    test_component.name = "component_1"
    test_pipeline_components = [test_component]

    test_username = "test_username"
    test_user = None

    test_analysis = None
    test_analysis_name = "test_analysis"

    def setUp(self):
        self.test_user = User(username=self.test_username)
        self.test_user.save()

        self.test_analysis = Analysis.objects.create(name=self.test_analysis_name,
                                                     owner=self.test_user)

        self.pipeline = Pipeline()
        self.pipeline.available_dependencies = self.test_available_dependencies
        self.pipeline.pipeline_components = self.test_pipeline_components
        self.pipeline.analysis = self.test_analysis

    def test_serialize_pipeline(self):
        json_dict = PipelineSerializer.serialize(self.pipeline)

        self.assertEqual(self.test_available_dependencies, json_dict['available_dependencies'])
        self.assertTrue(self.test_component.name in json_dict['pipeline_components'])
        self.assertEqual(self.test_analysis.id, json_dict['analysis'])

    def test_deserialize_pipeline(self):
        json_dict = {
            'available_dependencies': self.test_available_dependencies,
            'pipeline_components': [self.test_component.name],
            'analysis': self.test_analysis.id,
        }

        mock_factory = mock.MagicMock()
        mock_factory.create_component = mock.MagicMock(return_value=self.test_component)

        deserialized_pipeline = PipelineSerializer.deserialize(json_dict, mock_factory)

        mock_factory.create_component.assert_called_with(self.test_component.name)
        self.assertEqual(self.test_available_dependencies, deserialized_pipeline.available_dependencies)
        self.assertEqual(self.test_pipeline_components, deserialized_pipeline.pipeline_components)
        self.assertEqual(self.test_analysis, deserialized_pipeline.analysis)


class PipelineTasksTestCase(TestCase):
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
    def test_run_gbk_component_pipeline(self):
        pipeline_name = "pipeline"
        genomes = Genome.objects.filter(owner=self.test_user)

        pipeline = Pipeline()
        pipeline.append_component(SetupGbkPipelineComponent())
        pipeline.create_database_entry(pipeline_name, genomes, self.test_user)

        celery_result = run_pipeline_wrapper(pipeline)
        result = celery_result.get()

        self.assertEqual(self.test_genome_1.gbk.path, result['gbk_paths'][self.test_genome_1.id])
        self.assertEqual(self.test_genome_2.gbk.path, result['gbk_paths'][self.test_genome_2.id])
        self.assertTrue('gbk_paths' in result['available_dependencies'])
        self.assertTrue('setup_gbk' in result['pipeline_components'])

        analysis = Analysis.objects.get(id=pipeline.analysis.id)
        analysis_component = AnalysisComponent.objects.get(type__name='setup_gbk', analysis=pipeline.analysis)

        self.assertIsNotNone(analysis.celery_task_id)
        self.assertIsNotNone(analysis_component.celery_task_id)
        self.assertNotEqual(analysis.celery_task_id, analysis_component.celery_task_id)

    def tearDown(self):
        for genome in Genome.objects.all():
            genome.delete()
