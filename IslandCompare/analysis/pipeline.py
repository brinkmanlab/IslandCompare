import abc
from analysis.models import AnalysisType, AnalysisComponent, Analysis
import logging
from datetime import datetime
import json


class PipelineComponent(abc.ABC):
    """
    Abstract pipeline component used by the pipeline.
    """
    logger = logging.getLogger(__name__)

    name = ""
    """The name of the component to be defined in the inheriting class."""
    dependencies = []
    """The required fields in the report for the component to work."""
    result_types = []
    """The fields added to the report by the component."""

    def __init__(self, param=None):
        """
        Constructor for the pipeline component
        :param param: Any parameters (dictionary) needed by the pipeline component
        """
        if param is None:
            self.param = dict()
        else:
            self.param = param
            self.logger.info("Param set to: \n{}".format(json.dumps(self.param)))

    def create_database_entry(self, analysis):
        """
        Creates the database entry.
        :param analysis:
        :return:
        """
        analysis_type, analysis_exists = AnalysisType.objects.get_or_create(name=self.name)
        AnalysisComponent.objects.create(type=analysis_type,
                                         analysis=analysis)
        self.logger.info("Created Analysis Component: {} for Analysis with Id: {}".format(analysis_type.name, analysis.id))

    def record_start_time(self, report):
        """
        Records the start time in the database.
        :param report:
        :return:
        """
        type = AnalysisType.objects.get(name=self.name)
        component = AnalysisComponent.objects.get(analysis__id__exact=report['analysis'],
                                                  type=type)
        component.start_time = datetime.now()
        component.save()

    def record_complete_time(self, report):
        """
        Records the complete time in the database
        :param report:
        :return:
        """
        type = AnalysisType.objects.get(name=self.name)
        component = AnalysisComponent.objects.get(analysis__id__exact=report['analysis'],
                                                  type=type)
        component.complete_time = datetime.now()
        component.save()

    def validate_dependencies(self, available_dependencies, raise_exception=False):
        """
        Validates the components dependencies against a list of available fields
        :param available_dependencies: a list of available fields
        :param raise_exception: raise a runtime error if set to true
        :return: true if dependencies satisfied, false otherwise
        """
        for dependency in self.dependencies:
            if dependency not in available_dependencies:
                if raise_exception:
                    raise(RuntimeError("Missing dependency: {}".format(dependency)))
                else:
                    return False
        return True

    def setup(self, report):
        """
        Method to contain any setup required by the component
        To be defined in inheriting class
        :param report:
        :return:
        """
        return

    @abc.abstractmethod
    def analysis(self, report):
        """
        Method that runs and adds analysis data to the report.
        Fields defined in result_types should be added to the report here
        :param report:
        :return:
        """
        return

    def cleanup(self):
        """
        Method to contain any cleanup required by the component
        :return:
        """
        return

    def run(self, report):
        """
        Runs the component
        :param report:
        :return:
        """
        self.logger.info("Begin Running Component: {} for Analysis with Id: {}".format(self.name, report['analysis']))
        self.record_start_time(report)
        self.setup(report)
        self.analysis(report)
        self.cleanup()
        self.record_complete_time(report)
        self.logger.info("End Running Component: {} for Analysis with Id: {}".format(self.name, report['analysis']))
        return report


class Pipeline(object):
    """
    The pipeline that runs pipeline components. Set to currently run components in FIFO.
    """
    logger = logging.getLogger(__name__)

    pipeline_components = None
    """A list of pipeline components"""
    available_dependencies = None
    """A list of available dependencies"""
    analysis = None
    """The analysis"""
    failed_components = None
    """A list of failed components"""

    def __init__(self):
        """
        Constructor
        """
        self.pipeline_components = []
        self.available_dependencies = []
        self.analysis = None
        self.failed_components = {}

    def create_database_entry(self, name, genomes, owner):
        """
        Creates a database entry for the pipeline
        :param name: name of the analysis
        :param genomes: the genomes in teh analysis
        :param owner: the owner of the analysis
        :return:
        """
        self.analysis = Analysis.objects.create(name=name,
                                                owner=owner)
        self.analysis.genomes = genomes
        self.analysis.save()
        self.logger.info("Created Database Entry with Id: {}".format(self.analysis.id))

        for component in self.pipeline_components:
            component.create_database_entry(self.analysis)

    def append_component(self, pipeline_component):
        """
        Add a pipeline component to the pipeline
        :param pipeline_component:
        :return:
        """
        pipeline_component.validate_dependencies(available_dependencies=self.available_dependencies,
                                                 raise_exception=True)
        self.pipeline_components.append(pipeline_component)
        self.available_dependencies = self.available_dependencies + pipeline_component.result_types


class PipelineSerializer(object):
    """
    Serializer for a pipeline and all components that it contains.
    """
    @staticmethod
    def serialize(pipeline):
        """
        Serializes the pipeline as a dict
        :param pipeline:
        :return: the serialized pipeline
        """
        serialized_pipeline_components = [component.name for component in pipeline.pipeline_components]
        serialized_pipeline_param = {component.name: component.param for component in pipeline.pipeline_components}
        serialized_dict = {
            'analysis': pipeline.analysis.id,
            'available_dependencies': pipeline.available_dependencies,
            'pipeline_components': serialized_pipeline_components,
            'component_param': serialized_pipeline_param,
            'failed_components': pipeline.failed_components
        }
        return serialized_dict

    @staticmethod
    def deserialize(serialized_dict, pipeline_component_factory):
        """
        Deserializes the pipeline into objects
        :param serialized_dict:
        :param pipeline_component_factory: a factory that takes in pipeline component names and returns their object
        :return: the pipeline object
        """
        pipeline = Pipeline()

        for component in serialized_dict['pipeline_components']:
            pipeline.pipeline_components.append(
                pipeline_component_factory.create_component(component, serialized_dict['component_param'][component])
            )
        pipeline.available_dependencies = serialized_dict['available_dependencies']
        pipeline.analysis = Analysis.objects.get(id=serialized_dict['analysis'])
        pipeline.failed_components = serialized_dict['failed_components']

        return pipeline
