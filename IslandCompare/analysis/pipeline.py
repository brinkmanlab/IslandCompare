import abc
from analysis.models import AnalysisType, AnalysisComponent, Analysis
import logging
from datetime import datetime
import json


class PipelineComponent(abc.ABC):
    logger = logging.getLogger(__name__)

    name = ""
    dependencies = []
    result_types = []

    def __init__(self, param=None):
        if param is None:
            self.param = dict()
        else:
            self.param = param
            self.logger.info("Param set to: \n{}".format(json.dumps(self.param)))

    def create_database_entry(self, analysis):
        analysis_type, analysis_exists = AnalysisType.objects.get_or_create(name=self.name)
        AnalysisComponent.objects.create(type=analysis_type,
                                         analysis=analysis)
        self.logger.info("Created Analysis Component: {} for Analysis with Id: {}".format(analysis_type.name, analysis.id))

    def record_start_time(self, report):
        type = AnalysisType.objects.get(name=self.name)
        component = AnalysisComponent.objects.get(analysis__id__exact=report['analysis'],
                                                  type=type)
        component.start_time = datetime.now()
        component.save()

    def record_complete_time(self, report):
        type = AnalysisType.objects.get(name=self.name)
        component = AnalysisComponent.objects.get(analysis__id__exact=report['analysis'],
                                                  type=type)
        component.complete_time = datetime.now()
        component.save()

    def validate_dependencies(self, available_dependencies, raise_exception=False):
        for dependency in self.dependencies:
            if dependency not in available_dependencies:
                if raise_exception:
                    raise(RuntimeError("Missing dependency: {}".format(dependency)))
                else:
                    return False
        return True

    def setup(self, report):
        return

    @abc.abstractmethod
    def analysis(self, report):
        return

    def cleanup(self):
        return

    def run(self, report):
        self.logger.info("Begin Running Component: {} for Analysis with Id: {}".format(self.name, report['analysis']))
        self.record_start_time(report)
        self.setup(report)
        self.analysis(report)
        self.cleanup()
        self.record_complete_time(report)
        self.logger.info("End Running Component: {} for Analysis with Id: {}".format(self.name, report['analysis']))
        return report


class Pipeline(object):
    logger = logging.getLogger(__name__)

    pipeline_components = None
    available_dependencies = None
    analysis = None

    def __init__(self):
        self.pipeline_components = []
        self.available_dependencies = []
        self.analysis = None

    def create_database_entry(self, name, genomes, owner):
        self.analysis = Analysis.objects.create(name=name,
                                                owner=owner)
        self.analysis.genomes = genomes
        self.analysis.save()
        self.logger.info("Created Database Entry with Id: {}".format(self.analysis.id))

        for component in self.pipeline_components:
            component.create_database_entry(self.analysis)

    def append_component(self, pipeline_component):
        pipeline_component.validate_dependencies(available_dependencies=self.available_dependencies,
                                                 raise_exception=True)
        self.pipeline_components.append(pipeline_component)
        self.available_dependencies = self.available_dependencies + pipeline_component.result_types


class PipelineSerializer(object):
    @staticmethod
    def serialize(pipeline):
        serialized_pipeline_components = [component.name for component in pipeline.pipeline_components]
        serialized_pipeline_param = {component.name: component.param for component in pipeline.pipeline_components}
        serialized_dict = {
            'analysis': pipeline.analysis.id,
            'available_dependencies': pipeline.available_dependencies,
            'pipeline_components': serialized_pipeline_components,
            'component_param': serialized_pipeline_param
        }
        return serialized_dict

    @staticmethod
    def deserialize(serialized_dict, pipeline_component_factory):
        pipeline = Pipeline()

        for component in serialized_dict['pipeline_components']:
            pipeline.pipeline_components.append(
                pipeline_component_factory.create_component(component, serialized_dict['component_param'][component])
            )
        pipeline.available_dependencies = serialized_dict['available_dependencies']
        pipeline.analysis = Analysis.objects.get(id=serialized_dict['analysis'])

        return pipeline
