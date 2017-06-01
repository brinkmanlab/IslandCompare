from __future__ import absolute_import, unicode_literals
from celery import shared_task
from analysis.models import AnalysisComponent
from analysis import components
from analysis.pipeline import PipelineSerializer
import logging


class PipelineComponentFactory(object):
    @staticmethod
    def create_component(name, param=None):
        if name == "setup_gbk":
            return components.SetupGbkPipelineComponent(param)
        if name == "gbk_metadata":
            return components.GbkMetadataComponent(param)
        if name == "parsnp":
            return components.ParsnpPipelineComponent(param)
        if name == "mauve":
            return components.MauvePipelineComponent(param)
        if name == "sigi":
            return components.SigiHMMPipelineComponent(param)
        if name == "islandpath":
            return components.IslandPathPipelineComponent(param)
        if name == "merge_gis":
            return components.MergeIslandsPipelineComponent(param)
        if name == "start_pipeline":
            return components.StartPipelineComponent(param)
        if name == "end_pipeline":
            return components.EndPipelineComponent(param)
        if name == "mash_mcl":
            return components.MashMclClusterPipelineComponent(param)
        if name == "user_newick":
            return components.UserNewickPipelineComponent(param)
        if name == "rgi":
            return components.RGIPipelineComponent(param)
        raise(RuntimeError("Given component does not exist: {}".format(name)))


logger = logging.getLogger(__name__)


def run_pipeline_wrapper(self, pipeline):
    task = run_pipeline.s(PipelineSerializer.serialize(pipeline)).apply_async()
    logger.info("Job scheduled with celery task id: {}".format(task.task_id))
    return task


@shared_task(bind=True)
def run_pipeline(self, serialized_pipeline):
    logger.info("Begin scheduling of pipeline components with celery task id: {}".format(self.request.id))
    pipeline = PipelineSerializer.deserialize(serialized_pipeline, PipelineComponentFactory())
    pipeline.analysis.celery_task_id = self.request.id
    pipeline.analysis.save()

    index = len(pipeline.pipeline_components) - 1
    next_component = None

    while index != 0:
        if next_component is None:
            next_component = run_pipeline_component.s(pipeline.analysis.id,
                                                      pipeline.pipeline_components[index].name,
                                                      pipeline.pipeline_components[index].param)
        else:
            current_component = run_pipeline_component.s(pipeline.analysis.id,
                                                         pipeline.pipeline_components[index].name,
                                                         pipeline.pipeline_components[index].param)
            current_component.link(next_component)
            next_component = current_component
        index -= 1

    task_chain = run_pipeline_component.s(serialized_pipeline,
                                          pipeline.analysis.id,
                                          pipeline.pipeline_components[index].name)
    task_chain.link(next_component)

    logger.info("End scheduling of pipeline components with celery task id: {}".format(self.request.id))
    return task_chain.apply_async()


@shared_task(bind=True)
def run_pipeline_component(self, report, analysis_id, pipeline_component_name, pipeline_components_param=None):
    pipeline_component = PipelineComponentFactory.create_component(pipeline_component_name, pipeline_components_param)
    logger.info("Attempting to run component: {} for analysis with id: {}".format(pipeline_component_name,
                                                                                   analysis_id))

    analysis_component = AnalysisComponent.objects.get(analysis__id=analysis_id,
                                                       type__name=pipeline_component.name)
    analysis_component.celery_task_id = self.request.id
    analysis_component.save()

    return pipeline_component.run(report)
