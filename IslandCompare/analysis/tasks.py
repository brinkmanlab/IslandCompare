from __future__ import absolute_import, unicode_literals
from celery import shared_task, chain
from analysis.models import AnalysisComponent
from analysis import components
from analysis.pipeline import PipelineSerializer
import logging


class PipelineComponentFactory(object):
    """
    Pipeline component factory that converts names into pipeline components
    """
    @staticmethod
    def create_component(name, param=None):
        """
        Converts a pipeline components name into the associated object
        :param name: The name of the component
        :param param: A dict containing any parameters for the component
        :return:
        """
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
        if name == "user_gi":
            return components.UserGIPipelineComponent(param)
        raise(RuntimeError("Given component does not exist: {}".format(name)))


logger = logging.getLogger(__name__)


def run_pipeline_wrapper(self, pipeline):
    """
    Wrapper used to run the pipeline in celery. Allows an analysis to have a celery task id.
    :param self:
    :param pipeline:
    :return:
    """
    task = run_pipeline.s(PipelineSerializer.serialize(pipeline)).apply_async()
    logger.info("Job scheduled with celery task id: {}".format(task.task_id))
    return task


@shared_task(bind=True)
def run_pipeline(self, serialized_pipeline):
    """
    Runs the pipeline
    :param self:
    :param serialized_pipeline:
    :return:
    """
    logger.info("Begin scheduling of pipeline components with celery task id: {}".format(self.request.id))
    pipeline = PipelineSerializer.deserialize(serialized_pipeline, PipelineComponentFactory())
    pipeline.analysis.celery_task_id = self.request.id
    pipeline.analysis.save()

    component_list = serialized_pipeline['pipeline_components']
    tasks = [run_pipeline_component.s(serialized_pipeline, pipeline.analysis.id, "start_pipeline")]
    for comp in component_list[1:]:
        if comp == "sigi":
            tasks.append(run_pipeline_component.s(pipeline.analysis.id, comp)
                         .on_error(sigi_error_handler.s(pipeline.analysis.id)))
        elif comp == "islandpath":
            tasks.append(run_pipeline_component.s(pipeline.analysis.id, comp)
                         .on_error(ipath_error_handler.s(pipeline.analysis.id)))
        else:
            tasks.append(run_pipeline_component.s(pipeline.analysis.id, comp))
    task_chain = chain(*tasks)

    logger.info("End scheduling of pipeline components with celery task id: {}".format(self.request.id))
    return task_chain.apply_async()


@shared_task(bind=True)
def run_pipeline_component(self, report, analysis_id, pipeline_component_name, pipeline_components_param=None):
    """
    Runs a pipeline component
    :param self:
    :param report:
    :param analysis_id:
    :param pipeline_component_name:
    :param pipeline_components_param:
    :return:
    """
    pipeline_component = PipelineComponentFactory.create_component(pipeline_component_name, pipeline_components_param)
    logger.info("Attempting to run component: {} for analysis with id: {}".format(pipeline_component_name, analysis_id))

    analysis_component = AnalysisComponent.objects.get(analysis__id=analysis_id, type__name=pipeline_component.name)
    analysis_component.celery_task_id = self.request.id
    analysis_component.save()

    return pipeline_component.run(report)


@shared_task()
def sigi_error_handler(context, exc, traceback, pipeline_id):
    logger.info("Sigi HMM Failed! analysis id: {}".format(pipeline_id))
    del context.args[0]["sigi_gis"]
    task_chain = chain(run_pipeline_component.s(context.args[0], pipeline_id, "islandpath")
                           .on_error(ipath_error_handler.s(pipeline_id)),
                       run_pipeline_component.s(pipeline_id, "merge_gis"),
                       run_pipeline_component.s(pipeline_id, "mash_mcl"),
                       run_pipeline_component.s(pipeline_id, "end_pipeline"))
    return task_chain.apply_async()

@shared_task()
def ipath_error_handler(context, exc, traceback, pipeline_id):
    logger.info("Islandpath Failed! analysis id: {}".format(pipeline_id))
    del context.args[0]["islandpath_gis"]
    if "sigi_gis" in context.args[0]:
        context.args[0]["merge_gis"] = context.args[0]["sigi_gis"]
        return chain(run_pipeline_component.s(context.args[0], pipeline_id, "mash_mcl"),
                     run_pipeline_component.s(pipeline_id, "end_pipeline")).apply_async()
    return run_pipeline_component.s(context.args[0], pipeline_id, "end_pipeline").apply_async()
