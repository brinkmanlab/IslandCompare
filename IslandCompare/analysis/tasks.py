from __future__ import absolute_import, unicode_literals
from celery import shared_task
from analysis.models import AnalysisComponent
from analysis import components
from analysis.pipeline import PipelineSerializer


class PipelineComponentFactory(object):
    @staticmethod
    def create_component(name):
        if name == "setup_gbk":
            return components.SetupGbkPipelineComponent()
        raise(RuntimeError("Given component does not exist: {}".format(name)))


def run_pipeline_wrapper(pipeline):
    return run_pipeline.s(PipelineSerializer.serialize(pipeline)).apply().get()


@shared_task(bind=True)
def run_pipeline(self, serialized_pipeline):
    pipeline = PipelineSerializer.deserialize(serialized_pipeline, PipelineComponentFactory())
    pipeline.analysis.celery_task_id = self.request.id
    pipeline.analysis.save()

    task_chain = None

    for component in pipeline.pipeline_components:
        if task_chain is None:
            task_chain = run_pipeline_component.s(serialized_pipeline, pipeline.analysis.id, component.name)
        else:
            task_chain.link(run_pipeline_component.s(pipeline.analysis.id, component.name))

    return task_chain.apply_async()


@shared_task(bind=True)
def run_pipeline_component(self, report, analysis_id, pipeline_component_name):
    pipeline_component = PipelineComponentFactory.create_component(pipeline_component_name)
    analysis_component = AnalysisComponent.objects.get(analysis__id=analysis_id,
                                                       type__name=pipeline_component.name)
    analysis_component.celery_task_id = self.request.id
    analysis_component.save()

    return pipeline_component.run(report)
