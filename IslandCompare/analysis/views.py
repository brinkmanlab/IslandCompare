import re
from rest_framework import generics
from rest_framework.permissions import IsAuthenticated
from rest_framework.views import APIView
from analysis.serializers import AnalysisSerializer, RunAnalysisSerializer, ReportCsvSerializer, \
    ReportVisualizationOverviewSerializer, ReportGeneCsvSerializer, AnalysisGenomicIslandSerializer
from analysis.models import Analysis, AnalysisComponent
from rest_framework.response import Response
from analysis.pipeline import Pipeline
from analysis import components
from analysis.tasks import run_pipeline_wrapper
from celery.result import AsyncResult
from rest_framework.parsers import FormParser, MultiPartParser
from genomes.models import GenomicIsland



class AnalysisListView(generics.ListAPIView):
    """
    List Analyses
    """
    permission_classes = [IsAuthenticated]
    serializer_class = AnalysisSerializer

    def get_queryset(self):
        queryset = Analysis.objects.filter(owner=self.request.user).order_by('-id')

        num_list = self.request.query_params.get('num', None)
        if num_list is not None:
            queryset = queryset[:int(num_list)]

        return queryset


class AnalysisRetrieveUpdateView(generics.RetrieveUpdateAPIView):
    """
    Retrieve or Update Analysis
    """
    permission_classes = [IsAuthenticated]
    serializer_class = AnalysisSerializer

    def get_queryset(self):
        return Analysis.objects.filter(owner=self.request.user)

class AnalysisDestroyView(generics.RetrieveDestroyAPIView):
    """
    Destroy Analysis
    """
    permission_classes = [IsAuthenticated]
    serializer_class = AnalysisSerializer

    def get_queryset(self):
        return Analysis.objects.filter(owner=self.request.user)

class AnalysisRunView(APIView):
    """
    Run an Analysis
    """
    permission_classes = [IsAuthenticated]
    run_pipeline_callback = run_pipeline_wrapper
    parser_classes = (MultiPartParser, FormParser,)

    def post(self, request):
        serializer = RunAnalysisSerializer(data=request.data, context={'request': request})
        serializer.is_valid(raise_exception=True)

        pipeline = Pipeline()
        pipeline.append_component(components.StartPipelineComponent())
        pipeline.append_component(components.SetupGbkPipelineComponent())
        pipeline.append_component(components.GbkMetadataComponent())
        pipeline.append_component(components.RGIPipelineComponent())

        if 'newick' in serializer.validated_data:
            serializer.validated_data['newick'].seek(0)
            user_newick_component = components.UserNewickPipelineComponent()
            user_newick_file_contents = serializer.validated_data['newick'].read().decode('utf-8')
            user_newick_file_contents = re.sub(r'(\.genbank|\.gbff)', ".gbk", user_newick_file_contents)
            user_newick_component.set_newick(user_newick_file_contents)
            pipeline.append_component(user_newick_component)
        else:
            pipeline.append_component(components.ParsnpPipelineComponent())

        pipeline.append_component(components.MauvePipelineComponent())

        if 'gi' in serializer.validated_data:
            serializer.validated_data['gi'].seek(0)
            user_gi_component = components.UserGIPipelineComponent()
            user_gi_file_contents = serializer.validated_data['gi'].read().decode('utf-8')
            user_gi_file_contents = re.sub(r'(\.genbank|\.gbff)', ".gbk", user_gi_file_contents)
            user_gi_component.set_gi(user_gi_file_contents)
            pipeline.append_component(user_gi_component)
        else:
            pipeline.append_component(components.SigiHMMPipelineComponent())
            pipeline.append_component(components.IslandPathPipelineComponent())
            pipeline.append_component(components.MergeIslandsPipelineComponent())
            pipeline.append_component(components.MashMclClusterPipelineComponent())

        pipeline.append_component(components.EndPipelineComponent())
        pipeline.create_database_entry(name=serializer.validated_data['name'],
                                       genomes=serializer.validated_data['genomes'],
                                       owner=self.request.user)

        self.run_pipeline_callback(pipeline)

        # Analysis in Pipeline object can now be stale so retrieve the data again
        analysis = Analysis.objects.get(id=pipeline.analysis.id)

        analysis_serializer = AnalysisSerializer(analysis)
        return Response(analysis_serializer.data)


class AnalysisResultsView(generics.RetrieveAPIView):
    """
    Retrieve the Visualization Overview JSON
    """
    permission_classes = [IsAuthenticated]
    serializer_class = AnalysisSerializer

    def get_queryset(self):
        return Analysis.objects.filter(owner=self.request.user)

    def retrieve(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=kwargs['pk'])
        task = AsyncResult(analysis.celery_task_id)

        if task.status == 'SUCCESS':
            end_pipeline = AnalysisComponent.objects.get(analysis=analysis,
                                                         type__name="end_pipeline")
            result = AsyncResult(end_pipeline.celery_task_id).get()
            return Response(ReportVisualizationOverviewSerializer(result).data)
        elif task.status == 'FAILURE':
            response = Response()
            response.status_code = 500
            response.data = {'content': "Job has errored"}
            return response
        else:
            response = Response()
            response.status_code = 202
            response.data = {'content': "Job has not completed"}
            return response


class ExportAnalysisResultView(generics.RetrieveAPIView):
    """
    Retrieve a CSV of the Analysis GIs
    """
    permission_classes = [IsAuthenticated]
    serializer_class = AnalysisSerializer

    def get_queryset(self):
        return Analysis.objects.filter(owner=self.request.user)

    def retrieve(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=kwargs['pk'])
        task = AsyncResult(analysis.celery_task_id)

        if task.status == 'SUCCESS':
            return Response(ReportCsvSerializer(analysis).data)
        elif task.status == 'FAILURE':
            response = Response()
            response.status_code = 500
            response.data = {'content': "Job has errored"}
            return response
        else:
            response = Response()
            response.status_code = 202
            response.data = {'content': "Job has not completed"}
            return response

class ExportAnalysisGenesView(generics.RetrieveAPIView):
    """
    Retrieve a CSV of the Analysis GI genes
    """
    permission_classes = [IsAuthenticated]

    def retrieve(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=kwargs['pk'])
        task = AsyncResult(analysis.celery_task_id)

        if task.status == 'SUCCESS':
            return Response(ReportGeneCsvSerializer(analysis).data)
        elif task.status == 'FAILURE':
            response = Response()
            response.status_code = 500
            response.data = {'content': "Job has errored"}
            return response
        else:
            response = Response()
            response.status_code = 202
            response.data = {'content': "Job has not completed"}
            return response

class AnalysisGenomicIslandRetrieveView(generics.RetrieveAPIView):
    """
    Retrieve GIs of the analysis
    """
    permission_classes = [IsAuthenticated]
    serializer_class = AnalysisGenomicIslandSerializer

    def get_queryset(self):
        return Analysis.objects.filter(owner=self.request.user)

    def retrieve(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=kwargs['pk'])
        gis = GenomicIsland.objects.filter(genome__in=analysis.genomes.all())

        serializer = self.serializer_class(gis, many=True)

        return Response(serializer.data)
