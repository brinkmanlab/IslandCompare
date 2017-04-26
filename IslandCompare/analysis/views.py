from rest_framework import generics
from rest_framework.permissions import IsAuthenticated
from rest_framework.views import APIView
from analysis.serializers import AnalysisSerializer, RunAnalysisSerializer, ReportCsvSerializer
from analysis.models import Analysis
from rest_framework.response import Response
from analysis.pipeline import Pipeline, PipelineSerializer
from analysis.tasks import run_pipeline_wrapper
from celery.result import AsyncResult

# Create your views here.


class AnalysisListView(generics.ListAPIView):
    """
    List Analyses
    """
    permission_classes = [IsAuthenticated]
    serializer_class = AnalysisSerializer

    def get_queryset(self):
        return Analysis.objects.filter(owner=self.request.user)


class AnalysisRetrieveUpdateView(generics.RetrieveUpdateAPIView):
    """
    Retrieve or Update Analysis
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

    def post(self, request):
        serializer = RunAnalysisSerializer(data=request.data, context={'request': request})
        serializer.is_valid(raise_exception=True)

        pipeline = Pipeline()
        pipeline.create_database_entry(name=serializer.validated_data['name'],
                                       genomes=serializer.validated_data['genomes'],
                                       owner=self.request.user)

        self.run_pipeline_callback(PipelineSerializer.serialize(pipeline))

        analysis_serializer = AnalysisSerializer(pipeline.analysis)
        return Response(analysis_serializer.data)


class ExportAnalysisResultView(generics.RetrieveAPIView):
    """
    Retrieve a CSV of the Analysis
    """
    permission_classes = [IsAuthenticated]
    serializer_class = AnalysisSerializer

    def get_queryset(self):
        return Analysis.objects.filter(owner=self.request.user)

    def retrieve(self, request, *args, **kwargs):
        analysis = Analysis.objects.get(id=kwargs['pk'])
        task = AsyncResult(analysis.celery_task_id)

        if task.status == 'SUCCESS':
            result = task.get()
            return Response(ReportCsvSerializer(result).data)
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
