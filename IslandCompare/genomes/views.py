from rest_framework import generics, response
from rest_framework.permissions import IsAuthenticated
from genomes.serializers import GenomeSerializer, GenomeUploadSerializer, GenomeGenesSerializer
from genomes.models import Genome
from rest_framework.parsers import MultiPartParser, FormParser
from rest_framework.response import Response
import re


# Create your views here.


class GenomeListView(generics.ListAPIView):
    """
    List Genomes.
    """
    permission_classes = [IsAuthenticated]
    serializer_class = GenomeSerializer

    def get_queryset(self):
        return Genome.objects.filter(owner=self.request.user)


class GenomeUploadView(generics.CreateAPIView):
    """
    Upload a Genome.
    """
    permission_classes = [IsAuthenticated]
    serializer_class = GenomeUploadSerializer
    parser_classes = (FormParser, MultiPartParser)

    def create(self, request, *args, **kwargs):
        files = self.request.FILES['genomes']
        files.name = re.sub(r'(\.genbank|\.gbff)$', ".gbk", files.name)
        genome_serializer = GenomeSerializer(
            data={'name': files.name, 'gbk': files},
            context={'request': request}
        )
        if genome_serializer.is_valid():
            genome_serializer.save()
            return response.Response(status=201)
        else:
            errors = ""
            for key in genome_serializer.errors:
                errors += key.upper() + ": " + ", ".join(genome_serializer.errors[key]) + "\n"
            return response.Response(errors, status=400)

    def get_queryset(self):
        return Genome.objects.filter(owner=self.request.user)


class GenomeRetrieveUpdateDestroyView(generics.RetrieveUpdateDestroyAPIView):
    """
    Retrieve, Update, or Destroy a Genome.
    """
    permission_classes = [IsAuthenticated]
    serializer_class = GenomeSerializer

    def get_queryset(self):
        return Genome.objects.filter(owner=self.request.user)


class GenomeGeneRetrieveView(generics.RetrieveAPIView):
    """
    Retrieve Gene Information for a Group of Genomes
    """
    permission_classes = [IsAuthenticated]
    serializer_class = GenomeGenesSerializer

    def get_queryset(self):
        return Genome.objects.filter(owner=self.request.user)

    def retrieve(self, request, *args, **kwargs):
        queryset = self.get_queryset()
        url_queries = self.request.query_params

        serializer = GenomeGenesSerializer(queryset.get(id=kwargs['pk']))

        if 'start' in url_queries and 'end' in url_queries:
            serializer.start_cut_off = int(url_queries['start'])
            serializer.end_cut_off = int(url_queries['end'])

        return Response(serializer.data)
