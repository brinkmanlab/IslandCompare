from rest_framework import generics, response
from rest_framework.permissions import IsAuthenticated
from genomes.serializers import GenomeSerializer, GenomeUploadSerializer, GeneSerializer, GenomicIslandSerializer
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
            return response.Response(data=genome_serializer.data['id'], status=201)
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

    def delete(self, request, *args, **kwargs):
        genome = self.get_object()
        for analysis in genome.analysis_set.all():
            analysis.delete()
        genome.delete()
        return response.Response(status=204)


class GenomeGeneRetrieveView(generics.RetrieveAPIView):
    """
    Retrieve Gene Information for a Genome
    """
    permission_classes = [IsAuthenticated]
    serializer_class = GeneSerializer

    def get_queryset(self):
        return Genome.objects.filter(owner=self.request.user)

    def retrieve(self, request, *args, **kwargs):
        queryset = self.get_queryset()
        url_queries = self.request.query_params

        genes = queryset.get(id=kwargs['pk']).gene_set.all()

        if 'start' in url_queries:
            genes = genes.filter(end__gte=url_queries['start'])
        if 'end' in url_queries:
            genes = genes.filter(start__lte=url_queries['end'])

        serializer = GeneSerializer(genes, many=True)

        return Response({'genes': serializer.data})

class GenomicIslandRetrieveView(generics.RetrieveAPIView):
    """
    Retrieve Genomic Islands for a Genome
    """
    permission_classes = [IsAuthenticated]
    serializer_class = GenomicIslandSerializer

    def get_queryset(self):
        return Genome.objects.filter(owner=self.request.user)

    def retrieve(self, request, *args, **kwargs):
        url_queries = self.request.query_params
        genome = Genome.objects.get(id=kwargs['pk'])
        gis = genome.genomicisland_set.all()

        if 'start' in url_queries:
            gis = gis.filter(end__gte=url_queries['start'])
        if 'end' in url_queries:
            gis = gis.filter(start__lte=url_queries['end'])

        serializer = self.serializer_class(gis, many=True)

        return Response(serializer.data)
