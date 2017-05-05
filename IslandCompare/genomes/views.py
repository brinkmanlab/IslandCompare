from rest_framework import generics, response
from rest_framework.permissions import IsAuthenticated
from genomes.serializers import GenomeSerializer, GenomeUploadSerializer
from genomes.models import Genome
from rest_framework.parsers import MultiPartParser, FormParser


# Create your views here.


class GenomeListView(generics.ListAPIView):
    """
    List a Genome.
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
        genome_serializer = GenomeSerializer(
            data={'name': files.name, 'gbk': files},
            context={'request': request}
        )
        if genome_serializer.is_valid():
            genome_serializer.save()
            return response.Response(status=201)
        else:
            return response.Response(status=400)

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


