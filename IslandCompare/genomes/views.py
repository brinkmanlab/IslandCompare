from rest_framework import generics
from rest_framework.permissions import IsAuthenticated
from genomes.serializers import GenomeSerializer
from genomes.models import Genome

# Create your views here.


class GenomeListCreateView(generics.ListCreateAPIView):
    """
    List or Create a Genome.
    """
    permission_classes = [IsAuthenticated]
    serializer_class = GenomeSerializer

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


