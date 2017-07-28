from django.conf.urls import url
from genomes import views

urlpatterns = [
    url(r'^$', views.GenomeListView.as_view(), name="genome"),
    url(r'^genes/(?P<pk>[0-9]+)$', views.GenomeGeneRetrieveView.as_view(), name="genome_genes"),
    url(r'^upload/$', views.GenomeUploadView.as_view(), name="genome_upload"),
    url(r'^details/(?P<pk>[0-9]+)$', views.GenomeRetrieveUpdateDestroyView.as_view(), name="genome_details"),
    url(r'^islands/(?P<pk>[0-9]+)$', views.GenomicIslandRetrieveView.as_view(), name="genomic_islands")
]
