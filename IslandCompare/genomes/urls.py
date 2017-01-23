from django.conf.urls import url
from genomes import views

urlpatterns = [
    url(r'^$', views.GenomeListCreateView.as_view(), name="genome"),
    url(r'^details/(?P<pk>[0-9]+)$', views.GenomeRetrieveUpdateDestroyView.as_view(), name="genome_details"),
]
