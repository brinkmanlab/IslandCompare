from django.conf.urls import url
from genomes import views

urlpatterns = [
    url(r'^$', views.GenomeListView.as_view(), name="genome"),
    url(r'^upload/$', views.GenomeUploadView.as_view(), name="genome_upload"),
    url(r'^details/(?P<pk>[0-9]+)$', views.GenomeRetrieveUpdateDestroyView.as_view(), name="genome_details"),
]
