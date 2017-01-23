from django.conf.urls import url
from analysis import views

urlpatterns = [
    url(r'^$', views.AnalysisListView.as_view(), name="analysis"),
    url(r'^details/(?P<pk>[0-9]+)$', views.AnalysisRetrieveUpdateView.as_view(), name="analysis_details"),
    url(r'^run', views.AnalysisRunView.as_view(), name="analysis_run"),
]
