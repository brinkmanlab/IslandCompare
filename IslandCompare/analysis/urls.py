from django.conf.urls import url
from analysis import views

urlpatterns = [
    url(r'^$', views.AnalysisListView.as_view(), name="analysis"),
    url(r'^details/(?P<pk>[0-9]+)$', views.AnalysisRetrieveUpdateView.as_view(), name="analysis_details"),
    url(r'^export/(?P<pk>[0-9]+)$', views.ExportAnalysisResultView.as_view(), name="analysis_export"),
    url(r'^results/(?P<pk>[0-9]+)$', views.AnalysisResultsView.as_view(), name="analysis_results"),
    url(r'^run/', views.AnalysisRunView.as_view(), name="analysis_run"),
    url(r'^delete/(?P<pk>[0-9]+)$', views.AnalysisDestroyView.as_view(), name="analysis_destroy"),
]
