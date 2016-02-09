"""IslandCompare URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/1.9/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  url(r'^$', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  url(r'^$', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.conf.urls import url, include
    2. Add a URL to urlpatterns:  url(r'^blog/', include('blog.urls'))
"""
from django.conf.urls import url
import views

urlpatterns = [
    url(r'^$', views.index),
    url(r'^login', views.signIn),
    url(r'^logout', views.logOut),
    url(r'^manage', views.genomeManage),
    url(r'^genomeUpload', views.uploadGenome),
    url(r'^getGenomes', views.getGenomes),
    url(r'^getJobs', views.getJobs),
    url(r'^submitJob', views.runComparison),
    url(r'^retrieveMauve', views.retrieveMauveFile),
    url(r'^getAlignment',views.getAlignment)
]