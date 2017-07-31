from django.db import models
from django.contrib.auth.models import User
from genomes.models import Genome
from django.utils import timezone


class Analysis(models.Model):
    id = models.AutoField(primary_key=True)
    celery_task_id = models.CharField(max_length=765)
    name = models.CharField(max_length=100)
    genomes = models.ManyToManyField(Genome)
    submit_time = models.DateTimeField(default=timezone.now)
    start_time = models.DateTimeField(null=True)
    complete_time = models.DateTimeField(null=True)
    owner = models.ForeignKey(User)
    clusters = models.TextField(blank=True)

    class Meta:
        unique_together = ('name', 'owner')


class AnalysisType(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100, unique=True)


class AnalysisComponent(models.Model):
    id = models.AutoField(primary_key=True)
    celery_task_id = models.CharField(max_length=765)
    type = models.ForeignKey(AnalysisType)
    start_time = models.DateTimeField(null=True)
    complete_time = models.DateTimeField(null=True)
    analysis = models.ForeignKey(Analysis)

    class Meta:
        unique_together = ('type', 'analysis')
