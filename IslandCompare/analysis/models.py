from django.db import models
from django.contrib.auth.models import User
from genomes.models import Genome
from django.utils import timezone

# Create your models here.


class Analysis(models.Model):
    id = models.AutoField(primary_key=True)
    celery_task_id = models.CharField(max_length=765, unique=True)
    name = models.CharField(max_length=100)
    genomes = models.ManyToManyField(Genome)
    submit_time = models.DateTimeField(default=timezone.now)
    start_time = models.DateTimeField(null=True)
    complete_time = models.DateTimeField(null=True)
    owner = models.ForeignKey(User)

    class Meta:
        unique_together = ('name', 'owner')


class AnalysisType(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100, unique=True)


class AnalysisComponent(models.Model):
    id = models.AutoField(primary_key=True)
    celery_task_id = models.CharField(max_length=765, unique=True)
    type = models.ForeignKey(AnalysisType)
    start_time = models.DateTimeField(null=True)
    complete_time = models.DateTimeField(null=True)
    analysis = models.ForeignKey(Analysis)

    class Meta:
        unique_together = ('type', 'analysis')
