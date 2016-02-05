from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User

# Create your models here.

class Genome(models.Model):
    id = models.AutoField(primary_key=True)
    uploadedName = models.CharField(max_length=100)
    uploader = models.ForeignKey(User)
    genbank = models.FileField(upload_to='gbk/', blank=True)
    embl = models.FileField(upload_to='embl/', blank=True)
    name = models.CharField(max_length=100)

class Job(models.Model):
    STATUS_CHOICES = (
        ('Q', 'On Queue'),
        ('R', 'Running'),
        ('C', 'Complete'),
        ('F', 'Failed'),
    )
    id = models.AutoField(primary_key=True)
    genomes = models.ManyToManyField(Genome)
    status = models.CharField(max_length=1,choices=STATUS_CHOICES)
    jobType = models.CharField(max_length=20)
