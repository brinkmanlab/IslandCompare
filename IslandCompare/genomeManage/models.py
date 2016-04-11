from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User
from django.db.models.signals import post_delete
from django.dispatch.dispatcher import receiver
import os
import shutil

# Create your models here.

class SigiHMMOutput(models.Model):
    id = models.AutoField(primary_key=True)
    embloutput = models.FileField(upload_to='sigi/', blank=True)
    gffoutput = models.FileField(upload_to='sigi/', blank=True)

class Genome(models.Model):
    id = models.AutoField(primary_key=True)
    givenName = models.CharField(max_length=100, null=True, unique=True)
    uploadedName = models.CharField(max_length=100)
    uploader = models.ForeignKey(User)
    length = models.BigIntegerField(blank=True, null=True)
    description = models.TextField(blank=True)
    genbank = models.FileField(upload_to='gbk/', blank=True)
    embl = models.FileField(upload_to='embl/', blank=True)
    fna = models.FileField(upload_to='fna/', blank=True)
    name = models.CharField(max_length=100)
    sigi = models.ForeignKey(SigiHMMOutput, null=True)
    class Meta:
        unique_together = ('uploader', 'uploadedName',)

class Job(models.Model):
    STATUS_CHOICES = (
        ('Q', 'On Queue'),
        ('R', 'Running'),
        ('C', 'Complete'),
        ('F', 'Failed'),
    )
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100, blank=True)
    genomes = models.ManyToManyField(Genome)
    status = models.CharField(max_length=1,choices=STATUS_CHOICES)
    jobType = models.CharField(max_length=20)
    owner = models.ForeignKey(User)
    submitTime = models.DateTimeField()
    completeTime = models.DateTimeField(null=True)
    optionalGIFile = models.FileField(upload_to='gi/', blank=True)

class MauveAlignment(models.Model):
    id = models.AutoField(primary_key=True)
    jobId = models.ForeignKey(Job)
    backboneFile = models.FileField(upload_to='mauve/',blank=True)

@receiver(post_delete, sender=MauveAlignment)
def mauveCleanUp(sender, instance, **kwargs):
    targetDirectory = os.path.dirname(instance.backboneFile.name)
    shutil.rmtree(targetDirectory)

class Parsnp(models.Model):
    id = models.AutoField(primary_key=True)
    jobId = models.ForeignKey(Job)
    treeFile = models.FileField(blank=True)

@receiver(post_delete, sender=Parsnp)
def parsnpCleanUp(sender, instance, **kwargs):
    targetDirectory = os.path.dirname(instance.treeFile.name)
    shutil.rmtree(targetDirectory)
