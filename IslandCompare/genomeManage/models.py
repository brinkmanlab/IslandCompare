from __future__ import unicode_literals

from django.db import models
from django.contrib.auth.models import User

# Create your models here.

class Genome(models.Model):
    id = models.AutoField(primary_key=True)
    uploadedName = models.CharField(max_length=100)
    uploader = models.ForeignKey(User)
    genbank = models.FileField(upload_to='gbk/')
    name = models.CharField(max_length=100)
