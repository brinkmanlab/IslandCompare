from django.db import models
from django.contrib.auth.models import User
from django.db.models.signals import post_delete
from django.dispatch import receiver


class Genome(models.Model):
    id = models.AutoField(primary_key=True)
    name = models.CharField(max_length=100)
    gbk = models.FileField(upload_to="gbk")
    owner = models.ForeignKey(User)

    class Meta:
        unique_together = ('name', 'owner')

@receiver(post_delete, sender=Genome)
def delete_gbk(sender, instance, **kwargs):
    if instance.gbk is not None:
        instance.gbk.delete(save=False)

class Gene(models.Model):
    id = models.AutoField(primary_key=True)
    type = models.CharField(max_length=4)
    gene = models.CharField(max_length=12)
    locus_tag = models.CharField(max_length=12)
    product = models.CharField(max_length=200)
    start = models.IntegerField()
    end = models.IntegerField()
    strand = models.SmallIntegerField()
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)

class GenomicIsland(models.Model):
    id = models.AutoField(primary_key=True)
    method = models.CharField(max_length=10)
    start = models.IntegerField()
    end = models.IntegerField()
    genome = models.ForeignKey(Genome, on_delete=models.CASCADE)

class UserGenomicIsland(GenomicIsland):
    analysis = models.ForeignKey("analysis.Analysis")
    color = models.CharField(max_length=30)
