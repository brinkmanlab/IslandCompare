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