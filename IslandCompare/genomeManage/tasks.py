from __future__ import absolute_import
from IslandCompare.celery import app
from celery import shared_task
from Bio import SeqIO
from genomeManage.models import Genome
from django.conf import settings

@shared_task
def parseGenbankFile(sequenceid):
    sequence=Genome.objects.get(id=sequenceid)
    for record in SeqIO.parse(open(settings.MEDIA_ROOT+"/"+sequence.genbank.name),"genbank"):
        sequence.name = record.id
        sequence.save()
        return None
