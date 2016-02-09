from __future__ import absolute_import
from IslandCompare.celery import app
from celery import shared_task
from Bio import SeqIO
from genomeManage.models import Genome, Job, MauveAlignment
from django.conf import settings
from genomeManage.libs import mauvewrapper

@shared_task
def parseGenbankFile(sequenceid):
    sequence=Genome.objects.get(id=sequenceid)
    for record in SeqIO.parse(open(settings.MEDIA_ROOT+"/"+sequence.genbank.name),"genbank"):
        sequence.name = record.id
        sequence.save()
        return None

@shared_task
def runMauveAlignment(jobId,sequenceIdList):
    currentJob = Job.objects.get(id=jobId)
    currentJob.status = 'R'
    currentJob.save()
    sequencePathList = []
    outputfilename = settings.MEDIA_ROOT+"/mauve/"+("-".join(sequenceIdList))
    for genomeid in sequenceIdList:
        currentGenome = Genome.objects.get(id=genomeid)
        sequencePathList.append(settings.MEDIA_ROOT+"/"+currentGenome.genbank.name)
    mauvewrapper.runMauve(sequencePathList,outputfilename)
    mauvealignmentjob = MauveAlignment.objects.get(jobId=currentJob)
    mauvealignmentjob.backboneFile = outputfilename+".backbone"
    mauvealignmentjob.save()
    currentJob.status = 'C'
    currentJob.save()


