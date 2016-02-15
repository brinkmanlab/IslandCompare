from __future__ import absolute_import
from IslandCompare.celery import app
from celery import shared_task
from Bio import SeqIO
from genomeManage.models import Genome, Job, MauveAlignment
from django.conf import settings
from genomeManage.libs import mauvewrapper
from Bio import SeqIO

@shared_task
def parseGenbankFile(sequenceid):
    # Parses a Genomes gbk file and updates the database with data from the file
    # Takes the genomes sequenceid (primary key) and returns None
    # Currently updates the name of the genome in the database
    # Also creates an embl file from the gbk file
    sequence=Genome.objects.get(id=sequenceid)
    for record in SeqIO.parse(open(settings.MEDIA_ROOT+"/"+sequence.genbank.name),"genbank"):
        sequence.name = record.id
        break
    SeqIO.convert(settings.MEDIA_ROOT+"/"+sequence.genbank.name, "genbank",
                  settings.MEDIA_ROOT+"/embl/"+sequence.name+".embl", "embl")
    sequence.embl = settings.MEDIA_ROOT+"/embl/"+sequence.name+".embl"
    sequence.save()

@shared_task
def runMauveAlignment(jobId,sequenceIdList):
    # Given a jobId and a list of genomeIds this will run Mauve on the input genomes gbk files
    # On start, job status will be updated to running in the database and will change on completion of function
    currentJob = Job.objects.get(id=jobId)
    currentJob.status = 'R'
    currentJob.save()
    sequencePathList = []
    outputfilename = settings.MEDIA_ROOT+"/mauve/"+("-".join(sequenceIdList))
    for genomeid in sequenceIdList:
        currentGenome = Genome.objects.get(id=genomeid)
        sequencePathList.append(settings.MEDIA_ROOT+"/"+currentGenome.genbank.name)

    try:
        mauvewrapper.runMauve(sequencePathList,outputfilename)
        mauvealignmentjob = MauveAlignment.objects.get(jobId=currentJob)
        mauvealignmentjob.backboneFile = outputfilename+".backbone"
        mauvealignmentjob.save()
        currentJob.status = 'C'
    except:
        currentJob.status = 'F'
    currentJob.save()


