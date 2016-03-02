from __future__ import absolute_import
from IslandCompare.celery import app
from celery import shared_task
from genomeManage.models import Genome, Job, MauveAlignment, SigiHMMOutput, Parsnp
from django.conf import settings
from genomeManage.libs import mauvewrapper, sigihmmwrapper, parsnpwrapper, fileconverter
from genomeManage.email import sendAnalysisCompleteEmail
from Bio import SeqIO
import os
import StringIO
from django.core.files.base import ContentFile

@shared_task
def parseGenbankFile(sequenceid):
    # Parses a Genomes gbk file and updates the database with data from the file
    # Takes the genomes sequenceid (primary key) and returns None
    # Currently updates the name of the genome in the database
    # Also creates an embl file and faa file from the gbk file
    sequence=Genome.objects.get(id=sequenceid)
    for record in SeqIO.parse(open(settings.MEDIA_ROOT+"/"+sequence.genbank.name),"genbank"):
        sequence.name = record.id
        sequence.length = len(record.seq)
        sequence.description = record.description
        break
    emblOutputHandle = StringIO.StringIO()
    faaOutputHandle = StringIO.StringIO()
    SeqIO.convert(settings.MEDIA_ROOT+"/"+sequence.genbank.name, "genbank",
                  emblOutputHandle, "embl")
    emblFileString = emblOutputHandle.getvalue()
    emblOutputHandle.close()
    sequence.embl.save(sequence.name+".embl", ContentFile(emblFileString))
    fileconverter.convertGbkToFna(settings.MEDIA_ROOT+"/"+sequence.genbank.name, faaOutputHandle)
    faaFileString = faaOutputHandle.getvalue()
    faaOutputHandle.close()
    sequence.faa.save(sequence.name+".fna", ContentFile(faaFileString))
    sequence.save()

@shared_task
def runAnalysisPipeline(jobId,sequenceIdList):
    # Runs mauve, sigihmm, and parsnp on the input sequence list
    # JobId is the job id, sequenceIdList is a list of genome ids
    # will update status jobId on completion of pipeline
    currentJob = Job.objects.get(id=jobId)
    currentJob.status = 'R'
    currentJob.save()

    try:
        for id in sequenceIdList:
            currentJob.genomes.add(Genome.objects.get(id=id))
            runSigiHMM(id)
        mauveJob = MauveAlignment(jobId=currentJob)
        mauveJob.save()
        runMauveAlignment(currentJob.id, sequenceIdList)
        parsnpJob = Parsnp(jobId=currentJob)
        parsnpJob.save()
        runParsnp(currentJob.id,sequenceIdList)
        sendAnalysisCompleteEmail(currentJob.owner.email,currentJob.id)
        currentJob.status = 'C'
    except:
        currentJob.status = 'F'

    currentJob.save()

@shared_task
def runMauveAlignment(jobId,sequenceIdList):
    # Given a jobId and a list of genomeIds this will run Mauve on the input genomes gbk files
    # On start, job status will be updated to running in the database and will change on completion of function
    currentJob = Job.objects.get(id=jobId)
    sequencePathList = []
    outputfilename = settings.MEDIA_ROOT+"/mauve/"+("-".join(sequenceIdList))
    for genomeid in sequenceIdList:
        currentGenome = Genome.objects.get(id=genomeid)
        sequencePathList.append(settings.MEDIA_ROOT+"/"+currentGenome.genbank.name)

    mauvewrapper.runMauve(sequencePathList,outputfilename)
    mauvealignmentjob = MauveAlignment.objects.get(jobId=currentJob)
    mauvealignmentjob.backboneFile = outputfilename+".backbone"
    mauvealignmentjob.save()

@shared_task
def runSigiHMM(sequenceId):
    # Given a genomeIds this will run SigiHMM on the input genome file
    currentGenome = Genome.objects.get(id=sequenceId)
    outputbasename = settings.MEDIA_ROOT+"/sigi/"+currentGenome.name
    sigihmmwrapper.runSigiHMM(settings.MEDIA_ROOT+"/"+currentGenome.embl.name,
                              outputbasename+".embl",outputbasename+".gff")
    sigi = SigiHMMOutput(embloutput=outputbasename+".embl",gffoutput=outputbasename+".gff")
    sigi.save()
    currentGenome.sigi = sigi
    currentGenome.save()

@shared_task
def runParsnp(jobId, sequenceIdList):
    # Given a jobId and sequenceIdList, this will create an output directory in the parsnp folder and
    # fill it with the output created by running parsnp
    # this will also update the parsnp job in the database to have the path to the tree file
    outputDir = settings.MEDIA_ROOT+"/parsnp/"+str(jobId)
    os.mkdir(outputDir)
    faaInputList = []
    for sequenceId in sequenceIdList:
        seq = Genome.objects.get(id=sequenceId)
        faaInputList.append(settings.MEDIA_ROOT+"/"+seq.faa.name)
    parsnpwrapper.runParsnp(faaInputList,outputDir)
    currentJob = Job.objects.get(id=jobId)
    parsnpjob = Parsnp.objects.get(jobId=currentJob)
    parsnpjob.treeFile = outputDir+"/parsnp.tree"
    parsnpjob.save()