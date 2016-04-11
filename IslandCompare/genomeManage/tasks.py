from __future__ import absolute_import
from IslandCompare.celery import app
from celery import shared_task, chord, group
from genomeManage.models import Genome, Job, MauveAlignment, SigiHMMOutput, Parsnp
from django.conf import settings
from genomeManage.libs import mauvewrapper, sigihmmwrapper, parsnpwrapper, fileconverter
from genomeManage.email import sendAnalysisCompleteEmail
from Bio import SeqIO
import os
import StringIO
from django.core.files.base import ContentFile
import datetime
import pytz
import errno
import logging

@shared_task
def parseGenbankFile(sequenceid):
    # Parses a Genomes gbk file and updates the database with data from the file
    # Takes the genomes sequenceid (primary key) and returns None
    # Currently updates the name of the genome in the database
    # Also creates an embl file and faa file from the gbk file
    sequence=Genome.objects.get(id=sequenceid)
    for record in SeqIO.parse(open(settings.MEDIA_ROOT+"/"+sequence.genbank.name),"genbank"):
        sequence.givenName = record.id +"-"+str(sequenceid)
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
    sequence.fna.save(".".join(sequence.uploadedName.split(".")[0:-1])+".fna", ContentFile(faaFileString))
    sequence.save()

@shared_task
def runAnalysisPipeline(jobId,sequenceIdList,userNewickPath=None, userGiPath=None):
    # Runs mauve, sigihmm, and parsnp on the input sequence list
    # JobId is the job id, sequenceIdList is a list of genome ids
    # will update status jobId on completion of pipeline
    currentJob = Job.objects.get(id=jobId)
    currentJob.status = 'R'
    currentJob.save()

    try:
        # build the group of jobs to be run in parallel
        jobBuilder = []
        for id in sequenceIdList:
            currentJob.genomes.add(Genome.objects.get(id=id))
            # add each SIGIHMM run to the job list if user has not provided them
            if userGiPath is None:
                jobBuilder.append(runSigiHMM.s(id))

        # Run parsnp on the genomes or accept users input file
        if userNewickPath is None:
            parsnpJob = Parsnp(jobId=currentJob)
            parsnpJob.save()
            # parsnp must complete before parallel mauve is run (sync)
            treeOutput = runParsnp(currentJob.id,sequenceIdList)
        else:
            # gets the parsnp file from the user provided path
            logging.debug("File uploaded: "+userNewickPath)
            treeOutput = parsnpwrapper.getOrderedLeavesWithGenome(userNewickPath, currentJob)

        # Run mauve on the genomes, using ordered tree as a guide for which genomes to align and how
        # to merge them together, add mauve to the group of jobs to be run in parallel
        mauveJob = MauveAlignment(jobId=currentJob)
        mauveJob.save()

        # Runs parallel mauve then stitches outputs together
        jobBuilder.append(runParallelMauveAlignment.s(treeOutput,currentJob.id))
        # Runs regular mauve
        # jobBuilder.append(runMauveAlignment.s(jobId,sequenceIdList))

        # run joblist in parallel and end pipeline
        chord(group(jobBuilder))(endAnalysisPipeline.si(currentJob.id))
    except:
        # Something happened, end pipeline and throw appropriate error
        endAnalysisPipeline(currentJob.id, complete=False)
        raise Exception("Error Occured While Running Analysis Pipeline")

@shared_task
def endAnalysisPipeline(jobId, complete=True):
    # Called at the end of analysis pipeline to set job status and send email to user
    currentJob = Job.objects.get(id=jobId)
    if complete:
        currentJob.status = 'C'
        sendAnalysisCompleteEmail(currentJob.owner.email,currentJob.id)
    else:
        currentJob.status= 'F'
    currentJob.completeTime = datetime.datetime.now(pytz.timezone('US/Pacific'))
    currentJob.save()


@shared_task
def runParallelMauveAlignment(orderedIdList,jobId):
    # breaks up an orderedIdList and runs mauve down the list in pairs
    # stitches the mauve outputs into 1 file on completion
    outputPathList = []
    # run mauve alignments in parallel
    mauveJobBuilder = []
    for sequenceCounter in range(len(orderedIdList)-1):
        mauveJobBuilder.append(runMauveAlignment.s(jobId,[orderedIdList[sequenceCounter],orderedIdList[sequenceCounter+1]],sequenceCounter))
        outputPathList.append(settings.MEDIA_ROOT+"/mauve/"+str(jobId)+"/"+str(sequenceCounter)+".backbone")

    # note must give order of sequences in outputFile to merge so it can merge the sequences correctly
    # where the array position = the sequence number and its value is the column it is in
    orderMauveList = []  # the array that contains the column position of the sequences

    # iterate through a sorted id list, when sorted it gives the sequence ids in order (in mauve file)
    for x in sorted(orderedIdList):
        # add index of the found value to the ordered list (this tells us column in output file)
        orderMauveList.append(orderedIdList.index(x))

    # run all alignments and merge when all alignments are completed (Keep in same order as phylogenetic tree)
    chord(group(mauveJobBuilder))(mergeMauveAlignments.si(jobId,outputPathList,None))

@shared_task
def mergeMauveAlignments(jobId,backbonepaths,orderList):
    # merges the mauve backbone outputs together into 1 file
    currentJob = Job.objects.get(id=jobId)
    outputFile = settings.MEDIA_ROOT+"/mauve/"+str(jobId)+"/"+"merged"+".backbone"

    mauvewrapper.combineMauveBackbones(backbonepaths,outputFile,orderList)

    mauvealignmentjob = MauveAlignment.objects.get(jobId=currentJob)
    mauvealignmentjob.backboneFile = outputFile
    mauvealignmentjob.save()

@shared_task
def runMauveAlignment(jobId,sequenceIdList,outputBaseName=None):
    # Given a jobId and a list of genomeIds this will run Mauve on the input genomes gbk files
    # On start, job status will be updated to running in the database and will change on completion of function
    currentJob = Job.objects.get(id=jobId)
    sequencePathList = []
    outputbase = outputBaseName

    if outputbase is None:
        outputbase = ("-".join(sequenceIdList))

    # Create the mauve job directory (Has to be done this way to avoid race condition)
    try:
        os.mkdir(settings.MEDIA_ROOT+"/mauve/"+str(jobId))
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/mauve/"+str(jobId)):
            pass

    outputfilename = settings.MEDIA_ROOT+"/mauve/"+str(jobId)+"/"+str(outputbase)

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
    outputbasename = settings.MEDIA_ROOT+"/sigi/"+currentGenome.name+sequenceId
    sigihmmwrapper.runSigiHMM(settings.MEDIA_ROOT+"/"+currentGenome.embl.name,
                              outputbasename+".embl",outputbasename+".gff")
    sigi = SigiHMMOutput(embloutput=outputbasename+".embl",gffoutput=outputbasename+".gff")
    sigi.save()
    currentGenome.sigi = sigi
    currentGenome.save()

@shared_task
def runParsnp(jobId, sequenceIdList, returnTree=True):
    # Given a jobId and sequenceIdList, this will create an output directory in the parsnp folder and
    # fill it with the output created by running parsnp
    # this will also update the parsnp job in the database to have the path to the tree file
    # returns an ordered parsnp tree on completion if returnTree=True, else returns the path to file
    outputDir = settings.MEDIA_ROOT+"/parsnp/"+str(jobId)
    os.mkdir(outputDir)
    fnaInputList = []
    for sequenceId in sequenceIdList:
        seq = Genome.objects.get(id=sequenceId)
        fnaInputList.append(settings.MEDIA_ROOT+"/"+seq.fna.name)
    parsnpwrapper.runParsnp(fnaInputList,outputDir)
    currentJob = Job.objects.get(id=jobId)
    parsnpjob = Parsnp.objects.get(jobId=currentJob)
    parsnpjob.treeFile = outputDir+"/parsnp.tree"
    parsnpjob.save()

    if returnTree:
        # get the left to right order of the outputted parsnp tree
        return parsnpwrapper.getOrderedLeavesWithGenome(outputDir+"/parsnp.tree",currentJob)

    else:
        return outputDir+"/parsnp.tree"
