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
from celery.task.sets import TaskSet
import time
import logging

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
    sequence.fna.save(sequence.name+".fna", ContentFile(faaFileString))
    sequence.save()

@shared_task
def runAnalysisPipeline(jobId,sequenceIdList):
    # Runs mauve, sigihmm, and parsnp on the input sequence list
    # JobId is the job id, sequenceIdList is a list of genome ids
    # will update status jobId on completion of pipeline
    # TODO this has potential to deadlock fix this
    currentJob = Job.objects.get(id=jobId)
    currentJob.status = 'R'
    currentJob.save()

    try:
        # Run SIGIHMM on all of the genomes and add the genomes to the job
        sigihmmJobBuilder = []
        for id in sequenceIdList:
            currentJob.genomes.add(Genome.objects.get(id=id))
            # async running of SIGIHMM, does not matter order of when this completes
            # check near end of function for job completion
            sigihmmJobBuilder.append(runSigiHMM.s(id))
        sigihmmjobs = TaskSet(tasks=sigihmmJobBuilder)
        sigihmmresult = sigihmmjobs.apply_async()

        # Run parsnp on the genomes
        parsnpJob = Parsnp(jobId=currentJob)
        parsnpJob.save()
        # parsnp must complete before parallel mauve is run
        treeOutput = runParsnp(currentJob.id,sequenceIdList)

        # get the left to right order of the outputted parsnp tree
        newick = parsnpwrapper.newickToArray(treeOutput)
        treeOrder = parsnpwrapper.getLeftToRightOrderTree(newick)

        logging.info("TreeOrder: ")
        logging.info(treeOrder)

        genomes = currentJob.genomes.all()
        genomeDict = {}

        for genome in genomes:
            genomeDict[".".join(os.path.basename(genome.fna.name).split(".")[0:-1])] = genome.id

        treeOrderedIds = []
        for name in treeOrder:
            treeOrderedIds.append(genomeDict[name])

        # Run mauve on the genomes, using ordered tree as a guide for which genomes to align and how
        # to merge them together
        mauveJob = MauveAlignment(jobId=currentJob)
        mauveJob.save()
        runParallelMauveAlignment(currentJob.id, treeOrderedIds)

        # Check if all needed jobs are completed
        # TODO This can cause deadlock
        while not sigihmmresult.ready():
            time.sleep(30) # wait 30 seconds before checking again

        # Check if jobs were successful
        if not sigihmmresult.successful():
            raise Exception("SIGIHMM Job Failed")

        sendAnalysisCompleteEmail(currentJob.owner.email,currentJob.id)
        currentJob.status = 'C'
    except:
        currentJob.status = 'F'
        raise
    finally:
        currentJob.completeTime = datetime.datetime.now(pytz.timezone('US/Pacific'))
        currentJob.save()

@shared_task
def runParallelMauveAlignment(jobId,orderedIdList):
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

    # run all alignments and merge when all alignments are completed
    chord(group(mauveJobBuilder))(mergeMauveAlignments.si(jobId,outputPathList,orderMauveList))

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
def runParsnp(jobId, sequenceIdList):
    # Given a jobId and sequenceIdList, this will create an output directory in the parsnp folder and
    # fill it with the output created by running parsnp
    # this will also update the parsnp job in the database to have the path to the tree file
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
    return outputDir+"/parsnp.tree"
