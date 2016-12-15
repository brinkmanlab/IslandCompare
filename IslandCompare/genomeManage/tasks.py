from __future__ import absolute_import
from IslandCompare.celery import app
from celery import shared_task, chord, group
from genomeManage.models import Genome, Job, MauveAlignment, SigiHMMOutput, Parsnp
from django.conf import settings
from genomeManage.libs import mauvewrapper, sigihmmwrapper, parsnpwrapper, fileconverter, clusterer, genomeparser, vsearchwrapper, mashwrapper, giparser
from genomeManage.email import sendAnalysisCompleteEmail
from Bio import SeqIO
import os
import StringIO
from django.core.files.base import ContentFile
import datetime
import pytz
import errno
import logging
import pickle
import numpy as np
from mcl.mcl_clustering import mcl

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
    sequence.embl.save(str(sequenceid)+".embl", ContentFile(emblFileString))
    fileconverter.convertGbkToFna(settings.MEDIA_ROOT+"/"+sequence.genbank.name, faaOutputHandle)
    faaFileString = faaOutputHandle.getvalue()
    faaOutputHandle.close()
    sequence.fna.save(str(sequenceid)+".fna", ContentFile(faaFileString))
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
            currentGenome = Genome.objects.get(id=id)
            currentJob.genomes.add(currentGenome)
            # add each SIGIHMM run to the job list if user has not provided them
            if userGiPath is None and currentGenome.sigi is None:
                jobBuilder.append(runSigiHMM.s(id))

        # Run parsnp on the genomes or accept users input file
        if userNewickPath is None:
            parsnpJob = Parsnp(jobId=currentJob, isUserProvided=False)
            parsnpJob.save()
            # parsnp must complete before parallel mauve is run (sync)
            treeOutput = runParsnp(currentJob.id,sequenceIdList)
        else:
            # gets the parsnp file from the user provided path
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
        chord(group(jobBuilder))(endAnalysisPipeline.si(currentJob.id, sequenceIdList))
    except Exception:
        # Something happened, end pipeline and throw appropriate error
        endAnalysisPipeline(currentJob.id, sequenceIdList, complete=False)
        raise

@shared_task
def endAnalysisPipeline(jobId, sequenceIdList, complete=True):
    # Called at the end of analysis pipeline to set job status and send email to user
    currentJob = Job.objects.get(id=jobId)

    # Only run clustering if user did not provide a GI file
    if currentJob.optionalGIFile.name == "":
        logging.info("GI file was not found, running Mash-MCL cluster")
        mclMashCluster(jobId, 0.9, sequenceIdList)
    else:
        if not giparser.determineIfColorsProvided(currentJob.optionalGIFile.name):
            logging.info("GI file does not have colors, running Mash-MCL cluster")
            mclMashClusterGi(jobId, currentJob.optionalGIFile.name, sequenceIdList)

    try:
        parsnpjob = Parsnp.objects.get(jobId=jobId)
    except:
        parsnpjob = None
    try:
        mauvejob = MauveAlignment.objects.get(jobId=jobId)
    except:
        mauvejob = None
    if complete and (parsnpjob.success != False) and (mauvejob.success != False):
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

    mauvealignmentjob = MauveAlignment.objects.get(jobId=currentJob)

    try:
        mauvewrapper.combineMauveBackbones(backbonepaths,outputFile,orderList)
        mauvealignmentjob.backboneFile = outputFile
        mauvealignmentjob.success = True
    except:
        mauvealignmentjob.success = False

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
        logging.info("Creating Directory: "+settings.MEDIA_ROOT+"/mauve/"+str(jobId))
        os.mkdir(settings.MEDIA_ROOT+"/mauve/"+str(jobId))
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/mauve/"+str(jobId)):
            pass

    outputfilename = settings.MEDIA_ROOT+"/mauve/"+str(jobId)+"/"+str(outputbase)

    for genomeid in sequenceIdList:
        currentGenome = Genome.objects.get(id=genomeid)
        sequencePathList.append(settings.MEDIA_ROOT+"/"+currentGenome.genbank.name)

    mauvealignmentjob = MauveAlignment.objects.get(jobId=currentJob)
    mauvealignmentjob.backboneFile = outputfilename+".backbone"

    try:
        mauvewrapper.runMauve(sequencePathList,outputfilename)
    except:
        mauvealignmentjob.success=False
        raise Exception("Mauve Failed")
    finally:
        mauvealignmentjob.save()

@shared_task
def runSigiHMM(sequenceId):
    # Given a genomeIds this will run SigiHMM on the input genome file
    currentGenome = Genome.objects.get(id=sequenceId)
    outputbasename = settings.MEDIA_ROOT+"/sigi/"+currentGenome.name+sequenceId
    sigi = SigiHMMOutput(embloutput=outputbasename+".embl",gffoutput=outputbasename+".gff")

    try:
        logging.info("Running SigiHMM on " + settings.MEDIA_ROOT+"/"+currentGenome.embl.name)
        sigihmmwrapper.runSigiHMM(settings.MEDIA_ROOT+"/"+currentGenome.embl.name,
                              outputbasename+".embl",outputbasename+".gff")
        sigi.success=True
    except:
        sigi.success=False
        raise Exception("Sigi-HMM Failed")
    finally:
        sigi.save()
        currentGenome.sigi = sigi
        currentGenome.save()

@shared_task
def runParsnp(jobId, sequenceIdList, returnTree=True):
    # Given a jobId and sequenceIdList, this will create an output directory in the parsnp folder and
    # fill it with the output created by running parsnp
    # this will also update the parsnp job in the database to have the path to the tree file and the status of the job
    # returns an ordered parsnp tree on completion if returnTree=True, else returns the path to file
    currentJob = Job.objects.get(id=jobId)
    parsnpjob = Parsnp.objects.get(jobId=currentJob)
    outputDir = settings.MEDIA_ROOT+"/parsnp/"+str(jobId)
    os.mkdir(outputDir)
    fnaInputList = []
    for sequenceId in sequenceIdList:
        seq = Genome.objects.get(id=sequenceId)
        fnaInputList.append(settings.MEDIA_ROOT+"/"+seq.fna.name)
    try:
        parsnpwrapper.runParsnp(fnaInputList,outputDir)
        parsnpjob.treeFile = outputDir+"/parsnp.tree"
        parsnpjob.success = True
    except:
        parsnpjob.success = False
        raise Exception("Running Parsnp Failed")
    finally:
        parsnpjob.save()

    if returnTree:
        # get the left to right order of the outputted parsnp tree
        return parsnpwrapper.getOrderedLeavesWithGenome(outputDir+"/parsnp.tree",currentJob)

    else:
        return outputDir+"/parsnp.tree"

@shared_task
def clusterGis(jobId, sequenceIdList):
    logging.info("Running ClusterGI")
    seqRecords = []

    logging.debug("SequenceId List: " + str(sequenceIdList))
    for sequenceId in sequenceIdList:
        seq = Genome.objects.get(id=sequenceId)
        sigiFile = seq.sigi.gffoutput.name
        genbankFile = settings.MEDIA_ROOT+"/"+seq.genbank.name
        counter = 0
        for entry in sigihmmwrapper.parseSigiGFF(sigiFile):
            logging.info("Adding entry: " + str(entry) + "-" + str(counter))
            entrySequence = genomeparser.getSubsequence(genbankFile, entry['start'], entry['end'], counter)
            seqRecords.append(entrySequence)
            counter += 1

    try:
        logging.info("Creating Directory: "+settings.MEDIA_ROOT+"/vsearch/"+str(jobId))
        os.mkdir(settings.MEDIA_ROOT+"/vsearch/"+str(jobId))
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/vsearch/"+str(jobId)):
            pass

    genomeparser.writeFastaFile(settings.MEDIA_ROOT+"/vsearch/"+str(jobId)+"/completeFile.fna", seqRecords)
    vsearchwrapper.cluster(settings.MEDIA_ROOT+"/vsearch/"+str(jobId)+"/completeFile.fna", settings.MEDIA_ROOT+"/vsearch/"+str(jobId)+"/output", 0.9)
    logging.info("Ending ClusterGI")

@shared_task
def greedyMashCluster(jobId, threshold, sequenceIdList):
    logging.info("Running Greedy Mash Cluster")
    logging.debug("SequenceId List: " + str(sequenceIdList))
    fnaIslandPathsList = []
    clusteredIslandList = []

    # Create the mash job directory (Has to be done this way to avoid race condition)
    try:
        logging.info("Creating Directory: "+settings.MEDIA_ROOT+"/mash/"+str(jobId))
        os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId))
        os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna")
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)):
            pass

    for sequenceId in sequenceIdList:
        try:
            os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)):
                pass

        seq = Genome.objects.get(id=sequenceId)
        sigiFile = seq.sigi.gffoutput.name
        genbankFile = settings.MEDIA_ROOT+"/"+seq.genbank.name
        counter = 0
        for entry in sigihmmwrapper.parseSigiGFF(sigiFile):
            logging.info("Adding entry: " + str(entry) + "-" + str(counter))
            entrySequence = genomeparser.getSubsequence(genbankFile, entry['start'], entry['end'], counter)
            genomeparser.writeFastaFile(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId+"/"+str(counter), [entrySequence])
            fnaIslandPathsList.append(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId+"/"+str(counter))
            counter += 1

    while len(fnaIslandPathsList) > 0:
        if os.path.isfile(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"compoundScratch.msh"):
            os.remove(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"compoundScratch.msh")
        mashwrapper.createCompoundSketch(fnaIslandPathsList, settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/compoundScratch")

        if os.path.isfile(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output"):
            os.remove(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")
        mashwrapper.calculateMashDistance(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/compoundScratch.msh", fnaIslandPathsList[0], settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")

        distances = mashwrapper.mashOutputFileParser(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")
        distancesWithinThreshold = [island['referenceId'] for island in filter(lambda x: float(x['mashDistance'])<threshold, distances)]
        clusteredIslandList.append(distancesWithinThreshold)

        distances = filter(lambda x: float(x['mashDistance'])>=threshold, distances)
        fnaIslandPathsList = [island['referenceId'] for island in distances]

    outputList = {}
    for sequenceId in sequenceIdList:
        outputList[sequenceId] = {}

    clusterCounter = 0
    for clusteredIslands in clusteredIslandList:
        for island in clusteredIslands:
            splitIsland = island.split('/')[-2:]
            currentSequenceId = splitIsland[0]
            islandId = splitIsland[1]
            outputList[currentSequenceId][islandId] = clusterCounter
        clusterCounter += 1
    outputList['numberClusters'] = clusterCounter-1

    pickle.dump(outputList, open(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/clusters.p", "wb"))

@shared_task
def mclMashClusterGi(jobId, giFile, sequenceIdList):
    fnaIslandPathsList = []
    distanceMatrix = []

    try:
        logging.info("Creating Directory: "+settings.MEDIA_ROOT+"/mash/"+str(jobId))
        os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId))
        os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna")
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)):
            pass

    giDict = giparser.parseGiFile(giFile)

    for sequenceId in sequenceIdList:
        try:
            os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)):
                pass

        seq = Genome.objects.get(id=sequenceId)

        logging.info("Name used to parse dict: " + seq.uploadedName)
        print(giDict)

        try:
            genomeIslands = giDict[seq.uploadedName]
        except:
            genomeIslands = []
        genbankFile = settings.MEDIA_ROOT+"/"+seq.genbank.name

        counter = 0
        for island in genomeIslands:
            entrySequence = genomeparser.getSubsequence(genbankFile, island['start'], island['end'], counter)
            genomeparser.writeFastaFile(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId+"/"+str(counter), [entrySequence])
            fnaIslandPathsList.append(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId+"/"+str(counter))
            counter += 1

    mashwrapper.createCompoundSketch(fnaIslandPathsList, settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/compoundScratch")

    for island in fnaIslandPathsList:
        logging.info("Processing island: " + island)
        if os.path.isfile(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output"):
            os.remove(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")
        mashwrapper.calculateMashDistance(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/compoundScratch.msh", island, settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")

        distanceMatrix.append([float(i['mashDistance']) for i in mashwrapper.mashOutputFileParser(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")])

    baseMatrix = []
    for i in range(len(fnaIslandPathsList)):
        nextRow = []
        for j in range(len(fnaIslandPathsList)):
            nextRow.append(1)
        baseMatrix.append(nextRow)
    numpyBaseMatrix = np.array(baseMatrix)
    numpyDistanceMatrix = np.array(distanceMatrix)

    mclAdjacencyMatrix = np.subtract(numpyBaseMatrix, numpyDistanceMatrix)
    np.set_printoptions(threshold='nan')

    M, clusters = mcl(mclAdjacencyMatrix)
    outputList = {}

    for sequenceId in sequenceIdList:
        outputList[str(sequenceId)] = {}
    islandIdList = [i for i in range(len(fnaIslandPathsList))]

    numberClusters = 0

    logging.info(outputList)

    while len(islandIdList) > 0:
        numberClusters += 1
        currentCluster = clusters[islandIdList[0]]
        for i in currentCluster:
            island = fnaIslandPathsList[i]
            splitIsland = island.split('/')[-2:]
            currentSequenceId = splitIsland[0]
            islandId = splitIsland[1]
            logging.info("Current Sequence Id: " + str(currentSequenceId))
            logging.info("Current Island Id: " + str(islandId))
            outputList[str(currentSequenceId)][str(islandId)] = numberClusters - 1
        remainingIslands = filter(lambda x: x not in currentCluster, islandIdList)
        islandIdList = remainingIslands

    outputList['numberClusters'] = numberClusters-1
    pickle.dump(outputList, open(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/clusters.p", "wb"))



@shared_task
def mclMashCluster(jobId, threshold, sequenceIdList):
    fnaIslandPathsList = []
    distanceMatrix = []

    # Create the mash job directory (Has to be done this way to avoid race condition)
    try:
        logging.info("Creating Directory: "+settings.MEDIA_ROOT+"/mash/"+str(jobId))
        os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId))
        os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna")
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)):
            pass

    for sequenceId in sequenceIdList:
        try:
            os.mkdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(settings.MEDIA_ROOT+"/mash/"+str(jobId)):
                pass

        seq = Genome.objects.get(id=sequenceId)
        sigiFile = seq.sigi.gffoutput.name
        genbankFile = settings.MEDIA_ROOT+"/"+seq.genbank.name
        counter = 0
        for entry in sigihmmwrapper.parseSigiGFF(sigiFile):
            logging.info("Adding entry: " + str(entry) + "-" + str(counter))
            entrySequence = genomeparser.getSubsequence(genbankFile, entry['start'], entry['end'], counter)
            genomeparser.writeFastaFile(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId+"/"+str(counter), [entrySequence])
            fnaIslandPathsList.append(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/fna/"+sequenceId+"/"+str(counter))
            counter += 1

    mashwrapper.createCompoundSketch(fnaIslandPathsList, settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/compoundScratch")

    for island in fnaIslandPathsList:
        logging.info("Processing island: " + island)
        if os.path.isfile(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output"):
            os.remove(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")
        mashwrapper.calculateMashDistance(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/compoundScratch.msh", island, settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")

        distanceMatrix.append([float(i['mashDistance']) for i in mashwrapper.mashOutputFileParser(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/output")])

    baseMatrix = []
    for i in range(len(fnaIslandPathsList)):
        nextRow = []
        for j in range(len(fnaIslandPathsList)):
            nextRow.append(1)
        baseMatrix.append(nextRow)
    numpyBaseMatrix = np.array(baseMatrix)
    numpyDistanceMatrix = np.array(distanceMatrix)

    mclAdjacencyMatrix = np.subtract(numpyBaseMatrix, numpyDistanceMatrix)
    np.set_printoptions(threshold='nan')

    M, clusters = mcl(mclAdjacencyMatrix)
    outputList = {}

    for sequenceId in sequenceIdList:
        outputList[str(sequenceId)] = {}
    islandIdList = [i for i in range(len(fnaIslandPathsList))]

    numberClusters = 0

    logging.info(outputList)

    while len(islandIdList) > 0:
        numberClusters += 1
        currentCluster = clusters[islandIdList[0]]
        for i in currentCluster:
            island = fnaIslandPathsList[i]
            splitIsland = island.split('/')[-2:]
            currentSequenceId = splitIsland[0]
            islandId = splitIsland[1]
            logging.info("Current Sequence Id: " + str(currentSequenceId))
            logging.info("Current Island Id: " + str(islandId))
            outputList[str(currentSequenceId)][str(islandId)] = numberClusters - 1
        remainingIslands = filter(lambda x: x not in currentCluster, islandIdList)
        islandIdList = remainingIslands

    outputList['numberClusters'] = numberClusters-1
    pickle.dump(outputList, open(settings.MEDIA_ROOT+"/mash/"+str(jobId)+"/clusters.p", "wb"))
