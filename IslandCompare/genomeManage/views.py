from django.shortcuts import render
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth import logout
from django.views.decorators.http import require_http_methods
from django.http import JsonResponse, HttpResponse
from models import Genome, Job, MauveAlignment, Parsnp
from django.forms.models import model_to_dict
from tasks import parseGenbankFile, runAnalysisPipeline
from django.contrib.auth.models import User
from libs import sigihmmwrapper, parsnpwrapper, gbkparser, mauvewrapper
from django.conf import settings
import logging
import datetime
import pytz
import os

HOMOLOGOUSREGIONDIFFERENCE = 500

# Create your views here.
def index(request):
    return render(request,"index.html")

def about(request):
    return render(request,"about.html")

def signIn(request):
    username = request.POST.get('username','')
    password = request.POST.get('password','')

    user = authenticate(username=username,password=password)
    if user is not None:
        login(request,user)
        return index(request)

    return render(request,"login.html")

def logOut(request):
    logout(request)
    return index(request)

@require_http_methods(["POST"])
def createUser(request):
    username = request.POST.get('username','')
    password = request.POST.get('password','')

    user = User.objects.create_user(username, username, password)
    user.save()
    return signIn(request)

@login_required(login_url='/login')
def genomeManage(request):
    # Returns the manage.html page, this page makes ajax calls to
    # getJobs(request) and getGenomes(request)
    return render(request,"manage.html")

@login_required(login_url='/login')
@require_http_methods(["POST"])
def uploadGenome(request):
    # Takes uploaded files and creates new Genome objects according to the contents of the file
    # If uploaded file not in .gbk, .gb, or .embl format than no Genome object will be created
    downloadedFiles = request.FILES.getlist('uploadedGenomes')
    for uploadedfile in downloadedFiles:
        if uploadedfile.name.endswith('.gbk') or uploadedfile.name.endswith('.gb'):
            genome = Genome(uploadedName=uploadedfile.name,uploader=request.user,genbank=uploadedfile)
            genome.save()
            parseGenbankFile(genome.id)
        elif uploadedfile.name.endswith('.embl'):
            genome = Genome(uploadedName=uploadedfile.name,uploader=request.user,embl=uploadedfile)
            genome.save()
    return index(request)

@login_required(login_url='/login')
def runComparison(request):
    # Runs Mauve, parsnp, and Sigi-HMM on the genomes given in the jobCheckList
    # jobCheckList is given as a list of Genome.id
    # Creates a Job object with status in Queue ('Q') at start
    sequencesChecked = request.POST.get("selectedSequences").split(',')
    jobName = request.POST.get("optionalJobName")

    currentJob = Job(status='Q', name=jobName, jobType='Analysis',
                     owner=request.user,submitTime=datetime.datetime.now(pytz.timezone('US/Pacific')))

    currentJob.save()
    runAnalysisPipeline.delay(currentJob.id,sequencesChecked)
    return getJobs(request)

@login_required(login_url='/login')
def retrieveMauveFile(request):
    # Given a MauveAlignment.id will return a mauve file associated to id to user
    jobid = request.GET.get('id')
    job = Job.objects.get(id=jobid)

    if job.owner != request.user:
        return HttpResponse('Unauthorized', status=401)

    mauvejob = MauveAlignment.objects.get(jobId=job)
    output = open(mauvejob.backboneFile.name,'r')
    lines = output.readlines()
    output.close()
    return HttpResponse(lines)

@login_required(login_url='/login')
def getAlignment(request):
    # Returns the alignment.html file to the user
    # Ajax calls exist within alignment.html to retrieve required information to generate visualization
    # makes calls to retrieveMauveFile(request) and retrieveGenomesInJob(request)
    jobId = request.GET.get('id')
    return render(request,"alignment.html",{'id':jobId})

# Methods below all return JSON

@login_required(login_url='/login')
def getGenomes(request):
    # Returns JSON of all genomes owned by the current user
    # Called by manage.html
    genomes = Genome.objects.filter(uploader=request.user)
    tableData = {}
    outputArray = []
    for genome in genomes:
        currentGenome = []
        currentGenome.append(genome.id)
        currentGenome.append(genome.name)
        currentGenome.append(genome.length)
        currentGenome.append(genome.description)
        currentGenome.append(genome.uploadedName)
        outputArray.append(currentGenome)
    tableData['data']=outputArray

    return JsonResponse(tableData, safe=False)

@login_required(login_url='/login')
def getJobs(request):
    # Returns JSON of all jobs owned by the current user
    # Called by manage.html
    jobs = Job.objects.filter(owner=request.user)
    tableData = dict()
    outputArray = list()
    for job in jobs:
        currentJob = []
        currentJob.append(job.id)
        currentJob.append(job.name)
        currentJob.append(job.jobType)
        currentJob.append(job.submitTime.strftime("%Y-%m-%d %H:%M:%S"))

        if job.completeTime is not None:
            currentJob.append(job.completeTime.strftime("%Y-%m-%d %H:%M:%S"))
        else:
            currentJob.append("Not Completed")

        currentJob.append(job.status)
        outputArray.append(currentJob)
    tableData['data']=outputArray

    return JsonResponse(tableData, safe=False)

@login_required(login_url='/login')
def getAlignmentJSON(request):
    # Returns all data in JSON format needed to construct an alignment on the client
    # Will only work properly when provided mauve file is in the same order as the phylogenetic tree
    jobid = request.GET.get('id','')
    job = Job.objects.get(id=jobid)

    if job.owner != request.user:
        return HttpResponse('Unauthorized', status=401)

    outputDict ={}

    # Get the phylogenetic tree in an array
    parsnpjob = Parsnp.objects.get(jobId=job)
    outputDict['tree']=parsnpwrapper.newickToArray(parsnpjob.treeFile.name)

    # Gets the leaves of the tree from left to right
    # Assume Mauve output is ordered from first genome in input file to last genome in input file
    # If this is the case than when mauve is run, input is ordered by genome id
    treeOrder = parsnpwrapper.getLeftToRightOrderTree(parsnpjob.treeFile.name)

    logging.info("Tree Order: ")
    logging.info(treeOrder)

    # Get all the genomes in a job
    genomes = job.genomes.all()
    allgenomes = []
    count = 0

    logGenomeList = []
    for genome in genomes:
        genomedata = dict()
        genomedata['id']=count
        genomedata['name']= ".".join(os.path.basename(genome.fna.name).split(".")[0:-1])
        genomedata['length'] = genome.length
        genomedata['gis'] = sigihmmwrapper.parseSigiGFF(genome.sigi.gffoutput.name)
        genomedata['genes'] = gbkparser.getGenesFromGbk(settings.MEDIA_ROOT+"/"+genome.genbank.name)
        allgenomes.append(genomedata)
        count += 1
        logGenomeList.append(genomedata['name'])
    logging.info("Genome List: ")
    logging.info(logGenomeList)

    # Order the genomes....can write a better algorithm here if needed
    OrderedGenomeList = []
    for genomename in treeOrder:
        for x in allgenomes:
            if genomename == x['name']:
                OrderedGenomeList.append(x)

    outputDict['genomes']=OrderedGenomeList

    # Only get homologous regions for sequences that are side by side on parsnp tree
    # This prepares an array containing these homologous regions
    mauvejob = MauveAlignment.objects.get(jobId=job)
    logging.info("Mauve File Being Parsed: "+mauvejob.backboneFile.name)
    outputDict['backbone'] = mauvewrapper.parseMauveBackbone(mauvejob.backboneFile.name)

    trimmedHomologousRegionsDict = {}
    for sequenceIndex in range(len(treeOrder)-1):
        topid = sequenceIndex
        bottomid = sequenceIndex+1

        logging.debug("Top Sequence: "+str(topid))
        logging.debug("Bottom Sequence: "+str(bottomid))

        sequenceRegions = []
        for region in outputDict['backbone']:
            topSequence = region[topid]
            bottomSequence = region[bottomid]

            logging.debug("Positions Top and Bottom Sequence: ")
            logging.debug((topSequence,bottomSequence))

            # Dont send regions with no homologous regions
            if not((int(topSequence[0])==0 and int(topSequence[1])==0) or (int(bottomSequence[0])==0 and int(bottomSequence[1])==0)):
                sequenceRegions.append([[int(topSequence[0]),int(topSequence[1])],[int(bottomSequence[0]),int(bottomSequence[1])]])
        # sort the sequence regions from left to right on top strand in preperation of aggregation
        sequenceRegions.sort(key=lambda x:int(x[0][0]))
        logging.debug("Current Region: \n")
        logging.debug(sequenceRegions)
        logging.debug("Size Sequence Regions: "+str(len(sequenceRegions)))

        currentRegion = 0
        logging.debug("Current Region: "+str(currentRegion))

        currentRegionValue = sequenceRegions[currentRegion]

        aggregateList = []

        # Merge homologous regions that are closer than (HOMOLOGOUSREGIONDIFFERENCE) together
        for regionIndex in range(1,len(sequenceRegions)):
            # potentially merge results if end of region 1 is close to start of region 2 (top strand)
            if abs(sequenceRegions[regionIndex][0][0] - currentRegionValue[0][1]) < HOMOLOGOUSREGIONDIFFERENCE:
                # cases to deal with, either inversion or not an inversion
                # check second strand for this condition

                # start with easiest case, no gap between regions on second strand and no inversion
                if sequenceRegions[regionIndex][1][0] >= 0 and currentRegionValue[1][1] >= 0 and abs(sequenceRegions[regionIndex][1][0] - currentRegionValue[1][1]) < HOMOLOGOUSREGIONDIFFERENCE:
                    currentRegionValue = [[currentRegionValue[0][0],sequenceRegions[regionIndex][0][1]],
                                          [currentRegionValue[1][0],sequenceRegions[regionIndex][1][1]]]

                # if both regions are inversions and have no gap between regions on second strand
                elif sequenceRegions[regionIndex][1][0]<0 and currentRegionValue[1][1]<0 and abs(sequenceRegions[regionIndex][1][1] - currentRegionValue[1][0]) < HOMOLOGOUSREGIONDIFFERENCE:
                    currentRegionValue = [[currentRegionValue[0][0],sequenceRegions[regionIndex][0][1]],
                                         [sequenceRegions[regionIndex][1][0],currentRegionValue[1][1]]]

                # if trying to match an inversion with a non-inversion then append and continue or
                # second strand gap larger than (HOMOLOGOUSREGIONDIFFERENCE)
                else:
                    aggregateList.append(currentRegionValue)
                    currentRegion = regionIndex
                    currentRegionValue = sequenceRegions[currentRegion]

            # if gap on strand 1 is too large to merge
            else:
                aggregateList.append(currentRegionValue)
                currentRegion = regionIndex
                currentRegionValue = sequenceRegions[currentRegion]

        aggregateList.append(currentRegionValue)
        trimmedHomologousRegionsDict[sequenceIndex]=aggregateList

    outputDict['backbone']=trimmedHomologousRegionsDict

    return JsonResponse(outputDict, safe=False)

@login_required(login_url='/login')
def retrieveGenomesInJob(request):
    # Returns JSON of all genomes used in a job
    # Called by alignment.html
    jobid = request.GET.get('jobid','')
    getGenomeDetails = request.GET.get('details', 1)
    job = Job.objects.get(id=jobid)
    genomes = job.genomes.all()
    data = []
    for genome in genomes:
        genomedata = model_to_dict(genome)
        del genomedata['genbank']
        del genomedata['embl']
        del genomedata['sigi']
        del genomedata['fna']
        if getGenomeDetails == 1:
            genomedata['gis'] = sigihmmwrapper.parseSigiGFF(genome.sigi.gffoutput.name)
            genomedata['genes'] = gbkparser.getGenesFromGbk(settings.MEDIA_ROOT+"/"+genome.genbank.name)
        data.append(genomedata)
    return JsonResponse(data, safe=False)

@login_required(login_url='/login')
def getTreeData(request):
    # Returns JSON of phylogeny data associated with an alignment (generated with parsnp)
    # Called by alignment.html
    jobid = request.GET.get('jobid','')
    parsnpjob = Parsnp.objects.get(jobId=jobid)
    treeFile = parsnpjob.treeFile.name
    data = parsnpwrapper.newickToArray(treeFile)
    return JsonResponse(data, safe=False)
