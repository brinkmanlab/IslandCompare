from django.db import IntegrityError
from django.shortcuts import render
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.views.decorators.http import require_http_methods
from django.http import JsonResponse, HttpResponse
from models import Genome, Job, MauveAlignment, Parsnp
from django.forms.models import model_to_dict
from tasks import parseGenbankFile, runAnalysisPipeline
from django.contrib.auth.models import User
from libs import sigihmmwrapper, parsnpwrapper, gbkparser, mauvewrapper, giparser, vsearchwrapper
from django.conf import settings
from django.core.files.base import ContentFile
from django.core.files.storage import default_storage
import logging
import datetime
import pytz
import os
import pickle

# Used to determine when to merge genomic islands predicted by SIGIHMM together.
# This will merge any genomic islands closer than HOMOLOGOUSREGIONDIFFERENCE together.
HOMOLOGOUSREGIONDIFFERENCE = settings.HOMOLOGOUSREGIONDIFFERENCE

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
        if uploadedfile.name.endswith('.gbk'):
            genome = Genome(uploadedName=uploadedfile.name,uploader=request.user,genbank=uploadedfile)
            try:
                genome.save()
                parseGenbankFile(genome.id)
            except IntegrityError:
                logging.error("File with same name has already been uploaded.")
                return HttpResponse('File with same name has already been uploaded', status=403)
        elif uploadedfile.name.endswith('.gb') or uploadedfile.name.endswith('.gbff') or uploadedfile.name.endswith('.genbank'):
            uploadedfile.name = os.path.splitext(uploadedfile.name)[0] + ".gbk"
            genome = Genome(uploadedName=uploadedfile.name,uploader=request.user,genbank=uploadedfile)
            try:
                genome.save()
                parseGenbankFile(genome.id)
            except IntegrityError:
                logging.error("File with same name has already been uploaded.")
                return HttpResponse('File with same name has already been uploaded', status=403)
        elif uploadedfile.name.endswith('.embl'):
            genome = Genome(uploadedName=uploadedfile.name,uploader=request.user,embl=uploadedfile)
            try:
                genome.save()
            except IntegrityError:
                logging.error("File with same name has already been uploaded.")
                return HttpResponse('File with same name has already been uploaded', status=403)
    return index(request)

@login_required(login_url='/login')
def runComparison(request):
    # Runs Mauve, parsnp, and Sigi-HMM on the genomes given in the jobCheckList
    # jobCheckList is given as a list of Genome.id
    # Creates a Job object with status in Queue ('Q') at start
    # if optional newick path is given, then use that newick file instead of generating own
    sequencesChecked = request.POST.get("selectedSequences").split(',')
    jobName = request.POST.get("optionalJobName")
    optionalNewick = request.FILES.get('newick',None)
    optionalGi = request.FILES.get('gi',None)
    jobGenomes = Genome.objects.filter(id__in=sequencesChecked)

    currentJob = Job(status='Q', name=jobName, jobType='Analysis',
                     owner=request.user,submitTime=datetime.datetime.now(pytz.timezone('US/Pacific')))
    currentJob.save()

    # Add the genomes to the job
    currentJob.genomes = jobGenomes
    currentJob.save()

    # If an optionalGI file is given, then save the gi file
    if optionalGi is not None:
        gifile = settings.MEDIA_ROOT+"/gi/"+str(currentJob.id)+".gis"
        default_storage.save(gifile, ContentFile(optionalGi.read()))
        currentJob.optionalGIFile = gifile
        currentJob.save()
        optionalGi = gifile

    # If an optional newick file is given, then save the newick file and generate mauve alignment using it.
    if optionalNewick is not None:
        outputDir = settings.MEDIA_ROOT+"/parsnp/"+str(currentJob.id)
        os.mkdir(outputDir)

        parsnpjob = Parsnp(jobId=currentJob, isUserProvided=True)
        default_storage.save(outputDir+"/parsnp.tree", ContentFile(optionalNewick.read()))
        parsnpjob.treeFile = outputDir+"/parsnp.tree"
        parsnpjob.save()
        runAnalysisPipeline.delay(currentJob.id,sequencesChecked,parsnpjob.treeFile.name, optionalGi)
    # If an optional newick file is not given, include generation of newick file into pipeline
    else:
        runAnalysisPipeline.delay(currentJob.id,sequencesChecked, None, optionalGi)

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
        currentGenome.append(genome.givenName)
        currentGenome.append(genome.name)
        currentGenome.append(genome.length)
        currentGenome.append(genome.description)
        currentGenome.append(genome.uploadedName)
        outputArray.append(currentGenome)
    tableData['data']=outputArray

    return JsonResponse(tableData, safe=False)

@login_required(login_url='/login')
@require_http_methods(["POST"])
def updateGenome(request):
    # Updates a genomes givenName
    success = False
    try:
        genomeId = request.POST.get("id")
        newName = request.POST.get("name")

        targetGenome = Genome.objects.get(id=genomeId)
        targetGenome.givenName = newName
        targetGenome.save()
        success = True
    except:
        raise Exception("Updating Genome Name Failed")
    finally:
        return JsonResponse({"success": success})

@login_required(login_url='/login')
@require_http_methods(['POST'])
def removeJob(request):
    # Removes a job from the database
    success = False
    try:
        jobId = request.POST.get("jobId")
        job = Job.objects.get(id=jobId)
        logging.info("Deleting Job with ID: "+jobId)
        job.delete()
        success = True
    except:
        raise Exception("Deleting Job Has Failed")
    finally:
        return JsonResponse({"success": success})

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

        allGenomes = job.genomes
        logging.info("Number of genomes in job: %s", allGenomes.count())

        parsnpStatus = Parsnp.objects.filter(jobId=job.id)
        if len(parsnpStatus) > 0:
            currentJob.append(parsnpStatus[0].success)
        else:
            currentJob.append(None)

        try:
            mauveStatus = MauveAlignment.objects.get(jobId=job.id)
            currentJob.append(mauveStatus.success)
        except:
            currentJob.append(None)

        # Get all genomes in a job and get the status of all their sigi files
        numberSuccessSigi = allGenomes.filter(sigi__success=True).count()
        logging.info("Number of successful Sigi jobs: %s", numberSuccessSigi)
        if numberSuccessSigi == allGenomes.count():
            currentJob.append(True)
        elif len(allGenomes.filter(sigi__success=False)) > 0:
            currentJob.append(False)
        else:
            currentJob.append(None)

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
    getGenes = request.GET.get('getGenes','0')

    if job.owner != request.user:
        return HttpResponse('Unauthorized', status=401)

    outputDict ={}

    # Get the phylogenetic tree in an array
    parsnpjob = Parsnp.objects.get(jobId=job)
    # client no longer needs this data as it takes the newick file directly
    # outputDict['tree']=parsnpwrapper.newickToArray(parsnpjob.treeFile.name)

    # Get the raw newick file and sends it to client
    with open (parsnpjob.treeFile.name) as newickfile:
        rawNewick = newickfile.read()
        outputDict['newick'] = rawNewick

    # Get all the genomes in a job
    genomes = job.genomes.all()
    allgenomes = []
    count = 0

    # create a gi dict is optional gi file was added by user
    giDict = None
    if job.optionalGIFile != "":
        giDict = giparser.parseGiFile(job.optionalGIFile.name)

    clusterInfo = pickle.load(open("/data/mash/"+jobid+"/clusters.p", "rb"))

    def get_spaced_colors(n):
        max_value = 16581375 #255**3
        interval = int(max_value / n)
        colors = [hex(I)[2:].zfill(6) for I in range(0, max_value, interval)]

        return ['#%02x%02x%02x' % (int(i[:2], 16), int(i[2:4], 16), int(i[4:], 16)) for i in colors]

    colorIndex = get_spaced_colors(clusterInfo['numberClusters'])
    logging.info("Number colors generated: " + str(len(colorIndex)))

    for genome in genomes:
        genomedata = dict()
        genomedata['id']=count
        genomedata['givenName'] = genome.givenName
        genomedata['name']= ".".join(os.path.basename(genome.fna.name).split(".")[0:-1])
        genomedata['length'] = genome.length
        genomedata['splitName'] = "".join(genome.givenName.split(".")[0:-1])

        # if optional gi file was uploaded use those values instead of ones from sigi
        if giDict is not None:
            genomedata['gis'] = giDict[genome.uploadedName]
        else:
            genomedata['gis'] = sigihmmwrapper.parseSigiGFF(genome.sigi.gffoutput.name)

            for i in range(len(genomedata['gis'])):
                color = colorIndex[int(clusterInfo[str(genome.id)][str(i)])]
                genomedata['gis'][i]['color'] = color

        # NOTE: loading genes into the json response takes the most amount of time, so only retrieve is asked to
        if int(getGenes) == 1:
            genomedata['genes'] = gbkparser.getGenesFromGbk(settings.MEDIA_ROOT+"/"+genome.genbank.name)
        else:
            genomedata['genes'] = []
        allgenomes.append(genomedata)
        count += 1

    # if logging is set to info, then print out genome list followed by a list of all the genomes
    if logging.getLogger().isEnabledFor(logging.INFO):
        logging.info("Genome List: ")
        logging.info([".".join(os.path.basename(loggenome.fna.name).split(".")[0:-1]) for loggenome in genomes])

    # Gets the leaves of the tree from left to right
    # Assume Mauve output is ordered from first genome in input file to last genome in input file
    # If this is the case than when mauve is run, input is ordered by genome id
    treeOrder = parsnpwrapper.getLeftToRightOrderTree(parsnpjob.treeFile.name)

    logging.info("Tree Order: ")
    logging.info(treeOrder)

    # Order the genomes....can write a better algorithm here if needed
    OrderedGenomeList = []
    parsnpEntry = job.parsnp_set.all()[0]
    for genomename in treeOrder:
        print(genomename)
        if not parsnpEntry.isUserProvided:
            for x in allgenomes:
                if genomename == x['name']:
                    OrderedGenomeList.append(x)
        else:
            for x in allgenomes:
                print(x['splitName'])
                if genomename == x['splitName']:
                    OrderedGenomeList.append(x)

    outputDict['genomes']=OrderedGenomeList

    # Only get homologous regions for sequences that are side by side on parsnp tree
    # This prepares an array containing these homologous regions
    mauvejob = MauveAlignment.objects.get(jobId=job)
    logging.info("Mauve File Being Parsed: "+mauvejob.backboneFile.name)
    outputDict['backbone'] = mauvewrapper.parseMauveBackbone(mauvejob.backboneFile.name)

    # aggregates homologous regions that are closer than HOMOLOGOUSREGIONSIZE together and remove stated nonhomolous regions
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
def getGenesInJob(request):
    # returns a json response of all the genes in a job
    jobid = request.GET.get('id','')
    job = Job.objects.get(id=jobid)

    outputDict = {}
    for genome in job.genomes.all():
        outputDict[".".join(os.path.basename(genome.fna.name).split(".")[0:-1])] =\
            gbkparser.getGenesFromGbk(settings.MEDIA_ROOT+"/"+genome.genbank.name)
    return JsonResponse(outputDict)


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
