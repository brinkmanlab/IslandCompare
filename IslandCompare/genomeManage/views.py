from django.shortcuts import render
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth import logout
from django.views.decorators.http import require_http_methods
from django.http import JsonResponse, HttpResponse
from models import Genome, Job, MauveAlignment
from django.forms.models import model_to_dict
from tasks import parseGenbankFile, runMauveAlignment, runSigiHMM, runParsnp
from django.contrib.auth.models import User
from libs import sigihmmwrapper

# Create your views here.
def index(request):
    return render(request,"index.html")

def signIn(request):
    username = request.POST.get('username','')
    password = request.POST.get('password','')

    user = authenticate(username=username,password=password)
    if user is not None:
        login(request,user)
        return genomeManage(request)

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
    downloadedFiles = request.FILES.getlist('genomeFiles')
    for uploadedfile in downloadedFiles:
        if uploadedfile.name.endswith('.gbk') or uploadedfile.name.endswith('.gb'):
            genome = Genome(uploadedName=uploadedfile.name,uploader=request.user,genbank=uploadedfile)
            genome.save()
            parseGenbankFile.delay(genome.id)
        elif uploadedfile.name.endswith('.embl'):
            genome = Genome(uploadedName=uploadedfile.name,uploader=request.user,embl=uploadedfile)
            genome.save()
    return index(request)

@login_required(login_url='/login')
def runComparison(request):
    # Runs Mauve and Sigi-HMM on the genomes given in the jobCheckList
    # jobCheckList is given as a list of Genome.id
    # Creates a Job object with status in Queue ('Q') at start
    sequencesChecked = request.POST.getlist('jobCheckList')
    currentJob = Job(status='Q',jobType='Mauve',owner=request.user)
    currentJob.save()
    for id in sequencesChecked:
        currentJob.genomes.add(Genome.objects.get(id=id))
        runSigiHMM.delay(id)
    mauveJob = MauveAlignment(jobId=currentJob)
    mauveJob.save()
    runMauveAlignment.delay(currentJob.id, sequencesChecked)
    runParsnp.delay(currentJob.id,sequencesChecked)
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
    data = []
    for genome in genomes:
        genomedata = model_to_dict(genome)
        del genomedata['genbank']
        del genomedata['embl']
        del genomedata['sigi']
        del genomedata['faa']
        data.append(genomedata)
    return JsonResponse(data, safe=False)

@login_required(login_url='/login')
def getJobs(request):
    # Returns JSON of all jobs owned by the current user
    # Called by manage.html
    jobs = Job.objects.filter(owner=request.user)
    data = []
    for job in jobs:
        jobdata = model_to_dict(job)
        del jobdata['genomes']
        del jobdata['owner']
        data.append(jobdata)
    return JsonResponse(data, safe=False)

@login_required(login_url='/login')
def retrieveGenomesInJob(request):
    # Returns JSON of all genomes used in a job
    # Called by alignment.html
    jobid = request.GET.get('jobid','')
    job = Job.objects.get(id=jobid)
    genomes = job.genomes.all()
    data = []
    for genome in genomes:
        genomedata = model_to_dict(genome)
        del genomedata['genbank']
        del genomedata['embl']
        del genomedata['sigi']
        del genomedata['faa']
        genomedata['gis'] = sigihmmwrapper.parseSigiGFF(genome.sigi.gffoutput.name)
        data.append(genomedata)
    return JsonResponse(data, safe=False)