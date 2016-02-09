from django.shortcuts import render
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth import logout
from django.views.decorators.http import require_http_methods
from django.http import JsonResponse
from models import Genome
from django.forms.models import model_to_dict
from tasks import parseGenbankFile

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

@login_required(login_url='/login')
def genomeManage(request):
    return render(request,"manage.html")

@login_required(login_url='/login')
@require_http_methods(["POST"])
def uploadGenome(request):
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

# Methods below all return JSON

@login_required(login_url='/login')
def getGenomes(request):
    genomes = Genome.objects.filter(uploader=request.user)
    data = []
    for genome in genomes:
        genomedata = model_to_dict(genome)
        del genomedata['genbank']
        del genomedata['embl']
        data.append(genomedata)
    return JsonResponse(data, safe=False)
