from django.shortcuts import render
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth import logout
from django.views.decorators.http import require_http_methods
from models import Genome

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
    for genbankFile in downloadedFiles:
        genome = Genome(uploadedName=genbankFile.name,uploader=request.user,genbank=genbankFile)
        genome.save()
    #TODO Make this redirect to genome Manage later, use ajax
    return index(request)

