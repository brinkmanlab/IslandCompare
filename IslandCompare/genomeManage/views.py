from django.shortcuts import render
from django.contrib.auth import authenticate, login
from django.contrib.auth.decorators import login_required
from django.contrib.auth import logout

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