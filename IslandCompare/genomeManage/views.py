from django.shortcuts import render
from django.contrib.auth import authenticate, login

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

def genomeManage(request):
    return render(request,"manage.html")