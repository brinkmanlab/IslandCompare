IslandCompare

master: [![Build Status](https://travis-ci.com/adrianclim/IslandCompare.svg?token=SoRFeR6YxfonSdfpVcpV&branch=master)](https://travis-ci.com/adrianclim/IslandCompare)
dev: 

**Setup Development Server**
Requirements:
    Vagrant

Steps:
    cd to directory with Vagrantfile
    vagrant up

To login to server:
    vagrant ssh

To run server:
    cd to directory with manage.py
    python manage runserver 0:8000
    The server will now be accessible at localhost:8000

To run worker threads {this runs stuff like Mauve}:
    cd into IslandCompare directory {the project not the app}
    celery -A IslandCompare worker -l info