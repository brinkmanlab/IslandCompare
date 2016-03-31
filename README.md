IslandCompare

master: [![Build Status](https://travis-ci.com/brinkmanlab/IslandCompare.svg?token=SoRFeR6YxfonSdfpVcpV&branch=master)](https://travis-ci.com/brinkmanlab/IslandCompare)
<br>
dev: [![Build Status](https://travis-ci.com/brinkmanlab/IslandCompare.svg?token=SoRFeR6YxfonSdfpVcpV&branch=dev)](https://travis-ci.com/brinkmanlab/IslandCompare)

**Setup Development Server**
Requirements:
    Vagrant

Steps:
    (THIS SHOULD NOT BE USED FOR PRODUCTION)<br>
    Would probably be ok to use locally though<br>
    cd to directory with Vagrantfile<br>
    vagrant up<br>

To login to server:<br>
    vagrant ssh<br>

To run server:<br>
    cd to directory with manage.py<br>
    python manage runserver 0:8000<br>
    The server will now be accessible at localhost:8000<br>

To run worker threads {this runs stuff like Mauve}:<br>
    cd into IslandCompare directory {the project not the app}<br>
    celery -A IslandCompare worker -l info<br>

**Setup Production Server**
To Do: Write this!, the above steps should not be used for a production server!

**Software Currently Being Used**<br>
1. Mauve<br>
2. SIGI-HMM (Colombo)<br>
3. Parsnp<br>

**Visualization JavaScript Libraries**<br>
4. d3.phylogram.js (slightly modified)<br>
5. newick.js <br>
6. d3.js <br>