# IslandCompare

Dev: [![Build Status](https://travis-ci.com/brinkmanlab/IslandCompare.svg?token=SoRFeR6YxfonSdfpVcpV&branch=dev)](https://travis-ci.com/brinkmanlab/IslandCompare)  
Master: [![Build Status](https://travis-ci.com/brinkmanlab/IslandCompare.svg?token=SoRFeR6YxfonSdfpVcpV&branch=master)](https://travis-ci.com/brinkmanlab/IslandCompare)

## Setup Development Server
**Requirements:**
- Vagrant

**Steps (DO NOT USE FOR PRODUCTION):**  
Run vagrant up  
``` bash
cd to directory with Vagrantfile
vagrant up
```
Login to server:  
``` bash
vagrant ssh
```
Run the Django Server:
``` bash
cd to directory with manage.py
python3 manage.py runserver 0:8000
```
The server will now be accessible by opening index.html in the web directory.

To run the worker (This runs stuff like Mauve):  
``` bash
cd into IslandCompare directory {the project not the app}
celery -A IslandCompare worker -l info
```

Alternatively, multiple workers can be added:
```bash
celery -A IslandCompare worker -l info --concurrency=n  (where n = number of workers)
```
## Setup Production Server
To Do: Write this!, the above steps should not be used for a production server!
- For the web client, the url to the REST api (found at web/dist/js/islandcompare-rest-urls.js) will need to be updated.
- The REST api and web client can be deployed separately.

## Software Currently Being Used
### Bioinformatic Software Used
1. Mauve
2. SIGI-HMM (Colombo)
3. Parsnp
4. Mash
5. IslandPath
### Frameworks and Libraries being Used
1. Django (https://www.djangoproject.com)
2. Django REST Framework(http://www.django-rest-framework.org)
### Visualization JavaScript Libraries Used
1. d3.phylogram.js (slightly modified)
2. newick.js 
3. d3.js

## Notes
- The IslandCompare directory contains the REST API. This can run without any client.
- The web directory contains the web page client which consumes the REST API.

## Testing
To run unit tests find the manage.py file and run the following command.
```bash
python3 manage.py test
```
To run integration tests find the manage.py file and run the following command.
```bash
python3 manage.py test --pattern="integration_tests.*"
```