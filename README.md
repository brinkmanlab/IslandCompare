#IslandCompare

master: [![Build Status](https://travis-ci.com/brinkmanlab/IslandCompare.svg?token=SoRFeR6YxfonSdfpVcpV&branch=master)](https://travis-ci.com/brinkmanlab/IslandCompare)  
dev: [![Build Status](https://travis-ci.com/brinkmanlab/IslandCompare.svg?token=SoRFeR6YxfonSdfpVcpV&branch=dev)](https://travis-ci.com/brinkmanlab/IslandCompare)

##Setup Development Server
**Requirements:**
- Vagrant

**Steps (DO NOT USE FOR PRODUCTION):**  
Run vagrant up  
``` ruby
cd to directory with Vagrantfile
vagrant up
```
Login to server:  
``` ruby
vagrant ssh
```
Run the Django Server:
``` ruby
cd to directory with manage.py
python manage.py runserver 0:8000
```
The server will now be accessible at localhost:8000  

To run the worker (This runs stuff like Mauve):  
``` ruby
cd into IslandCompare directory {the project not the app}
celery -A IslandCompare worker -l info
```

Alternatively, multiple workers can be added:
```ruby
celery -A IslandCompare worker -l info --concurrency=n  (where n = number of workers)
```
##Setup Production Server
To Do: Write this!, the above steps should not be used for a production server!

##Software Currently Being Used
###Bioinformatic Software Used
1. Mauve
2. SIGI-HMM (Colombo)
3. Parsnp

###Visualization JavaScript Libraries Used
1. d3.phylogram.js (slightly modified)
2. newick.js 
3. d3.js 


## Notes
A different database can be used (currently postgres) in production

Note: Mauve outputs a file called .bbols in the same directory as manage.py.... I'm not sure if I can change this path.
Also Parsnp outputs a file called allmums.out also in the same directory as manage.py

Also: I believe workers need to have access to the same MEDIA_ROOT folder<br>
workers will need to have access to the apps (mauve, sigihmm, parsnp)
