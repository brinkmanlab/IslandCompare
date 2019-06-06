# IslandCompare

## Installation
There are two ways this can be deployed, docker-compose or into an existing galaxy installation.
Currently the dependant tools have not been deployed on the Galaxy tool shed, additional steps will have to be taken to install the workflow.
See the Galaxy documentation on how to do that.
All tools are in https://github.com/brinkmanlab/galaxy-tools

### Docker compose
Docker recipies have been created to customise galaxy in a docker container with the workflow installed.
Install docker, docker-compose, and clone this repo. In the directory this file is in run `docker-compose up -d`.
It will launch and you can then browse to [http://localhost:80/](http://localhost:80/) to use the workflow. 
The workflow will be listed in the "Shared Workflows" section of galaxy.

### Existing Galaxy installation
Install [Ephemeris](https://ephemeris.readthedocs.io/en/latest/installation.html) on your local computer and clone this repo.
Aquire an API key from your Galaxy instance.
In the root of the cloned repo run:

```sh
workflow-to-tools -w ./workflows/IslandCompare.ga -o ./workflows/IslandCompare.deps.yml -l "IslandCompare"
shed-tools install -v -g '_url to your galaxy instance_' -a _your api key_ -t ./workflows/IslandCompare.deps.yml
workflow-install -v -g '_url to your galaxy instance_' -a _your api key_ --publish_workflows -w ./workflows/IslandCompare.ga
```

The workflow will now be installed to your Galaxy instance.

# Notes
See workflow_notes for details of the configuration settings of each tool in the workflow.
