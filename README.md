# IslandCompare

Genomic island prediction software developed to facilitate the analysis of microbial 
population datasets. IslandCompare is designed to process sets of microbial genomes and present genomic island content with 
an interactive visual designed to enable exploration of cross-genome genomic island content.

IslandCompare exists as nothing more than a [Galaxy](http://github.com/galaxyproject/galaxy) workflow JSON file and a client side only web UI that invokes the workflow via Galaxies API.

## Installation

Galaxy must be patched to allow registration of users via the API. See [this patch](user_create.patch).

The visualization plugin must be installed into Galaxy. See the [multivis](http://github.com/brinkmanlab/multivis) repo for more information.

### Manual
- Install and configure a [Galaxy](http://github.com/galaxyproject/galaxy) instance. The minimum required version is Galaxy 19.09.
- Download the [workflow](workflow/workflows/IslandCompare_unpacked.ga) and import it via Galaxies web interface.
- Publicly share the workflow via the workflow settings.
- Manually install all tools. See http://github.com/brinkmanlab/galaxy-tools for instructions.

### Automated
- Install and configure a [Galaxy](http://github.com/galaxyproject/galaxy) instance. The minimum required version is Galaxy 19.09.
- Git clone http://github.com/brinkmanlab/galaxy-tools into the Galaxy server and add the `galaxy-tools/tool_conf.xml` full path to the Galaxy config `tool_config_file:` list.
- Install [Ephemeris](https://ephemeris.readthedocs.io/en/latest/installation.html) on your local computer and clone this repo.
- Also clone http://github.com/brinkmanlab/galaxy-tools locally.
- Aquire an admin API key from your Galaxy instance, not the master API key.

In the root of the cloned IslandCompare repo run:

```sh
SERVER='http://yourgalaxy.com' KEY='adminapikey' WORKFLOWS=./workflow/workflows TOOLCONF=/path/to/galaxy-tools/tool_conf.xml ./deployment/post-start-actions.sh
```

The workflow will now be installed to your Galaxy instance along with all the tool dependencies. While installing tool dependencies errors may appear, wait some time and rerun the above command to ensure everything was installed.

### Front-end
See ./ui/README.md for instructions to build the IslandCompare website.

## Notes
- See workflow_notes for details of the configuration settings of each tool in the workflow.
- ./.github/workflows/nodejs.yml specifies to the GitHub CI how to automatically deploy the front-end to preconfigured urls. Changes pushed to `./ui/` will automatically be built and deployed.
