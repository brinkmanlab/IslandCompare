# IslandCompare

Genomic island prediction software developed to facilitate the analysis of microbial 
population datasets. IslandCompare is designed to process sets of microbial genomes and present genomic island content with 
an interactive visual to enable exploration of cross-genome genomic island content.

IslandCompare exists as nothing more than a [Galaxy](http://github.com/galaxyproject/galaxy) workflow JSON file and a 
client side only web UI that invokes the workflow via Galaxies API. A [command line interface](https://github.com/brinkmanlab/islandcompare-cli/) 
is also available that will talk to Galaxies API, invoking the workflow.

IslandCompare operates on Genbank or EMBLE formatted data. It will attempt to stitch together draft genomes as some tools do not work
with multi-contig datasets. It will also accept a premade Newick formatted file rather than generate a phylogenetic tree. The resulting output
includes a GFF3 file containing all of the results, along with the generated newick file, any stitched datasets, and a GFF3 file containing only
the genomic islands.

## Use
IslandCompare is publicly hosted for your use at https://islandcompare.ca. There you can upload data, run analysis, and visualize the result.

If you prefer to deploy your own instance of IslandCompare, a containerized deployment of Galaxy is available along with
scripts to automatically deploy the IslandCompare workflow and dependencies. See the following section for more information.
The primary intended means of interacting with a local deployment of IslandCompare is via the [command line interface](https://github.com/brinkmanlab/islandcompare-cli/).

## Installation

### Automated
Automated deployments were prepared using Terraform. See [./deployment/README.md](deployment/README.md) for more information
regarding deployment and running an analysis.

### Manual
- Install and configure a [Galaxy](http://github.com/galaxyproject/galaxy) instance. The minimum required version is Galaxy 20.09.
- Download the [workflow](workflow/workflows/IslandCompare.ga) and import it via Galaxies web interface.
- Publicly share the workflow via the workflow settings.
- Manually install all tools. See http://github.com/brinkmanlab/galaxy-tools for instructions.
- The visualization plugin must be installed into Galaxy. See the [multivis](http://github.com/brinkmanlab/multivis) repo for more information.

### Front-end
See [./ui/README.md](ui/README.md) for instructions to build the IslandCompare website.

## Project layout

* `./ui` - Web front end source code. See [./ui/README.md](ui/README.md).
* `./docs` - High level documentation of the workflow
* `./workflow/*.rule` - Rules used in the `Apply rules to collection` tool throughout the workflow
* `./workflow/prepare_worflow` - Script that automates cleaning up and copying into the repository downloaded workflows from Galaxy
* `./workflow/workflow_notes` - Text file documenting various settings of the tools throughout the workflow for reproduction purposes
* `./workflow/scripts/*.gawk` - All scripts used in the `awkscript` tool throughout the workflow
* `./workflow/workflows/*.ga` - Galaxy workflows and subworkflows

### Deployment

Terraform is used to deploy the various resources needed to run Galaxy to the cloud provider of choice.

* `./destinations` - Terraform modules responsible for deployment into the various providers
* `./deployment` - Usage examples for the destination modules

## Notes
- See workflow_notes for details of the configuration settings of each tool in the workflow.
- ./.github/workflows/nodejs.yml specifies to the GitHub CI how to automatically deploy the front-end to preconfigured urls. Changes pushed to `./ui/` will automatically be built and deployed.
