# IslandCompare

Genomic island prediction software developed to facilitate the analysis of microbial 
population datasets. IslandCompare is designed to process sets of microbial genomes and present genomic island content with 
an interactive visual to enable exploration of cross-genome genomic island content.

IslandCompare exists as nothing more than a [Galaxy](http://github.com/galaxyproject/galaxy) workflow JSON file and a client side only web UI that invokes the workflow via Galaxies API.

## Installation

### Manual
- Install and configure a [Galaxy](http://github.com/galaxyproject/galaxy) instance. The minimum required version is Galaxy 20.09.
- Download the [workflow](workflow/workflows/IslandCompare.ga) and import it via Galaxies web interface.
- Publicly share the workflow via the workflow settings.
- Manually install all tools. See http://github.com/brinkmanlab/galaxy-tools for instructions.
- The visualization plugin must be installed into Galaxy. See the [multivis](http://github.com/brinkmanlab/multivis) repo for more information.

### Automated
Automated deployments were prepared using Terraform. See [./deployment/README.md](deployment/README.md).

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
