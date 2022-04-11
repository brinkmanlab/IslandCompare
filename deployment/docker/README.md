# IslandCompare Docker example deployment

If you are running Docker on OSX (Mac), first see the [related subheading below](#osx-peculiarities). Ensure that you are running the most recent version of docker and that docker can
be [run without root privileges](https://docs.docker.com/engine/install/linux-postinstall/). Change the current working directory to `./docker`.
Modify `./changeme.auto.tfvars` with your custom values. You must at least set the `docker_gid` variable to a group id with write access
to `/var/run/docker.sock`. Run `stat /var/run/docker.sock` to show the owning group id.

Run the following to start an instance on your local computer using docker:

```shell script
terraform init
./deploy.sh
```

To run an analysis, download the [IslandCompare-CLI](https://raw.githubusercontent.com/brinkmanlab/islandcompare-cli/master/islandcompare.py) tool.
Terraform will have generated `./env.sh`; run `source ./env.sh` to configure the CLI tool. See
the [CLI ReadMe](https://github.com/brinkmanlab/islandcompare-cli/blob/master/README.md) for instructions on its use. The instructions will include
command line arguments for configuring host and api key, which can be ommitted as `./env.sh` configured them globally.

To visualize your data you can access it locally at http://localhost:8000/plugins/visualizations/islandcompare/static/index.html

Browse to http://localhost:8000/ to access the Galaxy deployment directly. When logging in, use the username "admin" and locate your admin password by running `terraform output -json` and identifying the value under "galaxy_admin_password".

To shut down this instance, run `./destroy.sh`. This will delete the instance, all of its data, and the container images. Docker may fail to unmount
CVMFS during shutdown; run `sudo fusermount -u ./microbedb/mount` if you encounter `transport endpoint is not connected` errors.

### OSX Peculiarities

OSX does not natively support Docker, it runs Docker within a Linux virtual machine. Due to this, mounting MicrobeDB via CVMFS following the instructions above will fail with an error.

To work around this CVMFS must be installed and configured manually. First ensure that [FUSE](http://osxfuse.github.io/) is enabled by
running `kextstat | grep -i fuse`. Download the [CVMFS package](https://ecsft.cern.ch/dist/cvmfs/cvmfs-2.8.0/cvmfs-2.8.0.pkg). Install the pkg and
reboot. Copy [../destinations/docker/cvmfs.config](../destinations/docker/cvmfs.config) to `/etc/cvmfs/default.local`.
Copy [./microbedb.brinkmanlab.ca.pub](./microbedb.brinkmanlab.ca.pub) to `/etc/cvmfs/keys/microbedb.brinkmanlab.ca.pub`. Ensure everything is
configured properly by running `sudo cvmfs_config chksetup`. You **MUST** mount the CVMFS repository under a shared folder as configured in your
Docker settings. By default, `/tmp` should be included as a shared folder, and you can mount the repository to `/tmp/microbedb`. Ensure `/tmp/microbedb`
exists and run `sudo mount -t cvmfs microbedb.brinkmanlab.ca /tmp/microbedb`.

Once CVMFS is configured and you can see the database files in the location you chose to mount the repository, add the following
to `changeme.auto.tfvars`:

```hcl
docker_gid = 0
microbedb_mount_path = "/tmp/microbedb"
docker_socket_path = "/run/host-services/docker.proxy.sock"
enable_CVMFS = false
```

The last modification you need to make is to allow group write access to the docker socket within containers. To do this, run the following:
```shell
docker run --rm -v /run/host-services/docker.proxy.sock:/run/host-services/docker.proxy.sock alpine chmod g+w /run/host-services/docker.proxy.sock
```

Mounting CVMFS and the above command need to be run any time you restart your system before you can run deploy.sh or submit analyses.
See https://github.com/docker/for-mac/issues/3431 for more information about the issue.

Appropriate resources [must be allocated](https://stackoverflow.com/a/50770267/15446750) to the Docker VM or OSX will randomly kill the application and tools during an analysis.
A minimum of 8GB of RAM and 16GB of swap space is required. It is recommended to provide as much swap space as possible to avoid out of memory issues.


<!-- BEGIN_TF_DOCS -->
## Providers

| Name | Version |
|------|---------|
| <a name="provider_local"></a> [local](#provider\_local) | 2.1.0 |

## Modules

| Name | Source | Version |
|------|--------|---------|
| <a name="module_admin_user"></a> [admin\_user](#module\_admin\_user) | ../../../galaxy-container/modules/bootstrap_admin | n/a |
| <a name="module_galaxy"></a> [galaxy](#module\_galaxy) | github.com/brinkmanlab/galaxy-container.git//destinations/docker | n/a |
| <a name="module_islandcompare"></a> [islandcompare](#module\_islandcompare) | ../../destinations/docker | n/a |

## Inputs

| Name | Description | Type | Default | Required |
|------|-------------|------|---------|:--------:|
| <a name="input_debug"></a> [debug](#input\_debug) | Enabling will put the deployment into a mode suitable for debugging | `bool` | n/a | yes |
| <a name="input_docker_gid"></a> [docker\_gid](#input\_docker\_gid) | GID with write permission to /var/run/docker.sock | `number` | n/a | yes |
| <a name="input_docker_socket_path"></a> [docker\_socket\_path](#input\_docker\_socket\_path) | Host path to docker socket | `string` | `"/var/run/docker.sock"` | no |
| <a name="input_email"></a> [email](#input\_email) | Email address to send automated emails from | `string` | n/a | yes |
| <a name="input_enable_CVMFS"></a> [enable\_CVMFS](#input\_enable\_CVMFS) | Automatically mount CVMFS from a preconfigured Docker container (Not available for OSX and Windows) | `bool` | `true` | no |
| <a name="input_host_port"></a> [host\_port](#input\_host\_port) | Host port to expose galaxy service | `number` | n/a | yes |
| <a name="input_instance"></a> [instance](#input\_instance) | Unique deployment instance identifier | `string` | n/a | yes |
| <a name="input_microbedb_mount_path"></a> [microbedb\_mount\_path](#input\_microbedb\_mount\_path) | Path on host to mount microbedb to share with jobs | `string` | `"./microbedb/mount"` | no |

## Outputs

| Name | Description |
|------|-------------|
| <a name="output_galaxy_admin_api_key"></a> [galaxy\_admin\_api\_key](#output\_galaxy\_admin\_api\_key) | n/a |
| <a name="output_galaxy_admin_password"></a> [galaxy\_admin\_password](#output\_galaxy\_admin\_password) | n/a |
| <a name="output_galaxy_endpoint"></a> [galaxy\_endpoint](#output\_galaxy\_endpoint) | n/a |
| <a name="output_galaxy_master_api_key"></a> [galaxy\_master\_api\_key](#output\_galaxy\_master\_api\_key) | n/a |

## Resources

| Name | Type |
|------|------|
| [local_file.env](https://registry.terraform.io/providers/hashicorp/local/latest/docs/resources/file) | resource |
<!-- END_TF_DOCS -->