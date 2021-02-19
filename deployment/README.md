# IslandCompare deployment

Example deployments are provided in this folder for various destinations. For production use, it is recommended to create your own deployment recipe
using terraform modules provided in `../desinations`. Terraform is the deployment manager software used for all deployment destinations.

To install Terraform, check that your systems package manager provides it or download it from [here](https://www.terraform.io/downloads.html).

## Run local

If you are running Docker on OSX (Mac), first see
the [related subheading below](https://github.com/brinkmanlab/IslandCompare/tree/master/deployment#osx-peculiarities). Ensure docker can
be [run without root privileges](https://docs.docker.com/engine/install/linux-postinstall/). Change the current working directory to `./docker`.
Modify `./changeme.auto.tfvars` with any custom values you like. You must at least set the `docker_gid` variable to a group id with write access
to `/var/run/docker.sock`. Run `stat /var/run/docker.sock` (or `stat -x /var/run/docker.sock` on OSX) to show the owning group id.

Run the following to start an instance on your local computer using docker:

```shell script
terraform init
./deploy.sh
```

To run an analysis, download the [IslandCompare-CLI](https://raw.githubusercontent.com/brinkmanlab/islandcompare-cli/master/islandcompare.py).
Terraform will have generated `./env.sh`, run `source ./env.sh` to configure the CLI tool. See
the [CLI ReadMe](https://github.com/brinkmanlab/islandcompare-cli/blob/master/README.md) for instructions on its use. The instructions will include
command line arguments for configuring host and api key, those can be ommitted as `./env.sh` configured them globally.

To visualize your data you can access it locally at http://localhost:8000/plugins/visualizations/islandcompare/static/index.html

Browse to http://localhost:8000/ to access the Galaxy deployment directly.

To shut down this instance, run `./destroy.sh`. This will delete the instance, all of its data, and the container images. Docker may fail to unmount
CVMFS during shutdown, run `sudo fusermount -u ./microbedb/mount` if you encounter `transport endpoint is not connected` errors.

### OSX Peculiarities

OSX does not natively support Docker, it runs Docker within a Linux virtual machine. This workaround means that support is limited to only the most
basic use case. While mounting MicrobeDB via CVMFS, it will fail with an error.

To work around this CVMFS must be installed and configured manually. Follow
the [CVMFS installation and configuration instructions](https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html#mac-os-x),
replacing `cvmfs-config.cern.ch` with `microbedb.brinkmanlab.ca` in the configuration. Refer
to [../destinations/docker/cvmfs.config](../destinations/docker/cvmfs.config) for the content of the `default.local` file in the instructions.
Use [./microbedb.brinkmanlab.ca.pub](deployment/microbedb.brinkmanlab.ca.pub) for the repository public key. You **MUST** mount the CVMFS repository
under a shared folder as configured in your Docker settings. By default `/tmp` should be included as a shared folder and you can mount the repository
to `/tmp/microbedb`.

Once CVMFS is configured and you can see the database files in the location you chose to mount the repository, add the following
to `changeme.auto.tfvars`:

```
microbedb_mount_path = "/tmp/microbedb"
enable_CVMFS = false
```

See https://github.com/docker/for-mac/issues/3431 for more information about the issue.

## Deploy to cloud

Several terraform destinations have been configured. Select one from the `./destinations/` folder that you wish to use.
Modify `./changeme.auto.tfvars` with any custom values you like. Ensure you are authenticated with your cloud provider and that the required
environment variables are set for the respective terraform provider. Review the relevant cloud provider section below for additional configuration.
Once fully prepared, run `./deploy.sh` to deploy the application to the cloud.

### AWS

Select the region deployed to by exporting `export AWS_DEFAULT_REGION='us-west-2'` or creating an aws provider configuration block in the terraform
definitions. See the [supported regions for EKS](https://docs.aws.amazon.com/general/latest/gr/eks.html) as not all regions support deployment. This
step is independent of the default region setting in the next step.

Install the [AWS CLI tool](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html)
and [aws-iam-authenticator](https://docs.aws.amazon.com/eks/latest/userguide/install-aws-iam-authenticator.html). Configure the aws cli tool by
running `aws configure` and fill in the requested info. Proceed with deployment.

Additionally:
Galaxy is deployed into an AWS EKS cluster. Run `aws-iam-authenticator token -i islandcompare --token-only` to get the required token for the
dashboard. Configure `kubectl` by running `aws eks --region us-west-2 update-kubeconfig --name islandcompare`. `--name` must match the value
of `instance` in changeme.auto.tfvars. Refer to the Kubernetes section for the remaining information.

### Azure

TODO

### Google Cloud

TODO

### OpenStack

TODO

### Kubernetes

All cloud deployments include a dashboard server that provides administrative control of the cluster. To access
it, [install kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/) and run `kubectl proxy` in a separate terminal.
Visit [here](http://localhost:8001/api/v1/namespaces/kube-system/services/https:dashboard-chart-kubernetes-dashboard:https/proxy/#/login) to access
the dashboard.

To check the state of the cluster run `kubectl describe node`.

### Existing Kubernetes cluster

Configure the Kubernetes terraform provider and deploy the `./destinations/k8s` module.

### Existing Nomad cluster

Configure the Nomad terraform provider and deploy the `./destinations/nomad` module.