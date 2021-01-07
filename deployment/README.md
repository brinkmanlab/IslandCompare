# IslandCompare deployment

Example deployments are provided in this folder for various destinations. For production use, it is
recommended to create your own deployment recipe using the terraform modules provided in `../desinations`. Terraform
is the deployment manager software used for all deployment destinations.

To install Terraform, check that your systems package manager provides it or download it from [here](https://www.terraform.io/downloads.html).

## Run local
Ensure docker can be [run without root privileges](https://docs.docker.com/engine/install/linux-postinstall/).
Change the current working directory to `./docker`.
Modify `./changeme.auto.tfvars` with any custom values you like.

Run the following to start an instance on your local computer using docker:
```shell script
terraform init
./deploy.sh
```

Browse to http://localhost:8000/ to access the Galaxy deployment.

To run an analysis download the [IslandCompare-CLI](https://raw.githubusercontent.com/brinkmanlab/islandcompare-cli/master/islandcompare.py).
Terraform will have generated `./env.sh`, run `source ./env.sh` to configure the CLI tool. See the [CLI ReadMe](https://github.com/brinkmanlab/islandcompare-cli/blob/master/README.md) for instructions on its use.
To visualize your data you can access it locally at http://localhost:8000/plugins/visualizations/islandcompare/static/index.html

To shut down this instance, run `./destroy.sh`. This will delete the instance, all of its data, and the container images.
Docker may fail to unmount CVMFS during shutdown, run `sudo fusermount -u ./microbedb/mount` if you encounter `transport endpoint is not connected` errors.

## Deploy to cloud

Several terraform destinations have been configured. Select one from the `./destinations/` folder that you wish to use.
Modify `./changeme.auto.tfvars` with any custom values you like. Ensure you are authenticated with your cloud provider
and that the required environment variables are set for the respective terraform provider.

### AWS

Install the [AWS CLI tool](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html) and [aws-iam-authenticator](https://docs.aws.amazon.com/eks/latest/userguide/install-aws-iam-authenticator.html).
Galaxy is deployed into a AWS EKS cluster. Run `aws-iam-authenticator token -i galaxy --token-only` to get the required token for the dashboard.

Configure `kubectl` by running `aws eks --region us-west-2 update-kubeconfig --name galaxy`.

Refer to the Kubernetes section for the remaining information.

### Azure
TODO

### Google Cloud
TODO

### OpenStack
TODO

### Kubernetes

All cloud deployments include a dashboard server that provides administrative control of the cluster.
To access it, [install kubectl](https://kubernetes.io/docs/tasks/tools/install-kubectl/) and run `kubectl proxy` in a separate terminal.
Visit [here](http://localhost:8001/api/v1/namespaces/kube-system/services/https:dashboard-chart-kubernetes-dashboard:https/proxy/#/login) to
access the dashboard.

To check the state of the cluster run `kubectl describe node`.

### Existing Kubernetes cluster

Configure the Kubernetes terraform provider and deploy the `./destinations/k8s` module.

### Existing Nomad cluster

Configure the Nomad terraform provider and deploy the `./destinations/nomad` module.