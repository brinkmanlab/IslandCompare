# IslandCompare deployment

Example deployments are provided in this folder for various destinations. For production use, it is recommended to create your own deployment recipe
using terraform modules provided in `../desinations`. Terraform is the deployment manager software used for all deployment destinations.

To install Terraform, check that your system's package manager provides it or download it [here](https://www.terraform.io/downloads.html).

## Run local

See [docker/](docker/) for instructions.

## Deploy to cloud

Several terraform destinations have been configured. Select one from the `./destinations/` folder that you plan to use.

### AWS

See [aws/](aws/) for instructions.

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