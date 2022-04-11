# IslandCompare AWS example deployment

Change the current working directory to `./deployment/aws`. Modify `./changeme.auto.tfvars` with any custom values you like, especially the region.
See the [supported regions for EKS](https://docs.aws.amazon.com/general/latest/gr/eks.html) as not all regions support deployment. This step is independent of the default region setting in the next step.

Install the [AWS CLI tool](https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html) and [aws-iam-authenticator](https://docs.aws.amazon.com/eks/latest/userguide/install-aws-iam-authenticator.html).
Configure the aws cli tool by running `aws configure` and fill in the requested info.

Once fully prepared, run `./deploy.sh` to deploy the application to the cloud.

Additionally:
Galaxy is deployed into an AWS EKS cluster. Run `aws-iam-authenticator token -i islandcompare --token-only` to get the required token for the
dashboard. Configure `kubectl` by running `aws eks --region us-west-2 update-kubeconfig --name islandcompare`. `--name` must match the value
of `instance` in changeme.auto.tfvars. Refer to the Kubernetes section for the remaining information.


<!-- BEGIN_TF_DOCS -->
## Providers

| Name | Version |
|------|---------|
| <a name="provider_aws"></a> [aws](#provider\_aws) | n/a |
| <a name="provider_kubernetes"></a> [kubernetes](#provider\_kubernetes) | n/a |

## Modules

| Name | Source | Version |
|------|--------|---------|
| <a name="module_admin_user"></a> [admin\_user](#module\_admin\_user) | github.com/brinkmanlab/galaxy-container.git//modules/bootstrap_admin | n/a |
| <a name="module_cloud"></a> [cloud](#module\_cloud) | github.com/brinkmanlab/cloud_recipes.git//aws | n/a |
| <a name="module_cvmfs"></a> [cvmfs](#module\_cvmfs) | github.com/brinkmanlab/cloud_recipes.git//util/k8s/cvmfs | n/a |
| <a name="module_galaxy"></a> [galaxy](#module\_galaxy) | github.com/brinkmanlab/galaxy-container.git//destinations/aws | n/a |
| <a name="module_islandcompare"></a> [islandcompare](#module\_islandcompare) | ../../destinations/aws | n/a |

## Inputs

| Name | Description | Type | Default | Required |
|------|-------------|------|---------|:--------:|
| <a name="input_debug"></a> [debug](#input\_debug) | Enabling will put the deployment into a mode suitable for debugging | `bool` | `false` | no |
| <a name="input_email"></a> [email](#input\_email) | Email address to send automated emails from | `string` | n/a | yes |
| <a name="input_instance"></a> [instance](#input\_instance) | Unique deployment instance identifier | `string` | n/a | yes |
| <a name="input_region"></a> [region](#input\_region) | AWS region to deploy into | `string` | n/a | yes |

## Outputs

| Name | Description |
|------|-------------|
| <a name="output_galaxy_admin_password"></a> [galaxy\_admin\_password](#output\_galaxy\_admin\_password) | n/a |
| <a name="output_galaxy_endpoint"></a> [galaxy\_endpoint](#output\_galaxy\_endpoint) | n/a |

## Resources

| Name | Type |
|------|------|
| [kubernetes_namespace.instance](https://registry.terraform.io/providers/hashicorp/kubernetes/latest/docs/resources/namespace) | resource |
| [kubernetes_persistent_volume_claim.microbedb](https://registry.terraform.io/providers/hashicorp/kubernetes/latest/docs/resources/persistent_volume_claim) | resource |
| [aws_eks_cluster.cluster](https://registry.terraform.io/providers/hashicorp/aws/latest/docs/data-sources/eks_cluster) | data source |
| [aws_eks_cluster_auth.cluster](https://registry.terraform.io/providers/hashicorp/aws/latest/docs/data-sources/eks_cluster_auth) | data source |
<!-- END_TF_DOCS -->