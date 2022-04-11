# AWS Deployment Module
<!-- BEGIN_TF_DOCS -->
## Providers

| Name | Version |
|------|---------|
| <a name="provider_aws"></a> [aws](#provider\_aws) | n/a |

## Modules

| Name | Source | Version |
|------|--------|---------|
| <a name="module_k8s"></a> [k8s](#module\_k8s) | ../k8s | n/a |

## Inputs

| Name | Description | Type | Default | Required |
|------|-------------|------|---------|:--------:|
| <a name="input_data_dir"></a> [data\_dir](#input\_data\_dir) | n/a | `string` | n/a | yes |
| <a name="input_debug"></a> [debug](#input\_debug) | n/a | `bool` | `false` | no |
| <a name="input_eks"></a> [eks](#input\_eks) | Instance of EKS module output state | `any` | n/a | yes |
| <a name="input_instance"></a> [instance](#input\_instance) | n/a | `string` | `""` | no |
| <a name="input_microbedb_path"></a> [microbedb\_path](#input\_microbedb\_path) | Path to microbedb.sqlite as mounted by Galaxy | `string` | `"/cvmfs/microbedb.brinkmanlab.ca/microbedb.sqlite"` | no |
| <a name="input_namespace"></a> [namespace](#input\_namespace) | Instance of kubernetes\_namespace to provision instance resources under | `any` | n/a | yes |
| <a name="input_nfs_server"></a> [nfs\_server](#input\_nfs\_server) | n/a | `string` | n/a | yes |
| <a name="input_user_data_volume_name"></a> [user\_data\_volume\_name](#input\_user\_data\_volume\_name) | n/a | `string` | n/a | yes |
| <a name="input_uwsgi_gid"></a> [uwsgi\_gid](#input\_uwsgi\_gid) | GID of Galaxy process | `number` | `null` | no |
| <a name="input_uwsgi_uid"></a> [uwsgi\_uid](#input\_uwsgi\_uid) | UID of Galaxy process | `number` | `null` | no |
| <a name="input_vpc"></a> [vpc](#input\_vpc) | Instance of VPC module output state | `any` | n/a | yes |

## Outputs

No outputs.

## Resources

| Name | Type |
|------|------|
| [aws_region.current](https://registry.terraform.io/providers/hashicorp/aws/latest/docs/data-sources/region) | data source |
<!-- END_TF_DOCS -->