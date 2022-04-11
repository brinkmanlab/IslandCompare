# Kubernetes Deployment Module

<!-- BEGIN_TF_DOCS -->
## Providers

No providers.

## Modules

| Name | Source | Version |
|------|--------|---------|
| <a name="module_galaxy"></a> [galaxy](#module\_galaxy) | ../galaxy | n/a |

## Inputs

| Name | Description | Type | Default | Required |
|------|-------------|------|---------|:--------:|
| <a name="input_data_dir"></a> [data\_dir](#input\_data\_dir) | n/a | `string` | n/a | yes |
| <a name="input_debug"></a> [debug](#input\_debug) | n/a | `bool` | `false` | no |
| <a name="input_instance"></a> [instance](#input\_instance) | n/a | `string` | `""` | no |
| <a name="input_microbedb_path"></a> [microbedb\_path](#input\_microbedb\_path) | Path to microbedb.sqlite as mounted by Galaxy | `string` | `"/cvmfs/microbedb.brinkmanlab.ca/microbedb.sqlite"` | no |
| <a name="input_namespace"></a> [namespace](#input\_namespace) | Instance of kubernetes\_namespace to provision instance resources under | `any` | n/a | yes |
| <a name="input_uwsgi_gid"></a> [uwsgi\_gid](#input\_uwsgi\_gid) | GID of Galaxy process | `number` | `null` | no |
| <a name="input_uwsgi_uid"></a> [uwsgi\_uid](#input\_uwsgi\_uid) | UID of Galaxy process | `number` | `null` | no |

## Outputs

No outputs.

## Resources

No resources.
<!-- END_TF_DOCS -->