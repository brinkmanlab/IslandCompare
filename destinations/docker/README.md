# Docker Deployment Module
<!-- BEGIN_TF_DOCS -->
## Providers

| Name | Version |
|------|---------|
| <a name="provider_docker"></a> [docker](#provider\_docker) | n/a |
| <a name="provider_null"></a> [null](#provider\_null) | n/a |

## Modules

| Name | Source | Version |
|------|--------|---------|
| <a name="module_galaxy"></a> [galaxy](#module\_galaxy) | ../galaxy | n/a |

## Inputs

| Name | Description | Type | Default | Required |
|------|-------------|------|---------|:--------:|
| <a name="input_data_dir"></a> [data\_dir](#input\_data\_dir) | n/a | `string` | n/a | yes |
| <a name="input_debug"></a> [debug](#input\_debug) | n/a | `bool` | `false` | no |
| <a name="input_enable_CVMFS"></a> [enable\_CVMFS](#input\_enable\_CVMFS) | Automatically mount CVMFS from a preconfigured Docker container (Not available for Windows) | `bool` | `true` | no |
| <a name="input_instance"></a> [instance](#input\_instance) | n/a | `string` | `""` | no |
| <a name="input_microbedb_key_path"></a> [microbedb\_key\_path](#input\_microbedb\_key\_path) | Path on host to microbedb CVMFS repository pub key | `string` | n/a | yes |
| <a name="input_microbedb_mount_path"></a> [microbedb\_mount\_path](#input\_microbedb\_mount\_path) | Path on host to mount microbedb to share with jobs | `string` | n/a | yes |
| <a name="input_microbedb_path"></a> [microbedb\_path](#input\_microbedb\_path) | Path to microbedb.sqlite as mounted by Galaxy | `string` | `"/cvmfs/microbedb.brinkmanlab.ca/microbedb.sqlite"` | no |
| <a name="input_uwsgi_gid"></a> [uwsgi\_gid](#input\_uwsgi\_gid) | GID of Galaxy process | `number` | `null` | no |
| <a name="input_uwsgi_uid"></a> [uwsgi\_uid](#input\_uwsgi\_uid) | UID of Galaxy process | `number` | `null` | no |

## Outputs

No outputs.

## Resources

| Name | Type |
|------|------|
| [docker_container.microbedb](https://registry.terraform.io/providers/kreuzwerker/docker/latest/docs/resources/container) | resource |
| [docker_image.microbedb](https://registry.terraform.io/providers/kreuzwerker/docker/latest/docs/resources/image) | resource |
| [null_data_source.wait_on_microbedb](https://registry.terraform.io/providers/hashicorp/null/latest/docs/data-sources/data_source) | data source |
<!-- END_TF_DOCS -->