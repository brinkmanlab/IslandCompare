# Galaxy Deployment Module
<!-- BEGIN_TF_DOCS -->
## Providers

| Name | Version |
|------|---------|
| <a name="provider_galaxy"></a> [galaxy](#provider\_galaxy) | n/a |

## Modules

No modules.

## Inputs

| Name | Description | Type | Default | Required |
|------|-------------|------|---------|:--------:|
| <a name="input_data_dir"></a> [data\_dir](#input\_data\_dir) | n/a | `string` | n/a | yes |
| <a name="input_debug"></a> [debug](#input\_debug) | n/a | `bool` | `false` | no |
| <a name="input_instance"></a> [instance](#input\_instance) | n/a | `string` | `""` | no |
| <a name="input_microbedb_path"></a> [microbedb\_path](#input\_microbedb\_path) | Path to microbedb.sqlite as mounted by Galaxy | `string` | `"/cvmfs/microbedb.brinkmanlab.ca/microbedb.sqlite"` | no |
| <a name="input_uwsgi_gid"></a> [uwsgi\_gid](#input\_uwsgi\_gid) | GID of Galaxy process | `number` | `null` | no |
| <a name="input_uwsgi_uid"></a> [uwsgi\_uid](#input\_uwsgi\_uid) | UID of Galaxy process | `number` | `null` | no |

## Outputs

No outputs.

## Resources

| Name | Type |
|------|------|
| [galaxy_history.data_managers](https://registry.terraform.io/providers/brinkmanlab/galaxy/latest/docs/resources/history) | resource |
| [galaxy_job.microbedb](https://registry.terraform.io/providers/brinkmanlab/galaxy/latest/docs/resources/job) | resource |
| [galaxy_job.rgi](https://registry.terraform.io/providers/brinkmanlab/galaxy/latest/docs/resources/job) | resource |
| [galaxy_repository.islandcompare](https://registry.terraform.io/providers/brinkmanlab/galaxy/latest/docs/resources/repository) | resource |
| [galaxy_repository.microbedb](https://registry.terraform.io/providers/brinkmanlab/galaxy/latest/docs/resources/repository) | resource |
| [galaxy_stored_workflow.islandcompare](https://registry.terraform.io/providers/brinkmanlab/galaxy/latest/docs/resources/stored_workflow) | resource |
| [galaxy_workflow_repositories.islandcompare](https://registry.terraform.io/providers/brinkmanlab/galaxy/latest/docs/data-sources/workflow_repositories) | data source |
<!-- END_TF_DOCS -->