data "null_data_source" "wait_on_microbedb" {
  count = var.enable_CVMFS ? 1 : 0
  depends_on = [docker_container.microbedb]
  inputs = {
    microbedb_path = var.microbedb_path
  }
}

module "galaxy" {
  source = "../galaxy"

  data_dir = var.data_dir
  instance = var.instance
  uwsgi_uid = var.uwsgi_uid
  uwsgi_gid = var.uwsgi_gid
  microbedb_path = var.enable_CVMFS ? data.null_data_source.wait_on_microbedb.0.outputs.microbedb_path : var.microbedb_path
}