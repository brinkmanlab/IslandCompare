locals {
  namespace = var.namespace
}

module "galaxy" {
  source = "../galaxy"

  data_dir = var.data_dir
  instance = var.instance
  uwsgi_uid = var.uwsgi_uid
  uwsgi_gid = var.uwsgi_gid
  microbedb_path = var.microbedb_path
}