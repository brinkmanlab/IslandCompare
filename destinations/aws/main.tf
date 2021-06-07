locals {
  namespace = var.namespace
}

data "aws_region" "current" {}

module "k8s" {
  source = "../k8s"
  data_dir = var.data_dir
  instance = local.instance
  uwsgi_uid = var.uwsgi_uid
  uwsgi_gid = var.uwsgi_gid
  namespace = local.namespace
  microbedb_path = var.microbedb_path
}