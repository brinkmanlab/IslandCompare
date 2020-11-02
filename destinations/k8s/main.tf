locals {
  namespace = var.namespace != null ? var.namespace : kubernetes_namespace.instance[0]
}

resource "kubernetes_namespace" "instance" {
  count = var.namespace == null ? 1 : 0
  metadata {
    name = local.instance
  }
}

module "galaxy" {
  source = "../galaxy"

  admin_api_key = var.admin_api_key
  data_dir = var.data_dir
  endpoint = var.endpoint
  instance = var.instance
  uwsgi_uid = var.uwsgi_uid
  uwsgi_gid = var.uwsgi_gid
  microbedb_path = var.microbedb_path
}