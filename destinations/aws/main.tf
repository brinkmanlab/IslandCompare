locals {
  namespace = var.namespace != null ? var.namespace : kubernetes_namespace.instance[0]
}

data "aws_region" "current" {}

resource "kubernetes_namespace" "instance" {
  count = var.namespace == null ? 1 : 0
  metadata {
    name = local.instance
  }
}

module "k8s" {
  source = "../k8s"
  data_dir = var.data_dir
  endpoint = var.endpoint
  instance = local.instance
  uwsgi_uid = var.uwsgi_uid
  uwsgi_gid = var.uwsgi_gid
  namespace = local.namespace
  admin_api_key = var.admin_api_key
}