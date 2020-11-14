module "galaxy" {
  source = "../galaxy"
  depends_on = [docker_container.microbedb]

  data_dir = var.data_dir
  instance = var.instance
  uwsgi_uid = var.uwsgi_uid
  uwsgi_gid = var.uwsgi_gid
  microbedb_path = var.microbedb_path
}