locals {
  ansible_galaxy = yamldecode(file("${path.root}/galaxy/vars.yml"))
  mail_name      = regex("(?m)^mail.*hostname=(?P<mail_name>[^ ]+)", file("${path.root}/galaxy/inventory.ini")).mail_name
  mail_port      = regex("(?m)^mail.*port=(?P<mail_port>[^ ]+)", file("${path.root}/inventory.ini")).mail_port
  name_suffix    = var.instance != "" ? "-${var.instance}" : ""
  debug          = true
  region = "us-west-2"
  instance = var.instance != "" ? "-${var.instance}" : "default"
  cloud = var.create_cloud ? module.cloud.0 : data.terraform_remote_state.cloud.0
}

provider "aws" {
  region  = local.region
  version = "~> 2.0"
}

module "cloud" {
  count = var.create_cloud ? 1 : 0
  source = "./cloud_recipies/aws"
  cluster_name = "IslandCompare"
  instance = var.instance
  debug = local.debug
}

data "terraform_remote_state" "cloud" {
  count = var.create_cloud ? 0 : 1
  backend = "remote"
  config {
    # TODO
  }
}

module "galaxy_aws" {
  source                  = "./galaxy/destinations/aws"
  uwsgi_port              = local.ansible_galaxy.uwsgi.port
  web_name                = local.ansible_galaxy.containers.web.name
  app_name                = local.ansible_galaxy.containers.app.name
  worker_name             = local.ansible_galaxy.containers.worker.name
  db_name                 = local.ansible_galaxy.containers.db.name
  db_data_volume_name     = local.ansible_galaxy.volumes.db_data.name
  galaxy_root_volume_name = local.ansible_galaxy.volumes.galaxy_root.name
  user_data_volume_name   = local.ansible_galaxy.volumes.user_data.name
  data_dir                = local.ansible_galaxy.paths.data
  root_dir                = local.ansible_galaxy.paths.root
  config_dir              = local.ansible_galaxy.paths.config
  galaxy_app_image        = "brinkmanlab/${local.ansible_galaxy.containers.app.name}"
  galaxy_web_image        = "brinkmanlab/${local.ansible_galaxy.containers.web.name}"
  instance                = var.instance
  galaxy_conf = {
    email_from          = var.email
    error_email_to      = var.email
    #require_login       = true
    #allow_user_creation = false
  }
  image_tag = var.image_tag
  mail_name = local.mail_name
  mail_port = local.mail_port
  email     = var.email
  debug     = local.debug
  eks = local.cloud.eks
  vpc = local.cloud.vpc
}

data "aws_eks_cluster" "cluster" {
  name = local.cloud.eks.cluster_id
}

data "aws_eks_cluster_auth" "cluster" {
  name = local.cloud.eks.cluster_id
}

provider "kubernetes" {
  host                   = data.aws_eks_cluster.cluster.endpoint
  cluster_ca_certificate = base64decode(data.aws_eks_cluster.cluster.certificate_authority.0.data)
  token                  = data.aws_eks_cluster_auth.cluster.token
  load_config_file       = false
}

resource "random_password" "admin_user" {
  length  = 16
  special = false
}

resource "null_resource" "admin_user" {
  depends_on = [module.galaxy_aws]
  provisioner "local-exec" {
    command = "curl -H 'Content-Type: application/json' -d '${jsonencode({
    username = "brinkmanlab"
    password = random_password.admin_user.result
    email    = "brinkman-ws@sfu.ca"
  })}' /api/users?key=${module.galaxy_aws.master_api_key}"
  }
}

data "external" "api_key" {
  program = ["curl", "--user", "brinkman-ws@sfu.ca:${random_password.admin_user.result}", "--", "${module.galaxy_aws.endpoint}/api/authenticate/baseauth"]
}

module "islandcompare_aws" {
  depends_on = [module.galaxy_aws]
  source = "./destinations/aws"
  data_dir = local.ansible_galaxy.paths.data
  nfs_server = module.galaxy_aws.nfs_server
  user_data_volume_name = local.ansible_galaxy.volumes.user_data.name
}