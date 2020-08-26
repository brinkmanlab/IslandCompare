locals {
  name_suffix = var.instance != "" ? "-${var.instance}" : ""
  debug       = true
  region      = "us-west-2"
  instance    = var.instance != "" ? "-${var.instance}" : "default"
  cloud       = module.cloud
}

provider "aws" {
  region = local.region
}

module "cloud" {
  source = "github.com/brinkmanlab/cloud_recipes.git//aws"
  cluster_name = "IslandCompare"
  instance = var.instance
  debug = local.debug
}

module "galaxy_aws" {
  source   = "github.com/brinkmanlab/galaxy-container.git//destinations/aws"
  instance = var.instance
  galaxy_conf = {
    email_from     = var.email
    error_email_to = var.email
    #require_login       = true
    #allow_user_creation = false
  }
  image_tag = var.image_tag
  email     = var.email
  debug     = local.debug
  eks       = local.cloud.eks
  vpc       = local.cloud.vpc
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
  depends_on            = [module.galaxy_aws]
  source                = "./destinations/aws"
  data_dir              = module.galaxy_aws.data_dir
  nfs_server            = module.galaxy_aws.nfs_server
  user_data_volume_name = module.galaxy_aws.user_data_volume_name
}