provider "aws" {
  region  = var.region
  version = "~> 2.0"
}

module "cloud" {
  source             = "github.com/brinkmanlab/cloud_recipes.git//aws" #?ref=v0.1.2"
  cluster_name       = var.instance
  autoscaler_version = "1.17.3"
  #docker_registry_proxies = {
  #  quay = {
  #    hostname = "quay.io"
  #    url = "https://quay.io"
  #    username = var.quay_io_user
  #    password = var.quay_io_password
  #  }
  #}
  debug = var.debug
}

data "aws_eks_cluster" "cluster" {
  name = module.cloud.eks.cluster_id
}

data "aws_eks_cluster_auth" "cluster" {
  name = module.cloud.eks.cluster_id
}

provider "kubernetes" {
  host                   = data.aws_eks_cluster.cluster.endpoint
  cluster_ca_certificate = base64decode(data.aws_eks_cluster.cluster.certificate_authority.0.data)
  token                  = data.aws_eks_cluster_auth.cluster.token
  load_config_file       = false
}

module "galaxy" {
  source   = "github.com/brinkmanlab/galaxy-container.git//destinations/aws"
  instance = var.instance
  galaxy_conf = {
    email_from     = var.email
    error_email_to = var.email
    require_login  = true
    #allow_user_creation = false
    #cleanup_job = "never"
    slow_query_log_threshold = 500
  }
  image_tag   = "dev"
  admin_users = [var.email]
  email       = var.email
  debug       = var.debug
  eks         = module.cloud.eks
  vpc         = module.cloud.vpc
  tool_mappings = {
    "progressivemauve": "c1m2"
    "rgi": "c4m8"
    "tp_awk_tool": "tiny"
    "tp_text_file_with_recurring_lines": "tiny"
    "sigihmm": "c4m8"
    "feature-merge": "c1m8"
    "islandpath": "c1m2"
    "mash": "c4m8"
    "parsnp": "c16m32"
    "mcl": "c4m8"
    "biopython-convert": "c1m2"
    "extract-tree-order": "tiny"
    "mauve-contig-mover": "c1m2"
    "mauve-contig-mover-stitch": "c1m1"
    "make-unique-id": "tiny"
    "sendmail": "local-tiny"
    "sha256sum": "tiny"
    "awkscript": "tiny"
  }
}

module "admin_user" {
  source         = "github.com/brinkmanlab/galaxy-container.git//modules/bootstrap_admin"
  email          = var.email
  galaxy_url     = "http://${module.galaxy.endpoint}"
  master_api_key = module.galaxy.master_api_key
  username       = "admin"
}

provider "galaxy" {
  host    = "http://${module.galaxy.endpoint}"
  api_key = module.admin_user.api_key
}

module "islandcompare" {
  depends_on            = [module.galaxy]
  source                = "github.com/brinkmanlab/IslandCompare.git//destinations/aws"
  instance              = var.instance
  data_dir              = module.galaxy.data_dir
  nfs_server            = module.galaxy.nfs_server
  user_data_volume_name = module.galaxy.user_data_volume_name
  debug = var.debug
}