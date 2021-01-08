locals {
  microbedb_mount_path = "/cvmfs/microbedb.brinkmanlab.ca"
}

provider "aws" {
  region  = var.region
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

resource "kubernetes_namespace" "instance" {
  depends_on = [module.cloud]
  metadata {
    name = var.instance
  }
}

module "cvmfs" {
  depends_on = [module.cloud]
  source = "github.com/brinkmanlab/cloud_recipes.git//util/k8s/cvmfs"
  cvmfs_keys = {
    "microbedb.brinkmanlab.ca" = file("../microbedb.brinkmanlab.ca.pub")
  }
  servers = ["http://stratum-1.sfu.brinkmanlab.ca/cvmfs/@fqrn@", "http://stratum-1.cedar.brinkmanlab.ca/cvmfs/@fqrn@"]
}

resource "kubernetes_persistent_volume_claim" "microbedb" {
  wait_until_bound = false
  metadata {
    name      = "microbedb"
    namespace = kubernetes_namespace.instance.metadata.0.name
  }
  spec {
    access_modes       = ["ReadOnlyMany"]
    storage_class_name = module.cvmfs.storageclasses["microbedb.brinkmanlab.ca"].metadata.0.name
    resources {
      requests = {
        storage = "1Gi"
      }
    }
  }
}

module "galaxy" {
  source   = "github.com/brinkmanlab/galaxy-container.git//destinations/aws"
  instance = var.instance
  namespace = kubernetes_namespace.instance
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
  extra_mounts = {
    "microbedb" = {
      claim_name = kubernetes_persistent_volume_claim.microbedb.metadata.0.name
      path = local.microbedb_mount_path
      read_only = true
    }
  }
  extra_job_mounts = ["${kubernetes_persistent_volume_claim.microbedb.metadata.0.name}:${local.microbedb_mount_path}"]
  visualizations = ["https://github.com/brinkmanlab/multiviz/releases/download/v1.0.0/IslandCompare.tar.gz"]
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
  source                = "../../destinations/aws"
  instance              = var.instance
  namespace             = module.galaxy.namespace
  data_dir              = module.galaxy.data_dir
  nfs_server            = module.galaxy.nfs_server
  user_data_volume_name = module.galaxy.user_data_volume_name
  microbedb_path = "${local.microbedb_mount_path}/microbedb.sqlite"
  eks = module.cloud.eks
  vpc = module.cloud.vpc
  debug = var.debug
}