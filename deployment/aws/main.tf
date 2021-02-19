locals {
  microbedb_mount_path = "/cvmfs/microbedb.brinkmanlab.ca"
}

provider "aws" {
  region = var.region
}

module "cloud" {
  source = "github.com/brinkmanlab/cloud_recipes.git//aws"
  #?ref=v0.1.2"
  cluster_name = var.instance
  autoscaler_version = "1.20.0"
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
  host = data.aws_eks_cluster.cluster.endpoint
  cluster_ca_certificate = base64decode(data.aws_eks_cluster.cluster.certificate_authority.0.data)
  token = data.aws_eks_cluster_auth.cluster.token
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
  servers = [
    "http://stratum-1.sfu.brinkmanlab.ca/cvmfs/@fqrn@",
    "http://stratum-1.cedar.brinkmanlab.ca/cvmfs/@fqrn@",
  ]
}

resource "kubernetes_persistent_volume_claim" "microbedb" {
  wait_until_bound = false
  metadata {
    name = "microbedb"
    namespace = kubernetes_namespace.instance.metadata.0.name
  }
  spec {
    access_modes = ["ReadOnlyMany"]
    storage_class_name = module.cvmfs.storageclasses["microbedb.brinkmanlab.ca"].metadata.0.name
    resources {
      requests = {
        storage = "1Gi"
      }
    }
  }
}

module "galaxy" {
  source = "github.com/brinkmanlab/galaxy-container.git//destinations/aws"
  instance = var.instance
  namespace = kubernetes_namespace.instance
  galaxy_conf = {
    email_from = var.email
    error_email_to = var.email
    require_login = true
    #allow_user_creation = false
    #cleanup_job = "never"
    slow_query_log_threshold = 500
    maximum_workflow_jobs_per_scheduling_iteration = 1100
    history_local_serial_workflow_scheduling = true
  }
  image_tag = "dev"
  admin_users = [
    var.email]
  email = var.email
  debug = var.debug
  eks = module.cloud.eks
  vpc = module.cloud.vpc
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
  visualizations = [
    "https://github.com/brinkmanlab/multiviz/releases/download/v1.0.0/IslandCompare.tar.gz"]
  # Deal with MCM locking up and k8s runner cant set max time except per runner
  plugins = <<EOF
<plugin id="k8s-mcm" type="runner" load="galaxy.jobs.runners.kubernetes:KubernetesJobRunner">
  <param id="k8s_persistent_volume_claims" from_environ="K8S_VOLUMES">galaxy-user-data:/data</param>
  <param id="k8s_namespace" from_environ="K8S_NAMESPACE">default</param>
  <param id="k8s_use_service_account">True</param>
  <param id="k8s_galaxy_instance_id" from_environ="K8S_NAMESPACE">default</param>
  <param id="k8s_affinity">nodeAffinity:
requiredDuringSchedulingIgnoredDuringExecution:
  nodeSelectorTerms:
  - matchExpressions:
    - key: WorkClass
      operator: In
      values:
      - compute
</param>
  <param id="k8s_cleanup_job" from_environ="K8S_CLEANUP_JOB">always</param>
  <param id="k8s_run_as_user_id" from_environ="K8S_RUN_AS_UID">1000</param>
  <param id="k8s_run_as_group_id" from_environ="K8S_RUN_AS_GID">1000</param>
  <param id="k8s_fs_group_id">1000</param>
  <param id="k8s_walltime_limit">14400</param>
</plugin>
  EOF
  job_destinations = <<EOF
<destination id="mcm" runner="k8s-mcm">
    <param id="requests_cpu">1</param>
    <param id="requests_memory">2G</param>
    <param id="limits_cpu">2</param>
    <param id="limits_memory">4G</param>
    <param id="enabled" from_environ="K8S_ENABLED">False</param>
    <param id="tmp_dir">True</param>
    <param id="docker_enabled">True</param>
    <param id="docker_repo_default" from_environ="K8S_DEFAULT_REPO">docker.io</param>
    <param id="docker_owner_default" from_environ="K8S_DEFAULT_OWNER">brinkmanlab</param>
    <param id="docker_image_default" from_environ="K8S_DEFAULT_IMAGE">galaxy-app</param>
    <param id="docker_tag_default" from_environ="K8S_DEFAULT_IMAGE_TAG">latest</param>
    <resubmit condition="(unknown_error or walltime_reached) and attempt &lt;= 4" delay="attempt * 15" />
    <resubmit condition="memory_limit_reached" destination="c1m4" />
</destination>
  EOF
}

module "admin_user" {
  source = "github.com/brinkmanlab/galaxy-container.git//modules/bootstrap_admin"
  email = var.email
  galaxy_url = "http://${module.galaxy.endpoint}"
  master_api_key = module.galaxy.master_api_key
  username = "admin"
}

provider "galaxy" {
  host = "http://${module.galaxy.endpoint}"
  apikey = module.admin_user.api_key
  wait_for_host = 90
}

module "islandcompare" {
  source = "../../destinations/aws"
  data_dir = module.galaxy.data_dir
  debug = var.debug
  eks = module.cloud.eks
  instance = var.instance
  microbedb_path = "${local.microbedb_mount_path}/microbedb.sqlite"
  namespace = module.galaxy.namespace
  nfs_server = module.galaxy.nfs_server
  user_data_volume_name = module.galaxy.user_data_volume_name
  vpc = module.cloud.vpc
}