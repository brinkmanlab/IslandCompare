provider "aws" {
  region  = var.region
  version = "~> 2.0"
}

module "galaxy" {
  source   = "github.com/brinkmanlab/galaxy-container.git//destinations/docker"
  instance = var.instance
  galaxy_conf = {
    email_from     = var.email
    error_email_to = var.email
    require_login  = true
    builds_file_path = "/data/microbedb/builds.txt"
    tool_data_table_config_path = "/data/microbedb/tool_data_table_conf.xml"
    #allow_user_creation = false
    #cleanup_job = "never"
    slow_query_log_threshold = 500
  }
  image_tag   = "dev"
  admin_users = [var.email]
  email       = var.email
  debug       = var.debug
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
  galaxy_url     = module.galaxy.endpoint
  master_api_key = module.galaxy.master_api_key
  username       = "admin"
}

provider "galaxy" {
  host    = module.galaxy.endpoint
  api_key = module.admin_user.api_key
}

module "islandcompare" {
  depends_on            = [module.galaxy]
  source                = "github.com/brinkmanlab/IslandCompare.git//destinations/docker"
  instance              = var.instance
  data_dir              = module.galaxy.data_dir
  nfs_server            = module.galaxy.nfs_server
  user_data_volume_name = module.galaxy.user_data_volume_name
  admin_api_key = module.admin_user.api_key
  endpoint = module.galaxy.endpoint
  debug = var.debug
}