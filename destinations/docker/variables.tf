variable "microbedb_mount_path" {
  type = string
  description = "Path on host to mount microbedb to share with jobs"
}

variable "microbedb_key_path" {
  type = string
  description = "Path on host to microbedb CVMFS repository pub key"
}

variable "enable_CVMFS" {
  type = bool
  description = "Automatically mount CVMFS from a preconfigured Docker container (Not available for OSX and Windows)"
  default = true
}