variable "microbedb_mount_path" {
  type = string
  description = "Path on host to mount microbedb to share with jobs"
}

variable "microbedb_key_path" {
  type = string
  description = "Path on host to microbedb CVMFS repository pub key"
}