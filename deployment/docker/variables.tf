variable "instance" {
  type = string
}

variable "debug" {
  type = bool
}

variable "email" {
  type = string
}

variable "microbedb_mount_path" {
  type = string
  description = "Path on host to mount microbedb to share with jobs"
  default = "./microbedb/mount"
}

variable "enable_CVMFS" {
  type = bool
  description = "Automatically mount CVMFS from a preconfigured Docker container (Not available for OSX and Windows)"
  default = true
}

variable "host_port" {
  type = number
  description = "Host port to expose galaxy service"
}

variable "docker_gid" {
  type = number
  description = "GID with write permission to /var/run/docker.sock"
}