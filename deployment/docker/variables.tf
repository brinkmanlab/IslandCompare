variable "instance" {
  type = string
  description = "Unique deployment instance identifier"
}

variable "debug" {
  type = bool
  description = "Enabling will put the deployment into a mode suitable for debugging"
}

variable "email" {
  type = string
  description = "Email address to send automated emails from"
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

variable "docker_socket_path" {
  type = string
  description = "Host path to docker socket"
  default = "/var/run/docker.sock"
}