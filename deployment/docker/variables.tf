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

variable "host_port" {
  type = number
  description = "Host port to expose galaxy service"
}