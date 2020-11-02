locals {
  instance = var.instance == "" ? "default" : var.instance
}

variable "instance" {
  type    = string
  default = ""
}

variable "uwsgi_uid" {
  type        = number
  default     = null
  description = "UID of Galaxy process"
}

variable "uwsgi_gid" {
  type        = number
  default     = null
  description = "GID of Galaxy process"
}

variable "admin_api_key" {
  type = string
}

variable "endpoint" {
  type = string
}

variable "data_dir" {
  type = string
}

variable "microbedb_path" {
  type = string
  default = "/cvmfs/microbedb.brinkmanlab.ca/microbedb.sqlite"
  description = "Path to microbedb.sqlite as mounted by Galaxy"
}