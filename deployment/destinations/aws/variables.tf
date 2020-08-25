locals {
  instance = var.instance == "" ? "default" : var.instance
}

variable "instance" {
  type = string
  default = ""
}

variable "data_dir" {
  type = string
}

variable "nfs_server" {
  type = string
}

variable "user_data_volume_name" {
  type = string
}