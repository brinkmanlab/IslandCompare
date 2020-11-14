variable "nfs_server" {
  type = string
}

variable "user_data_volume_name" {
  type = string
}

variable "eks" {
  description = "Instance of EKS module output state"
}

variable "vpc" {
  description = "Instance of VPC module output state"
}

variable "namespace" {
  description = "Instance of kubernetes_namespace to provision instance resources under"
}