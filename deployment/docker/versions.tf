terraform {
  required_providers {
    aws = {
      source = "hashicorp/aws"
    }
    kubernetes = {
      source = "hashicorp/kubernetes"
    }
    galaxy = {
      source = "brinkmanlab/galaxy"
    }
  }
  required_version = ">= 0.13"
}
