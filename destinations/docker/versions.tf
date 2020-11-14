terraform {
  required_providers {
    external = {
      source = "hashicorp/external"
    }
    null = {
      source = "hashicorp/null"
    }
    random = {
      source = "hashicorp/random"
    }
    docker = {
      source = "terraform-providers/docker"
    }
  }
  required_version = ">= 0.13"
}