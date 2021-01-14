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
      source = "kreuzwerker/docker"
    }
  }
  required_version = ">= 0.13"
}