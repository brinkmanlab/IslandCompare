terraform {
  required_providers {
    null = {
      source = "hashicorp/null"
    }
    random = {
      source = "hashicorp/random"
    }
    galaxy = {
      source = "brinkmanlab/galaxy"
    }
  }
  required_version = ">= 0.13"
}