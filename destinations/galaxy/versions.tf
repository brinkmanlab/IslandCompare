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
      version = "0.1.0"
    }
  }
  required_version = ">= 0.13"
}