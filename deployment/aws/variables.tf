variable "instance" {
  type = string
  description = "Unique deployment instance identifier"
}

variable "debug" {
  type = bool
  default = false
  description = "Enabling will put the deployment into a mode suitable for debugging"
}

variable "email" {
  type = string
  description = "Email address to send automated emails from"
}

variable "region" {
  type = string
  description = "AWS region to deploy into"
}