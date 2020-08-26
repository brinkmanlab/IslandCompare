variable "instance" {
  type        = string
  default     = "islandcompare"
  description = "Specify a unique instance name for this deployment"
}

variable "image_tag" {
  type    = string
  default = "latest"
}

variable "email" {
  type        = string
  default     = "nolan_w@sfu.ca"
  description = "Email address to send automated emails from"
}

variable "create_cloud" {
  type        = bool
  default     = true
  description = "Create underlying cloud resources"
}