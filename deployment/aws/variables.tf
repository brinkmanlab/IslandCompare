variable "instance" {
  type = string
}

variable "debug" {
  type = bool
}

variable "email" {
  type = string
}

variable "region" {
  type = string
}

variable "base_url" {
  type        = string
  description = "The externally visible URL for accessing this instance of IRIDA. This key is used by the e-mailer when sending out e-mail notifications (password resets, for example) and embeds this URL directly in the body of the e-mail."
}