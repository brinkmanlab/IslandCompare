output "galaxy_admin_password" {
  value = module.admin_user.password
  sensitive = true
}

output "galaxy_admin_api_key" {
  value = module.admin_user.api_key
  sensitive = true
}

output "galaxy_master_api_key" {
  value = module.galaxy.master_api_key
  sensitive = true
}

output "galaxy_endpoint" {
  value = module.galaxy.endpoint
}