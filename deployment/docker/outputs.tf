output "galaxy_admin_password" {
  value = module.admin_user.password
}

output "galaxy_endpoint" {
  value = module.galaxy.endpoint
}