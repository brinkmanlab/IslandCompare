output "admin_password" {
  value = random_password.admin_user.result
}