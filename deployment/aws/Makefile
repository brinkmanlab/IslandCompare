# This file is only necessary due to https://github.com/hashicorp/terraform/issues/4149

.PHONY: apply
apply:
  terraform apply -target=module.cloud -auto-approve
  terraform apply -target=module.galaxy -auto-approve
  terraform apply -target=module.admin_user -auto-approve
  terraform apply -auto-approve

.PHONY: destroy
destroy:
  terraform destroy -target=module.islandcompare -auto-approve
  terraform destroy -target=module.admin_user -auto-approve
  terraform destroy -target=module.galaxy -auto-approve
  terraform destroy -auto-approve
