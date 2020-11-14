# This file is only necessary due to https://github.com/hashicorp/terraform/issues/4149

.PHONY: apply
apply:
	terraform apply -target=module.galaxy -auto-approve
	terraform apply -target=module.admin_user -auto-approve
	# -parallelism=1 due to https://github.com/galaxyproject/galaxy/issues/10651
	terraform apply -target=module.islandcompare -auto-approve -parallelism=1
	terraform apply -auto-approve
	# Run make destroy to shutdown and delete everything

.PHONY: destroy
destroy:
	terraform destroy -target=module.islandcompare -auto-approve
	terraform destroy -target=module.admin_user -auto-approve
	terraform destroy -target=module.galaxy -auto-approve
	terraform destroy -auto-approve
