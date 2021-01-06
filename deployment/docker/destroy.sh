#!/usr/bin/env bash
# This file is only necessary due to https://github.com/hashicorp/terraform/issues/4149

terraform destroy -target=module.islandcompare -auto-approve
terraform destroy -target=module.admin_user -auto-approve
terraform destroy -target=module.galaxy -auto-approve
terraform destroy -auto-approve
docker ps -a
docker volume list
docker network list

echo 'If docker listed any relevant resources to islandcompare try rerunning ./destroy.sh. If it still fails to delete the resources then you will need to manually delete them.'
echo 'If you manually delete resources you will need to delete your terraform.tfstate file before running ./deploy.sh again.'