#!/usr/bin/env bash
# This file is only necessary due to https://github.com/hashicorp/terraform/issues/4149

terraform destroy -target=module.islandcompare -auto-approve
terraform destroy -target=module.admin_user -auto-approve
terraform destroy -target=module.galaxy -auto-approve
terraform destroy -target=kubernetes_persistent_volume_claim.microbedb -auto-approve
terraform destroy -target=kubernetes_namespace.instance -auto-approve
terraform destroy -auto-approve