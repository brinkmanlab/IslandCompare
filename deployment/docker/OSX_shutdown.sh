#!/usr/bin/env bash

terraform destroy -target=module.islandcompare.docker_container.microbedb -auto-approve