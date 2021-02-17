#!/usr/bin/env bash

docker run --rm --privileged --pid=host justincormack/nsenter1 -- mount --make-shared /tmp
terraform apply -target=module.islandcompare.docker_container.microbedb -auto-approve