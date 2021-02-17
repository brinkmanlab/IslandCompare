#!/usr/bin/env bash

# https://github.com/docker/for-mac/issues/3431

docker run --rm --privileged --pid=host justincormack/nsenter1 -- mount --verbose --make-shared /tmp
terraform apply -target=module.islandcompare.docker_container.microbedb -auto-approve