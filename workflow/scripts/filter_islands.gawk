#!/usr/bin/env gawk -f
# Filter islands based on size specified in 'minimum_size' environment variable
BEGIN { FS=OFS="\t"; }
/^#/ || ($5-$4 >= ENVIRON["minimum_size"]) { print }
