# Expects input from "Collection Information" tool
# Output id of any collection length > 1
BEGIN { FS="\t"; }
NR >= 1 && $3 > 1 { print $1 }
