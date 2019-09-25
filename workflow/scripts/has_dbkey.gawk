# Output collection id's which have an assigned genome in the metadata input
BEGIN { FS=" "; i=1; }
$1 == "dbkey:" && $2 != "?" { print i++; nextfile }


