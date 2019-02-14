BEGIN { print "##gff-version 3"; }
match($0, /^([^\t]+)\t(.+)$/, a) { files[a[1]] = a[2]; next }
/^\(.*;$/ { for (i in files) gsub(i, files[i]); print "##newick: " $0; next }
match($0, /^[- ]*([^ ]+)[- ]*([0-9]+)/, a) { print "##sequence-region " a[1] " 1 " a[2]; next }
