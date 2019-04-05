{ i=0; while (match(substr($0, i+=RSTART+RLENGTH), /'([^']+)'/, a)) { sub(/\.ref$/, "", a[1]); print a[1]; }}
