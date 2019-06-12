match($0, /^>([^ ]+)/, a) { file=FILENAME; sub(/.*\//, "", file); print file, a[1]; nextfile }
