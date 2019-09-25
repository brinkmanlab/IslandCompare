BEGIN { OFS=FS="\t" }
tool_input == 0 && match($0, /^>([^ ]+)/, a) {
    names[a[1]] = tool_input_id;
    nextfile;
}
tool_input == 1 && ($1 in names) {
    print names[$1];
    delete names[$1];
}
END {
    # Output elements not found in alignment order
    for (name in names) print names[name];
}
