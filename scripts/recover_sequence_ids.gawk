# Replace all sequence IDs in second tool input with first sequence ID in first tool input
BEGIN { FS=OFS="\t"; }
tool_input == 0 && /^[^#]/ { id = $1; nextfile }
tool_input == 1 && /^[^#]/ { $1 = id }
tool_input == 1 { print }
