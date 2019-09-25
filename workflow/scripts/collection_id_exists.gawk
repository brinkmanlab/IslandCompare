# If tool_input_id begins with query, output nothing, else output result
match(tool_input_id, "^"ENVIRON["query"]":") { found = 1; exit }
{ nextfile }
END { if (!found) print ENVIRON["result"] }
