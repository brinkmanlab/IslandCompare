# Count number of contigs and output ids of drafts
BEGINFILE { i=0; }
/^[^#]/ { ++i; }
i > 1 { print tool_input_id; i=0; nextfile }
