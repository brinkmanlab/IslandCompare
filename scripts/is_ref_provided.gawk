# If the reference is set, print all input ids
BEGIN { if (ENVIRON["reference"] == "none" || ENVIRON["reference"] == "") exit 0; }
{ print tool_input_id; nextfile }
