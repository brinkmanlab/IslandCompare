#!/usr/bin/env -S gawk -f 
tool_input==0 { newick_found=1 }
tool_input==1 && !newick_found { print tool_input_id; nextfile }
