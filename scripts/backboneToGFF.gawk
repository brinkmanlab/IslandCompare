#Read sequence names
tool_input == 1 {
    file=FILENAME; 
    sub(/.*\//, "", file); 
    datasets[file] = tool_input_id;
    nextfile
}

#Read sequence order from XFMA
tool_input==2 && match($1, /^#Sequence([0-9]+)File/, a) { ordinals[a[1]-1] = datasets[$2]; next}
tool_input==2 && /^#BackboneFile/ { nextfile }

#Output alignment
tool_input==3 && FNR>1 {
    seq=0;
    for (i=1; i < NF; i+=2) {
        if ($i != "0" && (abs($(i+1) - $i) > ENVIRON["minimum_homologous_region"])) {
            base_seq=seq;
            for (j=i+2; j < NF; j+=2) {
                ++seq;
                if ($j != "0" && (abs($(j+1) - $j) > ENVIRON["minimum_homologous_region"])) {
                    print ordinals[base_seq], "progressiveMauve", "match", abs($i), abs($(i+1)), ".", int_to_strand($i), ".", "Target=" ordinals[seq] " " abs($j) " " abs($(j+1)) " " int_to_strand($j);
                }
            }
            break;
        }
        ++seq;
    }
    next
}
