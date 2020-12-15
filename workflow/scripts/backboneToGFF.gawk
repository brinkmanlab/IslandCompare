# Convert progressiveMauves backbone output to GFF3 alignment records
function abs(v) { return v < 0 ? -v : v }
function int_to_strand(v) { return v < 0 ? "-" : "+" }

# Read sequence names
tool_input == 0 && match($0, /^>([^ ]+)/, a) {
    file=FILENAME;
    sub(/.*\//, "", file);
    datasets[file] = a[1];
    nextfile
}

# Read sequence order from XFMA
tool_input==1 && match($1, /^#Sequence([0-9]+)File/, a) { ordinals[a[1]-1] = datasets[$2]; next}
tool_input==1 && /^#BackboneFile/ { nextfile }

# Output alignment
tool_input==2 && FNR>1 {
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
