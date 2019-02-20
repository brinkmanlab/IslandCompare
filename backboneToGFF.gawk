function abs(v) { return v < 0 ? -v : v }
function strand(v) { return v < 0 ? "-" : "+" }
NR>1 {
    for (i=1; i < NF; i+=2) {
        if ($i != "0") {
            for (j=i+2; j < NF; j+=2) {
                if ($j != "0") {
                    print "seq" i-1, "progressiveMauve", "match", abs($i), abs($(i+1)), ".", strand($i), ".", "Target=seq" j-1 " " abs($j) " " abs($(j+1)) " " strand($j)
                }
            }
            break;
        }
    }
}
