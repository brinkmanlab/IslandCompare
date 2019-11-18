# Collects outputs from IslandCompare workflow and combines into GFF
# Why not Python? AWK is far more concise and the language is far more stable for this application.

# See ./color.gawk for a writeup on this function
function generate_distinct_color(n){
    n -= 1; # Clusters are 1-indexed, shift to 0-index
    if (n > 0) {
        #calculate log2 n TODO implement using integer log2 as it is more efficient
        log2n=log(n)/log(2);
        floor=int(log2n);
        ceil=(log2n>floor?floor+1:floor);
        h=((n - 2^floor) * 2 + 1) / (2^ceil);
        h *= 0.9; # Shift color range to exclude circular h values
    } else {
        h=0;
    }
    s = 1;
    v = 1;
    # The following maps into a hsv color gradient
    # Taken from https://gist.github.com/mjackson/5311256
    i = int(h * 6);
    f = h * 6 - i;
    p = v * (1 - s);
    q = v * (1 - f * s);
    t = v * (1 - (1 - f) * s);
    switch (i % 6) {
        case 0: r = v; g = t; b = p; break;
        case 1: r = q; g = v; b = p; break;
        case 2: r = p; g = v; b = t; break;
        case 3: r = p; g = q; b = v; break;
        case 4: r = t; g = p; b = v; break;
        case 5: r = v; g = p; b = q; break;
    }
    return sprintf("#%02X%02X%02X", r * 255, g * 255, b * 255)
}

# Expected inputs
# 0: MCL Clustering matrix
# 1: Newick tree
# 2: Genomic islands + user annotations GFF
# 3: RGI annotation summary
# 4: Alignment GFF

BEGIN {
    OFS=FS="\t";
    print "##gff-version 3";
}

#Read clusters
tool_input==0 { for (i=1; i<=NF; i++) { clusters[$i] = FNR;} next}

#Output newick and column header
tool_input==1 { 
    print "##newick: "$0;
    print "#Sequence ID", "Feature Source", "Feature Type", "Start", "End", "Score", "Strand", "Phase", "Attributes";
    nextfile 
}

#Output user annotations
tool_input==2 && $3 != "genomic_island" && /^[^#]/ { print; }

#Output islands with cluster and color
tool_input==2 && $3 == "genomic_island" {
    i=$1 ":" $4 "-" $5;
    if (i in clusters) {
        if (length($9)>0 && substr($9, length($9)-1) != ";") $9 = $9 ";";
        cluster = clusters[i];
        $9 = $9 "cluster=" cluster ";color=" generate_distinct_color(cluster);
    }
    print;
    next
}

#Output RGI
tool_input==3 && FNR==1 { split( $0, tags ); next }
tool_input==3 { match($1, /ID=[^_]*([^;]+)/, a); print gensub(a[1]" *$", "", 1, $2),"RGI-CARD","gene",$3,$4,$8,$5,".","Name="$9";Alias=ARO:"$11";"tags[6]"="$6";"tags[7]"="$7";"tags[9]"="$9";"tags[10]"="$10";"tags[12]"="$12";"tags[13]"="$13";"tags[14]"="$14";"tags[16]"="$16";"tags[17]"="$17";"tags[15]"="$15";"; next}

#Output alignment
tool_input==4 { print; }

