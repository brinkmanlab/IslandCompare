# Collects outputs from IslandCompare workflow and combines into GFF
# Why not Python? AWK is far more concise and the language is far more stable for this application.

function gff_encode(v) {
    gsub(/\n {21}/, " ", v);
    gsub(/\t/, "%09", v);
    gsub(/,/, "%2C", v);
    gsub(/=/, "%3D", v);
    gsub(/;/, "%3B", v);
    return v;
}

function unescape_posix(a) {
    gsub("{amp}", "&", a);
    gsub("{slash}", "/", a);
    gsub("{baskslash}", "\\", a);
    gsub("{question}", "\\?", a);
    gsub("{percent}", "%", a);
    gsub("{star}", "\\*", a);
    gsub("{colon}", ":", a);
    gsub("{pipe}", "\\|", a);
    gsub("{dblquot}", "\"", a);
    gsub("{lt}", "<", a);
    gsub("{gt}", ">", a);
    gsub("{dot}", "\\.", a);
    gsub("{space}", " ", a);
    gsub("\t", "{tab}", a);
    gsub("{quot}", "'", a);
    gsub("{esc}", "{", a); #Must be last
    return a;
}

function generate_distinct_color(n){
    #TODO implement color generator
    return (n in colors)?colors[n]:"#000000"
}

function abs(v) { return v < 0 ? -v : v }
function int_to_strand(v) { return v < 0 ? "-" : "+" }

# Expected inputs
# 1: Sequence lengths
# 2: MCL Clustering matrix
# 3: Newick tree
# 4: Genomic island GFF file + user annotations
# 5: RGI annotation summary
# 6: Alignment GFF

BEGIN {
    last_tool_input = tool_input;
    tool_input_index = 0;
    min_cluster_size = ENVIRON["min_cluster_size"] ? ENVIRON["min_cluster_size"] : 2;
    split("#e6194b #3cb44b #ffe119 #4363d8 #f58231 #911eb4 #46f0f0 #f032e6 #bcf60c #fabebe #008080 #e6beff #9a6324 #fffac8 #800000 #aaffc3 #808000 #ffd8b1 #000075 #808080 #ffffff #000000", colors, " ");
    print "##gff-version 3";
}

BEGINFILE {
    if (last_tool_input != tool_input) { tool_input_index = 0; } # Catch input change and reset index
    else ++tool_input_index;
    last_tool_input = tool_input;
}

#Output sequence lengths
tool_input==1 { print "##sequence-region " $1 " 1 " $2; next }

#Read clusters
tool_input==2 && (NF >= min_cluster_size) { for (i=1; i<=NF; i++) { clusters[$i] = FNR;} next}

#Output newick
tool_input==3 { print "##newick: "unescape_posix(gensub(/'([^']+)\.ref'/, "'\\1'", "g", $0)); nextfile }

#Output islands with cluster and color
#TODO add offsets from stitching
tool_input==4 && /^[^#]/ {
    i=$1 ":" ($4-1) "-" $5;
    if (i in clusters) {
        if (length($9)>0 && substr($9, length($9)-1) != ";") $9 = $9 ";";
        cluster = clusters[i];
        $9 = $9 "cluster=" cluster ";color=" generate_distinct_color(cluster);
    }
    print;
    next
}

#Output RGI
tool_input==5 && FNR==1 { split( $0, tags ); next }
tool_input==5 { match($1, /ID=[^_]*([^;]+)/, a); print gensub(a[1]" *$", "", 1, $2),"RGI-CARD","gene",$3,$4,$8,$5,".","Name="$9";Alias=ARO:"$11";"tags[6]"="$6";"tags[7]"="$7";"tags[9]"="$9";"tags[10]"="$10";"tags[12]"="$12";"tags[13]"="$13";"tags[14]"="$14";"tags[16]"="$16";"tags[17]"="$17";"tags[15]"="$15";"; next}

#Output alignment
tool_input==6 { print; }

