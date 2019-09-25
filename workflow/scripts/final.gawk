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

function abs(v) { return v < 0 ? -v : v }
function int_to_strand(v) { return v < 0 ? "-" : "+" }

# Expected inputs
# 1: Fasta's used for mauve alignment
# 2: Genbank files to read sequence lengths
# 3: MCL Clustering matrix
# 4: Newick tree
# 5: Genomic island GFF file
# 6: RGI annotation summary
# 7: mauve xfma
# 8: mauve alignment tables
# 9: Genbank files again, to output gene annotations

BEGIN {
    last_tool_input = tool_input;
    tool_input_index = 0;
    min_cluster_size = ENVIRON["min_cluster_size"] ? ENVIRON["min_cluster_size"] : 2;
    print "##gff-version 3";
}

BEGINFILE { 
    if (tool_input==1 || tool_input==8) { FS="\n {21}/"; RS="\n {5}\\<"; } else { FS=OFS="\t"; RS="\n" };
    if (last_tool_input != tool_input) { tool_input_index = 0; } # Catch input change and reset index
    else ++tool_input_index;
    last_tool_input = tool_input;
}

#Read dataset -> sequence id
tool_input==0 && match($0, /^>([^ ]+)/, a) { file=FILENAME; sub(/.*\//, "", file); datasets[file] = a[1]; nextfile }

#Output sequence lengths
tool_input==1 && FNR==1 && match($0, /\nACCESSION +([^\n]+)/, a) { sequence = a[1]; next }
tool_input==1 && match($1, /source +([0-9]+)[^0-9]+([0-9]+)/, a) { print "##sequence-region " sequence " " a[1] " " a[2]; nextfile }

#Read clusters
tool_input==2 && (NF >= min_cluster_size) { for (i=1; i<=NF; i++) { clusters[$i] = FNR;} next}

#Output newick
tool_input==3 { print "##newick: "unescape_posix(gensub(/'([^']+)\.ref'/, "'\\1'", "g", $0)); nextfile }

#Output islands with cluster and color
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

#Read sequence order from XFMA
tool_input==6 && match($1, /^#Sequence([0-9]+)File/, a) { ordinals[tool_input_index, a[1]-1] = datasets[$2]; next}
tool_input==6 && /^#BackboneFile/ { nextfile }

#Output alignment
tool_input==7 && FNR>1 {
    seq=0;
    for (i=1; i < NF; i+=2) {
        if ($i != "0" && (abs($(i+1) - $i) > ENVIRON["minimum_homologous_region"])) {
            base_seq=seq;
            for (j=i+2; j < NF; j+=2) {
                ++seq;
                if ($j != "0" && (abs($(j+1) - $j) > ENVIRON["minimum_homologous_region"])) {
                    print ordinals[tool_input_index, base_seq], "progressiveMauve", "match", abs($i), abs($(i+1)), ".", int_to_strand($i), ".", "Target=" ordinals[tool_input_index, seq] " " abs($j) " " abs($(j+1)) " " int_to_strand($j);
                }
            }
            break;
        }
        ++seq;
    }
    next
}

#Append gene annotations from GBK input files
tool_input==8 { record_type="gene"; locus_tag=""; gene=""; product=""; codon_start="."; }
tool_input==8 && FNR==1 && match($0, /\nACCESSION +([^\n]+)/, a) { sequence = a[1]; next }
tool_input==8 && match($0, /\/gene="([^"]+)/, a) { gene=gff_encode(a[1]); }
tool_input==8 && match($0, /\/pseudo/) { record_type="pseudogene"; }
tool_input==8 && match($0, /\/locus_tag="([^"]+)/, a) { locus_tag=gff_encode(a[1]); }
tool_input==8 && match($0, /\/product="([^"]+)/, a) { product=gff_encode(a[1]); }
tool_input==8 && match($0, /\/codon_start=([0-9])/, a) { codon_start=a[1]-1; }
tool_input==8 && match($1, /(\w+) +(complement\()?([0-9]+)[^0-9]+([0-9]+)/, a) {
    type=a[1]; start=a[3]; stop=a[4]; strand=(a[2]==""?"+":"-");
    print sequence, "Genbank", record_type, start, stop, ".", strand, codon_start, "Name="gene";Type="type";locus_tag="locus_tag";product="product;
    next
}

