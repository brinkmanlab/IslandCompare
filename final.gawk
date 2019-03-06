function encode(v) {
    gsub(/\n {21}/, " ", v);
    gsub(/\t/, "%09", v);
    gsub(/,/, "%2C", v);
    gsub(/=/, "%3D", v);
    gsub(/;/, "%3B", v);
    return v;
}

function generate_distinct_color(n){
  return (n in colors)?colors[n]:"#000000"
}

function abs(v) { return v < 0 ? -v : v }
function int_to_strand(v) { return v < 0 ? "-" : "+" }

BEGINFILE { if (tool_input==2 || tool_input==9) { FS="\n {21}/"; RS="\n {5}\\<"; } else { FS=OFS="\t"; RS="\n" } }

BEGIN {
    min_cluster_size = ENVIRON['min_cluster_size'] ? ENVIRON['min_cluster_size'] : 2;
    split("#e6194b #3cb44b #ffe119 #4363d8 #f58231 #911eb4 #46f0f0 #f032e6 #bcf60c #fabebe #008080 #e6beff #9a6324 #fffac8 #800000 #aaffc3 #808000 #ffd8b1 #000075 #808080 #ffffff #000000", colors, " ");
    print "##gff-version 3";
}

#Read dataset -> sequence id
tool_input==1 && match($0, /^>([^ ]+)/, a) { file=FILENAME; sub(/.*\//, "", file); datasets[file] = a[1]; nextfile }

#Output sequence lengths
tool_input==2 && FNR==1 && match($0, /\nACCESSION +([^\n]+)/, a) { sequence = a[1]; next }
tool_input==2 && match($1, /source +([0-9]+)[^0-9]+([0-9]+)/, a) { print "##sequence-region " sequence " " a[1] " " a[2]; nextfile }

#Read clusters
tool_input==3 && (NF >= min_cluster_size) { for (i=1; i<=NF; i++) { clusters[$i] = FNR;} next}

#Output newick
# First gensub deals with ParSNP stuffing '.ref' at the end of the random dataset chosen
tool_input==4 { $0 = gensub(/(dataset_[0-9]+\.dat)\.ref/, "\\1", "g"); for (i in datasets) gsub(i, datasets[i], $0); print "##newick: " $0; nextfile }

#Output islands with cluster and color
tool_input==5 && /^[^#]/ {
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
tool_input==6 && FNR==1 { split( $0, tags ); next }
tool_input==6 { match($1, /ID=[^_]*([^;]+)/, a); print gensub(a[1]" *$", "", 1, $2),"RGI-CARD","gene",$3,$4,$8,$5,".","Name="$9";Alias=ARO:"$11";"tags[6]"="$6";"tags[7]"="$7";"tags[9]"="$9";"tags[10]"="$10";"tags[12]"="$12";"tags[13]"="$13";"tags[14]"="$14";"tags[16]"="$16";"tags[17]"="$17";"tags[15]"="$15";"; next}

#Read sequence order from XFMA
tool_input==7 && match($1, /^#Sequence([0-9]+)File/, a) { ordinals[a[1]-1] = datasets[$2]; next}
tool_input==7 && /^#BackboneFile/ { nextfile }

#Output alignment
tool_input==8 && FNR>1 {
    seq=0;
    for (i=1; i < NF; i+=2) {
        if ($i != "0" && (abs($(i+1) - $i) > ENVIRON['minimum_homologous_region'])) {
            base_seq=seq;
            for (j=i+2; j < NF; j+=2) {
                ++seq;
                if ($j != "0" && (abs($(j+1) - $j) > ENVIRON['minimum_homologous_region'])) {
                    print ordinals[base_seq], "progressiveMauve", "match", abs($i), abs($(i+1)), ".", int_to_strand($i), ".", "Target=" ordinals[seq] " " abs($j) " " abs($(j+1)) " " int_to_strand($j);
                }
            }
            break;
        }
        ++seq;
    }
    next
}

#Append gene annotations from GBK input files
tool_input==9 { locus_tag=""; gene=""; product=""; codon_start="."; }
tool_input==9 && FNR==1 && match($0, /\nACCESSION +([^\n]+)/, a) { sequence = a[1]; next }
tool_input==9 && match($0, /\/gene="([^"]+)/, a) { gene=encode(a[1]); }
tool_input==9 && match($0, /\/locus_tag="([^"]+)/, a) { locus_tag=encode(a[1]); }
tool_input==9 && match($0, /\/product="([^"]+)/, a) { product=encode(a[1]); }
tool_input==9 && match($0, /\/codon_start=([0-9])/, a) { codon_start=a[1]-1; }
tool_input==9 && match($1, /(\w+) +(complement\()?([0-9]+)[^0-9]+([0-9]+)/, a) {
    type=a[1]; start=a[3]; stop=a[4]; strand=(a[2]==""?"+":"-");
    print sequence, "Genbank", "gene", start, stop, ".", strand, codon_start, "Name="gene";Type="type";locus_tag="locus_tag";product="product;
    next
}

