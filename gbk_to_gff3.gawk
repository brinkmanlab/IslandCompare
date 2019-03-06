function encode(v) {
    gsub(/\n {21}/, " ", v);
    gsub(/\t/, "%09", v);
    gsub(/,/, "%2C", v);
    gsub(/=/, "%3D", v);
    gsub(/;/, "%3B", v);
    return v;
}

BEGIN { FS="\n {21}/"; RS="\n {5}\\<"; print "##gff-version 3"; }
FNR==1 && match($0, /\nACCESSION +([^\n]+)/, a) { sequence = a[1]; next }
match($1, /source +([0-9]+)[^0-9]+([0-9]+)/, a) { print "##sequence-region " sequence " " a[1] " " a[2]; next }

{ locus_tag=""; gene=""; product=""; codon_start="."; }
match($0, /\/gene="([^"]+)/, a) { gene=encode(a[1]); }
match($0, /\/locus_tag="([^"]+)/, a) { locus_tag=encode(a[1]); }
match($0, /\/product="([^"]+)/, a) { product=encode(a[1]); }
match($0, /\/codon_start=([0-9])/, a) { codon_start=a[1]-1; }
match($1, /(\w+) +(complement\()?([0-9]+)[^0-9]+([0-9]+)/, a) {
    type=a[1]; start=a[3]; stop=a[4]; strand=(a[2]==""?"+":"-");
    print sequence, "Genbank", "gene", start, stop, ".", strand, codon_start, "Name="gene";Type="type";locus_tag="locus_tag";product="product;
    next
}
