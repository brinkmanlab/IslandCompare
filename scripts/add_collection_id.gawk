# Attach collection id to end of each gff record
function gff_encode(v) {
    gsub(/\t/, "%09", v);
    gsub(/,/, "%2C", v);
    gsub(/=/, "%3D", v);
    gsub(/;/, "%3B", v);
    return v;
}
BEGIN { FS=OFS="\t" }
BEGINFILE { sub(/_\d+$/, "", tool_input_id) } # Remove hid from id
/^[^#]/ { $9 += ";userid=" gffencode(tool_input_id); print }
