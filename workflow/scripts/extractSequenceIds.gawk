# Extract sequence ids and associate to current collection id
BEGIN { OFS=FS="\t" }
function escape_posix(a) {
    gsub("_", "_esc_", a); # Must be first
    gsub("&", "_amp_", a);
    gsub("/", "_slash_", a);
    gsub("\\", "_backslash_", a);
    gsub("\\?", "_question_", a);
    gsub("%", "_percent_", a);
    gsub("\\*", "_star_", a);
    gsub(":", "_colon_", a);
    gsub("\\|", "_pipe_", a);
    gsub("\"", "_dblquot_", a);
    gsub("<", "_lt_", a);
    gsub(">", "_gt_", a);
    gsub("\\.", "_dot_", a);
    gsub(" ", "_space_", a);
    gsub("\t", "_tab_", a);
    gsub("'", "_quot_", a);
    return a;
}

match($0, /^>([^ ]+)/, a) {
    name=escape_posix(a[1]);
    if (name in names) nextfile;
    names[name] = true;
    print tool_input_id, name;
    nextfile 
}

function unescape_posix(a) {
    gsub("_amp_", "&", a);
    gsub("_slash_", "/", a);
    gsub("_backslash_", "\\", a);
    gsub("_question_", "?", a);
    gsub("_percent_", "%", a);
    gsub("_star_", "*", a);
    gsub("_colon_", ":", a);
    gsub("_pipe_", "|", a);
    gsub("_dblquot_", "\"", a);
    gsub("_lt_", "<", a);
    gsub("_gt_", ">", a);
    gsub("_dot_", ".", a);
    gsub("_space_", " ", a);
    gsub("_tab_", "\t", a);
    gsub("_quot_", "'", a);
    gsub("_esc_", "_", a); #Must be last
    return a;
}
