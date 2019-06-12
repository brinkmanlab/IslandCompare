# Extract sequence ids and associate to current collection id
function escape_posix(a) {
    gsub("{", "{esc}", a); # Must be first
    gsub("&", "{amp}", a);
    gsub("/", "{slash}", a);
    gsub("\\", "{baskslash}", a);
    gsub("\\?", "{question}", a);
    gsub("%", "{percent}", a);
    gsub("\\*", "{star}", a);
    gsub(":", "{colon}", a);
    gsub("\\|", "{pipe}", a);
    gsub("\"", "{dblquot}", a);
    gsub("<", "{lt}", a);
    gsub(">", "{gt}", a);
    gsub("\\.", "{dot}", a);
    gsub(" ", "{space}", a);
    gsub("\t", "{tab}", a);
    gsub("'", "{quot}", a);
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
}
