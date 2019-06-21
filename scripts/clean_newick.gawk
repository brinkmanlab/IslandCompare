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

{ print unescape_posix(gensub(/'([^']+)\.ref'/, "'\\1'", "g", $0)); nextfile }
