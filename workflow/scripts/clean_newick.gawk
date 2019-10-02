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

{ print unescape_posix(gensub(/'([^']+)\.ref'/, "'\\1'", "g", $0)); nextfile }
