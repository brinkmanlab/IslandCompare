BEGIN {
    split("/,\\,?,%,*,:,|,\",<,>,., ,'", bad_chars, ",");
}

match($0, /^>([^ ]+)/, a) {
    name=a[1];
    for (c in bad_chars) gsub("\\"bad_chars[c], "", name);
    print name; 
    nextfile 
}
