#!/usr/bin/env -S gawk -f 
# Stitch input multi-fastas in order, padding with 'N'.

BEGIN {
    #Schedule inputs for two passes
    ARGV[ARGC] = "pass=sequence";
    for (i=1; i<ARGC; ++i)
        ARGV[ARGC+i] = ARGV[i]
    ARGC = (ARGC-1) * 2 + 2;

    pass = "annotations";
    FS=" "; # Use for parsing fasta header
    pos=1; # 1-index of stitched sequence
}

BEGINFILE {
    first=1;
    quality=0;
    if (!stitched_id) {
        stitched_id = FILENAME;
        # get basename of FILENAME
        sub(/^.*\//, "", stitched_id);
        #sub(/\.[^.]*$/, "", stitched_id);
        # spaces have special meaning in fasta
        gsub(" ", "_", stitched_id);
    }
    
    OFS="\t";
    ORS="\n";
    if (pass == "annotations") {
        print "##gff-version 3";
        OFS="\t";
        ORS="\n";
    } else if (pass == "sequence") {
        print "##FASTA";
        print ">" stitched_id;
        ORS="";
        OFS="";
    }
}

# Reset quality state
/^>/ {
    quality=0;
}

# Store sequence information
pass == "annotations" && /^>/ {
    if (current_id) {
        # Output record
        print stitched_id, "stitch.gawk", "contig", start, pos-1, ".", ".", ".", "ID="current_id;
    }
    current_id=substr($1, 2);
    if (first) {
        first=0;
    } else {
        pos += ENVIRON["padding_length"];
    }
    start = pos;
    next;
}

pass == "annotations" && !quality {
    pos += length();
}

# Output sequence
pass == "sequence" && /^>/ {
    if (first) {
        first=0;
    } else {
        for (i = 0; i < ENVIRON["padding_length"]; ++i)
            if (pass == "sequence")
                print "N";
    }
    next;
}

# Default print if not quality lines
pass == "sequence" && !quality {
    print;
}

# Allow for fastq input
/^+/ { quality=1; next }

ENDFILE {
    pass = next_pass[pass]
    if (pass) {
        # Set next file to this one for second pass

        
    }
}
