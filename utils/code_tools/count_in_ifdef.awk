# this file counts occurances of the #if !defined(REMOVE_TRACEMANAGER) in a file
# usage in MS6_9 was
# epd@leconte:~/qmcpack$ rg REMOVE_TRACEMANAGER src | awk -e '!seen[$0]++ {print $1}' | sed -E 's/:#if//' | sed -E 's/:#cmakedefine//' | awk -e '! /CMakeLists.txt/ {print}'  | xargs gawk -f utils/code_tools/count_in_ifdef.awk | awk '{total += $2; print} END { printf "total lines: %d\n", total }'

BEGIN { print "Filename : lines for TM : if !defined(REMOVE_TRACEMANAGER) blocks | block sizes";
    file_count = 0;
}

BEGINFILE {
    delete trace_man;
    start = 0;
    end = 0;
    start_comment = 0;
    end_comment = 0;
    i = 0;
    in_trace_block = 0;
    in_comment_block = 0;
    total_trace = 0;
    delete trace_lines;
    file_count++;
}
/^[ \t]*#if[\t ]+!defined\(REMOVE_TRACEMANAGER\)/ {
    in_trace_block = 1;
    start = FNR;
}

in_trace_block == 1 && in_comment_block == 1 && /\*\// {
    in_comment_block = 0;
    end_comment = FNR;
    trace_lines[i] -= end_comment - start_comment;
}

in_trace_block == 1 && (/^[ \t]*\/\*\*/ && (! /\*\//)) {
    in_comment_block = 1;
    start_comment = FNR;
}

in_trace_block == 1 && (/^[ \t]*\/\// || /^[ \t]*$/) {
    trace_lines[i] -= 1;
}   

in_trace_block == 1 && (/^[ \t]*#else/ || /^[ \t]*#elif/ || /^[ \t]*#endif/) {
    in_trace_block = 0;
    end = FNR;
    total_trace += trace_lines[i] + (end - start);
    trace_lines[i++] = int( end - start );
}



ENDFILE {
    num_ifdefs = length(trace_lines)
    printf "%s %d : %d : ", FILENAME, total_trace, num_ifdefs;
	for( i = 0; i < length(trace_lines); i++) {
	    printf " %d ", trace_lines[i];
	}
    printf "\n"
}

END {
        printf "over %d files\n", file_count
}
