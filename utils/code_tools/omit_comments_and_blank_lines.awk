# this file counts occurances of the #if !defined(REMOVE_TRACEMANAGER) in a file
# usage in MS6_9 was
# epd@leconte:~/qmcpack$ rg REMOVE_TRACEMANAGER src | awk -e '!seen[$0]++ {print $1}' | sed -E 's/:#if//' | sed -E 's/:#cmakedefine//' | awk -e '! /CMakeLists.txt/ {print}'  | xargs gawk -f utils/code_tools/count_in_ifdef.awk | awk '{total += $2; print} END { printf "total lines: %d\n", total }'

BEGINFILE {
    start_comment = 0;
    end_comment = 0;
    lines = 0
}

in_comment_block == 1 && /\*\// {
    in_comment_block = 0;
    end_comment = FNR;
}

/^[ \t]*\/\*\*/ && (! /\*\//) {
    in_comment_block = 1;
    start_comment = FNR;
}

in_comment_block != 1 && (! (/^[ \t]*\/\// || /^[ \t]*$/)) {
  lines++
}   


ENDFILE {
    printf "%s %d : ", FILENAME, lines;
    printf "\n"
}
