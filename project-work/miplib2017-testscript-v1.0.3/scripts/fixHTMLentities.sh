#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# import some useful functions that are reused by other scripts
. $(dirname "${BASH_SOURCE[0]}")/run_functions.sh

if [[ $# -eq 0 ]]
then
    echo "Usage: $0 files" >&2
    echo "Replaces all occurrences of HTML entities with the unescaped ones."
    echo "This operation is executed in place, so the files are changed in the process."
    echo
    exit 1
fi

for FILE in $@
do
    unescapeHTMLentities ${FILE}
done
