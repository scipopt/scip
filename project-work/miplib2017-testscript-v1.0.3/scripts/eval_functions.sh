#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# Reads the value of the key ${2} from file ${1}.
# Expected format:
# @{KEY} {VALUE}
function readMetaData() {
    METAFILE=${1}
    KEY=${2}
    
    echo $(cat ${METAFILE} | sed -n "s/^@${KEY}\s\+\(.\+\)$/\1/p" )
}
