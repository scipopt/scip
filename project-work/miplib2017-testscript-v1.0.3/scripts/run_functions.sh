#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

# formates the time given in ${1} from seconds into HH:MM:SS
function formatTime() {
    TMP=$((${1}))
    RESULT=""
    DIVISORS=(60 60 24)
    for((i=0; i<=2; i++))
    do
        printf -v RESULT "%02d${RESULT}" $(expr ${TMP} % ${DIVISORS[i]})
        # separate the numbers by colons except for the last (HH hours)
        if test $i -lt 2
        then
            RESULT=":${RESULT}"
        fi
        TMP=`expr ${TMP} / ${DIVISORS[i]}`
    done
    if test ${TMP} -gt 0
    then
        RESULT=${TMP}-${RESULT}
    fi
    echo $RESULT
}

# function to strip version of, e.g., scip-3.2... to only scip and scipampl.* to scipampl
function stripversion() {
    NAMENOPATH=$(basename $1)
    # by '%%', Trim the longest match from the end
    NAMENOVERSION=${NAMENOPATH%%-*}
    NAMENOVERSION=${NAMENOVERSION%%\.*}
    echo ${NAMENOVERSION}
}

# function to capitalize a whole string
function capitalize {
    echo "${1}" | tr '[:lower:]' '[:upper:]'
}

# removes the the path and the'.mps.gz' extension from the instance file in ${1}
function stripPathAndInstanceExtension() {
    echo $(basename ${1} ".mps.gz")
}

# removes the path leading to the file in ${1}
function stripPath() {
    echo $(basename ${1})
}

# Permutes the instance under ${1} with the permutation seed ${2}, using the binary ${3}/permute
# The ${4}/permute has to be SCIP version 4.0 (or higher), as relevant settings were renamed.
function permuteInstance() {
    INSTANCEPATH=${1}
    PERMUTATION=${2}
    BINPATH=${3}

    if [[ ! -e ${BINPATH}/permute ]]
    then
        echo "ERROR: permutations can only be run with link 'permute' set in 'bin' folder; see bin/README"
        exit -1
    fi

    INSTANCEDIR=$(dirname ${INSTANCEPATH})
    INSTANCENAME=$(basename ${INSTANCEPATH} ".mps.gz")

    INSTANCE=${INSTANCEDIR}/${INSTANCENAME}\#p${PERMUTATION}.mps

    if [[ ! -e ${INSTANCE}.gz ]]
    then
        bash -c " ${BINPATH}/permute -c  ' set randomization permutationseed ${PERMUTATION} set randomization advanced permutevars TRUE read ${INSTANCEPATH} write problem ${INSTANCE} q' "
        bash -c " gzip ${INSTANCE}"
    fi &> /dev/null
    INSTANCE=${INSTANCE}.gz
    echo ${INSTANCE}
}

# Replaces the "{INSTANCENAME}" with ${2} in the pattern ${1}.
function replaceInPattern() {
    local PATTERN=${1}
    local INSTANCENAME=${2}
    
    local RESULT=$(echo ${PATTERN} | sed "s/{INSTANCENAME}/${INSTANCENAME}/g")
    echo ${RESULT}
}

# Create the meta file ${1} used by ipet.
function writeMetaFile() {
    METAFILE=${1}
    PERMUTATION=${2}
    TESTNAME=${3}
    BINNAME=${4}
    TIMELIMIT=${5}
    MEMLIMITMB=${6}
    THREADS=${7}
    QUEUE=${8}
    EXCLUSIVE=${9}
    
    echo "@Permutation ${PERMUTATION}" > ${METAFILE}
    echo "@TestName ${TESTNAME}" >> ${METAFILE}
    echo "@Solver $(stripversion ${BINNAME})" >> ${METAFILE}
    echo "@BinName ${BINNAME}" >> ${METAFILE}
    echo "@TimeLimit ${TIMELIMIT}" >> ${METAFILE}
    echo "@MemLimit ${MEMLIMITMB}" >> ${METAFILE}
    echo "@Threads ${THREADS}" >> ${METAFILE}
    if [[ -n "${QUEUE}" ]]
    then
        echo "@Queue ${QUEUE}" >> ${METAFILE}
    fi
    if [[ -n "${EXCLUSIVE}" ]]
    then
        echo "@Exclusive ${EXCLUSIVE}" >> ${METAFILE}
    fi
}

# Removes HTML entities from the given file ${1} IN PLACE.
function unescapeHTMLentities() {
    sed -i 's/&amp;/\&/g; s/&lt;/\</g; s/&gt;/\>/g; s/&quot;/\"/g; s/#&#39;/\'"'"'/g; s/&ldquo;/\"/g; s/&rdquo;/\"/g;' ${1}
}
