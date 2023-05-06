#!/usr/bin/env bash
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2017            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#* Note that this script uses the mpsfile option to read MPS files. This is  *
#* only available in SAS/OR 15.1 or later. For earlier versions use          *
#* run_sas_legacy.sh which utilizes the sas2sasd macro.                      *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SOLVER=${1}
BINNAME=${2}
NAME=${3}
TIMELIMIT=${4}
MEMLIMIT=${5}   # The memory limit (in MB)
SOLFILE=${6}
THREADS=${7}
MIPGAP=${8}

TMPFILE=check.${SOLVER}.tmp

echo > ${TMPFILE}
echo > ${SOLFILE}

echo "options nodate nonumber pagesize=MAX fullstimer;title1;"                                      >> ${TMPFILE}
echo "filename inzip ZIP \"${NAME}\" GZIP;"                                                         >> ${TMPFILE}
echo "proc optmilp mpsfile=inzip primalout=outfile timetype=real maxtime=${TIMELIMIT} loglevel=3 ABSOBJGAP=${MIPGAP} RELOBJGAP=${MIPGAP};" >> ${TMPFILE}
if [[ ${THREADS} -gt 0 ]] # SAS cannot handle 0 as default number of threads
then
    echo "performance nthreads=${THREADS};"                                                         >> ${TMPFILE}
fi
echo "run;"                                                                                         >> ${TMPFILE}
echo "%put &_OROPTMILP_;"                                                                           >> ${TMPFILE}
echo "%let solution_status = %scan(%substr(&_OROPTMILP_, %index(&_OROPTMILP_, SOLUTION_STATUS)), 2, \" =\");" >> ${TMPFILE}
echo "%let objective = %scan(%substr(&_OROPTMILP_, %index(&_OROPTMILP_, OBJECTIVE)), 2, \" =\");"   >> ${TMPFILE}
if [[ -e $SOLFILE ]]
then
echo "data _null_ ;"                                                                                >> ${TMPFILE}
echo "  set outfile;"                                                                               >> ${TMPFILE}
echo "  format _VALUE_ BEST32.;"							            >> ${TMPFILE}
echo "  file \"${SOLFILE}\";"                                                                       >> ${TMPFILE}
echo "  if _n_=1 then"                                                                              >> ${TMPFILE}
echo "      do;"                                                                                    >> ${TMPFILE}
echo "          if index(\"&solution_status\",\"NOSOL\") > 0 then stop;"                            >> ${TMPFILE}
echo "          if \"&solution_status\" ne \"INFEASIBLE\" then"                                     >> ${TMPFILE}
echo "          	put \"=obj= &objective\";"                                                  >> ${TMPFILE}
echo "          else"                                                                               >> ${TMPFILE}
echo "            do;"                                                                              >> ${TMPFILE}
echo "              put \"=infeas=\";"                                                              >> ${TMPFILE}
echo "              stop;"                                                                          >> ${TMPFILE}
echo "            end;"                                                                             >> ${TMPFILE}
echo "      end;"                                                                                   >> ${TMPFILE}
echo "  put _VAR_  _VALUE_;"                                                                        >> ${TMPFILE}
echo "run;"                                                                                         >> ${TMPFILE}
fi
echo "endsas;"                                                                                      >> ${TMPFILE}

$BINNAME -STDIO -nodms -memsize ${MEMLIMIT}M -autoexec ${TMPFILE}
