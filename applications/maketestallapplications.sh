#!/usr/bin/env bash
#
# Run all tests of applications, arguments are passed to the make test command.
# Parameter "-q" turns on the quiet mode, i.e., does not output the logging of the programs.
#
echo "Running all tests on applications."

# set this if you want to stop execution of this script when detecting a fail in the applications tests
: ${STOPONFAIL:=no}

# stops the script on error
set -e

APPLICATIONS=$(for f in *;do if [[ -d $f  ]]; then echo $f;fi; done)

LPSOLVERS=(spx2)
OPTS=(dbg)

# parse command line
MAKEARGS=""
QUIET=0
LIBTYPE="static"
LIBEXT="a"
for i in $@
do
   if test "$i" = "-q"
   then
      echo "Quiet mode."
      QUIET=1
   else
      # test whether we want to build shared libraries (need bash "testing" w.r.t. regexp)
      if [[ "$i" =~ SHARED[\s]*=true ]]
      then
         LIBTYPE="shared"
         LIBEXT="so"
      fi
      MAKEARGS="${MAKEARGS} $i"
   fi
done

# determine architecture
ARCH=`uname -m | \
   sed \
   -e 's/sun../sparc/' \
   -e 's/i.86/x86/' \
   -e 's/i86pc/x86/' \
   -e 's/[0-9]86/x86/' \
   -e 's/amd64/x86_64/' \
   -e 's/IP../mips/' \
   -e 's/9000..../hppa/' \
   -e 's/Power\ Macintosh/ppc/' \
   -e 's/00........../pwr4/'`
OSTYPE=`uname -s | tr '[:upper:]' '[:lower:]' | \
   sed \
   -e 's/cygwin.*/cygwin/' \
   -e 's/irix../irix/' \
   -e 's/windows.*/windows/' \
   -e 's/mingw.*/mingw/'`

APPLICATIONLOG=${PWD}/applicationtestsummary.log

# prepare log file
echo "" > ${APPLICATIONLOG}

# If script is exiting on error, then write this to log
trap "echo Last command failed >> ${APPLICATIONLOG}" ERR

# If script is exiting on exit, then output a summary to the log
trap "
echo
echo
echo ===== Summary =====

cat ${APPLICATIONLOG}
" EXIT

# pretest
for OPT in ${OPTS[@]}
do
   for LPS in ${LPSOLVERS[@]}
   do
      LPILIB=../lib/${LIBTYPE}/liblpi${LPS}.${OSTYPE}.${ARCH}.gnu.${OPT}.${LIBEXT}
      if test ! -e ${LPILIB}
      then
         echo "Error: "${LPILIB}" does not exist, please compile SCIP with OPT="${OPT}" and LPS="${LPS}"." >> ${APPLICATIONLOG}
         echo "Error: "${LPILIB}" does not exist, please compile SCIP with OPT="${OPT}" and LPS="${LPS}"."
         exit 1
      fi
      SCIPLIB=../lib/${LIBTYPE}/libscipbase.${OSTYPE}.${ARCH}.gnu.${OPT}.${LIBEXT}
      if test ! -e ${SCIPLIB}
      then
         echo "Error: "${SCIPLIB}" does not exist, please compile SCIP with OPT="${OPT}" and LPS="${LPS}"." >> ${APPLICATIONLOG}
         echo "Error: "${SCIPLIB}" does not exist, please compile SCIP with OPT="${OPT}" and LPS="${LPS}"."
         exit 1
      fi
   done
done

# run tests
for APPLICATION in ${APPLICATIONS}
do
   # See issues #1100 and #1169
   if test ${APPLICATION} = "PolySCIP" || test ${APPLICATION} = "STP";
   then
      continue
   fi
   echo
   echo ===== ${APPLICATION} ===== >> ${APPLICATIONLOG}
   echo ===== ${APPLICATION} =====
   echo
   pushd ${APPLICATION}
   for OPT in ${OPTS[@]}
   do
      for LPS in ${LPSOLVERS[@]}
      do
         echo make OPT=${OPT} LPS=${LPS} ${MAKEARGS} >> ${APPLICATIONLOG}
         make OPT=${OPT} LPS=${LPS} ${MAKEARGS}
         echo
         echo make OPT=${OPT} LPS=${LPS} ${MAKEARGS} test >> ${APPLICATIONLOG}

         if test ${QUIET} = 1
         then
            make OPT=${OPT} LPS=${LPS} ${MAKEARGS} test > /dev/null
         else
            make OPT=${OPT} LPS=${LPS} ${MAKEARGS} test
         fi

         # find most recently changed result file and display it
         if test -d check/results
         then
            RESFILE=`find check/results/check*.res -type f -printf '%T@ %p\n' | sort -n | tail -1 | cut -f2- -d" "`
            if test -e ${RESFILE}
            then
               cat ${RESFILE} >> ${APPLICATIONLOG}

               # exit immediately if there was a fail if STOPONFAIL is yes (||: to avoid error if grep output is empty)
               GREPFAILS=`grep "fail" ${RESFILE}` || :
               if test "${GREPFAILS}" != "" -a "${STOPONFAIL}" = "yes"
               then
                  echo -e "Testing "${APPLICATION}" failed:\n${GREPFAILS}\nsee ${RESFILE} in ${APPLICATION} directory for more details."
                  exit 1
               fi
            fi
         fi
         echo >> ${APPLICATIONLOG}
         echo
      done
   done
   popd > /dev/null
done
