#!/usr/bin/env bash
#
# Run all tests of examples, arguments are passed to the make test command.
# Parameter "-q" turns on the quiet mode, i.e., does not output the logging of the programs.
#

# stop on error
set -e

EXAMPLES=$(for f in *;do if [[ -d $f  ]]; then echo $f;fi; done)

LPSOLVERS=(spx2)
OPTS=(dbg)

echo "Running all tests on examples."

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
      MAKEARGS="$MAKEARGS $i"
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

EXAMPLELOG=${PWD}/exampletestsummary.log

# prepare log file
rm -f $EXAMPLELOG
touch $EXAMPLELOG

# if something is failing, then write Failed to log
trap "echo Last command failed >> $EXAMPLELOG" ERR

# print summary on exit and remove log
trap "
echo
echo
echo ===== Summary =====

cat $EXAMPLELOG
rm -f $EXAMPLELOG
" EXIT

# pretest
for OPT in ${OPTS[@]}
do
   for LPS in ${LPSOLVERS[@]}
   do
      LPILIB=../lib/$LIBTYPE/liblpi$LPS.$OSTYPE.$ARCH.gnu.$OPT.$LIBEXT
      if test ! -e $LPILIB
      then
         echo "Error: "$LPILIB" does not exist, please compile SCIP with OPT="$OPT" and LPS="$LPS"." | tee -a ../applicationtestsummary.log
         exit 1
      fi
      SCIPLIB=../lib/$LIBTYPE/libscip.$OSTYPE.$ARCH.gnu.$OPT.$LIBEXT
      if test ! -e $SCIPLIB
      then
         echo "Error: "$SCIPLIB" does not exist, please compile SCIP with OPT="$OPT" and LPS="$LPS"." | tee -a ../applicationtestsummary.log
         exit 1
      fi
   done
done

# run tests
for EXAMPLE in $EXAMPLES
do
   echo
   echo
   echo ===== $EXAMPLE ===== | tee -a $EXAMPLELOG
   echo
   pushd $EXAMPLE > /dev/null
   for OPT in ${OPTS[@]}
   do
      for LPS in ${LPSOLVERS[@]}
      do
         echo make OPT=$OPT LPS=$LPS $MAKEARGS | tee -a $EXAMPLELOG
         make OPT=$OPT LPS=$LPS $MAKEARGS

         echo
         echo make OPT=$OPT LPS=$LPS $MAKEARGS test | tee -a $EXAMPLELOG
         if test $QUIET = 1
         then
            make OPT=$OPT LPS=$LPS $MAKEARGS test > /dev/null
         else
            make OPT=$OPT LPS=$LPS $MAKEARGS test
         fi

         # find most recently changed result file and display it ("|| :" to ignore error)
         RESFILE=`ls -tr check/results/*.res 2>/dev/null | tail -1` || :
         if [ -n "$RESFILE" ] && [ -e "$RESFILE" ]
         then
            cat $RESFILE >> $EXAMPLELOG
         fi
         echo | tee -a $EXAMPLELOG
      done
   done
   popd > /dev/null
done
