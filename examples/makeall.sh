#!/usr/bin/env bash

# stop on error
set -e

EXAMPLES=(Binpacking CallableLibrary Eventhdlr GMI LOP MIPSolver Queens Relaxator SCFLP TSP VRP)
LPSOLVERS=(spx2 cpx none)
OPTS=(opt dbg)
LIBTYPE=(static shared)

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

for EXAMPLE in ${EXAMPLES[@]}
do
   echo
   echo
   echo ===== $EXAMPLE =====
   echo
   pushd $EXAMPLE > /dev/null
   echo
   for OPT in ${OPTS[@]}
   do
      for LPS in ${LPSOLVERS[@]}
      do
         for TYPE in ${LIBTYPE[@]}
         do
            if test "$TYPE" = "shared"
            then
               SHAREDVAL="true"
               LIBEXT="so"
            else
               SHAREDVAL="false"
               LIBEXT="a"
            fi

            SCIPLIB=../../lib/$TYPE/libscip.$OSTYPE.$ARCH.gnu.$OPT.$LIBEXT
            LPILIB=../../lib/$TYPE/liblpi${LPS}.$OSTYPE.$ARCH.gnu.$OPT.$LIBEXT
            if test -e $SCIPLIB && test -e $LPILIB
            then
               echo make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL clean
               make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL clean
               echo
               echo make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL
               make OPT=$OPT LPS=$LPS SHARED=$SHAREDVAL
            else
               echo "$SCIPLIB or $LPILIB do not exist - skipping combination ("$OPT", "$LPS", "$TYPE")"
            fi
            echo
         done
      done
   done
   popd > /dev/null
done
