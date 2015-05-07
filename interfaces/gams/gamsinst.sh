#!/bin/bash
# Copyright (C) 2011 GAMS Development and others
# All Rights Reserved.
# This file is distributed under the Eclipse Public License.
#
# Author: Michael Bussieck, Stefan Vigerske

gamspath="$1"

if test -z "$gamspath" ; then
  echo "Need to get path of GAMS system as first argument."
  exit 1
fi

# TODO adapt for DOS ?
#gmscmp=${gamspath}/gmscmpNT.txt
gmscmp="${gamspath}/gmscmpun.txt"

#libname="libGamsCoin.dll";
lib=`pwd`"/lib/libgamsscip.so";

gmscmporig="${gmscmp}.orig"
gmscmpbak="${gmscmp}.bak"

if ! test -r "$gmscmp" ;
then
   echo "File $gmscmp not found or not readable, cannot edit."
   exit 1
fi

if ! test -e "$lib"
then
   echo "Solver library $lib not found, cannot install."
   exit 1
fi

echo "Adding or updating entry for SCIPDEV into $gmscmp (pointing to $lib)."

# keep backup of original gmscmpun.txt file
if ! test -r "$gmscmporig"
then
   cp "$gmscmp" "$gmscmporig"
fi

# keep backup of current gmscmpun.txt file
cp -f "$gmscmp" "$gmscmpbak"

awk -vlib="$lib" '
BEGIN {
   fileType      = 111; 
   dictType      = 0; 
   licCodes      = "0001020304"; 
   defaultOkFlag = 1;
   hiddenFlag    = 0;
   # TODO adapt for DOS
   #scriptCmd  = "gmsgennt.cmd";
   #execCmd    = "gmsgennx.exe";
   scriptCmd = "gmsgenus.run";
   execCmd   = "gmsgenux.out";

   written["SCIPDEV"] = 0;
   libid["SCIPDEV"] = "scp";
   dicttype["SCIPDEV"] = 5;
   modeltypes["SCIPDEV"] = "RMIP MIP QCP RMIQCP NLP DNLP RMINLP CNS MIQCP MINLP";

   startBlock = 0;
}

function writeConfig(solverID) {
   print solverID, fileType, dicttype[solverID], licCodes, defaultOkFlag, hiddenFlag, "2", modeltypes[solverID];
   print scriptCmd;
   print execCmd;
   print lib, libid[solverID], "1";
   written[solverID] = 1;
}

(/^*/ || /^ *$/) { print $0 }

/^DEFAULTS/ {
   for( solverID in written )
      if( !written[solverID] )
      {
         writeConfig(solverID)
         print "";
      }
   print;
   next;
}

!(/^*/ || /^ *$/) {
   if( startBlock < 0 )
   {
      startBlock++;
      next;
   }
   if( $1 in written && !written[$1] )
   {
      writeConfig($1)
      startBlock = -($7+1);
      next;
   }
   print;
}
' $gmscmpbak > $gmscmp


#echo "Installing $lib in $gamspath"

#if test -e "${gamspath}/$libname" ;
#then
#   rm -f "${gamspath}/$libname"
#fi
#cp "${libdir}/$libname" "${gamspath}/$libname"
#ln -s "${libdir}/$libname" "${gamspath}/$libname"

# hide libstdc++ and libgfortran - this was only necessary with GAMS < 24.3
#if test -e "${gamspath}/libstdc++.so.6" ; then
#  echo "Moving ${gamspath}/libstdc++.so.6 to ${gamspath}/libstdc++.so.6.hide"
#  mv "${gamspath}/libstdc++.so.6" "${gamspath}/libstdc++.so.6.hide"
#fi
#if test -e "${gamspath}/libgfortran.so.3" ; then
#  echo "Moving ${gamspath}/libgfortran.so.3 to ${gamspath}/libgfortran.so.3.hide"
#  mv "${gamspath}/libgfortran.so.3" "${gamspath}/libgfortran.so.3.hide"
#fi
