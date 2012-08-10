#!/bin/sh

set -e

wgetcmd=wget
wgetcount=`which wget 2>/dev/null | wc -l`
if test ! $wgetcount = 1; then
  echo "Utility wget not found in your PATH."
  if test `uname` = Darwin; then
    wgetcmd=ftp
    echo "Using ftp command instead."
  elif test `uname` = FreeBSD; then
    wgetcmd=fetch
    echo "Using fetch command instead."
  else
    exit -1
  fi
fi

echo " "
echo "Running script for downloading the source code for the ASL"
echo " "

rm -f solvers.tar
echo "Downloading the source code from www.netlib.org..."
$wgetcmd ftp://www.netlib.org/ampl/solvers.tar

rm -rf solvers

echo "Unpacking the source code..."
tar xf solvers.tar

echo "Uncompressing the source files..."
gunzip -fr solvers

#echo "Adding No_dtoa to CFLAGS..."
#cd solvers
#mv makefile.u makefile.u.orig
#sed -e 's/CFLAGS = /CFLAGS = -DNo_dtoa /g' makefile.u.orig > makefile.u 
#cd ..

echo "Deleting the tar file..."
rm solvers.tar

echo " "
echo "Done downloading the source code for ASL."
echo " "

