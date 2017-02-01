#!/bin/bash

cd ..

# build a fresh version of SCIP
make -j clean
make -j

# run scip with some commands for the shell tutorial
bin/scip < doc/inc/shelltutorial/commands | tee doc/inc/shelltutorial/shelltutorialraw.txt

# cleanup of files created by the SCIP commands
rm stein27.lp
rm stein27.sol
rm _settingsfilename_
rm settingsfile.set

cd doc

# modify the raw log file by adding doxygen snippet marker via the python script
python inc/shelltutorial/insertsnippetstutorial.py

# finally build the scip documentation
doxygen scip.dxy
