#!/bin/bash

# change into SCIP root directory so that we can run bin/scip instead of ../bin/scip
cd ..

### START SHELL TUTORIAL

# build a fresh version of SCIP
make -j clean
make -j

# run scip with some commands for the shell tutorial
bin/scip < doc/inc/shelltutorial/commands | tee doc/inc/shelltutorial/shelltutorialraw.tmp

# cleanup of files created by the SCIP commands
rm stein27.lp stein27.sol settingsfile.set

cd doc

# modify the raw log file by adding doxygen snippet marker via the python script
python inc/shelltutorial/insertsnippetstutorial.py

### FINISHED SHELL TUTORIAL

### START FAQ GENERATION

cd inc/faq
python parser.py ./ localdoxysubstitutions && php localfaq.php > faq.inc
cd ../../

### FINISHED FAQ GENERATION

### START PARAMETER FILE CREATION

cd ..

bin/scip -c "set default set save doc/inc/parameters.set quit"

cd doc

### FINISHED FAQ GENERATION

# finally build the scip documentation
doxygen scip.dxy
