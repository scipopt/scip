#!/bin/bash

#
# generate doxygen documentation for SCIP
# requires python and php in PATH
#
# Optionally, a custom .dxy file can be passed for the doxygen configuration
#

: ${DOCLOG=${PWD}/builddoc.log}

echo "" > ${DOCLOG}
echo "running builddoc.sh in $(pwd)"

# stop on error
set -e

if [ -e MathJax ]; then
   echo "Found MathJax in doc folder."
   export MATHJAX_RELPATH="../MathJax"
else
   echo "No local MathJax found, will build doc with online version."
fi

if [ "$1" != "" ]
then
    DOXYFILE=$1
else
    DOXYFILE=scip.dxy
fi

if [ "$HTML_FILE_EXTENSION" = "" ]
then
    HTML_FILE_EXTENSION=shtml
fi

### START SHELL TUTORIAL

# change into SCIP root directory so that we can run bin/scip instead of ../bin/scip
cd ..

# build a fresh version of SCIP
# make -j clean
make -j ZIMPL=false >> ${DOCLOG}

# run scip with some commands for the shell tutorial
bin/scip < doc/inc/shelltutorial/commands | tee doc/inc/shelltutorial/shelltutorialraw.tmp >> ${DOCLOG}

# cleanup of files created by the SCIP commands
rm stein27.lp stein27.sol settingsfile.set

cd doc

# modify the raw log file by adding doxygen snippet marker via the python script
python3 inc/shelltutorial/insertsnippetstutorial.py

### FINISHED SHELL TUTORIAL

### START FAQ GENERATION

cd inc/faq

python3 parser.py --linkext $HTML_FILE_EXTENSION  && php localfaq.php > faq.inc
cd ../../

### FINISHED FAQ GENERATION

### START PARAMETER FILE CREATION

cd ..

bin/scip -c "set default set save doc/inc/parameters.set quit"

bin/scip -c "read doc/inc/simpleinstance/simple.lp optimize quit" > doc/inc/simpleinstance/output.log

cd doc

### FINISHED FAQ GENERATION
# finally build the scip documentation
doxygen $DOXYFILE
