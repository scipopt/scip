#!/usr/bin/env bash
# 
# This scripts generates all the parameter settings which can be create via
# the interactive shell of SCIP.
#

SETTINGS=(aggressive fast off)
PLUGINS=(heuristics presolving separating)
EMPHASES=(counter cpsolver easycip feasibility hardlp optimality)

if test ! -e settings
then
    mkdir settings
fi

for PLUGIN in ${PLUGINS[@]}
do
    for SETTING in ${SETTINGS[@]}
    do
	SETTINGNAME=settings/$PLUGIN\_$SETTING.set
	bin/scip -c "set default set $PLUGIN emphasis $SETTING set diff $SETTINGNAME quit"
    done
done

for EMPHASIS in ${EMPHASES[@]}
do
    SETTINGNAME=settings/$EMPHASIS.set
    bin/scip -c "set default set emphasis $EMPHASIS set save $SETTINGNAME quit"
done
