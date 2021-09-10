#!/usr/bin/env bash
#
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            *
#*                            fuer Informationstechnik Berlin                *
#*                                                                           *
#*  SCIP is distributed under the terms of the ZIB Academic License.         *
#*                                                                           *
#*  You should have received a copy of the ZIB Academic License              *
#*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: glue.sh,v 1.1 2010/10/26 03:52:21 bzfviger Exp $
#
# This script glues the instances listed in short.test into a single .cip file.
# First, each instance is translated into CIP format using generic variables names.
# Next, the variables and constraints are extracted, the names are prefixed with a unique identifier,
# and the objective coefficients are transformed for a minimization problem.
#
# Indicator instances are skipped, since the indicator-reader reports an error when parsing instances with generic names.
# Semicontinuous instances are skipped, since bounddisjunction constraints cannot be parsed yet.
#
# TODO: Compute the optimal value of the glued problem from the objective values specified in short.solu.
# TODO: norm objective functions from each instance to, e.g., 1.0

rm -f glue_vars.cip glue_conss.cip

j=0

for i in `cat ../check/testset/short.test` ;
do

# skip instances that cannot be parsed from cip yet or are infeasible
case $i in
  instances/Indicator/* ) continue ;;
  instances/MINLP* ) continue ;;
  instances/Pseudo* ) continue ;;
  instances/PseudoBoolean/normalized-t2001.13queen13.1111218308.opb ) continue ;;
  instances/MIP/stein27_inf.lp ) continue;
esac

name=`basename $i`
(( j++ ))

echo "$j: $name"

../bin/scip -q -c "r ../check/$i w genprob glue_tmp.cip q"

awk -v prefix=i$j '
BEGIN {
    inconss = 0;
  }
/  Sense/ {
    if ($3 == "minimize") sense = 1.0 ; else sense = -1.0;
  }
/CONSTRAINTS/ {
    inconss = 1; next;
  }
/END/ {
    exit;
  }
  { if (NR <= 7) next; }
  { gsub("<x", "<" prefix "_x");
    gsub("<~x", "<~" prefix "_x");
    gsub("<c", "<" prefix "_c");
    if( inconss ) {
      print >> "glue_conss.cip" ;
    }
    else {
      if (sense == -1.0)
        if (gsub("obj=-", "obj=") == 0)
         gsub("obj=", "obj=-")
      print >> "glue_vars.cip";
    }
  }
' glue_tmp.cip >> glue_vars.cip

done

echo "STATISTICS" > glue.cip
echo "OBJECTIVE" >> glue.cip
echo "  Sense : minimize" >> glue.cip
echo "VARIABLES" >> glue.cip
cat glue_vars.cip >> glue.cip
echo "CONSTRAINTS" >> glue.cip
cat glue_conss.cip >> glue.cip
echo "END" >> glue.cip

rm -f glue_{tmp,vars,conss}.cip

echo "presolving/components/maxintvars = 10000" > glue.set
echo "presolving/components/nodelimit = 10000000" >> glue.set

echo "Finished. Created file glue.cip and SCIP settings file glue.set."
