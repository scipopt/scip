#!/bin/gawk -f
#
BEGIN {
   lastprob = "";
}
#
# problem name
#
/^@01/ { 
   probfile = $2;
}
/^=ready=/ {
   lastprob = probfile;
}
END {
   printf(lastprob);
}
