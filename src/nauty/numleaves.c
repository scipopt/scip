/* This is a sample function to show how copyg can be compiled
   as a filter.  This one defines the filter to be TRUE if
   the graph has between Qlo and Qhi leaves.  Outdegree=1 counts
   as a leaf in the case of digraphs.

   To compile this into a program "numleaves", use this:
     gcc -o numleaves -O3 -DFILTER=numleaves copyg.c numleaves.c nauty.a
   or similarly with a different compiler.

   At execution time, the interval [Qlo,Qhi] is taken from the -Q switch.
   A missing lower bound is indicated by Qlo = -NOLIMIT and a missing
   upper bound by Qhi = NOLIMIT.
*/

#include "gtools.h"

boolean
numleaves(graph *g, boolean digraph, long Qlo, long Qhi, int m, int n)
{
    int i,leaves;
    set *gi;

  /* If m=1, it would be more efficient to use the POPCOUNT
     macro in this funtion without calling setsize().  But
     we won't do that as it is just an example. */

    leaves = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (setsize(gi,m) == 1) ++leaves;

    return leaves >= Qlo && (Qhi == NOLIMIT || leaves <= Qhi);
}
