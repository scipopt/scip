/* This is a sample of how to call geng as a procedure rather than
 * running it as a separate process.  The basic idea is to construct
 * an argument list for geng's main() function.  At compile time,
 * assign a name to the macros OUTPROC and GENG_MAIN.  A typical
 * Unix-style compilation command would be:
     gcc -o callgeng -O3 -DMAXN=WORDSIZE -DOUTPROC=myoutproc -DGENG_MAIN=geng_main \
       callgeng.c geng.c nauty.a
 * You will get best performance if you add machine-specific optimization arguments
 * such as -march=native.
 */

#include "gtools.h"

static unsigned long counter;

void
OUTPROC(FILE *outfile, graph *g, int n)
{
 /* This will be called for each graph. */

    ++counter;
}

int
GENG_MAIN(int argc, char *argv[]);

int
main(int argc, char *argv[])
{
    int n,geng_argc;
    char nstring[6]; /* string to put n into */
    char *geng_argv[6]; /* At least one bigger than geng_argc. */

  /* Set up geng argument list.  The 0-th argument is the command name.
   * There must be a NULL at the end.  This example is for connected
   * bipartite graphs of maximum degree at most 4 for 3-10 vertices. */

    for (n = 3; n <= 10; ++n)
    {
        sprintf(nstring,"%d",n);
        geng_argv[0] = "geng";
        geng_argv[1] = "-q";    /* stops geng writing stuff to stderr */
        geng_argv[2] = "-cb";
        geng_argv[3] = "-D4";
        geng_argv[4] = nstring;
        geng_argv[5] = NULL;
        geng_argc = 5;

        counter = 0;
        GENG_MAIN(geng_argc,geng_argv);
        printf("Number of graphs with %d vertices = %lu.\n",n,counter);
    }

    exit(0);
}
