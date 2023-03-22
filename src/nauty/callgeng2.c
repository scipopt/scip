/* This is a sample of how to call geng as a procedure in multiple threads.
 * Usage: callgeng2 [-N#] geng_args
 *   The argument of -N is the number of threads (default 10, maximum 100).
 * It must be linked with a version of geng and nauty that has been compiled
 * with the thread-local attribute.
 * The way to compile this program per the Posix standard and recent gcc compilers
 * is like this:
  gcc -o callgeng2 -O3 -march=native -pthread -DMAXN=WORDSIZE -DWORDSIZE=32 \
        -DUSE_TLS -DGENG_MAIN=geng_main -DSUMMARY=geng_summary  \
        callgeng2.c geng.c nautyW1.a
 * You might also need to add -lpthread or -lpthreads to the compilation.
 * On some systems these instructions will need modification.
 *
 * If output is not turned off using -u, graphs are written to stdout. Each graph
 * is written using a single fwrite() call, which is atomic on Linux for such
 * short strings but there is no guarantee that the output won't get mixed together
 * on other systems. Since this is only a short demonstration, I won't try to
 * fix that, but if you know thread programming you can figure out a way to have
 * a single thread do all the writing (use the hook OUTPROC, see callgeng.c).
 * Meanwhile, if you just want to count and not write, you don't need to worry
 * about that.
 */

#include "gtools.h"
#include <pthread.h>

#define MAXTHREADS 100
#define DEFAULTTHREADS 10
#define MAXGENGARGS 20

static nauty_counter count[MAXTHREADS];  /* Thread i only writes to count[i]. */
static TLS_ATTR int threadnumber;
typedef struct params { int res,mod; char **args; } params;

void SUMMARY(nauty_counter num, double cpu)
/* This will be called from each geng thread as it finishes. */
{
   count[threadnumber] = num;
}

/* Declare the main procedure of geng (usually main()) */
int GENG_MAIN(int argc, char *argv[]);

static void*
runit(void *threadarg)          /* Main routine for one thread */
{
    params *par;
    int geng_argc;
    char resmod[30]; /* string to put res/mod into */
    char *geng_argv[MAXGENGARGS+2]; /* At least one bigger than geng_argc. */
    char **args;

    par = (params*)(threadarg);
    threadnumber = par->res;
    args = par->args;

  /* Set up geng argument list.  The 0-th argument is the command name.
   * There must be a NULL at the end. */

    sprintf(resmod,"%d/%d",par->res,par->mod);

    geng_argv[0] = "callgeng2";
    geng_argv[1] = "-q";    /* stops geng writing stuff to stderr */
    geng_argc = 2;
    while (*args != NULL) geng_argv[geng_argc++] = *(args++);
    geng_argv[geng_argc++] = resmod;
    geng_argv[geng_argc] = NULL;

    GENG_MAIN(geng_argc,geng_argv);

    return NULL;
} 

int
main(int argc, char *argv[])
{
    int n,i,ret;
    char **gengargs;
    pthread_t thread[MAXTHREADS];
    params par[MAXTHREADS];
    nauty_counter totalcount;
    int threads;
    double t0,t1;
    boolean badargs;

    badargs = FALSE;

    if (argc < 2 || argc > MAXGENGARGS-2) badargs = TRUE;

  /* gengargs will be where the geng arguments start */

    if (!badargs)
    {
        if (sscanf(argv[1],"-N%d",&threads) == 1)
            gengargs = &argv[2];
        else
        {
            threads = DEFAULTTHREADS;
            gengargs = &argv[1];
        }
    }

    if (badargs || threads < 1 || threads > MAXTHREADS)
    {
        fprintf(stderr,"Usage: %s [-N#] geng_args...\n",argv[0]);
        exit(0);
    }

    nauty_check(WORDSIZE,1,MAXN,NAUTYVERSIONID);

    t0 = CPUTIME;

    for (i = 0; i < threads; ++i)
    {
        par[i].res = i;
        par[i].mod = threads;
        par[i].args = gengargs;
        if ((ret = pthread_create(&thread[i],NULL,runit,&par[i])) != 0)
        {
            fprintf(stderr,">E Thread creation failed, code=%d\n",ret);
            exit(1);
        }
    }

    for (i = 0; i < threads; ++i)
    {
        if ((ret = pthread_join(thread[i],NULL)) != 0)
        {
            fprintf(stderr,">E Thread joining failed, code=%d\n",ret);
            exit(1);
        }
    }

    totalcount = 0;

    for (i = 0; i < threads; ++i) totalcount += count[i];

    t1 = CPUTIME;

    fprintf(stderr,">Z " COUNTER_FMT " graphs made in %.2f seconds.\n",
                   totalcount,t1-t0);

    exit(0);
}
