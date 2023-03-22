/* copyg.c version 2.5; B D McKay, June 2020 */

#ifdef FILTER
#define HELPUSECMD
#define USAGE "[-v] [-gszfp#:#qhx] [infile [outfile]]"
#else
#define USAGE "copyg [-gszfp#:#qhx] [infile [outfile]]"
#endif

#define HELPTEXT \
"  Copy a file of graphs with possible format conversion.\n\
\n\
     -g  Use graph6 format for output\n\
     -s  Use sparse6 format for output\n\
     -z  Use digraph6 format for output\n\
     -i  Use incremental sparse6 format for output\n\
     In the absence of -g, -s, -z or -i, the format\n\
     depends on the header or, if none, the first input line.\n\
     As an exception, digraphs are always written in digraph6.\n\
\n\
     -p# -p#:#  \n\
        Specify range of input lines (first is 1)\n\
        This can fail if the input has incremental lines.\n\
     -f  With -p, assume input lines of fixed length\n\
        (ignored if header or first line has sparse6 format).\n\
     -I#  Have at most this number of incremental steps\n\
        in a row.  Implies -i. \n\
\n\
     -h  Write a header.\n\
     -x  Don't write a header.\n\
     In the absence of -h and -x, a header is written if\n\
     there is one in the input.\n\
\n\
     -q  Suppress auxiliary output.\n"

#ifdef FILTER
#define HELPTEXT1 HELPTEXT
#define HELPTEXT2 \
"\n\
   -Q# -Q#:# a number or range to pass to the filter.\n\
\n\
   Only inputs satisfying the filter, or not satisfying the\n\
   filter if -v is given, are passed through.\n"
#endif

/***********************************************************************/

/* How to use this as a skeleton to construct a filter.

   Suppose we want to make a filter "numleaves" that filters graphs
   according to their number of leaves.

   (1) In a separate file numleaves.c (present in the package), write
   a procedure with prototype matching
      boolean
      numleaves(graph *g, boolean digraph, long Qlo, long Qhi, int m, int n);
   It should return a nonzero value if g has between Qlo and Qhi leaves,
   and return 0 otherwise. The values of Qlo and Qhi are obtained from
   the -Q switch, with missing lower bound -NOLIMIT and missing upper
   bound NOLIMIT.

   (2) Compile numleaves.c together with this program with the preprocessor
   variable FILTER defined as "numleaves" (without the quotes). Link with
   nauty.a. For the gcc, icc, and some other compilers, you can do it like
   this:
     gcc -o numleaves -O3 -DFILTER=numleaves copyg.c numleaves.c nauty.a

   (3) That's it.  Now "numleaves -Q2:7" will copy through graphs with 2-7
   leaves, while "numleaves -v -Q7" will copy through graphs whose number
   of leaves is not 7.

   Within the procedure numleaves(), the index number of the current graph
   is available in the external nauty_counter variable nin.  The first
   graph has index 1.

   If you also define SUMMARY to be a procedure with prototype
   void name(void), it is called before the final summary (unless
   -q is given).  You can't use SUMMARY without FILTER.
*/

/***********************************************************************/

#include "gtools.h"

#ifdef FILTER
extern boolean
FILTER(graph *g, boolean digraph, long Qlo, long Qhi, int m, int n);
#ifdef SUMMARY
extern void SUMMARY(void);
#endif
#endif

nauty_counter nin;

int
main(int argc, char *argv[])
{
    graph *g,*gprev,*gbase;
    int m,n,codetype;
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    int outcode;
    nauty_counter nout;
    int argnum,j,nprev,mprev;
    char *arg,sw;
    boolean sswitch,fswitch,pswitch,qswitch,gswitch;
    boolean hswitch,xswitch,iswitch,Iswitch,zswitch;
    boolean badargs,digraph,vswitch,Qswitch,kept;
    long Qlo,Qhi;
    long pval1,pval2,maxin,refresh,inclines;
    double t0,t1;

    HELP; PUTVERSION;

    iswitch = Iswitch = sswitch = fswitch = pswitch = FALSE;
    gswitch = qswitch = xswitch = hswitch = zswitch = FALSE;
    vswitch = Qswitch = FALSE;
    Qlo = Qhi = 0;  /* avoid some warnings */
    infilename = outfilename = NULL;

    argnum = 0;
    badargs = FALSE;
    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('g',gswitch)
                else SWBOOLEAN('z',zswitch)
                else SWBOOLEAN('q',qswitch)
                else SWBOOLEAN('f',fswitch)
                else SWBOOLEAN('h',hswitch)
                else SWBOOLEAN('i',iswitch)
                else SWBOOLEAN('x',xswitch)
#ifdef FILTER
                else SWBOOLEAN('v',vswitch)
                else SWRANGE('Q',":-",Qswitch,Qlo,Qhi,"filter -Q")
#endif
                else SWLONG('I',Iswitch,refresh,"copyg -I")
                else SWRANGE('p',":-",pswitch,pval1,pval2,"copyg -p")
                else badargs = TRUE;
            }
        }
        else
        {
            ++argnum;
            if      (argnum == 1) infilename = arg;
            else if (argnum == 2) outfilename = arg;
            else                  badargs = TRUE;
        }
    }

    if (Iswitch) iswitch = TRUE;
    if ((sswitch!=0) + (gswitch!=0) + (iswitch!=0) + (zswitch!=0) > 1) 
        gt_abort(">E copyg: -s, -g, -z and -i/-I are incompatible\n");
    if (hswitch && xswitch) 
        gt_abort(">E copyg: -h and -x are incompatible\n");

    if (badargs || argnum > 2)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!qswitch)
    {
        fprintf(stderr,">A %s",argv[0]);
        if (sswitch || gswitch || iswitch || Iswitch || fswitch
            || zswitch || pswitch || xswitch || hswitch || vswitch || Qswitch)
            fprintf(stderr," -");
        if (sswitch) fprintf(stderr,"s");
        if (gswitch) fprintf(stderr,"g");
        if (zswitch) fprintf(stderr,"z");
        if (hswitch) fprintf(stderr,"h");
        if (Iswitch) fprintf(stderr,"I%ld",refresh); 
        else if (iswitch) fprintf(stderr,"i");
        if (xswitch) fprintf(stderr,"x");
        if (fswitch) fprintf(stderr,"f");
        if (pswitch) writerange(stderr,'p',pval1,pval2);
        if (vswitch) fprintf(stderr,"v");
        if (Qswitch) writerange(stderr,'Q',Qlo,Qhi);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,fswitch,
                                        pswitch ? pval1 : 1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

    if (!outfilename || outfilename[0] == '-')
    {
        outfilename = "stdout";
        outfile = stdout;
    }
    else if ((outfile = fopen(outfilename,"w")) == NULL)
    {
        fprintf(stderr,"Can't open output file %s\n",outfilename);
        gt_abort(NULL);
    }

    if (gswitch)                  outcode = GRAPH6;
    else if (sswitch || iswitch)  outcode = SPARSE6;
    else if (zswitch)             outcode = DIGRAPH6;
    else if ((codetype&GRAPH6))   outcode = GRAPH6;
    else if ((codetype&SPARSE6))  outcode = SPARSE6;
    else if ((codetype&DIGRAPH6)) outcode = DIGRAPH6;
    else
    {
        outcode = GRAPH6;
        fprintf(stderr,
             ">W copyg doesn't handle this graph format, writing graph6.\n");
    }

    if (hswitch || (!xswitch && (codetype&HAS_HEADER)))
    {
        if (outcode == SPARSE6)       writeline(outfile,SPARSE6_HEADER);
        else if (outcode == DIGRAPH6) writeline(outfile,DIGRAPH6_HEADER);
        else                          writeline(outfile,GRAPH6_HEADER);
    }

    nin = nout = 0;
    if (!pswitch || pval2 == NOLIMIT) maxin = NOLIMIT;
    else if (pval1 < 1) maxin = pval2;
    else                maxin = pval2 - pval1 + 1;

#ifdef FILTER
    if (!Qswitch) { Qlo = -NOLIMIT; Qhi = NOLIMIT; }
#endif

    gprev = NULL;
    nprev = mprev = 1;
    inclines = 0;
    t0 = CPUTIME;
    while (nin < maxin || maxin == NOLIMIT)
    {
        if ((g = readgg_inc(infile,NULL,0,&m,&n,gprev,mprev,nprev,&digraph))
                                                      == NULL) break;
        ++nin;
        kept = TRUE;
#ifdef FILTER
        if ((!vswitch && !FILTER(g,digraph,Qlo,Qhi,m,n)) ||
                           (vswitch && FILTER(g,digraph,Qlo,Qhi,m,n))) 
            kept = FALSE;
        else
            ++nout;
#endif

        if (!kept)
        { }
        else if (iswitch)  
        {
            if (digraph) gt_abort(
                ">Z incremental sparse6 is incompatible with digraphs\n");
            gbase = gprev;
            if (nprev != n || mprev != m) gbase = NULL;
            if (Iswitch > 0 && inclines == refresh) gbase = NULL;
            if (gbase == NULL) inclines = 0; else ++inclines;
            writeis6(outfile,g,gbase,m,n);
        }
        else if (outcode == readg_code)   writelast(outfile);
        else if (digraph)             writed6(outfile,g,m,n);
        else if (outcode == SPARSE6)  writes6(outfile,g,m,n);
        else if (outcode == DIGRAPH6) writed6(outfile,g,m,n);
        else                          writeg6(outfile,g,m,n);

        if (gprev) FREES(gprev);
        gprev = g;
        nprev = n;
        mprev = m;
    }
    t1 = CPUTIME;

    if (!qswitch)
#ifdef FILTER
#ifdef SUMMARY
        SUMMARY();
#endif
        fprintf(stderr,">Z  " COUNTER_FMT " graphs read from %s, "
                              COUNTER_FMT " written to %s; %.2f sec\n",
                       nin,infilename,nout,outfilename,t1-t0);
#else
        fprintf(stderr,">Z  " COUNTER_FMT " graphs copied from %s to %s\n",
                       nin,infilename,outfilename);
#endif

    exit(0);
}
