/* converseg.c  version 3.0; B D McKay, Dec 2020. */

#define USAGE "converseg [-q] [-a|-c] [infile [outfile]]"

#define HELPTEXT \
" Take the converse digraphs of a file of directed graphs.\n\
\n\
    The output file has a header if and only if the input file does.\n\
    Undirected graphs are passed through without change, while\n\
    directed graphs are written in digraph6 format.\n\
\n\
    -a  Also output the original graph (before the converse)\n\
    -c  Output only self-converse digraphs\n\
\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 
#include "nautinv.h"

/**************************************************************************/

static void
conv(graph *g, int m, int n)
/* Replace g by its converse */
{
    int i,j;
    graph *gi,*gj;

    for (i = 0, gi = g; i < n; ++i, gi += m)
        for (j = i+1, gj = gi+m; j < n; ++j, gj += m)
            if ((ISELEMENT(gi,j)!=0) + (ISELEMENT(gj,i)!=0) == 1)
            {
                FLIPELEMENT(gi,j);
                FLIPELEMENT(gj,i);
            }
}

/**************************************************************************/

static boolean
isselfconverse(graph *g, int m, int n)
/* Test if self-converse */
{
    int i,j,deg;
    size_t ij;
    DYNALLSTAT(graph,gc,gc_sz);
    DYNALLSTAT(graph,h,h_sz);
    DYNALLSTAT(graph,hc,hc_sz);
    nauty_counter gcode,gccode;
    graph *gi;
    static const int fuzz[] = {037541,017543,061532,005257,026416};
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_DIGRAPH(options);
    statsblk stats;

    DYNALLOC2(graph,gc,gc_sz,m,n,">E converseg");
    for (ij = 0; ij < n*(size_t)m; ++ij) gc[ij] = g[ij];
    conv(gc,m,n);

    gcode = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
        deg = 0;
        for (j = 0; j < m; ++j) deg += POPCOUNT(gi[j]);
        gcode += (deg ^ fuzz[deg%5]);
    }

    gccode = 0;
    for (i = 0, gi = gc; i < n; ++i, gi += m)
    {
        deg = 0;
        for (j = 0; j < m; ++j) deg += POPCOUNT(gi[j]);
        gccode += (deg ^ fuzz[deg%5]);
    }
    
    if (gcode != gccode) return FALSE;

    DYNALLOC2(graph,h,h_sz,m,n,">E converseg");
    DYNALLOC2(graph,hc,hc_sz,m,n,">E converseg");
    DYNALLOC1(int,lab,lab_sz,n,">E converseg");
    DYNALLOC1(int,ptn,ptn_sz,n,">E converseg");
    DYNALLOC1(int,orbits,orbits_sz,n,">E converseg");

    options.getcanon = TRUE;
    densenauty(g,lab,ptn,orbits,&options,&stats,m,n,h);
    densenauty(gc,lab,ptn,orbits,&options,&stats,m,n,hc);

    for (ij = 0; ij < n*(size_t)m; ++ij)
        if (h[ij] != hc[ij]) return FALSE;

    return TRUE;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet;
    boolean digraph,also,self;
    int j,m,n,argnum;
    int codetype,outcode;
    graph *g;
    nauty_counter nin,nout;
    char *arg,sw;
    double t;

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    badargs = FALSE;
    quiet = also = self = FALSE;

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
                     SWBOOLEAN('a',also)
                else SWBOOLEAN('c',self)
                else SWBOOLEAN('q',quiet)
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

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!quiet)
    {
        fprintf(stderr,">A converseg");
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        if (also || self)
        {
            fprintf(stderr,"-");
            if (also) fprintf(stderr,"a");
            if (self) fprintf(stderr,"c");
        }
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (also && self)
        gt_abort(">E converseg: -a and -c are incompatible\n");

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
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

    if      (codetype&SPARSE6)  outcode = SPARSE6;
    else if (codetype&DIGRAPH6) outcode = DIGRAPH6;
    else                        outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6)       writeline(outfile,SPARSE6_HEADER);
        else if (outcode == DIGRAPH6) writeline(outfile,DIGRAPH6_HEADER);
        else                          writeline(outfile,GRAPH6_HEADER);
    }

    gtools_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = nout = 0;
    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;

        if (!digraph)
        {
            writelast(outfile);
            ++nout;
        }
        else if (self)
        {
            if (isselfconverse(g,m,n))
            {
                writelast(outfile);
                ++nout;
            }
        }
        else
        {
            if (also)
            {
                writed6(outfile,g,m,n);
                ++nout;
            }
            conv(g,m,n);
            writed6(outfile,g,m,n);
            ++nout;
        }
        FREES(g);
    }
    t = CPUTIME - t;

    if (!quiet)
    {
        if (self)
            fprintf(stderr,">Z  " COUNTER_FMT 
                " digraphs read from %s; " COUNTER_FMT
                " self-converse written to %s in %3.2f sec.\n",
                nin,infilename,nout,outfilename,t);
        else
            fprintf(stderr,">Z  " COUNTER_FMT 
                " graphs converted from %s to %s in %3.2f sec.\n",
                nin,infilename,outfilename,t);
    }

    exit(0);
}
