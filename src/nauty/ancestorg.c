/* ancestorg.c  version 1.0; B D McKay, October 2014. */

#define USAGE "ancestorg [-q] [-g#:#|-g#] [infile [outfile]]"

#define HELPTEXT \
" The g-th generation ancestor of a graph is the graph obtained by removing\n\
  the final g vertices. The 0-th generation ancestor is the graph itself.\n\
  For each input graph, write the ancestors whose generation is given by the\n\
  g argument.  No zero-sized graphs are written.\n\
  Output is always in graph6 format.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -g# -g#:#  Specify a generation or range of generations (default: all)\n\
    -q  Suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h" 
static nauty_counter nout;

/**************************************************************************/

static void
writegenerations(FILE *f, graph *g, int m, int n, long mingen, long maxgen,
                 boolean digraph, int outcode)
/* write all the specified generations (mingen and maxgen are already sensible */
/* g is destroyed */
{
    int i,gen;
    set *gi;

    for (gen = 0; gen <= maxgen; ++gen)
    {
        if (gen >= mingen) 
        {
            if (digraph)                 writed6(f,g,m,n);
            else if (outcode == SPARSE6) writes6(f,g,m,n);
            else                         writeg6(f,g,m,n);
            ++nout;
        }
        for (i = 0, gi = g; i < n-1; ++i, gi += m) DELELEMENT(gi,n-1);
        --n;
    }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet,gswitch,digraph;
    int i,j,m,n,v,argnum;
    int n1,m1;
    int codetype,outcode;
    graph *g;
    nauty_counter nin;
    char *arg,sw;
    setword *gv;
    long mingen,maxgen;
    double t;

    HELP;

    infilename = outfilename = NULL;
    badargs = FALSE;
    gswitch = quiet = FALSE;

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
                    SWBOOLEAN('q',quiet)
                else SWRANGE('g',":-",gswitch,mingen,maxgen,">E ancestorg -g")
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

    if (!gswitch)
    {
        mingen = 0;
        maxgen = NAUTY_INFINITY;
    }
    if (mingen < 0) mingen = 0;

    if (!quiet)
    {
        fprintf(stderr,">A ancestorg");
        if (gswitch) fprintf(stderr," -g%ld:%ld",mingen,maxgen);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

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

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = nout = 0;
    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;

        writegenerations(outfile,g,m,n,mingen,(maxgen >= n ? n-1 : maxgen),
                                                              digraph,outcode);
        FREES(g);
    }
    t = CPUTIME - t;

    if (!quiet)
        fprintf(stderr,
             ">Z  " COUNTER_FMT " graphs read from %s, "
                    COUNTER_FMT " written to %s; %3.2f sec.\n",
             nin,infilename,nout,outfilename,t);

    exit(0);
}
