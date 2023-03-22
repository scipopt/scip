/* assembleg.c  version 1.1; B D McKay, Sep 2018. */

#define USAGE "assembleg -n#|-n#:# [-i#|i#:#] [-L] [-q] [infile [outfile]]"

#define HELPTEXT \
" Assemble input graphs as components of output graphs.\n\
\n\
    The output file has no header.\n\
    If the input has any directed graphs, all outputs are directed.\n\
    Otherwise, the output format is determined by the header\n\
       or first input.\n\
    The input graphs had better all fit into memory at once,\n\
       unless -L is given, in which case only the graphs of at\n\
       most half the output size are stored at once.\n\
    The output graphs will be non-isomorphic if the input\n\
       graphs are connected and non-isomorphic.\n\
\n\
    -n# -n#:#  Give range of output sizes (compulsory)\n\
    -i# -i#:#  Give range of input sizes to use\n\
    -L  Assume all input graphs strictly larger than maxn/2\n\
         vertices follow any smaller graphs in the input,\n\
         where maxn is the largest size specified by -n.\n\
         This can greatly reduce memory consumption.\n\
    -c  Also write graphs consisting of a single input\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

static int ninputs;             /* Number of inputs */
static nauty_counter nout;      /* Number of outputs */
typedef graph *graphptr;
static graphptr *gin;           /* Their contents */
static int *size;               /* Their sizes */
static int outcode;

/**************************************************************************/

static void
insertg(graph *g, int ng, graph *h, int nh, int n)
/* Insert nh-vertex graph starting at vertex ng in graph g.
   n is the total size available.
   It is assumed that g is empty apart from 0..ng-1. */
{
    int i,j,m,mh;
    set *gi,*hi;

    m = SETWORDSNEEDED(n);
    mh = SETWORDSNEEDED(nh);

    for (i = 0, hi = h, gi = g + ng*m; i < nh;
                            ++i, hi += mh, gi += m)
    {
        for (j = -1; (j = nextelement(hi,mh,j)) >= 0; )
            ADDELEMENT(gi,ng+j);
    }
}

static void
removeg(graph *g, int ng, int nh, int n)
/* Remove a subgraph that was in position ng..ng+nh-1. */
{
    set *gi;
    int i,j,m;

    m = SETWORDSNEEDED(n);

    for (i = ng, gi = g + ng*m; i < ng+nh; ++i, gi += m)
    {
        for (j = ng; j < ng+nh; ++j)
            DELELEMENT(gi,j);
    }
}

/**************************************************************************/

#define SORT_NAME sortbysize
#define SORT_OF_SORT 2
#define SORT_TYPE1 int
#define SORT_TYPE2 graphptr
#include "sorttemplates.c"

static void
readinputs(FILE *f, int imin, int imax)
/* Read inputs and sort by size */
{
    size_t tablesize;
    graph *g;
    int m,n;
    boolean digraph;

    if ((gin = malloc(sizeof(graphptr)*10000)) == NULL ||
        (size = malloc(sizeof(int)*10000)) == NULL)
            gt_abort(">E malloc failed in readinputs()\n");
    tablesize = 10000;
    
    ninputs = 0;

    while (TRUE)
    {
        if ((g = readgg(f,NULL,0,&m,&n,&digraph)) == NULL) break;
        if (digraph) outcode = DIGRAPH6;
        if (n < imin || n > imax) continue;

        if (ninputs == tablesize)
        {
            tablesize = 3*tablesize/2 + 10000;
            if ((gin = realloc(gin,sizeof(graphptr)*tablesize)) == NULL ||
                (size = realloc(size,sizeof(int)*tablesize)) == NULL)
                    gt_abort(">E realloc failed in readinputs()\n");
        }

        gin[ninputs] = g;
        size[ninputs] = n;
        ++ninputs;
    }

    if (ninputs < 0 || ninputs > tablesize)
        gt_abort(">E Some overflow problem in readinputs()\n");

    sortbysize(size,gin,ninputs);
}

static void
readsomeinputs(FILE *f, int imin, int imax, int maxsize,
                                                  graph **g, int *n)
/* Read inputs and sort by size.  Stop when either EOF or a
   graph bigger than maxsize is read. */
{
    size_t tablesize;
    int m;
    boolean digraph;

    if ((gin = malloc(sizeof(graphptr)*10000)) == NULL ||
        (size = malloc(sizeof(int)*10000)) == NULL)
            gt_abort(">E malloc failed in readsomeinputs()\n");
    tablesize = 10000;

    ninputs = 0;

    while (TRUE)
    {
        if ((*g = readgg(f,NULL,0,&m,n,&digraph)) == NULL) break;
        if (digraph) outcode = DIGRAPH6;
        if (*n > maxsize) break;

        if (*n < imin || *n > imax) continue;

        if (ninputs == tablesize)
        {
            tablesize += 10000;
            if ((gin = realloc(gin,sizeof(graphptr)*tablesize)) == NULL ||
                (size = realloc(size,sizeof(int)*tablesize)) == NULL)
                    gt_abort(">E realloc failed in readsomeinputs()\n");
        }

        gin[ninputs] = *g;
        size[ninputs] = *n;
        ++ninputs;
    }

    if (ninputs < 0 || ninputs > tablesize)
        gt_abort(">E Some overflow problem in readinputs()\n");

    sortbysize(size,gin,ninputs);
}

/**************************************************************************/

static void
assemble(graph *g, int nmin, int nmax, int sofar, int lastpos,
          boolean writeconn, FILE *outfile)
/* Recursively add one more graph.
   sofar = how many vertices so far. */
{
    int pos,newsize;

    for (pos = lastpos; pos < ninputs; ++pos)
    {
        newsize = sofar + size[pos];
        if (newsize > nmax) break;

        insertg(g,sofar,gin[pos],size[pos],nmax);
        if (newsize >= nmin && (sofar > 0 || writeconn))
        {
            if (outcode == DIGRAPH6)
                writed6(outfile,g,SETWORDSNEEDED(nmax),newsize);
            else if (outcode == GRAPH6)
                writeg6(outfile,g,SETWORDSNEEDED(nmax),newsize);
            else
                writes6(outfile,g,SETWORDSNEEDED(nmax),newsize);
            ++nout;
        }
        assemble(g,nmin,nmax,newsize,pos,writeconn,outfile);
        removeg(g,sofar,size[pos],nmax);
    }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet;
    boolean nswitch,cswitch,iswitch,Lswitch;
    boolean digraph;
    int j,mm,n,argnum;
    int codetype;
    graph *gread,*gout;
    char *arg,sw;
    double t;
    long nmin,nmax,mmax,imin,imax;

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    badargs = FALSE;
    iswitch = nswitch = cswitch = quiet = Lswitch = FALSE;

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
                else SWBOOLEAN('c',cswitch)
                else SWBOOLEAN('L',Lswitch)
                else SWRANGE('n',":-",nswitch,nmin,nmax,"assembleg -n")
                else SWRANGE('i',":-",iswitch,imin,imax,"assembleg -i")
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

    if (!nswitch) gt_abort(">E assembleg: -n is compulsory\n");

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (nmin <= 0) nmin = 1;
    if (nmin > NAUTY_INFINITY-2) nmin = NAUTY_INFINITY-2;

    if (!quiet)
    {
        fprintf(stderr,">A assembleg -");
        if (nmin == nmax) fprintf(stderr,"n%ld",nmin);
        else fprintf(stderr,"n%ld:%ld",nmin,nmax);
        if (iswitch)
        {
            if (imin == imax) fprintf(stderr,"i%ld",imin);
            else fprintf(stderr,"i%ld:%ld",imin,imax);
        }
        if (cswitch) fprintf(stderr,"c");
        if (Lswitch) fprintf(stderr,"L");
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (!iswitch || imin <= 0) imin = 1;
    if (!iswitch || imax > nmax) imax = nmax;
    if (!cswitch && imax == nmax) --imax;

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

    gtools_check(WORDSIZE,1,1,NAUTYVERSIONID);

    t = CPUTIME;

    mmax = SETWORDSNEEDED(nmax);

    if ((gout = malloc(mmax*sizeof(graph)*nmax)) == NULL)
        gt_abort(">E assembleg: malloc() failed in main()\n");

    nout = 0;

    if (Lswitch)
    {
        readsomeinputs(infile,(int)imin,(int)imax,(int)nmax/2,&gread,&n);

        EMPTYSET(gout,mmax*(size_t)nmax);
        assemble(gout,(int)nmin,(int)nmax,0,0,cswitch,outfile);

        while (gread)
        {
            if (n >= imin && n <= imax)
            { 
                EMPTYSET(gout,mmax*(size_t)nmax);
                insertg(gout,0,gread,n,(int)nmax);
    
                if (n >= nmin && n <= nmax && cswitch)
                {
                    if (outcode == DIGRAPH6)
                        writed6(outfile,gout,mmax,n);
                    else if (outcode == GRAPH6)
                        writeg6(outfile,gout,mmax,n);
                    else
                        writes6(outfile,gout,mmax,n);
                    ++nout;
                }

                assemble(gout,(int)nmin,(int)nmax,n,0,cswitch,outfile);
            }
            FREES(gread);

            if ((gread = readgg(infile,NULL,0,&mm,&n,&digraph)) == NULL) break;
            if (digraph) outcode = DIGRAPH6;
            if (n <= nmax/2) gt_abort(">E assembleg -L : inputs in bad order\n");
        }
    }
    else
    {
        readinputs(infile,(int)imin,(int)imax);

        EMPTYSET(gout,mmax*(size_t)nmax);
        assemble(gout,(int)nmin,(int)nmax,0,0,cswitch,outfile);
    }
 
    t = CPUTIME - t;

    if (!quiet)
        fprintf(stderr,">Z %d graphs read from %s; " COUNTER_FMT 
                " graphs written to %s in %3.2f sec.\n",
                ninputs,infilename,nout,outfilename,t);

    exit(0);
}
