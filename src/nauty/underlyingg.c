/* underlyingg.c  version 1.0; B D McKay, Dec 2017. */

#define USAGE "underlyingg [-q] [infile [outfile]]"

#define HELPTEXT \
" Take the underlying undirected graphs of a file of graphs.\n\
\n\
    The output file has no header.\n\
    Undirected graphs are passed through without change, while\n\
    Underlying graphs of digraphs are written in sparse6 format.\n\
\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

/**************************************************************************/

static void
underlying(graph *g, int m, int n)
/* Replace g by its underlying graph */
{
    int i,j;
    graph *gi,*gj;

    for (i = 0, gi = g; i < n; ++i, gi += m)
        for (j = i+1, gj = gi+m; j < n; ++j, gj += m)
            if (ISELEMENT(gi,j) || ISELEMENT(gj,i))
            {
                ADDELEMENT(gi,j);
                ADDELEMENT(gj,i);
            }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet;
    boolean digraph;
    int j,m,n,argnum;
    int codetype,outcode;
    graph *g;
    nauty_counter nin;
    char *arg,sw;
    double t;

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    badargs = FALSE;
    quiet = FALSE;

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
        fprintf(stderr,">A underlyingg");
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

    gtools_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = 0;
    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;

        if (!digraph)
            writelast(outfile);
        else
        {
            underlying(g,m,n);
            writes6(outfile,g,m,n);
        }
        FREES(g);
    }
    t = CPUTIME - t;

    if (!quiet)
        fprintf(stderr,">Z  " COUNTER_FMT 
                " graphs converted from %s to %s in %3.2f sec.\n",
                nin,infilename,outfilename,t);

    exit(0);
}
