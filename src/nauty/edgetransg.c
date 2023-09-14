/* edgetransg.c  version 1.0; B D McKay, May 2017. */

#define USAGE "edgetransg [-t] [-q] [infile [outfile]]"

#define HELPTEXT \
" Select undirected graphs according to group action on vertices, edges and arcs.\n\
  Digraphs are not supported yet.\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -v  require vertex-transitive\n\
    -V  require not vertex-transitive\n\
    -e  require edge-transitive\n\
    -E  require not edge-transitive\n\
    -a  require arc-transitive\n\
    -A  require not arc-transitive\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 
#include "nautinv.h"
#include <limits.h>

/**************************************************************************/

typedef struct arc_st {int from,to;} arc;

DYNALLSTAT(arc,arclist,arclist_sz);
DYNALLSTAT(int,arcperm,arcperm_sz);
DYNALLSTAT(int,arcorbs,arcorbs_sz);
static int narcs,numarcorbs;

static int
findarc(arc *a, int na, int from, int to)
/* Find position of from->to.  (Must be present.) */
{
    int lo,mid,hi;

    lo = 0;
    hi = na - 1;

    while (lo <= hi)
    {
        mid = (lo + hi) / 2;
        if (a[mid].from == from)
        {
            if (a[mid].to == to)
                return mid;
            else if (a[mid].to > to)
                hi = mid - 1;
            else
                lo = mid + 1;
        }
        else if (a[mid].from > from)
            hi = mid - 1;
        else
            lo = mid + 1;
    }
    gt_abort(">E findarc error\n");
    return 0; /* never happens */
}

/**************************************************************************/

void
arcorbits(int ngens, int *p, int *orbs, int numorbs, int stab, int n)
{
    int i;

    if (ngens == 1) for (i = 0; i < narcs; ++i) arcorbs[i] = i;

    for (i = 0; i < narcs; ++i)
        arcperm[i] = findarc(arclist,narcs,p[arclist[i].from],p[arclist[i].to]);
    numarcorbs = orbjoin(arcorbs,arcperm,narcs);
}

static int
callnauty(sparsegraph *sg)
/* Call nauty and return the number of orbits.
   For non-trivial group, also fill in arc orbits */
{
    int n,i,j,k;
    DYNALLSTAT(int,lab,lab_sz);
    DYNALLSTAT(int,ptn,ptn_sz);
    DYNALLSTAT(int,orbits,orbits_sz);
    static DEFAULTOPTIONS_SPARSEGRAPH(options);
    statsblk stats;
    size_t *v;
    int *e,*d;

    n = sg->nv;

    DYNALLOC1(int,lab,lab_sz,n,"callnauty");
    DYNALLOC1(int,ptn,ptn_sz,n,"callnauty");
    DYNALLOC1(int,orbits,orbits_sz,n,"callnauty");

    if (sg->nde >= INT_MAX) gt_abort(">E too many arcs\n");
    narcs = (int)sg->nde;

    DYNALLOC1(arc,arclist,arclist_sz,narcs,"callnauty");
    DYNALLOC1(int,arcperm,arcperm_sz,narcs,"callnauty");
    DYNALLOC1(int,arcorbs,arcorbs_sz,narcs,"callnauty");

    sortlists_sg(sg);
 
    SG_VDE(sg,v,d,e);

    k = 0;
    for (i = 0; i < n; ++i)
    {
        for (j = v[i]; j < v[i]+d[i]; ++j)
        {
            arclist[k].from = i;
            arclist[k].to = e[j];
            ++k;
        }
    }
    if (k != narcs) gt_abort(">E narcs wrong\n");

    options.userautomproc = arcorbits;
    options.invarproc = distances_sg;
    options.mininvarlevel = 1;
    options.maxinvarlevel = 1;
    options.invararg = 2;

    sparsenauty(sg,lab,ptn,orbits,&options,&stats,NULL);

    return stats.numorbits;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet;
    boolean vswitch,Vswitch;
    boolean eswitch,Eswitch;
    boolean aswitch,Aswitch;
    boolean isvt,iset,isat;
    int i,j,norbs,argnum;
    int codetype,outcode;
    sparsegraph g,h;
    nauty_counter nin,nout;
    char *arg,sw;
    double t;

    HELP; PUTVERSION;

    SG_INIT(g);
    SG_INIT(h);

    infilename = outfilename = NULL;
    quiet = FALSE;
    vswitch = Vswitch = FALSE;
    eswitch = Eswitch = FALSE;
    aswitch = Aswitch = FALSE;

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
                     SWBOOLEAN('v',vswitch)
                else SWBOOLEAN('V',Vswitch)
                else SWBOOLEAN('e',eswitch)
                else SWBOOLEAN('E',Eswitch)
                else SWBOOLEAN('a',aswitch)
                else SWBOOLEAN('A',Aswitch)
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

#if 0
    if (!quiet)
    {
        fprintf(stderr,">A edgetransg");
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }
#endif

    if ((vswitch && Vswitch) || (eswitch && Eswitch) || (aswitch && Aswitch))
        gt_abort(">E edgetransg: contradictory switches\n");

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

    NODIGRAPHSYET(codetype);

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

    if (codetype&SPARSE6) outcode = SPARSE6;
    else                  outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
        else                    writeline(outfile,GRAPH6_HEADER);
    }

    nin = nout = 0;
    t = CPUTIME;
    while (read_sg(infile,&g))
    {
        ++nin;

        if (g.nde == 0 || g.nv == 1)
            isvt = iset = isat = TRUE;
        else
        {
            norbs = callnauty(&g);
            isvt = (norbs == 1);

            if (norbs == g.nv)
                iset = isat = FALSE;
            else
            {
                isat = (numarcorbs == 1);
                if (isat)
                    iset = TRUE;
                else
                {
                    for (i = 0; i < narcs; ++i)
                        arcperm[i] = findarc(arclist,narcs,arclist[i].to,arclist[i].from);
                    iset = orbjoin(arcorbs,arcperm,narcs) == 1;
                }
            }
        }

        if ( ((vswitch && isvt) || (Vswitch && !isvt) || (!vswitch && !Vswitch))
          && ((eswitch && iset) || (Eswitch && !iset) || (!eswitch && !Eswitch))
          && ((aswitch && isat) || (Aswitch && !isat) || (!aswitch && !Aswitch)) )
        {
            ++nout;
            writelast(outfile);
        }
    }
    t = CPUTIME - t;

    if (!quiet)
    {
        fprintf(stderr,">Z " COUNTER_FMT
                " graphs read from %s, " COUNTER_FMT 
                " written to %s in %3.2f sec.\n",
                nin,infilename,nout,outfilename,t);
    }

    exit(0);
}
