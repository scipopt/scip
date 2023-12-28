/* subdivideg.c  version 1.1; B D McKay, Oct 2017. */

#define USAGE "subdivideg [-k#] [-q] [infile [outfile]]"

#define HELPTEXT \
" Make the subdivision graphs of a file of graphs.\n\
    -k#  Subdivide each edge by # new vertices (default 1)\n\
\n\
    The output file has a header if and only if the input file does.\n\
\n\
    -q  Suppress auxiliary information.\n"

/*************************************************************************/

#include "gtools.h" 

/**************************************************************************/

static void
subdivisiongraph(sparsegraph *g, int k, sparsegraph *h)
/* h := subdivision graph of g, k new vertices per edge ;
   version for undirected graphs. */
{
    DYNALLSTAT(size_t,eno,eno_sz);   /* edge number */
    int *ge,*gd,*he,*hd;
    size_t *gv,*hv;
    int gnv,hnv;
    size_t i,j,l,gnde,hnde,num;
    size_t hi,lo,mid,w;

    sortlists_sg(g);
    if (k == 0)
    {
        copy_sg(g,h);
        return;
    }

    SG_VDE(g,gv,gd,ge);
    gnv = g->nv;
    gnde = g->nde;
    DYNALLOC1(size_t,eno,eno_sz,gnde,"subdivideg");

    hnv = gnv + k*(gnde/2);
    if ((hnv > 0 && hnv <= 0)
             || (gnde > 0 && ((size_t)(hnv-gnv))/(gnde/2) != k))
        gt_abort(">E subdivideg: output graph too large\n");
    hnde = gnde * (k+1);
    if (hnde/(k+1) != gnde)
        gt_abort(">E subdivideg: output graph too large\n");

    num = 0;
    for (i = 0; i < gnv; ++i)
    {
        for (j = gv[i]; j < gv[i]+gd[i]; ++j)
        {
            if (ge[j] == i)
                gt_abort(">E subdivideg can't handle undirected loops\n");
            else if (ge[j] > i)
                eno[j] = num++;
            else
            {
                lo = gv[ge[j]];
                hi = lo + gd[ge[j]] - 1;
                while (lo <= hi)
                {
                    mid = lo + (hi-lo)/2;
                    if (ge[mid] == i) break;
                    else if  (ge[mid] < i) lo = mid+1;
                    else hi = mid-1;
                }
                if (lo > hi)
                    gt_abort(">E subdivideg : binary search failed\n");
                eno[j] = eno[mid];
            }
        }
    }

    SG_ALLOC(*h,hnv,hnde,"subdivideg");
    h->nv = hnv;
    h->nde = hnde;
    SG_VDE(h,hv,hd,he);

    for (i = 0; i < gnv; ++i)
    {
        hd[i] = gd[i];
        hv[i] = gv[i];
    }
    for (i = gnv; i < hnv; ++i)
    {
        hd[i] = 2;
        hv[i] = gnde + 2*(i-gnv);
    }

    for (i = 0; i < gnv; ++i)
    {
        for (j = gv[i]; j < gv[i]+gd[i]; ++j)
            if (ge[j] > i)
            {
                w = gnv + k*eno[j];
                he[j] = w;
                he[hv[w]] = i;
                for (l = 1; l < k; ++l)
                {
                    he[hv[w]+1] = w+1;
                    he[hv[w+1]] = w;
                    ++w;
                }
            }
            else
            {
                w = gnv + k*eno[j] + k - 1;
                he[j] = w;
                he[hv[w]+1] = i;
            }
    }
}

/**************************************************************************/

static void
subdivisiondigraph(sparsegraph *g, int k, sparsegraph *h)
/* h := subdivision graph of g, k new vertices per edge ;
   version for directed graphs. */
{
    int *ge,*gd,*he,*hd;
    size_t *gv,*hv,hvend;
    int gnv,hnv;
    size_t i,j,l,t,v,gnde,hnde;
    size_t hi,lo,mid,w;

    sortlists_sg(g);
    if (k == 0)
    {
        copy_sg(g,h);
        return;
    }

    SG_VDE(g,gv,gd,ge);
    gnv = g->nv;
    gnde = g->nde;

    if (gnv + k*gnde > NAUTY_INFINITY-2)
        gt_abort(">E result would be too large\n");
    hnv = gnv + k*gnde;
    hnde = gnde + k*gnde;

    SG_ALLOC(*h,hnv,hnde,"subdivisiong");
    h->nv = hnv;
    h->nde = hnde;
    SG_VDE(h,hv,hd,he);

    for (i = 0; i < gnv; ++i)
    {
        hv[i] = gv[i];
        hd[i] = gd[i];
    }
    hvend = hv[gnv-1] + hd[gnv-1];
    v = gnv;

    for (i = 0; i < gnv; ++i)
    {
        for (j = gv[i]; j < gv[i]+gd[i]; ++j)
        {
            l = ge[j];
            he[j] = v;
            for (t = 0; t < k; ++t)
            {
                hd[v+t] = 1;
                hv[v+t] = hvend+t;
                he[hv[v+t]] = v+t+1;
            }
            he[hv[v+k-1]] = l;
            v += t;
            hvend += t;
        }
    }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs,quiet,kswitch,digraph;
    int j,argnum,kvalue,nloops;
    int codetype,outcode;
    SG_DECL(g); SG_DECL(h);
    nauty_counter nin;
    char *arg,sw;
    double t;

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    quiet = kswitch = FALSE;
    kvalue = 1;

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
                else SWINT('k',kswitch,kvalue,"subdivideg -k")
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
        fprintf(stderr,">A subdivideg");
        if (kswitch) fprintf(stderr," -k%d",kvalue);
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

    if (codetype&SPARSE6) outcode = SPARSE6;
    else                  outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6) writeline(outfile,SPARSE6_HEADER);
        else                    writeline(outfile,GRAPH6_HEADER);
    }

    nin = 0;
    t = CPUTIME;
    while (read_sgg_loops(infile,&g,&nloops,&digraph))
    {
        ++nin;

        if (digraph) subdivisiondigraph(&g,kvalue,&h);
        else         subdivisiongraph(&g,kvalue,&h);

        if (digraph)                 writed6_sg(outfile,&h);
        else if (outcode == SPARSE6) writes6_sg(outfile,&h);
        else                         writeg6_sg(outfile,&h);
    }
    t = CPUTIME - t;

    if (!quiet)
    {
        fprintf(stderr,">Z " COUNTER_FMT
                       " graphs converted from %s to %s in %3.2f sec.\n",
                       nin,infilename,outfilename,t);
    }

    exit(0);
}
