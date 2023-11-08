/* delptg.c  version 3.3; B D McKay, Jan 2022. */

#define USAGE \
 "delptg [-lq] [-a|-b] [-d#|-d#:#] [-v#|-v#:#] [-r#] [-n#] [-m#|-i] [infile [outfile]]"

#define HELPTEXT \
" Delete some vertices from a file of graphs.\n\
\n\
    The output file has a header if and only if the input file does.\n\
    No isomorph reduction is done.\n\
\n\
    -l  Canonically label outputs\n\
    -d# -d#:# Only remove vertices with original degree in the given range\n\
        For digraphs, the out-degree is used.\n\
    -n# The number of vertices to delete (default 1).\n\
    -v# -v#:# Vertex number or numbers that it is allowed to delete\n\
              (first vertex is number 0).\n\
    -m# Lower bound on minimum degree of output graphs.\n\
    -r# Choose # random sets of points (not necessarily different)\n\
    -a  The deleted points must be adjacent.\n\
    -b  The deleted points must be non-adjacent.\n\
    -i  Leave deleted vertices as isolates, not compatible with -m.\n\
    No empty graphs are output. No warning is issued if\n\
           -d, -v -n, -m together imply no graphs are output.\n\
    For digraphs, out-degree is used for -d and -m.\n\
    -q  Suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h" 

static FILE *outfile;
static nauty_counter nout;
static int outcode;
static boolean digraph,dolabel;

/**************************************************************************/

static void
writeone(graph *g, int m, int n, int *del, int ndel,
                             int outmindeg, boolean isolates)
/* Delete the stated vertices and write it is mindeg is high enough.
 * The vertices to delete are del[0] < del[1] < ... < del[ndel-1]. */
{
    int i,j,k,nx,mx,deg;
    graph *gi,*gxi,*gq;
    ssize_t ii;
#if MAXN
    graph gx[MAXN*MAXM];
    graph hx[MAXN*MAXM];
    int lab[MAXN];
#else
    DYNALLSTAT(graph,gx,gx_sz);
    DYNALLSTAT(graph,hx,hx_sz);
    DYNALLSTAT(int,lab,lab_sz);
#endif

    if (isolates)
    {
        nx = n;
        mx = m;
    }
    else
    {
        nx = n - ndel;    /* Always positive */
        mx = SETWORDSNEEDED(nx);
    }

#if !MAXN
    DYNALLOC2(graph,gx,gx_sz,mx,nx,"delptg");
    if (dolabel) DYNALLOC2(graph,hx,hx_sz,mx,nx,"delptg");
    DYNALLOC1(int,lab,lab_sz,nx,"delptg");
#endif

    if (isolates)
    {
        for (ii = mx*(ssize_t)nx; --ii >= 0; )
            gx[ii] = g[ii];
        for (j = 0; j < ndel; ++j)
        {
            k = del[j];
            EMPTYSET(GRAPHROW(gx,k,mx),mx);
            for (i = 0, gxi = gx; i < nx; ++i, gxi += mx)
                DELELEMENT(gxi,k);
        }
    }
    else
    {
        j = 0;
        for (i = 0; i < del[0]; ++i) lab[j++] = i;
        for (k = 1; k < ndel; ++k)
            for (i = del[k-1]+1; i < del[k]; ++i) lab[j++] = i;
        for (i = del[ndel-1]+1; i < n; ++i) lab[j++] = i;
    
        EMPTYSET(gx,nx*(size_t)mx);
    
        for (i = 0, gxi = (set*)gx; i < nx; ++i, gxi += mx)
        {
            gi = GRAPHROW(g,lab[i],m);
            for (j = 0; j < nx; ++j)
            {
                k = lab[j];
                if (ISELEMENT(gi,k)) ADDELEMENT(gxi,j);
            }
        }
    }

    if (outmindeg > 0)
    {
        for (i = 0, gxi = (set*)gx; i < nx; ++i, gxi += mx)
        {
            deg = 0;
            for (j = 0; j < mx; ++j) deg += POPCOUNT(gxi[j]);
            if (deg < outmindeg) return;
        }
    }

    if (dolabel)
    {
        fcanonise(gx,mx,nx,hx,NULL,digraph);
        gq = hx;
    }
    else
        gq = gx;

    if (outcode == DIGRAPH6 || digraph) writed6(outfile,gq,mx,nx);
    else if (outcode == SPARSE6)        writes6(outfile,gq,mx,nx);
    else                                writeg6(outfile,gq,mx,nx);
    ++nout;
}

/**************************************************************************/

static void
delrandom(int ndel, int *del, graph *g, int m, int n, int outmindeg,
                int *okverts, int nvok, boolean isolates)
/* Delete ndel vertices chosen at random from okverts[0..nvok-1] */
{
    int i,k;

    for (i = 0; i < nvok; ++i) del[i] = okverts[i];

    if (2*ndel < nvok)
    {                   /* Choose ndel */
        k = 0;
        do
        {
            i = KRAN(nvok);
            if (del[i] >= 0)
            {
                del[i] = -1 - del[i];
                ++k;
            }
        } while (k < ndel);
        k = 0;
        for (i = 0; i < nvok; ++i)
            if (del[i] < 0) del[k++] = -1 - del[i];
    }
    else
    {       /* Remove nvok-ndel */
        k = 0;
        do
        {
            i = KRAN(nvok);
            if (del[i] >= 0)
            {
                del[i] = -1;
                ++k;
            }
        } while (k < nvok-ndel);
        k = 0;
        for (i = 0; i < nvok; ++i)
            if (del[i] >= 0) del[k++] = del[i];
    }

    writeone(g,m,n,del,ndel,outmindeg,isolates);
}

/**************************************************************************/

static void
search(int level, int ndel, int *del, graph *g, int m, int n, int outmindeg,
        int *okverts, int lastok, int nvok, boolean isolates)
/* level vertices have been chosen, last one was okverts[lastok] */
{
    int i;

    if (level == ndel)
    {
        writeone(g,m,n,del,ndel,outmindeg,isolates);
        return;
    }

    for (i = lastok+1; i <= level + nvok - ndel; ++i) 
    {
        del[level] = okverts[i];
        search(level+1,ndel,del,g,m,n,
                  outmindeg,okverts,i,nvok,isolates);
    }
}

/**************************************************************************/

static void
search_a(int level, int ndel, int *del, graph *g, int m, int n, set *avail,
         int outmindeg, int *okverts, int lastok, int nvok, boolean isolates)
/* level vertices have been chosen, last one was okverts[lastok] */
{
    int i,j,v;
    set *gv;

    if (level == ndel)
    {
        writeone(g,m,n,del,ndel,outmindeg,isolates);
        return;
    }

    for (i = lastok+1; i <= level + nvok - ndel; ++i) 
    {
        v = okverts[i];
        if (ISELEMENT(avail,v))
        {
            gv = g + m*(size_t)v;
            if (level < ndel-1)
                for (j = 0; j < m; ++j) avail[m+j] = avail[j] & gv[j];
            del[level] = v;
            search_a(level+1,ndel,del,g,m,n,
                       avail+m,outmindeg,okverts,i,nvok,isolates);
        }
    }
}

/**************************************************************************/

static void
search_b(int level, int ndel, int *del, graph *g, int m, int n, set *avail,
         int outmindeg, int *okverts, int lastok, int nvok, boolean isolates)
/* level vertices have been chosen, last one was okverts[lastok] */
{
    int i,j,v;
    set *gv;

    if (level == ndel)
    {
        writeone(g,m,n,del,ndel,outmindeg,isolates);
        return;
    }

    for (i = lastok+1; i <= level + nvok - ndel; ++i) 
    {
        v = okverts[i];
        if (ISELEMENT(avail,v))
        {
            gv = g + m*(size_t)v;
            if (level < ndel-1)
                for (j = 0; j < m; ++j) avail[m+j] = avail[j] & ~gv[j];
            del[level] = v;
            search_b(level+1,ndel,del,g,m,n,
                      avail+m,outmindeg,okverts,i,nvok,isolates);
        }
    }
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile;
    boolean badargs,quiet,dswitch,nswitch,vswitch,mswitch;
    boolean adj,nonadj,isolates,delrand;
    int i,j,m,n,v,argnum;
    int ndel,outmindeg;
    int codetype;
    graph *g;
    nauty_counter nin;
    char *arg,sw;
    setword *gv;
    long mindeg,maxdeg,irand,numrand;
    long minv,maxv;
    int iminv,imaxv,actmaxv,nvok;
    int degv;
    double t;
#if MAXN
    int okverts[MAXN];
    int del[MAXN];
#else
    DYNALLSTAT(int,okverts,okverts_sz);
    DYNALLSTAT(int,del,del_sz);
#endif
    DYNALLSTAT(set,avail,avail_sz);

    HELP; PUTVERSION;

    infilename = outfilename = NULL;
    badargs = vswitch = mswitch = FALSE;
    dswitch = nswitch = quiet = FALSE;
    adj = nonadj = isolates = delrand = FALSE;

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
                     SWBOOLEAN('l',dolabel)
                else SWBOOLEAN('q',quiet)
                else SWBOOLEAN('a',adj)
                else SWBOOLEAN('b',nonadj)
                else SWBOOLEAN('i',isolates)
                else SWINT('n',nswitch,ndel,">E delptg -n")
                else SWRANGE('d',":-",dswitch,mindeg,maxdeg,">E delptg -d")
                else SWRANGE('v',":-",vswitch,minv,maxv,">E delptg -v")
                else SWINT('m',mswitch,outmindeg,">E delptg -m")
                else SWLONG('r',delrand,numrand,">E delptg -r")
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

    if (nswitch && ndel < 1) gt_abort(">E delptg: bad argument for -n\n");
    if (adj && nonadj) gt_abort(">E delptg: -a and -b are incompatible\n");
    if (mswitch && isolates)
        gt_abort(">E delptg: -m and -i are incompatible\n");
    if (delrand && (adj || nonadj))
        gt_abort(">E delptg: -r is incompatible with -a and -b\n");

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!quiet)
    {
        fprintf(stderr,">A delptg");
        if (dolabel||adj||nonadj||isolates)
            fprintf(stderr," -%s%s%s%s",(dolabel?"l":""),
                (adj?"a":""),(nonadj?"b":""),(isolates?"i":""));
        if (dswitch) fprintf(stderr," -d%ld:%ld",mindeg,maxdeg);
        if (nswitch) fprintf(stderr," -n%d",ndel);
        if (vswitch) fprintf(stderr," -v%ld:%ld",minv,maxv);
        if (mswitch) fprintf(stderr," -m%d",outmindeg);
        if (delrand) fprintf(stderr," -r%ld",numrand);
        if (argnum > 0) fprintf(stderr," %s",infilename);
        if (argnum > 1) fprintf(stderr," %s",outfilename);
        fprintf(stderr,"\n");
        fflush(stderr);
    }

    if (vswitch)
    {
        iminv = (minv < 0 ? 0 : (int) minv);
        imaxv = (maxv >= NAUTY_INFINITY ? NAUTY_INFINITY-1 : maxv);
    }
    else
    {
        iminv = 0;
        imaxv = NAUTY_INFINITY-1;
    }

    if (!mswitch) outmindeg = 0;

    if (!dswitch)
    {
        mindeg = 0;
        maxdeg = NAUTY_INFINITY;
    }

    if (!nswitch) ndel = 1;
    if (ndel == 1) adj = nonadj = FALSE;

    if (dolabel) nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    nin = nout = 0;
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

    if (codetype&SPARSE6)       outcode = SPARSE6;
    else if (codetype&DIGRAPH6) outcode = DIGRAPH6;
    else                        outcode = GRAPH6;

    if (codetype&HAS_HEADER)
    {
        if (outcode == SPARSE6)       writeline(outfile,SPARSE6_HEADER);
        else if (outcode == DIGRAPH6) writeline(outfile,DIGRAPH6_HEADER);
        else                          writeline(outfile,GRAPH6_HEADER);
    }

    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++nin;

        if (n <= ndel)
        {
            FREES(g);
            continue;
        }

        actmaxv = (imaxv < n ? imaxv : n-1);
        if (actmaxv - iminv + 1 < ndel)
        {
            FREES(g);
            continue;
        }

#if !MAXN
        DYNALLOC1(int,okverts,okverts_sz,n,"delptg");
        DYNALLOC1(int,del,del_sz,n,"delptg");
#endif
        if (adj || nonadj)
        {
            DYNALLOC2(set,avail,avail_sz,m,ndel-1,"delptg");
            EMPTYSET(avail,m);
        }

        nvok = 0;
        for (v = 0, gv = g; v < n && v <= actmaxv; ++v, gv += m)
        {
            if (v < iminv) continue;
            degv = 0;
            for (i = 0; i < m; ++i) degv += POPCOUNT(gv[i]);
            if (degv < mindeg || degv > maxdeg) continue;
            okverts[nvok++] = v;
        }

        if (nvok < ndel)
        {
            FREES(g);
            continue;
        }

        if (adj || nonadj)
        {
            for (i = 0; i < nvok; ++i)
                ADDELEMENT(avail,okverts[i]);
        }

        if (delrand) INITRANBYTIME;

        if (delrand)
        {
            for (irand = 0; irand < numrand; ++irand)
                delrandom(ndel,del,g,m,n,outmindeg,okverts,nvok,isolates);
        }
        else if (adj)
           search_a(0,ndel,del,g,m,n,avail,outmindeg,okverts,-1,nvok,isolates);
        else if (nonadj)
           search_b(0,ndel,del,g,m,n,avail,outmindeg,okverts,-1,nvok,isolates);
        else
           search(0,ndel,del,g,m,n,outmindeg,okverts,-1,nvok,isolates);

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
