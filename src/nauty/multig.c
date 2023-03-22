/* multig.c version 2.1; B D McKay, March 4, 2022 */

#define USAGE \
"multig [-q] [-V] [-u|-T|-G|-A|-B] [-e#|-e#:#] \n" \
"       [-m#] [-f#] [-D#|-r#|-l#] [infile [outfile]]"

#define HELPTEXT \
" Read undirected loop-free graphs and replace their edges with multiple\n\
  edges in all possible ways (multiplicity at least 1).\n\
  Isomorphic multigraphs derived from the same input are suppressed.\n\
  If the input graphs are non-isomorphic then the output graphs are also.\n\
\n\
    -e# | -e#:#  specify a value or range of the total number of edges\n\
                 counting multiplicities\n\
    -m# maximum edge multiplicity (minimum is 1)\n\
    -D# upper bound on maximum degree\n\
    -r# make regular of specified degree (incompatible with -l, -D, -e)\n\
    -l# make regular multigraphs with multiloops, degree #\n\
                 (incompatible with -r, -D, -e)\n\
    -V  read the T format as produced by vcolg and obey the vertex colours\n\
        in computing the automorphism group. If -T or -G is used as the\n\
        output format, a list of the input colours is included.\n\
    Either -l, -r, -D, -e or -m with a finite maximum must be given\n\
    -f# Use the group that fixes the first # vertices setwise\n\
    -T  use a simple text output format (nv ne {v1 v2 mult})\n\
    -G  like -T but includes group size as third item (if less than 10^10)\n\
          The group size does not include exchange of isolated vertices.\n\
    -A  write as the upper triangle of an adjacency matrix, row by row,\n\
        including the diagonal, and preceded by the number of vertices\n\
    -B  write as an integer matrix preceded by the number of rows and\n\
        number of columns, where -f determines the number of rows\n\
    -u  no output, just count them\n\
    -q  suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h"
#include "naugroup.h"

nauty_counter mg_nin,mg_nout,mg_skipped;
FILE *outfile;

#define MAXNV 128 
#define MAXNE 1024

static int v0[MAXNE+MAXNV],v1[MAXNE+MAXNV];
static int edgeno[MAXNV][MAXNV];
static int lastlev[MAXNE];
static int colour[MAXNV];

static int ix[MAXNE+MAXNV],nix;
static boolean first;
static int lastreject[MAXNV];
static boolean lastrejok;
static unsigned long groupsize;
static unsigned long newgroupsize;
static boolean Vswitch,Gswitch,Tswitch,Aswitch,Bswitch;
static int Brows;

#define GROUPTEST_NOT 
#ifdef GROUPTEST
static long long totallab;
#endif

#define PATHCOUNTS_NOT
#ifdef PATHCOUNTS
static long long count0,count1,count2,count3,count4,count5;
static nauty_counter oldlo;
#endif

/* INPUTGRAPH/INPUTGRAPHC feature
 *
 * If INPUTGRAPH is defined, it must expand as the name of a procedure
 * with prototype like  int INPUTGRAPH(graph *g, int m, int n).
 * This procedure will be called for each input graph before it is
 * processed. The graph will be skipped if the value 0 is returned.
 *
 * If INPUTGRAPHC is defined, it must expand as the name of a procedure
 * with prototype like  int INPUTGRAPH(graph *g, int *col, int m, int n).
 * This procedure will be called for each input graph before it is
 * processed. col[i] is the colour of vertex i as read using -V,
 * otherwise it is 0.
 * The graph will be skipped if the value 0 is returned.
 *
 * At most one of INPUTGRAPH and INPUTGRAPHC can be defined.
 */

/* OUTPROC/OUTPROCC feature
 *
 * If OUTPROC is defined at compile time, and -u is not used, the
 * procedure OUTPROC is called for each graph.  This must be linked
 * by the compiler.  The arguments are
 * f = open output file (FILE*)
 * n = number of vertices (int)
 * ne = number of edges (int)
 * gp = unsigned long group size ignoring isolated vertices
 *                                       (note: may have overflowed)
 * v0[*], v1[*], ix[*] = int arrays.  The edges are
 * v0[i]-v1[i] with multiplicity ix[i] for i=0..ne-1.   ix[i]>0 always.
 *
 *
 * If OUTPROCC is defined at compile time, and -u is not used, the
 * procedure OUTPROCC is called for each graph.  This must be linked
 * by the compiler.  The arguments are
 * f = open output file (FILE*)
 * n = number of vertices (int)
 * ne = number of edges (int)
 * gp = unsigned long group size ignoring isolated vertices
 *                                      (note: may have overflowed)
 * col[*], v0[*], v1[*], ix[*] = int arrays. The colour of vertex
 * i is col[i] (this is the value read in with -V and is 0 in the
 * absence of -V. The edges are v0[i]-v1[i] with multiplicity ix[i]
 * for i=0..ne-1.   ix[i]>0 always.
 * 
 * At most one of OUTPROC and OUTPROCC can be defined. */

/* SUMMARY feature
 *
 * If SUMMARY is defined, it must expand as the name of a procedure
 * with prototype  void SUMMARY(void).  It is called at the end before
 * the normal summary (which can be suppressed with -q).  The numbers of
 * graphs read and multigraphs produced are available in the global variables
 * mg_nin and mg_nout (type nauty_counter).
 */

#if defined(INPUTGRAPH) && defined(INPUTGRAPHC)
#error -- At most one of INPUTGRAPH and INPUTGRAPHC can be defined.
#endif

#ifdef INPUTGRAPH
extern int INPUTGRAPH(graph*,int,int);
#endif

#ifdef INPUTGRAPHC
extern int INPUTGRAPHC(graph*,int*,int,int);
#endif

#if defined(OUTPROC) && defined(OUTPROCC)
#error -- At most one of OUTPROC and OUTPROCC can be defined.
#endif

#ifdef OUTPROC
extern void OUTPROC(FILE*,int,int,unsigned long,int*,int*,int*);
#endif

#ifdef OUTPROCC
extern void OUTPROCC(FILE*,int,int,unsigned long,int*,int*,int*,int*);
#endif

#ifdef SUMMARY
extern void SUMMARY(void);
#endif

/**************************************************************************/

#if 0
static void
writeautom(int *p, int n)
/* Called by allgroup as a test. */
{
    int i;

    for (i = 0; i < n; ++i) printf(" %2d",p[i]);
    printf("\n");
}
#endif

/**************************************************************************/

static boolean
ismax(int *p, int n)
/* test if x^p <= x */
{
    int i,k;

    for (i = 0; i < nix; ++i)
    {
        k = edgeno[p[v1[i]]][p[v0[i]]];

        if (ix[k] > ix[i]) return FALSE;
        else if (ix[k] < ix[i]) return TRUE;
    }

    ++newgroupsize;
    return TRUE;
}

/**************************************************************************/

static void
testmax(int *p, int n, int *abort)
/* Called by allgroup2. */
{
    int i;

    if (first)
    {                       /* only the identity */
        first = FALSE;
        return;
    }

    if (!ismax(p,n))
    {
        *abort = 1;
        for (i = 0; i < n; ++i) lastreject[i] = p[i];
        lastrejok = TRUE;
    }
}

/**************************************************************************/

static void
printam(FILE *f, int n, int ne, int *ix)
/* Write adjacency matrix formats */
{
    int i,j;

    if (Aswitch)
    {
        fprintf(f,"%d ",n);
        for (i = 0; i < n; ++i)
            for (j = i; j < n; ++j)
                fprintf(f," %d",(edgeno[i][j]>=0 ? ix[edgeno[i][j]]: 0));
        fprintf(f,"\n");
    }
    else
    {
        if (Brows <= 0 || Brows > n)
        {
            fprintf(stderr,">E multig: impossible matrix size for output\n");
            exit(1);
        }
        fprintf(f,"%d %d",Brows,n-Brows);

        for (i = 0; i < Brows; ++i)
        {
            fprintf(f," ");
            for (j = Brows; j < n; ++j)
                fprintf(f," %d",(edgeno[i][j]>=0 ? ix[edgeno[i][j]]: 0));
        }
        fprintf(f,"\n");
    }
}

/**************************************************************************/

static void
trythisone(grouprec *group,
       boolean lswitch, int *deg, int maxdeg, int ne, int n)
/* Try one solution, accept if minimal. */
{
    int i,j,ne2;
    boolean accept;
    char line[20*(MAXNV+MAXNE)+50],s[20],*p;

    nix = ne;
    newgroupsize = 1;

    if (!group || groupsize == 1)
        accept = TRUE;
    else if (lastrejok && !ismax(lastreject,n))
        accept = FALSE;
    else if (lastrejok && groupsize == 2)
        accept = TRUE;
    else
    {
        newgroupsize = 1;
        first = TRUE;

        if (allgroup2(group,testmax) == 0)
            accept = TRUE;
        else
            accept = FALSE;
    }

    if (accept)
    {
#ifdef GROUPTEST
        if (groupsize % newgroupsize != 0)
                    gt_abort("group size error\n");
        totallab += groupsize/newgroupsize;
#endif

        ++mg_nout;

        if (outfile)
        {
            ne2 = ne;
            if (lswitch)
                for (i = 0; i < n; ++i)
                    if (deg[i] < maxdeg)
                    {
                        v0[ne2] = v1[ne2] = i;
                        edgeno[i][i] = ne2;
                        ix[ne2] = (maxdeg-deg[i])/2;
                        ++ne2;
                    }
#ifdef OUTPROC
            OUTPROC(outfile,n,ne2,newgroupsize,v0,v1,ix);
#else
#ifdef OUTPROCC
            OUTPROCC(outfile,n,ne2,newgroupsize,colour,v0,v1,ix);
#else
            if (Aswitch || Bswitch)
                printam(outfile,n,ne2,ix);
            else
            {
#define PUTINT(xx) { unsigned long ul = (xx); char *sp; \
 if (ul == 0) *(p++) = '0'; \
 else { sp = s; while (ul) { *(sp++) = (ul % 10) + '0'; ul /= 10; } \
        while (sp > s) { *(p++) = *(--sp); } }}
#define SPC *(p++) = ' '

                p = line;
                PUTINT(n); SPC; PUTINT(ne2);
                if (Gswitch) { SPC; PUTINT(newgroupsize); }
                if (Vswitch)
                    for (i = 0; i < n; ++i) { SPC; PUTINT(colour[i]); }
                SPC;
                for (i = 0; i < ne2; ++i)
                {  SPC; PUTINT(v0[i]);
                   SPC; PUTINT(v1[i]);
                   SPC; PUTINT(ix[i]);
                }
                *(p++) = '\n';
                *(p++) = '\0';
                fputs(line,outfile);
            }
#endif
#endif
            for (i = ne2; --i >= ne; )
            {
                j = v0[i];
                edgeno[j][j] = -1;
            }
        }
        return;
    }
    else
        return;
}

/**************************************************************************/

static void
scan(int level, int ne, long minedges, long maxedges, long sofar,
    long maxmult, grouprec *group, int n)
/* Recursive scan for default case */
{
    int left;
    long min,max,k;

    if (level == ne)
    {
        trythisone(group,FALSE,NULL,0,ne,n);
        return;
    }

    left = ne - level - 1;
    min = minedges - sofar - maxmult*left;
    if (min < 1) min = 1;
        max = maxedges - sofar - left;
    if (max > maxmult) max = maxmult;

    for (k = min; k <= max; ++k)
    {
        ix[level] = k;
        scan(level+1,ne,minedges,maxedges,sofar+k,maxmult,group,n);
    }

    return;
}

/**************************************************************************/

static void
scan_md(int level, int ne, long minedges, long maxedges, long sofar,
    long maxmult, grouprec *group, int n, int *deg, int maxdeg)
/* Recursive scan, maxdeg version */
{
    int left;
    int min,max,k;
    int x1,x2;

    if (level == ne)
    {
        trythisone(group,FALSE,deg,maxdeg,ne,n);
        return;
    }

    x1 = v0[level];
    x2 = v1[level];
    left = ne - level - 1;
    min = minedges - sofar - maxmult*left;
    if (min < 1) min = 1;
    max = maxedges - sofar - left;
    if (max > maxmult) max = maxmult;
    if (deg[x1] + max - 1 > maxdeg) max = maxdeg - deg[x1] + 1;
    if (deg[x2] + max - 1 > maxdeg) max = maxdeg - deg[x2] + 1;

    for (k = min; k <= max; ++k)
    {
        ix[level] = k;
        deg[x1] += k-1; deg[x2] += k-1;
        scan_md(level+1,ne,minedges,maxedges,sofar+k,maxmult,group,
                n,deg,maxdeg);
        deg[x1] -= k-1; deg[x2] -= k-1;
    }

    return;
}

/**************************************************************************/

static void
scan_lp(int level, int ne, long minedges, long maxedges, long sofar,
    long maxmult, grouprec *group, int n, int *deg, int maxdeg)
/* Recursive scan, regular-with-loops version. */
{
    int left;
    long min,max,k;
    int x1,x2;
    boolean odd,even;

    if (level == ne)
    {
        trythisone(group,TRUE,deg,maxdeg,ne,n);
        return;
    }

    x1 = v0[level];
    x2 = v1[level];
    left = ne - level - 1;
    min = minedges - sofar - maxmult*left;
    if (min < 1) min = 1;
    max = maxedges - sofar - left;
    if (max > maxmult) max = maxmult;
    if (deg[x1] + max - 1 > maxdeg) max = maxdeg - deg[x1] + 1;
    if (deg[x2] + max - 1 > maxdeg) max = maxdeg - deg[x2] + 1;

    odd = even = FALSE;
    if (lastlev[x1] == level)
    {
        if (((maxdeg-deg[x1])&1) == 1) even = TRUE; else odd = TRUE;
    }
    if (lastlev[x2] == level)
    {
        if (((maxdeg-deg[x2])&1) == 1) even = TRUE; else odd = TRUE;
    }
    if (even && odd) return;

    for (k = min; k <= max; ++k)
    {
        if (even && (k&1) == 1) continue;
        if (odd && (k&1) == 0) continue;

        ix[level] = k;
        deg[x1] += k-1; deg[x2] += k-1;
        scan_lp(level+1,ne,minedges,maxedges,sofar+k,maxmult,group,
                n,deg,maxdeg);
        deg[x1] -= k-1; deg[x2] -= k-1;
    }

    return;
}

/**************************************************************************/

static void
scan_reg(int level, int ne, long minedges, long maxedges, long sofar,
    long maxmult, grouprec *group, int n, int *delta, int *def, int maxdeg)
/* Recursive scan, regular version. */
{
    int left;
    long min,max,k;
    int x1,x2;

    if (level == ne)
    {
        trythisone(group,FALSE,NULL,maxdeg,ne,n);
        return;
    }

    x1 = v0[level];
    x2 = v1[level];
    left = ne - level - 1;
    min = minedges - sofar - maxmult*left;
    if (min < 1) min = 1;
    max = maxedges - sofar - left;
    if (max > maxmult) max = maxmult;
    if (max > def[x1] + 1) max = def[x1] + 1;
    if (max > def[x2] + 1) max = def[x2] + 1;

    if (min < def[x2] + 1 - delta[x1]) min = def[x2] + 1 - delta[x1];
    if (min < def[x1] + 1 - delta[x2]) min = def[x1] + 1 - delta[x2];

    if (lastlev[x1] == level && min < def[x1] + 1) min = def[x1] + 1;
    if (lastlev[x2] == level && min < def[x2] + 1) min = def[x2] + 1;

    for (k = min; k <= max; ++k)
    {
        ix[level] = k;
        delta[x1] += k-1 - def[x2];
        delta[x2] += k-1 - def[x1];
        def[x1] -= k-1;
        def[x2] -= k-1;
        scan_reg(level+1,ne,minedges,maxedges,sofar+k,maxmult,group,
                n,delta,def,maxdeg);
        def[x1] += k-1;
        def[x2] += k-1;
        delta[x1] -= k-1 - def[x2];
        delta[x2] -= k-1 - def[x1];
    }

    return;
}

/**************************************************************************/

static void
multi(graph *g, int nfixed, long minedges, long maxedges, long maxmult,
      int maxdeg, boolean lswitch, int m, int n)
{
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[2*MAXNV];
    grouprec *group;
    int ne;
    int i,j,k,j0,j1,thisdeg,maxd,x0,x1;
    set *gi;
    int lab[MAXNV],ptn[MAXNV],orbits[MAXNV],deg[MAXNV];
    int delta[MAXNV],def[MAXNV];
    set active[(MAXNV+WORDSIZE-1)/WORDSIZE];
    boolean isreg,hasloops;

#ifdef PATHCOUNTS
    ++count0;
#endif

    j0 = -1;  /* last vertex with degree 0 */
    j1 = n;   /* first vertex with degree > 0 */
 
    ne = 0;
    maxd = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
        thisdeg = 0;
        for (j = 0; j < m; ++j) thisdeg += POPCOUNT(gi[j]);
        deg[i] = thisdeg;
        if (thisdeg > maxd) maxd = thisdeg;
        if (thisdeg == 0) lab[++j0] = i;
        else              lab[--j1] = i;
        ne += thisdeg;
    }
    ne /= 2;

    if (maxdeg >= 0 && maxd > maxdeg) return;

#ifdef PATHCOUNTS
    ++count1;
#endif

    if (Aswitch || Bswitch)
        for (i = 0; i < n; ++i)
        for (j = 0; j < n; ++j)
            edgeno[i][j] = -1;

    if (ne == 0 && minedges <= 0
                && (!lswitch || (lswitch && (maxdeg&1) == 0)))
    {
        trythisone(NULL,lswitch,deg,maxdeg,0,n);
        return;
    }

#ifdef PATHCOUNTS
    ++count2;
#endif

    k = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
        for (j = i; (j = nextelement(gi,m,j)) >= 0; )
        {
            v0[k] = i;
            v1[k] = j;
            edgeno[i][j] = edgeno[j][i] = k;
            lastlev[i] = lastlev[j] = k;
            ++k;
        }
    }

    isreg = !lswitch && (maxdeg >= 0 && 2*minedges == n*(long)maxdeg);
        /* Case of regular multigraphs */

    if (isreg)  /* regular case */
    /* Condition: def(v) <= total def of neighbours */
    {
        for (i = 0; i < n; ++i)
        {
            def[i] = maxdeg - deg[i];
            delta[i] = -def[i];
        }

        for (i = 0; i < k; ++i)
        {
            x0 = v0[i]; x1 = v1[i];
            delta[x0] += def[x1];
            delta[x1] += def[x0];
        }

        for (i = 0; i < n; ++i) if (delta[i] < 0) return;
    }

    if ((isreg || lswitch) && (maxdeg & n & 1) == 1) return;
    if (isreg && j0 >= 0 && maxdeg > 0) return;
    if (lswitch && j0 >= 0 && (maxdeg&1) == 1) return;

#ifdef PATHCOUNTS
    ++count3;
#endif

    if (maxedges == NOLIMIT)
    {
        if (maxmult == NOLIMIT) maxedges = maxdeg*n/2;
        else                    maxedges = ne*maxmult;
    }
    if (maxmult == NOLIMIT) maxmult = maxedges - ne + 1;
    if (maxdeg >= 0 && maxmult > maxdeg) maxmult = maxdeg;
    if (maxedges < ne || ne*maxmult < minedges) return;

#ifdef PATHCOUNTS
    ++count4;
#endif

    if (n > MAXNV || ne > MAXNE)
        gt_abort(">E multig: MAXNV or MAXNE exceeded\n");

    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    for (i = 0; i < n; ++i) ptn[i] = 1;
    ptn[n-1] = 0;
    EMPTYSET(active,m);
    if (j0 != n-1) ADDELEMENT(active,j0+1);

    for (i = 0; i <= j0; ++i) ptn[i] = 0;

    for (i = j0+1; i < n; ++i)
    if (lab[i] < nfixed) break;

    if (i != j0+1 && i != n)
    {
        ptn[i-1] = 0;
        ADDELEMENT(active,i);
    }

    hasloops = FALSE;
    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (ISELEMENT(gi,i))
        {
            hasloops = TRUE;
            break;
        }

    options.defaultptn = FALSE;
    options.digraph = hasloops;
    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;

    nauty(g,lab,ptn,active,orbits,&options,&stats,workspace,2*MAXNV,m,n,NULL);

    if (stats.grpsize2 == 0)
        groupsize = stats.grpsize1 + 0.1;
    else
        groupsize = 0;

    group = groupptr(FALSE);
    makecosetreps(group);

    lastrejok = FALSE;

    if (isreg)
        scan_reg(0,ne,minedges,maxedges,0,maxmult,group,n,delta,def,maxdeg);
    else if (lswitch)
        scan_lp(0,ne,minedges,maxedges,0,maxmult,group,n,deg,maxdeg);
    else if (maxdeg >= 0)
        scan_md(0,ne,minedges,maxedges,0,maxmult,group,n,deg,maxdeg);
    else
        scan(0,ne,minedges,maxedges,0,maxmult,group,n);
}

/**************************************************************************/

static void
vmulti(graph *g, int nfixed, long minedges, long maxedges, long maxmult,
      int maxdeg, boolean lswitch, int ne, int m, int n)
/* In this version, v0[0..ne-1], v1[0..ne-1], edgeno[*][*]
    and colour[0..ne-1] are already filled in. */
{
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[2*MAXNV];
    grouprec *group;
    int i,j,k,thisdeg,maxd,x0,x1;
    set *gi;
    int lab[MAXNV],ptn[MAXNV],orbits[MAXNV],deg[MAXNV];
    int delta[MAXNV],def[MAXNV];
    set active[SETWORDSNEEDED(MAXNV)];
    int weight[MAXNV];
    boolean isreg,hasloops;

#ifdef PATHCOUNTS
    ++count0;
#endif

    k = -1;
    maxd = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
    {
        thisdeg = 0;
        for (j = 0; j < m; ++j) thisdeg += POPCOUNT(gi[j]);
        deg[i] = thisdeg;
        if (thisdeg > maxd) maxd = thisdeg;
        if (thisdeg == 0)    weight[i] = -i-1;
        else if (i < nfixed) weight[i] = 2*colour[i];
        else                 weight[i] = 2*colour[i]+1;
    }

    if (maxdeg >= 0 && maxd > maxdeg) return;

#ifdef PATHCOUNTS
    ++count1;
#endif

    if (ne == 0 && minedges <= 0
                && (!lswitch || (lswitch && (maxdeg&1) == 0)))
    {
        trythisone(NULL,lswitch,deg,maxdeg,0,n);
        return;
    }

#ifdef PATHCOUNTS
    ++count2;
#endif

    isreg = !lswitch && (maxdeg >= 0 && 2*minedges == n*(long)maxdeg);
        /* Case of regular multigraphs */

    if (isreg)  /* regular case */
    /* Condition: def(v) <= total def of neighbours */
    {
        for (i = 0; i < n; ++i)
        {
            def[i] = maxdeg - deg[i];
            delta[i] = -def[i];
        }

        for (i = 0; i < k; ++i)
        {
            x0 = v0[i]; x1 = v1[i];
            delta[x0] += def[x1];
            delta[x1] += def[x0];
        }

        for (i = 0; i < n; ++i) if (delta[i] < 0) return;
    }

    if ((isreg || lswitch) && (maxdeg & n & 1) == 1) return;
    // if (isreg && j0 >= 0 && maxdeg > 0) return;             DECODE
    // if (lswitch && j0 >= 0 && (maxdeg&1) == 1) return;

#ifdef PATHCOUNTS
    ++count3;
#endif

    if (maxedges == NOLIMIT)
    {
        if (maxmult == NOLIMIT) maxedges = maxdeg*n/2;
        else                    maxedges = ne*maxmult;
    }
    if (maxmult == NOLIMIT) maxmult = maxedges - ne + 1;
    if (maxdeg >= 0 && maxmult > maxdeg) maxmult = maxdeg;
    if (maxedges < ne || ne*maxmult < minedges) return;

#ifdef PATHCOUNTS
    ++count4;
#endif

    if (n > MAXNV || ne > MAXNE)
        gt_abort(">E multig: MAXNV or MAXNE exceeded\n");

    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    setlabptn(weight,lab,ptn,n);
    EMPTYSET(active,m);
    ADDELEMENT(active,0);

    for (i = 1; i < n; ++i)
        if (ptn[i-1] == 0) ADDELEMENT(active,i);

    hasloops = FALSE;
    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (ISELEMENT(gi,i))
        {
            hasloops = TRUE;
            break;
        }
 
    options.defaultptn = FALSE;
    options.digraph = hasloops;
    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;

    nauty(g,lab,ptn,active,orbits,&options,&stats,workspace,2*MAXNV,m,n,NULL);

    if (stats.grpsize2 == 0)
        groupsize = stats.grpsize1 + 0.1;
    else
        groupsize = 0;

    group = groupptr(FALSE);
    makecosetreps(group);

    lastrejok = FALSE;

    if (isreg)
        scan_reg(0,ne,minedges,maxedges,0,maxmult,group,n,delta,def,maxdeg);
    else if (lswitch)
        scan_lp(0,ne,minedges,maxedges,0,maxmult,group,n,deg,maxdeg);
    else if (maxdeg >= 0)
        scan_md(0,ne,minedges,maxedges,0,maxmult,group,n,deg,maxdeg);
    else
        scan(0,ne,minedges,maxedges,0,maxmult,group,n);
}

/**************************************************************************/

static boolean
readvcol(FILE *f, int *m, int *n, int *ne, graph *g)
/* Read the T style output as from vcolg.
   The edges go into v0[*],v1[*],eno[*][*] and the vertex colours into colour[*].
   Return FALSE iff the input is empty.
*/
{
    int i,j,nn,nne,mm,x,y;
    
    if (fscanf(f,"%d",&nn) != 1) return FALSE;
    if (nn > MAXNV) gt_abort(">E multig : too many vertices\n");

    if (fscanf(f,"%d",&nne) != 1)
        gt_abort(">E multig : incomplete input (1)\n");
    else if (nne > MAXNE)
        gt_abort(">E multig : too many edges\n");

    *n = nn;
    *ne = nne;
    *m = mm = SETWORDSNEEDED(nn);

    for (i = 0; i < nn; ++i)
        if (fscanf(f,"%d",&colour[i]) != 1)
            gt_abort(">E multig : incomplete input (2)\n");

    EMPTYGRAPH(g,mm,nn);

    if (Aswitch || Bswitch)
        for (i = 0; i < nn; ++i)
        for (j = 0; j < nn; ++j)
            edgeno[i][j] = -1;

    for (i = 0; i < nne; ++i)
    {
        if (fscanf(f,"%d%d",&v0[i],&v1[i]) != 2)
            gt_abort(">E multig : incomplete input (3)\n");
        ADDELEMENT(g+mm*v0[i],v1[i]);
        ADDELEMENT(g+mm*v1[i],v0[i]);
        edgeno[v0[i]][v1[i]] = edgeno[v1[i]][v0[i]] = i;
    }

    return TRUE;
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    graph *g;
    graph gg[MAXNV*SETWORDSNEEDED(MAXNV)];
    int i,m,n,ne,codetype;
    int argnum,j,nfixed,maxdeg,regdeg,ldeg;
    char *arg,sw;
    boolean badargs;
    boolean fswitch,uswitch,eswitch,qswitch,mswitch,Dswitch;
    boolean lswitch,rswitch;
    long minedges,maxedges,maxmult;
    double t;
    char *infilename,*outfilename;
    FILE *infile;
    char msg[201];
    int msglen;

    HELP; PUTVERSION;

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    rswitch = fswitch = Tswitch = Gswitch = Vswitch = FALSE;
    uswitch = eswitch = mswitch = qswitch = FALSE;
    lswitch = Aswitch = Bswitch = Dswitch = FALSE;
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
                     SWLONG('m',mswitch,maxmult,"multig -m")
                else SWBOOLEAN('q',qswitch)
                else SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('T',Tswitch)
                else SWBOOLEAN('G',Gswitch)
                else SWBOOLEAN('A',Aswitch)
                else SWBOOLEAN('B',Bswitch)
                else SWBOOLEAN('V',Vswitch)
                else SWINT('f',fswitch,nfixed,"multig -f")
                else SWINT('D',Dswitch,maxdeg,"multig -D")
                else SWINT('r',rswitch,regdeg,"multig -r")
                else SWINT('l',lswitch,ldeg,"multig -l")
                else SWRANGE('e',":-",eswitch,minedges,maxedges,"multig -e")
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

    if (badargs || argnum > 2)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if ((Gswitch!=0) + (Tswitch!=0) + (uswitch!=0)
          + (Aswitch!=0) + (Bswitch!=0) >= 2)
        gt_abort(">E multig: -G, -T, -A, -B and -u are incompatible\n");

#ifndef OUTPROC
    if (!Tswitch && !Gswitch && !Aswitch && !Bswitch && !uswitch)
        gt_abort(">E multig: must use -A, -B, -T, -G or -u\n");
#endif

    if (rswitch && (Dswitch || eswitch))
        gt_abort(">E multig: -r is incompatible with -D and -e\n");

    if (lswitch && (rswitch || Dswitch || eswitch))
        gt_abort(">E multig: -l is incompatible with -r, -D and -e\n");

    if (!eswitch)
    {
        minedges = 0;
        maxedges = NOLIMIT;
    }
    if (!mswitch) maxmult = NOLIMIT;
    if (!fswitch) nfixed = 0;

    if (Bswitch && nfixed == 0)
        gt_abort(">E multig: -B requires -f# with #>0\n");
    if (fswitch) Brows = nfixed;

    if (maxedges >= NOLIMIT && maxmult >= NOLIMIT
                  && !Dswitch && !rswitch && !lswitch)
        gt_abort(
      ">E multig: either -D or -e or -m or -r must impose a real limit\n");

    if (!qswitch)
    {
        msg[0] = '\0';
        CATMSG0(">A multig");
        if (eswitch || mswitch || uswitch || (fswitch && nfixed > 0)
              || lswitch || rswitch || Dswitch || Tswitch
              || Gswitch || Aswitch || Bswitch || Vswitch)
            CATMSG0(" -");
        if (mswitch) CATMSG1("m%ld",maxmult);
        if (Vswitch) CATMSG0("V");
        if (uswitch) CATMSG0("u");
        if (Tswitch) CATMSG0("T");
        if (Gswitch) CATMSG0("G");
        if (Aswitch) CATMSG0("A");
        if (Bswitch) CATMSG0("B");
        if (fswitch) CATMSG1("f%d",nfixed);
        if (eswitch) CATMSG2("e%ld:%ld",minedges,maxedges);
        if (Dswitch) CATMSG1("D%d",maxdeg);
        if (rswitch) CATMSG1("r%d",regdeg);
        if (lswitch) CATMSG1("l%d",ldeg);
        msglen = strlen(msg);
        if (argnum > 0) msglen += strlen(infilename);
        if (argnum > 1) msglen += strlen(outfilename);
        if (msglen >= 196)
        {
            fputs(msg,stderr);
            if (argnum > 0) fprintf(stderr," %s",infilename);
            if (argnum > 1) fprintf(stderr," %s",outfilename);
            fprintf(stderr,"\n");
        }
        else
        {
            if (argnum > 0) CATMSG1(" %s",infilename);
            if (argnum > 1) CATMSG1(" %s",outfilename);
            CATMSG0("\n");
            fputs(msg,stderr);
        }
        fflush(stderr);
    }

    if (rswitch)
    {
        eswitch = Dswitch = TRUE;
        maxdeg = regdeg;
    }

    if (lswitch)
    {
        eswitch = Dswitch = TRUE;
        maxdeg = ldeg;
    }

    if (infilename && infilename[0] == '-') infilename = NULL;

    if (Vswitch)
    {
        if (!infilename)
            infile = stdin;
        else if ((infile = fopen(infilename,"r")) == NULL)
            gt_abort(">E multig -V : cannot open input file\n");
        codetype = 0;
    }
    else
    {
        infile = opengraphfile(infilename,&codetype,FALSE,1);
        if (!infile) exit(1);
    }
    if (!infilename) infilename = "stdin";

    NODIGRAPHSYET(codetype);

    if (uswitch)
        outfile = NULL;
    else
    {
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
    }

    mg_nin = mg_nout = mg_skipped = 0;

    if (Vswitch) g = gg;

    t = CPUTIME;
    while (TRUE)
    {
        if (Vswitch)
        {
            if (!readvcol(infile,&m,&n,&ne,g)) break;
        }
        else
        {
            if ((g = readg(infile,NULL,0,&m,&n)) == NULL) break;
            if (n > MAXNV) gt_abort(">E multig: increase MAXNV\n");
            for (i = 0; i < n; ++i) colour[i] = 0;
        }
        ++mg_nin;
#ifdef PATHCOUNTS
        oldlo = mg_nout;
#endif

#ifdef INPUTGRAPH
        if (!INPUTGRAPH(g,m,n))
        {
           ++mg_skipped;
           if (!Vswitch) FREES(g);
           continue;
        }
#endif

#ifdef INPUTGRAPHC
        if (!INPUTGRAPHC(g,colour,m,n))
        {
           ++mg_skipped;
           if (!Vswitch) FREES(g);
           continue;
        }
#endif
 
        if (rswitch)
        {
            minedges = ((long)n * (long)regdeg + 1) / 2;
            maxedges = ((long)n * (long)regdeg) / 2;
        }
        if (lswitch)
        {
            maxedges = ((long)n * (long)ldeg) / 2;
            if ((ldeg & 1) == 1) minedges = (n + 1) / 2;
            else                 minedges = 0;
        }

        if (Vswitch) vmulti(g,nfixed,minedges,maxedges,maxmult,
                                     (Dswitch?maxdeg:-1),lswitch,ne,m,n);
        else         multi(g,nfixed,minedges,maxedges,maxmult,
                                     (Dswitch?maxdeg:-1),lswitch,m,n);
#ifdef PATHCOUNTS
        if (mg_nout != oldlo) ++count5;
#endif
        if (!uswitch && ferror(outfile))
            gt_abort(">E multig output error\n");
        if (!Vswitch) FREES(g);
    }
    t = CPUTIME - t;

#ifdef SUMMARY
    SUMMARY();
#endif

    if (!qswitch)
    {
        fprintf(stderr,">Z " COUNTER_FMT " graphs read from %s; ",
                        mg_nin,infilename);
#ifdef INPUTGRAPH
        fprintf(stderr,COUNTER_FMT " skipped; ",mg_skipped);
#endif
        PRINT_COUNTER(stderr,mg_nout);
        if (!uswitch)
            fprintf(stderr," multigraphs written to %s",outfilename);
        else
            fprintf(stderr," multigraphs generated");
        fprintf(stderr,"; %.2f sec\n",t);
    }

#ifdef GROUPTEST
    fprintf(stderr,"Group test = %lld\n",totallab);
#endif

#ifdef PATHCOUNTS
    fprintf(stderr,"Counts: %lld %lld %lld %lld %lld %lld\n",
            count0,count1,count2,count3,count4,count5);
#endif

    exit(0);
}
