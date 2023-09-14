/* vcolg.c version 3.1; B D McKay, Apr 24, 2021 */

#define USAGE \
"vcolg [-q] [-u|-T|-o] [-e#|-e#:#] [-m#] [-c#,..,#] [-f#] [infile [outfile]]"

#define HELPTEXT \
"  Read graphs or digraphs and colour their vertices in\n\
  all possible ways with colours 0,1,2,... .\n\
  Isomorphic graphs derived from the same input are suppressed.\n\
  If the input graphs are non-isomorphic then the output graphs are also.\n\
\n\
    -e# | -e#:#  specify a value or range of the total value of the colours\n\
    -m# number of available colours (default 2 if -c not given)\n\
    -c#,..,#  specify the maximum number of vertices of each colour\n\
        The total must at least equal the number of vertices in the input.\n\
    -d#,..,#  minimum vertex degree for each colour (out-degree for digraphs)\n\
    -D#,..,#  maximum vertex degree for each colour (out-degree for digraphs)\n\
         -d and -D can have fewer colours than -m/-c but not more\n\
    -f# Use the group that fixes the first # vertices setwise\n\
    -T  Use a simple text output format (nv ne {col} {v1 v2})\n\
    -o  Use sparse6 (undirected) or digraph6 (directed) for output,\n\
          provided m=2 and the inputs have no loops.\n\
    -u  no output, just count them\n\
    -q  suppress auxiliary information\n"

/*************************************************************************/

#include "gtools.h"
#include "naugroup.h"
#include "nautinv.h"
#include "gutils.h"

nauty_counter vc_nin,vc_nout;
FILE *outfile;

#define MAXNV 128   /* Maximum number of vertices */

static int col[MAXNV];
static boolean first;
static int lastreject[MAXNV];
static boolean lastrejok;
static unsigned long groupsize;
static unsigned long newgroupsize;
static boolean Tswitch,oswitch;

static int fail_level;

#define GROUPTEST_NOT 
#ifdef GROUPTEST
static long long totallab;
#endif

/* If OUTPROC is defined at compile time, and -u is not used, the
 * procedure OUTPROC is called for each graph.  This must be linked
 * by the compiler.  The arguments are
 * f = open output file
 * g = the input graph
 * col[0..n-1] = the colours
 * m,n = usual nauty meanings
 */

/* SUMMARY feature
 *
 * If SUMMARY is defined, it must expand as the name of a procedure
 * with prototype  void SUMMARY(void).  It is called at the end before
 * the normal summary (which can be suppressed with -q).  The numbers of
 * graphs read and coloured graphs produced are available in the global
 * variables vc_nin and vc_nout (type nauty_counter).
 */

#ifdef OUTPROC
extern void OUTPROC(FILE*,graph*,int*,int,int);
#endif

#ifdef SUMMARY
extern void SUMMARY(void);
#endif

/**************************************************************************/

static int
ismax(int *p, int n)
/* test if col^p <= col */
{
    int i,k;
    int fail;

    fail = 0;
    for (i = 0; i < n; ++i)
    {
        k = p[i];
        if (k > fail) fail = k;
        if (col[k] > col[i])
        {
            fail_level = fail;
            return FALSE;
        }
        else if (col[k] < col[i]) return TRUE;
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
writeone(FILE *outfile, graph *g, int *col, boolean digraph, int m, int n)
{
    int i;
    set *gi;

    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (col[i] == 1) ADDELEMENT(gi,i);

    if (digraph)
        writed6(outfile,g,m,n);
    else
        writes6(outfile,g,m,n);

    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (col[i] == 1) DELELEMENT(gi,i);
}

static int
trythisone(grouprec *group, graph *g, boolean digraph, int m, int n)
/* Try one solution, accept if maximal. */
/* Return value is level to return to. */
{
    int i,j;
    boolean accept;
    graph *gi;
    size_t ne;
    char s[20],*p;
    DYNALLSTAT(char,line,line_sz);

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

        ++vc_nout;

        if (oswitch)
            writeone(outfile,g,col,digraph,m,n);
        else if (outfile)
        {
#ifdef OUTPROC
            OUTPROC(outfile,g,col,m,n);
#else
            ne = 0;
            for (gi = g + m*(size_t)n; --gi >= g; )
                ne += POPCOUNT(*gi);
            if (!digraph)
            {
                for (i = 0, gi = g; i < n; ++i, gi += m)
                    if (ISELEMENT(gi,i)) ++ne;
                ne /= 2;
            }
#define PUTINT(xx) { unsigned long ul = (xx); char *sp; \
 if (ul == 0) *(p++) = '0'; \
 else { sp = s; while (ul) { *(sp++) = (ul % 10) + '0'; ul /= 10; } \
        while (sp > s) { *(p++) = *(--sp); } }}
#define SPC *(p++) = ' '

            DYNALLOC1(char,line,line_sz,20*(n+ne)+50,"vcolg output");
            p = line;
            PUTINT(n); SPC; PUTINT(ne);
            for (i = 0; i < n; ++i) { SPC; PUTINT(col[i]); }
            SPC;
            for (i = 0, gi = g; i < n; ++i, gi += m)
            {
                for (j = (digraph?-1:i-1); (j = nextelement(gi,m,j)) >= 0; )
                {
                    SPC; PUTINT(i); SPC; PUTINT(j);
                }
            }
            *(p++) = '\n';
            *(p++) = '\0';
            fputs(line,outfile);
#endif
        }
        return n-1;
    }
    else
        return fail_level-1;
}

/**************************************************************************/

static int
scan(int level, graph *g, boolean digraph, int *prev, long mincols,
     long maxcols, long sofar, long *colcount, long *mindeg,
     long *maxdeg, long *deg, int numcols, grouprec *group, int m, int n)
/* Recursive scan for default case */
/* Returned value is level to return to. */
{
    int left;
    long min,max,k,ret;

    if (level == n)
        return trythisone(group,g,digraph,m,n);

    left = n - level - 1;
    min = mincols - sofar - numcols*left;
    if (min < 0) min = 0;
    max = maxcols - sofar;
    if (max >= numcols) max = numcols - 1;
    if (prev[level] >= 0 && col[prev[level]] < max)
        max = col[prev[level]];

    for (k = min; k <= max; ++k)
    {
        if (colcount[k] <= 0) continue;
        if (mindeg[k] > deg[level] || maxdeg[k] < deg[level]) continue;
        --colcount[k];
        col[level] = k;
        ret = scan(level+1,g,digraph,prev,mincols,maxcols,
                       sofar+k,colcount,mindeg,maxdeg,deg,numcols,group,m,n);
        ++colcount[k];
        if (ret < level) return ret;
    }

    return level-1;
}

/**************************************************************************/

static void
colourgraph(graph *g, int nfixed, long mincols, long maxcols,
         long *colcount, long *mindeg, long *maxdeg, long *deg,
         int numcols, int m, int n)
{
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[8*MAXNV];
    grouprec *group;
    int i,j,k,nloops;
    set *gi,*gj;
    int lab[MAXNV],ptn[MAXNV],orbits[MAXNV];
    boolean loop[MAXNV];
    int prev[MAXNV]; /* If >= 0, earlier point that must have greater colour */
    int weight[MAXNV];
    int region,start,stop;

    if (n > MAXNV) gt_abort(">E vcolg: MAXNV exceeded\n");
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    nloops = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (ISELEMENT(gi,i))
        {
            DELELEMENT(gi,i);
            loop[i] = TRUE;
            ++nloops;
        }
        else
            loop[i] = FALSE;

    for (region = 0; region < 2; ++region)
    {
        if (region == 0)
        {
            if (nfixed == 0) continue;
            start = 0;
            stop = nfixed;
            if (stop > n) stop = n;
        }
        else
        {
            if (nfixed >= n) continue;
            start = nfixed;
            stop = n;
        }
        
        for (i = start, gi = g + m*(size_t)start; i < stop; ++i, gi += m)
        {
            /* Find most recent equivalent j. */
            for (j = i-1, gj = gi-m; j >= start; --j, gj -= m)
            {
                if (loop[j] != loop[i]) continue;
                for (k = 0; k < m; ++k) if (gi[k] != gj[k]) break;
                if (k < m)
                {
                    FLIPELEMENT(gi,i); FLIPELEMENT(gj,j);
                    for (k = 0; k < m; ++k) if (gi[k] != gj[k]) break;
                    FLIPELEMENT(gi,i); FLIPELEMENT(gj,j);
                }
                if (k == m) break;
            }
            if (j >= start)
            {
                prev[i] = j;
                weight[i] = weight[j] + 1;
            }
            else
            {
                prev[i] = -1;
                weight[i] = 0;
            }
        }
    }

    if (n == 0)
    {
        scan(0,g,FALSE,prev,mincols,maxcols,0,colcount,
                               NULL,NULL,NULL,numcols,FALSE,m,n);
        return;
    }

    for (i = nfixed; i < n; ++i) weight[i] += nfixed;

    if (maxcols == NOLIMIT || maxcols > n*numcols) maxcols = n*numcols;
    if (n*numcols < mincols) return;

    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;
    options.defaultptn = FALSE;
    options.digraph = (nloops > 0);

    setlabptn(weight,lab,ptn,n);

    if (nloops > 0)
        for (i = 0, gi = g; i < n; ++i, gi += m)
            if (loop[i]) ADDELEMENT(gi,i);
 
    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,8*MAXNV,m,n,NULL);

    if (stats.grpsize2 == 0)
        groupsize = stats.grpsize1 + 0.1;
    else
        groupsize = 0;

    group = groupptr(FALSE);
    makecosetreps(group);

    if (stats.numorbits < n)
    {
        j = n;
        for (i = 0; i < n; ++i)
            if (orbits[i] < i && orbits[i] < j) j = orbits[i];

        for (i = j + 1; i < n; ++i)
            if (orbits[i] == j) prev[i] = j;
    }

    lastrejok = FALSE;
    for (i = 0; i < n; ++i) col[i] = 0;

    scan(0,g,FALSE,prev,mincols,maxcols,0,colcount,
                           mindeg,maxdeg,deg,numcols,group,m,n);
}

/**************************************************************************/

static void
colourdigraph(graph *g, int nfixed, long mincols, long maxcols,
         long *colcount, long *mindeg, long *maxdeg, long *deg,
         int numcols, int m, int n)
{
    static DEFAULTOPTIONS_GRAPH(options);
    statsblk stats;
    setword workspace[8*MAXNV];
    grouprec *group;
    int i,j,k,nloops;
    size_t ii;
    set *gi,*gj,*gci,*gcj;
    int lab[MAXNV],ptn[MAXNV],orbits[MAXNV];
    boolean loop[MAXNV];
    int prev[MAXNV]; /* If >= 0, earlier point that must have greater colour */
    int weight[MAXNV];
    int region,start,stop;
    DYNALLSTAT(graph,gconv,gconv_sz);

    if (n > MAXNV) gt_abort(">E vcolg: MAXNV exceeded\n");
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    DYNALLOC2(graph,gconv,gconv_sz,n,m,"colourdigraph");

    nloops = 0;
    for (i = 0, gi = g; i < n; ++i, gi += m)
        if (ISELEMENT(gi,i))
        {
            DELELEMENT(gi,i);
            loop[i] = TRUE;
            ++nloops;
        }
        else
            loop[i] = FALSE;

    for (ii = 0; ii < m*(size_t)n; ++ii) gconv[ii] = g[ii];
    converse(gconv,m,n);

    for (region = 0; region < 2; ++region)
    {
        if (region == 0)
        {
            if (nfixed == 0) continue;
            start = 0;
            stop = nfixed;
            if (stop > n) stop = n;
        }
        else
        {
            if (nfixed >= n) continue;
            start = nfixed;
            stop = n;
        }
        
        for (i = start,
                    gi = g + m*(size_t)start, gci = gconv + m*(size_t)start;
             i < stop; ++i, gi += m, gci += m)
        {
            /* Find most recent equivalent j. */
            for (j = i-1, gj = gi-m, gcj = gci-m; j >= start;
                                                   --j, gj -= m, gcj -= m)
            {
                if (loop[j] != loop[i]
                       || ISELEMENT(gi,j) != ISELEMENT(gj,i)) continue;
                for (k = 0; k < m; ++k)
                     if (gi[k] != gj[k] || gci[k] != gcj[k]) break;
                if (k < m)
                {
                    FLIPELEMENT(gi,i); FLIPELEMENT(gj,j);
                    FLIPELEMENT(gci,i); FLIPELEMENT(gcj,j);
                    for (k = 0; k < m; ++k)
                        if (gi[k] != gj[k] || gci[k] != gcj[k]) break;
                    FLIPELEMENT(gci,i); FLIPELEMENT(gcj,j);
                    FLIPELEMENT(gi,i); FLIPELEMENT(gj,j);
                }
                if (k == m) break;
            }
            if (j >= start)
            {
                prev[i] = j;
                weight[i] = weight[j] + 1;
            }
            else
            {
                prev[i] = -1;
                weight[i] = 0;
            }
        }
    }

    for (i = nfixed; i < n; ++i) weight[i] += nfixed;

    if (maxcols == NOLIMIT || maxcols > n*(long)numcols)
        maxcols = n*(long)numcols;
    if (n*(long)numcols < mincols) return;

    if (n == 0)
    {
        scan(0,g,TRUE,prev,mincols,maxcols,0,colcount,
             mindeg,maxdeg,deg,numcols,FALSE,m,n);
        return;
    }

    options.userautomproc = groupautomproc;
    options.userlevelproc = grouplevelproc;
    options.defaultptn = FALSE;
    options.digraph = TRUE;
    options.invarproc = adjacencies;
    options.maxinvarlevel = n;

    setlabptn(weight,lab,ptn,n);

    if (nloops > 0)
        for (i = 0, gi = g; i < n; ++i, gi += m)
            if (loop[i]) ADDELEMENT(gi,i);
 
    nauty(g,lab,ptn,NULL,orbits,&options,&stats,workspace,8*MAXNV,m,n,NULL);

    if (stats.grpsize2 == 0)
        groupsize = stats.grpsize1 + 0.1;
    else
        groupsize = 0;

    group = groupptr(FALSE);
    makecosetreps(group);

    if (stats.numorbits < n)
    {
        j = n;
        for (i = 0; i < n; ++i)
            if (orbits[i] < i && orbits[i] < j) j = orbits[i];

        for (i = j + 1; i < n; ++i)
            if (orbits[i] == j) prev[i] = j;
    }

    lastrejok = FALSE;
    for (i = 0; i < n; ++i) col[i] = 0;

    scan(0,g,TRUE,prev,mincols,maxcols,0,colcount,
                                mindeg,maxdeg,deg,numcols,group,m,n);
}

/**************************************************************************/

int
main(int argc, char *argv[])
{
    graph *g,*gi;
    int m,n,codetype;
    int argnum,i,j,nfixed;
    char *arg,sw;
    boolean badargs,digraph,cswitch,dswitch,Dswitch;
    boolean fswitch,uswitch,eswitch,qswitch,mswitch;
    long mincols,maxcols,totcols;
    int numcols;
    double t;
    char *infilename,*outfilename;
    FILE *infile;
    char msg[201];
    int msglen,collen,dlen,Dlen;
    long colcount[MAXNV],mindeg[MAXNV],maxdeg[MAXNV],deg[MAXNV];

    HELP; PUTVERSION;

    nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

    fswitch = Tswitch = oswitch = cswitch = dswitch = FALSE;
    uswitch = eswitch = mswitch = qswitch = Dswitch = FALSE;
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
                     SWINT('m',mswitch,numcols,"vcolg -m")
                else SWBOOLEAN('q',qswitch)
                else SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('T',Tswitch)
                else SWBOOLEAN('o',oswitch)
                else SWSEQUENCEMIN('c',",",cswitch,colcount,1,MAXNV,collen,"vcolg -c")
                else SWSEQUENCEMIN('d',",",dswitch,mindeg,1,MAXNV,dlen,"vcolg -d")
                else SWSEQUENCEMIN('D',",",Dswitch,maxdeg,1,MAXNV,Dlen,"vcolg -D")
                else SWINT('f',fswitch,nfixed,"vcolg -f")
                else SWRANGE('e',":-",eswitch,mincols,maxcols,"vcolg -e")
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

    if ((Tswitch!=0) + (oswitch!=0) + (uswitch!=0) >= 2)
        gt_abort(">E vcolg: -T, -o and -u are incompatible\n");

#ifndef OUTPROC
    if (!Tswitch && !oswitch && !uswitch)
        gt_abort(">E vcolg: must use -T, -o or -u\n");
#endif

    if (!mswitch) numcols = 2;

    if (oswitch && numcols != 2)
        gt_abort(">E vcolg: -o is only allowed for 2 colours\n");

    if (!eswitch)
    {
        mincols = 0;
        maxcols = NOLIMIT;
    }
    if (!fswitch) nfixed = 0;

    if (cswitch && mswitch && numcols != collen)
        gt_abort(">E vcolg: -m and -c disagree on number of colours\n");

    if (cswitch)
    {
        numcols = collen;
        totcols = 0;
        for (i = 0; i < numcols; ++i)
            if (colcount[i] < 0)
                gt_abort(">E vcolg: negative counts not allowed for -c\n");
            else
            {
                totcols += colcount[i];
                if (totcols < 0) { totcols = NOLIMIT; break; } /* catch overflow */
            }
    }
    else
        for (i = 0; i < numcols; ++i) colcount[i] = NOLIMIT;

    if (dswitch && dlen > numcols)
        gt_abort(">E vcolg: -d has too many colours\n");
    if (Dswitch && Dlen > numcols)
        gt_abort(">E vcolg: -D has too many colours\n");

    for (i = (dswitch ? dlen : 0); i < numcols; ++i)
        mindeg[i] = 0;
    for (i = (Dswitch ? Dlen : 0); i < numcols; ++i)
        maxdeg[i] = NOLIMIT;
    for (i = 0; i < numcols; ++i)
        if (mindeg[i] > maxdeg[i])
            gt_abort(">E vcolg : contradictory bound from -d/-D\n");

    if (cswitch && !qswitch)
    { 
        fprintf(stderr,">c"); for (i = 0; i < numcols; ++i)
            fprintf(stderr," %ld",colcount[i]);
        fprintf(stderr,"\n");
    }

    if (dswitch && !qswitch)
    { 
        fprintf(stderr,">d");
        for (i = 0; i < numcols; ++i)
            fprintf(stderr," %ld",mindeg[i]);
        fprintf(stderr,"\n");
    }

    if (Dswitch && !qswitch)
    { 
        fprintf(stderr,">D");
       
        for (i = 0; i < numcols; ++i)
            if (maxdeg[i] == NOLIMIT) fprintf(stderr," -");
            else                      fprintf(stderr," %ld",maxdeg[i]);
        fprintf(stderr,"\n");
    }

    if (!qswitch)
    {
        msg[0] = '\0';
        CATMSG0(">A vcolg");
        if (eswitch || mswitch || uswitch || (fswitch && nfixed > 0)
              || Tswitch || oswitch)
            CATMSG0(" -");
        if (mswitch) CATMSG1("m%d",numcols);
        if (uswitch) CATMSG0("u");
        if (Tswitch) CATMSG0("T");
        if (oswitch) CATMSG0("o");
        if (fswitch) CATMSG1("f%d",nfixed);
        if (eswitch) CATMSG2("e%ld:%ld",mincols,maxcols);
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

    if (infilename && infilename[0] == '-') infilename = NULL;
    infile = opengraphfile(infilename,&codetype,FALSE,1);
    if (!infile) exit(1);
    if (!infilename) infilename = "stdin";

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

    vc_nin = vc_nout = 0;

    t = CPUTIME;
    while (TRUE)
    {
        if ((g = readgg(infile,NULL,0,&m,&n,&digraph)) == NULL) break;
        ++vc_nin;

        if (cswitch && n > totcols)
            gt_abort(">E vcolg: not enough colours for input\n");

        if (oswitch && loopcount(g,m,n) > 0)
            gt_abort(">E vcolg: loops in input are not allowed for -o\n");

        for (i = 0, gi = g; i < n; ++i, gi += m)
        {
            deg[i] = 0;
            for (j = 0; j < m; ++j) deg[i] += POPCOUNT(gi[j]);
        }

        if (!digraph)
            colourgraph(g,nfixed,mincols,maxcols,colcount,
                                   mindeg,maxdeg,deg,numcols,m,n);
        else
            colourdigraph(g,nfixed,mincols,maxcols,colcount,
                                   mindeg,maxdeg,deg,numcols,m,n);

        if (!uswitch && ferror(outfile)) gt_abort(">E vcolg output error\n");
        FREES(g);
    }
    t = CPUTIME - t;

#ifdef SUMMARY
    SUMMARY();
#endif

    if (!qswitch)
    {
        fprintf(stderr,">Z ");
        PRINT_COUNTER(stderr,vc_nin);
        fprintf(stderr," graphs read from %s",infilename);
        fprintf(stderr,"; ");
        PRINT_COUNTER(stderr,vc_nout);
        if (!uswitch)
            fprintf(stderr," coloured graphs written to %s",outfilename);
        else
            fprintf(stderr," coloured graphs generated");
        fprintf(stderr,"; %.2f sec\n",t);
    }

#ifdef GROUPTEST
    fprintf(stderr,"Group test = %lld\n",totallab);
#endif

    exit(0);
}
