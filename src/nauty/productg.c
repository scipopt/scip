/* productg.c  version 1.2; B D McKay, October 2022. */

#define USAGE "productg [-c|-l|-L|-k|-t|-a#] [infile [outfile]]"

#define HELPTEXT \
" Read two graphs in graph6/sparse6 format and write their product\n\
  in sparse6 or dreadnaut format.\n\
  -d : Dreadnaut format (default is sparse6)\n\
\n\
  -c : Cartesian product\n\
  -l : Lexicographic product (G1[G2])\n\
  -L : Lexicographic product (G2[G1])\n\
  -t : Tensor (Kronecker, categorical) product\n\
  -k : Complete (strong, normal) product\n\
  -a# : general case (argument is a 3-digit octal number)\n\
        Add these values giving the condition for an edge:\n\
   Code:    400  200   100   040 020 010      004   002    001\n\
   Graph1: same same  same   adj adj adj    nonadj nonadj nonadj\n\
   Graph2: same  adj nonadj same adj nonadj  same   adj   nonadj\n"

/*************************************************************************/

#undef MAXN
#define MAXN 0
#include "gtools.h" 
#include "gutils.h"

#define ISADJ1(v1,w1) ISELEMENT(g1+m1*(v1),w1)
#define ISADJ2(v2,w2) ISELEMENT(g2+m2*(v2),w2)
#define ISEQ(x,y) ((x)==(y))
#define N(x,y) ((x)*(long)n2+(y))   /* Vertex number of (x,y) */

#define VCODE1(v1,w1) (ISEQ(v1,w1) ? 0 : ISADJ1(v1,w1) ? 1 : 2)
#define VCODE2(v2,w2) (ISEQ(v2,w2) ? 0 : ISADJ2(v2,w2) ? 1 : 2)
#define ECODE(v1,v2,w1,w2) (3*VCODE1(v1,w1)+VCODE2(v2,w2))
#define ADJ(v1,v2,w1,w2) (adjcode & (0400 >> ECODE(v1,v2,w1,w2)))

/* 
 *  Code:   0400  0200  0100  0040  0020  0010  0004  0002  0001
 *  Graph1:  eq    eq    eq    adj   adj   adj   non   non   non 
 *  Graph2:  eq    adj   non    eq   adj   non   eq    adj   non
 *
 *  Cartesian product = 0200+0040 = 0240
 *  Tensor product = 0020
 *  Lexicographic product = 0200+0070=0270  or  0040+0222 = 0262
 *  Complete product = 0260
 */

#define SWOCT(c,bool,val,id) if (sw==c) {bool=TRUE;arg_oct(&arg,&val,id);}

/**************************************************************************/

static int
longoctvalue(char **ps, long *l)
{
    boolean neg,pos;
    long sofar,last;
    char *s;

    s = *ps;
    pos = neg = FALSE;
    if (*s == '-')
    {
        neg = TRUE;
        ++s;
    }
    else if (*s == '+')
    {
        pos = TRUE;
        ++s;
    }

    if (*s < '0' || *s > '7') 
    {
        *ps = s;
        return (pos || neg) ? ARG_ILLEGAL : ARG_MISSING;
    }

    sofar = 0;

    for (; *s >= '0' && *s <= '7'; ++s)
    {
        last = sofar;
        sofar = sofar * 8 + (*s - '0');
        if (sofar < last || sofar > MAXARG)
        {
            *ps = s;
            return ARG_TOOBIG;
        }
    }
    *ps = s;
    *l = neg ? -sofar : sofar;
    return ARG_OK;
}

/*************************************************************************/

void
arg_oct(char **ps, int *val, char *id)
{
    int code;
    long longval;

    code = longoctvalue(ps,&longval);
    *val = longval;
    if (code == ARG_MISSING || code == ARG_ILLEGAL)
    {
        fprintf(stderr,">E %s: missing argument value\n",id);
        gt_abort(NULL);
    }
    else if (code == ARG_TOOBIG || *val != longval)
    {
        fprintf(stderr,">E %s: argument value too large\n",id);
        gt_abort(NULL);
    }
}

/************************************************************************/

int
main(int argc, char *argv[])
{
    char *infilename,*outfilename;
    FILE *infile,*outfile;
    boolean badargs;
    int i,j,m1,n1,m2,n2,argnum;
    int codetype;
    graph *g1,*g2;
        char *arg,sw;
    int nv,v1,v2,w1,w2;
    boolean aswitch,cswitch,lswitch,Lswitch,tswitch,kswitch;
    boolean digraph,dreadnaut;
    int adjcode;
    SG_DECL(sg);
    int *d,*e,xx,yy;
    size_t *v,twone;

    aswitch = cswitch = lswitch = Lswitch = tswitch = kswitch = FALSE;
    dreadnaut = FALSE;

    HELP;

    infilename = outfilename = NULL;
    badargs = FALSE;

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
                    SWBOOLEAN('c',cswitch)
                else SWBOOLEAN('t',tswitch)
                else SWBOOLEAN('l',lswitch)
                else SWBOOLEAN('L',Lswitch)
                else SWBOOLEAN('k',kswitch)
                else SWBOOLEAN('d',dreadnaut)
                else SWOCT('a',aswitch,adjcode,"productg -a")
                else
                   badargs = TRUE;
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

    if ((aswitch!=0) + (cswitch!=0) + (lswitch!=0)  + (kswitch!=0) 
             + (Lswitch!=0) + (tswitch!=0) + (aswitch!=0) != 1)
    {
        fprintf(stderr,
                ">E productg: exactly one of -k,-c,-l,-L,-t,-a# is needed\n");
        exit(1);
    }

    if (cswitch) adjcode = 0240;
    if (tswitch) adjcode = 0020;
    if (lswitch) adjcode = 0270;
    if (Lswitch) adjcode = 0262;
    if (kswitch) adjcode = 0260;

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

    if ((g1 = readgg(infile,NULL,0,&m1,&n1,&digraph)) == NULL)
        gt_abort(">E first graph not found\n");
    if (digraph) gt_abort(">E productg does not support digraphs yet.\n");

    if ((g2 = readgg(infile,NULL,0,&m2,&n2,&digraph)) == NULL)
        gt_abort(">E second graph not found\n");
    if (digraph) gt_abort(">E productg does not support digraphs yet.\n");

    if ((long)n1 * (long)n2 > NAUTY_INFINITY-2)
        gt_abort(">E product would be too large\n");

    if (dreadnaut)
    {
        fprintf(outfile,"n=%ld $=0 g\n",n1*(long)n2);

        for (v1 = 0; v1 < n1; ++v1)
        for (v2 = 0; v2 < n2; ++v2)
        {
            fprintf(outfile,"%ld:",N(v1,v2));
            for (w1 = v1; w1 < n1; ++w1)
            for (w2 = 0; w2 < n2; ++w2)
            {
                if ((w1 > v1 || w2 > v2) && ADJ(v1,v2,w1,w2))
                    fprintf(outfile," %ld",N(w1,w2));
            }
            fprintf(outfile,"\n");
        }
        fprintf(outfile,". $$\n");
    }
    else
    {
        nv = n1*n2;
        SG_ALLOC(sg,nv,0,"productg");
        v = sg.v;
        d = sg.d;
        for (i = 0; i < nv; ++i) d[i] = 0;
        for (v1 = 0; v1 < n1; ++v1)
        for (v2 = 0; v2 < n2; ++v2)
        {
            xx = N(v1,v2);
            for (w1 = v1; w1 < n1; ++w1)
            for (w2 = 0; w2 < n2; ++w2)
            {
                if ((w1 > v1 || w2 > v2) && ADJ(v1,v2,w1,w2))
                {
                    yy = N(w1,w2);
                    ++d[xx];
                    ++d[yy];
                }
            }
        }
        twone = 0;
        for (i = 0; i < nv; ++i) twone += d[i];
        SG_ALLOC(sg,nv,twone,"productg");
        if (v != sg.v || d != sg.d)
            gt_abort(">E productg: SG_ALLOC problem\n");
        e = sg.e;

        v[0] = 0;
        for (i = 1; i < nv; ++i) v[i] = v[i-1] + d[i-1];
        for (i = 0; i < nv; ++i) d[i] = 0;

        for (v1 = 0; v1 < n1; ++v1)
        for (v2 = 0; v2 < n2; ++v2)
        {
            xx = N(v1,v2);
            for (w1 = v1; w1 < n1; ++w1)
            for (w2 = 0; w2 < n2; ++w2)
            {
                if ((w1 > v1 || w2 > v2) && ADJ(v1,v2,w1,w2))
                {
                    yy = N(w1,w2);
                    e[v[xx]+d[xx]++] = yy;
                    e[v[yy]+d[yy]++] = xx;
                }
            }
        }

        sg.nv = nv;
        sg.nde = twone;
        writes6_sg(outfile,&sg);
    }

    fprintf(stderr,">Z Wrote graph of order %d to %s\n",
                       n1*n2,outfilename);

    exit(0);
}
