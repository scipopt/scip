/* gentree version 1.3; Brendan McKay Oct 2022 */
/* This program is a wrapper for the program FreeTrees.c written
 * by Gang Li & Frank Ruskey.  See below for their original
 * comments. */
/* TODO: splitlevinc */

#define USAGE \
 "gentreeg [-D#] [-Z#:#] [-ulps] [-q] n [res/mod] [file]"

#define HELPTEXT \
" Generate (unrooted) trees.\n\
\n\
  n, n1:n2 : the number of vertices or a range\n\
              Outputs are in order of the number of vertices.\n\
   res/mod : only generate subset res out of subsets 0..mod-1\n\
\n\
     -D#   : an upper bound for the maximum degree\n\
     -Z#:# : bounds for the diameter\n\
\n\
     -s    : use sparse6 output (default)\n\
     -p    : write a parent array\n\
     -l    : write a level array \n\
     -u    : do not output any graphs, just generate and count them\n\
\n\
     -q    : suppress auxiliary output\n\
\n\
  See program text for much more information.\n"

/**************************************************************
Counts for n=1..45:

 1: 1
 2: 1
 3: 1
 4: 1
 5: 2
 6: 3
 7: 6
 8: 11
 9: 23
10: 47
11: 106
12: 235
13: 551
14: 1301
15: 3159
16: 7741
17: 19320
18: 48629
19: 123867
20: 317955
21: 823065
22: 2144505
23: 5623756
24: 14828074
25: 39299897
26: 104636890
27: 279793450
28: 751065460
29: 2023443032
30: 5469566585
31: 14830871802
32: 40330829030
33: 109972410221
34: 300628862480
35: 823779631721
36: 2262366343746
37: 6226306037178
38: 17169677490714
39: 47436313524262
40: 131290543779126
41: 363990257783343
42: 1010748076717151
43: 2810986483493475
44: 7828986221515605
45: 21835027912963086
********************************/

/* Comments on original program by original authors */
/*==============================================================*/
/* program: freetree.c                                          */
/* purpose: generating all free trees                           */
/* input  : n -- number of nodes                                */
/*          m -- max degree                                     */
/*          lb,ub -- lower and upper bound on diameter          */
/*          res/mod -- splitting into parts                     */
/* output : listing of free trees in relex order                */ 
/* date   : September 1995, updated Sept 2000                   */
/* programmers: Gang Li & Frank Ruskey                          */
/* algorithm: From the paper: G. Li and F. Ruskey,  "The        */
/*    Advantages of Forward Thinking in Generating Rooted and   */
/*    Free Trees",  10th Annual ACM-SIAM Symposium on Discrete  */
/*    Algorithms (SODA), (1999) S939-940.  See the web page at  */
/*    http://www.theory.csc.UVic.CA/~fruskey/                   */
/*      Publications/RootedFreeTree.html                        */
/* more info: See                                               */
/*    http://www.theory.csc.UVic.CA/~cos/inf/FreeTrees.html     */
/*==============================================================*/

#define MAXN  128  /* max size of the tree.
                      Check MAXOUTLEN if more than 1000 */
#include "gtools.h"

#define NAUTY_PGM  4   /* 1 = geng, 2 = genbg, 3 = gentourng, 4 = gentreeg */

/*****************************************************************

OUTPROC feature.

   By defining the C preprocessor variable OUTPROC at compile time
   (for Unix the syntax is  -DOUTPROC=procname  on the cc command),
   gentreeg can be made to call
       void OUTPROC(FILE *f, int *par, int n)
   instead of the usual output routines.

   par[1..n] is the parent array of the tree.
   The edges in 1..n labelling are {j,par[j]} for j=2..n.
   WARNING: The vertices are 1..n, not 0..n-1, and also par[]
   uses elements 1..n, not 0..n-1.

   Your procedure can be in a separate file so long as it is linked
   with gentreeg. 

PRUNE feature.

   By defining the C preprocessor variable PRUNE at compile time,
   gentreeg can be made to call
        int PRUNE(int *par, int n)
   for each tree.  The meaning of par[] is above and n is the number
   of vertices.  If a non-zero value is returned, the tree is rejected:
   it is not included in the count, and output routines including
   OUTPROC are not called.

SUMMARY

   If the C preprocessor variable SUMMARY is defined at compile
   time, the procedure/macro SUMMARY(nauty_counter nout, double cpu)
   is called just before the program exits. The purpose is to allow
   reporting of statistics collected by PRUNE or OUTPROC. The values
   nout and cpu are the output count and cpu time as reported on the
   >Z line. Output should be written to stderr.

CALLING FROM A PROGRAM

   It is possible to call gentreeg from another program instead of
   using it as a stand-alone program. The main requirement is to
   change the name of the main program to be other than "main". This
   is done by defining the preprocessor variable GENTREEG_MAIN. You
   might also like to define OUTPROC to be the name of a procedure
   to receive the graphs. To call the program you need to define an
   argument list argv[] consistent with the usual one; don't forget
   that argv[0] is the command name and not the first argument. The
   value of argc is the number of strings in argv[]; that is, one
   more than the number of arguments.

EFFICIENCY

   Note that the trees are generated much faster than they can be
   written.  A substantial time saving can be achieved using one
   of the above methods instead of writing the trees and then
   reading them into another program.
*/

static int
 par[MAXN+1],                 /* parent position of i  */
 maxchild,                    /* max number of children */
 chi[MAXN+1],                 /* number of children of a node */
 nextp[MAXN+1],               /* next good pos to add nodes */
 rChi[MAXN+1],                /* the right most child of node i */
 ub;                          /* upper bound on something */

static nauty_counter nout;                 /* number of trees */
static FILE *outfile;

static int splitlevel,splitcount,mod,res;   /* -s res/mod */
static boolean irred; /* homeomorphically irreducible */

static int nv;  /* number of vertices */
static int mindiam; 
static int maxdeg; 
static int maxdiam;

/* MAXOUTLEN must be at least the longest output line length */
#define MAXOUTLEN (10 + 4*MAXN)
static char outstring[MAXOUTLEN];

void (*outproc)(FILE *f, int vpar[], int n);
#ifdef OUTPROC
extern void OUTPROC(FILE *f, int vpar[], int n);
#endif
#ifdef PRUNE
extern int PRUNE(int vpar[], int n);
#endif

void WritePar(FILE *f, int vpar[], int n)
/* Write parent array */
{
   int i,j;
   char *pout,*p,one[8];
   size_t len;

   pout = outstring;

   for (i=1;i<=nv;++i) {
       j = vpar[i];
       if (j == 0) *(pout++) = '0';
       else {
          for (p = one; j > 0; j /= 10)
              *(p++) = '0' + j%10;
          while (--p >= one) *(pout++) = *p;
       }
       if (i < nv) *(pout++) = ' ';
       else       *(pout++) = '\n';
   }
   len = pout - outstring;

   if (fwrite(outstring,sizeof(char),len,f) != len)
        gt_abort(">E gentreeg: fwrite() failed\n");
}

static void
WriteLev(FILE *f, int vpar[], int n)
/* Write levels array */
{
   int i,j;
   char *pout,*p,one[8];
   size_t len;
   int lev[MAXN+1];

   lev[1]=0;
   for ( i=2; i<=nv; ++i ) lev[i] = lev[vpar[i]]+1;

   pout = outstring;

   for (i=1;i<=nv;++i) {
       j = lev[i];
       if (j == 0) *(pout++) = '0';
       else {
          for (p = one; j > 0; j /= 10)
              *(p++) = '0' + j%10;
          while (--p >= one) *(pout++) = *p;
       }
       if (i < nv) *(pout++) = ' ';
       else       *(pout++) = '\n';
   }
   len = pout - outstring;

   if (fwrite(outstring,sizeof(char),len,f) != len)
        gt_abort(">E gentreeg: fwrite() failed\n");
}

static void
DontWrite(FILE *f, int vpar[], int n)
/* Null print routine */
{
}

static void
WriteS6(FILE *f, int vpar[], int n)
/* Write in sparse6 format */
{
    char *pout;
    int nb,i,j,lastj,x,k,r,rr,topbit;
    size_t len;

    pout = outstring;
    *pout++ = ':';
    encodegraphsize(n,&pout);

    for (i = n-1, nb = 0; i != 0 ; i >>= 1, ++nb) {}
    topbit = 1 << (nb-1);
    k = 6;
    x = 0;

    lastj = 0;
    for (j = 1; j < n; ++j)
    {
        i = vpar[j+1] - 1;
        if (j == lastj)
        {
            x <<= 1;
            if (--k == 0)
            {
                *pout++ = 63 + x;
                k = 6;
                x = 0;
            }
        }
        else
        {
            x = (x << 1) | 1;
            if (--k == 0)
            {
                *pout++ = 63 + x;
                k = 6;
                x = 0;
            }
            if (j > lastj+1)
            {
                for (r = 0, rr = j; r < nb; ++r, rr <<= 1)
                {
                    if (rr & topbit) x = (x << 1) | 1;
                    else             x <<= 1;
                    if (--k == 0)
                    {
                        *pout++ = 63 + x;
                        k = 6;
                        x = 0;
                    }
                }
                x <<= 1;
                if (--k == 0)
                {
                    *pout++ = 63 + x;
                    k = 6;
                    x = 0;
                }
            }
            lastj = j;
        }
        for (r = 0, rr = i; r < nb; ++r, rr <<= 1)
        {
            if (rr & topbit) x = (x << 1) | 1;
            else             x <<= 1;
            if (--k == 0)
            {
                *pout++ = 63 + x;
                k = 6;
                x = 0;
            }
        }
    }

    if (k != 6) *pout++ = 63 + ((x << k) | ((1 << k) - 1));

    *pout++ = '\n';
    len = pout - outstring;

    if (fwrite(outstring,sizeof(char),len,f) != len)
        gt_abort(">E gentreeg: fwrite() failed\n");
}

static boolean
ishi(int *par, int n)
/* Test if it has a vertex of degree 2 */
{
    int degm1[MAXN+1];  /* Degrees minus 1 */
    int i;

    degm1[1] = -1;
    for (i = 2; i <= n; ++i) degm1[i] = 0;
    for (i = 2; i <= n; ++i) ++degm1[par[i]];
    for (i = 1; i <= n; ++i) if (degm1[i] == 1) return TRUE;

    return FALSE;
}

static void
WriteIt(int level)
{
   if (level < splitlevel && res != 0) return;

   if (irred && ishi(par,nv)) return;
#ifdef PRUNE
   if (PRUNE(par,nv)) return;
#endif

   ++nout;
   (*outproc)(outfile,par,nv);
}

static boolean
good( int p, int h, int t ) {
  if (p==2 && mindiam<=2 && t==0) return TRUE;
  if (t == 1) {
     if (2*h>=mindiam+1 && 2*h <= maxdiam+1) {
        if ((p-1)*2 >= nv) return TRUE;  
        else if (p - h-1 == 1) {
          if (par[p]> 2) return TRUE;
        } else
          if ((p - h-1 >=2) && ((par[h+2]>2) || (par[h+3]>2))) return TRUE;
     }
  } else 
  if (nv-p >= h && 2*h>=mindiam) { 
     if (maxdiam==nv-1 && nv%2==0) return 2*h<=maxdiam+1; 
     else return 2*h <= maxdiam;
  }
  return FALSE;
} /* good */   


static void
Gen( int level, int p, int s, int cL, int h, int l, int n, int f, int g )
/* The main generation procedure. */
{
  int hh,flag,entry,temp;

  if (level == splitlevel)
  {
      if (splitcount-- != 0) return;
      splitcount = mod - 1;
  }

  if (p > n) 
     if (f == 0) {
        if (good(p-1,h,0)) Gen(level+1, p, 2, p-2,h,n,nv,1,0); 
        if (good(p-1,h,1)) Gen(level+1, p, 2, p-3,h,n,nv,1,1);
     } else { WriteIt(level);}
  else {
     if (cL == 0) {
        if ( p< ub+2 ) par[p] = p-1;
        else {
           Gen(level+1, p, p-1, 1, h, l, n, f, g );
           return; 
        }
     } else  
        if (par[p-cL] < s) par[p]=par[p-cL];       
        else {
           par[p] = cL + par[p-cL];
           if (g==1) 
           {
               if (((l-1)*2 < n) && (p-cL<=l) && (
                   ((p-cL+1<l) &&  (par[p-cL+1]==2)  
                   && (p-cL+2<=l) && (par[p-cL+2]==2))     /*case 1*/
                   || ((p-cL+1==l) && (par[p-cL+1]==2))    /*case 2*/
                   || (p-cL+1>l))) {                       /*case 3*/
                  s= par[p]; cL= p-s;
                  par[p] = par[par[p]];
               } else 
                  if (par[p-cL]==2) par[p]=1;  
           }
        }
        if (s!=0 || p<=ub+1) {
           chi[par[p]] = chi[par[p]] + 1;
           temp = rChi[par[p]]; rChi[par[p]] = p;
           if (chi[par[p]] <= ((par[p]==1)?maxchild:maxchild-1)) {
              if (chi[par[p]] < (par[p]==1?maxchild:maxchild-1)) nextp[p] = par[p];
              else nextp[p] = nextp[par[p]];
              Gen(level+1, p+1, s, cL,h,l,n,f,g ); 
           }
           chi[par[p]] = chi[par[p]] - 1;
           rChi[par[p]] = temp;
        }
        if (s==0 && 2*(p-2)<mindiam) return;

        nextp[p] = nextp[par[p]];
        entry = nextp[p];
        flag= 0; hh=1;
        while ((((f==0) && (entry>=2)) ||
               ((f==1) && (entry>=1))) && (flag==0)) {
           if (s==0) h = p-2;
           if (p<=l+h-g) hh = 0; 
           if ((f==0) || (hh==1)) {
              /*s = par[p]; par[p] = par[s];*/
              par[p] = entry;

              chi[entry] = chi[entry] + 1;
              temp = rChi[par[p]];  rChi[par[p]] = p;
              if (chi[entry] >= (entry==1?maxchild:maxchild-1)) nextp[p] = nextp[entry];
              if (f == 0) Gen(level+1,p+1,temp,p-temp,h,0,nv-h+1,f,g);
              else if (hh == 1) Gen(level+1,p+1,temp,p-temp,h,l,n,f,g);
              chi[entry] = chi[entry] - 1;
              rChi[par[p]] = temp;
              entry = nextp[entry];
              nextp[p] = entry;
           } else flag= 1;
        }
        if (f == 0) {
           if (good(p-1,h,0)) Gen(level+1,p,2,p-2,h,p-1,nv,1,0);
           if (good(p-1,h,1)) Gen(level+1,p,2,p-3,h,p-1,nv,1,1);
        }
  }
} /* Gen */

/**************************************************************************/
/**************************************************************************/

#ifdef GENTREEG_MAIN
int
GENTREEG_MAIN(int argc, char *argv[])
#else
int
main(int argc, char *argv[])
#endif
{
    char *arg;
    boolean badargs,uswitch,sswitch,pswitch,lswitch,Zswitch;
    boolean qswitch,Dswitch,gotmr,gotf;
    long Z1,Z2,Z1arg,Z2arg;
    char *outfilename,sw;
    int i,j,argnum;
    int splitlevinc;
    int minnv,maxnv;
    double t1,t2;
    int maxdegarg;
    char msg[201];

    HELP; PUTVERSION;

    badargs = FALSE;
    uswitch = sswitch = pswitch = lswitch = Zswitch = FALSE;
    qswitch = Dswitch = gotmr = gotf = irred = FALSE;
    outfilename = NULL;

    splitlevinc = 0;
    
    argnum = 0;
    for (j = 1; !badargs && j < argc; ++j)
    {
        arg = argv[j];
        if (arg[0] == '-' && arg[1] != '\0')
        {
            ++arg;
            while (*arg != '\0')
            {
                sw = *arg++;
                     SWBOOLEAN('u',uswitch)
                else SWBOOLEAN('s',sswitch)
                else SWBOOLEAN('p',pswitch)
                else SWBOOLEAN('l',lswitch)
                else SWBOOLEAN('i',irred)
                else SWRANGE('Z',":-",Zswitch,Z1arg,Z2arg,"gentreeg -Z")
                else SWBOOLEAN('q',qswitch)
                else SWINT('D',Dswitch,maxdegarg,"gentreeg -D")
#ifdef PLUGIN_SWITCHES
PLUGIN_SWITCHES
#endif
                else badargs = TRUE;
            }
        }
        else if (arg[0] == '-' && arg[1] == '\0')
            gotf = TRUE;
        else
        {
            if (argnum == 0)
            {
                if (sscanf(arg,"%d:%d",&minnv,&maxnv) == 2)
                {}
                else if (sscanf(arg,"%d",&minnv) != 1)
                    badargs = TRUE;
                else
                    maxnv = minnv;
                ++argnum;
            }
            else if (gotf)
                badargs = TRUE;
            else
            {
                if (!gotmr)
                {
                    if (sscanf(arg,"%d/%d",&res,&mod) == 2)
                    { 
                        gotmr = TRUE; 
                        continue; 
                    }
                }
                if (!gotf)
                {
                    outfilename = arg;
                    gotf = TRUE;
                    continue;
                }
            }
        }
    }

    gtools_check(WORDSIZE,SETWORDSNEEDED(maxnv),maxnv,NAUTYVERSIONID);

    if (argnum == 0)
        badargs = TRUE;
    else if (minnv < 1 || maxnv > MAXN || minnv > maxnv)
    {
        fprintf(stderr,">E gentreeg: n must be in the range 1..%d\n",MAXN);
        exit(1);
    }

    if (badargs)
    {
        fprintf(stderr,">E Usage: %s\n",USAGE);
        GETHELP;
        exit(1);
    }

    if (!gotmr)
    {
        mod = 1;
        res = 0;
    }

    if (!Dswitch || maxdegarg >= maxnv) maxdegarg = maxnv - 1;
    if (!Zswitch || Z2arg >= maxnv) Z2arg = maxnv - 1;
    if (!Zswitch || Z1arg <= 0) Z1arg = (minnv == 1 ? 0 : minnv == 2 ? 1 : 2);

    if (Z2arg < Z1arg || Z2arg < (minnv == 1 ? 0 : minnv == 2 ? 1 : 2))
        gt_abort(">E gentreeg: impossible diameter bounds\n");
    if (maxdegarg < (maxnv == 1 ? 0 : maxnv == 2 ? 1 : 2))
        gt_abort(">E gentreeg: impossible degree bound\n");

    if (res < 0 || res >= mod)
        gt_abort(">E gentreeg: must have 0 <= res < mod\n");

    if ((lswitch!=0) + (pswitch!=0) + (sswitch!=0) + (uswitch!=0) > 1)
        gt_abort(">E gentreeg: -uslp are incompatible\n");

#ifdef OUTPROC
    outproc = OUTPROC;
#else
    if      (uswitch)  outproc = DontWrite;
    else if (lswitch)  outproc = WriteLev;
    else if (pswitch)  outproc = WritePar;
    else               outproc = WriteS6;
#endif

#ifdef PLUGIN_INIT
PLUGIN_INIT
#endif

    if (qswitch)
        outfile = stdout;
    else if (!gotf || outfilename == NULL)
    {
        outfilename = "stdout";
        outfile = stdout;
    }
    else if ((outfile = fopen(outfilename,"w")) == NULL)
    {
        fprintf(stderr,
              ">E gentree: can't open %s for writing\n",outfilename);
        gt_abort(NULL);
    }

    if (!qswitch)
    {
        msg[0] = '\0';
        if (strlen(argv[0]) > 75)
            fprintf(stderr,">A %s",argv[0]);
        else
            CATMSG1(">A %s",argv[0]);
        if (irred) CATMSG0(" -i");
       
        CATMSG2(" Z=%d:%d",mindiam,maxdiam);
        if (minnv == maxnv)
            CATMSG2(" D=%d n=%d",maxdeg,nv);
        else
            CATMSG3(" D=%d n=%d:%d",maxdeg,minnv,maxnv);
        if (mod > 1) CATMSG2(" class=%d/%d",res,mod);
        CATMSG0("\n");
        fputs(msg,stderr);
        fflush(stderr);
    }

    t1 = CPUTIME;

    nout = 0;

    for (nv = minnv; nv <= maxnv; ++nv)
    {
        maxdeg = maxdegarg;
        if (!Dswitch || maxdeg >= nv) maxdeg = nv - 1;
     
        Z1 = Z1arg; Z2 = Z2arg;
        if (!Zswitch || Z2 >= nv) Z2 = nv - 1;
        if (!Zswitch || Z1 < 0) Z1 = (nv == 1 ? 0 : nv == 2 ? 1 : 2);
        if (Z1 < 1 && nv == 2) Z1 = 1;
        if (Z1 < 2 && nv > 2) Z1 = 2;
        mindiam = (int)Z1;
        maxdiam = (int)Z2;

        if (nv == 1)
        {
            if (res == 0 && mindiam <= 0)
            {
                par[1] = 0;
                WriteIt(0);
            }
        }
        else if (nv == 2)
        {
            if (res == 0 && maxdeg >= 1 && mindiam <= 1 && maxdiam >= 2)
            {
                par[1] = 0;
                par[2] = 1;
                WriteIt(0);
            }
        }
        else
        {
            if (mod > 1 && nv > 3)
            {
                splitlevel = nv/2;
                splitcount = res;
            }
            else
            {
                splitlevel = -1;
                mod = 1;
            }
            if (maxdeg > 0) maxchild = maxdeg; 
            else maxchild = -1;

            ub = (maxdiam+1)/2;
            for(i=1;i<=nv;i++) chi[i]=0;

            par[1] = 0; par[2] = 1;
            nextp[1]=0; nextp[2]=1;chi[1]=1;
            Gen( 0, 3, 0, 0, ub, 0, nv-ub+1, 0, 0 ); 
        }
    }

    t2 = CPUTIME;

#ifdef SUMMARY
    SUMMARY(nout,t2-t1);
#endif

    if (!qswitch)
    {
        fprintf(stderr,">Z " COUNTER_FMT 
                " trees generated in %3.2f sec\n",nout,t2-t1);
    }

#ifdef GENTREEG_MAIN
    return 0;
#else
    exit(0);
#endif
}
