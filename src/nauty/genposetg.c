/* poset generator version 1.0, October 8, 2022. */

/* Compile: cc -O3 -march=native genposetg.c nautyS1.a
   If -DCOUNT is added, only counting is allowed and writing is disabled.

   With the t switch, write topologically labelled posets to stdout
   in digraph6 format.  All edges x-->y have x<y.
   
   With the o switch, write posets in digraph6 format without
   topologically sorting. This is faster but maybe the outputs will
   be harder for you to process.

   This version is limited to 16 vertices and some expertise is needed
   to make it work for larger sizes. You can't just increase WORDSIZE!
   However, since it already takes 5 years or so just to count up to
   16 vertices, you might consider whether the next size (about 100 times
   as hard) is plausible.

   Counts:
       1  1
       2  2
       3  5
       4  16
       5  63
       6  318
       7  2045
       8  16999
       9  183231
      10  2567284
      11  46749427
      12  1104891746
      13  33823827452
      14  1338193159771
      15  68275077901156
      16  4483130665195087

   The next two values are unknown, but are probably close to 4e17 and 4e19.

   PLUGIN facility. You can specify a file to be included with code to
   restrict the output. Create a file, say "posetplugin.c" and specify it
   on the genposetg compilation command like this:  -DPLUGIN='"posetplugin.c"'
   (including all the quotes). This file can contain
   1. Declaration of static variables and procedures.
   2. Optional macros POSET_PRUNE0, POSET_PRUNE1, and POSET_SUMMARY .
    POSET_PRUNE0(pos,n) is called in the output procedure for each poset.
        The poset is pos[0..n-1] in internal labelling.
    POSET_PRUNE1(pos,n) is the same except now the poset is in topological order.
        In each case you can discard the poset by executing "return".
        Note that POSET_PRUNE1 is only invoked if you use the 't' switch.
        You can use both POSET_PRUNE0 and POSET_PRUNE1 if you want.
    POSET_SUMMARY is placed just before the final summary line and is intended
        for writing statistics that the other macros have collected.

   For example, this plugin will restrict outputs to those posets with
   exactly two sources and two sinks.
-----------------------
#define POSET_PRUNE0(pos,n) { int i,sinks=0; graph gin=0; \
    for (i = 0; i < n; ++i) { if (pos[i] == 0) ++sinks; gin |= pos[i]; } \
    if (sinks != 2 || n - POPCOUNT(gin) != 2) return; }

#define POSET_SUMMARY fprintf(stderr,"That does it.\n")
-----------------------

   Authors:  Gunnar Brinkmann and Brendan McKay.
*/

#if defined(WORDSIZE) && WORDSIZE != 16
#error "Only WORDSIZE=16 is supported for this program."
#endif

#define WORDSIZE 16
#define MAXN WORDSIZE
#define MAXGRAPH ((1<<WORDSIZE)-1)
#define MAXCHECK 6 /* Konstante, ab wann versucht wird, Nachbarschaften zu faerben, um .. */
#define CHECKSIZE (MAXCHECK>MAXN ? MAXCHECK : MAXN)

#include "gtools.h"
#include <ctype.h>

/*-------------------------------------------*/
/* fuer nauty: */
static int generators[MAXN][MAXN];
static int number_of_generators;
static setword transposition[MAXN];   /* If transposition[i] is non-zero,
       the generator is a transposition and the two moved elements
       are given in transposition[i].  If it is zero, the generator
       may or may not be a transposition and the generator is given
       in full in generator[i].  PERMSET(x,i,y) sets y to the image
       of x under transposition/generators[i]).  x and y are setwords,
       not pointers. */
#define PERMSET(x,i,y) \
  if (transposition[i]) {setword temp = (x)&transposition[i]; \
    if (temp && temp!=transposition[i]) y=(x)^transposition[i]; else y=(x);}\
  else permset(&x,&y,1,generators[i]);


static int lab[MAXN][MAXN], ptn[MAXN][MAXN], orbits[MAXN];
/* lab[i] und ptn[i] beschreiben die durch die Automorphismen
   vorgegebene Faerbung der Knoten soweit sie bisher (d.h. auf
   level[i] festgelegt werden kann.  Sie wird immer neu berechnet,
   sobald ein neues level angefangen wird. Die alte wird dann
   zwischengespeichert. Die Knoten auf dem aktuell bearbeiteten Level
   bekommen immer die gleiche Farbe. */
static DEFAULTOPTIONS(options);
static DEFAULTOPTIONS(options_canon);
static DEFAULTOPTIONS(options_final);
static statsblk stats ;
static setword workspace[100*MAXN];
/*-------------------------------------------*/

static long long int nautycalls_canon=0, nautycalls_group=0;
static long long int anzahl_posets=0, written_posets=0;
static long long int punkte_letztes_level[MAXN+1]={0};
static long long int anzahl_level[MAXN+1]={0};
static long long int nauty_count1=0,nauty_count2=0,nauty_count3=0,nauty_count4=0,nauty_count5=0;
static long long int minmax_counter[MAXN+1][MAXN+1]={{0}};
static long long int *minmaxline;
static long long int relations[((MAXN*(MAXN-1)))/2+1]={0};

static int punktzahl, maxpunktzahl, mod=0, rest=0, splitpunktzahl=0, splitzaehler= -1;
static int topsort;  /* 0 for internal order, 1 for topological order */

/* die Levelnumerierung startet bei 0 -- das Ende von Listen wird durch eine 0 markiert. */
static graph subsets[MAXN][MAXGRAPH+1]; /* subsets[i]: Die Liste der direkten oberen Nachbarn auf level i */
static int anzahl_subsets[MAXN]; /* die Laenge dieser Liste, d.h. die Anzahl der Moeglichkeiten
                             fuer DIREKTE obere Nachbarn*/

static graph nbs[MAXN][MAXGRAPH+1]; /* nbs[i-1]: Die Liste der moeglichen Nachbarschaften fuer Knoten auf
                              level i */
static int anzahl_nbs[MAXN]; /* die Laenge dieser Liste */

static graph all_ancestors[MAXGRAPH+1]={0}; /* all_ancestors[i] sind alle Vorgaenger der Menge "i" als
                                      binaere Menge (oder Graph) interpretiert */

static graph nonequiv_nbs[MAXN][MAXGRAPH+1]; 
                           /* hier werden die wirklich verschiedenen nbs gespeichert.
                              Haengt nicht vom level ab, sondern von der Punktzahl */
static int maxorbsize_level[MAXN]; /* maxorbsize_level[i]: Wir gross ist der groesste Orbit der 
                               Automorphismengruppe fuer das Poset eingeschraenkt auf level 
                               0...i auf den Nachbarschaften ? */
static graph orbitnumber[MAXN][MAXGRAPH]; /* hier werden die Orbitindizees gespeichert, wenn ein
                                      neues level begonnen wird */
static int first_point_level[MAXN]; /* Die Nummer des kleinsten Knotens auf dem Level */

static graph *jump[MAXN][MAXN+1]; /* ist i auf level j, so ist jump[i][j] die erste nbs, die
                              j nicht als kleinstes Element hat */

#define MAXMARK (INT_MAX-2)
static int markvalue_ = MAXMARK;
static int marks_[MAXGRAPH+1];
#define RESETMARKS {int mki_; if ((++markvalue_) >= MAXMARK) \
       { markvalue_ = 1; for (mki_=0;mki_<=MAXGRAPH;++mki_) marks_[mki_]=0;}}
#define UNMARK(x) (marks_[x] = 0)
#define ISMARKED(x) (marks_[x] == markvalue_)
#define UNMARKED(x) (marks_[x] != markvalue_)
#define MARK(x) (marks_[x] = markvalue_)

#define FOREACH(i,set) for ((i)=elementlist[set]; *(i)!=MAXN; (i)++)

static int points_on_level[MAXN]={0}; /* Die Anzahl der Knoten auf dem Level */
static graph poset[MAXN]={0};
static graph undirected_poset[MAXN]={0};
static graph all_points_up_to[MAXN]={0}; /* all_points_up_to[i] ist die Bitmap, die alle Knoten 
                                     bis level i enthaelt */
static graph all_points_on[MAXN]={0}; /* all_points_on[i] ist die Bitmap, die alle Knoten 
                                  auf level i enthaelt */


/* BDM char next_element[MAXGRAPH+1][MAXN+1]; * die Funktion nextelement() als Vektor */
static char elementlist[MAXGRAPH+1][MAXN+1]; /* Die Liste der Elemente */
/* BDM char first_bit[MAXGRAPH+1];  * das Makro FIRSTBIT als Vektor */

static int output=0;

#ifdef PLUGIN
#include PLUGIN
#endif

/**********************************************************************************/

static void writeposet(graph *g, int n)
/* Write a poset to stdout, optionally in upper-triangular form */
{
    graph h[MAXN];
    int p[MAXN],pinv[MAXN];
    int i,j,k;
    setword left,new,w;

#ifdef POSET_PRUNE0
    POSET_PRUNE0(g,n);
#endif

    if (!topsort)
    {
        writed6(stdout,g,1,n);
#if defined(POSET_PRUNE0) || defined(POSET_PRUNE1)
        written_posets++;
#endif
        return;
    }

    for (i = 0; i < n; ++i) h[i] = g[i];

    left = ALLMASK(n);
    j = 0;
    while (left)
    {
        w = 0;
        for (i = 0; i < n; ++i) w |= h[i];
        new = left & ~w;
        left &= ~new;
        while (new)
        {
            TAKEBIT(k,new);
            p[j++] = k;
            h[k] = 0;
        }
    }
    
    if (j != n) gt_abort(">E topsort error\n");

    for (i = 0; i < n; ++i) pinv[p[i]] = i;

    for (i = 0; i < n; ++i)
    {
        left = g[p[i]];
        w = 0;
        while (left)
        {
            TAKEBIT(k,left);
            w |= bit[pinv[k]];
        }
        h[i] = w;
    }

#ifdef POSET_PRUNE1
    POSET_PRUNE1(h,n);
#endif
#if defined(POSET_PRUNE0) || defined(POSET_PRUNE1)
    written_posets++;
#endif

    writed6(stdout,h,1,n);
}

/**********************************************************************************/

static void sammle_permutationen(int count, int perm[], int orbits[],
                          int numorbits, int stabvertex, int n)
{
  memcpy(generators+number_of_generators,perm,sizeof(int)*n);
  transposition[number_of_generators] = 0;

  number_of_generators++;
}
/**********************************************************************************/

static void init_nauty()
/* initialises the nauty variables */
{
 /* tc_level = 0 is not the default in the most recent
    editions of nauty.  However, it is better for very small
    graphs so we set it here. */

options.getcanon=FALSE;
options.userautomproc = sammle_permutationen;
options.defaultptn=FALSE;
options.writeautoms=FALSE;
options.writemarkers=FALSE;
options.digraph=FALSE;
options.tc_level = 0;

options_canon.getcanon=TRUE;
options_canon.userautomproc = sammle_permutationen;
options_canon.defaultptn=FALSE;
options_canon.writeautoms=FALSE;
options_canon.writemarkers=FALSE;
options_canon.digraph=FALSE;
options_canon.tc_level = 0;

options_final.getcanon=TRUE;
options_final.defaultptn=FALSE;
options_final.writeautoms=FALSE;
options_final.writemarkers=FALSE;
options_final.digraph=FALSE;
options_final.tc_level = 0;
}


/**********************************************************************************/


static boolean refinex(graph *g, int *lab, int *ptn, int *numcells,
                int *count, set *active, setword good_vertices, int n)
/*  custom version of refine1() which can exit quickly if required.

  If good_vertices is not 0, a FALSE return is made if it can be determined
  that the last cell containing a member of good_vertices also contains n-1.
  In that case (FALSE return) the refinement process may be incomplete.

  In all other cases, the refinement is completed and TRUE is returned. */
{
        register int i,c1,c2,labc1;
        register setword x,lact;
        int split1,split2,cell1,cell2;
        int cnt,bmin,bmax;
        set *gptr;
        setword workset;
        int workperm[MAXN];
        int bucket[MAXN+2];

        if (n == 1)  return TRUE;

        lact = *active;

        split1 = -1;
        while (*numcells < n && lact)
        {
            x = lact & (-lact);
            split1 = FIRSTBIT(x);
            lact ^= x;
            
            for (split2 = split1; ptn[split2] > 0; ++split2)
            {}
            if (split1 == split2)       /* trivial splitting cell */
            {
                gptr = GRAPHROW(g,lab[split1],1);
                for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
                {
                    for (cell2 = cell1; ptn[cell2] > 0; ++cell2)
                    {}
                    if (cell1 == cell2)
                        continue;
                    c1 = cell1;
                    c2 = cell2;
                    while (c1 <= c2)
                    {
                        labc1 = lab[c1];
                        if (ISELEMENT1(gptr,labc1))
                            ++c1;
                        else
                        {
                            lab[c1] = lab[c2];
                            lab[c2] = labc1;
                            --c2;
                        }
                    }
                    if (c2 >= cell1 && c1 <= cell2)
                    {
                        ptn[c2] = 0;
                        ++*numcells;
                        lact |= bit[c1];
                    }
                }
            }

            else        /* nontrivial splitting cell */
            {
                workset = 0;
                for (i = split1; i <= split2; ++i)
                    workset |= bit[lab[i]];

                for (cell1 = 0; cell1 < n; cell1 = cell2 + 1)
                {
                    for (cell2 = cell1; ptn[cell2] > 0; ++cell2)
                    {}
                    if (cell1 == cell2)
                        continue;
                    i = cell1;
                    if ((x = workset & g[lab[i]]))     /* not == */
                        cnt = POPCOUNT(x);
                    else
                        cnt = 0;
                    count[i] = bmin = bmax = cnt;
                    bucket[cnt] = 1;
                    while (++i <= cell2)
                    {
                        if ((x = workset & g[lab[i]])) /* not == */
                            cnt = POPCOUNT(x);
                        else
                            cnt = 0;
                        while (bmin > cnt)
                            bucket[--bmin] = 0;
                        while (bmax < cnt)
                            bucket[++bmax] = 0;
                        ++bucket[cnt];
                        count[i] = cnt;
                    }
                    if (bmin == bmax)
                    {
                        continue;
                    }
                    c1 = cell1;
                    for (i = bmin; i <= bmax; ++i)
                        if (bucket[i])
                        {
                            c2 = c1 + bucket[i];
                            bucket[i] = c1;
                            if (c1 != cell1)
                            {
                                lact |= bit[c1];
                                ++*numcells;
                            }
                            if (c2 <= cell2)
                                ptn[c2-1] = 0;
                            c1 = c2;
                        }
                    for (i = cell1; i <= cell2; ++i)
                        workperm[bucket[count[i]]++] = lab[i];
                    for (i = cell1; i <= cell2; ++i)
                        lab[i] = workperm[i];
                }
            }

            if (good_vertices)
            {
                for (i = n-1; (bit[lab[i]]&good_vertices) == 0 ; --i) {}
        
                while (ptn[i] > 0) ++i;

                for (; lab[i] != n-1;)
                {
                    --i;
                    if (i < 0 || ptn[i] == 0) return FALSE;
                }
            }
        }

        return TRUE;
}


/*******************************************************************************
***/

static boolean call_nauty(int *lab, int *ptn, int n, int *orbits,
                                   setword good_vertices, optionblk *options)
/* Call nauty on undirected_poset, maybe just pretending.  On exit,
     stats.numorbits gives the number of orbits.  Other fields of stats
     are not necessarily set.
   Only n <= WORDSIZE is supported!

   If options->getcanon==TRUE:
       If the canonically last orbit containing v in good_vertices (must exist!)
          also contains vertex n-1, then lab[] will be in canonical order and
          the value TRUE is returned.  Also orbits[] is correct (if not NULL)
          and generators are stored in transposition[]/generators[] if
          options->userautomproc!=NULL.
       Else the value FALSE is returned.  In this case no other information
          about the group or labelling can be assumed correct.
   If options->getcanon==FALSE:
       good_vertices is ignored.  lab[] is not necessarily in canonical order.
       The value TRUE is returned.  Also orbits[] is correct (if not NULL)
       and generators are stored in transposition[]/generators[] if
       options->userautomproc!=NULL.
*/
{
    setword active;
    int tmp_perm[MAXN];
    int tmp_orbits[MAXN],*orbs;
    int i,j,tpos,tpos2,numcells;
    int a,b,c,d,which,min;
    graph canong[MAXN];
    boolean cells_are_orbits;

    /* First call refinex().  This will perform a refinement of the 
    partition while at the same time looking for the possibility of
    a quick reject or accept. */

    active = 0;
    numcells = 0;
    for (i = 0; i < n; ++i)
    {
        ++numcells;
        active |= bit[i];
        while (ptn[i]) ++i;
    }

    if (!refinex(undirected_poset,lab,ptn,&numcells,tmp_perm,&active,
        (options->getcanon?good_vertices:0),n)) return FALSE;

    /* The partition returned by refine1() is known from theory to be the
    orbits partition in many cases.  nauty contains a procedure for this.
    It includes the case where there are n, n-1 or n-2 cells and in those
    cases we can easily find the group as well as just the orbits. */

    cells_are_orbits = cheapautom(ptn,0,FALSE,n);

    if (numcells == n)
    {
        ++nauty_count1;
        stats.numorbits = n;
        if (orbits) for (i = 0; i < n; ++i) orbits[i] = i;

        /* lab is already in canonical order */

        return TRUE;
    }
    else if (numcells == n-1)
    {
        /* This is a transposition (ab). */

        ++nauty_count2;
        for (tpos = 0; ptn[tpos] == 0; ++tpos) {}
        a = lab[tpos]; b = lab[tpos+1];
        stats.numorbits = n-1;

        if (options->userautomproc)   /* Make the generator (ab). */
            transposition[number_of_generators++] = bit[a]|bit[b];

        if (orbits)
        {
            for (i = 0; i < n; ++i) orbits[i] = i;
            if (a > b) orbits[a] = b; else orbits[b] = a;
        }

        /* lab is already in canonical order */

        return TRUE;
    }
    else if (numcells == n-2)
    {
        /* There are three cases here.  
           which=1 : (ab) and (bc)
           which=2 : (ab) and (cd)
           which=3 : (ab)(cd) */

        ++nauty_count3;
        stats.numorbits = n-2;

        for (tpos = 0; ptn[tpos] == 0; ++tpos) {}
        a = lab[tpos]; b = lab[tpos+1];
        if (ptn[tpos+1] > 0)
        {
            which = 1;
            c = lab[tpos+2];
        }
        else
        {
            for (tpos2 = tpos+2; ptn[tpos2] == 0; ++tpos2) {}
            c = lab[tpos2]; d = lab[tpos2+1];
            if ((undirected_poset[a]&bit[c]) == (undirected_poset[b]&bit[c]))
                which = 2;
            else
                which = 3;
        }

        if (options->userautomproc) 
        {
            if (which == 1)   /* (ab) and (bc) */
            {
                transposition[number_of_generators++] = bit[a]|bit[b];
                transposition[number_of_generators++] = bit[b]|bit[c];
            }
            else if (which == 2)  /* (ab) and (cd) */
            {   
                transposition[number_of_generators++] = bit[a]|bit[b];
                transposition[number_of_generators++] = bit[c]|bit[d];
            }
            else  /* (ab)(cd) */
            {
                for (i = 0; i < n; ++i) tmp_perm[i] = i;
                tmp_perm[a] = b; tmp_perm[b] = a;
                tmp_perm[c] = d; tmp_perm[d] = c;
                (options->userautomproc)(0,tmp_perm,NULL,0,0,n);
            }
        }

        if (orbits)
        {
            for (i = 0; i < n; ++i) orbits[i] = i;
            if (which == 1)  /* {a,b,c} */
            {
                min = (a < b ? a : b); min = (min < c ? min : c);
                orbits[a] = orbits[b] = orbits[c] = min;
            }
            else    /* {a,b} and {c,d} */
            {
                orbits[a] = orbits[b] = (a < b ? a : b);
                orbits[c] = orbits[d] = (c < d ? c : d);
            }
        }

        if (options->getcanon && which == 3)
        {
            /* lab is already in canonical order except possibly in the
               case which == 3.  In that case c and d might need swapping. */

            if (undirected_poset[a] & bit[c])
            { 
                lab[tpos2] = d; lab[tpos2+1] = c;
            }
        }

        return TRUE;
    }
    else if (cells_are_orbits)
    {
        ++nauty_count4;
        stats.numorbits = numcells;

        /* We know the canonicity test (if any) has been passed.  However,
        we might need the group and the orbits.  If only the orbits, we
        will make them from the cells.  Otherwise we need nauty. */

        if (options->userautomproc)
        {
            orbs = (orbits ? orbits : tmp_orbits);
            active = 0;   /* No need to refine further */
            nauty(undirected_poset,lab,ptn,&active,orbs,options,&stats,
                                               workspace,2,1,n,canong);
            return TRUE;
        }

        if (orbits)
        {
            for (i = 0; i < n; ++i)
            {
                min = lab[i];
                for (j = i; ptn[j] > 0; )
                {
                    ++j;
                    if (lab[j] < min) min = lab[j];
                }
                
                /* Now i..j is the cell and min is its least vertex */
                for ( ; i <= j; ++i) orbits[lab[i]] = min;
            }
        }

        return TRUE;
    }
    else
    {
        ++nauty_count5;
        orbs = (orbits ? orbits : tmp_orbits);
        active = 0;   /* No need to refine further */
        nauty(undirected_poset,lab,ptn,&active,orbs,options,&stats,
                                               workspace,2,1,n,canong);

        if (options->getcanon)
        {
            for (i = n-1; (bit[lab[i]]&good_vertices) == 0 ; --i) {}
            return orbs[lab[i]] == orbs[n-1];
        }

        return TRUE;
    }
}



/**********************SCHREIBENBS*************SCHREIB_ALLE_NBS******************************************/

static void schreibenbs(graph nbs[])
{
  int i;

  fprintf(stderr,"nbs: ");
  for (i=-1; (i=nextelement(nbs,1,i))>=0; ) fprintf(stderr,"%d ",i);
  fprintf(stderr,"\n");
}

static void schreib_alle_nbs(graph nbs[])
{

  fprintf(stderr,"Die nbsliste:\n"); 
  for ( ;*nbs!=0; nbs++) schreibenbs(nbs);
  fprintf(stderr,"\n");
}

/*******************************DO_SUBSETS***************************************/

static void do_subsets(int first, int end, graph **list, graph buffer)

/* really fills in the entries (started by all_subsets) */

{
  graph oldbuffer;

  oldbuffer=buffer;

  /* Diese Routine fuellt jeweils das naechste Element ein */

  for ( ; first < end; first++)
    { 
      ADDELEMENT(&buffer,first);
      do_subsets(first+1, end, list, buffer);
      buffer=oldbuffer;
    }

  ADDELEMENT(&buffer,end); 
  **list=buffer; 
  (*list)++; 

  **list=oldbuffer; /* oldbuffer ist hier nur beim ersten Level leer -- da wird
                       die Abschlussnull gesetzt */
  (*list)++;

  return;
}

/*******************************ALL_SUBSETS***************************************/

static void all_subsets(int start, int end, graph list[])

/* Computes all nonempty subsets of the set {start, .... ,end} in
   lexicographically decreasing order and writes them into list[]. In
   fact the lexicographic value is just the opposite of the numerical
   one, that is: {0} is larger than {1}, {1,3} larger than {2,3}
   larger than {2} etc.

 */

{

  graph *runlist;

  runlist=list;
  if (end<start) { fprintf(stderr,"Error in all_subsets: end= %d start=%d\n",end, start);
                   exit(1); }

  do_subsets(start,end,&runlist,(graph)0);
}

/***********************************USAGE****************************************/

static void usage(char name[])

{ fprintf(stderr,"\nUsage: %s n [o|t] [q] [m x y] where n <= 16 is the number of points\n",name );
  fprintf(stderr,"   Generate the Hasse diagrams of the posets with n points\n");
  fprintf(stderr,"   o  causes digraph6 output in arbitrary order to be written to stdout\n");
  fprintf(stderr,"   t  causes digraph6 output in topological order to be written to stdout\n");
  fprintf(stderr,"   q  supresses statistics except for the final count\n");
  fprintf(stderr,"   m x y  with 0 <= x < y divides the generation\n"
              "        into y parts and writes only part x.\n");
  exit(0); } 

/**********************************SCHREIBEPOSET***********************************/

static void schreibeposet()
{
  int i,j;
fprintf(stderr,"\n\nDas Poset:\n");
  for (i=0; i<punktzahl; i++)
    { fprintf(stderr,"%d:",i); 
    for (j=-1; (j=nextelement(poset+i,1,j))>=0 ; ) fprintf(stderr," %d",j);
    fprintf(stderr,"\n");
    }


}

/**********************************SCHREIBEPOSET***********************************/

static void schreibeposet2()
{
  int i,j;
fprintf(stderr,"\n\nDas ungerichtete Poset:\n");
  for (i=0; i<punktzahl; i++)
    { fprintf(stderr,"%d:",i); 
    for (j=-1; (j=nextelement(undirected_poset+i,1,j))>=0 ; ) fprintf(stderr," %d",j);
    fprintf(stderr,"\n");
    }


}

/**********************************AUFSCHREIBEN***********************************/

static void aufschreiben()
{
#if 0
 graph puffer;
      puffer=punktzahl;
      fwrite(&puffer,sizeof(graph),1,stdout);
      fwrite(poset,sizeof(graph),punktzahl,stdout);
#endif
      writeposet(poset,punktzahl);
}

/**********************************AUFSCHREIBEN_2***********************************/

static void aufschreiben_2(int last_level, graph maxpoints,int relationscounter)

     /* Hier werden die Graphen konstruiert, die mit einer Kette von
        jeweils nur einem Knoten pro Level anfangen. Sie entstehen aus
        denen mit weniger Punkten, indem das gesamte erste Level mit
        einem zusaetzlichen Punkt verbunden wird und danach eine Kette
        entsprechender Laenge angehaengt wird. */

{
#if 0
 graph puffer;
#endif
 graph poset2[MAXN];
 char *i2;
 int i;

  if (mod && (punktzahl<=splitpunktzahl)) { splitzaehler++;
                                            if (splitzaehler==rest) splitzaehler-=mod;
                                            else return; }

  minmax_counter[1][POPCOUNT(maxpoints)]++;

  anzahl_posets++;
  i=maxpunktzahl-punktzahl;
  relations[relationscounter+(punktzahl*i)+((i*(i-1))/2)]++;

  punkte_letztes_level[points_on_level[last_level]]++;
  anzahl_level[last_level+(maxpunktzahl-punktzahl)]++;

  if (output)
    {
#if 0
      puffer=maxpunktzahl;
      fwrite(&puffer,sizeof(graph),1,stdout);
#endif
      memcpy(poset2,poset,sizeof(graph)*punktzahl);
      poset2[maxpunktzahl-1]=0;
      for (i=punktzahl; i<(maxpunktzahl-1); i++) poset2[i]=bit[i+1];
      FOREACH(i2,all_points_up_to[0]) ADDELEMENT(poset2+(*i2),punktzahl);
#if 0
      fwrite(poset2,sizeof(graph),maxpunktzahl,stdout);
#endif
      writeposet(poset2,maxpunktzahl);
    }
}

/**********************************AUFSCHREIBEN_3***********************************/

static void aufschreiben_3()

     /* Nur fuer den Sonderfall, dass das gesamte poset eine Kette ist */

{
#if 0
 graph puffer;
#endif
 graph poset2[MAXN];
 int i;

  anzahl_posets++;
  relations[(maxpunktzahl*(maxpunktzahl-1))/2]++;

  punkte_letztes_level[1]++;
  anzahl_level[maxpunktzahl-1]++;
  minmax_counter[1][1]++;

  if (output)
    {
#if 0
      puffer=maxpunktzahl;
      fwrite(&puffer,sizeof(graph),1,stdout);
#endif
      poset2[maxpunktzahl-1]=0;
      for (i=0; i<(maxpunktzahl-1); i++) poset2[i]=bit[i+1];
#if 0
      fwrite(poset2,sizeof(graph),maxpunktzahl,stdout);
#endif
      writeposet(poset2,maxpunktzahl);
    }

}

/*********************************compute_orbits**************************************/

static int compute_orbits(graph *all_nbs,graph *nonequiv_nbs, int bestimme_nummern, int level,
                    int minorbitzahl, int *maxorbsize)

/* Berechnet die Orbits und schreibt ein Element jedes Orbits in nonequiv_nbs -- und zum
   Abschluss eine 0 

   Gibt die Anzahl der Orbits zurueck, wenn neue Nummern berechnet
   werden (d.h. noch keiner auf dem Orbit) und die Anzahl der Orbits mit
   STRIKT groesserer Orbitzahl als minorbitzahl sonst.

*/

{
static graph locallist[MAXGRAPH+1];
static graph *listend, *listrun, image;
static  int i, orbsize;
static graph orbitcounter;
static graph *orbitnummer;


 RESETMARKS;
 orbitnummer=orbitnumber[level-1];

 if (bestimme_nummern) /* orbitzahlen belegen -- dann muss auch maxorbsize belegt werden */
   { orbitcounter=0;
   *maxorbsize=1;
   for ( ; *all_nbs != 0; all_nbs++)
     if (UNMARKED(*all_nbs))
       { *nonequiv_nbs=*all_nbs; nonequiv_nbs++;
       listend=listrun=locallist;
       *listend=*all_nbs; listend++;
       /* Breitensuche: */
       orbitcounter++;
       MARK(*listrun); orbitnummer[*listrun]=orbitcounter;
       while (listrun != listend)
         { 
           for (i=0; i<number_of_generators; i++)
             { PERMSET(*listrun,i,image); /*permset(listrun,&image,1,generators[i]); */
               if (UNMARKED(image))
                 { orbitnummer[image]=orbitcounter; MARK(image); *listend=image; listend++; }
             }
           listrun++;
         }
       orbsize=listend-locallist;
       if (orbsize> *maxorbsize) *maxorbsize=orbsize;
       }
   }/* ende orbitzahlen neu belegen */
else 
  { /* nur nbs mit geeigneter Orbitzahl rausschreiben */
    orbitcounter=0;
 for ( ; *all_nbs != 0; all_nbs++)
   if ((orbitnummer[*all_nbs]>=minorbitzahl) && UNMARKED(*all_nbs))
     { *nonequiv_nbs=*all_nbs; nonequiv_nbs++;
     if (orbitnummer[*all_nbs]>minorbitzahl) orbitcounter++;
       listend=listrun=locallist;
       *listend=*all_nbs; listend++;
       /* Breitensuche: */
       MARK(*listrun);
       while (listrun != listend)
         { 
           for (i=0; i<number_of_generators; i++)
             { PERMSET(*listrun,i,image); /*permset(listrun,&image,1,generators[i]); */
               if (UNMARKED(image))
                 { MARK(image); *listend=image; listend++; }
             }
           listrun++;
         }
     }
  }/* ende points_on_level !=0 */

 *nonequiv_nbs=0;
 return (int)orbitcounter;


}





/**********************************COMPARE_SETS******************************************/

static int compare_sets(const void *a, const void *b)
{ return ((int)(orbits[*((int*)a)]-(int)orbits[*((int*)b)])); }


/**********************************ADD_LAST_POINT******************************************/

/* Fuegt den letzten Punkt hinzu -- wie add_point, nur dass Sachen, die der Weiterfuehrung
   dienen nicht mehr gemacht werden. Darf nur fuer den letzten Punkt aufgerufen werden ! */

static void add_last_point (int level, int do_nauty, graph *last_nbs, graph maxpoints, int relationscounter)

     /* do_nauty ist genau dann 0, wenn nauty gerade ausgefuehrt wurde, um die
        Kanonizitaet des letzten Knotens zu testen, orbits, generators, etc 
        also noch aktuell sind */
{
  graph *local_nonequiv_nbs; /* hier werden die wirklich verschiedenen nbs gespeichert.
                                Bei trivialer Gruppe gleich nbs[level-1] */
  graph *previous_nbs;
  int *local_lab;
  int *local_ptn;
  int nautylab[MAXN], nautyptn[MAXN];
  graph *local_orbitnumber;
  int i, j, l, test, test2, merke;
  int newlevel=0;
  graph dummy, dummy2, *run, *ende, *local_subsets, *local_nbs;
  setword good_vertices;
  boolean accept;
  int op_triv=0; /* Ist die Operation der Gruppe auf den Nachbarschaften leicht als
                    trivial zu erkennen ? */
  unsigned long int faerbung[CHECKSIZE]; /* Fuer eine heuristische Ueberpruefung, ob die Gruppe
                                           trivial operiert */
  unsigned long int merkef;
  graph which_nbs[MAXN];
  int nicht_kanonisch, levelm1, anzahl_save_orbits;
 char *i2;

levelm1=level-1;


if (points_on_level[level]==0)
  {
    if (level>1) all_points_on[levelm1]=(all_points_up_to[levelm1]&(~all_points_up_to[level-2]));

    
    /* Jetzt die direkten oberen Nachbarn: */
    
    anzahl_subsets[levelm1]=(1<<points_on_level[levelm1])-1;
    local_subsets=subsets[levelm1];
    all_subsets(first_point_level[levelm1],punktzahl-1,local_subsets);
    
    
    
    /* Jetzt die Vorgaenger fuer diese neuen Mengen berechnen */
    for (i=anzahl_subsets[levelm1]-1; i>=0; i--)
      { 
        if (POPCOUNT(local_subsets[i])>1) /* die sind schon festgelegt und wuerden sogar zerstoert, wenn
                                             sie behandelt wuerden */
          { dummy=local_subsets[i];
          all_ancestors[dummy]=0;
          FOREACH(i2,dummy) all_ancestors[dummy] |= all_ancestors[bit[*i2]];
          }
      }
    
    
    /* Jetzt mit den vorigen Nachbarschaften kombinieren: */
    
    if (level==1) {local_nbs=local_subsets; 
                   memcpy(nbs[0],local_subsets,sizeof(graph)*anzahl_subsets[0]);
                   anzahl_nbs[levelm1]=anzahl_subsets[levelm1]; }
    else
      { local_nbs=run=nbs[levelm1]; 
      for (i=0; (dummy=local_subsets[i])!=0; i++)
        { *run=dummy;  run++;
        dummy2=all_ancestors[dummy];
        for (l=level-2; (l>=0) && (dummy2!=all_points_up_to[l]); l--)
          if ((~dummy2) & all_points_on[l]) /* since all nbs on this level contain at least
                                             one point from this level, there should be no
                                             in dummy2 */
            { merke=first_point_level[l];
            for (previous_nbs=nbs[l] ; *previous_nbs !=0; merke++ )
              { 
              if (ISELEMENT(&dummy2,merke)) previous_nbs=jump[l][merke];
              else
                { 
                  ende=jump[l][merke];
                  for (;*previous_nbs!=0;previous_nbs++)
                    if ((dummy2 & *previous_nbs)==0) /* Nur wenn der Schnitt leer ist, ist es
                                                        eine moegliche Nachbarschaft */
                      { *run=(dummy | *previous_nbs);  
                      all_ancestors[*run] = dummy2 | all_ancestors[*previous_nbs];
                      run++; }
                }
              } /*ende for ueber alle vorherigen nbs */
            }
        }
      *run=0;
      anzahl_nbs[levelm1]=run-local_nbs;
      /* Jetzt ein paar Heuristiken, bei denen man erkennt, dass die Gruppe trivial operiert */

    if ((anzahl_nbs[levelm1]==1) 
        || ((anzahl_nbs[levelm1]==2) && (POPCOUNT(local_nbs[0])!=POPCOUNT(local_nbs[1]))))
      { maxorbsize_level[levelm1]=1; 
        op_triv=1; }
      else
        if (anzahl_nbs[levelm1] < MAXCHECK)
          { local_orbitnumber=orbitnumber[0];
          for (j=anzahl_nbs[levelm1]-1; j>=0; j--) 
            faerbung[j]=POPCOUNT(local_nbs[j]&all_points_on[levelm1])
              +(local_orbitnumber[local_nbs[j]&all_points_up_to[0]]<<1);
          for (i=level-2; i>0; i--)
            { local_orbitnumber=orbitnumber[i];
            for (j=anzahl_nbs[levelm1]-1; j>=0; j--)
              faerbung[j]+=
                (local_orbitnumber[local_nbs[j]&all_points_on[i]]<<(i+1));
            }
          /* wenn alle Farben verschieden sind operiert die Gruppe mit Sicherheit trivial */
          for (j=anzahl_nbs[levelm1]-1,test=1; j>0; j--)
            for (i=0; i<j; i++) if (faerbung[i]==faerbung[j]) test=0;
          if (test)
            { 
              maxorbsize_level[levelm1]=1;
              op_triv=1; }
          } /* Ende 2. Heuristik */
    }
  local_lab=lab[levelm1]; local_ptn=ptn[levelm1];
  } /* Ende neues level */
 else { local_nbs= nbs[levelm1];  /* existieren schon */ local_lab=lab[level]; local_ptn=ptn[level]; }



 if ((((points_on_level[level]==0) && !op_triv) || (maxorbsize_level[levelm1]!=1)) && do_nauty)
   { 
     number_of_generators=0;
     if ((points_on_level[level]==0) && (maxorbsize_level[level-2]==1)) /* fuer level==1 wird das nicht 
                                                                           aufgerufen, da do_nauty=0 */
     {
       for (i=0; i<=first_point_level[levelm1];i++) orbits[i]=i;
       for ( ; i<=punktzahl-1; i++)
         { if (poset[i]==poset[i-1])
           { 
             transposition[number_of_generators++] = bit[i-1]|bit[i];
             orbits[i]=orbits[i-1];
           }
         else orbits[i]=i;
         }
       if (number_of_generators==0) stats.numorbits=punktzahl; else stats.numorbits=0; /* Hack */
     }
     else {
       memcpy(nautylab,local_lab,sizeof(int)*punktzahl);
       memcpy(nautyptn,local_ptn,sizeof(int)*punktzahl);
       call_nauty(nautylab,nautyptn,punktzahl,orbits,0,&options);
     }
}


if (points_on_level[level]==0) 
  { newlevel=1; 
  first_point_level[level]=punktzahl;
  if (op_triv || (stats.numorbits==punktzahl)) maxorbsize_level[levelm1]=1; 
  else maxorbsize_level[levelm1]=INT_MAX; /* erstmal */
 }

  
/* Jetzt die wirklich verschiedenen Nachbarschaften berechnen */

 maxpoints |= bit[punktzahl];

/* Nun zum Hinzufuegen des Punktes: */

 if (points_on_level[level]==0)
 {

   if (maxorbsize_level[levelm1]==1) { local_nonequiv_nbs=local_nbs; 
                                       local_orbitnumber=orbitnumber[levelm1];
                                       anzahl_save_orbits=anzahl_nbs[levelm1];
                                      }
   else { 
     local_nonequiv_nbs=nonequiv_nbs[punktzahl]; /* punktzahl, bevor noch einer hinzugefuegt wird. 
                                                    Hier wird nur Speicher bereitgestellt. */
     anzahl_save_orbits=
       compute_orbits(local_nbs,local_nonequiv_nbs,1,level,0,maxorbsize_level+levelm1);
   }


   anzahl_posets+=anzahl_save_orbits;
   punkte_letztes_level[points_on_level[level]+1]+=anzahl_save_orbits;
   anzahl_level[level]+=anzahl_save_orbits;
   for (run=local_nonequiv_nbs; *run!=0; run++)
     { dummy=all_ancestors[*run];
     relations[relationscounter+POPCOUNT(dummy)]++;
     minmaxline[POPCOUNT(maxpoints & ~dummy)]++;
#ifndef COUNT
     if (output)
       { poset[punktzahl]=*run; punktzahl++; aufschreiben(); punktzahl--; }
#endif
     }
   } /* ende noch kein Punkt auf dem Level */
   
 else /* d.h. es ist mindestens ein Punkt bereits auf dem Level */
 {
   local_orbitnumber=orbitnumber[levelm1];

   if (maxorbsize_level[levelm1]==1) { local_nonequiv_nbs=last_nbs;
                                       i= last_nbs-nbs[levelm1];
                                       anzahl_save_orbits=anzahl_nbs[levelm1]-i;
 }
   else { 
     local_nonequiv_nbs=nonequiv_nbs[punktzahl]; /* punktzahl, bevor noch einer hinzugefuegt wird 
                                                    hier wird nur Speicher bereitgestellt. */
     anzahl_save_orbits=
       compute_orbits(local_nbs,local_nonequiv_nbs,0,level,local_orbitnumber[poset[punktzahl-1]],0);
   }

 /* Jetzt wirklich den Punkt hinzufuegen */

 local_lab[punktzahl]=punktzahl; local_ptn[punktzahl]=0;
 local_ptn[punktzahl-1]=1; /* nauty wird hier nur gerufen, wenn die orbitzahlen der beiden
                              letzten gleich sind */

 punktzahl++;
 points_on_level[level]++;

 anzahl_posets+=anzahl_save_orbits;
 punkte_letztes_level[points_on_level[level]]+=anzahl_save_orbits;
 anzahl_level[level]+=anzahl_save_orbits;
 /* die mit strikt groesserer Orbitnummer koennen schon mal gesammelt hinzugezaehlt werden */


 for (run=local_nonequiv_nbs; *run!=0; run++)
   { poset[punktzahl-1]=*run;

   if ((maxorbsize_level[levelm1]==1) || 
       (local_orbitnumber[*run]>local_orbitnumber[poset[punktzahl-2]]))
     { dummy=all_ancestors[*run];
       relations[relationscounter+POPCOUNT(dummy)]++;
       minmaxline[POPCOUNT(maxpoints & ~dummy)]++;
#ifndef COUNT
     if (output)
       { aufschreiben(); }
#endif
     }
   else /* zunaechst noch zwei Heuristiken: Wenn alle mit der gleichen orbitnummer identische
           Nachbarschaften haben (d.h. nicht verschiedene, die nur im gleichen Orbit liegen), 
           so ist es auch OK */
     { j=local_orbitnumber[*run];
     for (i=punktzahl-2, test=1; (local_orbitnumber[poset[i]]==j) 
            && (i>=first_point_level[level]);i--) if (poset[i]!=*run) test=0;
     if (test)
       { anzahl_posets++;
       punkte_letztes_level[points_on_level[level]]++;
       anzahl_level[level]++;
       dummy=all_ancestors[*run];
       relations[relationscounter+POPCOUNT(dummy)]++;
       minmaxline[POPCOUNT(maxpoints & ~dummy)]++;
#ifndef COUNT
     if (output)
       { aufschreiben(); }
#endif
       }
     else /* noch ein Sonderfall: Wenn die maximale Orbitgroesse 2 ist und noch kein Knoten
             mit Nachbarn eines anderen Orbits auf dem Level, so ist der Knoten kanonisch,
             der zu der nbs mit mehr Nachbarn auf dem Level adjazent ist. Bei Gleichheit
             sind die Knoten aequivalent */
       if ((maxorbsize_level[levelm1]==2) && 
           (local_orbitnumber[*run]==local_orbitnumber[poset[first_point_level[level]]]))
         { j=local_orbitnumber[*run];
         dummy=*run; test=1; test2=0;
         for (i=punktzahl-2; (local_orbitnumber[poset[i]]==j) && (i>=first_point_level[level]);i--) 
           if (poset[i]==dummy) test++; else test2++;
         if (test>=test2)
           { anzahl_posets++;
           punkte_letztes_level[points_on_level[level]]++;
           anzahl_level[level]++;
           dummy=all_ancestors[*run];
           relations[relationscounter+POPCOUNT(dummy)]++;
           minmaxline[POPCOUNT(maxpoints & ~dummy)]++;
#ifndef COUNT
           if (output)
             { aufschreiben(); }
#endif
           }
         }
   else /* sonst muss vielleicht doch nauty gefragt werden, ob der letzte Punkt kanonisch ist */
     { 
       /* zunaechst aber noch eine Heuristik: Den Knoten werden Farben abhaengig von ihrer Nachbarschaft
          zugeordnet. Sie bestehen zum einen aus der Anzahl der Knoten mit ihrer Nachbarschaft und
          zum anderen aus dem Schnitt mit den anderen nbs der Knoten des gleichen Levels. Der
          kanonische muss auf jeden Fall die gleiche Orbitnummer haben wie der letzte. */

       RESETMARKS;
       for (i=points_on_level[level]; i>0; i--) { faerbung[i]=0; which_nbs[i]=0; }
       faerbung[0]=1; which_nbs[0]=*run;
       MARK(*run);
       l=local_orbitnumber[*run];
       for (i=punktzahl-2; (i>=first_point_level[level]) && (local_orbitnumber[poset[i]]==l); i--)
         { for (j=0; (which_nbs[j]!=0) && (which_nbs[j]!=poset[i]); j++);
         if (which_nbs[j]==0) { which_nbs[j]=poset[i]; faerbung[j]=1; } else faerbung[j]++;
         }
       test=1; merkef=faerbung[0]; nicht_kanonisch=0;
       for (j=1; (which_nbs[j]!=0); j++) { if (faerbung[j]==merkef) { test=0; /* evtl. nicht zu entscheiden */
                                                                     MARK(which_nbs[j]); }
                                           else if (faerbung[j]>merkef) nicht_kanonisch=1; }

       if ((test==0) && (nicht_kanonisch==0)) /* dann probieren wir eine etwas kompliziertere
                                                 Faerbung */
         { 
         test=1;
         for (i=punktzahl-2, dummy=0; (i>=first_point_level[level]); i--)
           if (UNMARKED(poset[i])) dummy|=poset[i];

         RESETMARKS;
         MARK(*run);
         dummy2=(*run & dummy);
         faerbung[0]+=(orbitnumber[0][dummy2 &(all_points_up_to[0])]);
         for (i=levelm1;i>0; i--) 
           faerbung[0]+=(orbitnumber[i][dummy2 &all_points_on[i]])<<i;
         /* Teilmengen der Punkte eines Levels haben immer eine Orbitnummer */



         for (j=1; (which_nbs[j]!=0); j++) 
           if (faerbung[j]==merkef) 
                { dummy2=(which_nbs[j] & dummy);
                  faerbung[j]+=(orbitnumber[0][dummy2 &(all_points_up_to[0])]);
                  for (i=levelm1;i>0; i--) 
                    faerbung[j]+=
                      (orbitnumber[i][dummy2 &all_points_on[i]])<<i;

                  if (faerbung[j]==faerbung[0]) { test=0; MARK(which_nbs[j]); }
                  else if (faerbung[j]>faerbung[0]) nicht_kanonisch=1; }
         }


       if ((test==0) && (nicht_kanonisch==0)) /* OK -- dann halt doch nauty... */
         {

           undirected_poset[punktzahl-1]=*run;
           FOREACH(i2,*run) ADDELEMENT(undirected_poset+(*i2),punktzahl-1);

           number_of_generators=0;
           memcpy(nautylab,local_lab,sizeof(int)*punktzahl);
           memcpy(nautyptn,local_ptn,sizeof(int)*punktzahl);
           good_vertices = 0;
           for (i = first_point_level[level]; i < punktzahl; ++i)
               if (ISMARKED(poset[i])) ADDELEMENT(&good_vertices,i);
           accept = call_nauty(nautylab,nautyptn,punktzahl,/*TMPorbits*/NULL,good_vertices,&options_final);
           FOREACH(i2,*run) DELELEMENT(undirected_poset+(*i2),punktzahl-1);

           undirected_poset[punktzahl-1]=0;
           j=local_orbitnumber[poset[punktzahl-1]]; /* nur zum Puffern */
           if (accept)
             { anzahl_posets++;
             punkte_letztes_level[points_on_level[level]]++;
             anzahl_level[level]++;
             dummy=all_ancestors[*run];
             relations[relationscounter+POPCOUNT(dummy)]++;
             minmaxline[POPCOUNT(maxpoints & ~dummy)]++;
#ifndef COUNT
             if (output)
               { aufschreiben(); }
#endif
             } /* ende akzeptiert mit nauty */

         } /* ende halt doch nauty */
       else /* konnte ohne nauty entschieden werden */
         if ((test==1) && (nicht_kanonisch==0))
             { anzahl_posets++;
             punkte_letztes_level[points_on_level[level]]++;
             anzahl_level[level]++;
             dummy=all_ancestors[*run];
             relations[relationscounter+POPCOUNT(dummy)]++;
             minmaxline[POPCOUNT(maxpoints & ~dummy)]++;
#ifndef COUNT
             if (output)
               { aufschreiben(); }
#endif
             }
     } /* else nauty... */
     } /* else mit Heuristik... */
   }
 DELELEMENT(all_points_up_to+level,punktzahl-1);
 points_on_level[level]--;
 punktzahl--;
 local_ptn[punktzahl-1]=0; 
 } /* ende mindestens ein Punkt bereits auf dem Level */

}







/**********************************ADD_POINT******************************************/

/* Fuegt einen weiteren Punkt hinzu */

static void add_point (int level, int do_nauty, graph *last_nbs, graph maxpoints, int relationscounter)

     /* do_nauty ist genau dann 0, wenn nauty gerade ausgefuehrt wurde, um die
        Kanonizitaet des letzten Knotens zu testen, orbits, generators, etc 
        also noch aktuell sind */
{
  graph *local_nonequiv_nbs; /* hier werden die wirklich verschiedenen nbs gespeichert.
                                Bei trivialer Gruppe gleich nbs[level-1] */
  graph *previous_nbs;
  int *local_lab;
  int *local_ptn;
  int nautylab[MAXN], nautyptn[MAXN];
  graph *local_orbitnumber;
  int i, j, l, test, test2, merke, nr;
  int newlevel=0;
  graph dummy, dummy2, dummymax, *run, *ende, *local_subsets, *local_nbs;
  setword good_vertices;
  boolean accept;

  int merkegenerators[MAXN][MAXN];
  int merke_numgen;
  int merkeorbits[MAXN];
  setword merketrans[MAXN];
  int merkestatsnumorbits;
  int op_triv=0; /* Ist die Operation der Gruppe auf den Nachbarschaften leicht als
                    trivial zu erkennen ? */
  unsigned long int faerbung[CHECKSIZE]; /* Fuer eine heuristische Ueberpruefung, ob die Gruppe
                                           trivial operiert */
  unsigned long int merkef;
  graph which_nbs[MAXN];
  int nicht_kanonisch, levelm1;
  char *i2;

levelm1=level-1;


  if (mod && (punktzahl==splitpunktzahl)) { splitzaehler++;
                                            if (splitzaehler==rest) splitzaehler-=mod;
                                            else return; }

if (points_on_level[level]==0)
  {
    if (level>1) all_points_on[levelm1]=(all_points_up_to[levelm1]&(~all_points_up_to[level-2]));
    
    /* Jetzt die direkten oberen Nachbarn: */
    
    anzahl_subsets[levelm1]=(1<<points_on_level[levelm1])-1;
    local_subsets=subsets[levelm1];
    all_subsets(first_point_level[levelm1],punktzahl-1,local_subsets);
    
    
    
    /* Jetzt die Vorgaenger fuer diese neuen Mengen berechnen */
    for (i=anzahl_subsets[levelm1]-1; i>=0; i--)
      { 
        if (POPCOUNT(local_subsets[i])>1) /* die sind schon festgelegt und wuerden sogar zerstoert, wenn
                                             sie behandelt wuerden */
          { dummy=local_subsets[i];
          all_ancestors[dummy]=0;
          FOREACH(i2,dummy) all_ancestors[dummy] |= all_ancestors[bit[*i2]];
          }
      }
    
    
    /* Jetzt mit den vorigen Nachbarschaften kombinieren: */
    
    if (level==1) {local_nbs=local_subsets; 
                   memcpy(nbs[0],local_subsets,sizeof(graph)*anzahl_subsets[0]);
                   anzahl_nbs[levelm1]=anzahl_subsets[levelm1]; }
    else
      { local_nbs=run=nbs[levelm1]; 
      for (i=0; (dummy=local_subsets[i])!=0; i++)
        { *run=dummy;  run++;
        dummy2=all_ancestors[dummy];
        for (l=level-2; (l>=0) && (dummy2!=all_points_up_to[l]); l--)
          if ((~dummy2) & all_points_on[l])
            { merke=first_point_level[l];
            for (previous_nbs=nbs[l] ; *previous_nbs !=0; merke++ )
              { 
              if (ISELEMENT(&dummy2,merke)) previous_nbs=jump[l][merke];
              else
                { 
                  ende=jump[l][merke];
                  for (;*previous_nbs!=0;previous_nbs++)
                    if ((dummy2 & *previous_nbs)==0) /* Nur wenn der Schnitt leer ist, ist es
                                                        eine moegliche Nachbarschaft */
                      { *run=(dummy | *previous_nbs);  
                      all_ancestors[*run] = dummy2 | all_ancestors[*previous_nbs];
                      run++; }
                }
              } /*ende for ueber alle vorherigen nbs */
            }
        }
      *run=0;
      anzahl_nbs[levelm1]=run-local_nbs;

      /* Jetzt ein paar Heuristiken, bei denen man sofort erkennt, dass die Gruppe trivial operiert */

    if ((anzahl_nbs[levelm1]==1) 
        || ((anzahl_nbs[levelm1]==2) && (POPCOUNT(local_nbs[0])!=POPCOUNT(local_nbs[1]))))
      { maxorbsize_level[levelm1]=1; 
        op_triv=1; }
      else
        if (anzahl_nbs[levelm1] < MAXCHECK)
          { local_orbitnumber=orbitnumber[0];
          for (j=anzahl_nbs[levelm1]-1; j>=0; j--) 
            faerbung[j]=POPCOUNT(local_nbs[j]&all_points_on[levelm1])
              +(local_orbitnumber[local_nbs[j]&all_points_up_to[0]]<<1);
          for (i=level-2; i>0; i--)
            { local_orbitnumber=orbitnumber[i];
            for (j=anzahl_nbs[levelm1]-1; j>=0; j--)
              faerbung[j]+=
                (local_orbitnumber[local_nbs[j]&all_points_on[i]]<<(i+1));
            }
          /* wenn alle Farben verschieden sind operiert die Gruppe mit Sicherheit trivial */
          for (j=anzahl_nbs[levelm1]-1,test=1; j>0; j--)
            for (i=0; i<j; i++) if (faerbung[i]==faerbung[j]) test=0;
          if (test)
            { 
              maxorbsize_level[levelm1]=1;
              op_triv=1; }
          } /* Ende 2. Heuristik */
      } /* ende level>1 */

    if (punktzahl<=maxpunktzahl-2)
    for (i=first_point_level[levelm1], run=local_nbs; i<punktzahl; i++)
      { while (ISELEMENT(run,i)) run++; 
        jump[levelm1][i]=run; }
    local_lab=lab[levelm1]; local_ptn=ptn[levelm1];
  } /* Ende neues level */
 else { local_nbs= nbs[levelm1];  /* existieren schon */ local_lab=lab[level]; local_ptn=ptn[level]; }



 if ((((points_on_level[level]==0) && !op_triv) || (maxorbsize_level[levelm1]!=1)) && do_nauty)
   { 
     number_of_generators=0;
     if ((points_on_level[level]==0) && (maxorbsize_level[level-2]==1)) /* fuer level==1 wird das nicht 
                                                                           aufgerufen, da do_nauty=0 */
     {
       for (i=0; i<=first_point_level[levelm1];i++) orbits[i]=i;
       for ( ; i<=punktzahl-1; i++)
         { if (poset[i]==poset[i-1])
           { 
             transposition[number_of_generators++] = bit[i-1]|bit[i];
             orbits[i]=orbits[i-1];
           }
         else orbits[i]=i;
         }
       if (number_of_generators==0) stats.numorbits=punktzahl; else stats.numorbits=0; /* Hack */
     }
     else {
       memcpy(nautylab,local_lab,sizeof(int)*punktzahl);
       memcpy(nautyptn,local_ptn,sizeof(int)*punktzahl);
       call_nauty(nautylab,nautyptn,punktzahl,orbits,0,&options);
     }
}


if (points_on_level[level]==0) 
  /* ein neues level -- lab kann neu festgelegt werden */
  { newlevel=1; 
    first_point_level[level]=punktzahl;
    local_lab=lab[level]; local_ptn=ptn[level];
    /* zunaechst lab: */
    if (op_triv || (stats.numorbits==punktzahl)) 
      { maxorbsize_level[levelm1]=1; 
      for (i=0; i<punktzahl; i++)
        { local_lab[i]=i; local_ptn[i]=0;}
      }
    else 
      { maxorbsize_level[levelm1]=INT_MAX; /* erstmal */
      for (i=0; i<punktzahl; i++) local_lab[i]=i;
      qsort((void *)local_lab, punktzahl, sizeof(int), compare_sets);
      local_ptn[punktzahl-1]=0;
      for (i=punktzahl-2; i>=0; i--)
        { if (orbits[local_lab[i+1]] == orbits[local_lab[i]]) local_ptn[i]=1; 
        else local_ptn[i]=0;
        }
      }
  }

  
/* Jetzt die wirklich verschiedenen Nachbarschaften berechnen */


/* Nun zum Hinzufuegen des Punktes: */

 maxpoints |= bit[punktzahl];


 if (points_on_level[level]==0)
 {

   if (maxorbsize_level[levelm1]==1) { local_nonequiv_nbs=local_nbs; 
                                       local_orbitnumber=orbitnumber[levelm1];
                                       for (i=1, run=local_nonequiv_nbs; *run != 0; run++, i++)
                                         local_orbitnumber[*run]=i;
                                      }
   else { 
     local_nonequiv_nbs=nonequiv_nbs[punktzahl]; /* punktzahl, bevor noch einer hinzugefuegt wird. 
                                                    Hier wird nur Speicher bereitgestellt. */
     compute_orbits(local_nbs,local_nonequiv_nbs,1,level,0,maxorbsize_level+levelm1);
     if (maxorbsize_level[levelm1]==1)
       { for (i=0; i<punktzahl; i++) local_ptn[i]=0; local_nonequiv_nbs=local_nbs;}
     /* SATZ: Der so gefaerbte Graph hat die gleichen Orbits auf den nbs auch fuer kommende
        Level und liefert ein korrektes Kanonizitaetskriterium */

   }




 /* Jetzt wirklich den Punkt hinzufuegen */

 local_lab[punktzahl]=punktzahl; local_ptn[punktzahl]=0;
 punktzahl++;
 points_on_level[level]++;
 ADDELEMENT(all_points_up_to+level,punktzahl-1);
 for (run=local_nonequiv_nbs; *run!=0; run++)
   { poset[punktzahl-1]=undirected_poset[punktzahl-1]=*run; 
   FOREACH(i2,*run) ADDELEMENT(undirected_poset+(*i2),punktzahl-1);
   dummy=all_ancestors[*run];
   all_ancestors[bit[punktzahl-1]]= (bit[punktzahl-1] | dummy);
   dummymax=maxpoints & (~dummy);
   nr=relationscounter+POPCOUNT(dummy);
   aufschreiben_2(level,dummymax,nr);
   if (punktzahl==(maxpunktzahl-1)) { add_last_point(level,1,run,dummymax,nr); 
                                      add_last_point(level+1,1,run,dummymax,nr); }
   else { add_point(level,1,run,dummymax,nr); 
          add_point(level+1,1,run,dummymax,nr); }

   FOREACH(i2,*run) DELELEMENT(undirected_poset+(*i2),punktzahl-1);
   }
 poset[punktzahl-1]=undirected_poset[punktzahl-1]=0;
 DELELEMENT(all_points_up_to+level,punktzahl-1);
 points_on_level[level]--;
 punktzahl--;
 } /* ende noch kein Punkt auf dem Level */

 else /* d.h. es ist mindestens ein Punkt bereits auf dem Level */
 {
   local_orbitnumber=orbitnumber[levelm1];

   if (maxorbsize_level[levelm1]==1) local_nonequiv_nbs=last_nbs;
   else { 
     local_nonequiv_nbs=nonequiv_nbs[punktzahl]; /* punktzahl, bevor noch einer hinzugefuegt wird 
                                                    hier wird nur Speicher bereitgestellt. */
     compute_orbits(local_nbs,local_nonequiv_nbs,0,level,local_orbitnumber[poset[punktzahl-1]],0);
   }


 /* Jetzt wirklich den Punkt hinzufuegen */

 local_lab[punktzahl]=punktzahl; local_ptn[punktzahl]=0;
 punktzahl++;
 points_on_level[level]++;
 ADDELEMENT(all_points_up_to+level,punktzahl-1);
 for (run=local_nonequiv_nbs; *run!=0; run++)
   { poset[punktzahl-1]=undirected_poset[punktzahl-1]=*run;
   FOREACH(i2,*run) ADDELEMENT(undirected_poset+(*i2),punktzahl-1);
   if (local_orbitnumber[poset[punktzahl-2]]==local_orbitnumber[*run]) 
                             local_ptn[punktzahl-2]=1; else local_ptn[punktzahl-2]=0;

   dummy=all_ancestors[*run];
   all_ancestors[bit[punktzahl-1]]= (bit[punktzahl-1] | dummy);
   dummymax=maxpoints & (~dummy);
   nr=relationscounter+POPCOUNT(dummy);

   if ((maxorbsize_level[levelm1]==1) || 
       (local_orbitnumber[*run]>local_orbitnumber[poset[punktzahl-2]]))
     {
     aufschreiben_2(level,dummymax,nr);
     if (punktzahl==(maxpunktzahl-1)) { add_last_point(level,1,run,dummymax,nr); 
                                        add_last_point(level+1,1,run,dummymax,nr); }
     else { add_point(level,1,run,dummymax,nr); 
            add_point(level+1,1,run,dummymax,nr); }
     }
   else /* zunaechst noch zwei Heuristiken: Wenn alle mit der gleichen orbitnummer identische
           Nachbarschaften haben (d.h. nicht verschiedene, die nur im gleichen Orbit liegen), 
           so ist es auch OK */
     { j=local_orbitnumber[*run];
     for (i=punktzahl-2, test=1; (local_orbitnumber[poset[i]]==j) 
            && (i>=first_point_level[level]);i--) if (poset[i]!=*run) test=0;
     if (test)
       {
       aufschreiben_2(level,dummymax,nr);
       if (punktzahl==(maxpunktzahl-1)) { add_last_point(level,1,run,dummymax,nr); 
                                          add_last_point(level+1,1,run,dummymax,nr); }
       else { add_point(level,1,run,dummymax,nr); 
              add_point(level+1,1,run,dummymax,nr); }
       }
     else /* noch ein Sonderfall: Wenn die maximale Orbitgroesse 2 ist und noch kein Knoten
             mit Nachbarn eines anderen Orbits auf dem Level, so ist der Knoten kanonisch,
             der zu der nbs mit mehr Nachbarn auf dem Level adjazent ist. Bei Gleichheit,
             sind die Knoten aequivalent */
       if ((maxorbsize_level[levelm1]==2) && 
           (local_orbitnumber[*run]==local_orbitnumber[poset[first_point_level[level]]]))
         { j=local_orbitnumber[*run];
         dummy=*run; test=1; test2=0;
         for (i=punktzahl-2; (local_orbitnumber[poset[i]]==j) && (i>=first_point_level[level]);i--) 
           if (poset[i]==dummy) test++; else test2++;
         if (test>=test2)
           {
             aufschreiben_2(level,dummymax,nr);
             if (punktzahl==(maxpunktzahl-1)) { add_last_point(level,1,run,dummymax,nr); 
                                                add_last_point(level+1,1,run,dummymax,nr); }
             else { add_point(level,1,run,dummymax,nr); 
                    add_point(level+1,1,run,dummymax,nr); }
           }
         }
   else /* sonst muss vielleicht doch nauty gefragt werden, ob der letzte Punkt kanonisch ist */
     { 
       /* zunaechst aber noch eine Heuristik: Den Knoten werden Farben abhaengig von ihrer Nachbarschaft
          zugeordnet. Sie bestehen zum einen aus der Anzahl der Knoten mit ihrer Nachbarschaft und
          zum anderen aus dem Schnitt mit den anderen nbs der Knoten des gleichen Levels. Der
          kanonische muss auf jeden Fall die gleiche Orbitnummer haben wie der letzte. */

       RESETMARKS;
       for (i=points_on_level[level]; i>0; i--) { faerbung[i]=0; which_nbs[i]=0; }
       faerbung[0]=1; which_nbs[0]=*run;
       MARK(*run);
       l=local_orbitnumber[*run];
       for (i=punktzahl-2; (i>=first_point_level[level]) && (local_orbitnumber[poset[i]]==l); i--)
         { for (j=0; (which_nbs[j]!=0) && (which_nbs[j]!=poset[i]); j++);
         if (which_nbs[j]==0) { which_nbs[j]=poset[i]; faerbung[j]=1; } else faerbung[j]++;
         }
       test=1; merkef=faerbung[0]; nicht_kanonisch=0;
       for (j=1; (which_nbs[j]!=0); j++) { if (faerbung[j]==merkef) { test=0; /* evtl. nicht zu entscheiden */
                                                                     MARK(which_nbs[j]); }
                                           else if (faerbung[j]>merkef) nicht_kanonisch=1; }


       if ((test==0) && (nicht_kanonisch==0)) /* dann probieren wir eine etwas kompliziertere
                                                 Faerbung */
         { 
         test=1;
         for (i=punktzahl-2, dummy=0; (i>=first_point_level[level]); i--)
           if (UNMARKED(poset[i])) dummy|=poset[i];

         RESETMARKS;
         MARK(*run);
         dummy2=(*run & dummy);
         faerbung[0]+=(orbitnumber[0][dummy2 &(all_points_up_to[0])]);
         for (i=levelm1;i>0; i--) 
           faerbung[0]+=(orbitnumber[i][dummy2 &all_points_on[i]])<<i;
         /* Teilmengen der Punkte eines Levels haben immer eine Orbitnummer */



         for (j=1; (which_nbs[j]!=0); j++) 
           if (faerbung[j]==merkef) 
                { dummy2=(which_nbs[j] & dummy);
                  faerbung[j]+=(orbitnumber[0][dummy2 &(all_points_up_to[0])]);
                  for (i=levelm1;i>0; i--) 
                    faerbung[j]+=
                      (orbitnumber[i][dummy2 &all_points_on[i]])<<i;

                  if (faerbung[j]==faerbung[0]) { test=0; MARK(which_nbs[j]); }
                  else if (faerbung[j]>faerbung[0]) nicht_kanonisch=1; }
         }


       if ((test==0) && (nicht_kanonisch==0)) /* OK -- dann halt doch nauty... */
         {
           number_of_generators=0;
           memcpy(nautylab,local_lab,sizeof(int)*punktzahl);
           memcpy(nautyptn,local_ptn,sizeof(int)*punktzahl);
           good_vertices = 0;
           for (i = first_point_level[level]; i < punktzahl; ++i)
               if (ISMARKED(poset[i])) ADDELEMENT(&good_vertices,i);
           accept = call_nauty(nautylab,nautyptn,punktzahl,orbits,good_vertices,&options_canon);
           j=local_orbitnumber[poset[punktzahl-1]]; /* nur zum Puffern */
           if (accept)
             {
               aufschreiben_2(level,dummymax,nr);
               if (punktzahl==(maxpunktzahl-1))
                 { add_last_point(level+1,0,run,dummymax,nr); 
                 /* hier wird nauty garantiert nicht noch einmal
                    aufgerufen, deshalb muss man sich die Sachen 
                    nicht merken */
                 add_last_point(level,0,run,dummymax,nr); }
               else
                 {
                   merkestatsnumorbits=stats.numorbits;
                   merke_numgen=number_of_generators;
                   memcpy(merkeorbits,orbits,sizeof(int)*punktzahl);
                   memcpy(merkegenerators,generators,sizeof(int)*number_of_generators*MAXN);
                   memcpy(merketrans,transposition,sizeof(setword)*number_of_generators);

                   add_point(level,0,run,dummymax,nr);
                   
                   stats.numorbits=merkestatsnumorbits;
                   number_of_generators=merke_numgen;
                   memcpy(orbits,merkeorbits,sizeof(int)*punktzahl);
                   memcpy(generators,merkegenerators,sizeof(int)*number_of_generators*MAXN);
                   memcpy(transposition,merketrans,sizeof(setword)*number_of_generators);

                   add_point(level+1,0,run,dummymax,nr);
                 }
             } /* ende akzeptiert mit nauty */

         } /* ende halt doch nauty */
       else /* konnte ohne nauty entschieden werden */
         if ((test==1) && (nicht_kanonisch==0))
             { 
               aufschreiben_2(level,dummymax,nr);
               if (punktzahl==(maxpunktzahl-1)) { add_last_point(level,1,run,dummymax,nr); 
                                                  add_last_point(level+1,1,run,dummymax,nr); }
               else { add_point(level,1,run,dummymax,nr); add_point(level+1,1,run,dummymax,nr); }
             }
     } /* else nauty... */
     } /* else mit Heuristik... */
   FOREACH(i2,*run) DELELEMENT(undirected_poset+(*i2),punktzahl-1);
   }
 poset[punktzahl-1]=undirected_poset[punktzahl-1]=0;
 DELELEMENT(all_points_up_to+level,punktzahl-1);
 points_on_level[level]--;
 punktzahl--;
 local_ptn[punktzahl-1]=0; 
 } /* ende mindestens ein Punkt bereits auf dem Level */

}

/***************************************MAIN*******************************************/

int main(int argc, char *argv[])

  { int i,j,k,quiet,pruned;
  graph dummy;
  double time0,time1;

  nauty_check(WORDSIZE,1,1,NAUTYVERSIONID);

  if (argc==1) usage(argv[0]); 
  if (argc == 2 && (strcmp(argv[1],"-help") == 0 || strcmp(argv[1],"--help") == 0))
      usage(argv[0]);

  if (!isdigit(argv[1][0])) usage(argv[0]);

  quiet = topsort = 0;
  for (i=2;i<argc;i++)
    {
       if (argv[i][0]=='o') { output=1;
#ifdef COUNT
    fprintf(stderr,"The program was compiled with -DCOUNT, so it cannot output posets !\n");
    exit(1);
#endif
    }
    else
    if (argv[i][0]=='t') { output=1; topsort=1;
#ifdef COUNT
    fprintf(stderr,"The program was compiled with -DCOUNT, so it cannot output posets !\n");
    exit(1);
#endif
    }
    else
    if (argv[i][0]=='q') quiet = 1;
    else 
    if ((argv[i][0]=='m') && isdigit(argv[i+1][0]) &&  isdigit(argv[i+2][0]))
      { i++; rest=atoi(argv[i]); i++; mod=atoi(argv[i]); }
    else usage(argv[0]); }

  if (mod && ((rest<0) || (rest>=mod))) usage(argv[0]);

  if (sizeof(unsigned long int)<4) { fprintf(stderr,">E sizeof(unsigned long int) should be at least 4\n");
                                     exit(1); }

  if (!quiet) { for (i=0;i<argc;i++) fprintf(stderr,"%s ",argv[i]); fprintf(stderr,"\n"); }

  time0 = CPUTIME;

for (i=0;i<MAXN;i++) orbitnumber[i][0]=0;

#if 0   // BDM
for (i=0;i<=MAXGRAPH;i++) 
  { dummy=i;
    first_bit[i]=FIRSTBIT(dummy);
  for (j=0;j<MAXN;j++) 
    { if ((k=nextelement(&dummy,1,j))>=0) next_element[i][j]=k;
    else next_element[i][j]=MAXN;
    if ((k=nextelement(&dummy,1,-1))>=0) next_element[i][MAXN]=k;
    else next_element[i][MAXN]=MAXN;
    }
    for (j=MAXN, k=0; (j=next_element[i][j])!=MAXN; ) { elementlist[i][k]=j; k++; }
    elementlist[i][k]=MAXN;
  }
#else
 for (i=0;i<=MAXGRAPH;i++)
 {
    dummy = i;
    k = 0;
    while (dummy)
    { TAKEBIT(elementlist[i][k],dummy); 
      k++;
    }
    elementlist[i][k]=MAXN;
 }
#endif

  maxpunktzahl=atoi(argv[1]);

  if (maxpunktzahl>MAXN) { fprintf(stderr,">E No more than %d points are supported.\n",MAXN);
                           fprintf(stderr,"Note that changing MAXN is a nontrivial task.\n");
                           exit(1); }

  if (mod && (maxpunktzahl<6)) { fprintf(stderr,">E Need at least 6 vertices for splitting\n");
                                 exit(1); }
  if (mod) { splitpunktzahl=maxpunktzahl-4; 
             if (splitpunktzahl>11) splitpunktzahl=11; }

  init_nauty();

  first_point_level[0]=0;
  for (points_on_level[0]=2; points_on_level[0]<maxpunktzahl; points_on_level[0]++)
    if ((mod==0) || (rest==0) || (points_on_level[0]<=splitpunktzahl))
    { punktzahl=points_on_level[0];
    maxorbsize_level[0]=punktzahl;
    all_points_up_to[0]=0;
    for (i=0; i<punktzahl; i++) { lab[0][i]=i; ptn[0][i]=1; all_ancestors[bit[i]]=bit[i]; 
                                  ADDELEMENT(all_points_up_to,i); orbits[i]=1;}
    all_points_on[0]=all_points_up_to[0];
      ptn[0][punktzahl-1]=0;
      number_of_generators=0;
      for (i=0; i<punktzahl-1; i++)
        { 

          transposition[number_of_generators++] = bit[i]|bit[i+1];
        }
      stats.numorbits=1;

      minmaxline=minmax_counter[punktzahl];
      aufschreiben_2(0,all_points_on[0],0);
      if (punktzahl==(maxpunktzahl-1)) add_last_point(1,0,0,all_points_on[0],0); 
      else add_point(1,0,0,all_points_on[0],0);
    }
  /* Die Faelle nur ein Level bzw ueberall nur ein Punkt extra: */
  if (rest==0) /* d.h. auch wenn garnicht gesplittet wird */
    { 
      if (maxpunktzahl>1) aufschreiben_3();

      for (i=0, all_points_up_to[0]=0; i<maxpunktzahl; i++) 
        { poset[i]=0; ADDELEMENT(all_points_up_to,i); }
      points_on_level[0]=punktzahl=maxpunktzahl; 
      minmaxline=minmax_counter[punktzahl];
      anzahl_posets++;
      punkte_letztes_level[maxpunktzahl]++;
      anzahl_level[0]++;
      relations[0]++;
      minmaxline[punktzahl]++;
#ifndef COUNT
             if (output)
               { aufschreiben(); }
#endif
    }

    time1 = CPUTIME;

if (!quiet)
{
fprintf(stderr,"\nPosets with  1 point on the last level:  %lld\n",punkte_letztes_level[1]);
  for (i=2; i<=maxpunktzahl; i++) fprintf(stderr,"Posets with %2d points on the last level: %lld\n",i,punkte_letztes_level[i]);
  fprintf(stderr,"\n");
  fprintf(stderr,"Posets with  1 level:  %lld\n",anzahl_level[0]);
  for (i=2; i<=maxpunktzahl; i++) fprintf(stderr,"Posets with %2d levels: %lld\n",i,anzahl_level[i-1]);
  fprintf(stderr,"nauty counts: %lld %lld %lld %lld %lld\n",nauty_count1,nauty_count2,nauty_count3,nauty_count4,nauty_count5);

  fprintf(stderr,"\n");
  fprintf(stderr,"Nontrivial relations | posets\n");
  for (i=0; i<=(maxpunktzahl*(maxpunktzahl-1))/2; i++)
    fprintf(stderr,"r(%d)=  %lld\n",i,relations[i]);

  fprintf(stderr,"\n");
  fprintf(stderr,"minimal points | maximal points | posets\n");
  for (i=1; i<=maxpunktzahl;i++)
    for (j=1; j<=maxpunktzahl;j++)
      fprintf(stderr,"p(%d,%d)=  %lld\n",i,j,minmax_counter[i][j]);
}

#ifdef POSET_SUMMARY
   POSET_SUMMARY;
#endif

#if defined(POSET_PRUNE0) || defined(POSET_PRUNE1)
  pruned=output;
#else
  pruned=0;
#endif
  if (pruned)
    fprintf(stderr,">Z %lld posets generated, %lld written; %.2f sec\n",
           anzahl_posets,written_posets,time1-time0);
  else
    fprintf(stderr,">Z %lld posets %s; %.2f sec\n",
           anzahl_posets,(output?"written":"generated"),time1-time0);

  exit(0);
}
