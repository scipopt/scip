/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_tworowcomb.c
 * @brief  derive variable bounds from convex combination of two rows
 * @author Dieter Weninger
 * @author Patrick Gemander
 */

//TODO remove these
//#define SCIP_DEBUG
//#define SCIP_DEBUG_BOUNDS
//#define SCIP_DEBUG_BREAKPOINTS
//#define SCIP_DEBUG_CLIQUE
//#define SCIP_TEST_CLIQUE
//#define SCIP_DEBUG_HASHING
//#define SCIP_MORE_DEBUG
//#define SCIP_DEBUG_SIGNS
//#define SCIP_DEBUG_SUBSCIP
//#define SCIP_HASHBLOCK_INFO

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_matrix.h"
#include "presol_tworowcomb.h"

// TODO Fix aligning
#define PRESOL_NAME            "tworowcomb"
#define PRESOL_DESC            "derive variable bounds from convex combination of two rows"
#define PRESOL_PRIORITY         -4000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS        0 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_MAXCONSIDEREDNONZEROS  100
#define DEFAULT_MAXRETRIEVEFAILS       1000
#define DEFAULT_MAXCOMBINEFAILS        1000
#define DEFAULT_MAXHASHFAC             10
#define DEFAULT_MAXPAIRFAC             1

#define MIN_DENOM 1.e-10

/*
 * Data structures
 */

enum signum {UP, DN, POS, NEG, CLQ};

/** presolver data */
struct SCIP_PresolData
{
   SCIP_Real maxpairfac;
   SCIP_Real maxhashfac;
   int maxretrievefails;
   int maxcombinefails;
   int maxconsiderednonzeros;
   int nchgbnds;
   int nuselessruns;
};


/** structure representing a pair of row indices; used for lookup in a hashtable */
struct RowPair
{
   int row1idx;
   int row2idx;
};

typedef struct RowPair ROWPAIR;


/*
 * Local methods
 */

/** encode contents of a rowpair as void* pointer */
static
void*
encodeRowPair(
   ROWPAIR*              rowpair             /**< pointer to rowpair */
   )
{
   return (void*)SCIPcombineTwoInt(rowpair->row1idx, rowpair->row2idx);
}

/** compute single int hashvalue for two ints */
static
int
hashIndexPair(
   int                   idx1,               /**< first integer index */
   int                   idx2                /**< second integer index */
   )
{
   uint32_t hash = SCIPhashOne(SCIPcombineTwoInt(idx1, idx2));
   return *((int*) &hash);
}

static SCIP_RETCODE addEntry
(
 SCIP* scip,                   /**< SCIP datastructure */
 int* pos,
 int* listsize,
 int** hashlist,
 int** rowidxlist,
 int hash,
 int rowidx
)
{
   if( (*pos) >= (*listsize) )
   {
      int newsize  = SCIPcalcMemGrowSize(scip, (*pos) + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, hashlist, (*listsize), newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, rowidxlist, (*listsize), newsize) );
      (*listsize) = newsize;
   }

   (*hashlist)[(*pos)] = hash;
   (*rowidxlist)[(*pos)] = rowidx;
   (*pos)++;

   return SCIP_OKAY;
}

static SCIP_RETCODE findNextBlock
(
   int*                 list,
   int                  len,
   int*                 start,
   int*                 end
)
{
   int i;
   (*start) = (*end);
   i = (*end) + 1;
   while( i < len && list[i] == list[i-1] )
      i++;

   (*end) = i;

   return SCIP_OKAY;
}

static int calcCliqueMaximums(
   SCIP*                scip,
   int*                 varinds,             /**< variable index array */
   int*                 cliquevarpos,        /**< positions of clique variables in index array */
   int                  cliquesize,          /**< size of current clique */
   SCIP_Real*           row1coefs,           /**< coefficients of first row */
   SCIP_Real*           row2coefs,           /**< coefficients of second row */
   enum signum*         signs,               /**< sign behaviour of variables */
   int*                 nbreakpoints,        /**< number of breakpoints between 0 and 1 */
   SCIP_Real*           breakpoints,         /**< variable breakpoints */
   int*                 cliquemaxinds        /**< array containing in which clique this variable has a maximum */
   )
{
   int i;
   int idx;
   int maxidx;
   int maxpos;
   int newmaxpos;
   int firstmaxpos;
   SCIP_Real lambda;
   SCIP_Real minlambda;
   SCIP_Real* gradients;
   SCIP_Real breakpointval;

#ifdef SCIP_DEBUG_CLIQUE
   SCIPdebugMsg(scip, "calculating maximums for clique %d\n", -cliquemaxinds[cliquevarpos[0]]);
#endif

   //TODO Should I add safeguards for out of bounds stuff?

   SCIP_CALL( SCIPallocBufferArray(scip, &gradients, cliquesize) );

   // calculate gradients
   for( i = 0; i < cliquesize; i++ )
   {
      gradients[i] = row1coefs[varinds[cliquevarpos[i]]] - row2coefs[varinds[cliquevarpos[i]]];
   }

   SCIPsortRealInt(gradients, cliquevarpos, cliquesize);

#ifdef SCIP_DEBUG_CLIQUE
   for( i = 0; i < cliquesize; i++ )
      SCIPdebugMsg(scip, "var_%d: %g + %g * lambda\n", i, row2coefs[varinds[cliquevarpos[i]]], gradients[i]);
#endif


   // find maximum for lambda = 0
   maxpos = 0;
   for( i = 1; i < cliquesize; i++ )
   {
      if( SCIPisGE(scip, row2coefs[varinds[cliquevarpos[i]]], row2coefs[varinds[cliquevarpos[maxpos]]]) )
         maxpos = i;
   }

   // variable is relevant only if its coefficient is non-negative, if all coefs are negative we set firstmaxpos to -1
   if( SCIPisPositive(scip, row2coefs[varinds[cliquevarpos[maxpos]]]) || (SCIPisZero(scip, row2coefs[varinds[cliquevarpos[maxpos]]]) && SCIPisPositive(scip, gradients[maxpos])) )
      firstmaxpos = maxpos;
   else
      firstmaxpos = -1;

#ifdef SCIP_DEBUG_CLIQUE
   SCIPdebugMsg(scip, "first maxvar: var_%d at coef(0) = %g, firstmaxpos = %d\n", maxpos, row2coefs[varinds[cliquevarpos[maxpos]]], firstmaxpos);
#endif

   // find next maximum
   minlambda = 0.0;
   if( firstmaxpos == -1 )
      maxidx = -1;
   else
      maxidx = varinds[cliquevarpos[firstmaxpos]];

   while( maxpos < cliquesize && SCIPisLT(scip, minlambda, 1.0) )
   {
      newmaxpos = -1;

      // all coefficients are negative
      if( maxidx == -1 )
      {
         // find variable with lowest non-negative lambda such that (row2coef - row1coef) * lambda + row2coef >= 0 */
         minlambda = 2.0;
         for( i = maxpos; i < cliquesize; i++ )
         {
            if( SCIPisPositive(scip, gradients[i]) )
            {
               idx = varinds[cliquevarpos[i]];
               lambda = -row2coefs[idx] / gradients[i];
               if( SCIPisLE(scip, lambda, minlambda) )
               {
                  minlambda = lambda;
                  newmaxpos = i;
               }
            }
         }
      }
      else
      {
         minlambda = 2.0;
         idx = varinds[cliquevarpos[maxpos]];
         for( i = maxpos + 1; i < cliquesize; i++ )
         {
            if( !SCIPisEQ(scip, gradients[i], gradients[maxpos]) )
            {
               lambda = (row2coefs[idx] - row2coefs[varinds[cliquevarpos[i]]]) / (gradients[i] - gradients[maxpos]);
               if( SCIPisLE(scip, lambda, minlambda) )
               {
                  minlambda = lambda;
                  newmaxpos = i;
               }
            }
         }
      }

      if( newmaxpos == -1 || SCIPisGE(scip, minlambda, 1.0) )
      {
         /* check whether last segment becomes negative */
         if( maxidx != -1 && SCIPisNegative(scip, row1coefs[varinds[cliquevarpos[maxpos]]]) )
         {
            assert(SCIPisPositive(scip, row2coefs[varinds[cliquevarpos[maxpos]]]));
            assert(cliquemaxinds[cliquevarpos[firstmaxpos]] < 0);
            assert(firstmaxpos != -1); // cheesy use of firstmaxpos as breakpoint dummy variable for the negative segment
            breakpoints[cliquevarpos[firstmaxpos]] = -row2coefs[maxidx] / gradients[maxpos]; // lambda where old maxcoef hits zero
            maxidx = -1;
            (*nbreakpoints)++;
         }
         else
            break;
      }
      else
      {
         breakpointval = row2coefs[varinds[cliquevarpos[newmaxpos]]] + minlambda * gradients[newmaxpos];
         assert(maxidx == -1 || SCIPisEQ(scip, breakpointval, row2coefs[maxidx] + minlambda * gradients[maxpos]));

         /* check if next segment can become negative */
         if( SCIPisNegative(scip, breakpointval) || (SCIPisZero(scip, breakpointval) && !SCIPisPositive(scip, gradients[newmaxpos])) )
         {
            assert(cliquemaxinds[cliquevarpos[firstmaxpos]] < 0);
            assert(firstmaxpos != -1); // cheesy use of firstmaxpos as breakpoint dummy variable for the negative segment
            breakpoints[cliquevarpos[firstmaxpos]] = -row2coefs[maxidx] / gradients[maxpos]; // lambda where old maxcoef hits zero
            maxidx = -1;
         }
         else
         {
            maxidx = varinds[cliquevarpos[newmaxpos]];
            breakpoints[cliquevarpos[newmaxpos]] = minlambda;
            cliquemaxinds[cliquevarpos[newmaxpos]] = -cliquemaxinds[cliquevarpos[newmaxpos]];
         }

         maxpos = newmaxpos;
         (*nbreakpoints)++;

#ifdef SCIP_DEBUG_CLIQUE
         if( maxidx == -1)
            SCIPdebugMsg(scip, "next maxvar: firstmaxvar as dummy for negative clique, cliquemaxinds[%d] = %d\n", varinds[cliquevarpos[firstmaxpos]], cliquemaxinds[varinds[cliquevarpos[firstmaxpos]]]);
         else
            SCIPdebugMsg(scip, "next maxvar: var_%d at coef(%g) = %g, cliquemaxinds[%d] = %d, \n", maxpos, minlambda, breakpointval, maxidx, cliquemaxinds[maxidx]);
#endif
      }
   }

   SCIPfreeBufferArray(scip, &gradients);

   if( firstmaxpos == -1 )
      return -1;
   else
      return varinds[cliquevarpos[firstmaxpos]];
}

/** try two-row combine for given rows */
static
SCIP_RETCODE combineRows
(
   SCIP*                scip,                /**< SCIP datastructure */
   SCIP_MATRIX*         matrix,              /**< the constraint matrix */
   int                  row1,
   int                  row2,
   SCIP_Bool            swaprow1,
   SCIP_Bool            swaprow2,
   SCIP_Real*           lbs,
   SCIP_Real*           ubs,
   SCIP_Bool*           success             /**< we return (success ||  found better bounds") */
)
{
   int i;
   int j;
   int ncols;
   int nvars;
   int* varinds;
   int nbreakpoints;
   SCIP_Real* breakpoints;
   int* cliquemaxinds;
   int idx;
   int idx1;
   int idx2;
   int row1len;
   int row2len;
   int* row1idxptr;
   int* row2idxptr;
   SCIP_Real* row1valptr;
   SCIP_Real* row2valptr;
   SCIP_Real* row1coefs;
   SCIP_Real* row2coefs;
   enum signum* signs;
   SCIP_Real  b1;
   SCIP_Real  b2;
   int ninfs;
   int l1infs;
   SCIP_Real  l1;
   SCIP_Real  l2;
   SCIP_Real* newlbs;
   SCIP_Real* newubs;
   SCIP_Real coef;
   int sign;
   SCIP_VAR* var;

   int shift;
   int nbinvars;
   SCIP_VAR** binvars;
   int* binvarpos;
   int* cliquepartition;
   int ncliques;
   int* currentmaxinds;
   int cliqueidx;

#ifdef SCIP_DEBUG_SUBSCIP
   SCIP* subscip;
   SCIP_VAR** subvars;
   SCIP_CONS* subcons1;
   SCIP_CONS* subcons2;
   SCIP_Real* subrow1coefs;
   SCIP_Real* subrow2coefs;
   SCIP_Real* subsciplbs;
   SCIP_Real* subscipubs;
#endif

   //TODO Check if I can ask row if its variable bounds are tightened with standard bound tightening
#ifdef SCIP_MORE_DEBUG
   SCIPdebugMsg(scip, "combining rows %d (%s) and %d (%s) with swaps %d/%d\n", row1, SCIPmatrixGetRowName(matrix, row1), row2, SCIPmatrixGetRowName(matrix, row2), swaprow1, swaprow2);
#endif

   ncols = SCIPmatrixGetNColumns(matrix);

   row1len = SCIPmatrixGetRowNNonzs(matrix, row1);
   row1idxptr = SCIPmatrixGetRowIdxPtr(matrix, row1);
   row1valptr = SCIPmatrixGetRowValPtr(matrix, row1);

   row2len = SCIPmatrixGetRowNNonzs(matrix, row2);
   row2idxptr = SCIPmatrixGetRowIdxPtr(matrix, row2);
   row2valptr = SCIPmatrixGetRowValPtr(matrix, row2);

   SCIP_CALL( SCIPallocBufferArray(scip, &row1coefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &row2coefs, ncols) );

   SCIPsortIntReal(row1idxptr, row1valptr, row1len);
   SCIPsortIntReal(row2idxptr, row2valptr, row2len);

   /* swap rows if necessary */
   if( swaprow1 )
   {
      for( i = 0; i < row1len; i++ )
         row1coefs[row1idxptr[i]] = -row1valptr[i];

      assert(!SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, row1)));
      b1 = -SCIPmatrixGetRowRhs(matrix, row1);
   }
   else
   {
      for( i = 0; i < row1len; i++ )
         row1coefs[row1idxptr[i]] = row1valptr[i];

      b1 = SCIPmatrixGetRowLhs(matrix, row1);
   }

   if( swaprow2 )
   {
      for( i = 0; i < row2len; i++ )
         row2coefs[row2idxptr[i]] = -row2valptr[i];

      assert(!SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, row2)));
      b2 = -SCIPmatrixGetRowRhs(matrix, row2);
   }
   else
   {
      for( i = 0; i < row2len; i++ )
         row2coefs[row2idxptr[i]] = row2valptr[i];

      b2 = SCIPmatrixGetRowLhs(matrix, row2);
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &signs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varinds, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &breakpoints, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cliquemaxinds, ncols) );

   /* calculate cancellation breakpoints and sign behaviour of non-binary variables */
   i = 0;
   j = 0;
   nvars = 0;
   nbreakpoints = 0;
   while( i < row1len && j < row2len )
   {
      assert(i+1 == row1len || row1idxptr[i] < row1idxptr[i+1]);
      assert(j+1 == row2len || row2idxptr[j] < row2idxptr[j+1]);

      idx1 = row1idxptr[i];
      idx2 = row2idxptr[j];

      /* We use 2.0 as default value for "no cancellation". For cancellations, this will be replaced by values in (0,1).
       * A value larger than 1.0 is used because we sort the array and want to put non-cancellations at the end. */
      breakpoints[nvars] = 2.0;
      // We use 0 as default value for "not in any clique"
      cliquemaxinds[nvars] = 0;

      if( idx1 == idx2 )
      {
         if( (SCIPisNegative(scip, row1coefs[idx1]) && SCIPisPositive(scip, row2coefs[idx2])) ||
             (SCIPisPositive(scip, row1coefs[idx1]) && SCIPisNegative(scip, row2coefs[idx2])) )
         {
            if( SCIPisNegative(scip, row2coefs[idx2]) )
               signs[idx1] = UP;
            else
               signs[idx1] = DN;

            breakpoints[nvars] = row2coefs[idx2]/(row2coefs[idx2] - row1coefs[idx1]);
            nbreakpoints++;
         }
         else if( SCIPisPositive(scip, row1coefs[idx1]) )
            signs[idx1] = POS;
         else
            signs[idx1] = NEG;

         varinds[nvars] = idx1;
         i++;
         j++;
      }
      else if( idx1 < idx2 )
      {
         if( SCIPisPositive(scip, row1coefs[idx1]) )
            signs[idx1] = POS;
         else
            signs[idx1] = NEG;

         /* We will access this entry later on, so we explicitly write a zero here */
         row2coefs[idx1] = 0.0;

         varinds[nvars] = idx1;
         i++;
      }
      else
      {
         assert(idx1 > idx2);
         if( SCIPisPositive(scip, row2coefs[idx2]) )
            signs[idx2] = POS;
         else
            signs[idx2] = NEG;

         /* We will access this entry later on, so we explicitly write a zero here */
         row1coefs[idx2] = 0.0;

         varinds[nvars] = idx2;
         j++;
      }
      nvars++;
   }

   while( i < row1len )
   {
      idx1 = row1idxptr[i];

      if( SCIPisPositive(scip, row1coefs[idx1]) )
         signs[idx1] = POS;
      else
         signs[idx1] = NEG;

      /* We will access this entry later on, so we explicitly write a zero here */
      row2coefs[idx1] = 0.0;

      varinds[nvars] = idx1;
      breakpoints[nvars] = 2.0;
      cliquemaxinds[nvars] = 0;

      nvars++;
      i++;
   }

   while( j < row2len )
   {
      idx2 = row2idxptr[j];

      if( SCIPisPositive(scip, row2coefs[idx2]) )
         signs[idx2] = POS;
      else
         signs[idx2] = NEG;

      /* We will access this entry later on, so we explicitly write a zero here */
      row1coefs[idx2] = 0.0;

      varinds[nvars] = idx2;
      breakpoints[nvars] = 2.0;
      cliquemaxinds[nvars] = 0;

      nvars++;
      j++;
   }

   /* calculate clique breakpoints */
   SCIP_CALL( SCIPallocBufferArray(scip, &binvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvarpos, nvars) );

   nbinvars = 0;
   for( i = 0; i < nvars; i++ )
   {
      var = SCIPmatrixGetVar(matrix, varinds[i]);
      if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
      {
         binvars[nbinvars] = var;
         binvarpos[nbinvars] = i;
         nbinvars++;
      }
#ifdef SCIP_DEBUG_SIGNS
      idx = varinds[i];
      if( signs[idx] == UP )
         assert(SCIPisGT(scip, row1coefs[idx], 0.0) && SCIPisLT(scip, row2coefs[idx], 0.0));
      else if( signs[idx] == DN )
         assert(SCIPisLT(scip, row1coefs[idx], 0.0) && SCIPisGT(scip, row2coefs[idx], 0.0));
      else if( signs[idx] == POS )
         assert(SCIPisGE(scip, row1coefs[idx], 0.0) && SCIPisGE(scip, row2coefs[idx], 0.0));
      else if( signs[idx] == NEG )
         assert(SCIPisLE(scip, row1coefs[idx], 0.0) && SCIPisLE(scip, row2coefs[idx], 0.0));
#endif
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &cliquepartition, nbinvars) );

#ifdef SCIP_TEST_CLIQUE
   if( nbinvars > 0 )
   {
      SCIP_Bool infeasible;
      int nbdchgs;
      SCIPaddClique(scip, binvars, NULL, nbinvars, FALSE, &infeasible, &nbdchgs);
   }
#endif

   SCIPcalcCliquePartition(scip, binvars, nbinvars, cliquepartition, &ncliques);
   SCIPsortIntInt(cliquepartition, binvarpos, nbinvars);

#ifdef SCIP_DEBUG_CLIQUE
   SCIPdebugMsg(scip, "ncliques = %d, nbinvars = %d\n", ncliques, nbinvars);
   for( i = 0; i < nbinvars; i++ )
      SCIPdebugMsg(scip, "%s (%d) is in clique %d\n", SCIPmatrixGetColName(matrix, varinds[binvarpos[i]]), varinds[binvarpos[i]], cliquepartition[i]);

   SCIPdebugMsg(scip, "breakpoints without cliques: %d\n", nbreakpoints);
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &currentmaxinds, ncliques+1) ); // The entry at index 0 is a dummy to prevent excessive index shifting in the code

   /* remove cliques containing a single variable and compute clique maximums */
   shift = 0;
   for( i = 0; i < nbinvars; )
   {
#ifdef SCIP_DEBUG_CLIQUE
      SCIPdebugMsg(scip, "checking clique of var_%d\n", i);
#endif
      if( i + 1 < nbinvars && cliquepartition[i] == cliquepartition[i+1] )
      {
         int currentclique = cliquepartition[i];
         for( j = i; j < nbinvars && currentclique == cliquepartition[j]; j++ )
         {
            signs[varinds[binvarpos[j]]] = CLQ;
            binvarpos[j-shift] = binvarpos[j];
            cliquepartition[j-shift] = cliquepartition[j] - shift + 1; // +1 ensures that all clique IDs are >0 such that we can use -ID to signal an all-negative clique
            cliquemaxinds[binvarpos[j]] = -cliquepartition[j-shift]; // negative sign implies that this variable does not assume the maximum in this clique, will be adjusted in calcCliqueMaximums
            // variables in cliques have either a clique-breakpoint or no breakpoint at all
            if( !SCIPisEQ(scip, breakpoints[binvarpos[j]], 2.0) )
            {
               breakpoints[binvarpos[j]] = 2.0;
               nbreakpoints--;
            }
         }

#ifdef SCIP_DEBUG_CLIQUE
         SCIPdebugMsg(scip, "breakpoints before checking maximums of clique %d: %d\n", cliquepartition[binvarpos[i]], nbreakpoints);
#endif
         // size of current clique equals j - i
         assert((shift + j) <= ncols);
         idx = cliquepartition[i-shift]; // we enumerated the relevant cliques only
         currentmaxinds[idx] = calcCliqueMaximums(scip, varinds, &binvarpos[i-shift], j - i, row1coefs, row2coefs, signs, &nbreakpoints, breakpoints, cliquemaxinds);
#ifdef SCIP_DEBUG_CLIQUE
         SCIPdebugMsg(scip, "breakpoints after checking maximums of clique %d: %d, firstmaxidx = %d\n", idx, nbreakpoints, currentmaxinds[idx]);
#endif

         i = j;
      }
      else
      {
         shift++;
         i++;
      }
   }
   ncliques -= shift;
   nbinvars -= shift;

#ifdef SCIP_DEBUG_CLIQUE
   SCIPdebugMsg(scip, "after clearing single variable cliques: ncliques = %d, nbinvars = %d, nbreakpoints = %d\n", ncliques, nbinvars, nbreakpoints);
   for( i = 0; i < nbinvars; i++ )
      SCIPdebugMsg(scip, "%s (%d) is in clique %d with breakpoint %g\n", SCIPmatrixGetColName(matrix, varinds[binvarpos[i]]), varinds[binvarpos[i]], cliquepartition[i], breakpoints[binvarpos[i]]);
#endif

   /* The obvious preconditions for bound tightenings are met, so we try to calculate new bounds. */
   if( nbreakpoints >= 1 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &newlbs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newubs, nvars) );
#ifdef SCIP_DEBUG_SUBSCIP
      SCIP_CALL( SCIPallocBufferArray(scip, &subsciplbs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipubs, nvars) );
#endif

      SCIPsortRealIntInt(breakpoints, varinds, cliquemaxinds, nvars);

      for( i = 0; i < nvars; i++)
      {
         idx = varinds[i];
         newlbs[i] = lbs[idx];
         newubs[i] = ubs[idx];
      }

#ifdef SCIP_DEBUG_SIGNS
      SCIPdebugMsg(scip, "b1 = %g, b2 = %g, nbreakpoints = %d\n", b1, b2, nbreakpoints);
      for( i = 0; i < nvars; i++ )
      {
         idx = varinds[i];
         if( signs[idx] == UP )
            SCIPdebugMsg(scip, "%g <= %s <= %g, UP: %g to %g, breakpoint at %g\n", newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
         else if( signs[idx] == DN )
            SCIPdebugMsg(scip, "%g <= %s <= %g, DN: %g to %g, breakpoint at %g\n", newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
         else if( signs[idx] == POS )
            SCIPdebugMsg(scip, "%g <= %s <= %g, POS: %g to %g, breakpoint at %g\n", newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
         else if( signs[idx] == NEG )
            SCIPdebugMsg(scip, "%g <= %s <= %g, NEG: %g to %g, breakpoint at %g\n", newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
         else if( signs[idx] == CLQ )
            SCIPdebugMsg(scip, "%g <= %s <= %g, CLQ: %g to %g, breakpoint at %g, cliquemaxind = %d\n", newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i], cliquemaxinds[i]);
      }
#endif

#ifdef SCIP_DEBUG_CLIQUE
      SCIPdebugMsg(scip, "currentmaxinds[0] = dummy\n");
      for( i = 1; i < ncliques+1; i++ )
         if( currentmaxinds[i] >= 0 )
            SCIPdebugMsg(scip, "currentmaxinds[%d] = %d (%s)\n", i, currentmaxinds[i], SCIPmatrixGetColName(matrix, currentmaxinds[i]));
         else
            SCIPdebugMsg(scip, "currentmaxinds[%d] = %d (no maximal variable)\n", i, currentmaxinds[i]);
#endif

      /* calculate activity contributions of each row */
      l1 = b1;
      l2 = b2;
      l1infs = 0;
      ninfs = 0;
      for( i = 0; i < nvars; i++ )
      {
         idx = varinds[i];
         if( signs[idx] == UP || signs[idx] == NEG )
         {
            assert(SCIPisNegative(scip, row2coefs[idx]) || (SCIPisZero(scip, row2coefs[idx]) && SCIPisNegative(scip, row1coefs[idx])));
            if( !SCIPisInfinity(scip, -lbs[idx]) )
            {
               l1 -= row1coefs[idx] * lbs[idx];
               l2 -= row2coefs[idx] * lbs[idx];
            }
            else if ( SCIPisZero(scip, row2coefs[idx]) )
               l1infs++;
            else
               ninfs++;
         }
         else if( signs[idx] == DN || signs[idx] == POS )
         {
            assert(SCIPisPositive(scip, row2coefs[idx]) || (SCIPisZero(scip, row2coefs[idx]) && SCIPisPositive(scip, row1coefs[idx])));
            if( !SCIPisInfinity(scip, ubs[idx]) )
            {
               l1 -= row1coefs[idx] * ubs[idx];
               l2 -= row2coefs[idx] * ubs[idx];
            }
            else if( SCIPisZero(scip, row2coefs[idx]) )
               l1infs++;
            else
               ninfs++;
         }
      }

      /* calculate activity contributions of cliques, ignore dummy entry at index 0 */
      for( i = 1; i < ncliques+1; i++ )
      {
         if( currentmaxinds[i] != -1 )
         {
            idx = currentmaxinds[i];
            assert(!SCIPisNegative(scip, row2coefs[idx]));
            l1 -= row1coefs[idx];
            l2 -= row2coefs[idx];
         }
      }

      /* TODO Implement a check that skips all further computation if ninfs >= maxcancel + 2
       * where maxcancel is the maximum number of variables that can be cancelled at once */

#ifdef SCIP_DEBUG_SUBSCIP
      SCIPdebugMsg(scip, "calculating bounds via subscip\n");
      /* calculate bounds by explicitly solving the LP-relaxation of the two constraints */
      SCIP_CALL( SCIPcreate(&subscip) );
      SCIP_CALL( SCIPcreateProbBasic(subscip, "subscip") );
      SCIPsetMessagehdlrQuiet(subscip, TRUE);
      SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subvars, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subrow1coefs, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subrow2coefs, nvars) );
      SCIP_CALL( SCIPsetIntParam(subscip, "presolving/tworowcomb/maxrounds", 0) );

      for( i = 0; i < nvars; i++ )
      {
         SCIP_Bool dummy;

         idx = varinds[i];
         var = SCIPmatrixGetVar(matrix, idx);

         assert(SCIPisLE(scip, lbs[idx], ubs[idx]));

         subrow1coefs[i] = row1coefs[idx];
         subrow2coefs[i] = row2coefs[idx];

         SCIP_CALL( SCIPgetVarCopy(scip, subscip, var, &subvars[i], NULL, NULL, FALSE, &dummy) );
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, subvars[i], lbs[idx]) );
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, subvars[i], ubs[idx]) );
         SCIPchgVarType(subscip, subvars[i], SCIP_VARTYPE_CONTINUOUS, &dummy);
         SCIPchgVarObj(subscip, subvars[i], 0.0);
      }

      SCIP_CALL( SCIPcreateConsLinear(subscip, &subcons1, SCIPmatrixGetRowName(matrix, row1), nvars, subvars,
         subrow1coefs, b1, SCIPinfinity(subscip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, subcons1) );
      SCIP_CALL( SCIPreleaseCons(subscip, &subcons1) );

      SCIP_CALL( SCIPcreateConsLinear(subscip, &subcons2, SCIPmatrixGetRowName(matrix, row2), nvars, subvars,
         subrow2coefs, b2, SCIPinfinity(subscip), TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(subscip, subcons2) );
      SCIP_CALL( SCIPreleaseCons(subscip, &subcons2) );

      //SCIP_CALL( SCIPwriteOrigProblem(subscip, "../subscip.lp", "lp", FALSE) );

      for( i = 0; i < nvars; i++ )
      {
         SCIPchgVarObj(subscip, subvars[i], 1.0);

         SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE);
         SCIPsolve(subscip);
         subscipubs[i] = SCIPgetPrimalbound(subscip);
         SCIP_CALL( SCIPfreeTransform(subscip) );

         SCIPsetObjsense(subscip, SCIP_OBJSENSE_MINIMIZE);
         SCIPsolve(subscip);
         subsciplbs[i] = SCIPgetPrimalbound(subscip);
         SCIP_CALL( SCIPfreeTransform(subscip) );

         SCIPchgVarObj(subscip, subvars[i], 0.0);
      }

      SCIPfreeBufferArray(scip, &subrow2coefs);
      SCIPfreeBufferArray(scip, &subrow1coefs);
      SCIPfreeBufferArray(scip, &subvars);
      SCIP_CALL( SCIPfree(&subscip) );
      SCIPdebugMsg(scip, "finished calculating bounds via subscip\n");
#endif
      /* Calculate bounds for lambda = 0 */
#ifdef SCIP_DEBUG_BREAKPOINTS
      SCIPdebugMsg(scip, "lambda = 0, l1 = %g, l2 = %g, ninfs = %d\n", l1, l2, ninfs);
#endif

      /* try strengthening the bound of the one variable which adds the infinity */
      if( ninfs == 1 )
      {
#ifdef SCIP_DEBUG_BOUNDS
         SCIP_Real oldlb;
         SCIP_Real oldub;
#endif
         for( i = 0; i < nvars; i++ )
         {
#ifdef SCIP_DEBUG_BOUNDS
            oldlb = newlbs[i];
            oldub = newubs[i];
#endif
            idx = varinds[i];
            if( SCIPisPositive(scip, row2coefs[idx]) && SCIPisInfinity(scip, ubs[idx]) )
               newlbs[i] = MAX(newlbs[i], l2 / row2coefs[idx]);
            else if ( SCIPisNegative(scip, row2coefs[idx]) && SCIPisInfinity(scip, -lbs[idx]) )
               newubs[i] = MIN(newubs[i], l2 / row2coefs[idx]);
#ifdef SCIP_DEBUG_BOUNDS
            if( !SCIPisEQ(scip, oldlb, newlbs[i]) || !SCIPisEQ(scip, oldub, newubs[i]) )
               SCIPdebugMsg(scip, "%g <= %g <= %s <= %g <= %g\n", oldlb, newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], oldub);
#endif
         }
      }

      /* try strenghtening all bounds */
      else if( ninfs == 0 )
      {
#ifdef SCIP_DEBUG_BOUNDS
         SCIP_Real oldlb;
         SCIP_Real oldub;
#endif
         for( i = 0; i < nvars; i++ )
         {
#ifdef SCIP_DEBUG_BOUNDS
            oldlb = newlbs[i];
            oldub = newubs[i];
#endif
            idx = varinds[i];
            if( SCIPisPositive(scip, row2coefs[idx]) )
            {
               if( signs[idx] == CLQ )
                  newlbs[i] = MAX(newlbs[i], l2 / row2coefs[idx]);
               else
                  newlbs[i] = MAX(newlbs[i], (l2 + row2coefs[idx] * ubs[idx]) / row2coefs[idx]);
            }
            else if( SCIPisNegative(scip, row2coefs[idx]) )
            {
               if( signs[idx] == CLQ )
               {
                  if( cliquemaxinds[i] > 0 )
                     cliqueidx = cliquemaxinds[i];
                  else
                     cliqueidx = -cliquemaxinds[i];

                  if(currentmaxinds[cliqueidx] == -1)
                     newubs[i] = MIN(newubs[i], l2 / row2coefs[idx]);
                  else
                  {
                     // never leads to division by zero as negative coefficients can never be equal to the cliquemaximum
                     assert(!SCIPisZero(scip, row2coefs[idx] - row2coefs[currentmaxinds[cliqueidx]]));
                     newubs[i] = MIN(newubs[i], l2 / (row2coefs[idx] - row2coefs[currentmaxinds[cliqueidx]]));
                  }
               }
               else
                  newubs[i] = MIN(newubs[i], (l2 + row2coefs[idx] * lbs[idx]) / row2coefs[idx]);
            }
#ifdef SCIP_DEBUG_BOUNDS
            if( !SCIPisEQ(scip, oldlb, newlbs[i]) || !SCIPisEQ(scip, oldub, newubs[i]) )
               SCIPdebugMsg(scip, "%g <= %g <= %s <= %g <= %g\n", oldlb, newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], oldub);
#endif
         }
      }
      ninfs += l1infs;

      i = 0;
      while( i < nbreakpoints )
      {
         int nnewinfs;
         SCIP_Real l1update;
         SCIP_Real l2update;
         SCIP_Bool updated;

         assert(SCIPisPositive(scip, breakpoints[i]));
         assert(SCIPisLT(scip, breakpoints[i], 1.0));

         /* determine number of infinities and compute update for l1 and l2 */
         shift = 0;
         nnewinfs = 0;
         l1update = 0.0;
         l2update = 0.0;
         updated = FALSE;
         for( j = i; !updated; j++ )
         {
            idx = varinds[j];
            assert(signs[idx] == UP || signs[idx] == DN || signs[idx] == CLQ);
            if( signs[idx] == CLQ )
            {
               if( cliquemaxinds[j] < 0 )
               {
                  assert(currentmaxinds[-cliquemaxinds[j]] != -1 );
                  l1update += row1coefs[currentmaxinds[-cliquemaxinds[j]]];
                  l2update += row2coefs[currentmaxinds[-cliquemaxinds[j]]];
                  currentmaxinds[-cliquemaxinds[j]] = -1;
               }
               else if( currentmaxinds[cliquemaxinds[j]] == -1 )
               {
                  assert(SCIPisPositive(scip, row1coefs[idx]));
                  assert(SCIPisNegative(scip, row2coefs[idx]));
                  l1update -= row1coefs[idx];
                  l2update -= row2coefs[idx];
                  currentmaxinds[cliquemaxinds[j]] = idx;
               }
               else
               {
                  l1update += row1coefs[currentmaxinds[cliquemaxinds[j]]] - row1coefs[idx];
                  l2update += row2coefs[currentmaxinds[cliquemaxinds[j]]] - row2coefs[idx];
                  currentmaxinds[cliquemaxinds[j]] = idx;
               }
            }
            else
            {
               if( signs[idx] == UP )
                  sign = 1;
               else
                  sign = -1;

               if( !SCIPisInfinity(scip, -lbs[idx]) )
               {
                  l1update += sign * row1coefs[idx] * lbs[idx];
                  l2update += sign * row2coefs[idx] * lbs[idx];
               }
               else
               {
                  if( signs[idx] == UP  )
                     ninfs--;
                  else
                     nnewinfs++;
               }

               if( !SCIPisInfinity(scip, ubs[idx]) )
               {
                  l1update -= sign * row1coefs[idx] * ubs[idx];
                  l2update -= sign * row2coefs[idx] * ubs[idx];
               }
               else
               {
                  if( signs[idx] == UP  )
                     nnewinfs++;
                  else
                     ninfs--;
               }

               if( signs[idx] == UP )
                  signs[idx] = POS;
               else if( signs[idx] == DN )
                  signs[idx] = NEG;
            }

            if( j+1 >= nbreakpoints || !SCIPisEQ(scip, breakpoints[j], breakpoints[j+1]) )
               updated = TRUE;

            shift++;
         }

#ifdef SCIP_DEBUG_BREAKPOINTS
         SCIPdebugMsg(scip, "lambda_%d = %g, l1 = %g, l2 = %g, ninfs = %d\n", i, breakpoints[i], l1, l2, ninfs);
#endif

         assert(ninfs >= 0);

         /* if more than one infinity destroys our bounds we cannot tighten anything */
         if( ninfs <= 1 )
         {
            /* check for bounds to be tightened */
            for( j = 0; j < nvars; j++ )
            {
#ifdef SCIP_DEBUG_BOUNDS
               SCIP_Real oldlb;
               SCIP_Real oldub;
#endif
               /* skip the next variables as they are cancelled anyway */
               if( j == i )
                  j += shift;

               /* catches the special case where the entire remaining constraint is cancelled */
               if( j >= nvars )
                  break;

#ifdef SCIP_DEBUG_BOUNDS
               oldlb = newlbs[j];
               oldub = newubs[j];
#endif

               idx = varinds[j];
               coef = breakpoints[i] * row1coefs[idx] + (1 - breakpoints[i]) * row2coefs[idx];
               assert(!SCIPisEQ(scip, breakpoints[i], 2.0));

               // skip if the coefficient is too close to zero as it becomes numerically unstable
               if( SCIPisZero(scip, coef) )
                  continue;

               if( signs[idx] == POS || signs[idx] == DN )
               {
                  if( ninfs == 0 )
                     newlbs[j] = MAX(newlbs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2 + coef * ubs[idx]) / coef);
                  else if( SCIPisInfinity(scip, ubs[idx]) )
                     newlbs[j] = MAX(newlbs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2) / coef);
               }
               else if ( signs[idx] == NEG || signs[idx] == UP )
               {
                  if( ninfs == 0 )
                     newubs[j] = MIN(newubs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2 + coef * lbs[idx]) / coef);
                  else if( SCIPisInfinity(scip, -lbs[idx]) )
                     newubs[j] = MIN(newubs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2) / coef);
               }
               else if( ninfs == 0 )
               {
                  assert(signs[idx] == CLQ);
                  if( SCIPisNegative(scip, coef) )
                  {
                     if( cliquemaxinds[j] > 0 )
                        cliqueidx = cliquemaxinds[j];
                     else
                        cliqueidx = -cliquemaxinds[j];

                     if(currentmaxinds[cliqueidx] == -1)
                        newubs[j] = MIN(newubs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2) / coef);
                     else
                        newubs[j] = MIN(newubs[j], (breakpoints[i] * (l1 + row1coefs[currentmaxinds[cliqueidx]]) + (1 - breakpoints[i]) * (l2 + row2coefs[currentmaxinds[cliqueidx]])) / coef);
                  }
                  else
                     newlbs[j] = MAX(newlbs[j], (breakpoints[i] * l1 + (1 - breakpoints[i]) * l2) / coef);
               }

#ifdef SCIP_DEBUG_BOUNDS
               if( !SCIPisEQ(scip, oldlb, newlbs[j]) || !SCIPisEQ(scip, oldub, newubs[j]) )
                  SCIPdebugMsg(scip, "%g <= %g <= %s <= %g <= %g\n", oldlb, newlbs[j], SCIPmatrixGetColName(matrix, idx), newubs[j], oldub);
#endif
            }
         }

#ifdef SCIP_DEBUG_CLIQUE
      SCIPdebugMsg(scip, "currentmaxinds[0] = dummy\n");
      for( j = 1; j < ncliques+1; j++ )
         if( currentmaxinds[j] >= 0 )
            SCIPdebugMsg(scip, "currentmaxinds[%d] = %d (%s)\n", j, currentmaxinds[j], SCIPmatrixGetColName(matrix, currentmaxinds[j]));
         else
            SCIPdebugMsg(scip, "currentmaxinds[%d] = %d (no maximal variable)\n", j, currentmaxinds[j]);
#endif

         i += shift;
         ninfs += nnewinfs;
         l1 += l1update;
         l2 += l2update;
      }

      /* check infinities in first row */
      ninfs = 0;
      for( i = 0; i < nvars; i++ )
      {
         idx = varinds[i];
         if( (SCIPisPositive(scip, row1coefs[idx]) && SCIPisInfinity(scip, ubs[idx]))
             || (SCIPisNegative(scip, row1coefs[idx]) && SCIPisInfinity(scip, -lbs[idx])) )
            ninfs++;
      }

      /* calculate bounds for lambda = 1 */
#ifdef SCIP_DEBUG_BREAKPOINTS
      SCIPdebugMsg(scip, "lambda = 1, l1 = %g, l2 = %g, ninfs = %d\n", l1, l2, ninfs);
#endif

      /* try strengthening the bound of the one variable which adds the infinity */
      if( ninfs == 1 )
      {
#ifdef SCIP_DEBUG_BOUNDS
         SCIP_Real oldlb;
         SCIP_Real oldub;
#endif
         for( i = 0; i < nvars; i++ )
         {
#ifdef SCIP_DEBUG_BOUNDS
            oldlb = newlbs[i];
            oldub = newubs[i];
#endif
            idx = varinds[i];
            if( SCIPisPositive(scip, row1coefs[idx]) && SCIPisInfinity(scip, ubs[idx]) )
               newlbs[i] = MAX(newlbs[i], l1 / row1coefs[idx]);
            else if ( SCIPisNegative(scip, row1coefs[idx]) && SCIPisInfinity(scip, -lbs[idx]) )
               newubs[i] = MIN(newubs[i], l1 / row1coefs[idx]);
#ifdef SCIP_DEBUG_BOUNDS
            if( !SCIPisEQ(scip, oldlb, newlbs[i]) || !SCIPisEQ(scip, oldub, newubs[i]) )
               SCIPdebugMsg(scip, "%g <= %g <= %s <= %g <= %g\n", oldlb, newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], oldub);
#endif
         }
      }
      /* try strengthening all bounds */
      else if( ninfs == 0 )
      {
#ifdef SCIP_DEBUG_BOUNDS
         SCIP_Real oldlb;
         SCIP_Real oldub;
#endif
         for( i = 0; i < nvars; i++ )
         {
#ifdef SCIP_DEBUG_BOUNDS
            oldlb = newlbs[i];
            oldub = newubs[i];
#endif
            idx = varinds[i];
            if( SCIPisPositive(scip, row1coefs[idx]) )
            {
               if( signs[idx] == CLQ )
                  newlbs[i] = MAX(newlbs[i], l1 / row1coefs[idx]);
               else
                  newlbs[i] = MAX(newlbs[i], (l1 + row1coefs[idx] * ubs[idx]) / row1coefs[idx]);
            }
            else if( SCIPisNegative(scip, row1coefs[idx]) )
            {
               if( signs[idx] == CLQ )
               {
                  if( cliquemaxinds[i] > 0 )
                     cliqueidx = cliquemaxinds[i];
                  else
                     cliqueidx = -cliquemaxinds[i];

                  if(currentmaxinds[cliqueidx] == -1)
                     newubs[i] = MIN(newubs[i], l1 / row1coefs[idx]);
                  else
                  {
                     // never leads to division by zero as negative coefficients can never be equal to the cliquemaximum
                     assert(!SCIPisZero(scip, row1coefs[idx] - row1coefs[currentmaxinds[cliqueidx]]));
                     newubs[i] = MIN(newubs[i], l1 / (row1coefs[idx] - row1coefs[currentmaxinds[cliqueidx]]));
                  }
               }
               else
                  newubs[i] = MIN(newubs[i], (l1 + row1coefs[idx] * lbs[idx]) / row1coefs[idx]);
            }
#ifdef SCIP_DEBUG_BOUNDS
            if( !SCIPisEQ(scip, oldlb, newlbs[i]) || !SCIPisEQ(scip, oldub, newubs[i]) )
               SCIPdebugMsg(scip, "%g <= %g <= %s <= %g <= %g\n", oldlb, newlbs[i], SCIPmatrixGetColName(matrix, idx), newubs[i], oldub);
#endif
         }
      }

      /* update bound arrays and determine success */
      for( i = 0; i < nvars; i++ )
      {
         idx = varinds[i];
         var = SCIPmatrixGetVar(matrix, idx);

         assert(SCIPisLE(scip, lbs[idx], newlbs[i]));
         assert(SCIPisGE(scip, ubs[idx], newubs[i]));

#ifdef SCIP_DEBUG_SUBSCIP
         assert(SCIPisEQ(scip, newlbs[i], subsciplbs[i]));
         assert(SCIPisEQ(scip, newubs[i], subscipubs[i]));
#endif
         if( SCIPisGT(scip, newlbs[i], lbs[idx]) || SCIPisLT(scip, newubs[i], ubs[idx]) )
         {
            (*success) = TRUE;

            if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
            {
               lbs[idx] = SCIPceil(scip, newlbs[i]);
               ubs[idx] = SCIPfloor(scip, newubs[i]);
            }
            else
            {
               lbs[idx] = newlbs[i];
               ubs[idx] = newubs[i];
            }
         }
      }

#ifdef SCIP_DEBUG_SUBSCIP
      SCIPfreeBufferArray(scip, &subscipubs);
      SCIPfreeBufferArray(scip, &subsciplbs);
#endif
      SCIPfreeBufferArray(scip, &newubs);
      SCIPfreeBufferArray(scip, &newlbs);
   }

   SCIPfreeBufferArray(scip, &currentmaxinds);
   SCIPfreeBufferArray(scip, &cliquepartition);
   SCIPfreeBufferArray(scip, &binvarpos);
   SCIPfreeBufferArray(scip, &binvars);
   SCIPfreeBufferArray(scip, &cliquemaxinds);
   SCIPfreeBufferArray(scip, &breakpoints);
   SCIPfreeBufferArray(scip, &varinds);
   SCIPfreeBufferArray(scip, &signs);
   SCIPfreeBufferArray(scip, &row2coefs);
   SCIPfreeBufferArray(scip, &row1coefs);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRESOLCOPY(presolCopyTworowcomb)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of tworowcomb presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolCopyTworowcomb NULL
#endif


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeTworowcomb)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** initialization method of presolver (called after problem was transformed) */
static
SCIP_DECL_PRESOLINIT(presolInitTworowcomb)
{
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   presoldata->nchgbnds = 0;
   presoldata->nuselessruns = 0;

   return SCIP_OKAY;
}

/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitTworowcomb)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of tworowcomb presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitTworowcomb NULL
#endif


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreTworowcomb)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of tworowcomb presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreTworowbound NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#if 0
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreTworowcomb)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of tworowcomb presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitpreTworowcomb NULL
#endif


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecTworowcomb)
{  /*lint --e{715}*/

   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_Bool infeasible;
   SCIP_PRESOLDATA* presoldata;

   int i;
   int j;
   int k;

   SCIPinfoMessage(scip, NULL, "starting tworowcomb\n");

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;
   infeasible = FALSE;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   if( presoldata->nuselessruns >= 5 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, TRUE, &initialized, &complete) );

   /* we only work on pure MIPs currently */
   if( initialized && complete )
   {
      int oldnchgbds;
      int oldnfixedvars;
      int ndelcons;
      SCIP_Longint nrows;
      SCIP_Longint ncols;
      SCIP_Real* oldlbs;
      SCIP_Real* oldubs;
      SCIP_Real* newlbs;
      SCIP_Real* newubs;
      int* rowidxptr;
      SCIP_Real* rowvalptr;
      SCIP_VAR* var;

      SCIP_Longint maxhashes;

      int maxlen;
      int pospp;
      int listsizepp;
      int posmm;
      int listsizemm;
      int pospm;
      int listsizepm;
      int posmp;
      int listsizemp;

      int* hashlistpp;
      int* hashlistmm;
      int* hashlistpm;
      int* hashlistmp;

      int* rowidxlistpp;
      int* rowidxlistmm;
      int* rowidxlistpm;
      int* rowidxlistmp;

      SCIP_Bool finiterhs;
      int block1start;
      int block1end;
      int block2start;
      int block2end;

      SCIP_HASHSET* pairhashset;

#ifdef SCIP_HASHBLOCK_INFO
      int ppblocks;
      int mmblocks;
      int pmblocks;
      int mpblocks;

      int checkedppblocks;
      int checkedmmblocks;
      int checkedpmblocks;
      int checkedmpblocks;

      int pppairs;
      int pmpairs;

      int checkedpppairs;
      int checkedpmpairs;

      ppblocks = 0;
      mmblocks = 0;
      pmblocks = 0;
      mpblocks = 0;
      checkedppblocks = 0;
      checkedmmblocks = 0;
      checkedpmblocks = 0;
      checkedmpblocks = 0;
      pppairs = 0;
      pmpairs = 0;
      checkedpppairs = 0;
      checkedpmpairs = 0;
#endif

      nrows = SCIPmatrixGetNRows(matrix);
      ncols = SCIPmatrixGetNColumns(matrix);

      if( nrows == 1 )
      {
         SCIPmatrixFree(scip, &matrix);
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistpp, nrows) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistmm, nrows) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistpm, nrows) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &hashlistmp, nrows) );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowidxlistpp, nrows) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowidxlistmm, nrows) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowidxlistpm, nrows) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &rowidxlistmp, nrows) );

      //TODO What about collisions?
      pospp = 0;
      posmm = 0;
      pospm = 0;
      posmp = 0;
      listsizepp = nrows;
      listsizemm = nrows;
      listsizepm = nrows;
      listsizemp = nrows;
      maxhashes = nrows * presoldata->maxhashfac;
      for( i = 0; i < nrows; i++)
      {
         if( pospp + posmm + pospm + posmp > maxhashes )
            break;

         rowvalptr = SCIPmatrixGetRowValPtr(matrix, i);
         rowidxptr = SCIPmatrixGetRowIdxPtr(matrix, i);
         finiterhs = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, i));
         maxlen = MIN(presoldata->maxconsiderednonzeros, SCIPmatrixGetRowNNonzs(matrix, i));
         for( j = 0; j < maxlen; j++)
         {
            for( k = j+1; k < maxlen; k++)
            {
               if( SCIPisPositive(scip, rowvalptr[j]) )
               {
                  if(SCIPisPositive(scip, rowvalptr[k]) )
                  {
                     addEntry(scip, &pospp, &listsizepp, &hashlistpp, &rowidxlistpp,
                        hashIndexPair(rowidxptr[j],rowidxptr[k]), i);
                     if( finiterhs )
                        addEntry(scip, &posmm, &listsizemm, &hashlistmm, &rowidxlistmm,
                           hashIndexPair(rowidxptr[j],rowidxptr[k]), i);
                  }
                  else
                  {
                     addEntry(scip, &pospm, &listsizepm, &hashlistpm, &rowidxlistpm,
                        hashIndexPair(rowidxptr[j],rowidxptr[k]), i);
                     if( finiterhs )
                        addEntry(scip, &posmp, &listsizemp, &hashlistmp, &rowidxlistmp,
                           hashIndexPair(rowidxptr[j],rowidxptr[k]), i);
                  }
               }
               else
               {
                  if(SCIPisPositive(scip, rowvalptr[k]) )
                  {
                     addEntry(scip, &posmp, &listsizemp, &hashlistmp, &rowidxlistmp,
                        hashIndexPair(rowidxptr[j],rowidxptr[k]), i);
                     if( finiterhs )
                        addEntry(scip, &pospm, &listsizepm, &hashlistpm, &rowidxlistpm,
                           hashIndexPair(rowidxptr[j],rowidxptr[k]), i);
                  }
                  else
                  {
                     addEntry(scip, &posmm, &listsizemm, &hashlistmm, &rowidxlistmm,
                        hashIndexPair(rowidxptr[j],rowidxptr[k]), i);
                     if( finiterhs )
                        addEntry(scip, &pospp, &listsizepp, &hashlistpp, &rowidxlistpp,
                           hashIndexPair(rowidxptr[j],rowidxptr[k]), i);
                  }
               }
            }
         }
      }

#ifdef SCIP_DEBUG_HASHING
      SCIPdebugMsg(scip, "pp\n");
      for( i = 0; i < pospp; i++)
        SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistpp[i], rowidxlistpp[i]);
      SCIPdebugMsg(scip, "mm\n");
      for( i = 0; i < posmm; i++)
        SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistmm[i], rowidxlistmm[i]);
      SCIPdebugMsg(scip, "pm\n");
      for( i = 0; i < pospm; i++)
        SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistpm[i], rowidxlistpm[i]);
      SCIPdebugMsg(scip, "mp\n");
      for( i = 0; i < posmp; i++)
        SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistmp[i], rowidxlistmp[i]);
#endif
      SCIPdebugMsg(scip, "hashlist sizes: pp %d, mm %d, pm %d, mp %d \n", pospp, posmm, pospm, posmp);

      SCIPsortIntInt(hashlistpp, rowidxlistpp, pospp);
      SCIPsortIntInt(hashlistmm, rowidxlistmm, posmm);
      SCIPsortIntInt(hashlistpm, rowidxlistpm, pospm);
      SCIPsortIntInt(hashlistmp, rowidxlistmp, posmp);

#ifdef SCIP_DEBUG_HASHING
      SCIPdebugMsg(scip, "sorted pp\n");
      for( i = 0; i < pospp; i++)
        SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistpp[i], rowidxlistpp[i]);
      SCIPdebugMsg(scip, "sorted mm\n");
      for( i = 0; i < posmm; i++)
        SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistmm[i], rowidxlistmm[i]);
      SCIPdebugMsg(scip, "sorted pm\n");
      for( i = 0; i < pospm; i++)
        SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistpm[i], rowidxlistpm[i]);
      SCIPdebugMsg(scip, "sorted mp\n");
      for( i = 0; i < posmp; i++)
        SCIPdebugMsg(scip, "%d: hash  = %d, rowidx = %d\n", i, hashlistmp[i], rowidxlistmp[i]);
#endif

      SCIP_CALL( SCIPallocBufferArray(scip, &oldlbs, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &oldubs, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newlbs, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &newubs, ncols) );

      for( i = 0; i < SCIPmatrixGetNColumns(matrix); i++ )
      {
         var = SCIPmatrixGetVar(matrix, i);
         oldlbs[i] = SCIPvarGetLbLocal(var);
         oldubs[i] = SCIPvarGetUbLocal(var);
         newlbs[i] = oldlbs[i];
         newubs[i] = oldubs[i];
      }

      SCIP_CALL( SCIPhashsetCreate(&pairhashset, SCIPblkmem(scip), 1) );

      /* Process pp and mm lists */
      if( pospp > 0 && posmm > 0 )
      {
         SCIP_Longint ncombines;
         SCIP_Longint maxcombines;
         SCIP_Bool finished;
         SCIP_Bool success;
         int combinefails;
         int retrievefails;
         ROWPAIR rowpair;

         finished = FALSE;
         i = 0;
         j = 0;
         block1start = 0;
         block1end = 0;
         block2start = 0;
         block2end = 0;
         maxcombines = nrows * presoldata->maxpairfac;
         ncombines = 0;
         combinefails = 0;
         retrievefails = 0;
         findNextBlock(hashlistpp, pospp, &block1start, &block1end);
         findNextBlock(hashlistmm, posmm, &block2start, &block2end);
         SCIPdebugMsg(scip, "processing pp and mm\n");
#ifdef SCIP_HASHBLOCK_INFO
         checkedppblocks++;
         checkedmmblocks++;
#endif
         while( !finished )
         {
            //SCIPdebugMsg(scip, "block1: %d -> %d, block2: %d -> %d, hashes: %d, %d\n", block1start, block1end, block2start, block2end, hashlistpp[block1start], hashlistmm[block2start]);
            if( hashlistpp[block1start] == hashlistmm[block2start] )
            {
#ifdef SCIP_HASHBLOCK_INFO
               checkedpppairs++;
#endif
               for( i = block1start; i < block1end; i++ )
               {
                  for( j = block2start; j < block2end; j++ )
                  {
                     if( rowidxlistpp[i] != rowidxlistmm[j] )
                     {
                        rowpair.row1idx = MIN(rowidxlistpp[i], rowidxlistmm[j]);
                        rowpair.row2idx = MAX(rowidxlistpp[i], rowidxlistmm[j]);
                        if( !SCIPhashsetExists(pairhashset, encodeRowPair(&rowpair)) )
                        {
                           assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row1idx)));
                           assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row2idx)));

                           success = FALSE;

                           combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                       FALSE, FALSE, newlbs, newubs, &success);
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx)) )
                              combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                          TRUE, FALSE, newlbs, newubs, &success);
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx)) )
                              combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                          FALSE, TRUE, newlbs, newubs, &success);
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx))
                                && !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx)) )
                              combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                          TRUE, TRUE, newlbs, newubs, &success);
                           if( success )
                              combinefails = 0;
                           else
                              combinefails++;

                           SCIP_CALL( SCIPhashsetInsert(pairhashset, SCIPblkmem(scip), encodeRowPair(&rowpair)) );
                           ncombines++;
                           if( ncombines >= maxcombines || combinefails >= presoldata->maxcombinefails )
                              finished = TRUE;

                           retrievefails = 0;
                        }
                        else if( retrievefails < presoldata->maxretrievefails )
                           retrievefails++;
                        else
                           finished = TRUE;
                     }
                     if( finished )
                        break;
                  }
                  if( finished )
                     break;
               }

               if( block1end < pospp && block2end < posmm )
               {
                  findNextBlock(hashlistpp, pospp, &block1start, &block1end);
                  findNextBlock(hashlistmm, posmm, &block2start, &block2end);
#ifdef SCIP_HASHBLOCK_INFO
                  checkedppblocks++;
                  checkedmmblocks++;
#endif

               }
               else
                  finished = TRUE;
            }
            else if( hashlistpp[block1start] < hashlistmm[block2start] && block1end < pospp )
            {
               findNextBlock(hashlistpp, pospp, &block1start, &block1end);
#ifdef SCIP_HASHBLOCK_INFO
               checkedppblocks++;
#endif

            }
            else if( hashlistpp[block1start] > hashlistmm[block2start] && block2end < posmm )
            {
               findNextBlock(hashlistmm, posmm, &block2start, &block2end);
#ifdef SCIP_HASHBLOCK_INFO
               checkedmmblocks++;
#endif
            }
            else
               finished = TRUE;
         }
      }

      /* Process pm and mp lists */
      if( pospm > 0 && posmp > 0 )
      {
         SCIP_Longint ncombines;
         SCIP_Longint maxcombines;
         SCIP_Bool finished;
         SCIP_Bool success;
         int combinefails;
         int retrievefails;
         ROWPAIR rowpair;

         finished = FALSE;
         i = 0;
         j = 0;
         block1start = 0;
         block1end = 0;
         block2start = 0;
         block2end = 0;
         maxcombines = nrows * presoldata->maxpairfac;
         ncombines = 0;
         combinefails = 0;
         retrievefails = 0;
         findNextBlock(hashlistpm, pospm, &block1start, &block1end);
         findNextBlock(hashlistmp, posmp, &block2start, &block2end);
         SCIPdebugMsg(scip, "processing pm and mp\n");
#ifdef SCIP_HASHBLOCK_INFO
         checkedpmblocks++;
         checkedmpblocks++;
#endif
         while( !finished )
         {
            if( hashlistpm[block1start] == hashlistmp[block2start] )
            {
#ifdef SCIP_HASHBLOCK_INFO
               checkedpmpairs++;
#endif
               for( i = block1start; i < block1end; i++ )
               {
                  for( j = block2start; j < block2end; j++ )
                  {
                     if( rowidxlistpm[i] != rowidxlistmp[j] )
                     {
                        rowpair.row1idx = MIN(rowidxlistpm[i], rowidxlistmp[j]);
                        rowpair.row2idx = MAX(rowidxlistpm[i], rowidxlistmp[j]);
                        if( ! SCIPhashsetExists(pairhashset, encodeRowPair(&rowpair)) )
                        {
                           assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row1idx)));
                           assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row2idx)));

                           success = FALSE;

                           combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                       FALSE, FALSE, newlbs, newubs, &success);
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx)) )
                              combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                          TRUE, FALSE, newlbs, newubs, &success);
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx)) )
                              combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                          FALSE, TRUE, newlbs, newubs, &success);
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx))
                                && !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx)) )
                              combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                          TRUE, TRUE, newlbs, newubs, &success);
                           if( success )
                              combinefails = 0;
                           else
                              combinefails++;

                           SCIP_CALL( SCIPhashsetInsert(pairhashset, SCIPblkmem(scip), encodeRowPair(&rowpair)) );
                           ncombines++;
                           if( ncombines >= maxcombines || combinefails >= presoldata->maxcombinefails )
                              finished = TRUE;

                           retrievefails = 0;
                        }
                        else if( retrievefails < presoldata->maxretrievefails )
                           retrievefails++;
                        else
                           finished = TRUE;
                     }
                     if( finished )
                        break;
                  }
                  if( finished )
                     break;
               }

               if( block1end < pospm && block2end < posmp )
               {
                  findNextBlock(hashlistpm, pospm, &block1start, &block1end);
                  findNextBlock(hashlistmp, posmp, &block2start, &block2end);
#ifdef SCIP_HASHBLOCK_INFO
                  checkedpmblocks++;
                  checkedmpblocks++;
#endif
               }
               else
                  finished = TRUE;
            }
            else if( hashlistpm[block1start] < hashlistmp[block2start] && block1end < pospm )
            {
               findNextBlock(hashlistpm, pospm, &block1start, &block1end);
#ifdef SCIP_HASHBLOCK_INFO
               checkedpmblocks++;
#endif
            }
            else if( hashlistpm[block1start] > hashlistmp[block2start] && block2end < posmp )
            {
               findNextBlock(hashlistmp, posmp, &block2start, &block2end);
#ifdef SCIP_HASHBLOCK_INFO
               checkedmpblocks++;
#endif
            }
            else
               finished = TRUE;
         }
      }

      /* Apply reductions */
      oldnchgbds = *nchgbds;
      oldnfixedvars = *nfixedvars;
      for( i = 0; i < SCIPmatrixGetNColumns(matrix); i++ )
      {
         SCIP_Bool bndwastightened;
         SCIP_Bool fixed;

         var = SCIPmatrixGetVar(matrix, i);

         assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
            || (SCIPisEQ(scip, newlbs[i], SCIPceil(scip, newlbs[i])) && (SCIPisEQ(scip, newubs[i], SCIPfloor(scip, newubs[i])))));

         if( SCIPisEQ(scip, newlbs[i], newubs[i]) )
         {
            SCIP_CALL( SCIPfixVar(scip, var, newlbs[i], &infeasible, &fixed) );

            if( infeasible )
            {
               SCIPdebugMessage(" -> infeasible fixing of variable %s\n", SCIPvarGetName(var));
               break;
            }

            if( fixed )
            {
               SCIPdebugMessage("variable %s fixed to %g\n", SCIPvarGetName(var), newlbs[i]);
               (*nfixedvars)++;
            }
         }

         if( SCIPisLT(scip, oldlbs[i], newlbs[i]) )
         {
            SCIP_CALL( SCIPtightenVarLb(scip, var, newlbs[i], FALSE, &infeasible, &bndwastightened) );

            if( infeasible )
            {
               SCIPdebugMessage(" -> infeasible lower bound tightening of variable %s\n", SCIPvarGetName(var));
               break;
            }

            if( bndwastightened )
            {
               SCIPdebugMessage("lower bound of %s changed from %g to %g\n", SCIPvarGetName(var), oldlbs[i], newlbs[i]);
               (*nchgbds)++;
            }
         }

         if( SCIPisGT(scip, oldubs[i], newubs[i]) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var, newubs[i], FALSE, &infeasible, &bndwastightened) );

            if( infeasible )
            {
               SCIPdebugMessage(" -> infeasible upper bound tightening\n");
               break;
            }

            if( bndwastightened )
            {
               SCIPdebugMessage("upper bound of %s changed from %g to %g\n", SCIPvarGetName(var), oldubs[i], newubs[i]);
               (*nchgbds)++;
            }
         }
      }

      /* check for redundant constraints */
      ndelcons = 0;
      for( i = 0; i < nrows; i++ )
      {
         SCIP_Real rowinf;
         SCIP_Real rowsup;

         rowidxptr = SCIPmatrixGetRowIdxPtr(matrix, i);
         rowvalptr = SCIPmatrixGetRowValPtr(matrix, i);

         rowinf = 0;
         for( j = 0; j < SCIPmatrixGetRowNNonzs(matrix, i); j++ )
         {
            k = rowidxptr[j];
            if( SCIPisPositive(scip, rowvalptr[j]) )
            {
               if( !SCIPisInfinity(scip, -newlbs[k]) )
                  rowinf += rowvalptr[j] * newlbs[k];
               else
               {
                  rowinf = -SCIPinfinity(scip);
                  break;
               }
            }
            else if( SCIPisNegative(scip, rowvalptr[j]) )
            {
               if( !SCIPisInfinity(scip, newubs[k]) )
                  rowinf += rowvalptr[j] * newubs[k];
               else
               {
                  rowinf = -SCIPinfinity(scip);
                  break;
               }
            }
         }

         rowsup = 0;
         for( j = 0; j < SCIPmatrixGetRowNNonzs(matrix, i); j++ )
         {
            k = rowidxptr[j];
            if( SCIPisPositive(scip, rowvalptr[j]) )
            {
               if( !SCIPisInfinity(scip, newubs[k]) )
                  rowsup += rowvalptr[j] * newubs[k];
               else
               {
                  rowsup = SCIPinfinity(scip);
                  break;
               }
            }
            else if( SCIPisNegative(scip, rowvalptr[j]) )
            {
               if( !SCIPisInfinity(scip, -newlbs[k]) )
                  rowsup += rowvalptr[j] * newlbs[k];
               else
               {
                  rowsup = SCIPinfinity(scip);
                  break;
               }
            }
         }

         if( !SCIPisInfinity(scip, -rowinf) )
         {
            if( SCIPisLE(scip, SCIPmatrixGetRowLhs(matrix, i), rowinf) && SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, i)) )
            {
               ndelcons++;
               SCIP_CALL( SCIPdelCons(scip, SCIPmatrixGetCons(matrix, i)) );
               SCIPdebugMsg(scip, "removing redundant cxonstraint %s\n", SCIPmatrixGetRowName(matrix, i));
            }
            if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, i)) && SCIPisGT(scip, rowinf, SCIPmatrixGetRowRhs(matrix, i)) )
            {
               SCIPdebugMsg(scip, "infeasibility detected in %s, rowinf = %g\n", SCIPmatrixGetRowName(matrix, i), rowinf);
               infeasible = TRUE;
            }
         }
         if( !SCIPisInfinity(scip, rowsup) )
         {
            if( SCIPisGT(scip, SCIPmatrixGetRowLhs(matrix, i), rowsup) )
            {
               SCIPdebugMsg(scip, "infeasibility detected in %s, rowsup = %g\n", SCIPmatrixGetRowName(matrix, i), rowsup);
               infeasible = TRUE;
            }
            // TODO we currently can't do this for ranged rows or equality constraints
            /* if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, i)) && SCIPisGE(scip, SCIPmatrixGetRowRhs(matrix, i), rowsup) ) */
            /* { */
            /*    ndelcons++; */
            /*    SCIP_CALL( SCIPdelCons(scip, SCIPmatrixGetCons(matrix, i)) ); */
            /*    SCIPdebugMsg(scip, "removing redundant cxonstraint %s\n", SCIPmatrixGetRowName(matrix, i)); */
            /* } */
         }
      }
      (*ndelconss) += ndelcons;

#ifdef SCIP_HASHBLOCK_INFO
      block1start = 0;
      block1end = 0;
      block2start = 0;
      block2end = 0;
      if( pospp > 0 )
      {
         findNextBlock(hashlistpp, pospp, &block1start, &block1end);
         ppblocks++;
      }
      if( posmm > 0 )
      {
         findNextBlock(hashlistmm, posmm, &block2start, &block2end);
         mmblocks++;
      }

      while( block1end < pospp || block2end < posmm )
      {
         if( pospp > block1start && posmm > block2start && hashlistpp[block1start] == hashlistmm[block2start] )
            pppairs++;

         if( block1end < pospp && (posmm == 0 || block2end >= posmm || hashlistpp[block1start] <= hashlistmm[block2start]) )
         {
            findNextBlock(hashlistpp, pospp, &block1start, &block1end);
            ppblocks++;
         }

         if( block2end < posmm && (pospp == 0 || block1end >= pospp || hashlistpp[block1start] > hashlistmm[block2start]) )
         {
            findNextBlock(hashlistmm, posmm, &block2start, &block2end);
            mmblocks++;
         }
      }
      block1start = 0;
      block1end = 0;
      block2start = 0;
      block2end = 0;
      if( pospm > 0 )
      {
         findNextBlock(hashlistpm, pospm, &block1start, &block1end);
         pmblocks++;
      }
      if( posmp > 0 )
      {
         findNextBlock(hashlistmp, posmp, &block2start, &block2end);
         mpblocks++;
      }

      while( block1end < pospm || block2end < posmp )
      {
         if( pospm > block1start && posmp > block2start && hashlistpm[block1start] == hashlistmp[block2start] )
            pmpairs++;

         if( block1end < pospm && (posmp == 0 || block2end >= posmp || hashlistpm[block1start] <= hashlistmp[block2start]) )
         {
            findNextBlock(hashlistpm, pospm, &block1start, &block1end);
            pmblocks++;
         }

         if( block2end < posmp && (pospm == 0 || block1end >= pospm || hashlistpm[block1start] > hashlistmp[block2start]) )
         {
            findNextBlock(hashlistmp, posmp, &block2start, &block2end);
            mpblocks++;
         }
      }
      SCIPinfoMessage(scip, NULL, "2RC: %d/%d pp-blocks, %d/%d mm-blocks, %d/%d pm-blocks, %d/%d mp-blocks, %d/%d pp-pairs, %d/%d pm-pairs\n", checkedppblocks, ppblocks, checkedmmblocks, mmblocks, checkedpmblocks, pmblocks, checkedmpblocks, mpblocks, checkedpppairs, pppairs, checkedpmpairs, pmpairs);
#endif


      /* set result */
      if( *nchgbds > oldnchgbds || *nfixedvars > oldnfixedvars || ndelcons > 0 )
      {
         *result = SCIP_SUCCESS;
         presoldata->nuselessruns = 0;
         SCIPinfoMessage(scip, NULL, "tworowcomb evaluated %d pairs to tighten %d bounds, fix %d variables and delete %d redundant constraints\n", SCIPhashsetGetNElements(pairhashset), *nchgbds - oldnchgbds, *nfixedvars - oldnfixedvars, ndelcons);
      }
      else if( infeasible )
      {
         *result = SCIP_CUTOFF;
         SCIPinfoMessage(scip, NULL, "tworowcomb detected infeasibility\n");
      }
      else
      {
         presoldata->nuselessruns++;
         SCIPinfoMessage(scip, NULL, "tworowcomb evaluated %d pairs without success\n", SCIPhashsetGetNElements(pairhashset));
      }

      SCIPhashsetFree(&pairhashset, SCIPblkmem(scip));
      SCIPfreeBufferArray(scip, &newubs);
      SCIPfreeBufferArray(scip, &newlbs);
      SCIPfreeBufferArray(scip, &oldubs);
      SCIPfreeBufferArray(scip, &oldlbs);
      SCIPfreeBlockMemoryArray(scip, &rowidxlistmp, listsizemp);
      SCIPfreeBlockMemoryArray(scip, &rowidxlistpm, listsizepm);
      SCIPfreeBlockMemoryArray(scip, &rowidxlistmm, listsizemm);
      SCIPfreeBlockMemoryArray(scip, &rowidxlistpp, listsizepp);
      SCIPfreeBlockMemoryArray(scip, &hashlistmp, listsizemp);
      SCIPfreeBlockMemoryArray(scip, &hashlistpm, listsizepm);
      SCIPfreeBlockMemoryArray(scip, &hashlistmm, listsizemm);
      SCIPfreeBlockMemoryArray(scip, &hashlistpp, listsizepp);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the tworowcomb presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolTworowcomb(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create tworowcomb presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );
   /* TODO: (optional) create presolver specific data here */

   presol = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecTworowcomb,
         presoldata) );

   assert(presol != NULL);

   //SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyTworowcomb) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeTworowcomb) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitTworowcomb) );
   //SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitTworowcomb) );
   //SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreTworowcomb) );
   //SCIP_CALL( SCIPsetPresolExitpre(scip, presol, presolExitpreTworowcomb) );

   /* add tworowcomb presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowcomb/maxconsiderednonzeros",
         "maximal number of considered non-zeros within one row (-1: no limit)",
         &presoldata->maxconsiderednonzeros, FALSE, DEFAULT_MAXCONSIDEREDNONZEROS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowcomb/maxretrievefails",
         "maximal number of consecutive useless hashtable retrieves",
         &presoldata->maxretrievefails, FALSE, DEFAULT_MAXRETRIEVEFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowcomb/maxcombinefails",
         "maximal number of consecutive useless row combines",
         &presoldata->maxcombinefails, FALSE, DEFAULT_MAXCOMBINEFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/tworowcomb/maxhashfac",
         "Maximum number of hashlist entries as multiple of number of rows in the problem (-1: no limit)",
         &presoldata->maxhashfac, FALSE, DEFAULT_MAXHASHFAC, -1, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/tworowcomb/maxpairfac",
         "Maximum number of processed row pairs as multiple of the number of rows in the problem (-1: no limit)",
         &presoldata->maxpairfac, FALSE, DEFAULT_MAXPAIRFAC, -1, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
