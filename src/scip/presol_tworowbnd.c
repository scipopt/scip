/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_tworowbnd.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief  do bound tightening by using two rows
 * @author Dieter Weninger
 * @author Patrick Gemander
 *
 * TODO update documentation
 * Perform bound tightening on two inequalities with some common variables.
 *
 * Let two constraints be given:
 * \f{eqnarray*}{
 *   A_{iR} x_R + A_{iS} x_S              \geq b_i\\
 *   A_{kR} x_R              + A_{kT} x_T \geq b_k
 * \f}
 * with \f$N\f$ the set of variable indexes, \f$R \subseteq N\f$, \f$S \subseteq N\f$, \f$T \subseteq N\f$,
 * \f$R \cap S = \emptyset\f$, \f$R \cap T = \emptyset\f$, \f$S \cap T = \emptyset\f$ and \f$i \not= k\f$.
 *
 * Solve the following two LPs
 * \f{eqnarray*}{
 *   L = \min \{ A_{kR} x_R : A_{iR} x_R + A_{iS} x_S \geq b_i \}\\
 *   U = \max \{ A_{kR} x_R : A_{iR} x_R + A_{iS} x_S \geq b_i \}
 * \f}
 * and use \f$L\f$ and \f$U\f$ for getting bounds on \f$x_T\f$.
 *
 * If \f$L + \mbox{infimum}(A_{kT}x_T) \geq b_k\f$, then the second constraint above is redundant.
 */

//TODO remove these
//#define SCIP_DEBUG
//#define SCIP_DEBUG_SINGLEROWLP
//#define SCIP_DEBUG_2RB
//#define SCIP_DEBUG_BOUNDS
//#define SCIP_DEBUG_HASHING

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_matrix.h"
#include "presol_tworowbnd.h"

// TODO Fix aligning
#define PRESOL_NAME            "tworowbnd"
#define PRESOL_DESC            "do bound tigthening by using two rows"
#define PRESOL_PRIORITY         -2000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS        0 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_MAXCONSIDEREDNONZEROS  100
#define DEFAULT_MAXRETRIEVEFAILS       1000
#define DEFAULT_MAXCOMBINEFAILS        1000
#define DEFAULT_MAXHASHFAC             10
#define DEFAULT_MAXPAIRFAC             1

/*
 * Data structures
 */

/** presolver data */
struct SCIP_PresolData
{
   int maxpairfac;
   int maxhashfac;
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

static
SCIP_RETCODE addEntry
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

static
void findNextBlock
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
}

/** Solve single-row LP of the form
 *  min c^T x, s.t. a^T x >= b
 *  TODO proper documentation
 */
static
SCIP_RETCODE solveSingleRowLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            a,                  /**< constraint coefficients */
   SCIP_Real             b,                  /**< right hand side */
   SCIP_Real*            c,                  /**< objective coefficients */
   SCIP_Real*            lbs,                /**< lower variable bounds */
   SCIP_Real*            ubs,                /**< upper variable bounds */
   int                   len,                /**< length of arrays */
   SCIP_Real*            obj,                /**< objective value of solution */
   SCIP_Bool*            solvable            /**< status whether LP was solvable */
   )
{
   int i;
   int k;
   int nvars;
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real mincost;
   SCIP_Real maxgain;

#ifdef SCIP_DEBUG_SINGLEROWLP
   SCIPdebugMsg(scip, "solving single row LP with %d variables\n", len);
#endif

   nvars = 0;
   (*obj) = 0;
   (*solvable) = TRUE;
   mincost = SCIPinfinity(scip);
   maxgain = 0;
   for( i = 0; i < len; i++)
   {
      /* Handle variables with zero weight */
      if( SCIPisZero(scip, a[i]) )
      {
         /* a[i] = 0, c[i] > 0 */
         if( SCIPisPositive(scip, c[i]) )
         {
            if( SCIPisInfinity(scip, -lbs[i]) )
            {
               (*solvable) = FALSE;
               return SCIP_OKAY;
            }
            else
               (*obj) += c[i] * lbs[i];
         }
         /* a[i] = 0, c[i] < 0 */
         else if( SCIPisNegative(scip, c[i]) )
         {
            if( SCIPisInfinity(scip, ubs[i]) )
            {
               (*solvable) = FALSE;
               return SCIP_OKAY;
            }
            else
               (*obj) += c[i] * ubs[i];
         }
         /* Note that variables with a[i] = 0, c[i] = 0 can be ignored */
         continue;
      }

      /* Handle free variables */
      if( SCIPisInfinity(scip, -lbs[i]) && SCIPisInfinity(scip, ubs[i]) )
      {
         /* The problem is unbounded */
         if( (SCIPisPositive(scip, c[i]) && SCIPisNegative(scip, a[i])) ||
             (SCIPisNegative(scip, c[i]) && SCIPisPositive(scip, a[i])) )
         {
            (*solvable) = FALSE;
            return SCIP_OKAY;
         }
         else
         {
            mincost = MIN(mincost, c[i]/a[i]);
            maxgain = MAX(maxgain, c[i]/a[i]);
         }
         continue;
      }

      /* Swap variable orientation if lower bound is infinite */
      if( SCIPisInfinity(scip, -lbs[i]) )
      {
         c[i] = -c[i];
         a[i] = -a[i];
         lb = -ubs[i];
         ub = -lbs[i];
      }
      else
      {
         lb = lbs[i];
         ub = ubs[i];
      }

      /* Handle variables with infinite upper bound */
      if( SCIPisInfinity(scip, ub) )
      {
         if( SCIPisPositive(scip, a[i]) )
         {
            /* a[i] > 0, c[i] >= 0 */
            if( !SCIPisNegative(scip, c[i]) )
            {
               mincost = MIN(mincost, c[i]/a[i]);
            }
            /* a[i] > 0, c[i] < 0 */
            else
            {
               (*solvable) = FALSE;
               return SCIP_OKAY;
            }
         }
         /* a[i] < 0, c[i] < 0 */
         else if( SCIPisNegative(scip, c[i]) )
         {
            maxgain = MAX(maxgain, c[i]/a[i]);
         }
         /* a[i] < 0, c[i] >= 0 results in dual fixing of this variable, which is included in the bound shift below */

         /* Shift lower bound to zero */
         if( !SCIPisZero(scip, lb) )
         {
            (*obj) += c[i] * lb;
            b -= a[i] * lb;
         }
         continue;
      }

      /* Handle fixed variables */
      if( SCIPisEQ(scip, lb, ub) )
      {
         (*obj) += c[i] * lb;
         b -= a[i] * lb;
         continue;
      }

      /* Dual fixing for variables with finite bounds */
      if( !SCIPisNegative(scip, c[i]) && SCIPisNegative(scip, a[i]) )
      {
         (*obj) += c[i] * lb;
         b -= a[i] * lb;
         continue;
      }
      else if( !SCIPisPositive(scip, c[i]) && SCIPisPositive(scip, a[i]) )
      {
         (*obj) += c[i] * ub;
         b -= a[i] * ub;
         continue;
      }

      assert(!SCIPisInfinity(scip, -lb));
      assert(!SCIPisInfinity(scip, ub));

      /** At this point the variable has finite bounds and a[i],c[i] are both positive or both negative.
       * Normalize variable such that
       *  1. x_i \in [0,1]
       *  2. a[i] > 0
       *  3. c[i] >= 0
       * and calculate its "unit price" c[i]/a[i]. */
      if( SCIPisNegative(scip, a[i]) )
      {
         c[i] = -c[i];
         a[i] = -a[i];
         lb = -ubs[i];
         ub = -lbs[i];
      }

      assert(SCIPisPositive(scip, a[i]) && SCIPisPositive(scip, c[i]));

      /* Adjust objective offset and b to shift lower bound to zero */
      (*obj) += c[i] * lb;
      b -= a[i] * lb;

      /* Calculate unit price */
      c[nvars] = c[i] / a[i];

      /* Normalize bound [0, ub] to [0,1] */
      a[nvars] = (ub - lb) * a[i];
      nvars++;
   }

#ifdef SCIP_DEBUG_SINGLEROWLP
   SCIPdebugMsg(scip, "After preprocessing: obj = %g, b = %g, nvars = %d, mincost = %g, maxgain = %g\n", (*obj), b, nvars, mincost, maxgain);
#endif

   /** Actual solving starts here.
    * If maxgain > 0 holds, we have a variable that can relax the constraint to an arbitrary degree while yielding
    * a certain profit per unit. This will be called downslack. If mincost < inf holds, we have a variable that can
    * always satisfy the constraint at a certain unit price. This will be called upslack. */

   /* Problem is unbounded since the downslack variable yields higher gains than the upslack variable costs */
   if( SCIPisLT(scip, mincost, maxgain) )
   {
      (*solvable) = FALSE;
      return SCIP_OKAY;
   }
   /* Solution is trivial as we have slack variables of equal price for both directions */
   else if( SCIPisEQ(scip, mincost, maxgain) )
   {
      /* Use all elements with cost smaller than maxgain */
      for( i = 0; i < nvars; i++ )
      {
         if( SCIPisLT(scip, c[i], maxgain) )
         {
            (*obj) += c[i] * a[i];
            b -= a[i];
         }
      }
      /* Use slack variable to satisfy constraint */
      (*obj) += mincost * b;
      return SCIP_OKAY;
   }
   /** mincost > maxgain
    *  In this case we need to solve the problem for the remaining variables with mincost > c[i] > maxgain.
    */
   else
   {
      /* Only keep variables that are cheaper than the upslack variable */
      if( !SCIPisInfinity(scip, mincost) )
      {
         k = 0;
         for( i = 0; i < nvars; i++ )
         {
            if( SCIPisLT(scip, c[i], mincost) )
            {
               c[k] = c[i];
               a[k] = a[i];
               k++;
            }
         }
         nvars = k;
      }

      /* Exploit all variables that are cheaper than the downslack variable */
      if( !SCIPisZero(scip, maxgain) )
      {
         k = 0;
         for( i = 0; i < nvars; i++ )
         {
            if( SCIPisLE(scip, c[i], maxgain) )
            {
               (*obj) += c[i] * a[i];
               b -= a[i];
            }
            else
            {
               c[k] = c[i];
               a[k] = a[i];
               k++;
            }
         }
         if( !SCIPisPositive(scip, b) )
         {
            (*obj) += maxgain * b;
            return SCIP_OKAY;
         }
         nvars = k;
      }

#ifdef SCIP_DEBUG_SINGLEROWLP
      SCIPdebugMsg(scip, "After exploiting slacks: obj = %g, nvars = %d\n", (*obj), nvars);
#endif

      /* If there are no variables left we can trivially put together a solution or determine infeasibility */
      if( nvars == 0 )
      {
         if( !SCIPisInfinity(scip, mincost) )
         {
            (*obj) += mincost * b;
            return SCIP_OKAY;
         }
         else
         {
            (*solvable) = FALSE;
            return SCIP_OKAY;
         }
      }
      /* Solve the remaining part of the problem */
      else
      {
         assert(nvars > 0);
         //TODO Check what happens if I can pack all elements
         //TODO test whether k == nvars holds iff all elements can be used
#ifdef SCIP_DEBUG_SINGLEROWLP
         for( i = 0; i < nvars; i++ )
            SCIPdebugMsg(scip, "c[%d] = %g, a[%d] = %g\n", i, c[i], i, a[i]);
#endif

         SCIPselectWeightedReal(c, a, b, nvars, &k);

#ifdef SCIP_DEBUG_SINGLEROWLP
         SCIPdebugMsg(scip, "k-mean = %g at index %d\n", c[k], k, b);
         for( i = 0; i < nvars; i++ )
            SCIPdebugMsg(scip, "c[%d] = %g, a[%d] = %g\n", i, c[i], i, a[i]);
#endif

         /* Finalize objective value of solution. First we use all elements cheaper than the k-median */
         for( i = 0; i < k; i++ )
         {
            (*obj) += c[i] * a[i];
            b -= a[i];
         }

#ifdef SCIP_DEBUG_SINGLEROWLP
         SCIPdebugMsg(scip, "LP is solved: b = %g\n", b);
#endif

         /* If the constraint is not yet satisfied, we have to fix that */
         if( SCIPisPositive(scip, b) )
         {
            /* There exists an element to satisfy the constraint */
            if( k < nvars )
            {
               (*obj) += c[k] * b;
               return SCIP_OKAY;
            }
            /* There is an upslack variable to satisfy the constraint */
            else if( !SCIPisInfinity(scip, mincost) )
            {
#ifdef SCIP_DEBUG_SINGLEROWLP
               SCIPdebugMsg(scip, "We use %g units of upslack to satisfy the constraint\n", b);
#endif
               (*obj) += mincost * b;
               return SCIP_OKAY;
            }
            /* We cannot satisfy the constraint so the problem is infeasible */
            else
            {
               (*solvable) = FALSE;
               return SCIP_OKAY;
            }
         }
         /* The constraint is already satisfied, i.e. b <= 0 */
         else
         {
            return SCIP_OKAY;
         }
      }
   }
}

static
SCIP_RETCODE transformAndSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   row1idx,
   int                   row2idx,
   SCIP_Bool             swaprow1,
   SCIP_Bool             swaprow2,
   SCIP_Real*            amin,
   SCIP_Real*            amax,
   SCIP_Real*            cmin,
   SCIP_Real*            cmax,
   SCIP_Bool*            cangetbnd,
   SCIP_Real*            lbs,
   SCIP_Real*            ubs,
   SCIP_Real*            newlbsmin,
   SCIP_Real*            newlbsmax,
   SCIP_Real*            newubsmin,
   SCIP_Real*            newubsmax,
   SCIP_Bool*            success,            /**< return (success || "found better bounds") */
   SCIP_Bool*            redundant,          /**< return whether first row is redundant */
   SCIP_Bool*            infeasible          /**< we return (infeasible || "detected infeasibility") */
   )
{
   int i;
   int j;
   int idx1;
   int idx2;
   int row1len;
   int row2len;
   int* row1idxptr;
   int* row2idxptr;
   SCIP_Real* row1valptr;
   SCIP_Real* row2valptr;
   int nvars;
   SCIP_Real minact;
   SCIP_Real maxact;
   int maxinfs;
   int mininfs;

   SCIP_Bool minsolvable;
   SCIP_Real minobj;
   SCIP_Bool maxsolvable;
   SCIP_Real maxobj;
   SCIP_Bool minswapsolvable;
   SCIP_Real minswapobj;
   SCIP_Bool maxswapsolvable;
   SCIP_Real maxswapobj;

   SCIP_Real newbnd;

   assert(!swaprow1 || (swaprow1 && !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, row1idx))));
   assert(!swaprow2 || (swaprow2 && !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, row2idx))));

   row1len = SCIPmatrixGetRowNNonzs(matrix, row1idx);
   row2len = SCIPmatrixGetRowNNonzs(matrix, row2idx);
   row1idxptr = SCIPmatrixGetRowIdxPtr(matrix, row1idx);
   row2idxptr = SCIPmatrixGetRowIdxPtr(matrix, row2idx);
   row1valptr = SCIPmatrixGetRowValPtr(matrix, row1idx);
   row2valptr = SCIPmatrixGetRowValPtr(matrix, row2idx);

   /* getting bounds for row1 */
   i = 0;
   j = 0;
   nvars = 0;
   mininfs = 0;
   maxinfs = 0;
   minact = 0;
   maxact = 0;
   while( i < row1len && j < row2len )
   {
      idx1 = row1idxptr[i];
      idx2 = row2idxptr[j];

      if( idx1 == idx2 )
      {
         cmin[nvars] = row1valptr[i];
         amin[nvars] = row2valptr[j];
         newlbsmin[nvars] = lbs[idx1];
         newubsmin[nvars] = ubs[idx1];
         cangetbnd[idx1] = FALSE;
         nvars++;
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "%g <= (%s) <= %g  has coefs %g and %g, %d LP vars\n", lbs[idx1], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx1)), ubs[idx1], row1valptr[i], row2valptr[j], nvars);
#endif
         i++;
         j++;
      }
      else if( idx1 < idx2 )
      {
         if( SCIPisPositive(scip, row1valptr[i]) )
         {
            if( SCIPisInfinity(scip, ubs[idx1]) )
               maxinfs++;
            else
               maxact -= row1valptr[i] * ubs[idx1];

            if( SCIPisInfinity(scip, -lbs[idx1]) )
               mininfs++;
            else
               minact -= row1valptr[i] * lbs[idx1];
         }
         else
         {
            if( SCIPisInfinity(scip, -lbs[idx1]) )
               maxinfs++;
            else
               maxact -= row1valptr[i] * lbs[idx1];

            if( SCIPisInfinity(scip, ubs[idx1]) )
               mininfs++;
            else
               minact -= row1valptr[i] * ubs[idx1];

            cangetbnd[idx1] = TRUE;
         }
         i++;
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "%g <= (%s) <= %g has coefs %g and 0.0, minact = %g, maxact = %g\n", lbs[idx1], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx1)), ubs[idx1], row1valptr[i], minact, maxact);
#endif
      }
      else
      {
         cmin[nvars] = 0.0;
         amin[nvars] = row2valptr[j];
         newlbsmin[nvars] = lbs[idx2];
         newubsmin[nvars] = ubs[idx2];
         cangetbnd[idx2] = FALSE;
         nvars++;
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "%g <= (%s) <= %g has coefs 0.0 and %g, %d LP vars\n", lbs[idx2], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx2)), ubs[idx2], row2valptr[j], nvars);
#endif
         j++;
      }
   }
   while( i < row1len )
   {
      idx1 = row1idxptr[i];
      if( SCIPisPositive(scip, row1valptr[i]) )
      {
         if( SCIPisInfinity(scip, ubs[idx1]) )
            maxinfs++;
         else
            maxact -= row1valptr[i] * ubs[idx1];

         if( SCIPisInfinity(scip, -lbs[idx1]) )
            mininfs++;
         else
            minact -= row1valptr[i] * lbs[idx1];
      }
      else
      {
         if( SCIPisInfinity(scip, -lbs[idx1]) )
            maxinfs++;
         else
            maxact -= row1valptr[i] * lbs[idx1];

         if( SCIPisInfinity(scip, ubs[idx1]) )
            mininfs++;
         else
            minact -= row1valptr[i] * ubs[idx1];
      }
      cangetbnd[idx1] = TRUE;
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "%g <= (%s) <= %g  has coefs %g and 0.0, minact = %g, maxact = %g\n", lbs[idx1], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx1)), ubs[idx1], row1valptr[i], minact, maxact);
#endif
      i++;
   }
   while( j < row2len )
   {
      idx2 = row2idxptr[j];
      cmin[nvars] = 0.0;
      amin[nvars] = row2valptr[j];
      newlbsmin[nvars] = lbs[idx2];
      newubsmin[nvars] = ubs[idx2];
      nvars++;
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "%g <= (%s) <= %g has coefs 0.0 and %g, %d LP vars\n", lbs[idx2], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx2)), ubs[idx2], row2valptr[j], nvars);
#endif
      j++;
   }

#ifdef SCIP_DEBUG_2RB
   SCIPdebugMsg(scip, "right hand sides: %g and %g\n", SCIPmatrixGetRowLhs(matrix, row1idx), SCIPmatrixGetRowLhs(matrix, row2idx));
#endif

   maxsolvable = FALSE;
   minsolvable = FALSE;
   maxswapsolvable = FALSE;
   minswapsolvable = FALSE;
   /* solve LPs */
   if( maxinfs <= 1 )
   {
      // TODO fix naming of amin/amax/... as *min is used to store the initial values and *max is used to make a copy for the solving process -> * and *copy would be better names
      for( i = 0; i < nvars; i++ )
      {
         amax[i] = amin[i];
         cmax[i] = -cmin[i];
         newlbsmax[i] = newlbsmin[i];
         newubsmax[i] = newubsmin[i];
      }
      SCIP_CALL( solveSingleRowLP(scip, amax, SCIPmatrixGetRowLhs(matrix, row2idx), cmax, newlbsmax, newubsmax, nvars, &maxobj, &maxsolvable) );
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "max-LP solved: obj = %g\n", maxobj);
#endif
   }

   if( mininfs == 0 || (mininfs == 1 && swaprow1) )
   {
      // copy stuff
      for( i = 0; i < nvars; i++ )
      {
         amax[i] = amin[i];
         cmax[i] = cmin[i];
         newlbsmax[i] = newlbsmin[i];
         newubsmax[i] = newubsmin[i];
      }
      SCIP_CALL( solveSingleRowLP(scip, amax, SCIPmatrixGetRowLhs(matrix, row2idx), cmax, newlbsmax, newubsmax, nvars, &minobj, &minsolvable) );
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "min-LP solved: obj = %g\n", minobj);
#endif
   }

   if( swaprow2 )
   {
      if( maxinfs <= 1 )
      {
         // copy stuff
         for( i = 0; i < nvars; i++ )
         {
            amax[i] = -amin[i];
            cmax[i] = -cmin[i];
            newlbsmax[i] = newlbsmin[i];
            newubsmax[i] = newubsmin[i];
         }
         SCIP_CALL( solveSingleRowLP(scip, amax, -SCIPmatrixGetRowRhs(matrix, row2idx), cmax, newlbsmax, newubsmax, nvars, &maxswapobj, &maxswapsolvable) );
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "maxswap-LP solved: obj = %g\n", maxswapobj);
#endif
      }

      if( mininfs == 0 || (mininfs == 1 && swaprow1) )
      {
         // copy stuff
         for( i = 0; i < nvars; i++ )
         {
            amax[i] = -amin[i];
            cmax[i] = cmin[i];
            newlbsmax[i] = newlbsmin[i];
            newubsmax[i] = newubsmin[i];
         }
         SCIP_CALL( solveSingleRowLP(scip, amax, -SCIPmatrixGetRowRhs(matrix, row2idx), cmax, newlbsmax, newubsmax, nvars, &minswapobj, &minswapsolvable) );
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "minswap-LP solved: obj = %g\n", minswapobj);
#endif
      }
   }

   // perform bound tightening, infeasibility checks and redundancy checks
   if( maxinfs <= 1 && (maxsolvable || maxswapsolvable) )
   {
      SCIP_Real activity;

      if( maxsolvable && maxswapsolvable )
         activity = MAX(maxobj, maxswapobj) + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/
      else if( maxsolvable )
         activity = maxobj + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/
      else
         activity = maxswapobj + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/

      // infeasibility check
      if( maxinfs == 0 && SCIPisPositive(scip, activity) )
      {
         (*infeasible) = TRUE;
         (*success) = TRUE;
         return SCIP_OKAY;
      }

      // strengthen bounds of all variables outside overlap
      else if( maxinfs == 0 )
      {
         for( i = 0; i < row1len; i++ )
         {
            idx1 = row1idxptr[i];
            if( cangetbnd[idx1] )
            {
               if( SCIPisPositive(scip, row1valptr[i]) )
               {
                  if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                     newbnd = SCIPceil(scip, (activity + row1valptr[i] * ubs[idx1]) / row1valptr[i]);
                  else
                     newbnd = (activity + row1valptr[i] * ubs[idx1]) / row1valptr[i];

                  if( SCIPisGT(scip, newbnd, lbs[idx1]) )
                  {
#ifdef SCIP_DEBUG_BOUNDS
                     SCIPdebugMsg(scip, "%g <= %g <= %s <= %g\n", lbs[idx1], newbnd, SCIPmatrixGetColName(matrix, idx1), ubs[idx1]);
#endif
                     lbs[idx1] = newbnd;
                     (*success) = TRUE;
                  }
               }
               else
               {
                  assert(SCIPisNegative(scip, row1valptr[i]));
                  if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                     newbnd = SCIPfloor(scip, (activity + row1valptr[i] * lbs[idx1]) / row1valptr[i]);
                  else
                     newbnd = (activity + row1valptr[i] * lbs[idx1]) / row1valptr[i];

                  if( SCIPisLT(scip, newbnd, ubs[idx1]) )
                  {
#ifdef SCIP_DEBUG_BOUNDS
                     SCIPdebugMsg(scip, "%g <= %s <= %g <= %g\n", lbs[idx1], SCIPmatrixGetColName(matrix, idx1), newbnd, ubs[idx1]);
#endif
                     ubs[idx1] = newbnd;
                     (*success) = TRUE;
                  }
               }
            }
         }
      }
      // strengthen bound of the single variable contributing the infinity
      else
      {
         assert(maxinfs == 1);
         for( i = 0; i < row1len; i++ )
         {
            idx1 = row1idxptr[i];
            if( cangetbnd[idx1] )
            {
               if( SCIPisPositive(scip, row1valptr[i]) && SCIPisInfinity(scip, ubs[idx1]) )
               {
                  if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                     newbnd = SCIPceil(scip, activity / row1valptr[i]);
                  else
                     newbnd = activity / row1valptr[i];

                  if( SCIPisGT(scip, newbnd, lbs[idx1]) )
                  {
#ifdef SCIP_DEBUG_BOUNDS
                     SCIPdebugMsg(scip, "%g <= %g <= %s <= %g\n", lbs[idx1], newbnd, SCIPmatrixGetColName(matrix, idx1), ubs[idx1]);
#endif
                     lbs[idx1] = newbnd;
                     (*success) = TRUE;
                  }
               }
               else if( SCIPisInfinity(scip, -lbs[idx1]) )
               {
                  assert(SCIPisNegative(scip, row1valptr[i]));
                  if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                      || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                     newbnd = SCIPfloor(scip, activity / row1valptr[i]);
                  else
                     newbnd = activity / row1valptr[i];

                  if( SCIPisLT(scip, newbnd, ubs[idx1]) )
                  {
#ifdef SCIP_DEBUG_BOUNDS
                     SCIPdebugMsg(scip, "%g <= %s <= %g <= %g\n", lbs[idx1], SCIPmatrixGetColName(matrix, idx1), newbnd, ubs[idx1]);
#endif
                     ubs[idx1] = newbnd;
                     (*success) = TRUE;
                  }
               }
            }
         }
      }
   }

   // redundancy check
   if( mininfs == 0 && !swaprow1 )
   {
      if( (minsolvable && SCIPisGT(scip, minobj, SCIPmatrixGetRowLhs(matrix, row1idx) + minact))
          || (minswapsolvable && SCIPisGT(scip, minswapobj, SCIPmatrixGetRowLhs(matrix, row1idx) + minact)) ) /*lint !e644*/
         (*redundant) = TRUE;
      else
         (*redundant) = FALSE;
   }
   else
      (*redundant) = FALSE;

   /* in this case the objective is swapped. therefore the minimum and the maximum of the support switch roles */
   if( swaprow1 )
   {
      // perform bound tightening, infeasibility checks and redundancy checks
      if( mininfs <= 1 && (minsolvable || minswapsolvable) )
      {
         SCIP_Real activity;

         if( minsolvable && minswapsolvable )
            activity = MAX(minobj, minswapobj) - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;
         else if( minsolvable )
            activity = minobj - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;
         else
            activity = minswapobj - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;

         // infeasibility check
         if( mininfs == 0 && SCIPisPositive(scip, activity) )
         {
            (*infeasible) = TRUE;
            (*success) = TRUE;
            return SCIP_OKAY;
         }
         // strengthen bounds of all variables outside overlap
         else if( mininfs == 0 )
         {
            for( i = 0; i < row1len; i++ )
            {
               idx1 = row1idxptr[i];
               if( cangetbnd[idx1] )
               {
                  if( SCIPisNegative(scip, row1valptr[i]) ) // since we look at the swapped case, this represents a positive coefficient
                  {
                     if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                        newbnd = SCIPceil(scip, (activity - row1valptr[i] * ubs[idx1]) / (-row1valptr[i]));
                     else
                        newbnd = (activity - row1valptr[i] * ubs[idx1]) / (-row1valptr[i]);

                     if( SCIPisGT(scip, newbnd, lbs[idx1]) )
                     {
#ifdef SCIP_DEBUG_BOUNDS
                        SCIPdebugMsg(scip, "%g <= %g <= %s <= %g\n", lbs[idx1], newbnd, SCIPmatrixGetColName(matrix, idx1), ubs[idx1]);
#endif
                        lbs[idx1] = newbnd;
                        (*success) = TRUE;
                     }
                  }
                  else
                  {
                     assert(SCIPisPositive(scip, row1valptr[i])); // since we look at the swapped case, this represents a negative coefficient
                     if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                        newbnd = SCIPfloor(scip, (activity - row1valptr[i] * lbs[idx1]) / (-row1valptr[i]));
                     else
                        newbnd = (activity - row1valptr[i] * lbs[idx1]) / (-row1valptr[i]);

                     if( SCIPisLT(scip, newbnd, ubs[idx1]) )
                     {
#ifdef SCIP_DEBUG_BOUNDS
                        SCIPdebugMsg(scip, "%g <= %s <= %g <= %g\n", lbs[idx1], SCIPmatrixGetColName(matrix, idx1), newbnd, ubs[idx1]);
#endif
                        ubs[idx1] = newbnd;
                        (*success) = TRUE;
                     }
                  }
               }
            }
         }
         // strengthen bound of the single variable contributing the infinity
         else
         {
            assert(mininfs == 1);
            for( i = 0; i < row1len; i++ )
            {
               idx1 = row1idxptr[i];
               if( cangetbnd[idx1] )
               {
                  if( SCIPisNegative(scip, row1valptr[i]) && SCIPisInfinity(scip, ubs[idx1]) ) // since we look at the swapped case, this represents a positive coefficient
                  {
                     if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                        newbnd = SCIPceil(scip, activity / (-row1valptr[i]));
                     else
                        newbnd = activity / (-row1valptr[i]);

                     if( SCIPisGT(scip, newbnd, lbs[idx1]) )
                     {
#ifdef SCIP_DEBUG_BOUNDS
                        SCIPdebugMsg(scip, "%g <= %g <= %s <= %g\n", lbs[idx1], newbnd, SCIPmatrixGetColName(matrix, idx1), ubs[idx1]);
#endif
                        lbs[idx1] = newbnd;
                        (*success) = TRUE;
                     }
                  }
                  else if( SCIPisInfinity(scip, -lbs[idx1]) )
                  {
                     assert(SCIPisPositive(scip, row1valptr[i])); // since we look at the swapped case, this represents a negative coefficient
                     if( SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_BINARY
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_INTEGER
                         || SCIPvarGetType(SCIPmatrixGetVar(matrix, idx1)) == SCIP_VARTYPE_IMPLINT )
                        newbnd = SCIPfloor(scip, activity / (-row1valptr[i]));
                     else
                        newbnd = activity / (-row1valptr[i]);

                     if( SCIPisLT(scip, newbnd, ubs[idx1]) )
                     {
#ifdef SCIP_DEBUG_BOUNDS
                        SCIPdebugMsg(scip, "%g <= %s <= %g <= %g\n", lbs[idx1], SCIPmatrixGetColName(matrix, idx1), newbnd, ubs[idx1]);
#endif
                        ubs[idx1] = newbnd;
                        (*success) = TRUE;
                     }
                  }
               }
            }
         }
      }
   }

   return SCIP_OKAY;
}


static
SCIP_RETCODE twoRowBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   row1,
   int                   row2,
   SCIP_Bool             swaprow1,
   SCIP_Bool             swaprow2,
   SCIP_Real*            lbs,
   SCIP_Real*            ubs,
   SCIP_Bool*            delcons,
   SCIP_Bool*            success             /**< return (success || "found better bounds") */
   )
{
   SCIP_Real* amin;
   SCIP_Real* amax;
   SCIP_Real* cmin;
   SCIP_Real* cmax;
   SCIP_Real* newlbsmin;
   SCIP_Real* newlbsmax;
   SCIP_Real* newubsmin;
   SCIP_Real* newubsmax;
   SCIP_Bool* cangetbnd;
   SCIP_Bool row1redundant;
   SCIP_Bool row2redundant;
   SCIP_Bool infeasible;

#ifdef SCIP_DEBUG_2RB
   SCIPdebugMsg(scip, "combining rows %d (%s) and %d (%s)\n", row1, SCIPmatrixGetRowName(matrix, row1), row2, SCIPmatrixGetRowName(matrix, row2));
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &amin, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &amax, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cmin, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cmax, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlbsmin, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlbsmax, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newubsmin, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newubsmax, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cangetbnd, SCIPmatrixGetNColumns(matrix)) );

   // Sort matrix rows
   SCIPsortIntReal(SCIPmatrixGetRowIdxPtr(matrix, row1), SCIPmatrixGetRowValPtr(matrix, row1), SCIPmatrixGetRowNNonzs(matrix, row1));
   SCIPsortIntReal(SCIPmatrixGetRowIdxPtr(matrix, row2), SCIPmatrixGetRowValPtr(matrix, row2), SCIPmatrixGetRowNNonzs(matrix, row2));

   // Use row2 to strengthen row1
   infeasible = FALSE;
   SCIP_CALL( transformAndSolve(scip, matrix, row1, row2, swaprow1, swaprow2, amin, amax, cmin, cmax, cangetbnd, lbs, ubs, newlbsmin, newlbsmax, newubsmin, newubsmax, success, &row1redundant, &infeasible) );

   if( row1redundant )
   {
      delcons[row1] = TRUE;
   }

   // Switch roles and use row1 to strengthen row2
   SCIP_CALL( transformAndSolve(scip, matrix, row2, row1, swaprow2, swaprow1, amin, amax, cmin, cmax, cangetbnd, lbs, ubs, newlbsmin, newlbsmax, newubsmin, newubsmax, success, &row2redundant, &infeasible) );

   if( row2redundant && !row1redundant )
   {
      delcons[row2] = TRUE;
   }

   SCIPfreeBufferArray(scip, &cangetbnd);
   SCIPfreeBufferArray(scip, &newubsmax);
   SCIPfreeBufferArray(scip, &newubsmin);
   SCIPfreeBufferArray(scip, &newlbsmax);
   SCIPfreeBufferArray(scip, &newlbsmin);
   SCIPfreeBufferArray(scip, &cmax);
   SCIPfreeBufferArray(scip, &cmin);
   SCIPfreeBufferArray(scip, &amax);
   SCIPfreeBufferArray(scip, &amin);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeTworowbnd)
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
SCIP_DECL_PRESOLINIT(presolInitTworowbnd)
{
   SCIP_PRESOLDATA* presoldata;

   presoldata = SCIPpresolGetData(presol);
   presoldata->nchgbnds = 0;
   presoldata->nuselessruns = 0;

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecTworowbnd)
{  /*lint --e{715}*/

   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;
   SCIP_Bool infeasible;
   SCIP_PRESOLDATA* presoldata;

   int i;
   int j;
   int k;

   SCIPinfoMessage(scip, NULL, "starting tworowbnd\n");

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
      int nrows;
      int ncols;
      SCIP_Real* oldlbs;
      SCIP_Real* oldubs;
      SCIP_Real* newlbs;
      SCIP_Real* newubs;
      SCIP_Bool* delcons;
      int* rowidxptr;
      SCIP_Real* rowvalptr;
      SCIP_VAR* var;

      int maxhashes;

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
      // prevent overflow issues
      if( nrows != 0 && maxhashes / nrows != presoldata->maxhashfac )
      {
         maxhashes = INT_MAX;
      }

      for( i = 0; i < nrows; i++)
      {
         if( pospp + posmm + pospm + posmp > maxhashes )
            break;

         rowvalptr = SCIPmatrixGetRowValPtr(matrix, i);
         rowidxptr = SCIPmatrixGetRowIdxPtr(matrix, i);
         finiterhs = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, i));
         maxlen = MIN(presoldata->maxconsiderednonzeros, SCIPmatrixGetRowNNonzs(matrix, i)); /*lint !e666*/
         for( j = 0; j < maxlen; j++)
         {
            for( k = j+1; k < maxlen; k++)
            {
               if( SCIPisPositive(scip, rowvalptr[j]) )
               {
                  if(SCIPisPositive(scip, rowvalptr[k]) )
                  {
                     // right-shift is required because we want to sort the hashes later on
                     SCIP_CALL( addEntry(scip, &pospp, &listsizepp, &hashlistpp, &rowidxlistpp,
                        hashIndexPair(rowidxptr[j],rowidxptr[k])>>1, i) ); /*lint !e702*/
                     if( finiterhs )
                        SCIP_CALL( addEntry(scip, &posmm, &listsizemm, &hashlistmm, &rowidxlistmm,
                           hashIndexPair(rowidxptr[j],rowidxptr[k])>>1, i) ); /*lint !e702*/
                  }
                  else
                  {
                     SCIP_CALL( addEntry(scip, &pospm, &listsizepm, &hashlistpm, &rowidxlistpm,
                        hashIndexPair(rowidxptr[j],rowidxptr[k])>>1, i) ); /*lint !e702*/
                     if( finiterhs )
                        SCIP_CALL( addEntry(scip, &posmp, &listsizemp, &hashlistmp, &rowidxlistmp,
                           hashIndexPair(rowidxptr[j],rowidxptr[k])>>1, i) ); /*lint !e702*/
                  }
               }
               else
               {
                  if(SCIPisPositive(scip, rowvalptr[k]) )
                  {
                     // right-shift is required because we want to sort the hashes later on
                     SCIP_CALL( addEntry(scip, &posmp, &listsizemp, &hashlistmp, &rowidxlistmp,
                        hashIndexPair(rowidxptr[j],rowidxptr[k])>>1, i) ); /*lint !e702*/
                     if( finiterhs )
                        SCIP_CALL( addEntry(scip, &pospm, &listsizepm, &hashlistpm, &rowidxlistpm,
                           hashIndexPair(rowidxptr[j],rowidxptr[k])>>1, i) ); /*lint !e702*/
                  }
                  else
                  {
                     SCIP_CALL( addEntry(scip, &posmm, &listsizemm, &hashlistmm, &rowidxlistmm,
                        hashIndexPair(rowidxptr[j],rowidxptr[k])>>1, i) ); /*lint !e702*/
                     if( finiterhs )
                        SCIP_CALL( addEntry(scip, &pospp, &listsizepp, &hashlistpp, &rowidxlistpp,
                           hashIndexPair(rowidxptr[j],rowidxptr[k])>>1, i) ); /*lint !e702*/
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

      SCIP_CALL( SCIPallocCleanBufferArray(scip, &delcons, nrows) );
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
         int ncombines;
         int maxcombines;
         SCIP_Bool finished;
         SCIP_Bool success;
         SCIP_Bool swaprow1;
         SCIP_Bool swaprow2;
         int combinefails;
         int retrievefails;
         ROWPAIR rowpair;

         finished = FALSE;
         block1start = 0;
         block1end = 0;
         block2start = 0;
         block2end = 0;
         maxcombines = nrows * presoldata->maxpairfac;
         // prevent overflow issues
         if( nrows != 0 && maxcombines / nrows != presoldata->maxpairfac )
         {
            maxcombines = INT_MAX;
         }
         ncombines = 0;
         combinefails = 0;
         retrievefails = 0;
         findNextBlock(hashlistpp, pospp, &block1start, &block1end);
         findNextBlock(hashlistmm, posmm, &block2start, &block2end);
         SCIPdebugMsg(scip, "processing pp and mm\n");
         while( !finished )
         {
            //SCIPdebugMsg(scip, "block1: %d -> %d, block2: %d -> %d, hashes: %d, %d\n", block1start, block1end, block2start, block2end, hashlistpp[block1start], hashlistmm[block2start]);
            if( hashlistpp[block1start] == hashlistmm[block2start] )
            {
               for( i = block1start; i < block1end; i++ )
               {
                  for( j = block2start; j < block2end; j++ )
                  {
                     if( rowidxlistpp[i] != rowidxlistmm[j] )
                     {
                        rowpair.row1idx = MIN(rowidxlistpp[i], rowidxlistmm[j]);
                        rowpair.row2idx = MAX(rowidxlistpp[i], rowidxlistmm[j]);
                        if( ! SCIPhashsetExists(pairhashset, encodeRowPair(&rowpair)) )
                        {
                           assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row1idx)));
                           assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row2idx)));

                           success = FALSE;

                           swaprow1 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx));
                           swaprow2 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx));

                           SCIP_CALL( twoRowBound(scip, matrix, rowpair.row1idx, rowpair.row2idx, swaprow1, swaprow2, newlbs, newubs, delcons, &success) );

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
               }
               else
                  finished = TRUE;
            }
            else if( hashlistpp[block1start] < hashlistmm[block2start] && block1end < pospp )
               findNextBlock(hashlistpp, pospp, &block1start, &block1end);
            else if( hashlistpp[block1start] > hashlistmm[block2start] && block2end < posmm )
               findNextBlock(hashlistmm, posmm, &block2start, &block2end);
            else
               finished = TRUE;
         }
      }

      /* Process pm and mp lists */
      if( pospm > 0 && posmp > 0 )
      {
         int ncombines;
         int maxcombines;
         SCIP_Bool finished;
         SCIP_Bool success;
         SCIP_Bool swaprow1;
         SCIP_Bool swaprow2;
         int combinefails;
         int retrievefails;
         ROWPAIR rowpair;

         finished = FALSE;
         block1start = 0;
         block1end = 0;
         block2start = 0;
         block2end = 0;
         maxcombines = nrows * presoldata->maxpairfac;
         // prevent overflow issues
         if( nrows != 0 && maxcombines / nrows != presoldata->maxpairfac )
         {
            maxcombines = INT_MAX;
         }
         ncombines = 0;
         combinefails = 0;
         retrievefails = 0;
         findNextBlock(hashlistpm, pospm, &block1start, &block1end);
         findNextBlock(hashlistmp, posmp, &block2start, &block2end);
         SCIPdebugMsg(scip, "processing pm and mp\n");
         while( !finished )
         {
            if( hashlistpm[block1start] == hashlistmp[block2start] )
            {
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

                           swaprow1 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx));
                           swaprow2 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx));

                           SCIP_CALL( twoRowBound(scip, matrix, rowpair.row1idx, rowpair.row2idx, swaprow1, swaprow2, newlbs, newubs, delcons, &success) );

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
               }
               else
                  finished = TRUE;
            }
            else if( hashlistpm[block1start] < hashlistmp[block2start] && block1end < pospm )
               findNextBlock(hashlistpm, pospm, &block1start, &block1end);
            else if( hashlistpm[block1start] > hashlistmp[block2start] && block2end < posmp )
               findNextBlock(hashlistmp, posmp, &block2start, &block2end);
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

      /* Remove redundant constraints */
      ndelcons = 0;
      for( i = 0; i < nrows; i++ )
      {
         if( delcons[i] )
         {
            ndelcons++;
            delcons[i] = FALSE; //Remove nonzero from clean buffer array
            SCIPdebugMsg(scip, "removing redundant constraint %s\n", SCIPmatrixGetRowName(matrix, i));
            SCIP_CALL( SCIPdelCons(scip, SCIPmatrixGetCons(matrix, i)) );
         }
      }
      (*ndelconss) += ndelcons;

      /* check for redundant constraints */
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
            if( SCIPisLE(scip, SCIPmatrixGetRowLhs(matrix, i), rowinf) )
            {
               //TODO remove constraint
               //SCIPdebugMsg(scip, "%s is redundant\n", SCIPmatrixGetRowName(matrix, i));
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
            if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, i)) && SCIPisGE(scip, SCIPmatrixGetRowRhs(matrix, i), rowsup) )
            {
               //TODO remove constraint
               //SCIPdebugMsg(scip, "%s is redundant\n", SCIPmatrixGetRowName(matrix, i));
            }
         }
      }

      /* set result */
      if( *nchgbds > oldnchgbds || *nfixedvars > oldnfixedvars || ndelcons > 0 )
      {
         *result = SCIP_SUCCESS;
         presoldata->nuselessruns = 0;
         SCIPinfoMessage(scip, NULL, "tworowbnd evaluated %d pairs to tighten %d bounds, fix %d variables and delete %d redundant constraints\n", SCIPhashsetGetNElements(pairhashset), *nchgbds - oldnchgbds, *nfixedvars - oldnfixedvars, ndelcons);
      }
      else if( infeasible )
      {
         *result = SCIP_CUTOFF;
         SCIPinfoMessage(scip, NULL, "tworowbnd detected infeasibility\n");
      }
      else
      {
         presoldata->nuselessruns++;
         SCIPinfoMessage(scip, NULL, "tworowbnd evaluated %d pairs without success\n", SCIPhashsetGetNElements(pairhashset));
      }

      SCIPhashsetFree(&pairhashset, SCIPblkmem(scip));
      SCIPfreeBufferArray(scip, &newubs);
      SCIPfreeBufferArray(scip, &newlbs);
      SCIPfreeBufferArray(scip, &oldubs);
      SCIPfreeBufferArray(scip, &oldlbs);
      SCIPfreeCleanBufferArray(scip, &delcons);
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

/** creates the tworowbndb presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolTworowbnd(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create tworowbnd presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );
   /* TODO: (optional) create presolver specific data here */

   presol = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_TIMING,
         presolExecTworowbnd,
         presoldata) );

   assert(presol != NULL);

   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeTworowbnd) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitTworowbnd) );

   /* add tworowbnd presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxconsiderednonzeros",
         "maximal number of considered non-zeros within one row (-1: no limit)",
         &presoldata->maxconsiderednonzeros, FALSE, DEFAULT_MAXCONSIDEREDNONZEROS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxretrievefails",
         "maximal number of consecutive useless hashtable retrieves",
         &presoldata->maxretrievefails, FALSE, DEFAULT_MAXRETRIEVEFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxcombinefails",
         "maximal number of consecutive useless row combines",
         &presoldata->maxcombinefails, FALSE, DEFAULT_MAXCOMBINEFAILS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxhashfac",
         "Maximum number of hashlist entries as multiple of number of rows in the problem (-1: no limit)",
         &presoldata->maxhashfac, FALSE, DEFAULT_MAXHASHFAC, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/tworowbnd/maxpairfac",
         "Maximum number of processed row pairs as multiple of the number of rows in the problem (-1: no limit)",
         &presoldata->maxpairfac, FALSE, DEFAULT_MAXPAIRFAC, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
