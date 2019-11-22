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
 * Perform bound tightening on two inequalities with some common variables.
 * Two possible methods are being used.
 *
 * 1. LP-bound
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
 *
 *
 * 2. ConvComb with clique-extension
 * Given two constraints
 * \begin{equation}
 *   \begin{array}{c}
 *     A_{r\cdot} x \geq b_r \\
 *     A_{s\cdot} x \geq b_s \\
 *     \ell \leq x \leq u \\
 *   \end{array}
 * \end{equation}
 * this method determines promising values for $\lambda in (0,1)$ and
 * applies feasibility-based bound tightening on the convex combinations
 *
 * (\lambda A_{r\cdot} + (1 - \lambda) A_{s\cdot}) x \geq \lambda b_r + (1 - \lambda) b_s$.
 *
 * Additionally, cliques drawn from the SCIPcliqueTable are used
 * to further strengthen the above bound tightening.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*
 * Additional debug defines in this presolver
 * SCIP_DEBUG_HASHING
 * SCIP_DEBUG_BREAKPOINTS
 * SCIP_DEBUG_CLIQUE
 * SCIP_DEBUG_BOUNDS
 * SCIP_DEBUG_SINGLEROWLP
 */

#include <assert.h>

#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"
#include "scip/pub_matrix.h"
#include "presol_tworowbnd.h"
#include <string.h>

#define PRESOL_NAME                    "tworowbnd"
#define PRESOL_DESC                    "do bound tigthening by using two rows"
#define PRESOL_PRIORITY                -2000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS               -1       /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING                  SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_ENABLECOPY             TRUE     /**< should tworowbnd presolver be copied to sub-SCIPs? */
#define DEFAULT_MAXCONSIDEREDNONZEROS  100      /**< maximal number of considered non-zeros within one row (-1: no limit) */
#define DEFAULT_MAXRETRIEVEFAILS       1000     /**< maximal number of consecutive useless hashtable retrieves */
#define DEFAULT_MAXCOMBINEFAILS        1000     /**< maximal number of consecutive useless row combines */
#define DEFAULT_MAXHASHFAC             10       /**< maximal number of hashlist entries as multiple of number of rows in the problem (-1: no limit) */
#define DEFAULT_MAXPAIRFAC             1        /**< maximal number of processed row pairs as multiple of the number of rows in the problem (-1: no limit) */
#define DEFAULT_LPBOUND                FALSE    /**< should the bounds be tightened via LP-bound? */
#define DEFAULT_CONVCOMB               TRUE     /**< should the bounds be tightened via feasibility based bound-tightening on convex combined constraints? */

/*
 * Data structures
 */

/** Signum for convex-combined variable coefficients (\lambda * A_{ri} + (1 - \lambda) * A_{si})
 *  UP  - Coefficient changes from negative to positive for increasing lambda
 *  DN  - Coefficient changes from positive to negative for increasing lambda
 *  POS - Coefficient is positive for all lambda in (0,1)
 *  NEG - Coefficient is negative for all lambda in (0,1)
 *  CLQ - Coefficient belongs to a variable which is part of a clique
 */
enum signum {UP, DN, POS, NEG, CLQ};

/** presolver data */
struct SCIP_PresolData
{
   int maxpairfac;            /**< maximal number of processed row pairs as multiple of the number of rows in the problem (-1: no limit) */
   int maxhashfac;            /**< maximal number of hashlist entries as multiple of number of rows in the problem (-1: no limit) */
   int maxretrievefails;      /**< maximal number of consecutive useless hashtable retrieves */
   int maxcombinefails;       /**< maximal number of consecutive useless row combines */
   int maxconsiderednonzeros; /**< maximal number of considered non-zeros within one row (-1: no limit) */
   int nchgbnds;              /**< number of variable bounds changed by this presolver */
   int nuselessruns;          /**< number of runs where this presolver did not apply any changes */
   SCIP_Bool lpbound;         /**< should the bounds be tightened via LP-bound? */
   SCIP_Bool convcomb;        /**< should the bounds be tightened via feasibility based bound-tightening on convex combined constraints? */
   SCIP_Bool enablecopy;      /**< should tworowbnd presolver be copied to sub-SCIPs? */
};

/** structure representing a pair of row indices; used for lookup in a hashtable */
struct RowPair
{
   int row1idx; /**< first row index */
   int row2idx; /**< second row index */
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

/** add hash/rowidx pair to hashlist/rowidxlist **/
static
SCIP_RETCODE addEntry
(
 SCIP* scip,                  /**< SCIP datastructure */
 int* pos,                    /**< position of last entry added */
 int* listsize,               /**< size of hashlist and rowidxlist */
 int** hashlist,              /**< block memory array containing hashes */
 int** rowidxlist,            /**< block memory array containing row indices */
 int hash,                    /**< hash to be inserted */
 int rowidx                   /**< row index to be inserted */
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

/** Within a sorted list, get next block with same value
 *  E.g. for [h1, h1, h1, h2, h2, h2, h2, h3,...] and end = 0
 *  returns start = 0, end = 3
 *  and on a second call with end = 3 on the same list
 *  returns start = 3, end = 7.
 **/
static
void findNextBlock
(
   int*                 list,    /**< list of integers */
   int                  len,     /**< length of list */
   int*                 start,   /**< variable to contain start index of found block */
   int*                 end      /**< variable to contain end index of found block */
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
 *  min c^T x
 *  s.t. a^T x >= b
 *  lbs <= x <= ubs
 *
 *  First, the problem is transformed such that
 *  SCIPselectWeightedReal() can be applied, which
 *  then solves the problem as a continuous knapsack
 *  in linear time.
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

/* Transforms rows as required, solves the individual single-row LPs and then tightens bounds */
static
SCIP_RETCODE transformAndSolve(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object, rows specified by row1idx/row2idx must be sorted */
   int                   row1idx,            /**< index of first row */
   int                   row2idx,            /**< index of second row */
   SCIP_Bool             swaprow1,           /**< should row1 <= rhs be used in addition to lhs <= row1 */
   SCIP_Bool             swaprow2,           /**< should row2 <= rhs be used in addition to lhs <= row2 */
   SCIP_Real*            aoriginal,          /**< buffer array for original constraint coefficients */
   SCIP_Real*            acopy,              /**< buffer array for coefficients adjusted to single-row LP to be solved */
   SCIP_Real*            coriginal,          /**< buffer array for original objective coefficients */
   SCIP_Real*            ccopy,              /**< buffer array for coefficients adjusted to single-row LP to be solved */
   SCIP_Bool*            cangetbnd,          /**< buffer array for flags of which variables a bound can be generated */
   SCIP_Real*            lbs,                /**< buffer array for lower bounds for single-row LP */
   SCIP_Real*            ubs,                /**< buffer array for upper bounds for single-row LP */
   SCIP_Real*            newlbsoriginal,     /**< buffer array for new lower bounds not adjusted to individual single-row LPs */
   SCIP_Real*            newlbscopy,         /**< buffer array for adjusted lower bounds */
   SCIP_Real*            newubsoriginal,     /**< buffer array for new upper bounds not adjusted to individual single-row LPs */
   SCIP_Real*            newubscopy,         /**< buffer array for adjusted upper bounds */
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
         coriginal[nvars] = row1valptr[i];
         aoriginal[nvars] = row2valptr[j];
         newlbsoriginal[nvars] = lbs[idx1];
         newubsoriginal[nvars] = ubs[idx1];
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
         coriginal[nvars] = 0.0;
         aoriginal[nvars] = row2valptr[j];
         newlbsoriginal[nvars] = lbs[idx2];
         newubsoriginal[nvars] = ubs[idx2];
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
      coriginal[nvars] = 0.0;
      aoriginal[nvars] = row2valptr[j];
      newlbsoriginal[nvars] = lbs[idx2];
      newubsoriginal[nvars] = ubs[idx2];
      nvars++;
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "%g <= (%s) <= %g has coefs 0.0 and %g, %d LP vars\n", lbs[idx2], SCIPvarGetName(SCIPmatrixGetVar(matrix, idx2)), ubs[idx2], row2valptr[j], nvars);
#endif
      j++;
   }

#ifdef SCIP_DEBUG_2RB
   SCIPdebugMsg(scip, "right hand sides: %g and %g\n", SCIPmatrixGetRowLhs(matrix, row1idx), SCIPmatrixGetRowLhs(matrix, row2idx));
#endif

   /* solve single-row LPs */
   maxsolvable = FALSE;
   minsolvable = FALSE;
   maxswapsolvable = FALSE;
   minswapsolvable = FALSE;
   /* maximize overlap in first row with lhs <= row2 as constraint */
   if( maxinfs <= 1 )
   {
      for( i = 0; i < nvars; i++ )
      {
         acopy[i] = aoriginal[i];
         ccopy[i] = -coriginal[i];
         newlbscopy[i] = newlbsoriginal[i];
         newubscopy[i] = newubsoriginal[i];
      }
      SCIP_CALL( solveSingleRowLP(scip, acopy, SCIPmatrixGetRowLhs(matrix, row2idx), ccopy, newlbscopy, newubscopy, nvars, &maxobj, &maxsolvable) );
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "max-LP solved: obj = %g\n", maxobj);
#endif
   }

   /* minimize overlap in first row with lhs <= row2 as constraint */
   if( mininfs == 0 || (mininfs == 1 && swaprow1) )
   {
      /* copy coefficients */
      for( i = 0; i < nvars; i++ )
      {
         acopy[i] = aoriginal[i];
         ccopy[i] = coriginal[i];
         newlbscopy[i] = newlbsoriginal[i];
         newubscopy[i] = newubsoriginal[i];
      }
      SCIP_CALL( solveSingleRowLP(scip, acopy, SCIPmatrixGetRowLhs(matrix, row2idx), ccopy, newlbscopy, newubscopy, nvars, &minobj, &minsolvable) );
#ifdef SCIP_DEBUG_2RB
      SCIPdebugMsg(scip, "min-LP solved: obj = %g\n", minobj);
#endif
   }

   if( swaprow2 )
   {
     /* maximize overlap in first row with row2 <= rhs as constraint */
      if( maxinfs <= 1 )
      {
         /* copy coefficients */
         for( i = 0; i < nvars; i++ )
         {
            acopy[i] = -aoriginal[i];
            ccopy[i] = -coriginal[i];
            newlbscopy[i] = newlbsoriginal[i];
            newubscopy[i] = newubsoriginal[i];
         }
         SCIP_CALL( solveSingleRowLP(scip, acopy, -SCIPmatrixGetRowRhs(matrix, row2idx), ccopy, newlbscopy, newubscopy, nvars, &maxswapobj, &maxswapsolvable) );
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "maxswap-LP solved: obj = %g\n", maxswapobj);
#endif
      }

      /* minimize overlap in first row with row2 <= rhs as constraint */
      if( mininfs == 0 || (mininfs == 1 && swaprow1) )
      {
         /* copy coefficients */
         for( i = 0; i < nvars; i++ )
         {
            acopy[i] = -aoriginal[i];
            ccopy[i] = coriginal[i];
            newlbscopy[i] = newlbsoriginal[i];
            newubscopy[i] = newubsoriginal[i];
         }
         SCIP_CALL( solveSingleRowLP(scip, acopy, -SCIPmatrixGetRowRhs(matrix, row2idx), ccopy, newlbscopy, newubscopy, nvars, &minswapobj, &minswapsolvable) );
#ifdef SCIP_DEBUG_2RB
         SCIPdebugMsg(scip, "minswap-LP solved: obj = %g\n", minswapobj);
#endif
      }
   }

   /* perform bound tightening, infeasibility checks and redundancy checks */
   if( maxinfs <= 1 && (maxsolvable || maxswapsolvable) )
   {
      SCIP_Real activity;

      if( maxsolvable && maxswapsolvable )
         activity = MAX(maxobj, maxswapobj) + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/
      else if( maxsolvable )
         activity = maxobj + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/
      else
         activity = maxswapobj + SCIPmatrixGetRowLhs(matrix, row1idx) + maxact; /*lint !e644*/

      /* infeasibility check */
      if( maxinfs == 0 && SCIPisPositive(scip, activity) )
      {
         (*infeasible) = TRUE;
         (*success) = TRUE;
         return SCIP_OKAY;
      }

      /* strengthen bounds of all variables outside overlap */
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
      /* strengthen bound of the single variable contributing the infinity */
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

   /* redundancy check */
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
      /* perform bound tightening, infeasibility checks and redundancy checks */
      if( mininfs <= 1 && (minsolvable || minswapsolvable) )
      {
         SCIP_Real activity;

         if( minsolvable && minswapsolvable )
            activity = MAX(minobj, minswapobj) - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;
         else if( minsolvable )
            activity = minobj - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;
         else
            activity = minswapobj - SCIPmatrixGetRowRhs(matrix, row1idx) - minact;

         /* infeasibility check */
         if( mininfs == 0 && SCIPisPositive(scip, activity) )
         {
            (*infeasible) = TRUE;
            (*success) = TRUE;
            return SCIP_OKAY;
         }
         /* strengthen bounds of all variables outside overlap */
         else if( mininfs == 0 )
         {
            for( i = 0; i < row1len; i++ )
            {
               idx1 = row1idxptr[i];
               if( cangetbnd[idx1] )
               {
                  if( SCIPisNegative(scip, row1valptr[i]) ) /* since we look at the swapped case, this represents a positive coefficient */
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
                     assert(SCIPisPositive(scip, row1valptr[i])); /* since we look at the swapped case, this represents a negative coefficient */
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
         /* strengthen bound of the single variable contributing the infinity */
         else
         {
            assert(mininfs == 1);
            for( i = 0; i < row1len; i++ )
            {
               idx1 = row1idxptr[i];
               if( cangetbnd[idx1] )
               {
                  if( SCIPisNegative(scip, row1valptr[i]) && SCIPisInfinity(scip, ubs[idx1]) ) /* since we look at the swapped case, this represents a positive coefficient */
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
                     assert(SCIPisPositive(scip, row1valptr[i])); /* since we look at the swapped case, this represents a negative coefficient */
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


/* Create required buffer arrays and apply LP-based bound tightening in both directions */
static
SCIP_RETCODE applyLPboundTightening(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_MATRIX*          matrix,             /**< constraint matrix object */
   int                   row1,               /**< index of first row */
   int                   row2,               /**< index of seond row */
   SCIP_Bool             swaprow1,           /**< should row1 <= rhs be used in addition to lhs <= row1 */
   SCIP_Bool             swaprow2,           /**< should row2 <= rhs be used in addition to lhs <= row2 */
   SCIP_Real*            lbs,                /**< lower variable bounds */
   SCIP_Real*            ubs,                /**< upper variable bounds */
   SCIP_Bool*            delcons,            /**< flags which constraints are redundant and can be removed */
   SCIP_Bool*            success             /**< return (success || "found better bounds") */
   )
{
   SCIP_Real* aoriginal;
   SCIP_Real* acopy;
   SCIP_Real* coriginal;
   SCIP_Real* ccopy;
   SCIP_Real* newlbsoriginal;
   SCIP_Real* newlbscopy;
   SCIP_Real* newubsoriginal;
   SCIP_Real* newubscopy;
   SCIP_Bool* cangetbnd;
   SCIP_Bool row1redundant;
   SCIP_Bool row2redundant;
   SCIP_Bool infeasible;

#ifdef SCIP_DEBUG_2RB
   SCIPdebugMsg(scip, "combining rows %d (%s) and %d (%s)\n", row1, SCIPmatrixGetRowName(matrix, row1), row2, SCIPmatrixGetRowName(matrix, row2));
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &aoriginal, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &acopy, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coriginal, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ccopy, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlbsoriginal, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newlbscopy, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newubsoriginal, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newubscopy, SCIPmatrixGetNColumns(matrix)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cangetbnd, SCIPmatrixGetNColumns(matrix)) );

   /* Sort matrix rows */
   SCIPsortIntReal(SCIPmatrixGetRowIdxPtr(matrix, row1), SCIPmatrixGetRowValPtr(matrix, row1), SCIPmatrixGetRowNNonzs(matrix, row1));
   SCIPsortIntReal(SCIPmatrixGetRowIdxPtr(matrix, row2), SCIPmatrixGetRowValPtr(matrix, row2), SCIPmatrixGetRowNNonzs(matrix, row2));

   /* Use row2 to strengthen row1 */
   infeasible = FALSE;
   SCIP_CALL( transformAndSolve(scip, matrix, row1, row2, swaprow1, swaprow2, aoriginal, acopy, coriginal, ccopy, cangetbnd, lbs, ubs, newlbsoriginal, newlbscopy, newubsoriginal, newubscopy, success, &row1redundant, &infeasible) );

   if( row1redundant )
   {
      delcons[row1] = TRUE;
   }

   /* Switch roles and use row1 to strengthen row2 */
   SCIP_CALL( transformAndSolve(scip, matrix, row2, row1, swaprow2, swaprow1, aoriginal, acopy, coriginal, ccopy, cangetbnd, lbs, ubs, newlbsoriginal, newlbscopy, newubsoriginal, newubscopy, success, &row2redundant, &infeasible) );

   if( row2redundant && !row1redundant )
   {
      delcons[row2] = TRUE;
   }

   SCIPfreeBufferArray(scip, &cangetbnd);
   SCIPfreeBufferArray(scip, &newubscopy);
   SCIPfreeBufferArray(scip, &newubsoriginal);
   SCIPfreeBufferArray(scip, &newlbscopy);
   SCIPfreeBufferArray(scip, &newlbsoriginal);
   SCIPfreeBufferArray(scip, &ccopy);
   SCIPfreeBufferArray(scip, &coriginal);
   SCIPfreeBufferArray(scip, &acopy);
   SCIPfreeBufferArray(scip, &aoriginal);

   return SCIP_OKAY;
}

/* calculate clique breakpoints and mark in which cliques the variables have a maximum */
static
int calcCliqueMaximums(
   SCIP*                scip,
   int*                 varinds,             /**< variable index array */
   int*                 cliquevarpos,        /**< positions of clique variables in index array */
   int                  cliquesize,          /**< size of current clique */
   SCIP_Real*           row1coefs,           /**< coefficients of first row */
   SCIP_Real*           row2coefs,           /**< coefficients of second row */
   int*                 nbreakpoints,        /**< number of breakpoints between 0 and 1 */
   SCIP_Real*           breakpoints,         /**< variable breakpoints */
   int*                 cliquemaxinds,       /**< array containing in which clique this variable has a maximum */
   SCIP_Real*           gradients            /**< buffer array the function can use for handling the gradients */
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
   SCIP_Real breakpointval;

#ifdef SCIP_DEBUG_CLIQUE
   SCIPdebugMsg(scip, "calculating maximums for clique %d\n", -cliquemaxinds[cliquevarpos[0]]);
#endif

   /* calculate gradients */
   for( i = 0; i < cliquesize; i++ )
   {
      gradients[i] = row1coefs[varinds[cliquevarpos[i]]] - row2coefs[varinds[cliquevarpos[i]]];
   }

   SCIPsortRealInt(gradients, cliquevarpos, cliquesize);

#ifdef SCIP_DEBUG_CLIQUE
   for( i = 0; i < cliquesize; i++ )
      SCIPdebugMsg(scip, "var_%d: %g + %g * lambda\n", i, row2coefs[varinds[cliquevarpos[i]]], gradients[i]);
#endif


   /* find maximum for lambda = 0 */
   maxpos = 0;
   for( i = 1; i < cliquesize; i++ )
   {
      if( SCIPisGE(scip, row2coefs[varinds[cliquevarpos[i]]], row2coefs[varinds[cliquevarpos[maxpos]]]) )
         maxpos = i;
   }

   /* variable is relevant only if its coefficient is non-negative, if all coefs are negative we set firstmaxpos to -1 */
   if( SCIPisPositive(scip, row2coefs[varinds[cliquevarpos[maxpos]]]) || (SCIPisZero(scip, row2coefs[varinds[cliquevarpos[maxpos]]]) && SCIPisPositive(scip, gradients[maxpos])) )
      firstmaxpos = maxpos;
   else
      firstmaxpos = -1;

#ifdef SCIP_DEBUG_CLIQUE
   SCIPdebugMsg(scip, "first maxvar: var_%d at coef(0) = %g, firstmaxpos = %d\n", maxpos, row2coefs[varinds[cliquevarpos[maxpos]]], firstmaxpos);
#endif

   /* find next maximum */
   minlambda = 0.0;
   if( firstmaxpos == -1 )
      maxidx = -1;
   else
      maxidx = varinds[cliquevarpos[firstmaxpos]];

   while( maxpos < cliquesize && SCIPisLT(scip, minlambda, 1.0) )
   {
      newmaxpos = -1;

      /* all coefficients are negative */
      if( maxidx == -1 )
      {
         /* find variable with lowest non-negative lambda such that (row2coef - row1coef) * lambda + row2coef >= 0 */
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
            breakpoints[cliquevarpos[firstmaxpos]] = -row2coefs[maxidx] / gradients[maxpos]; /* lambda where old maxcoef hits zero */
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
            assert(firstmaxpos != -1); /* use firstmaxpos as breakpoint dummy variable for the negative segment */
            breakpoints[cliquevarpos[firstmaxpos]] = -row2coefs[maxidx] / gradients[maxpos]; /* lambda where old maxcoef hits zero */
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
   int                  row1,                /**< index of first row */
   int                  row2,                /**< index of second row */
   SCIP_Bool            swaprow1,            /**< should row1 <= rhs be used instead of lhs <= row1 */
   SCIP_Bool            swaprow2,            /**< should row2 <= rhs be used instead of lhs <= row2 */
   SCIP_Real*           lbs,                 /**< lower variable bounds */
   SCIP_Real*           ubs,                 /**< upper variable bounds */
   SCIP_Bool*           success              /**< we return (success ||  found better bounds") */
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
   SCIP_Real* gradients;

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

#ifdef SCIP_DEBUG
   /* swap = 0 means that lhs <= cons is used for bound tightening */
   SCIPdebugMsg(scip, "combining rows %d (%s) and %d (%s) with swaps %d/%d\n", row1, SCIPmatrixGetRowName(matrix, row1), row2, SCIPmatrixGetRowName(matrix, row2), swaprow1, swaprow2);
   SCIPinfoMessage(scip, NULL, "Cons1: \n");
   SCIPprintCons(scip, SCIPmatrixGetCons(matrix, row1), NULL);
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "Cons2: \n");
   SCIPprintCons(scip, SCIPmatrixGetCons(matrix, row2), NULL);
   SCIPinfoMessage(scip, NULL, "\n");
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
      /* We use 0 as default value for "not in any clique" */
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

#ifdef SCIP_DEBUG_BREAKPOINTS
   SCIPdebugMsg(scip, "breakpoints and signs before considering cliques:\n");
   for( i = 0; i < nvars; i++ )
   {
      idx = varinds[i];
      if( signs[idx] == UP )
         SCIPdebugMsg(scip, "%g <= %s <= %g, UP: %g to %g, breakpoint at %g\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
      else if( signs[idx] == DN )
         SCIPdebugMsg(scip, "%g <= %s <= %g, DN: %g to %g, breakpoint at %g\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
      else if( signs[idx] == POS )
         SCIPdebugMsg(scip, "%g <= %s <= %g, POS: %g to %g, breakpoint at %g\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
      else if( signs[idx] == NEG )
         SCIPdebugMsg(scip, "%g <= %s <= %g, NEG: %g to %g, breakpoint at %g\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
   }
#endif

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
#ifdef SCIP_DEBUG
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
   SCIP_CALL( SCIPallocBufferArray(scip, &gradients, nbinvars) );

   SCIP_CALL( SCIPcalcCliquePartition(scip, binvars, nbinvars, cliquepartition, &ncliques) );
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
            cliquepartition[j-shift] = cliquepartition[j] - shift + 1; /* +1 ensures that all clique IDs are >0 such that we can use -ID to signal an all-negative clique */
            cliquemaxinds[binvarpos[j]] = -cliquepartition[j-shift]; /* negative sign implies that this variable does not assume the maximum in this clique, will be adjusted in calcCliqueMaximums */
            /* variables in cliques have either a clique-breakpoint or no breakpoint at all */
            if( !SCIPisEQ(scip, breakpoints[binvarpos[j]], 2.0) )
            {
               breakpoints[binvarpos[j]] = 2.0;
               nbreakpoints--;
            }
         }

#ifdef SCIP_DEBUG_CLIQUE
         SCIPdebugMsg(scip, "breakpoints before checking maximums of clique %d: %d\n", cliquepartition[binvarpos[i]], nbreakpoints);
#endif
         /* size of current clique equals j - i */
         idx = cliquepartition[i-shift]; // we enumerated the relevant cliques only
         currentmaxinds[idx] = calcCliqueMaximums(scip, varinds, &binvarpos[i-shift], j - i, row1coefs, row2coefs, &nbreakpoints, breakpoints, cliquemaxinds, gradients);
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

#ifdef SCIP_DEBUG_CLIQUE
   nbinvars -= shift;
   SCIPdebugMsg(scip, "after clearing single variable cliques: ncliques = %d, nbinvars = %d, nbreakpoints = %d\n", ncliques, nbinvars, nbreakpoints);
   for( i = 0; i < nbinvars; i++ )
      SCIPdebugMsg(scip, "%s (%d) is in clique %d with breakpoint %g\n", SCIPmatrixGetColName(matrix, varinds[binvarpos[i]]), varinds[binvarpos[i]], cliquepartition[i], breakpoints[binvarpos[i]]);
#endif

#ifdef SCIP_DEBUG_BREAKPOINTS
   SCIPdebugMsg(scip, "breakpoints and signs after considering cliques\n");
   SCIPdebugMsg(scip, "b1 = %g, b2 = %g, nbreakpoints = %d\n", b1, b2, nbreakpoints);
   for( i = 0; i < nvars; i++ )
   {
      idx = varinds[i];
      if( signs[idx] == UP )
         SCIPdebugMsg(scip, "%g <= %s <= %g, UP: %g to %g, breakpoint at %g\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
      else if( signs[idx] == DN )
         SCIPdebugMsg(scip, "%g <= %s <= %g, DN: %g to %g, breakpoint at %g\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
      else if( signs[idx] == POS )
         SCIPdebugMsg(scip, "%g <= %s <= %g, POS: %g to %g, breakpoint at %g\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
      else if( signs[idx] == NEG )
         SCIPdebugMsg(scip, "%g <= %s <= %g, NEG: %g to %g, breakpoint at %g\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i]);
      else if( signs[idx] == CLQ )
         SCIPdebugMsg(scip, "%g <= %s <= %g, CLQ: %g to %g, breakpoint at %g, cliquemaxind = %d\n", lbs[i], SCIPmatrixGetColName(matrix, idx), ubs[i], row2coefs[idx], row1coefs[idx], breakpoints[i], cliquemaxinds[i]);
   }
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
      SCIP_CALL( SCIPsetIntParam(subscip, "presolving/tworowbnd/maxrounds", 0) );

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

#ifdef SCIP_DEBUG_CLIQUE
      SCIPdebugMsg(scip, "currentmaxinds[0] = dummy\n");
      for( i = 1; i < ncliques+1; i++ )
         if( currentmaxinds[i] >= 0 )
            SCIPdebugMsg(scip, "currentmaxinds[%d] = %d (%s)\n", i, currentmaxinds[i], SCIPmatrixGetColName(matrix, currentmaxinds[i]));
         else
            SCIPdebugMsg(scip, "currentmaxinds[%d] = %d (no maximal variable)\n", i, currentmaxinds[i]);
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

                  if( currentmaxinds[cliqueidx] == -1 )
                     newubs[i] = MIN(newubs[i], l2 / row2coefs[idx]);
                  else if( SCIPisLT(scip, (l2 + row2coefs[currentmaxinds[cliqueidx]]) / row2coefs[idx], 1.0) )
                     newubs[i] = MIN(newubs[i], 0.0);
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

      /* iterate over all breakpoints and apply bound tightening on each resulting convex combination */
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
         j = i;
         while( !updated )
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
            j++;
         }

#ifdef SCIP_DEBUG_BREAKPOINTS
         SCIPdebugMsg(scip, "lambda_%d = %g, l1 = %g, l2 = %g, ninfs = %d\n", i, breakpoints[i], l1, l2, ninfs);
#endif

#ifdef SCIP_DEBUG_CLIQUE
      SCIPdebugMsg(scip, "currentmaxinds[0] = dummy\n");
      for( j = 1; j < ncliques+1; j++ )
         if( currentmaxinds[j] >= 0 )
            SCIPdebugMsg(scip, "currentmaxinds[%d] = %d (%s)\n", j, currentmaxinds[j], SCIPmatrixGetColName(matrix, currentmaxinds[j]));
         else
            SCIPdebugMsg(scip, "currentmaxinds[%d] = %d (no maximal variable)\n", j, currentmaxinds[j]);
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
               /* catch the special case where the entire remaining constraint is cancelled */
               if( j >= nvars )
                  break;

#ifdef SCIP_DEBUG_BOUNDS
               oldlb = newlbs[j];
               oldub = newubs[j];
#endif

               idx = varinds[j];
               coef = breakpoints[i] * row1coefs[idx] + (1 - breakpoints[i]) * row2coefs[idx];
               assert(!SCIPisEQ(scip, breakpoints[i], 2.0));

               /* skip if the coefficient is too close to zero as it becomes numerically unstable */
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
                     else if(SCIPisLT(scip, (breakpoints[i] * (l1 + row1coefs[currentmaxinds[cliqueidx]]) + (1 - breakpoints[i]) * (l2 + row2coefs[currentmaxinds[cliqueidx]])) / coef, 1.0) )
                        newubs[j] = MIN(newubs[j], 0.0);
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
                     /* never leads to division by zero as negative coefficients can never be equal to the cliquemaximum */
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

#ifdef SCIP_DEBUG_BOUNDS
      SCIPdebugMsg(scip, "Updating bounds:\n");
#endif
      /* update bound arrays and determine success */
      for( i = 0; i < nvars; i++ )
      {
         idx = varinds[i];
         var = SCIPmatrixGetVar(matrix, idx);

         SCIPdebugMsg(scip, "%g <= %g <= %s <= %g <= %g\n", lbs[idx], newlbs[i], SCIPvarGetName(var), newubs[i], ubs[idx]);
         assert(SCIPisLE(scip, lbs[idx], newlbs[i]));
         assert(SCIPisGE(scip, ubs[idx], newubs[i]));

#ifdef SCIP_DEBUG_SUBSCIP
         /* TODO Adjust the subscip solve to consider the clique extension
          * assert(SCIPisEQ(scip, newlbs[i], subsciplbs[i]));
          * assert(SCIPisEQ(scip, newubs[i], subscipubs[i])); */
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
   SCIPfreeBufferArray(scip, &gradients);
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

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyTworowbnd)
{
   SCIP_PRESOLDATA* presoldata;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver if copying is enabled */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   if( presoldata->enablecopy )
   {
      SCIP_CALL( SCIPincludePresolTworowbnd(scip) );
   }

   return SCIP_OKAY;
}

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

   int i;
   int j;
   int k;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;
   infeasible = FALSE;

   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   if( presoldata->nuselessruns >= 5 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, FALSE, &initialized, &complete) );

   if( !initialized )
      return SCIP_OKAY;

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   if( nrows <= 1 )
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

   pospp = 0;
   posmm = 0;
   pospm = 0;
   posmp = 0;
   listsizepp = nrows;
   listsizemm = nrows;
   listsizepm = nrows;
   listsizemp = nrows;
   maxhashes = presoldata->maxhashfac == -1 ? SCIP_LONGINT_MAX : (((SCIP_Longint)nrows) * presoldata->maxhashfac);

   /* skim through the problem and create hashlists for combination candidates */
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
                  /* right-shift is required because we want to sort the hashes later on */
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
                  /* right-shift is required because we want to sort the hashes later on */
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

   /* Process pp and mm hashlists */
   if( pospp > 0 && posmm > 0 )
   {
      SCIP_Longint maxcombines;
      int ncombines;
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
      maxcombines = presoldata->maxpairfac == -1 ? SCIP_LONGINT_MAX : (((SCIP_Longint)nrows) * presoldata->maxpairfac);

      ncombines = 0;
      combinefails = 0;
      retrievefails = 0;
      findNextBlock(hashlistpp, pospp, &block1start, &block1end);
      findNextBlock(hashlistmm, posmm, &block2start, &block2end);
      SCIPdebugMsg(scip, "processing pp and mm\n");
      while( !finished )
      {
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
                     if( !SCIPhashsetExists(pairhashset, encodeRowPair(&rowpair)) )
                     {
                        assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row1idx)));
                        assert(!SCIPisInfinity(scip, -SCIPmatrixGetRowLhs(matrix, rowpair.row2idx)));

                        success = FALSE;

                        /* apply lp-based bound tightening */
                        if( presoldata->lpbound )
                        {
                           swaprow1 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx));
                           swaprow2 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx));

                           SCIP_CALL( applyLPboundTightening(scip, matrix, rowpair.row1idx, rowpair.row2idx, swaprow1, swaprow2, newlbs, newubs, delcons, &success) );
                        }

                        /* apply bound tightening on convex-combined rows */
                        if( presoldata->convcomb )
                        {
                           SCIP_CALL( combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                              FALSE, FALSE, newlbs, newubs, &success) );
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx)) )
                              SCIP_CALL( combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                 TRUE, FALSE, newlbs, newubs, &success) );
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx)) )
                              SCIP_CALL( combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                 FALSE, TRUE, newlbs, newubs, &success) );
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx))
                                 && !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx)) )
                              SCIP_CALL( combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                 TRUE, TRUE, newlbs, newubs, &success) );
                        }
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

   /* Process pm and mp hashlists */
   if( pospm > 0 && posmp > 0 )
   {
      SCIP_Longint maxcombines;
      int ncombines;
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
      maxcombines = presoldata->maxpairfac == -1 ? SCIP_LONGINT_MAX : (((SCIP_Longint)nrows) * presoldata->maxpairfac);

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

                        /* apply lp-based bound tightening */
                        if( presoldata->lpbound )
                        {
                           swaprow1 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx));
                           swaprow2 = !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx));

                           SCIP_CALL( applyLPboundTightening(scip, matrix, rowpair.row1idx, rowpair.row2idx, swaprow1, swaprow2, newlbs, newubs, delcons, &success) );
                        }

                        /* apply bound tightening on convex-combined rows */
                        if( presoldata->convcomb )
                        {
                           SCIP_CALL( combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                              FALSE, FALSE, newlbs, newubs, &success) );
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx)) )
                              SCIP_CALL( combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                 TRUE, FALSE, newlbs, newubs, &success) );
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx)) )
                              SCIP_CALL( combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                 FALSE, TRUE, newlbs, newubs, &success) );
                           if( !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row1idx))
                                 && !SCIPisInfinity(scip, SCIPmatrixGetRowRhs(matrix, rowpair.row2idx)) )
                              SCIP_CALL( combineRows(scip, matrix, rowpair.row1idx, rowpair.row2idx,
                                 TRUE, TRUE, newlbs, newubs, &success) );
                        }
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
            SCIPdebugMessage(" -> infeasible lower bound tightening: %s >= %g\n", SCIPvarGetName(var), newlbs[i]);
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
            SCIPdebugMessage(" -> infeasible upper bound tightening: %s <= %g\n", SCIPvarGetName(var), newubs[i]);
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
         delcons[i] = FALSE; /* Remove nonzero from clean buffer array */
         SCIPdebugMsg(scip, "removing redundant constraint %s\n", SCIPmatrixGetRowName(matrix, i));
         SCIP_CALL( SCIPdelCons(scip, SCIPmatrixGetCons(matrix, i)) );
      }
      else
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
         }
      }
   }
   (*ndelconss) += ndelcons;

   /* set result */
   if( *nchgbds > oldnchgbds || *nfixedvars > oldnfixedvars || ndelcons > 0 )
   {
      *result = SCIP_SUCCESS;
      presoldata->nuselessruns = 0;
   }
   else if( infeasible )
   {
      *result = SCIP_CUTOFF;
   }
   else
   {
      presoldata->nuselessruns++;
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

   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyTworowbnd) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeTworowbnd) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitTworowbnd) );

   /* add tworowbnd presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/tworowbnd/enablecopy",
         "should tworowbnd presolver be copied to sub-SCIPs?",
         &presoldata->enablecopy, TRUE, DEFAULT_ENABLECOPY, NULL, NULL) );
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
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/tworowbnd/lpbound",
         "Should the bounds be tightened via LP-bound?",
         &presoldata->lpbound, FALSE, DEFAULT_LPBOUND, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/tworowbnd/convcomb",
         "Should the bounds be tightened via feasibility based bound-tightening on convex combined constraints?",
         &presoldata->convcomb, FALSE, DEFAULT_CONVCOMB, NULL, NULL) );

   return SCIP_OKAY;
}
