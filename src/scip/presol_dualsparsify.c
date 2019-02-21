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

/**@file   presol_dualsparsify.c
 * @brief  cancel non-zeros of the constraint matrix based on the column
 * @author Dieter Weninger
 * @author Robert Lion Gottwald
 * @author Ambros Gleixner
 * @author Weikun Chen
 *
 * TODO: 1. allow to add constraints @done
 *       2. checkholes: @done
 *       3. use implied bounds information @done
 *       4. knapsack constraint (unknown error) @done
 *       5. update locked number control @done
 *       6. currently, it seems impossible to aggregated a binary variable in SCIP
 *          even after aggregation, the new variable is binary. This is because 
 *          the aggregated variable is marked as multi-aggregated variable but the 
 *          SCIPvarCompareActiveAndNegated function assert that this is wrong.
 *       7. update the fillins
 *
 * This presolver attempts to cancel non-zero entries of the constraint
 * matrix by adding scaled variables to other variables.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "scip/cons_linear.h"
#include "scip/presol_dualsparsify.h"
#include "scip/pub_cons.h"
#include "scip/pub_matrix.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_presol.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_presol.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"
#include <string.h>

#define PRESOL_NAME            "dualsparsify"
#define PRESOL_DESC            "eliminate non-zero coefficients"

#define PRESOL_PRIORITY           -240000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               -1    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_ENABLECOPY           TRUE    /**< should dualsparsify presolver be copied to sub-SCIPs? */
#define DEFAULT_CANCELLINEAR         TRUE    /**< should we cancel nonzeros in constraints of the linear constraint handler? */
#define DEFAULT_PRESERVEINTCOEFS    FALSE    /**< should we forbid cancellations that destroy integer coefficients? */
#define DEFAULT_MAX_CONT_FILLIN         1    /**< default value for the maximal fillin for continuous variables */
#define DEFAULT_MAX_BIN_FILLIN          1    /**< default value for the maximal fillin for binary variables */
#define DEFAULT_MAX_INT_FILLIN          1    /**< default value for the maximal fillin for integer variables (including binary) */
#define DEFAULT_MAXNONZEROS            -1    /**< maximal support of one equality to be used for cancelling (-1: no limit) */
#define DEFAULT_MAXCONSIDEREDNONZEROS  70    /**< maximal number of considered non-zeros within one row (-1: no limit) */
#define DEFAULT_ROWSORT               'd'    /**< order in which to process inequalities ('n'o sorting, 'i'ncreasing nonzeros, 'd'ecreasing nonzeros) */
#define DEFAULT_MAXRETRIEVEFAC      100.0    /**< limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints */
#define DEFAULT_WAITINGFAC            2.0    /**< number of calls to wait until next execution as a multiple of the number of useless calls */
#define DEFAULT_ISLOCKEDLIMIT       FALSE    /**< should locked limited be used to stop some aggregation */
#define DEFAULT_MINCANCELNNZSADDCONS  100    /**< minimalcancelledd nonzeros within one column if we need to add a constraint to the problem */

#define MAXSCALE                   1000.0    /**< maximal allowed scale for cancelling non-zeros */
#define MINSCALE                    0.001    /**< minimal allowed scale for cancelling non-zeros */


#define DEFAULT_MINACCEPTCANCELNNZS    10    /**< minimal cancel nonzeros when accepting to aggregate the (nonfree) variable */
#define DEFAULT_MINCONSIDEREDNNZS      10    /**< minimal number of considered non-zeros within one column */
#define DEFAULT_MAXCOMPAREDEVERYPAIR  200    /**< maximal number on the implementaion of doing reduction on every pair of variables */
#define DEFAULT_MAX_FILLINRATE       0.02    /**< cancel the variable if nfillins/ncancels is less than or equal this rate */
#define DEFAULT_MAX_FILLINRATE_FREE  0.01    /**< cancel the free variable if nfillins/ncancels is less than or equal this rate */

/*
 * Data structures
 */


/** presolver data */
struct SCIP_PresolData
{
   int                   ncancels;           /**< total number of canceled nonzeros (net value, i.e., removed minus added nonzeros) */
   int                   nfillin;            /**< total number of added nonzeros */
   int                   nfailures;          /**< number of calls to presolver without success */
   int                   nwaitingcalls;      /**< number of presolver calls until next real execution */
   int                   maxcontfillin;      /**< maximal fillin for continuous variables */
   int                   maxintfillin;       /**< maximal fillin for integer variables*/
   int                   maxbinfillin;       /**< maximal fillin for binary variables */
   int                   maxnonzeros;        /**< maximal support of one equality to be used for cancelling (-1: no limit) */
   int                   maxconsiderednonzeros;/**< maximal number of considered non-zeros within one row (-1: no limit) */
   SCIP_Real             maxretrievefac;     /**< limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints */
   SCIP_Real             waitingfac;         /**< number of calls to wait until next execution as a multiple of the number of useless calls */
   char                  rowsort;            /**< order in which to process inequalities ('n'o sorting, 'i'ncreasing nonzeros, 'd'ecreasing nonzeros) */
   SCIP_Bool             enablecopy;         /**< should dualsparsify presolver be copied to sub-SCIPs? */
   SCIP_Bool             cancellinear;       /**< should we cancel nonzeros in constraints of the linear constraint handler? */
   SCIP_Bool             preserveintcoefs;   /**< should we forbid cancellations that destroy integer coefficients? */
   SCIP_Bool             islockedlimit;      /**< should locked limit be used to store some aggregation */
   int                   naggregated;
   int                   mincancelnnzsaddcons;/**< minimal cancelled nonzeros within one column if we need to add a constraint to the problem */

   int                   minacceptcancelnnzs;/**< minimal cancel nonzeros when accepting to aggregate the (nonfree) variable */
   int                   minconsiderednnzs;  /**< minimal number of considered non-zeros within one column */
   int                   maxcompareeverypair;/**< maximal number on the implementaion of doing reduction on every pair of variables */
   SCIP_Real             maxfillinrate;      /**< cancel the variable if nfillins/ncancels is less than or equal this rate */
   SCIP_Real             maxfillinratefree;  /**< cancel the free variable if nfillins/ncancels is less than or equal this rate */
};

/** structure representing a pair of constraints in a cols; used for lookup in a hashtable */
struct ColConsPair
{
   int colindex;
   int consindex1;
   int consindex2;
   SCIP_Real conscoef1;
   SCIP_Real conscoef2;
};

typedef struct ColConsPair COLCONSPAIR;

/*
 * Local methods
 */

/** returns TRUE iff both keys are equal */
static
SCIP_DECL_HASHKEYEQ(consPairsEqual)
{  /*lint --e{715}*/
   SCIP* scip;
   COLCONSPAIR* conspair1;
   COLCONSPAIR* conspair2;
   SCIP_Real ratio1;
   SCIP_Real ratio2;

   scip = (SCIP*) userptr;
   conspair1 = (COLCONSPAIR*) key1;
   conspair2 = (COLCONSPAIR*) key2;

   if( conspair1->consindex1 != conspair2->consindex1 )
      return FALSE;

   if( conspair1->consindex2 != conspair2->consindex2 )
      return FALSE;

   ratio1 = conspair1->conscoef2 / conspair1->conscoef1;
   ratio2 = conspair2->conscoef2 / conspair2->conscoef1;

   return SCIPisEQ(scip, ratio1, ratio2);
}

/** returns the hash value of the key */
static
SCIP_DECL_HASHKEYVAL(consPairHashval)
{  /*lint --e{715}*/
   COLCONSPAIR* conspair;

   conspair = (COLCONSPAIR*) key;

   return SCIPhashTwo(SCIPcombineTwoInt(conspair->consindex1, conspair->consindex2),
                      SCIPrealHashCode(conspair->conscoef2 / conspair->conscoef1));
}

/** calculate max activity of one row without one column */
static
SCIP_Real getMaxActivitySingleRowWithoutCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col                 /**< column index */
   )
{
   int c;
   SCIP_Real val;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   SCIP_Real maxactivity;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(matrix != NULL);

   maxactivity = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      c = *rowpnt;
      val = *valpnt;

      if( c == col )
         continue;

      if( val > 0.0 )
      {
         ub = SCIPmatrixGetColUb(matrix, c);
         assert(!SCIPisInfinity(scip, ub));
         maxactivity += val * ub;
      }
      else if( val < 0.0 )
      {
         lb = SCIPmatrixGetColLb(matrix, c);
         assert(!SCIPisInfinity(scip, -lb));
         maxactivity += val * lb;
      }
   }

   return maxactivity;
}

/** calculate min activity of one row without one column */
static
SCIP_Real getMinActivitySingleRowWithoutCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col                 /**< column index */
   )
{
   int c;
   SCIP_Real val;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   SCIP_Real minactivity;
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(matrix != NULL);

   minactivity = 0;

   rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
   rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
   valpnt = SCIPmatrixGetRowValPtr(matrix, row);

   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      c = *rowpnt;
      val = *valpnt;

      if( c == col )
         continue;

      if( val > 0.0 )
      {
         lb = SCIPmatrixGetColLb(matrix, c);
         assert(!SCIPisInfinity(scip, -lb));
         minactivity += val * lb;
      }
      else if( val < 0.0 )
      {
         ub = SCIPmatrixGetColUb(matrix, c);
         assert(!SCIPisInfinity(scip, ub));
         minactivity += val * ub;
      }
   }

   return minactivity;
}

/** get min/max residual activity without the specified column */
static
void getMinMaxActivityResiduals(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this variable in this row */
   SCIP_Real*            minresactivity,     /**< minimum residual activity of this row */
   SCIP_Real*            maxresactivity,     /**< maximum residual activity of this row */
   SCIP_Bool*            isminsettoinfinity, /**< flag indicating if minresactiviy is set to infinity */
   SCIP_Bool*            ismaxsettoinfinity  /**< flag indicating if maxresactiviy is set to infinity */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   int nmaxactneginf;
   int nmaxactposinf;
   int nminactneginf;
   int nminactposinf;
   SCIP_Real maxactivity;
   SCIP_Real minactivity;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);
   assert(isminsettoinfinity != NULL);
   assert(ismaxsettoinfinity != NULL);

   lb = SCIPmatrixGetColLb(matrix, col);
   ub = SCIPmatrixGetColUb(matrix, col);

   *isminsettoinfinity = FALSE;
   *ismaxsettoinfinity = FALSE;

   nmaxactneginf = SCIPmatrixGetRowNMaxActNegInf(matrix, row);
   nmaxactposinf = SCIPmatrixGetRowNMaxActPosInf(matrix, row);
   nminactneginf = SCIPmatrixGetRowNMinActNegInf(matrix, row);
   nminactposinf = SCIPmatrixGetRowNMinActPosInf(matrix, row);

   maxactivity = SCIPmatrixGetRowMaxActivity(matrix, row);
   minactivity = SCIPmatrixGetRowMinActivity(matrix, row);

   if( val >= 0.0 )
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(nmaxactposinf >= 1);
         if( nmaxactposinf == 1 && nmaxactneginf == 0 )
            *maxresactivity = getMaxActivitySingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nmaxactneginf + nmaxactposinf) > 0 )
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
         else
            *maxresactivity = maxactivity - val * ub;
      }

      if( SCIPisInfinity(scip, -lb) )
      {
         assert(nminactneginf >= 1);
         if( nminactneginf == 1 && nminactposinf == 0 )
            *minresactivity = getMinActivitySingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nminactneginf + nminactposinf) > 0 )
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
         else
            *minresactivity = minactivity - val * lb;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(nmaxactneginf >= 1);
         if( nmaxactneginf == 1 && nmaxactposinf == 0 )
            *maxresactivity = getMaxActivitySingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nmaxactneginf + nmaxactposinf) > 0 )
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
         else
            *maxresactivity = maxactivity - val * lb;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         assert(nminactposinf >= 1);
         if( nminactposinf == 1 && nminactneginf == 0 )
            *minresactivity = getMinActivitySingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (nminactneginf + nminactposinf) > 0 )
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
         else
            *minresactivity = minactivity - val * ub;
      }
   }
}

/** calculate the upper bound of one variable from one row */
static
void getVarUpperBoundOfRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index of variable */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this column in this row */
   SCIP_Real*            rowub,              /**< upper bound of row */
   SCIP_Bool*            ubfound             /**< flag indicating that an upper bound was calculated */
   )
{
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(rowub != NULL);
   assert(ubfound != NULL);

   *rowub = SCIPinfinity(scip);
   *ubfound = FALSE;

   getMinMaxActivityResiduals(scip, matrix, col, row, val, &minresactivity, &maxresactivity,
         &isminsettoinfinity, &ismaxsettoinfinity);

   lhs = SCIPmatrixGetRowLhs(matrix, row);
   rhs = SCIPmatrixGetRowRhs(matrix, row);

   if( val > 0.0 )
   {
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) )
      {
         *rowub = (rhs - minresactivity)/val;
         *ubfound = TRUE;
      }
   }
   else
   {
      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) )
      {
         *rowub = (lhs - maxresactivity)/val;
         *ubfound = TRUE;
      }
   }
}


/** calculate the lower bound of one variable from one row */
static
void getVarLowerBoundOfRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index of variable */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this column in this row */
   SCIP_Real*            rowlb,              /**< lower bound of row */
   SCIP_Bool*            lbfound             /**< flag indicating that an lower bound was calculated */
   )
{
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(rowlb != NULL);
   assert(lbfound != NULL);

   *rowlb = -SCIPinfinity(scip);
   *lbfound = FALSE;

   getMinMaxActivityResiduals(scip, matrix, col, row, val, &minresactivity, &maxresactivity,
         &isminsettoinfinity, &ismaxsettoinfinity);

   lhs = SCIPmatrixGetRowLhs(matrix, row);
   rhs = SCIPmatrixGetRowRhs(matrix, row);

   if( val < 0.0 )
   {
      if( !isminsettoinfinity && !SCIPisInfinity(scip, rhs) )
      {
         *rowlb = (rhs - minresactivity)/val;
         *lbfound = TRUE;
      }
   }
   else
   {
      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -lhs) )
      {
         *rowlb = (lhs - maxresactivity)/val;
         *lbfound = TRUE;
      }
   }
}


/** verify whether variable upper bound is implied */
static
SCIP_Bool isUpperBoundImplied(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col                 /**< column index for implied free test */
   )
{
   SCIP_Real varub;
   SCIP_Real impliedub;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;
   SCIP_Bool ubimplied;

   assert(scip != NULL);
   assert(matrix != NULL);

   varub = SCIPmatrixGetColUb(matrix, col);
   if( SCIPisInfinity(scip, varub) )
      return TRUE;

   ubimplied = FALSE;
   impliedub = SCIPinfinity(scip);

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for( ; (colpnt < colend); colpnt++, valpnt++ )
   {
      SCIP_Real rowub;
      SCIP_Bool ubfound;

      getVarUpperBoundOfRow(scip, matrix, col, *colpnt, *valpnt, &rowub, &ubfound);

      if( ubfound && (rowub < impliedub) )
         impliedub = rowub;
   }

   if( SCIPisFeasLE(scip, impliedub, varub) )
      ubimplied = TRUE;

   return ubimplied;
}

/** verify whether variable lower bound is implied */
static
SCIP_Bool isLowerBoundImplied(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_MATRIX*          matrix,             /**< matrix containing the constraints */
   int                   col                 /**< column index for implied free test */
   )
{
   SCIP_Real varlb;
   SCIP_Real impliedlb;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;
   SCIP_Bool lbimplied;

   assert(scip != NULL);
   assert(matrix != NULL);

   varlb = SCIPmatrixGetColLb(matrix, col);
   if( SCIPisInfinity(scip, -varlb) )
      return TRUE;

   lbimplied = FALSE;
   impliedlb = -SCIPinfinity(scip);

   colpnt = SCIPmatrixGetColIdxPtr(matrix, col);
   colend = colpnt + SCIPmatrixGetColNNonzs(matrix, col);
   valpnt = SCIPmatrixGetColValPtr(matrix, col);

   for( ; (colpnt < colend); colpnt++, valpnt++ )
   {
      SCIP_Real rowlb;
      SCIP_Bool lbfound;

      getVarLowerBoundOfRow(scip, matrix, col, *colpnt, *valpnt, &rowlb, &lbfound);

      if( lbfound && (rowlb > impliedlb) )
         impliedlb = rowlb;
   }

   if( SCIPisFeasGE(scip, impliedlb, varlb) )
      lbimplied = TRUE;

   return lbimplied;
}


/* y = weght1*var[colidx1] + var[colidx2] */
static
SCIP_RETCODE aggregation(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_MATRIX*          matrix,             /**< the constraint matrix */
   SCIP_VAR**            vars,               /**< the current variables*/
   int                   colidx1,            /**< one of the indexes of column to try non-zero cancellation for */
   int                   colidx2,            /**< one of the indexes of column to try non-zero cancellation for */
   SCIP_Bool             isimpliedfree,      /**< is the aggregated variable implied free? */
   SCIP_Real             weight1            /**< weight variable one in the aggregated expression */
      )
{
   SCIP_VAR* tmpvars[2];
   SCIP_Real coefs[2];
   SCIP_VAR* aggregatedvar;
   SCIP_VAR* newvar;
   SCIP_CONS* newcons;
   SCIP_Real newlb;
   SCIP_Real newub;
   char newvarname[SCIP_MAXSTRLEN];
   char newconsname[SCIP_MAXSTRLEN];
   SCIP_Bool infeasible;
   SCIP_Bool aggregated;
   SCIP_VARTYPE newvartype;
   SCIP_Real constant;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert( !SCIPisZero(scip, weight1) );

   presoldata->naggregated += 1;
   aggregatedvar = vars[colidx2];

   (void) SCIPsnprintf(newvarname, SCIP_MAXSTRLEN, "dualsparsifyvar_%d", presoldata->naggregated);

   constant = 0.0;

   if( SCIPvarIsIntegral(vars[colidx1]) && SCIPvarIsIntegral(vars[colidx2]) )
      printf("aggregated integral variables\n");

   if( weight1 > 0 )
   {
      if(SCIPisInfinity(scip, -SCIPvarGetLbGlobal(vars[colidx1])) || SCIPisInfinity(scip, -SCIPvarGetLbGlobal(vars[colidx2])))
         newlb = -SCIPinfinity(scip);
      else
         newlb = weight1*SCIPvarGetLbGlobal(vars[colidx1]) + SCIPvarGetLbGlobal(vars[colidx2]);

      if(SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[colidx1])) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[colidx2])))
         newub = SCIPinfinity(scip);
      else
         newub = weight1*SCIPvarGetUbGlobal(vars[colidx1]) + SCIPvarGetUbGlobal(vars[colidx2]);
   }

   else
   {
      if(SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[colidx1])) || SCIPisInfinity(scip, -SCIPvarGetLbGlobal(vars[colidx2])))
         newlb = -SCIPinfinity(scip);
      else
         newlb = weight1*SCIPvarGetUbGlobal(vars[colidx1]) + SCIPvarGetLbGlobal(vars[colidx2]);

      if(SCIPisInfinity(scip, SCIPvarGetLbGlobal(vars[colidx1])) || SCIPisInfinity(scip, SCIPvarGetUbGlobal(vars[colidx2])))
         newub = SCIPinfinity(scip);
      else
         newub = weight1*SCIPvarGetLbGlobal(vars[colidx1]) + SCIPvarGetUbGlobal(vars[colidx2]);
   }
   if( SCIPvarIsIntegral(aggregatedvar) )
   {
      newvartype = (SCIPvarGetType(aggregatedvar) == SCIP_VARTYPE_IMPLINT) ? SCIP_VARTYPE_IMPLINT : SCIP_VARTYPE_INTEGER;
   }
   else
   {
      newvartype = SCIP_VARTYPE_CONTINUOUS;
   }

   lhs = SCIPvarGetLbGlobal(vars[colidx2]);
   rhs = SCIPvarGetUbGlobal(vars[colidx2]);


   SCIP_CALL( SCIPcreateVar(scip, &newvar, newvarname, newlb, newub, 0.0, newvartype,
            SCIPvarIsInitial(aggregatedvar), SCIPvarIsRemovable(aggregatedvar), NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, newvar) );
//   printf("%d, %d, %s, %s, %d, %8.4f\n", SCIPvarGetType(vars[colidx1]), SCIPvarGetType(vars[colidx2]), SCIPvarGetName(vars[colidx1]), SCIPvarGetName(vars[colidx2])
//        ,newvartype, constant );

   tmpvars[0] = vars[colidx1];
   tmpvars[1] = newvar;
   coefs[0] = -weight1;
   coefs[1] = 1;

   SCIP_CALL( SCIPmultiaggregateVar(scip, aggregatedvar, 2, tmpvars, coefs, constant, &infeasible, &aggregated) );
   assert(!infeasible);
   assert(aggregated);

   vars[colidx2] = newvar;

   if( !isimpliedfree )
   {
      (void) SCIPsnprintf(newconsname, SCIP_MAXSTRLEN, "dualsparsifycons_%d", presoldata->naggregated);

      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, newconsname, 2, tmpvars, coefs,
               lhs, rhs, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIPdebugPrintCons(scip, newcons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   }

   SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
}

/** try non-zero cancellation for given row */
static
SCIP_RETCODE cancelCol(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_MATRIX*          matrix,             /**< the constraint matrix */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_HASHTABLE*       pairtable,          /**< the hashtable containing ROWVARPAIR's of equations */
   SCIP_Bool*            isimpliedfrees,     /**< array to indicates whether it is impliedfree or not */
   SCIP_VAR**            vars,               /**< array to store the current variables */
   SCIP_Bool*            isblockedvar,       /**< array to indicates whether it is blocked or not */
   int                   colidx,             /**< index of row to try non-zero cancellation for */
   int                   maxcontfillin,      /**< maximal fill-in allowed for continuous variables */
   int                   maxintfillin,       /**< maximal fill-in allowed for integral variables */
   int                   maxbinfillin,       /**< maximal fill-in allowed for binary variables */
   int                   maxconsiderednonzeros, /**< maximal number of non-zeros to consider for cancellation */
   SCIP_Bool             preserveintcoefs,   /**< only perform non-zero cancellation if integrality of coefficients is preserved? */
   SCIP_Longint*         nuseless,           /**< pointer to update number of useless hashtable retrieves */
   int*                  nchgcoefs,          /**< pointer to update number of changed coefficients */
   int*                  ncanceled,          /**< pointer to update number of canceled nonzeros */
   int*                  nfillin,            /**< pointer to update the produced fill-in */
   SCIP_Bool             isaddedcons         /**< whether a linear constraint required to added to keep the validity */
   )
{
   int* cancelcolinds;
   SCIP_Real* cancelcolvals;
   SCIP_Real cancellb;
   SCIP_Real cancelub;
   int bestcand;
   int bestnfillin;
   SCIP_Real bestscale;
   SCIP_Real bestcancelrate;
   int* tmpinds;
   SCIP_Real* scores;
   SCIP_Real* tmpvals;
   int cancelcollen;
   int* colidxptr;
   SCIP_Real* colvalptr;
   int nchgcoef;
   int nretrieves;
   SCIP_Real mincancelrate;
   SCIP_Bool colisimpl;
   SCIP_VAR* cancelvar;
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   SCIP_Bool cancelisbin;
   SCIP_Real ncols;

   colisimpl = isimpliedfrees[colidx];
   ncols = SCIPmatrixGetNColumns(matrix);

   cancelcollen = SCIPmatrixGetColNNonzs(matrix, colidx);
   colidxptr = SCIPmatrixGetColIdxPtr(matrix, colidx);
   colvalptr = SCIPmatrixGetColValPtr(matrix, colidx);

   cancelvar = vars[colidx];

   mincancelrate = 0.0;

   SCIP_CALL( SCIPduplicateBufferArray(scip, &cancelcolinds, colidxptr, cancelcollen) );
   SCIP_CALL( SCIPduplicateBufferArray(scip, &cancelcolvals, colvalptr, cancelcollen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpinds, cancelcollen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvals, cancelcollen) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scores, cancelcollen) );

   cancellb = SCIPvarGetLbGlobal(cancelvar);
   cancelub = SCIPvarGetUbGlobal(cancelvar);
   cancelisbin = ((SCIPvarGetType(cancelvar) == SCIP_VARTYPE_BINARY) ||
                     (SCIPvarIsIntegral(cancelvar) && SCIPisZero(scip, cancellb) && SCIPisEQ(scip, cancelub, 1.0)));

   nchgcoef = 0;
   nretrieves = 0;
   while( TRUE ) /*lint !e716 */
   {
      int i;
      int j;
      COLCONSPAIR colconspair;
      int maxlen;

      bestcand = -1;
      bestnfillin = 0;
      bestscale = 1.0;
      bestcancelrate = 0.0;

      for( i = 0; i < cancelcollen; ++i )
      {
         tmpinds[i] = i;
         scores[i] = -SCIPmatrixGetRowNNonzs(matrix, cancelcolinds[i]) - 1.0*cancelcolinds[i]/(ncols);
      }

      SCIPsortRealInt(scores, tmpinds, cancelcollen);

      maxlen = cancelcollen;
      if( maxconsiderednonzeros >= 0 )
         maxlen = MIN(cancelcollen, maxconsiderednonzeros);

      for( i = 0; i < maxlen; ++i )
      {
         for( j = i + 1; j < maxlen; ++j )
         {
            int a,b;
            int ncancel;
            int ncontfillin;
            int nintfillin;
            int nbinfillin;
            int ntotfillin;
            int implcollen;
            COLCONSPAIR* implcolconspair;
            SCIP_Real* implcolvals;
            int* implcolinds;
            SCIP_VAR* implcolvar;
            SCIP_Real implcollb;
            SCIP_Real implcolub;
            SCIP_Bool implcolisbin;
            SCIP_Real scale;
            SCIP_Real cancelrate;
            int i1,i2;
            SCIP_Bool abortpair;
            SCIP_Real rowlhs;
            SCIP_Real rowrhs;

            i1 = tmpinds[i];
            i2 = tmpinds[j];
         
            
            assert(cancelcolinds[i] < cancelcolinds[j]);

            if( cancelcolinds[i1] < cancelcolinds[i2] )
            {
               colconspair.consindex1 = cancelcolinds[i1];
               colconspair.consindex2 = cancelcolinds[i2];
               colconspair.conscoef1 = cancelcolvals[i1];
               colconspair.conscoef2 = cancelcolvals[i2];
            }
            else
            {
               colconspair.consindex1 = cancelcolinds[i2];
               colconspair.consindex2 = cancelcolinds[i1];
               colconspair.conscoef1 = cancelcolvals[i2];
               colconspair.conscoef2 = cancelcolvals[i1];
            }

            implcolconspair = (COLCONSPAIR*)SCIPhashtableRetrieve(pairtable, (void*) &colconspair);
            nretrieves++;


            if( implcolconspair == NULL || implcolconspair->colindex == colidx || isblockedvar[implcolconspair->colindex] )
               continue;

            /* if the column we want to cancel is an implied free column, we will only use implied free columns
             * for canceling with less non-zeros and if the number of non-zeros is equal we use the
             * colindex as tie-breaker to avoid cyclic non-zero cancellation
             */
            implcollen = SCIPmatrixGetColNNonzs(matrix, implcolconspair->colindex);
            if( colisimpl && (cancelcollen < implcollen || (cancelcollen == implcollen && colidx < implcolconspair->colindex)) )
               continue;

            implcolvals = SCIPmatrixGetColValPtr(matrix, implcolconspair->colindex);
            implcolinds = SCIPmatrixGetColIdxPtr(matrix, implcolconspair->colindex);
            implcolvar = vars[implcolconspair->colindex];
            implcollb = SCIPvarGetLbGlobal(implcolvar);
            implcolub = SCIPvarGetUbGlobal(implcolvar);
            implcolisbin = ((SCIPvarGetType(implcolvar) == SCIP_VARTYPE_BINARY) || 
                              (SCIPvarIsIntegral(implcolvar) && SCIPisZero(scip, implcollb) && SCIPisEQ(scip, implcolub, 1.0)));
 
            scale = -colconspair.conscoef1 / implcolconspair->conscoef1;

            if( REALABS(scale) > MAXSCALE )
               continue;

            /* @todo maybe do more reduction if knspsack constraint handler supports downgrading constraint,
             * i.e., converting into a linear constraint
             */
            if( implcolisbin )
               continue;
            else if( SCIPvarIsIntegral(implcolvar) )
            {
               if( SCIPvarIsIntegral(cancelvar) )
               {
                  if( (SCIPvarGetType(implcolvar) != SCIP_VARTYPE_IMPLINT) && (SCIPvarGetType(cancelvar) == SCIP_VARTYPE_IMPLINT) )
                     continue;
                  //if( !SCIPisIntegral(scip, scale) || isHoleExist(scip, vars, implcolconspair->colindex, colidx, -scale) )
                  if( !SCIPisIntegral(scip, scale) )
                     continue;
               }
               else
                  continue;
            }

            a = 0;
            b = 0;
            ncancel = 0;

            ncontfillin = 0;
            nintfillin = 0;
            nbinfillin = 0;
            abortpair = FALSE;

            while( a < cancelcollen && b < implcollen )
            {
               if( cancelcolinds[a] == implcolinds[b] )
               {
                  SCIP_Real newcoef;

                  newcoef = cancelcolvals[a] + scale * implcolvals[b];

                  /* check if coefficient is cancelled */
                  if( SCIPisZero(scip, newcoef) )
                  {
                     ++ncancel;
                  }
                  /* otherwise, check if integral coefficients are preserved if the column is integral */
                  else if( (preserveintcoefs && SCIPvarIsIntegral(cancelvar) &&
                            SCIPisIntegral(scip, cancelcolvals[a]) && !SCIPisIntegral(scip, newcoef)) )
                  { 
                     abortpair = TRUE;
                     break;
                  }
                  /* finally, check if locks could be modified in a bad way due to flipped signs */
                  else if( COPYSIGN(1.0, newcoef) != COPYSIGN(1.0, cancelcolvals[a]) ) /*lint !e777*/
                  {
                     /* do not flip signs for non-canceled coefficients if this adds a lock to a variable that had at most one lock
                      * in that direction before, except if the other direction gets unlocked
                      */
                     rowrhs = SCIPmatrixGetRowRhs(matrix, cancelcolinds[a]);
                     rowlhs = SCIPmatrixGetRowLhs(matrix, cancelcolinds[a]);
                     if( (cancelcolvals[a] > 0.0 && ! SCIPisInfinity(scip, rowrhs)) ||
                         (cancelcolvals[a] < 0.0 && ! SCIPisInfinity(scip, -rowlhs)) )
                     {
                        /* if we get into this case the variable had a positive coefficient in a <= constraint or a negative
                         * coefficient in a >= constraint, e.g. an uplock. If this was the only uplock we do not abort their
                         * cancelling, otherwise we abort if we had a single or no downlock and add one
                         */
                        //TODO: we may change this 
                        if( presoldata->islockedlimit && (SCIPmatrixGetColNUplocks(matrix, colidx) > 1 &&
                            SCIPmatrixGetColNDownlocks(matrix, colidx) <= 1) )
                        {
                           abortpair = TRUE;
                           break;
                        }
                     }

                     if( (cancelcolvals[a] < 0.0 && ! SCIPisInfinity(scip, rowrhs)) ||
                         (cancelcolvals[a] > 0.0 && ! SCIPisInfinity(scip, -rowlhs)) )
                     {
                        /* symmetric case where the variable had a downlock */
                        if( presoldata->islockedlimit && (SCIPmatrixGetColNDownlocks(matrix, colidx) > 1 &&
                            SCIPmatrixGetColNUplocks(matrix, colidx) <= 1) )
                        {
                           abortpair = TRUE;
                           break;
                        }
                     }
                  }

                  ++a;
                  ++b;
               }
               else if( cancelcolinds[a] < implcolinds[b] )
               {
                  ++a;
               }
               else
               {
                  SCIP_Real newcoef;

                  //TODO: be more careful about this
                  rowrhs = SCIPmatrixGetRowRhs(matrix, implcolinds[b]);
                  rowlhs = SCIPmatrixGetRowLhs(matrix, implcolinds[b]);

                  newcoef = scale * implcolvals[b];

                  if( (newcoef > 0.0 && ! SCIPisInfinity(scip, rowrhs)) ||
                      (newcoef < 0.0 && ! SCIPisInfinity(scip, -rowlhs)) )
                  {
                     if( presoldata->islockedlimit && SCIPmatrixGetColNUplocks(matrix, colidx) <= 1 )
                     {
                        abortpair = TRUE;
                        ++b;
                        break;
                     }
                  }

                  if( (newcoef < 0.0 && ! SCIPisInfinity(scip, rowrhs)) ||
                      (newcoef > 0.0 && ! SCIPisInfinity(scip, -rowlhs)) )
                  {
                     if( presoldata->islockedlimit && SCIPmatrixGetColNDownlocks(matrix, colidx) <= 1 )
                     {
                        abortpair = TRUE;
                        ++b;
                        break;
                     }
                  }

                  ++b;

                  if( SCIPvarIsIntegral(cancelvar) )
                  {
                     if( SCIPvarIsBinary(cancelvar) && ++nbinfillin > maxbinfillin )
                     {
                        abortpair = TRUE;
                        break;
                     }

                     if( ++nintfillin > maxintfillin )
                     {
                        abortpair = TRUE;
                        break;
                     }
                  }
                  else
                  {
                     if( ++ncontfillin > maxcontfillin )
                     {
                        abortpair = TRUE;
                        break;
                     }
                  }
                  if( abortpair )
                  {
                     //printf( "fail due to too many fillins\n" );
                  }
               }    
            }

            if( abortpair )
               continue;

            //TODO: update the fillin here
            while( b < implcollen )
            {
               ++b;
               conshdlr = SCIPconsGetHdlr(SCIPmatrixGetCons(matrix, implcolinds[b]));
               conshdlrname = SCIPconshdlrGetName(conshdlr);
               if( SCIPvarIsIntegral(cancelvar) )
               {
                  if( SCIPvarIsBinary(cancelvar) && ++nbinfillin > maxbinfillin )
                     break;
                  if( ++nintfillin > maxintfillin )
                     break;
               }
               else
               {
                  if( ++ncontfillin > maxcontfillin )
                     break;
               }
            }


            if( ncontfillin > maxcontfillin || nbinfillin > maxbinfillin || nintfillin > maxintfillin )
            {
               //printf( "fail due to too many fillins\n" );
               continue;
            }

            ntotfillin = nintfillin + ncontfillin;

            if( ntotfillin >= ncancel )
               continue;

            cancelrate = (ncancel - ntotfillin) / (SCIP_Real) implcollen;
            
            /* if a linear constraint is needed to keep the validity, we require a large nonzero cancellation */
            if( isaddedcons && (ncancel - ntotfillin < presoldata->mincancelnnzsaddcons) )
               continue;

         

            if( cancelrate < mincancelrate )
               continue;

            if( cancelrate > bestcancelrate )
            {
               bestnfillin = ntotfillin;
               bestcand = implcolconspair->colindex;
               bestscale = scale;
               bestcancelrate = cancelrate;

               /* stop looking if the current candidate does not create any fill-in or alter coefficients */
               if( cancelrate == 1.0 )
                  break;
            }

            /* we accept the best candidate immediately if it does not create any fill-in or alter coefficients */
            if( bestcand != -1 && bestcancelrate == 1.0 )
               break;
         }
      }

      if( bestcand != -1 )
      {
         int a;
         int b;
         SCIP_Real* implcolvals;
         int* implcolinds;
         int implcollen;
         SCIP_Real implcollb;
         SCIP_Real implcolub;
         SCIP_VAR* implcolvar;
         SCIP_Bool implcolisbin;
         int tmpcollen;
   
         implcolvar = vars[bestcand];
         implcolvals = SCIPmatrixGetColValPtr(matrix, bestcand);
         implcolinds = SCIPmatrixGetColIdxPtr(matrix, bestcand);
         implcollen = SCIPmatrixGetColNNonzs(matrix, bestcand);
         implcollb = SCIPvarGetLbGlobal(implcolvar);
         implcolub = SCIPvarGetUbGlobal(implcolvar);
         implcolisbin = ((SCIPvarGetType(implcolvar) == SCIP_VARTYPE_BINARY) || 
                           (SCIPvarIsIntegral(implcolvar) && SCIPisZero(scip, implcollb) && SCIPisEQ(scip, implcolub, 1.0)));

         a = 0;
         b = 0;
         tmpcollen = 0;


         while( a < cancelcollen && b < implcollen )
         {
            conshdlr = SCIPconsGetHdlr(SCIPmatrixGetCons(matrix, implcolinds[b]));
            conshdlrname = SCIPconshdlrGetName(conshdlr);

            if( cancelcolinds[a] == implcolinds[b] )
            {
               SCIP_Real val = cancelcolvals[a] + bestscale * implcolvals[b];

               if( !SCIPisZero(scip, val) )
               {
                  tmpinds[tmpcollen] = cancelcolinds[a];
                  tmpvals[tmpcollen] = val;
                  ++tmpcollen;
               }
               ++nchgcoef;

               ++a;
               ++b;
            }
            else if( cancelcolinds[a] < implcolinds[b] )
            {
               tmpinds[tmpcollen] = cancelcolinds[a];
               tmpvals[tmpcollen] = cancelcolvals[a];
               ++tmpcollen;
               ++a;
            }
            else
            {
               tmpinds[tmpcollen] = implcolinds[b];
               tmpvals[tmpcollen] = implcolvals[b] * bestscale;
               ++nchgcoef;
               ++tmpcollen;
               ++b;
            }
         }

         while( a < cancelcollen )
         {
            tmpinds[tmpcollen] = cancelcolinds[a];
            tmpvals[tmpcollen] = cancelcolvals[a];
            ++tmpcollen;
            ++a;
         }

         while( b < implcollen )
         {
            conshdlr = SCIPconsGetHdlr(SCIPmatrixGetCons(matrix, implcolinds[b]));
            conshdlrname = SCIPconshdlrGetName(conshdlr);
            tmpinds[tmpcollen] = implcolinds[b];
            tmpvals[tmpcollen] = implcolvals[b] * bestscale;
            ++nchgcoef;
            ++tmpcollen;
            ++b;
         }

         /* update fill-in counter */
         *nfillin += bestnfillin;

         /* swap the temporary arrays so that the cancelrowinds and cancelrowvals arrays, contain the new
          * changed row, and the tmpinds and tmpvals arrays can be overwritten in the next iteration
          */
         SCIPswapPointers((void**) &tmpinds, (void**) &cancelcolinds);
         SCIPswapPointers((void**) &tmpvals, (void**) &cancelcolvals);
         cancelcollen = tmpcollen;
         SCIP_CALL( aggregation(scip, presoldata, matrix, vars, colidx, bestcand, !isaddedcons, -bestscale) );
      }
      else
         break;
   }


   if( nchgcoef != 0 )
   {

      /* update counters */
      *nchgcoefs += nchgcoef;
      *ncanceled += SCIPmatrixGetColNNonzs(matrix, colidx) - cancelcollen;

      isblockedvar[colidx] = TRUE;

      /* if successful, decrease the useless hashtable retrieves counter; the rationale here is that we want to keep
       * going if, after many useless calls that almost exceeded the budget, we finally reach a useful section; but we
       * don't allow a negative build-up for the case that the useful section is all at the beginning and we just want
       * to quit quickly afterwards
       */
      *nuseless -= nretrieves;
      *nuseless = MAX(*nuseless, 0);
   }
   else
   {
      /* if not successful, increase useless hashtable retrieves counter */
      *nuseless += nretrieves;
   }

   SCIPfreeBufferArray(scip, &scores);
   SCIPfreeBufferArray(scip, &tmpvals);
   SCIPfreeBufferArray(scip, &tmpinds);
   SCIPfreeBufferArray(scip, &cancelcolvals);
   SCIPfreeBufferArray(scip, &cancelcolinds);

   return SCIP_OKAY;
}

/** updates failure counter after one execution */
static
void updateFailureStatistic(
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_Bool             success             /**< was this execution successful? */
   )
{
   assert(presoldata != NULL);

   if( success )
   {
      presoldata->nfailures = 0;
      presoldata->nwaitingcalls = 0;
   }
   else
   {
      presoldata->nfailures++;
      presoldata->nwaitingcalls = (int)(presoldata->waitingfac*(SCIP_Real)presoldata->nfailures);
   }
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDualsparsify)
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
      SCIP_CALL( SCIPincludePresolDualsparsify(scip) );
   }

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualsparsify)
{  /*lint --e{715}*/
   SCIP_MATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;
   int ncols;
   int nrows;
//   int r;
   int c;
   int i;
   int j;
   int ncancels;
   int nfillins;
   int nchgcoef;
   int ncancel;
   int* perm;
   SCIP_Real* scores;
   int* colidxsorted;
   int* colsparsity;
   COLCONSPAIR* conspairs;
   int nconspairs;
   int conspairssize;
   SCIP_HASHTABLE* pairtable;
   SCIP_PRESOLDATA* presoldata;
   SCIP_Longint maxuseless;
   SCIP_Longint nuseless;
   SCIP_Bool* isimpliedfrees;
   int nimpliedfrees;
   SCIP_Bool* isblockedvar;
   SCIP_VAR** vars;
   int oldnchgcoefs;
   int numcancel;
   int nfillin;

   int *processedvarsidx;
   int *processedvarsnnz;
   SCIP_Bool *processedvarssign;
   SCIP_Bool *processedvarsisfree;
   int nprocessedvarsidx;

   nprocessedvarsidx = 0;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);

   if( presoldata->nwaitingcalls > 0 )
   {
      presoldata->nwaitingcalls--;
      SCIPdebugMsg(scip, "skipping dualsparsify: nfailures=%d, nwaitingcalls=%d\n", presoldata->nfailures,
         presoldata->nwaitingcalls);
      return SCIP_OKAY;
   }
   SCIPdebugMsg(scip, "starting dualsparsify. . .\n");
   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized && complete )
   {
      ncols = SCIPmatrixGetNColumns(matrix);
      nimpliedfrees = 0;

      /* sort column by row indices */
      for( i = 0; i < ncols; i++ )
      {
         int* colpnt = SCIPmatrixGetColIdxPtr(matrix, i);
         SCIP_Real* valpnt = SCIPmatrixGetColValPtr(matrix, i);
         SCIPsortIntReal(colpnt, valpnt, SCIPmatrixGetColNNonzs(matrix, i));
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &scores, SCIPmatrixGetNRows(matrix)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &perm, SCIPmatrixGetNRows(matrix)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &isimpliedfrees, SCIPmatrixGetNColumns(matrix)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, SCIPmatrixGetNColumns(matrix)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &isblockedvar, SCIPmatrixGetNColumns(matrix)) );

      /* loop over all columns and create column pairs */
      conspairssize = 0;
      nconspairs = 0;
      conspairs = NULL;
      SCIP_CALL( SCIPhashtableCreate(&pairtable, SCIPblkmem(scip), 1, SCIPhashGetKeyStandard, consPairsEqual, consPairHashval, (void*) scip) );

      /* collect implied free variables and their number of non-zeros */
      for( c = 0; c < ncols; c++ )
      {
         int nnonz;
         int lbimplied;
         int ubimplied;

         nnonz = SCIPmatrixGetColNNonzs(matrix, c);
         lbimplied = isLowerBoundImplied(scip, matrix, c);
         ubimplied = isUpperBoundImplied(scip, matrix, c);
         vars[c] = SCIPmatrixGetVar(matrix, c);
         isimpliedfrees[c] = FALSE;

         if( lbimplied && ubimplied )
         {
            nimpliedfrees += 1;
            isimpliedfrees[c] = TRUE;
         }
         isblockedvar[c] = FALSE;

         /* only consider implied free variables
          * skip singleton variables, because either the constraint is redundant
          * or the variables can be cancelled by variables substitution
          */
         if( nnonz >= 2 && (lbimplied && ubimplied) )
         {
            int* colinds;
            SCIP_Real* colvals;
            int npairs;
            int failshift;

            colinds = SCIPmatrixGetColIdxPtr(matrix, c);
            colvals = SCIPmatrixGetColValPtr(matrix, c);

            for( i = 0; i < nnonz; ++i )
            {
               perm[i] = i;
               scores[i] = -SCIPmatrixGetRowNNonzs(matrix, colinds[i]) - 1.0*colinds[i]/ncols;
            }

            SCIPsortRealInt(scores, perm, nnonz);

            if( presoldata->maxconsiderednonzeros >= 0 )
               nnonz = MIN(nnonz, presoldata->maxconsiderednonzeros);

            npairs = (nnonz * (nnonz - 1)) / 2;
            if( nconspairs + npairs > conspairssize )
            {
               int newsize = SCIPcalcMemGrowSize(scip, nconspairs + npairs);
               SCIP_CALL( SCIPreallocBufferArray(scip, &conspairs, newsize) );
               conspairssize = newsize;
            }

            /* if we are called after one or more failures, i.e., executions without finding cancellations, then we
             * shift the section of nonzeros considered; in the case that the maxconsiderednonzeros limit is hit, this
             * results in different variable pairs being tried and avoids trying the same useless cancellations
             * repeatedly
             */
            failshift = presoldata->nfailures*presoldata->maxconsiderednonzeros;

            for( i = 0; i < nnonz; ++i )
            {
               for( j = i + 1; j < nnonz; ++j )
               {
                  int i1;
                  int i2;

                  assert(nconspairs < conspairssize);
                  assert(conspairs != NULL);

                  i1 = perm[(i + failshift) % nnonz];
                  i2 = perm[(j + failshift) % nnonz];
                  conspairs[nconspairs].colindex = c;

                 if( colinds[i1] < colinds[i2])
                  {
                     conspairs[nconspairs].consindex1 = colinds[i1];
                     conspairs[nconspairs].consindex2 = colinds[i2];
                     conspairs[nconspairs].conscoef1 = colvals[i1];
                     conspairs[nconspairs].conscoef2 = colvals[i2];
                  }
                  else
                  {
                     conspairs[nconspairs].consindex1 = colinds[i2];
                     conspairs[nconspairs].consindex2 = colinds[i1];
                     conspairs[nconspairs].conscoef1 = colvals[i2];
                     conspairs[nconspairs].conscoef2 = colvals[i1];
                  }
                  ++nconspairs;
               }
            }
         }
      }

      {
      if( SCIPpresolGetNCalls(presol) <= 0 && nimpliedfrees != 0 )
         printf( "The number of implied free variables %d\n", nimpliedfrees );
      }

      /* insert conspairs into hash table */
      for( c = 0; c < nconspairs; ++c )
      {
         SCIP_Bool insert;
         COLCONSPAIR* otherconspair;

         assert(conspairs != NULL);

         insert = TRUE;

         /* check if this pair is already contained in the hash table;
          * The loop is required due to the non-transitivity of the hash functions
          */
         while( (otherconspair = (COLCONSPAIR*)SCIPhashtableRetrieve(pairtable, (void*) &conspairs[c])) != NULL )
         {
            /* if the previous variable pair has fewer or the same number of non-zeros in the attached row
             * we keep that pair and skip this one
             */
            if( SCIPmatrixGetColNNonzs(matrix, otherconspair->colindex) <= SCIPmatrixGetColNNonzs(matrix, conspairs[c].colindex) )
            {
               insert = FALSE;
               break;
            }

            /* this pairs row has fewer non-zeros, so remove the other pair from the hash table and loop */
            SCIP_CALL( SCIPhashtableRemove(pairtable, (void*) otherconspair) );
         }

         if( insert )
         {
            SCIP_CALL( SCIPhashtableInsert(pairtable, (void*) &conspairs[c]) );
         }
      }

      /* sort cols according to decreasingly sparsity */
      SCIP_CALL( SCIPallocBufferArray(scip, &colidxsorted, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colsparsity, ncols) );
      for( c = 0; c < ncols; ++c )
         colidxsorted[c] = c;
      for( c = 0; c < ncols; ++c )
         colsparsity[c] = -SCIPmatrixGetColNNonzs(matrix, c);
      SCIPsortIntInt(colsparsity, colidxsorted, ncols);


      /* loop over the columns and cancel non-zeros until maximum number of retrieves is reached */
      maxuseless = (SCIP_Longint)(presoldata->maxretrievefac * (SCIP_Real)ncols);
      nuseless = 0;
      oldnchgcoefs = *nchgcoefs;
      for( c = 0; c < ncols && nuseless <= maxuseless; c++ )
      {
         int colidx;

         colidx = colidxsorted[c];

         if( isblockedvar[colidx] )
            continue;


         /* since the function parameters for the max fillin are unsigned we do not need to handle the
          * unlimited (-1) case due to implicit conversion rules */
         SCIP_CALL( cancelCol(scip, matrix, presoldata, pairtable, isimpliedfrees, vars, isblockedvar, colidx, \
               presoldata->maxcontfillin == -1 ? INT_MAX : presoldata->maxcontfillin, \
               presoldata->maxintfillin == -1 ? INT_MAX : presoldata->maxintfillin, \
               presoldata->maxbinfillin == -1 ? INT_MAX : presoldata->maxbinfillin, \
               presoldata->maxconsiderednonzeros, presoldata->preserveintcoefs, \
               &nuseless, nchgcoefs, &numcancel, &nfillin, FALSE) );
      }

      if( numcancel > 0 )
      {
         *result = SCIP_SUCCESS;
      }
      else
      {
         SCIPhashtableRemoveAll(pairtable);
         nconspairs = 0;

         /* collect large nonzero entries variables and their number of non-zeros */
         for( c = 0; c < ncols; c++ )
         {
            int nnonz;
            nnonz = SCIPmatrixGetColNNonzs(matrix, c);
            vars[c] = SCIPmatrixGetVar(matrix, c);

            isblockedvar[c] = FALSE;

            /* only consider large nonzero entries and nonimplied free variables
             * skip singleton variables, because either the constraint is redundant
             * or the variables can be cancelled by variables substitution
             */
            if( nnonz >= presoldata->mincancelnnzsaddcons && !isimpliedfrees[c] )
            {
               int* colinds;
               SCIP_Real* colvals;
               int npairs;
               int failshift;

               colinds = SCIPmatrixGetColIdxPtr(matrix, c);
               colvals = SCIPmatrixGetColValPtr(matrix, c);

               for( i = 0; i < nnonz; ++i )
               {
                  perm[i] = i;
                  scores[i] = -SCIPmatrixGetRowNNonzs(matrix, colinds[i]) - 1.0*colinds[i]/ncols;
               }

               SCIPsortRealInt(scores, perm, nnonz);

               if( presoldata->maxconsiderednonzeros >= 0 )
                  nnonz = MIN(nnonz, presoldata->maxconsiderednonzeros);

               npairs = (nnonz * (nnonz - 1)) / 2;
               if( nconspairs + npairs > conspairssize )
               {
                  int newsize = SCIPcalcMemGrowSize(scip, nconspairs + npairs);
                  SCIP_CALL( SCIPreallocBufferArray(scip, &conspairs, newsize) );
                  conspairssize = newsize;
               }

               /* if we are called after one or more failures, i.e., executions without finding cancellations, then we
                * shift the section of nonzeros considered; in the case that the maxconsiderednonzeros limit is hit, this
                * results in different variable pairs being tried and avoids trying the same useless cancellations
                * repeatedly
                */
               failshift = presoldata->nfailures*presoldata->maxconsiderednonzeros;

               for( i = 0; i < nnonz; ++i )
               {
                  for( j = i + 1; j < nnonz; ++j )
                  {
                     int i1;
                     int i2;

                     assert(nconspairs < conspairssize);
                     assert(conspairs != NULL);

                     i1 = perm[(i + failshift) % nnonz];
                     i2 = perm[(j + failshift) % nnonz];
                     conspairs[nconspairs].colindex = c;

                    if( colinds[i1] < colinds[i2])
                     {
                        conspairs[nconspairs].consindex1 = colinds[i1];
                        conspairs[nconspairs].consindex2 = colinds[i2];
                        conspairs[nconspairs].conscoef1 = colvals[i1];
                        conspairs[nconspairs].conscoef2 = colvals[i2];
                     }
                     else
                     {
                        conspairs[nconspairs].consindex1 = colinds[i2];
                        conspairs[nconspairs].consindex2 = colinds[i1];
                        conspairs[nconspairs].conscoef1 = colvals[i2];
                        conspairs[nconspairs].conscoef2 = colvals[i1];
                     }
                     ++nconspairs;
                  }
               }
            }
         }

         /* insert conspairs into hash table */
         for( c = 0; c < nconspairs; ++c )
         {
            SCIP_Bool insert;
            COLCONSPAIR* otherconspair;

            assert(conspairs != NULL);

            insert = TRUE;

            /* check if this pair is already contained in the hash table;
             * The loop is required due to the non-transitivity of the hash functions
             */
            while( (otherconspair = (COLCONSPAIR*)SCIPhashtableRetrieve(pairtable, (void*) &conspairs[c])) != NULL )
            {
               /* if the previous variable pair has fewer or the same number of non-zeros in the attached row
                * we keep that pair and skip this one
                */
               if( SCIPmatrixGetColNNonzs(matrix, otherconspair->colindex) <= SCIPmatrixGetColNNonzs(matrix, conspairs[c].colindex) )
               {
                  insert = FALSE;
                  break;
               }

               /* this pairs row has fewer non-zeros, so remove the other pair from the hash table and loop */
               SCIP_CALL( SCIPhashtableRemove(pairtable, (void*) otherconspair) );
            }

            if( insert )
            {
               SCIP_CALL( SCIPhashtableInsert(pairtable, (void*) &conspairs[c]) );
            }
         }

         /* sort rows according to decreasingly sparsity */
         SCIP_CALL( SCIPallocBufferArray(scip, &colidxsorted, ncols) );
         SCIP_CALL( SCIPallocBufferArray(scip, &colsparsity, ncols) );
         for( c = 0; c < ncols; ++c )
            colidxsorted[c] = c;
         for( c = 0; c < ncols; ++c )
            colsparsity[c] = -SCIPmatrixGetColNNonzs(matrix, c);
         SCIPsortIntInt(colsparsity, colidxsorted, ncols);


         /* loop over the columns and cancel non-zeros until maximum number of retrieves is reached */
         maxuseless = (SCIP_Longint)(presoldata->maxretrievefac * (SCIP_Real)ncols);
         nuseless = 0;
         oldnchgcoefs = *nchgcoefs;
         for( c = 0; c < ncols && nuseless <= maxuseless; c++ )
         {
            int colidx;
            int nnonz;

            colidx = colidxsorted[c];
            nnonz = SCIPmatrixGetColNNonzs(matrix, colidx);

            if( isblockedvar[colidx] || nnonz < 100 )
               continue;



            /* since the function parameters for the max fillin are unsigned we do not need to handle the
             * unlimited (-1) case due to implicit conversion rules */
            SCIP_CALL( cancelCol(scip, matrix, presoldata, pairtable, isimpliedfrees, vars, isblockedvar, colidx, \
                  presoldata->maxcontfillin == -1 ? INT_MAX : presoldata->maxcontfillin, \
                  presoldata->maxintfillin == -1 ? INT_MAX : presoldata->maxintfillin, \
                  presoldata->maxbinfillin == -1 ? INT_MAX : presoldata->maxbinfillin, \
                  presoldata->maxconsiderednonzeros, presoldata->preserveintcoefs, \
                  &nuseless, nchgcoefs, &numcancel, &nfillin, TRUE) );
         }

         if( numcancel > 0 )
         {
            *result = SCIP_SUCCESS;
         }
      }



      updateFailureStatistic(presoldata, numcancel > 0);

      SCIPfreeBufferArray(scip, &colsparsity);
      SCIPfreeBufferArray(scip, &colidxsorted);

      SCIPhashtableFree(&pairtable);
      SCIPfreeBufferArrayNull(scip, &conspairs);

      SCIPfreeBufferArray(scip, &isblockedvar);
      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &isimpliedfrees);
      SCIPfreeBufferArray(scip, &perm);
      SCIPfreeBufferArray(scip, &scores);
   }
   /* if matrix construction fails once, we do not ever want to be called again */
   else
   {
      updateFailureStatistic(presoldata, FALSE);
      presoldata->nwaitingcalls = INT_MAX;
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeDualsparsify)
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
SCIP_DECL_PRESOLINIT(presolInitDualsparsify)
{
   SCIP_PRESOLDATA* presoldata;

   /* set the counters in the init (and not in the initpre) callback such that they persist across restarts */
   presoldata = SCIPpresolGetData(presol);
   presoldata->ncancels = 0;
   presoldata->nfillin = 0;
   presoldata->nfailures = 0;
   presoldata->nwaitingcalls = 0;
   presoldata->naggregated = 0;

   return SCIP_OKAY;
}

/** creates the dualsparsify presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualsparsify(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create dualsparsify presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecDualsparsify, presoldata) );

   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDualsparsify) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeDualsparsify) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitDualsparsify) );


   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/minacceptcancelnnzs",
         "minimal cancel nonzeros when accepting to aggregate variables",
         &presoldata->minacceptcancelnnzs, FALSE, DEFAULT_MINACCEPTCANCELNNZS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/minconsiderednnzs",
         "minimal number of considered non-zeros within one column",
         &presoldata->minconsiderednnzs, FALSE, DEFAULT_MINCONSIDEREDNNZS, 0, INT_MAX, NULL, NULL) );


   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/mincancelnnzsaddcons",
         "minimal cancelled non-zeros within one column if we need to add a constraint to the problem",
         &presoldata->mincancelnnzsaddcons, FALSE, DEFAULT_MINCANCELNNZSADDCONS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxcompareeverypair",
         "maximal number on the implementaion of doing reduction on every pair of variables",
         &presoldata->maxcompareeverypair, FALSE, DEFAULT_MAXCOMPAREDEVERYPAIR, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/dualsparsify/maxfillinrate",
         "cancel the variable if nfillins/ncancels is less than or equal this rate",
         &presoldata->maxfillinrate, TRUE, DEFAULT_MAX_FILLINRATE, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/dualsparsify/maxfillinratefree",
         "cancel the free variable if nfillins/ncancels is less than or equal this rate",
         &presoldata->maxfillinratefree, TRUE, DEFAULT_MAX_FILLINRATE, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxconsiderednonzeros",
         "maximal number of considered non-zeros within one row (-1: no limit)",
         &presoldata->maxconsiderednonzeros, TRUE, DEFAULT_MAXCONSIDEREDNONZEROS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxcontfillin",
         "maximal fillin for continuous variables (-1: unlimited)",
         &presoldata->maxcontfillin, FALSE, DEFAULT_MAX_CONT_FILLIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxbinfillin",
         "maximal fillin for binary variables (-1: unlimited)",
         &presoldata->maxbinfillin, FALSE, DEFAULT_MAX_BIN_FILLIN, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxintfillin",
         "maximal fillin for integer variables including binaries (-1: unlimited)",
         &presoldata->maxintfillin, FALSE, DEFAULT_MAX_INT_FILLIN, -1, INT_MAX, NULL, NULL) );
#if 0
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxnonzeros",
         "maximal support of one equality to be used for cancelling (-1: no limit)",
         &presoldata->maxnonzeros, TRUE, DEFAULT_MAXNONZEROS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/dualsparsify/maxconsiderednonzeros",
         "maximal number of considered non-zeros within one row (-1: no limit)",
         &presoldata->maxconsiderednonzeros, TRUE, DEFAULT_MAXCONSIDEREDNONZEROS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam(scip,
         "presolving/dualsparsify/rowsort",
         "order in which to process inequalities ('n'o sorting, 'i'ncreasing nonzeros, 'd'ecreasing nonzeros)",
         &presoldata->rowsort, TRUE, DEFAULT_ROWSORT, "nid", NULL, NULL) );
#endif
   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/dualsparsify/maxretrievefac",
         "limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints",
         &presoldata->maxretrievefac, TRUE, DEFAULT_MAXRETRIEVEFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/dualsparsify/waitingfac",
         "number of calls to wait until next execution as a multiple of the number of useless calls",
         &presoldata->waitingfac, TRUE, DEFAULT_WAITINGFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/dualsparsify/preserveintcoefs",
         "should we forbid cancellations that destroy integer coefficients?",
         &presoldata->preserveintcoefs, TRUE, DEFAULT_PRESERVEINTCOEFS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/dualsparsify/islockedlimit",
         "should locked limit be used to stop some aggregation",
         &presoldata->islockedlimit, TRUE, DEFAULT_ISLOCKEDLIMIT, NULL, NULL) );

   return SCIP_OKAY;
}
