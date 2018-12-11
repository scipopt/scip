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
 * @brief  cancel non-zeros of the constraint matrix
 * @author Dieter Weninger
 * @author Robert Lion Gottwald
 * @author Ambros Gleixner
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

#define PRESOL_PRIORITY            -24000    /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS               -1    /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_ENABLECOPY           TRUE    /**< should dualsparsify presolver be copied to sub-SCIPs? */
#define DEFAULT_CANCELLINEAR         TRUE    /**< should we cancel nonzeros in constraints of the linear constraint handler? */
#define DEFAULT_PRESERVEINTCOEFS     TRUE    /**< should we forbid cancellations that destroy integer coefficients? */
#define DEFAULT_MAX_CONT_FILLIN         0    /**< default value for the maximal fillin for continuous variables */
#define DEFAULT_MAX_BIN_FILLIN          0    /**< default value for the maximal fillin for binary variables */
#define DEFAULT_MAX_INT_FILLIN          0    /**< default value for the maximal fillin for integer variables (including binary) */
#define DEFAULT_MAXNONZEROS            -1    /**< maximal support of one equality to be used for cancelling (-1: no limit) */
#define DEFAULT_MAXCONSIDEREDNONZEROS  70    /**< maximal number of considered non-zeros within one row (-1: no limit) */
#define DEFAULT_ROWSORT               'd'    /**< order in which to process inequalities ('n'o sorting, 'i'ncreasing nonzeros, 'd'ecreasing nonzeros) */
#define DEFAULT_MAXRETRIEVEFAC      100.0    /**< limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints */
#define DEFAULT_WAITINGFAC            2.0    /**< number of calls to wait until next execution as a multiple of the number of useless calls */

#define MAXSCALE                   1000.0    /**< maximal allowed scale for cancelling non-zeros */
#define MINSCALE                    0.001    /**< minimal allowed scale for cancelling non-zeros */


#define DEFAULT_MINACCEPTCANCELNNZS    10    /**< minimal cancel nonzeros when accepting to aggregate the (nonfree) variable */
#define DEFAULT_MINCONSIDEREDNNZS      10    /**< minimal number of considered non-zeros within one column */
#define DEFAULT_MAXCOMPAREDEVERYPAIR  200    /**< maximal number on the implementaion of doing reduction on every pair of variables */
#define DEFAULT_MAX_FILLINRATE       0.05    /**< cancel the variable if nfillins/ncancels is less than or equal this rate */
#define DEFAULT_MAX_FILLINRATE_FREE   0.1    /**< cancel the free variable if nfillins/ncancels is less than or equal this rate */

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
   int                   naggregated;

   int                   minacceptcancelnnzs;/**< minimal cancel nonzeros when accepting to aggregate the (nonfree) variable */
   int                   minconsiderednnzs;  /**< minimal number of considered non-zeros within one column */
   int                   maxcompareeverypair;/**< maximal number on the implementaion of doing reduction on every pair of variables */
   SCIP_Real             maxfillinrate;      /**< cancel the variable if nfillins/ncancels is less than or equal this rate */
   SCIP_Real             maxfillinratefree;  /**< cancel the free variable if nfillins/ncancels is less than or equal this rate */
};

/*
 * Local methods
 */

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


/** verify which variable upper bounds are implied */
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

   varub = SCIPmatrixGetColUb(matrix, col);

   if( !SCIPisInfinity(scip, varub) && SCIPisFeasLE(scip, impliedub, varub) )
      ubimplied = TRUE;

   return ubimplied;
}

/** verify which variable lower bounds are implied */
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

   varlb = SCIPmatrixGetColLb(matrix, col);

   if( !SCIPisInfinity(scip, -varlb) && SCIPisFeasGE(scip, impliedlb, varlb) )
      lbimplied = TRUE;

   return lbimplied;
}


/* add variable colidx1 to variable colidx2 */
static
SCIP_RETCODE aggregation(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_MATRIX*          matrix,             /**< the constraint matrix */
   int                   colidx1,            /**< one of the indexes of column to try non-zero cancellation for */
   int                   colidx2,            /**< one of the indexes of column to try non-zero cancellation for */
   SCIP_Bool             isimpliedfree,      /**< is the aggregated variable implied free? */
   SCIP_Real             weight1             /**< weight variable one in the aggregated expression */
      )
{
   SCIP_VAR* vars[2];
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

   assert( !SCIPisZero(scip, weight1) );

   presoldata->naggregated += 1;
   aggregatedvar = SCIPmatrixGetVar(matrix, colidx2);

   (void) SCIPsnprintf(newvarname, SCIP_MAXSTRLEN, "dualsparsifyvar_%d", presoldata->naggregated);
   /* TODO: we need to consider the infinite bounds here */
   if( weight1 > 0 )
   {
      if(SCIPisInfinity(scip, -SCIPmatrixGetColLb(matrix, colidx1)) || SCIPisInfinity(scip, -SCIPmatrixGetColLb(matrix, colidx2)))
         newlb = -SCIPinfinity(scip);
      else
         newlb = weight1*SCIPmatrixGetColLb(matrix, colidx1) + SCIPmatrixGetColLb(matrix, colidx2);

      if(SCIPisInfinity(scip, SCIPmatrixGetColUb(matrix, colidx1)) || SCIPisInfinity(scip, SCIPmatrixGetColUb(matrix, colidx2)))
         newub = SCIPinfinity(scip);
      else
         newub = weight1*SCIPmatrixGetColUb(matrix, colidx1) + SCIPmatrixGetColUb(matrix, colidx2);
   }

   else
   {
      if(SCIPisInfinity(scip, SCIPmatrixGetColUb(matrix, colidx1)) || SCIPisInfinity(scip, -SCIPmatrixGetColLb(matrix, colidx2)))
         newlb = -SCIPinfinity(scip);
      else
         newlb = weight1*SCIPmatrixGetColUb(matrix, colidx1) + SCIPmatrixGetColLb(matrix, colidx2);

      if(SCIPisInfinity(scip, SCIPmatrixGetColLb(matrix, colidx1)) || SCIPisInfinity(scip, SCIPmatrixGetColUb(matrix, colidx2)))
         newub = SCIPinfinity(scip);
      else
         newub = weight1*SCIPmatrixGetColLb(matrix, colidx1) + SCIPmatrixGetColUb(matrix, colidx2);
   }

   SCIP_CALL( SCIPcreateVar(scip, &newvar, newvarname, newlb, newub, 0.0, SCIP_VARTYPE_CONTINUOUS,
            SCIPvarIsInitial(aggregatedvar), SCIPvarIsRemovable(aggregatedvar), NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(scip, newvar) );

   vars[0] = SCIPmatrixGetVar(matrix, colidx1);
   vars[1] = newvar;
   coefs[0] = -weight1;
   coefs[1] = 1;

   SCIP_CALL( SCIPmultiaggregateVar(scip, aggregatedvar, 2, vars, coefs, 0.0, &infeasible, &aggregated) );
   assert(!infeasible);
   assert(aggregated);

   if( !isimpliedfree )
   {
      (void) SCIPsnprintf(newconsname, SCIP_MAXSTRLEN, "dualsparsifycons_%d", presoldata->naggregated);

      SCIP_CALL( SCIPcreateConsLinear(scip, &newcons, newconsname, 2, vars, coefs,
               SCIPmatrixGetColLb(matrix, colidx2), SCIPmatrixGetColUb(matrix, colidx2),
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, newcons) );
      SCIPdebugPrintCons(scip, newcons, NULL);

      SCIP_CALL( SCIPreleaseCons(scip, &newcons) );
   }
   SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
}

static
SCIP_RETCODE cancelCol(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   SCIP_MATRIX*          matrix,             /**< the constraint matrix */
   int                   colidx1,            /**< one of the indexes of column to try non-zero cancellation for */
   int                   colidx2,            /**< one of the indexes of column to try non-zero cancellation for */
   SCIP_Bool             isimpliedfree1,     /**< is the variable one implied free */
   SCIP_Bool             isimpliedfree2,     /**< is the variable two implied free */
   SCIP_Bool*            success,            /**< pointer to store whether the aggregations succeed or not */
   SCIP_Real*            ratios,             /**< ratio of the vectors*/
   int*                  nratios,            /**< number of different ratios*/
   int*                  nchgcoefs,          /**< pointer to update number of changed coefficients */
   int*                  ncanceled,          /**< pointer to update number of canceled nonzeros */
   int*                  nfillin,            /**< pointer to update the produced fill-in */
   int*                  naddconss         /**< pointer to update the number of added constraint */
      )
{
   int i;
   int j;
   int varlen1;
   int varlen2;
   int *inds1;
   int *inds2;
   int nnz1;
   int nnz2;
   SCIP_VAR* var1;
   SCIP_VAR* var2;
   SCIP_Real *vals1;
   SCIP_Real *vals2;

   varlen1 =  SCIPmatrixGetColNNonzs(matrix, colidx1);
   varlen2 =  SCIPmatrixGetColNNonzs(matrix, colidx2);
   assert(varlen1 >= presoldata->minconsiderednnzs);
   assert(varlen2 >= presoldata->minconsiderednnzs);

   var1 = SCIPmatrixGetVar(matrix, colidx1);
   var2 = SCIPmatrixGetVar(matrix, colidx2);
   inds1 = SCIPmatrixGetColIdxPtr(matrix, colidx1);
   inds2 = SCIPmatrixGetColIdxPtr(matrix, colidx2);
   vals1 = SCIPmatrixGetColValPtr(matrix, colidx1);
   vals2 = SCIPmatrixGetColValPtr(matrix, colidx2);
   
   i = 0;
   j = 0;
   nnz1 = 0;
   nnz2 = 0;
   *nratios = 0;
   *success = FALSE;
   while(i < varlen1 && j < varlen2)
   {
      if(inds1[i] == inds2[j])
      {
         if( SCIPisZero(scip, vals1[i]) )
         {
            nnz2 += 1;
         }
         else if( SCIPisZero(scip, vals2[j]) )
         {
            nnz1 += 1;
         }
         else
         {
            ratios[*nratios] = vals1[i]/vals2[j];
            (*nratios)++;
         }
         i++;
         j++;
      }
      else if(inds1[i] < inds2[j])
      {
         i++;
         if( !SCIPisZero(scip, vals1[i]) )
            nnz1 += 1;
      }
      else
      {
         j++;
         if( !SCIPisZero(scip, vals2[j]) )
            nnz2 += 1;
      }
   }
   if( i < varlen1 )
      nnz1 += varlen1 - i;
   if( j < varlen2 )
      nnz2 += varlen2 - j;

   if( nnz1<*nratios || nnz2<*nratios )
   {
      SCIP_Real maxratio;
      SCIP_Real curratio;
      int nmaxratio;
      int ncurratio;
      int tmp;
      SCIP_Bool tmp_sign;
      int maxfillin1;
      int maxfillin2;
      int mincancel1;
      int mincancel2;
      SCIP_Bool isaccepted1;
      SCIP_Bool isaccepted2;

      SCIP_Bool isexecute;

      nmaxratio = 0;

      SCIPsortDownReal(ratios, *nratios);
      curratio = ratios[0];
      ncurratio = 1;
      for( i=1; i<*nratios; i++ )
      {
         if( SCIPisEQ(scip, curratio, ratios[i]) )
         {
            ncurratio++;
         }
         else
         {
            if( ncurratio > nmaxratio )
            {
               maxratio = curratio;
               nmaxratio = ncurratio;
            }
            else
            {
               curratio = ratios[i];
               ncurratio = 1;
            }
         }
      }

      if( ncurratio > nmaxratio )
      {
         maxratio = curratio;
         nmaxratio = ncurratio;
      }

      maxfillin1 = isimpliedfree1 ? floor(presoldata->maxfillinratefree*nmaxratio) : floor(presoldata->maxfillinrate*nmaxratio);
      maxfillin2 = isimpliedfree2 ? floor(presoldata->maxfillinratefree*nmaxratio) : floor(presoldata->maxfillinrate*nmaxratio);
      mincancel1 = isimpliedfree1 ? 0 : presoldata->minacceptcancelnnzs;
      mincancel2 = isimpliedfree2 ? 0 : presoldata->minacceptcancelnnzs;
      isaccepted1 = (maxfillin1 >= nnz1) && (nmaxratio >= nnz1 + mincancel1) && (REALABS(1.0/maxratio) > MINSCALE)
         && (REALABS(1.0/maxratio) < MAXSCALE);
      isaccepted2 = (maxfillin2 >= nnz2) && (nmaxratio >= nnz2 + mincancel2) && (REALABS(maxratio) > MINSCALE)
         && (REALABS(maxratio) < MAXSCALE);
      isexecute = TRUE;
 
      if( isaccepted1 && isaccepted2 )
      {
         if((isimpliedfree1 && isimpliedfree2) || (!isimpliedfree1 && !isimpliedfree2))
         {
            tmp_sign = nnz1 < nnz2 ? TRUE : FALSE;
         }
         else if( isimpliedfree1 )
         {
            tmp_sign = TRUE;
         }
         else
         {
            tmp_sign = FALSE;
         }
      }
      else if( isaccepted1 )
      {
         tmp_sign = TRUE;
      }
      else if( isaccepted2 )
      {
         tmp_sign = FALSE;
      }
      else
      {
         isexecute = FALSE;
      }

      if( isexecute )
      {
         if( tmp_sign )
         {
            *success = TRUE;

            SCIPdebugMsg(scip, "use dualsparsify method on varaibles %s and %s\n",
                  SCIPmatrixGetColName(matrix, colidx1), SCIPmatrixGetColName(matrix, colidx2));

            /* aggregate variable one */
            SCIP_CALL( aggregation(scip, presoldata, matrix, colidx2, colidx1, isimpliedfree1, 1.0/maxratio) );

            if( !isimpliedfree1 )
               *naddconss += 1;

            *ncanceled = *ncanceled + nmaxratio - nnz1;
            *nchgcoefs = *nchgcoefs + varlen2;
         }
         else
         {
            *success = TRUE;

            SCIPdebugMsg(scip, "use dualsparsify method on varaibles %s and %s\n",
                  SCIPmatrixGetColName(matrix, colidx1), SCIPmatrixGetColName(matrix, colidx2));

            /* aggregate variable two */
            SCIP_CALL( aggregation(scip, presoldata, matrix, colidx1, colidx2, isimpliedfree2, maxratio) );

            if( !isimpliedfree2 )
               *naddconss += 1;

            *ncanceled = *ncanceled + nmaxratio - nnz2;
            *nchgcoefs = *nchgcoefs + varlen1;
         }
      }
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
   int r;
   int i;
   int j;
   int ncancels;
   int nfillins;
   int nchgcoef;
   int ncancel;
   int* locks;
   int* rowidxsorted;
   int* rowsparsity;
   int nvarpairs;
   int varpairssize;
   SCIP_PRESOLDATA* presoldata;
   SCIP_Longint maxuseless;
   SCIP_Longint nuseless;

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
      SCIP_Real *ratios;
      int nratios;

      ncols = SCIPmatrixGetNColumns(matrix);
      nrows = SCIPmatrixGetNRows(matrix);
      nratios = 0;
      ncancels = 0;
      nchgcoef = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &processedvarsidx, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &processedvarsnnz, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &processedvarssign, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &processedvarsisfree, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &ratios, nrows) );

      for(i=0; i<ncols; i++)
      {
         int nnonz;
         SCIP_Bool lbimplied;
         SCIP_Bool ubimplied;

         nnonz = SCIPmatrixGetColNNonzs(matrix, i);
         if(nnonz > presoldata->minconsiderednnzs && SCIPvarGetType(SCIPmatrixGetVar(matrix, i)) == SCIP_VARTYPE_CONTINUOUS
               && !SCIPdoNotMultaggrVar(scip, SCIPmatrixGetVar(matrix, i)))
         {
            processedvarsidx[nprocessedvarsidx] = i;
            processedvarsnnz[nprocessedvarsidx] = nnonz;
            processedvarssign[nprocessedvarsidx] = TRUE;
            processedvarsisfree[nprocessedvarsidx] = FALSE;

            lbimplied = isLowerBoundImplied(scip, matrix, i);
            ubimplied = isUpperBoundImplied(scip, matrix, i);

            if( lbimplied && ubimplied )
            {
               processedvarsisfree[nprocessedvarsidx] = TRUE;
               SCIPdebugMsg(scip, "variable %s is implied free variable\n", SCIPmatrixGetColName(matrix, i));
            }
            nprocessedvarsidx++;
         }
      }

      for(i=0; i<nprocessedvarsidx; i++)
      {
         int* colpnt = SCIPmatrixGetColIdxPtr(matrix, processedvarsidx[i]);
         SCIP_Real* valpnt = SCIPmatrixGetColValPtr(matrix, processedvarsidx[i]);
         SCIPsortIntReal(colpnt, valpnt, SCIPmatrixGetColNNonzs(matrix, processedvarsidx[i]));
      }

      SCIPsortIntIntInt(processedvarsnnz, processedvarsidx, processedvarsisfree, nprocessedvarsidx);

      /* compare every pair of variables if the number of considered variables is small enough */
      if( nprocessedvarsidx > presoldata->maxcompareeverypair )
      {
         for(i=0; i<nprocessedvarsidx; i++)
         {
            /* only consider the variable that has not aggregated or used to aggregated */
            if( processedvarssign[i] )
            {
               for(j=i+1; j<nprocessedvarsidx; j++)
               {
                  SCIP_Bool success;

                  if( !processedvarssign[j] )
                     continue;

                  nchgcoef = 0;
                  ncancel = 0;
                  cancelCol(scip, presoldata, matrix, processedvarsidx[i], processedvarsidx[j], processedvarsisfree[i], processedvarsisfree[j],
                        &success, ratios, &nratios, &nchgcoef, &ncancel, &nfillins, naddconss);

                  if( success )
                  {
                     ncancels += ncancel;
                     *nchgcoefs += nchgcoef;
                     *naggrvars += 1;
                  //   *naddconss += 1;
                     processedvarssign[i] = FALSE;
                     processedvarssign[j] = FALSE;
                     break;
                  }
               }
            }
         }
      }
      /* compare only the neighbor pair of variables if the number of considered variables is too large */
      else
      {
         i=0;
         while(i < nprocessedvarsidx - 1)
         {
            SCIP_Bool success;

            nchgcoef = 0;
            ncancel = 0;
            cancelCol(scip, presoldata, matrix, processedvarsidx[i], processedvarsidx[i+1], processedvarsisfree[i], processedvarsisfree[i+1],
                  &success, ratios, &nratios, &nchgcoef, &ncancel, &nfillins, naddconss);

            if( success )
            {
               ncancels += ncancel;
               *nchgcoefs += nchgcoef;
               *naggrvars += 1;
         //      *naddconss += 1;
               i += 2;
            }
            else
               i++;
         }
      }

      presoldata->ncancels += ncancels;

      if( ncancels > 0 )
      {
         *result = SCIP_SUCCESS;
      }
      /* if matrix construction fails once, we do not ever want to be called again */
      else
      {
         presoldata->nwaitingcalls = INT_MAX;
      }

      SCIPfreeBufferArray(scip, &ratios);
      SCIPfreeBufferArray(scip, &processedvarsidx);
      SCIPfreeBufferArray(scip, &processedvarsnnz);
      SCIPfreeBufferArray(scip, &processedvarssign);
      SCIPfreeBufferArray(scip, &processedvarsisfree);
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
#if 0
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

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/dualsparsify/maxretrievefac",
         "limit on the number of useless vs. useful hashtable retrieves as a multiple of the number of constraints",
         &presoldata->maxretrievefac, TRUE, DEFAULT_MAXRETRIEVEFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/dualsparsify/waitingfac",
         "number of calls to wait until next execution as a multiple of the number of useless calls",
         &presoldata->waitingfac, TRUE, DEFAULT_WAITINGFAC, 0.0, SCIP_REAL_MAX, NULL, NULL) );
#endif

   return SCIP_OKAY;
}
