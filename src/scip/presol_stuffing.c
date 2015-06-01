/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_stuffing.c
 * @brief  fix singleton continuous variables
 * @author Dieter Weninger
 *
 * Investigate singleton continuous variables if one can be fixed at a bound.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>


#include "scip/pub_matrix.h"
#include "presol_stuffing.h"

#define PRESOL_NAME            "stuffing"
#define PRESOL_DESC            "fix redundant singleton continuous variables"
#define PRESOL_PRIORITY             -100     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

/** type of fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,         /**< fix variable at lower bound */
   NOFIX   =  0,         /**< do not fix variable */
   FIXATUB =  1          /**< fix variable at upper bound */
};
typedef enum Fixingdirection FIXINGDIRECTION;

/*
 * Local methods
 */

/** try to fix singleton continuous variables */
static
SCIP_RETCODE singletonColumnStuffing(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIPMILPMATRIX*       matrix,             /**< matrix containing the constraints */
   FIXINGDIRECTION*      varstofix,          /**< array holding fixing information */
   int*                  nfixings            /**< number of possible fixings */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* colratios;
   SCIP_Real* colcoeffs;
   SCIP_Bool* rowprocessed;
   int* rowpnt;
   int* rowend;
   int* colindices;
   int* dummy;
   SCIP_Bool* swapped;
   SCIP_Real constant1;
   SCIP_Real constant2;
   SCIP_Real coef;
   SCIP_Real value;
   SCIP_Real boundoffset;
   SCIP_Bool tryfixing;
   int idx;
   int col;
   int row;
   int fillcnt;
   int k;
   int nrows;
   int ncols;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(varstofix != NULL);
   assert(nfixings != NULL);

   nrows = SCIPmatrixGetNRows(matrix);
   ncols = SCIPmatrixGetNColumns(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &colindices, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colratios, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colcoeffs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &dummy, ncols) );

   SCIP_CALL( SCIPallocBufferArray(scip, &swapped, ncols) );
   BMSclearMemoryArray(swapped, ncols);

   SCIP_CALL( SCIPallocBufferArray(scip, &rowprocessed, nrows) );
   BMSclearMemoryArray(rowprocessed, nrows);

   for( col = 0; col < ncols; col++ )
   {
      /* consider only at rows with minimal one continuous singleton column */
      if( SCIPmatrixGetColNNonzs(matrix, col) == 1 &&
          SCIPvarGetType(SCIPmatrixGetVar(matrix, col)) == SCIP_VARTYPE_CONTINUOUS )
      {
         row = *(SCIPmatrixGetColIdxPtr(matrix, col));
         if( rowprocessed[row] )
            continue;

         rowprocessed[row] = TRUE;

         if( SCIPmatrixIsRowRhsInfinity(matrix, row) )
         {
            /* consider >= relation with obj > 0 and coef > 0 */
            fillcnt = 0;
            tryfixing = TRUE;
            constant1 = 0.0;
            constant2 = 0.0;

            rowpnt = SCIPmatrixGetRowIdxPtr(matrix, row);
            rowend = rowpnt + SCIPmatrixGetRowNNonzs(matrix, row);
            valpnt = SCIPmatrixGetRowValPtr(matrix, row);

            for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
            {
               SCIP_VAR* var;

               coef = *valpnt;
               idx = *rowpnt;
               var = SCIPmatrixGetVar(matrix, idx);

               if( SCIPmatrixGetColNNonzs(matrix, idx) == 1 && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               {
                  if( SCIPvarGetObj(var) > 0 && coef > 0 )
                  {
                     constant1 += coef * SCIPvarGetLbGlobal(var);
                     constant2 += coef * SCIPvarGetLbGlobal(var);
                     colratios[fillcnt] = SCIPvarGetObj(var) / coef;
                     colindices[fillcnt] = idx;
                     colcoeffs[fillcnt] = coef;
                     fillcnt++;
                  }
                  else if( SCIPvarGetObj(var) < 0 && coef < 0 )
                  {
                     /* swap column and bounds */
                     swapped[idx] = TRUE;
                     constant1 += -coef * -SCIPvarGetUbGlobal(var);
                     constant2 += -coef * -SCIPvarGetUbGlobal(var);
                     colratios[fillcnt] = -SCIPvarGetObj(var) / -coef;
                     colindices[fillcnt] = idx;
                     colcoeffs[fillcnt] = -coef;
                     fillcnt++;
                  }
                  else if( SCIPvarGetObj(var) > 0 && coef < 0 )
                  {
                     if( SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)) )
                     {
                        /* unbounded */
                        tryfixing = FALSE;
                        break;
                     }
                     else
                     {
                        constant1 += coef * SCIPvarGetLbGlobal(var);
                        constant2 += coef * SCIPvarGetLbGlobal(var);
                     }
                  }
                  else
                  {
                     /* obj < 0 and coef > 0 */
                     if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
                     {
                        /* unbounded */
                        tryfixing = FALSE;
                        break;
                     }
                     else
                     {
                        constant1 += coef * SCIPvarGetUbGlobal(var);
                        constant2 += coef * SCIPvarGetUbGlobal(var);
                     }
                  }
               }
               else
               {
                  /* discrete variables and non-singleton continuous variables */
                  if( coef > 0 )
                  {
                     constant1 += coef * SCIPvarGetUbGlobal(var);
                     constant2 += coef * SCIPvarGetLbGlobal(var);
                  }
                  else
                  {
                     constant1 += coef * SCIPvarGetLbGlobal(var);
                     constant2 += coef * SCIPvarGetUbGlobal(var);
                  }
               }
            }

            if( tryfixing )
            {
               SCIPsortRealRealIntInt(colratios, colcoeffs, colindices, dummy, fillcnt);

               /* verify which variable can be fixed */
               for( k = 0; k < fillcnt; k++ )
               {
                  SCIP_VAR* var;

                  idx = colindices[k];
                  var = SCIPmatrixGetVar(matrix, idx);

                  assert(SCIPmatrixGetColNNonzs(matrix, idx) == 1 && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

                  if( swapped[idx] )
                  {
                     value = colcoeffs[k] * (-SCIPvarGetLbGlobal(var));
                     boundoffset = colcoeffs[k] * (-SCIPvarGetUbGlobal(var));
                  }
                  else
                  {
                     value = colcoeffs[k] * SCIPvarGetUbGlobal(var);
                     boundoffset = colcoeffs[k] * SCIPvarGetLbGlobal(var);
                  }

                  if( SCIPisLE(scip, value, SCIPmatrixGetRowLhs(matrix, row) - constant1 + boundoffset) )
                  {
                     if( swapped[idx] )
                        varstofix[idx] = FIXATLB;
                     else
                        varstofix[idx] = FIXATUB;
                     (*nfixings)++;
                  }
                  else if( SCIPisLE(scip, SCIPmatrixGetRowLhs(matrix, row), constant2) )
                  {
                     if( swapped[idx] )
                        varstofix[idx] = FIXATUB;
                     else
                        varstofix[idx] = FIXATLB;
                     (*nfixings)++;
                  }

                  constant1 += value;
                  constant1 -= boundoffset;

                  constant2 += value;
                  constant2 -= boundoffset;
               }
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &rowprocessed);
   SCIPfreeBufferArray(scip, &swapped);
   SCIPfreeBufferArray(scip, &dummy);
   SCIPfreeBufferArray(scip, &colcoeffs);
   SCIPfreeBufferArray(scip, &colratios);
   SCIPfreeBufferArray(scip, &colindices);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyStuffing)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolStuffing(scip) );

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecStuffing)
{  /*lint --e{715}*/
   SCIPMILPMATRIX* matrix;
   SCIP_Bool initialized;
   SCIP_Bool complete;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   if( SCIPgetNContVars(scip) == 0 || SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   matrix = NULL;
   SCIP_CALL( SCIPmatrixCreate(scip, &matrix, &initialized, &complete) );

   if( initialized )
   {
      FIXINGDIRECTION* varstofix;
      int nfixings;
      int ncols;

      nfixings = 0;
      ncols = SCIPmatrixGetNColumns(matrix);

      SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, ncols) );
      BMSclearMemoryArray(varstofix, ncols);

      SCIP_CALL( singletonColumnStuffing(scip, matrix, varstofix, &nfixings) );

      if( nfixings > 0 )
      {
         int v;
         int oldnfixedvars;

         oldnfixedvars = *nfixedvars;

         /* look for fixable variables */
         for( v = ncols - 1; v >= 0; --v )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            SCIP_VAR* var;

            var = SCIPmatrixGetVar(matrix, v);

            if( SCIPvarGetNLocksUp(var) != SCIPmatrixGetColNUplocks(matrix, v) ||
               SCIPvarGetNLocksDown(var) != SCIPmatrixGetColNDownlocks(matrix, v) )
            {
               /* no fixing, locks not consistent */
               continue;
            }

            if( varstofix[v] == FIXATLB )
            {
               SCIP_Real lb;

               lb = SCIPvarGetLbGlobal(var);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

               /* fix at lower bound */
               SCIP_CALL( SCIPfixVar(scip, var, lb, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage(" -> infeasible fixing\n");
                  *result = SCIP_CUTOFF;

                  break;
               }
               assert(fixed);
               (*nfixedvars)++;
            }
            else if( varstofix[v] == FIXATUB )
            {
               SCIP_Real ub;

               ub = SCIPvarGetUbGlobal(var);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

               /* fix at upper bound */
               SCIP_CALL( SCIPfixVar(scip, var, ub, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage(" -> infeasible fixing\n");
                  *result = SCIP_CUTOFF;

                  break;
               }
               assert(fixed);

               (*nfixedvars)++;
            }
         }

         if( *result != SCIP_CUTOFF && *nfixedvars > oldnfixedvars )
            *result = SCIP_SUCCESS;
      }

      SCIPfreeBufferArray(scip, &varstofix);
   }

   SCIPmatrixFree(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the stuffing presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolStuffing(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecStuffing, NULL) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyStuffing) );

   return SCIP_OKAY;
}

/*lint --e{749}*/
