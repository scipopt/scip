/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    presol_domcol.c
 * @ingroup PRESOLVERS
 * @brief   dominated column presolver
 * @author  Dieter Weninger
 * @author  Gerald Gamrath
 *
 * This presolver looks for dominance relations between variable pairs.
 * From a dominance relation and certain bound/clique-constellations
 * variable fixings mostly at the lower bound of the dominated variable can be derived.
 *
 * @todo also run on general CIPs, if the number of locks of the investigated variables
 *       comes only from (upgraded) linear constraints
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

/* includes necessary for matrix data structure */
#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"

#include "presol_domcol.h"

#define PRESOL_NAME            "domcol"
#define PRESOL_DESC            "dominated column presolver"
#define PRESOL_PRIORITY         20000000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY                TRUE     /**< should presolver be delayed, if other presolvers found reductions? */

/*
 * Data structures
 */

/** type of fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,
   NOFIX   =  0,
   FIXATUB =  1
};
typedef enum Fixingdirection FIXINGDIRECTION;


/********************************************************************/
/************** matrix data structure and functions *****************/

/** constraint matrix data structure in column and row major format */
struct ConstraintMatrix
{
   SCIP_Real*            colmatval;          /**< coefficients in column major format */
   int*                  colmatind;          /**< row indexes in column major format */
   int*                  colmatbeg;          /**< column storage offset */
   int*                  colmatcnt;          /**< number of row entries per column */
   int                   ncols;              /**< complete number of columns */
   SCIP_Real*            lb;                 /**< lower bound per variable */
   SCIP_Real*            ub;                 /**< upper bound per variable */

   SCIP_VAR**            vars;               /**< variables pointer */

   SCIP_Real*            rowmatval;          /**< coefficients in row major format */
   int*                  rowmatind;          /**< column indexed in row major format */
   int*                  rowmatbeg;          /**< row storage offset */
   int*                  rowmatcnt;          /**< number of column entries per row */
   int                   nrows;              /**< complete number of rows */
   SCIP_Real*            lhs;                /**< left hand side per row */
   SCIP_Real*            rhs;                /**< right hand side per row */
   int                   nnonzs;             /**< sparsity counter */
   SCIP_Real*            minactivity;        /**< min activity per row */
   SCIP_Real*            maxactivity;        /**< max activity per row */
   int*                  minactivityneginf;  /**< min activity negative infinity counter */
   int*                  minactivityposinf;  /**< min activity positive infinity counter */
   int*                  maxactivityneginf;  /**< max activity negative infinity counter */
   int*                  maxactivityposinf;  /**< max activity positive infinity counter */
};
typedef struct ConstraintMatrix CONSTRAINTMATRIX;

/** transforms given variables, scalars, and constant to the corresponding active variables, scalars, and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< vars array to get active variables for */
   SCIP_Real**           scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant            /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c  */
   )
{
   int requiredsize;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(scalars != NULL);
   assert(*vars != NULL);
   assert(*scalars != NULL);
   assert(nvars != NULL);
   assert(constant != NULL);

   SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *nvars, constant, &requiredsize, TRUE) );

   if( requiredsize > *nvars )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, vars, requiredsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, scalars, requiredsize) );

      /* call function a second time with enough memory */
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, requiredsize, constant, &requiredsize, TRUE) );
      assert(requiredsize <= *nvars);
   }

   return SCIP_OKAY;
}

/** add one row to the constraint matrix */
static
SCIP_RETCODE addRow(
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix */
   SCIP_VAR**            vars,               /**< variables of this row */
   SCIP_Real*            vals,               /**< coefficients of this row */
   int                   nvars,              /**< number of variables of this row */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   int j;
   int probindex;
   int rowidx;

   assert(vars != NULL);
   assert(vals != NULL);

   rowidx = matrix->nrows;

   matrix->lhs[rowidx] = lhs;
   matrix->rhs[rowidx] = rhs;

   matrix->rowmatbeg[rowidx] = matrix->nnonzs;

   for( j = 0; j < nvars; j++ )
   {
      matrix->rowmatval[matrix->nnonzs] = vals[j];
      probindex = SCIPvarGetProbindex(vars[j]);
      assert(matrix->vars[probindex] == vars[j]);

      assert(0 <= probindex && probindex < matrix->ncols);
      matrix->rowmatind[matrix->nnonzs] = probindex;
      matrix->nnonzs = matrix->nnonzs + 1;
   }

   matrix->rowmatcnt[rowidx] = matrix->nnonzs - matrix->rowmatbeg[rowidx];

   ++(matrix->nrows);

   return SCIP_OKAY;
}

/** add one constraint to matrix */
static
SCIP_RETCODE addConstraint(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix */
   SCIP_VAR**            vars,               /**< variables of this constraint */
   SCIP_Real*            vals,               /**< variable coefficients of this constraint */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs                 /**< right hand side */
   )
{
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   SCIP_Real activeconstant;
   int nactivevars;
   int v;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(vars != NULL || nvars == 0);
   assert(SCIPisLE(scip, lhs, rhs));

   /* constraint is redundant */
   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY;

   /* we do not add empty constraints to the matrix */
   if( nvars == 0 )
      return SCIP_OKAY;

   activevars = NULL;
   activevals = NULL;
   nactivevars = nvars;
   activeconstant = 0.0;

   /* duplicate variable and value array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars ) );
   if( vals != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars ) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; v++ )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, &activevars, &activevals, &nactivevars, &activeconstant) );

   /* adapt left and right hand side */
   if( !SCIPisInfinity(scip, -lhs) )
      lhs -= activeconstant;
   if( !SCIPisInfinity(scip, rhs) )
      rhs -= activeconstant;

   /* add single row to matrix */
   if( nactivevars > 0 )
   {
      SCIP_CALL( addRow(matrix, activevars, activevals, nactivevars, lhs, rhs) );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevals);
   SCIPfreeBufferArray(scip, &activevars);

   return SCIP_OKAY;
}

/** transform row major format into column major format */
static
SCIP_RETCODE setColumnMajorFormat(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix              /**< constraint matrix */
   )
{
   int colidx;
   int i;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   int* fillidx;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(matrix->colmatval != NULL);
   assert(matrix->colmatind != NULL);
   assert(matrix->colmatbeg != NULL);
   assert(matrix->colmatcnt != NULL);
   assert(matrix->rowmatval != NULL);
   assert(matrix->rowmatind != NULL);
   assert(matrix->rowmatbeg != NULL);
   assert(matrix->rowmatcnt != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &fillidx, matrix->ncols) );
   BMSclearMemoryArray(fillidx, matrix->ncols);
   BMSclearMemoryArray(matrix->colmatcnt, matrix->ncols);

   for( i = 0; i < matrix->nrows; i++ )
   {
      rowpnt = matrix->rowmatind + matrix->rowmatbeg[i];
      rowend = rowpnt + matrix->rowmatcnt[i];
      for( ; rowpnt < rowend; rowpnt++ )
      {
         colidx = *rowpnt;
         (matrix->colmatcnt[colidx])++;
      }
   }

   matrix->colmatbeg[0] = 0;
   for( i = 0; i < matrix->ncols-1; i++ )
   {
      matrix->colmatbeg[i+1] = matrix->colmatbeg[i] + matrix->colmatcnt[i];
   }

   for( i = 0; i < matrix->nrows; i++ )
   {
      rowpnt = matrix->rowmatind + matrix->rowmatbeg[i];
      rowend = rowpnt + matrix->rowmatcnt[i];
      valpnt = matrix->rowmatval + matrix->rowmatbeg[i];
      for(; rowpnt < rowend; rowpnt++, valpnt++)
      {
         assert(*rowpnt < matrix->ncols);
         colidx = *rowpnt;
         matrix->colmatval[matrix->colmatbeg[colidx] + fillidx[colidx]] = *valpnt;
         matrix->colmatind[matrix->colmatbeg[colidx] + fillidx[colidx]] = i;
         fillidx[colidx]++;
      }
   }

   SCIPfreeBufferArray(scip, &fillidx);

   return SCIP_OKAY;
}

/** calculate min/max activity per row */
static
SCIP_RETCODE calcActivityBounds(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix              /**< constraint matrix */
   )
{
   SCIP_Real val;
   SCIP_Bool mininfinite;
   SCIP_Bool maxinfinite;
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   int col;
   int row;

   assert(scip != NULL);
   assert(matrix != NULL);

   for( row = 0; row < matrix->nrows; row++ )
   {
      mininfinite = FALSE;
      maxinfinite = FALSE;
      matrix->minactivity[row] = 0;
      matrix->maxactivity[row] = 0;
      matrix->minactivityneginf[row] = 0;
      matrix->minactivityposinf[row] = 0;
      matrix->maxactivityneginf[row] = 0;
      matrix->maxactivityposinf[row] = 0;

      rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
      rowend = rowpnt + matrix->rowmatcnt[row];
      valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

      for(; (rowpnt < rowend) && (!mininfinite || !maxinfinite); rowpnt++, valpnt++)
      {
         /* get column index */
         col = *rowpnt;

         /* get variable coefficient */
         val = *valpnt;

         if( matrix->ncols <= col)
         {
            assert(0);
         }

         if( val >= 0.0 )
         {
            if(SCIPisInfinity(scip, matrix->lb[col]))
               matrix->minactivityposinf[row]++;

            if(SCIPisInfinity(scip, -matrix->lb[col]))
               matrix->minactivityneginf[row]++;

            if(SCIPisInfinity(scip, -matrix->ub[col]))
               matrix->maxactivityneginf[row]++;

            if(SCIPisInfinity(scip, matrix->ub[col]))
               matrix->maxactivityposinf[row]++;

            mininfinite = mininfinite || SCIPisInfinity(scip, -matrix->lb[col]);
            maxinfinite = maxinfinite || SCIPisInfinity(scip, matrix->ub[col]);

            if( !mininfinite )
               matrix->minactivity[row] += val * matrix->lb[col];
            if( !maxinfinite )
               matrix->maxactivity[row] += val * matrix->ub[col];
         }
         else
         {
            if(SCIPisInfinity(scip, -matrix->ub[col]))
               matrix->minactivityposinf[row]++;

            if(SCIPisInfinity(scip, matrix->ub[col]))
               matrix->minactivityneginf[row]++;

            if(SCIPisInfinity(scip, matrix->lb[col]))
               matrix->maxactivityneginf[row]++;

            if(SCIPisInfinity(scip, -matrix->lb[col]))
               matrix->maxactivityposinf[row]++;

            mininfinite = mininfinite || SCIPisInfinity(scip, matrix->ub[col]);
            maxinfinite = maxinfinite || SCIPisInfinity(scip, -matrix->lb[col]);
            if( !mininfinite )
               matrix->minactivity[row] += val * matrix->ub[col];
            if( !maxinfinite )
               matrix->maxactivity[row] += val * matrix->lb[col];
         }
      }

      if( mininfinite )
         matrix->minactivity[row] = -SCIPinfinity(scip);
      if( maxinfinite )
         matrix->maxactivity[row] = SCIPinfinity(scip);
   }

   return SCIP_OKAY;
}

/** initialize matrix */
static
SCIP_RETCODE initMatrix(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX**    matrixptr,          /**< pointer to constraint matrix object to be initialized */
   SCIP_Bool*            initialized         /**< was the initialization successful? */
   )
{
   CONSTRAINTMATRIX* matrix;
   SCIP_CONSHDLR** conshdlrs;
   const char* conshdlrname;
   SCIP_Bool stopped;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   int nconshdlrs;
   int nconss;
   int nnonzstmp;
   int nvars;
   int c;
   int i;
   int v;

   nnonzstmp = 0;

   assert(scip != NULL);
   assert(matrixptr != NULL);
   assert(initialized != NULL);

   *initialized = FALSE;

   /* return if no variables or constraints are present */
   if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
      return SCIP_OKAY;

   /* loop over all constraint handlers and collect the number of checked constraints */
   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);
   nconss = 0;

   for( i = 0; i < nconshdlrs; ++i )
   {
      int nconshdlrconss;

      nconshdlrconss = SCIPconshdlrGetNCheckConss(conshdlrs[i]);

      if( nconshdlrconss > 0 )
      {
         conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);

         if( (strcmp(conshdlrname, "linear") != 0) && (strcmp(conshdlrname, "setppc") != 0)
            && (strcmp(conshdlrname, "logicor") != 0) && (strcmp(conshdlrname, "knapsack") != 0)
            && (strcmp(conshdlrname, "varbound") != 0) )
         {
            SCIPdebugMessage("unsupported constraint type <%s>: aborting domcol presolver\n", conshdlrname);
            break;
         }
      }

      nconss += nconshdlrconss;
   }

   /* do nothing if we have unsupported constraint types or no checked constraints */
   if( i < nconshdlrs || nconss == 0 )
      return SCIP_OKAY;


   stopped = FALSE;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* approximate number of nonzeros by taking for each variable the number of up- and downlocks;
    * this counts nonzeros in equalities twice, but can be at most two times as high as the exact number
    */
   for( i = nvars - 1; i >= 0; --i )
   {
      nnonzstmp += SCIPvarGetNLocksDown(vars[i]);
      nnonzstmp += SCIPvarGetNLocksUp(vars[i]);
   }

   /* do nothing if we have no entries */
   if( nnonzstmp == 0 )
      return SCIP_OKAY;

   /* build the matrix structure */
   SCIP_CALL( SCIPallocBuffer(scip, matrixptr) );
   matrix = *matrixptr;

   /* copy vars array and set number of variables */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &matrix->vars, SCIPgetVars(scip), nvars) );
   matrix->ncols = nvars;

   matrix->nrows = 0;
   matrix->nnonzs = 0;

   /* allocate memory for columns */
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatbeg, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatcnt, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lb, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->ub, matrix->ncols) );

   for( v = 0; v < matrix->ncols; v++ )
   {
      var = matrix->vars[v];
      assert(var != NULL);

      matrix->lb[v] = SCIPvarGetLbGlobal(var);
      matrix->ub[v] = SCIPvarGetUbGlobal(var);
   }

   /* allocate memory for rows */
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatbeg, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatcnt, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lhs, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rhs, nconss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivity, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivity, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivityneginf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivityposinf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivityneginf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivityposinf, nconss) );

   /* loop a second time over constraints handlers and add constraints to the matrix */
   for( i = 0; i < nconshdlrs; ++i )
   {
      SCIP_CONS** conshdlrconss;
      int nconshdlrconss;

      if( SCIPisStopped(scip) )
      {
         stopped = TRUE;
         break;
      }

      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);
      conshdlrconss = SCIPconshdlrGetCheckConss(conshdlrs[i]);
      nconshdlrconss = SCIPconshdlrGetNCheckConss(conshdlrs[i]);

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
         {
            cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsLinear(scip, cons),
                  SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
                  SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons)) );
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
         {
            SCIP_Real lhs;
            SCIP_Real rhs;

            cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            switch( SCIPgetTypeSetppc(scip, cons) )
            {
            case SCIP_SETPPCTYPE_PARTITIONING :
               lhs = 1.0;
               rhs = 1.0;
               break;
            case SCIP_SETPPCTYPE_PACKING :
               lhs = -SCIPinfinity(scip);
               rhs = 1.0;
               break;
            case SCIP_SETPPCTYPE_COVERING :
               lhs = 1.0;
               rhs = SCIPinfinity(scip);
               break;
            default:
               return SCIP_ERROR;
            }

            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsSetppc(scip, cons), NULL,
                  SCIPgetNVarsSetppc(scip, cons), lhs, rhs) );
         }
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
         {
            cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsLogicor(scip, cons),
               NULL, SCIPgetNVarsLogicor(scip, cons), 1.0, SCIPinfinity(scip)) );
         }
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         if( nconshdlrconss > 0 )
         {
            SCIP_Real* consvals;
            int valssize;

            valssize = 100;
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, valssize) );

            for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
            {
               SCIP_Longint* weights;

               cons = conshdlrconss[c];
               assert(SCIPconsIsTransformed(cons));

               weights = SCIPgetWeightsKnapsack(scip, cons);
               nvars = SCIPgetNVarsKnapsack(scip, cons);

               if( nvars > valssize )
               {
                  valssize = (int) (1.5 * nvars);
                  SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, valssize) );
               }

               for( v = 0; v < nvars; v++ )
                  consvals[v] = (SCIP_Real)weights[v];

               SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsKnapsack(scip, cons), consvals,
                     SCIPgetNVarsKnapsack(scip, cons), -SCIPinfinity(scip),
                     (SCIP_Real)SCIPgetCapacityKnapsack(scip, cons)) );
            }

            SCIPfreeBufferArray(scip, &consvals);
         }
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         if( nconshdlrconss > 0 )
         {
            SCIP_VAR** consvars;
            SCIP_Real* consvals;

            SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );
            consvals[0] = 1.0;

            for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
            {
               cons = conshdlrconss[c];
               assert(SCIPconsIsTransformed(cons));

               consvars[0] = SCIPgetVarVarbound(scip, cons);
               consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

               consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

               SCIP_CALL( addConstraint(scip, matrix, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
                     SCIPgetRhsVarbound(scip, cons)) );
            }

            SCIPfreeBufferArray(scip, &consvals);
            SCIPfreeBufferArray(scip, &consvars);
         }
      }
#ifndef NDEBUG
      else
      {
         assert(nconshdlrconss == 0);
      }
#endif
   }
   assert(matrix->nrows <= nconss);
   assert(matrix->nnonzs <= nnonzstmp);

   if( !stopped )
   {
      /* calculate row activity bounds */
      SCIP_CALL( calcActivityBounds(scip, matrix) );

      /* transform row major format into column major format */
      SCIP_CALL( setColumnMajorFormat(scip, matrix) );

      *initialized = TRUE;
   }

   return SCIP_OKAY;
}


/** frees the constraint matrix */
static
void freeMatrix(
   SCIP*                 scip,               /**< current SCIP instance */
   CONSTRAINTMATRIX**    matrix              /**< constraint matrix object */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);

   if( (*matrix) != NULL )
   {
      assert((*matrix)->colmatval != NULL);
      assert((*matrix)->colmatind != NULL);
      assert((*matrix)->colmatbeg != NULL);
      assert((*matrix)->colmatcnt != NULL);
      assert((*matrix)->lb != NULL);
      assert((*matrix)->ub != NULL);

      assert((*matrix)->rowmatval != NULL);
      assert((*matrix)->rowmatind != NULL);
      assert((*matrix)->rowmatbeg != NULL);
      assert((*matrix)->rowmatcnt != NULL);
      assert((*matrix)->lhs != NULL);
      assert((*matrix)->rhs != NULL);

      SCIPfreeBufferArray(scip, &((*matrix)->maxactivityposinf));
      SCIPfreeBufferArray(scip, &((*matrix)->maxactivityneginf));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivityposinf));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivityneginf));
      SCIPfreeBufferArray(scip, &((*matrix)->maxactivity));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivity));

      SCIPfreeBufferArray(scip, &((*matrix)->rhs));
      SCIPfreeBufferArray(scip, &((*matrix)->lhs));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatcnt));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatbeg));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatval));

      SCIPfreeBufferArray(scip, &((*matrix)->ub));
      SCIPfreeBufferArray(scip, &((*matrix)->lb));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatcnt));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatbeg));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->colmatval));

      (*matrix)->nrows = 0;
      (*matrix)->ncols = 0;
      (*matrix)->nnonzs = 0;

      SCIPfreeBufferArrayNull(scip, &((*matrix)->vars));

      SCIPfreeBuffer(scip, matrix);
   }
}

/************** end matrix data structure and functions *************/
/********************************************************************/


/*
 * Local methods
 */

/** get minimum/maximum residual activity without the specified variable */
static
void getActivityResiduals(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   SCIP_VAR*             var,                /**< variable for activity residual calculation */
   SCIP_Real             val,                /**< coefficient of this variable in this row */
   int                   row,                /**< row index */
   SCIP_Real*            minresactivity,     /**< minimum residual activity of this row */
   SCIP_Real*            maxresactivity      /**< maximum residual activity of this row */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(var != NULL);
   assert(row < matrix->nrows);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);

   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, lb) )
      {
         assert(matrix->minactivityposinf[row] >= 1);
         if( matrix->minactivityposinf[row] >= 2 )
            *minresactivity = SCIPinfinity(scip);
         else if( matrix->minactivityneginf[row] >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = matrix->minactivity[row];
      }
      else if( SCIPisInfinity(scip, -lb) )
      {
         assert(matrix->minactivityneginf[row] >= 1);
         if( matrix->minactivityposinf[row] >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( matrix->minactivityneginf[row] >= 2 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = matrix->minactivity[row];
      }
      else
      {
         if( matrix->minactivityposinf[row] >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( matrix->minactivityneginf[row] >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = matrix->minactivity[row] - val * lb;
      }

      if( SCIPisInfinity(scip, -ub) )
      {
         assert(matrix->maxactivityneginf[row] >= 1);
         if( matrix->maxactivityneginf[row] >= 2 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( matrix->maxactivityposinf[row] >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = matrix->maxactivity[row];
      }
      else if( SCIPisInfinity(scip, ub) )
      {
         assert(matrix->maxactivityposinf[row] >= 1);
         if( matrix->maxactivityneginf[row] >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( matrix->maxactivityposinf[row] >= 2 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = matrix->maxactivity[row];
      }
      else
      {
         if( matrix->maxactivityneginf[row] >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( matrix->maxactivityposinf[row] >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = matrix->maxactivity[row] - val * ub;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, -ub) )
      {
         assert(matrix->minactivityposinf[row] >= 1);
         if( matrix->minactivityposinf[row] >= 2 )
            *minresactivity = SCIPinfinity(scip);
         else if( matrix->minactivityneginf[row] >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = matrix->minactivity[row];
      }
      else if( SCIPisInfinity(scip, ub) )
      {
         assert(matrix->minactivityneginf[row] >= 1);
         if( matrix->minactivityposinf[row] >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( matrix->minactivityneginf[row] >= 2 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = matrix->minactivity[row];
      }
      else
      {
         if( matrix->minactivityposinf[row] >= 1 )
            *minresactivity = SCIPinfinity(scip);
         else if( matrix->minactivityneginf[row] >= 1 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = matrix->minactivity[row] - val * ub;
      }

      if( SCIPisInfinity(scip, lb) )
      {
         assert(matrix->maxactivityneginf[row] >= 1);
         if( matrix->maxactivityneginf[row] >= 2 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( matrix->maxactivityposinf[row] >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = matrix->maxactivity[row];
      }
      else if( SCIPisInfinity(scip, -lb) )
      {
         assert(matrix->maxactivityposinf[row] >= 1);
         if( matrix->maxactivityneginf[row] >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( matrix->maxactivityposinf[row] >= 2 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = matrix->maxactivity[row];
      }
      else
      {
         if( matrix->maxactivityneginf[row] >= 1 )
            *maxresactivity = -SCIPinfinity(scip);
         else if( matrix->maxactivityposinf[row] >= 1 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = matrix->maxactivity[row] - val * lb;
      }
   }
}

/** calculate bounds of the dominating variable by rowbound analysis.
 *  we use a special kind of predictive rowbound analysis by first setting the dominated
 *  variable to their lower bound.
 */
static
SCIP_RETCODE calcVarBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   row,                /**< current row index */
   int                   coldominating,      /**< column index of dominating variable */
   SCIP_Real             valdominating,      /**< row coefficient of dominating variable */
   int                   coldominated,       /**< column index of dominated variable */
   SCIP_Real             valdominated,       /**< row coefficient of dominated variable */
   SCIP_Bool*            ubcalculated,       /**< was a upper bound calculated? */
   SCIP_Real*            calculatedub,       /**< predicted upper bound */
   SCIP_Bool*            wclbcalculated,     /**< was a lower worst case bound calculated? */
   SCIP_Real*            calculatedwclb      /**< predicted worst case lower bound */
   )
{
   SCIP_VAR* vardominating;
   SCIP_VAR* vardominated;
   SCIP_Real lbdominated;
   SCIP_Real ubdominated;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);
   assert(0 <= coldominating && coldominating < matrix->ncols);
   assert(0 <= coldominated && coldominated < matrix->ncols);

   assert(ubcalculated != NULL);
   assert(calculatedub != NULL);
   assert(wclbcalculated != NULL);
   assert(calculatedwclb != NULL);

   assert(!SCIPisZero(scip, valdominating));
   assert(matrix->vars[coldominating] != NULL);

   vardominating = matrix->vars[coldominating];
   vardominated = matrix->vars[coldominated];

   *ubcalculated = FALSE;
   *wclbcalculated = FALSE;

   /* no rowbound analysis for multiaggregated variables */
   if( SCIPvarGetStatus(vardominating) == SCIP_VARSTATUS_MULTAGGR ||
      SCIPvarGetStatus(vardominated) == SCIP_VARSTATUS_MULTAGGR )
   {
      return SCIP_OKAY;
   }

   lhs = matrix->lhs[row];
   rhs = matrix->rhs[row];
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));

   /* get minimum/maximum activity of this row without the dominating variable */
   getActivityResiduals(scip, matrix, vardominating, valdominating, row, &minresactivity, &maxresactivity);

   assert(!SCIPisInfinity(scip, minresactivity));
   assert(!SCIPisInfinity(scip, -maxresactivity));

   lbdominated = SCIPvarGetLbLocal(vardominated);
   ubdominated = SCIPvarGetUbLocal(vardominated);

   *calculatedub = SCIPinfinity(scip);
   *calculatedwclb = -SCIPinfinity(scip);

   /* predictive rowbound analysis */

   if( valdominating > 0.0 )
   {
      /* upper bound calculation */
      if( !SCIPisInfinity(scip, -minresactivity) && !SCIPisInfinity(scip, rhs) )
      {
         /* <= */
         *ubcalculated = TRUE;
         if( SCIPisZero(scip, valdominated) )
            *calculatedub = (rhs - minresactivity)/valdominating;
         else
            if( valdominated < 0.0 )
               *calculatedub = (rhs - (minresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
            else
               *calculatedub = (rhs - (minresactivity /* - (valdominated * lbdominated) + (valdominated * lbdominated)*/))/valdominating;
      }

      /* worst case calculation for lower bound */
      if( !SCIPisInfinity(scip, -minresactivity) && !SCIPisInfinity(scip, -lhs) )
      {
         /* >= */
         *wclbcalculated = TRUE;
         if( SCIPisZero(scip, valdominated) )
            *calculatedwclb = (lhs - minresactivity)/valdominating;
         else
            if( valdominated < 0.0 )
               *calculatedwclb = (lhs - (minresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
            else
               *calculatedwclb = (lhs - minresactivity)/valdominating;
      }
   }
   else
   {
      /* worst case calculation for lower bound */
      if( !SCIPisInfinity(scip, maxresactivity) && !SCIPisInfinity(scip, rhs) )
      {
         /* <= */
         *wclbcalculated = TRUE;
         if( SCIPisZero(scip, valdominated) )
            *calculatedwclb = (rhs - maxresactivity)/valdominating;
         else
            if( valdominated < 0.0 )
               *calculatedwclb = (rhs - maxresactivity)/valdominating;
            else
               *calculatedwclb = (rhs - (maxresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
      }

      /* upper bound calculation */
      if( !SCIPisInfinity(scip, maxresactivity) && !SCIPisInfinity(scip, -lhs) )
      {
         /* >= */
         *ubcalculated = TRUE;
         if( SCIPisZero(scip, valdominated) )
            *calculatedub = (lhs - maxresactivity)/valdominating;
         else
            if( valdominated < 0.0 )
               *calculatedub = (lhs - (maxresactivity /* - (valdominated * lbdominated) + (valdominated * lbdominated)*/))/valdominating;
            else
               *calculatedub = (lhs - (maxresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
      }
   }

   return SCIP_OKAY;
}

/** try to find new variable bounds and update them whether they are better then the old bounds */
static
SCIP_RETCODE updateBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   row,                /**< row index */
   int                   col1,               /**< dominating variable index */
   SCIP_Real             val1,               /**< dominating variable coefficient */
   int                   col2,               /**< dominated variable index */
   SCIP_Real             val2,               /**< dominated variable coefficient */
   SCIP_Real*            upperbound,         /**< predicted upper bound */
   SCIP_Real*            wclowerbound        /**< predicted worst case lower bound */
   )
{
   SCIP_Bool ubcalculated;
   SCIP_Bool wclbcalculated;
   SCIP_Real newub;
   SCIP_Real newwclb;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(row < matrix->nrows);
   assert(col1 < matrix->ncols);
   assert(col2 < matrix->ncols);

   /* do predictive rowbound analysis */
   SCIP_CALL( calcVarBounds(scip, matrix, row, col1, val1, col2, val2,
         &ubcalculated, &newub, &wclbcalculated, &newwclb) );

   /* update bounds in case if they are better */
   if( ubcalculated )
   {
      if( newub < *upperbound )
         *upperbound = newub;
   }
   if( wclbcalculated )
   {
      if( newwclb > *wclowerbound )
         *wclowerbound = newwclb;
   }

   return SCIP_OKAY;
}

/** try to find possible variable fixings */
static
void findFixings(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_VAR*             dominatingvar,      /**< dominating variable */
   int                   dominatingidx,      /**< column index of the dominating variable */
   SCIP_Real             dominatingub,       /**< predicted upper bound of the dominating variable */
   SCIP_Real             dominatingwclb,     /**< predicted worst case lower bound of the dominating variable */
   SCIP_VAR*             dominatedvar,       /**< dominated variable */
   int                   dominatedidx,       /**< column index of the dominated variable */
   FIXINGDIRECTION*      varstofix,          /**< array holding fixing information */
   SCIP_Bool             onlybinvars,        /**< flag indicating only binary variables are present */
   int*                  npossiblefixings,   /**< counter for possible fixings */
   int*                  ncliquepreventions, /**< counter for clique preventions */
   int*                  nboundpreventions   /**> counter for bound preventions */
   )
{
   if( onlybinvars )
   {
      if( SCIPvarsHaveCommonClique(dominatingvar, TRUE, dominatedvar, TRUE, TRUE) &&
         (!SCIPvarsHaveCommonClique(dominatingvar, TRUE, dominatedvar, FALSE, TRUE) ||
            !SCIPvarsHaveCommonClique(dominatingvar, FALSE, dominatedvar, FALSE, TRUE) ) )
      {
         /* we have a (1->1)-clique with dominance relation (x->y) (x dominates y)
          * from dominance relation we know (1->0) is better then (0->1)
          * it follows, only (1->0) or (0->0) are possible => y=0
          */
         if( varstofix[dominatedidx] == NOFIX )
         {
            varstofix[dominatedidx] = FIXATLB;
            (*npossiblefixings)++;
         }
      }
      else if( SCIPvarsHaveCommonClique(dominatingvar, FALSE, dominatedvar, FALSE, TRUE) &&
         (!SCIPvarsHaveCommonClique(dominatingvar, TRUE, dominatedvar, TRUE, TRUE) ||
            !SCIPvarsHaveCommonClique(dominatingvar, TRUE, dominatedvar, FALSE, TRUE) ) )
      {
         /* we have a (0->0)-clique with dominance relation x->y (x dominates y)
          * from dominance relation we know (1->0) is better then (0->1)
          * it follows only (1->0) or (1->1) are possible => x=1
          */
         if( varstofix[dominatingidx] == NOFIX )
         {
            varstofix[dominatingidx] = FIXATUB;
            (*npossiblefixings)++;
         }
      }
      else
      {
         (*ncliquepreventions)++;
      }
   }
   else
   {
      if( SCIPisPositive(scip, SCIPvarGetObj(dominatingvar)) )
      {
         assert(SCIPisPositive(scip, SCIPvarGetObj(dominatedvar)));
         if( !SCIPisInfinity(scip, -dominatingwclb) &&
            SCIPisLE(scip, dominatingwclb, SCIPvarGetUbLocal(dominatingvar)) )
         {
            /* we have a x->y dominance relation with a positive obj coefficient
             * of the dominating variable x, thus the obj coefficient of the
             * dominated variable y is positive too. we need to secure feasibility
             * by testing if the predicted lower worst case bound is less equal the current upper bound.
             */
            if( varstofix[dominatedidx] == NOFIX )
            {
               varstofix[dominatedidx] = FIXATLB;
               (*npossiblefixings)++;
            }
         }
         else
         {
            (*nboundpreventions)++;
         }
      }
      else /* SCIPvarGetObj(dominatingvar) <= 0.0 */
      {
         if( !SCIPisInfinity(scip, dominatingub) &&
            SCIPisLE(scip, dominatingub, SCIPvarGetUbLocal(dominatingvar)) )
         {
            /* we have a x->y dominance relation with a negative or zero obj coefficient
             * of the dominating variable x, thus the obj coefficient of the
             * dominated variable y is positive, negative or zero. in all cases we have to look
             * if the predicted upper bound of the dominating variable is great enough.
             */
            if( varstofix[dominatedidx] == NOFIX )
            {
               varstofix[dominatedidx] = FIXATLB;
               (*npossiblefixings)++;
            }
         }
         else
         {
            (*nboundpreventions)++;
         }
      }
   }
}

/** find dominance relation between variable pairs */
static
SCIP_RETCODE findDominancePairs(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int*                  searchcols,         /**< indexes of variables for pair comparisons */
   int                   searchsize,         /**< number of variables for pair comparisons */
   SCIP_Bool             onlybinvars,        /**< flag indicating searchcols has only binary variable indexes */
   FIXINGDIRECTION*      varstofix,          /**< array holding information for later upper/lower bound fixing */
   int*                  npossiblefixings,   /**< found number of possible fixings */
   int*                  ndomrelations,      /**< found number of dominance relations */
   int*                  ncliquepreventions, /**< number of clique preventions for doing a variable fixing */
   int*                  nboundpreventions   /**< number of bound preventions for doing a variable fixing */
   )
{
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   SCIP_Real tmpupperboundcol1;
   SCIP_Real tmpupperboundcol2;
   SCIP_Real tmpwclowerboundcol1;
   SCIP_Real tmpwclowerboundcol2;
   int* rows1;
   int* rows2;
   int nrows1;
   int nrows2;
   SCIP_Bool col1domcol2;
   SCIP_Bool col2domcol1;
   int cnt1;
   int cnt2;
   int col1;
   int col2;
   int r1;
   int r2;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(searchcols != NULL);
   assert(varstofix != NULL);
   assert(npossiblefixings != NULL);
   assert(ndomrelations != NULL);
   assert(ncliquepreventions != NULL);
   assert(nboundpreventions != NULL);

   /* pair comparisons */
   for( cnt1 = 0; cnt1 < searchsize; cnt1++ )
   {
      for( cnt2 = cnt1+1; cnt2 < searchsize; cnt2++ )
      {
         /* get indexes of this variable pair */
         col1 = searchcols[cnt1];
         col2 = searchcols[cnt2];

         /* we always have minimize as obj sense */

         /* column 1 dominating column 2 */
         col1domcol2 = (SCIPvarGetObj(matrix->vars[col1]) <= SCIPvarGetObj(matrix->vars[col2]));

         /* column 2 dominating column 1 */
         col2domcol1 = (SCIPvarGetObj(matrix->vars[col2]) <= SCIPvarGetObj(matrix->vars[col1]));

         /* search only if nothing was found yet */
         col1domcol2 = col1domcol2 && (varstofix[col2] == NOFIX);
         col2domcol1 = col2domcol1 && (varstofix[col1] == NOFIX);

         /* we search only for a dominance relation if the lower bounds are not negative */
         if( !onlybinvars )
         {
            if( SCIPisLT(scip, SCIPvarGetLbLocal(matrix->vars[col1]), 0.0) ||
               SCIPisLT(scip, SCIPvarGetLbLocal(matrix->vars[col2]), 0.0) )
            {
               col1domcol2 = FALSE;
               col2domcol1 = FALSE;
            }
         }

         if( !col1domcol2 && !col2domcol1 )
            continue;

         /* get the data for both columns */
         vals1 = matrix->colmatval + matrix->colmatbeg[col1];
         rows1 = matrix->colmatind + matrix->colmatbeg[col1];
         nrows1 = matrix->colmatcnt[col1];
         r1 = 0;
         vals2 = matrix->colmatval + matrix->colmatbeg[col2];
         rows2 = matrix->colmatind + matrix->colmatbeg[col2];
         nrows2 = matrix->colmatcnt[col2];
         r2 = 0;

         /* do we have a obj constant? */
         if( nrows1 == 0 || nrows2 == 0 )
         {
            col1domcol2 = FALSE;
            col2domcol1 = FALSE;
            continue;
         }

         /* initialize temporary bounds */
         tmpupperboundcol1 = SCIPinfinity(scip);
         tmpupperboundcol2 = tmpupperboundcol1;
         tmpwclowerboundcol1 = -SCIPinfinity(scip);
         tmpwclowerboundcol2 = tmpwclowerboundcol1;

         /* compare rows of this column pair */
         while( (col1domcol2 || col2domcol1) && (r1 < nrows1 || r2 < nrows2))
         {
            assert((r1 >= nrows1-1) || (rows1[r1] < rows1[r1+1]));
            assert((r2 >= nrows2-1) || (rows2[r2] < rows2[r2+1]));

            /* there is a nonredundant row containing column 1 but not column 2 */
            if( r1 < nrows1 && (r2 == nrows2 || rows1[r1] < rows2[r2]) )
            {
               /* dominance depends on the type of relation */
               if( !SCIPisInfinity(scip, -matrix->lhs[rows1[r1]]) &&
                  !SCIPisInfinity(scip, matrix->rhs[rows1[r1]]) )
               {
                  /* no dominance relation for equations or ranged rows */
                  col2domcol1 = FALSE;
                  col1domcol2 = FALSE;
               }
               else if( SCIPisInfinity(scip, -matrix->lhs[rows1[r1]]) &&
                  !SCIPisInfinity(scip, matrix->rhs[rows1[r1]]) )
               {
                  /* <= relation, smaller coefficients dominate larger ones */
                  if( vals1[r1] < 0.0 )
                     col2domcol1 = FALSE;
                  else if( vals1[r1] > 0.0 )
                     col1domcol2 = FALSE;
               }
               else if( !SCIPisInfinity(scip, -matrix->lhs[rows1[r1]]) &&
                  SCIPisInfinity(scip, matrix->rhs[rows1[r1]]) )
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals1[r1] > 0.0 )
                     col2domcol1 = FALSE;
                  else if( vals1[r1] < 0.0 )
                     col1domcol2 = FALSE;
               }
               else
               {
                  col2domcol1 = FALSE;
                  col1domcol2 = FALSE;
               }

               /* update bounds only for dominance relation of non binary variables */
               if( col1domcol2 && !onlybinvars )
               {
                  SCIP_CALL( updateBounds(scip, matrix, rows1[r1], col1, vals1[r1], col2, 0.0,
                        &tmpupperboundcol1, &tmpwclowerboundcol1) );
               }

               r1++;
            }
            /* there is a nonredundant row containing column 2, but not column 1 */
            else if( r2 < nrows2 && (r1 == nrows1 || rows1[r1] > rows2[r2]) )
            {
               /* dominance depends on the type of relation */
               if( !SCIPisInfinity(scip, -matrix->lhs[rows2[r2]]) &&
                  !SCIPisInfinity(scip, matrix->rhs[rows2[r2]]) )
               {
                  /* no dominance relation for equations or ranged rows */
                  col2domcol1 = FALSE;
                  col1domcol2 = FALSE;
               }
               else if( SCIPisInfinity(scip, -matrix->lhs[rows2[r2]]) &&
                  !SCIPisInfinity(scip, matrix->rhs[rows2[r2]]) )
               {
                  /* <= relation, smaller coefficients dominate larger ones */
                  if( vals2[r2] > 0.0 )
                     col2domcol1 = FALSE;
                  else if( vals2[r2] < 0.0 )
                     col1domcol2 = FALSE;
               }
               else if( !SCIPisInfinity(scip, -matrix->lhs[rows2[r2]]) &&
                  SCIPisInfinity(scip, matrix->rhs[rows2[r2]]) )
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals2[r2] < 0.0 )
                     col2domcol1 = FALSE;
                  else if( vals2[r2] > 0.0 )
                     col1domcol2 = FALSE;
               }
               else
               {
                  col2domcol1 = FALSE;
                  col1domcol2 = FALSE;
               }

               /* update bounds only for dominance relation of non binary variables */
               if( col2domcol1 && !onlybinvars )
               {
                  SCIP_CALL( updateBounds(scip, matrix, rows2[r2], col2, vals2[r2], col1, 0.0,
                        &tmpupperboundcol2, &tmpwclowerboundcol2) );
               }
               r2++;
            }
            /* if both columns appear in a common row, compare the coefficients */
            else
            {
               assert(r1 < nrows1 && r2 < nrows2);
               assert(rows1[r1] == rows2[r2]);

               /* dominance depends on the type of inequality */
               if( !SCIPisInfinity(scip, -matrix->lhs[rows1[r1]]) &&
                  !SCIPisInfinity(scip, matrix->rhs[rows1[r1]]) )
               {
                  /* no dominance relation if coefficients differ for equations or ranged rows */
                  if( !SCIPisEQ(scip, vals1[r1], vals2[r2]) )
                  {
                     col2domcol1 = FALSE;
                     col1domcol2 = FALSE;
                  }
               }
               else if(  SCIPisInfinity(scip, -matrix->lhs[rows1[r1]]) &&
                  !SCIPisInfinity(scip, matrix->rhs[rows1[r1]]) )
               {
                  /* <= relation, smaller coefficients dominate larger ones */
                  if( vals1[r1] < vals2[r2] )
                     col2domcol1 = FALSE;
                  else if( vals1[r1] > vals2[r2] )
                     col1domcol2 = FALSE;
               }
               else if( !SCIPisInfinity(scip, -matrix->lhs[rows1[r1]]) &&
                  SCIPisInfinity(scip, matrix->rhs[rows1[r1]]) )
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals1[r1] > vals2[r2] )
                     col2domcol1 = FALSE;
                  else if( vals1[r1] < vals2[r2] )
                     col1domcol2 = FALSE;
               }
               else
               {
                  col1domcol2 = FALSE;
                  col2domcol1 = FALSE;
               }

               /* update bounds only for dominance relation of non binary variables */
               if( col1domcol2 && !onlybinvars )
               {
                  SCIP_CALL( updateBounds(scip, matrix, rows1[r1], col1, vals1[r1], col2, vals2[r2],
                        &tmpupperboundcol1, &tmpwclowerboundcol1) );
               }
               else if( col2domcol1 && !onlybinvars )
               {
                  SCIP_CALL( updateBounds(scip, matrix, rows2[r2], col2, vals2[r2], col1, vals1[r1],
                        &tmpupperboundcol2, &tmpwclowerboundcol2) );
               }

               r1++;
               r2++;
            }
         }

         /* a column is only dominated, if there are no more rows in which it is contained */
         col1domcol2 = col1domcol2 && r2 == nrows2;
         col2domcol1 = col2domcol1 && r1 == nrows1;

         if( !col1domcol2 && !col2domcol1 )
            continue;

         /* no dominance relation for left equations or ranged rows */
         while( r1 < nrows1 )
         {
            if( !SCIPisInfinity(scip, -matrix->lhs[rows1[r1]]) &&
               !SCIPisInfinity(scip, matrix->rhs[rows1[r1]]) )
            {
               col2domcol1 = FALSE;
               col1domcol2 = FALSE;
               break;
            }
            r1++;
         }
         while( r2 < nrows2 )
         {
            if( !SCIPisInfinity(scip, -matrix->lhs[rows2[r2]]) &&
               !SCIPisInfinity(scip, matrix->rhs[rows2[r2]]) )
            {
               col2domcol1 = FALSE;
               col1domcol2 = FALSE;
               break;
            }
            r2++;
         }

         if( col1domcol2 || col2domcol1 )
            (*ndomrelations)++;

         /* use dominance relation and cliquen/bound-information to find variable fixings */
         if( col1domcol2 )
         {
            findFixings(scip, matrix->vars[col1], col1, tmpupperboundcol1, tmpwclowerboundcol1, matrix->vars[col2],
               col2, varstofix, onlybinvars, npossiblefixings, ncliquepreventions, nboundpreventions);
         }
         else if( col2domcol1 )
         {
            findFixings(scip, matrix->vars[col2], col2, tmpupperboundcol2, tmpwclowerboundcol2, matrix->vars[col1],
               col1, varstofix, onlybinvars, npossiblefixings, ncliquepreventions, nboundpreventions);
         }
      }
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDomcol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDomcol(scip) );

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDomcol)
{  /*lint --e{715}*/
   CONSTRAINTMATRIX* matrix;
   SCIP_Bool initialized;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   /* do no dominated column presolving in case of probing and nonlinear processing
    * @todo SCIPisNLPEnabled() always returns FALSE during presolve, since the necessary flag is set after presolve (in exitpre, currently)
    */
   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* initialize constraint matrix */
   matrix = NULL;
   initialized = FALSE;
   SCIP_CALL( initMatrix(scip, &matrix, &initialized) );

   if( initialized )
   {
      int npossiblefixings;
      int ncliquepreventions;
      int nboundpreventions;
      int ndomrelations;
      int v;
      int r;
      FIXINGDIRECTION* varstofix;
      int* varsprocessed;
      int nconvarsfixed;
      int nintvarsfixed;
      int nbinvarsfixed;
      int nvars;
      int nrows;
      int* rowidxsorted;
      int* rowsparsity;
      int varcount;
      int* consearchcols;
      int* intsearchcols;
      int* binsearchcols;
      int nconfill;
      int nintfill;
      int nbinfill;
      SCIP_Bool onlybinvars;

      assert(SCIPgetNVars(scip) == matrix->ncols);

      npossiblefixings = 0;
      ncliquepreventions = 0;
      nboundpreventions = 0;
      ndomrelations = 0;
      nconvarsfixed = 0;
      nintvarsfixed = 0;
      nbinvarsfixed = 0;
      nvars = matrix->ncols;
      nrows = matrix->nrows;

      SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varsprocessed, nvars) );
      BMSclearMemoryArray(varstofix, nvars);
      BMSclearMemoryArray(varsprocessed, nvars);

      SCIP_CALL( SCIPallocBufferArray(scip, &consearchcols, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &intsearchcols, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &binsearchcols, nvars) );

      SCIP_CALL( SCIPallocBufferArray(scip, &rowidxsorted, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowsparsity, nrows) );
      for( r = 0; r < nrows; ++r )
      {
         rowidxsorted[r] = r;
         rowsparsity[r] = matrix->rowmatcnt[r];
      }

      /* sort rows per sparsity */
      SCIPsortIntInt(rowsparsity, rowidxsorted, nrows);

      /* verify if we only have binary variables */
      onlybinvars = (SCIPgetNIntVars(scip) + SCIPgetNImplVars(scip) + SCIPgetNContVars(scip)) == 0;

      /* search for dominance relations by row-sparsity order */
      varcount = 0;
      for( r = 0; r < nrows; ++r )
      {
         int rowidx;
         int* rowpnt;
         int* rowend;

         /* break if the time limit was reached; since the check is expensive, we only check all 1000 constraints */
         if( (r % 1000 == 0) && SCIPisStopped(scip) )
            break;

         rowidx = rowidxsorted[r];
         rowpnt = matrix->rowmatind + matrix->rowmatbeg[rowidx];
         rowend = rowpnt + matrix->rowmatcnt[rowidx];

         nconfill = 0;
         nintfill = 0;
         nbinfill = 0;

         if( !onlybinvars )
         {
            for( ; rowpnt < rowend; rowpnt++ )
            {
               if( varsprocessed[*rowpnt] == FALSE )
               {
                  int varidx;
                  varidx = *rowpnt;

                  /* higher variable types dominate smaller ones: bin <- int <- impl <- cont
                   * we search only for dominance relations between the same variable type
                   */
                  if( SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_CONTINUOUS )
                  {
                     consearchcols[nconfill++] = varidx;
                  }
                  else if( SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_INTEGER ||
                     SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_IMPLINT )
                  {
                     intsearchcols[nintfill++] = varidx;
                  }
                  else if( SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_BINARY )
                  {
                     binsearchcols[nbinfill++] = varidx;
                  }
               }
            }

            /* search for dominance relations between continuous variables */
            if( nconfill > 1 )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, consearchcols, nconfill, FALSE,
                     varstofix, &npossiblefixings, &ndomrelations, &ncliquepreventions, &nboundpreventions) );
            }
            for( v = 0; v < nconfill; ++v )
            {
               varsprocessed[consearchcols[v]] = TRUE;
            }
            varcount += nconfill;

            /* search for dominance relations between integer and impl-integer variables */
            if( nintfill > 1 )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, intsearchcols, nintfill, FALSE,
                     varstofix, &npossiblefixings, &ndomrelations, &ncliquepreventions, &nboundpreventions) );
            }
            for( v = 0; v < nintfill; ++v )
            {
               varsprocessed[intsearchcols[v]] = TRUE;
            }
            varcount += nintfill;

            /* search for dominance relations between binary variables */
            if( nbinfill > 1 )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, binsearchcols, nbinfill, TRUE,
                     varstofix, &npossiblefixings, &ndomrelations, &ncliquepreventions, &nboundpreventions) );
            }
            for( v = 0; v < nbinfill; ++v )
            {
               varsprocessed[binsearchcols[v]] = TRUE;
            }
            varcount += nbinfill;
         }
         else
         {
            /* we only have binary variables */
            for( ; rowpnt < rowend; rowpnt++ )
            {
               if( varsprocessed[*rowpnt] == FALSE )
               {
                  binsearchcols[nbinfill++] = *rowpnt;
               }
            }

            /* search for dominance relations between binary variables */
            if( nbinfill > 1 )
            {
               SCIP_CALL( findDominancePairs(scip, matrix, binsearchcols, nbinfill, TRUE,
                     varstofix, &npossiblefixings, &ndomrelations, &ncliquepreventions, &nboundpreventions) );
            }
            for( v = 0; v < nbinfill; ++v )
            {
               varsprocessed[binsearchcols[v]] = TRUE;
            }
            varcount += nbinfill;
         }

         if( varcount >= nvars )
         {
            break;
         }
      }

      if( npossiblefixings > 0 )
      {
         /* look for fixable variables */
         for( v = matrix->ncols-1; v >= 0; --v )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            SCIP_VAR* var;

            if( varstofix[v] == FIXATLB )
            {
               SCIP_Real lb;

               var = matrix->vars[v];
               lb = SCIPvarGetLbLocal(var);

               if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
               {
                  lb = SCIPfeasCeil(scip, lb);
               }

               /* fix at lower bound */
               SCIP_CALL( SCIPfixVar(scip, var, lb, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage(" -> infeasible fixing\n");
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               assert(fixed);
               (*nfixedvars)++;
               *result = SCIP_SUCCESS;

               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               {
                  nconvarsfixed++;
               }
               else if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
               {
                  nintvarsfixed++;
               }
               else
               {
                  nbinvarsfixed++;
               }
            }
            else if( varstofix[v] == FIXATUB )
            {
               SCIP_Real ub;

               var = matrix->vars[v];
               ub = SCIPvarGetUbLocal(var);

               if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
               {
                  ub = SCIPfeasFloor(scip, ub);
               }

               /* fix at upper bound */
               SCIP_CALL( SCIPfixVar(scip, var, ub, &infeasible, &fixed) );
               if( infeasible )
               {
                  SCIPdebugMessage(" -> infeasible fixing\n");
                  *result = SCIP_CUTOFF;
                  return SCIP_OKAY;
               }
               assert(fixed);
               (*nfixedvars)++;
               *result = SCIP_SUCCESS;

               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               {
                  nconvarsfixed++;
               }
               else if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
               {
                  nintvarsfixed++;
               }
               else
               {
                  nbinvarsfixed++;
               }
            }
         }
      }

      SCIPfreeBufferArray(scip, &rowsparsity);
      SCIPfreeBufferArray(scip, &rowidxsorted);
      SCIPfreeBufferArray(scip, &binsearchcols);
      SCIPfreeBufferArray(scip, &intsearchcols);
      SCIPfreeBufferArray(scip, &consearchcols);
      SCIPfreeBufferArray(scip, &varsprocessed);
      SCIPfreeBufferArray(scip, &varstofix);

      if( (nconvarsfixed + nintvarsfixed + nbinvarsfixed) > 0 )
      {
         SCIPdebugMessage("### %d vars [%d dom] (%d clqprev, %d bndprev) ===>>> fixed [c=%d, z=%d, b=%d]\n",
            matrix->ncols, ndomrelations, ncliquepreventions, nboundpreventions, nconvarsfixed, nintvarsfixed, nbinvarsfixed);
      }
      else
      {
          SCIPdebugMessage("### %d vars [%d dom] (%d clqprev, %d bndprev)\n",
            matrix->ncols, ndomrelations, ncliquepreventions, nboundpreventions);
      }
   }

   freeMatrix(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the domcol presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDomcol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_DELAY, presolExecDomcol, NULL) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDomcol) );

   return SCIP_OKAY;
}
