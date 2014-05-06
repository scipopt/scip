/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_implfree.c
 * @brief  exploit implied free variables for multi-aggregation
 * @author Dieter Weninger
 *
 * This presolver tries to find implied free variables within equalities and
 * exploits this variables for multi-aggregation.
 *
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

#include "presol_implfree.h"

#define PRESOL_NAME            "implfree"
#define PRESOL_DESC            "exploit implied free variables for multi-aggregation"
#define PRESOL_PRIORITY         12000000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY                TRUE     /**< should presolver be delayed, if other presolvers found reductions? */



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
   int*                  nuplocks;           /**< number of up locks per variable */
   int*                  ndownlocks;         /**< number of down locks per variable */

   SCIP_VAR**            vars;               /**< variables pointer */

   SCIP_Real*            rowmatval;          /**< coefficients in row major format */
   int*                  rowmatind;          /**< column indexed in row major format */
   int*                  rowmatbeg;          /**< row storage offset */
   int*                  rowmatcnt;          /**< number of column entries per row */

#ifdef SCIP_DEBUG
   const char**          rowname;            /**< name of row */
#endif
   int                   nrows;              /**< complete number of rows */
   SCIP_Real*            lhs;                /**< left hand side per row */
   SCIP_Real*            rhs;                /**< right hand side per row */

   SCIP_CONS**           cons;               /**< constraints pointer */

   SCIP_Bool*            isrhsinfinite;      /**< is right hand side infinity */
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
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix */
#ifdef SCIP_DEBUG
   const char*           name,               /**< name of constraint corresponding to row */
#endif
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
   SCIP_Real factor;
   SCIP_Bool rangedorequality;

   assert(vars != NULL);
   assert(vals != NULL);

   rowidx = matrix->nrows;
   rangedorequality = FALSE;

   if( SCIPisInfinity(scip, -lhs) )
   {
      factor = -1.0;
      matrix->lhs[rowidx] = -rhs;
      matrix->rhs[rowidx] = SCIPinfinity(scip);
      matrix->isrhsinfinite[rowidx] = TRUE;
   }
   else
   {
      factor = 1.0;
      matrix->lhs[rowidx] = lhs;
      matrix->rhs[rowidx] = rhs;
      matrix->isrhsinfinite[rowidx] = SCIPisInfinity(scip, matrix->rhs[rowidx]);

      if( !SCIPisInfinity(scip, rhs) )
         rangedorequality = TRUE;
   }
   assert(!SCIPisInfinity(scip, -matrix->lhs[rowidx]));


#ifdef SCIP_DEBUG
   matrix->rowname[rowidx] = name;
#endif

   matrix->rowmatbeg[rowidx] = matrix->nnonzs;

   /* = or ranged */
   if( rangedorequality )
   {
      assert(factor > 0);

      for( j = 0; j < nvars; j++ )
      {
         matrix->rowmatval[matrix->nnonzs] = factor * vals[j];
         probindex = SCIPvarGetProbindex(vars[j]);
         assert(matrix->vars[probindex] == vars[j]);

         matrix->nuplocks[probindex]++;
         matrix->ndownlocks[probindex]++;

         assert(0 <= probindex && probindex < matrix->ncols);
         matrix->rowmatind[matrix->nnonzs] = probindex;
         (matrix->nnonzs)++;
      }
   }
   /* >= or <= */
   else
   {
      for( j = 0; j < nvars; j++ )
      {
         /* due to the factor, <= constraints will be transfered to >= */
         matrix->rowmatval[matrix->nnonzs] = factor * vals[j];
         probindex = SCIPvarGetProbindex(vars[j]);
         assert(matrix->vars[probindex] == vars[j]);

         if( matrix->rowmatval[matrix->nnonzs] > 0 )
            matrix->ndownlocks[probindex]++;
         else
         {
            assert(matrix->rowmatval[matrix->nnonzs] < 0);
            matrix->nuplocks[probindex]++;
         }

         assert(0 <= probindex && probindex < matrix->ncols);
         matrix->rowmatind[matrix->nnonzs] = probindex;
         (matrix->nnonzs)++;
      }
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
#ifdef SCIP_DEBUG
   const char*           name,               /**< name of constraint corresponding to row */
#endif
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
#ifdef SCIP_DEBUG
      SCIP_CALL( addRow(scip, matrix, name, activevars, activevals, nactivevars, lhs, rhs) );
#else
      SCIP_CALL( addRow(scip, matrix, activevars, activevals, nactivevars, lhs, rhs) );
#endif
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

      for( ; rowpnt < rowend; rowpnt++, valpnt++ )
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
   int* rowpnt;
   int* rowend;
   SCIP_Real* valpnt;
   int col;
   int row;

   assert(scip != NULL);
   assert(matrix != NULL);

   for( row = 0; row < matrix->nrows; row++ )
   {
      matrix->minactivity[row] = 0;
      matrix->maxactivity[row] = 0;
      matrix->minactivityneginf[row] = 0;
      matrix->minactivityposinf[row] = 0;
      matrix->maxactivityneginf[row] = 0;
      matrix->maxactivityposinf[row] = 0;

      rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
      rowend = rowpnt + matrix->rowmatcnt[row];
      valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

      for( ; rowpnt < rowend; rowpnt++, valpnt++ )
      {
         /* get column index */
         col = *rowpnt;

         /* get variable coefficient */
         val = *valpnt;
         assert(!SCIPisZero(scip, val));

         assert(matrix->ncols > col);

         assert(!SCIPisInfinity(scip, matrix->lb[col]));
         assert(!SCIPisInfinity(scip, -matrix->ub[col]));

         /* positive coefficient */
         if( val > 0.0 )
         {
            if( SCIPisInfinity(scip, matrix->ub[col]) )
               matrix->maxactivityposinf[row]++;
            else
               matrix->maxactivity[row] += val * matrix->ub[col];

            if( SCIPisInfinity(scip, -matrix->lb[col]) )
               matrix->minactivityneginf[row]++;
            else
               matrix->minactivity[row] += val * matrix->lb[col];
         }
         /* negative coefficient */
         else
         {
            if( SCIPisInfinity(scip, -matrix->lb[col]) )
               matrix->maxactivityneginf[row]++;
            else
               matrix->maxactivity[row] += val * matrix->lb[col];

            if( SCIPisInfinity(scip, matrix->ub[col]) )
               matrix->minactivityposinf[row]++;
            else
               matrix->minactivity[row] += val * matrix->ub[col];
         }
      }
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
   int cnt;

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
            SCIPdebugMessage("unsupported constraint type <%s>: aborting implfree presolver\n", conshdlrname);
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
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->nuplocks, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->ndownlocks, matrix->ncols) );

   BMSclearMemoryArray(matrix->nuplocks, matrix->ncols);
   BMSclearMemoryArray(matrix->ndownlocks, matrix->ncols);

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

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->cons, nconss) );

   /* allocate memory for status of finiteness of row sides */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &matrix->isrhsinfinite, nconss) );

#ifdef SCIP_DEBUG
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowname, nconss) );
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivity, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivity, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivityneginf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivityposinf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivityneginf, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivityposinf, nconss) );

   cnt = 0;

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
            assert(cnt < nconss);
            matrix->cons[cnt++] = cons;

#ifdef SCIP_DEBUG
            SCIP_CALL( addConstraint(scip, matrix, SCIPconsGetName(cons), SCIPgetVarsLinear(scip, cons),
                  SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
                  SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons)) );
#else
            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsLinear(scip, cons),
                  SCIPgetValsLinear(scip, cons), SCIPgetNVarsLinear(scip, cons),
                  SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons)) );
#endif
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
            assert(cnt < nconss);
            matrix->cons[cnt++] = cons;

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

#ifdef SCIP_DEBUG
            SCIP_CALL( addConstraint(scip, matrix, SCIPconsGetName(cons), SCIPgetVarsSetppc(scip, cons), NULL,
                  SCIPgetNVarsSetppc(scip, cons), lhs, rhs) );
#else
            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsSetppc(scip, cons), NULL,
                  SCIPgetNVarsSetppc(scip, cons), lhs, rhs) );
#endif
         }
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         for( c = 0; c < nconshdlrconss && (c % 1000 != 0 || !SCIPisStopped(scip)); ++c )
         {
            cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));
            assert(cnt < nconss);
            matrix->cons[cnt++] = cons;

#ifdef SCIP_DEBUG
            SCIP_CALL( addConstraint(scip, matrix, SCIPconsGetName(cons), SCIPgetVarsLogicor(scip, cons),
               NULL, SCIPgetNVarsLogicor(scip, cons), 1.0, SCIPinfinity(scip)) );
#else
            SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsLogicor(scip, cons),
                  NULL, SCIPgetNVarsLogicor(scip, cons), 1.0, SCIPinfinity(scip)) );
#endif
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
               assert(cnt < nconss);
               matrix->cons[cnt++] = cons;

               weights = SCIPgetWeightsKnapsack(scip, cons);
               nvars = SCIPgetNVarsKnapsack(scip, cons);

               if( nvars > valssize )
               {
                  valssize = (int) (1.5 * nvars);
                  SCIP_CALL( SCIPreallocBufferArray(scip, &consvals, valssize) );
               }

               for( v = 0; v < nvars; v++ )
                  consvals[v] = (SCIP_Real)weights[v];

#ifdef SCIP_DEBUG
               SCIP_CALL( addConstraint(scip, matrix, SCIPconsGetName(cons), SCIPgetVarsKnapsack(scip, cons), consvals,
                     SCIPgetNVarsKnapsack(scip, cons), -SCIPinfinity(scip),
                     (SCIP_Real)SCIPgetCapacityKnapsack(scip, cons)) );
#else
               SCIP_CALL( addConstraint(scip, matrix, SCIPgetVarsKnapsack(scip, cons), consvals,
                     SCIPgetNVarsKnapsack(scip, cons), -SCIPinfinity(scip),
                     (SCIP_Real)SCIPgetCapacityKnapsack(scip, cons)) );
#endif
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
               assert(cnt < nconss);
               matrix->cons[cnt++] = cons;

               consvars[0] = SCIPgetVarVarbound(scip, cons);
               consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

               consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

#ifdef SCIP_DEBUG
               SCIP_CALL( addConstraint(scip, matrix, SCIPconsGetName(cons), consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
                     SCIPgetRhsVarbound(scip, cons)) );
#else
               SCIP_CALL( addConstraint(scip, matrix, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
                     SCIPgetRhsVarbound(scip, cons)) );
#endif
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
      assert((*matrix)->nuplocks != NULL);
      assert((*matrix)->ndownlocks != NULL);

      assert((*matrix)->rowmatval != NULL);
      assert((*matrix)->rowmatind != NULL);
      assert((*matrix)->rowmatbeg != NULL);
      assert((*matrix)->rowmatcnt != NULL);
      assert((*matrix)->lhs != NULL);
      assert((*matrix)->rhs != NULL);
#ifdef SCIP_DEBUG
      assert((*matrix)->rowname != NULL);
#endif

      SCIPfreeBufferArray(scip, &((*matrix)->maxactivityposinf));
      SCIPfreeBufferArray(scip, &((*matrix)->maxactivityneginf));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivityposinf));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivityneginf));
      SCIPfreeBufferArray(scip, &((*matrix)->maxactivity));
      SCIPfreeBufferArray(scip, &((*matrix)->minactivity));

#ifdef SCIP_DEBUG
      SCIPfreeBufferArray(scip, &((*matrix)->rowname));
#endif

      SCIPfreeMemoryArray(scip, &((*matrix)->isrhsinfinite));

      SCIPfreeBufferArray(scip, &((*matrix)->cons));

      SCIPfreeBufferArray(scip, &((*matrix)->rhs));
      SCIPfreeBufferArray(scip, &((*matrix)->lhs));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatcnt));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatbeg));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatind));
      SCIPfreeBufferArray(scip, &((*matrix)->rowmatval));

      SCIPfreeBufferArray(scip, &((*matrix)->ndownlocks));
      SCIPfreeBufferArray(scip, &((*matrix)->nuplocks));
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


/** calculate max activity of one row without one variable */
static
SCIP_Real getMaxActSingleRowWithoutCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
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

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);
   assert(0 <= col && col < matrix->ncols);

   maxactivity = 0;

   rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
   rowend = rowpnt + matrix->rowmatcnt[row];
   valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      c = *rowpnt;
      val = *valpnt;

      if( c == col )
         continue;

      assert(!SCIPisInfinity(scip, matrix->lb[c]));
      assert(!SCIPisInfinity(scip, -matrix->ub[c]));

      if( val > 0.0 )
      {
         assert(!SCIPisInfinity(scip, matrix->ub[c]));
         maxactivity += val * matrix->ub[c];
      }
      else if( val < 0.0 )
      {
         assert(!SCIPisInfinity(scip, -matrix->lb[c]));
         maxactivity += val * matrix->lb[c];
      }
   }

   return maxactivity;
}

/** calculate min activity of one row without one variable */
static
SCIP_Real getMinActSingleRowWithoutCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
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

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);
   assert(0 <= col && col < matrix->ncols);

   minactivity = 0;

   rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
   rowend = rowpnt + matrix->rowmatcnt[row];
   valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      c = *rowpnt;
      val = *valpnt;

      if( c == col )
         continue;

      assert(!SCIPisInfinity(scip, matrix->lb[c]));
      assert(!SCIPisInfinity(scip, -matrix->ub[c]));

      if( val > 0.0 )
      {
         assert(!SCIPisInfinity(scip,-matrix->lb[c]));
         minactivity += val * matrix->lb[c];
      }
      else if( val < 0.0 )
      {
         assert(!SCIPisInfinity(scip,matrix->ub[c]));
         minactivity += val * matrix->ub[c];
      }
   }

   return minactivity;
}

/** get minimum/maximum residual activity without the specified variable */
static
SCIP_RETCODE getActivityResiduals(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
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

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(col < matrix->ncols);
   assert(row < matrix->nrows);
   assert(minresactivity != NULL);
   assert(maxresactivity != NULL);
   assert(isminsettoinfinity != NULL);
   assert(ismaxsettoinfinity != NULL);

   lb = SCIPvarGetLbGlobal(matrix->vars[col]);
   ub = SCIPvarGetUbGlobal(matrix->vars[col]);

   assert(!SCIPisInfinity(scip, lb));
   assert(!SCIPisInfinity(scip, -ub));

   *isminsettoinfinity = FALSE;
   *ismaxsettoinfinity = FALSE;

   if( val >= 0.0 )
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(matrix->maxactivityposinf[row] >= 1);
         assert(matrix->maxactivityneginf[row] >= 0);
         if( matrix->maxactivityposinf[row] == 1 && matrix->maxactivityneginf[row] == 0 )
            *maxresactivity = getMaxActSingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (matrix->maxactivityneginf[row] + matrix->maxactivityposinf[row]) > 0 )
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
         else
            *maxresactivity = matrix->maxactivity[row] - val * ub;
      }

      if( SCIPisInfinity(scip, -lb) )
      {
         assert(matrix->minactivityneginf[row] >= 1);
         assert(matrix->minactivityposinf[row] >= 0);
         if( matrix->minactivityneginf[row] == 1 && matrix->minactivityposinf[row] == 0 )
            *minresactivity = getMinActSingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (matrix->minactivityneginf[row] + matrix->minactivityposinf[row]) > 0 )
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
         else
            *minresactivity = matrix->minactivity[row] - val * lb;
      }
   }
   else
   {
      if( SCIPisInfinity(scip, -lb) )
      {
         assert(matrix->maxactivityneginf[row] >= 1);
         assert(matrix->maxactivityposinf[row] >= 0);
         if( matrix->maxactivityneginf[row] == 1 && matrix->maxactivityposinf[row] == 0 )
            *maxresactivity = getMaxActSingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (matrix->maxactivityneginf[row] + matrix->maxactivityposinf[row]) > 0 )
         {
            *maxresactivity = SCIPinfinity(scip);
            *ismaxsettoinfinity = TRUE;
         }
         else
            *maxresactivity = matrix->maxactivity[row] - val * lb;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         assert(matrix->minactivityposinf[row] >= 1);
         assert(matrix->minactivityneginf[row] >= 0);
         if( matrix->minactivityposinf[row] == 1 && matrix->minactivityneginf[row] == 0 )
            *minresactivity = getMinActSingleRowWithoutCol(scip, matrix, row, col);
         else
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
      }
      else
      {
         if( (matrix->minactivityneginf[row] + matrix->minactivityposinf[row]) > 0 )
         {
            *minresactivity = -SCIPinfinity(scip);
            *isminsettoinfinity = TRUE;
         }
         else
            *minresactivity = matrix->minactivity[row] - val * ub;
      }
   }

   return SCIP_OKAY;
}

/** calculate the bounds of this variable from one row */
static
SCIP_RETCODE getVarBoundsOfRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index of variable */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< coefficient of this column in this row */
   SCIP_Real*            rowlb,              /**< lower bound of row */
   SCIP_Bool*            lbfound,            /**< flag indicating that a lower bound was calculated */
   SCIP_Real*            rowub,              /**< upper bound of row */
   SCIP_Bool*            ubfound             /**< flag indicating that a upper bound was calculated */
   )
{
   SCIP_Bool isminsettoinfinity;
   SCIP_Bool ismaxsettoinfinity;
   SCIP_Real minresactivity;
   SCIP_Real maxresactivity;

   assert(rowlb != NULL);
   assert(lbfound != NULL);
   assert(rowub != NULL);
   assert(ubfound != NULL);

   *rowlb = -SCIPinfinity(scip);
   *rowub = SCIPinfinity(scip);
   *lbfound = FALSE;
   *ubfound = FALSE;

   SCIP_CALL( getActivityResiduals(scip, matrix, col, row, val, &minresactivity, &maxresactivity,
         &isminsettoinfinity, &ismaxsettoinfinity) );

   if( val > 0.0 )
   {
      if( !isminsettoinfinity && !SCIPisInfinity(scip, matrix->rhs[row]) )
      {
         *rowub = (matrix->rhs[row] - minresactivity)/val;
         *ubfound = TRUE;
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -matrix->lhs[row]) )
      {
         *rowlb = (matrix->lhs[row] - maxresactivity)/val;
         *lbfound = TRUE;
      }
   }
   else
   {
      if( !isminsettoinfinity && !SCIPisInfinity(scip, matrix->rhs[row]) )
      {
         *rowlb = (matrix->rhs[row] - minresactivity)/val;
         *lbfound = TRUE;
      }

      if( !ismaxsettoinfinity && !SCIPisInfinity(scip, -matrix->lhs[row]) )
      {
         *rowub = (matrix->lhs[row] - maxresactivity)/val;
         *ubfound = TRUE;
      }
   }

   return SCIP_OKAY;
}

/** return true if the bounds of the variable are implied by another constraint */
static
SCIP_RETCODE isVarImpliedFree(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index for implied free test */
   int                   row,                /**< constraint planned for multi-aggregation */
   SCIP_Bool*            impliedfree         /**< flag indicating if this variable is an implied free variable */
   )
{
   SCIP_Real varub;
   SCIP_Real varlb;
   SCIP_Real impliedub;
   SCIP_Real impliedlb;
   int* colpnt;
   int* colend;
   SCIP_Real* valpnt;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);
   assert(0 <= row && row < matrix->nrows);
   assert(impliedfree != NULL);

   *impliedfree = FALSE;

   varub = SCIPvarGetUbGlobal(matrix->vars[col]);
   varlb = SCIPvarGetLbGlobal(matrix->vars[col]);
   impliedub = SCIPinfinity(scip);
   impliedlb = -SCIPinfinity(scip);

   colpnt = matrix->colmatind + matrix->colmatbeg[col];
   colend = colpnt + matrix->colmatcnt[col];
   valpnt = matrix->colmatval + matrix->colmatbeg[col];

   for( ; (colpnt < colend); colpnt++, valpnt++ )
   {
      SCIP_Real rowlb;
      SCIP_Real rowub;
      SCIP_Bool lbfound;
      SCIP_Bool ubfound;

      SCIP_CALL( getVarBoundsOfRow(scip,matrix,col,*colpnt,*valpnt,&rowlb,&lbfound,&rowub,&ubfound) );

      if( lbfound && rowlb > impliedlb )
         impliedlb = rowlb;

      if( ubfound && rowub < impliedub )
         impliedub = rowub;
   }

   /* @todo: we need to consider infinity bounds */
   if( SCIPisFeasLE(scip,impliedub,varub) && SCIPisFeasLE(scip,varlb,impliedlb) )
      *impliedfree = TRUE;

   return SCIP_OKAY;
}

/** calculate the amount of fill-in getting from multi-aggregation */
static
SCIP_Real getFillIn(
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   col,                /**< column index */
   int                   row                 /**< row index */
   )
{
   SCIP_Real fillin;
   int nonzerosnew;
   int nonzerosold;

   assert(matrix != NULL);
   assert(matrix->colmatcnt[col] > 0);
   assert(matrix->rowmatcnt[row] > 0);

   nonzerosold = matrix->rowmatcnt[row] + matrix->colmatcnt[col] - 1;
   nonzerosnew = (matrix->colmatcnt[col] - 1) * (matrix->rowmatcnt[row] - 1);

   fillin = (SCIP_Real)nonzerosnew / (SCIP_Real)nonzerosold;

   return fillin;
}

/** use a simple heuristic if the presence of the variable is advantageous for multi-aggregation */
static
SCIP_RETCODE advConsPresence(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Bool*            beneficial          /**< flag indicating if constraint presence is advantageous */
   )
{
   int* colpnt;
   int* colend;
   int equalitycnt;

   assert(0 <= col && col < matrix->ncols);
   assert(0 <= row && row < matrix->nrows);
   assert(beneficial != NULL);

   colpnt = matrix->colmatind + matrix->colmatbeg[col];
   colend = colpnt + matrix->colmatcnt[col];
   *beneficial = TRUE;
   equalitycnt = 0;

   /* currently we allow the presence of one equality and no ranged rows */
   /* @todo: verify cases with ranged rows and more than one equality */
   for( ; (colpnt < colend); colpnt++ )
   {
      if( !SCIPisInfinity(scip,-matrix->lhs[*colpnt]) )
      {
         if( !SCIPisInfinity(scip,matrix->rhs[*colpnt]) )
         {
            if( SCIPisEQ(scip,matrix->lhs[*colpnt],matrix->rhs[*colpnt]) )
            {
               /* equality */
               equalitycnt++;
               if( equalitycnt > 1 )
               {
                  *beneficial = FALSE;
                  break;
               }
            }
            else
            {
               /* ranged */
               *beneficial = FALSE;
               break;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** identify the candidates for multi-aggregation */
static
SCIP_RETCODE getMultiaggDelcons(
   SCIP*                 scip,               /**< SCIP data structure */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object */
   SCIP_Bool*            multiaggvars,       /**< array indicating multi-aggregatable variables */
   int*                  nummultiaggvars,    /**< number of multi-aggregatable variables */
   int*                  multiaggequalities  /**< array holding equality row indices for deletion */
   )
{
   int r;
   SCIP_Real bestfillin;
   int bestvaridx;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(multiaggvars != NULL);
   assert(nummultiaggvars != NULL);
   assert(multiaggequalities != NULL);

   for( r = 0; r < matrix->nrows; r++ )
   {
      /* we consider only long equalities */
      if( matrix->rowmatcnt[r] > 2 )
      {
         if( SCIPisEQ(scip, matrix->lhs[r], matrix->rhs[r]) )
         {
            int* rowpnt;
            int* rowend;

            rowpnt = matrix->rowmatind + matrix->rowmatbeg[r];
            rowend = rowpnt + matrix->rowmatcnt[r];

            bestfillin = 1.0;
            bestvaridx = -1;
            for( ; rowpnt < rowend; rowpnt++ )
            {
               SCIP_VAR* var;
               SCIP_Bool goodconspresence;

               var = matrix->vars[*rowpnt];

               /* @todo: add methods for discrete variables too */
               if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
                  continue;

               /* search for a continuous variable which is implied free, produces less fill-in
                * and has an advantageous constraint presence
                */
               /* @todo: consider numerical stability ! */

               SCIP_CALL( advConsPresence(scip,matrix,*rowpnt,r,&goodconspresence) );

               if( goodconspresence )
               {
                  SCIP_Bool impliedfree;

                  SCIP_CALL( isVarImpliedFree(scip,matrix,*rowpnt,r,&impliedfree) );

                  if( impliedfree )
                  {
                     SCIP_Real fillin;
                     fillin = getFillIn(matrix,*rowpnt,r);
                     if( fillin < bestfillin )
                     {
                        bestfillin = fillin;
                        bestvaridx = *rowpnt;
                     }
                  }
               }
            }
            if( bestvaridx > -1 && multiaggvars[bestvaridx] != TRUE )
            {
               assert(bestvaridx < matrix->ncols);
               multiaggvars[bestvaridx] = TRUE;
               multiaggequalities[bestvaridx] = r;
               (*nummultiaggvars)++;
            }
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
SCIP_DECL_PRESOLCOPY(presolCopyImplfree)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolImplfree(scip) );

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecImplfree)
{  /*lint --e{715}*/
   CONSTRAINTMATRIX* matrix;
   SCIP_Bool initialized;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   /* do we need this ? */
   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* initialize constraint matrix */
   matrix = NULL;
   initialized = FALSE;
   SCIP_CALL( initMatrix(scip, &matrix, &initialized) );

   if( initialized )
   {
      SCIP_Bool* isvartomultiagg;
      int* multiaggequalities;
      int nummultiaggvars;

      assert(SCIPgetNVars(scip) == matrix->ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &isvartomultiagg, matrix->ncols) );
      BMSclearMemoryArray(isvartomultiagg, matrix->ncols);

      SCIP_CALL( SCIPallocBufferArray(scip, &multiaggequalities, matrix->ncols) );
      BMSclearMemoryArray(multiaggequalities, matrix->ncols);

      nummultiaggvars = 0;

      SCIP_CALL( getMultiaggDelcons(scip,matrix,isvartomultiagg,&nummultiaggvars,multiaggequalities) );

      if( nummultiaggvars > 0 )
      {
         int v;
         for( v = 0; v < matrix->ncols; v++ )
         {
            if( isvartomultiagg[v] )
            {
               SCIP_VAR* multiaggvar;
               SCIP_VAR** vars;
               SCIP_Real* scalars;
               SCIP_Real multiaggcoef;
               SCIP_Real aggrconst;
               SCIP_Bool infeasible;
               SCIP_Bool aggregated;
               SCIP_CONS* multiaggcons;
               int row;
               int cnt;
               int* rowpnt;
               int* rowend;
               SCIP_Real* valpnt;
               int nvars;

               assert(SCIPvarGetType(matrix->vars[v]) == SCIP_VARTYPE_CONTINUOUS);

               multiaggvar = matrix->vars[v];
               row = multiaggequalities[v];
               assert(row < matrix->nrows);
               multiaggcons = matrix->cons[row];

               /* get the number of variables without the multiagg variable */
               nvars = matrix->rowmatcnt[row] - 1;

               /* allocate temporary memory */
               SCIP_CALL( SCIPallocBufferArray(scip, &scalars, nvars) );
               SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

               /* get multi-agg variable coefficient and vars */
               rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
               rowend = rowpnt + matrix->rowmatcnt[row];
               valpnt = matrix->rowmatval + matrix->rowmatbeg[row];
               cnt = 0;
               multiaggcoef = 0.0;
               for( ; rowpnt < rowend; rowpnt++, valpnt++ )
               {
                  if( *rowpnt == v )
                  {
                     multiaggcoef = *valpnt;
                     continue;
                  }

                  vars[cnt++] = matrix->vars[*rowpnt];
               }

               /* avoid division by zero */
               if( SCIPisEQ(scip,multiaggcoef,0.0) )
                  continue;

               /* it is possible to have an implied equality where we have
                * to distinguished to cases
                */
               assert(SCIPisInfinity(scip,-matrix->lhs[row]) || SCIPisInfinity(scip,matrix->rhs[row]) || SCIPisEQ(scip,matrix->rhs[row],matrix->lhs[row]));
               if( !SCIPisInfinity(scip,matrix->rhs[row]) )
                  aggrconst = matrix->rhs[row]/multiaggcoef;
               else
                  aggrconst = matrix->lhs[row]/multiaggcoef;

               /* calculate scalars */
               rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
               rowend = rowpnt + matrix->rowmatcnt[row];
               valpnt = matrix->rowmatval + matrix->rowmatbeg[row];
               cnt = 0;
               SCIPdebugMessage("constraint <%s>: multi-aggregate <%s> ==", SCIPconsGetName(multiaggcons), SCIPvarGetName(multiaggvar));
               for( ; rowpnt < rowend; rowpnt++, valpnt++ )
               {
                  if( *rowpnt == v )
                     continue;

                  scalars[cnt] = -(*valpnt)/multiaggcoef;
                  SCIPdebugPrintf(" %+.15g<%s>", scalars[cnt], SCIPvarGetName(matrix->vars[*rowpnt]));
                  cnt++;
               }
               SCIPdebugPrintf(" %+.15g, bounds of <%s>: [%.15g,%.15g]\n",
                  aggrconst, SCIPvarGetName(multiaggvar), SCIPvarGetLbGlobal(multiaggvar), SCIPvarGetUbGlobal(multiaggvar));

               /* perform the multi-aggregation */
               SCIP_CALL( SCIPmultiaggregateVar(scip, multiaggvar, nvars, vars, scalars, aggrconst, &infeasible, &aggregated) );
               assert(aggregated);

               SCIPfreeBufferArray(scip, &vars);
               SCIPfreeBufferArray(scip, &scalars);

               /* check for infeasible aggregation */
               if( infeasible )
               {
                  SCIPdebugMessage("constraint <%s>: infeasible multi-aggregation\n", SCIPconsGetName(multiaggcons));
                  return SCIP_OKAY;
               }

               (*naggrvars)++;

               /* delete constraint */
               SCIP_CALL( SCIPdelCons(scip,multiaggcons) );
               (*ndelconss)++;
            }
         }
      }

      SCIPfreeBufferArray(scip, &multiaggequalities);
      SCIPfreeBufferArray(scip, &isvartomultiagg);
   }

   freeMatrix(scip, &matrix);

   return SCIP_OKAY;
}

/*
 * presolver specific interface methods
 */

/** creates the implied free presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolImplfree(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOL* presol;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_DELAY, presolExecImplfree, NULL) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyImplfree) );

   return SCIP_OKAY;
}
