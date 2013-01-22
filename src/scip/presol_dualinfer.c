/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*#define SCIP_DEBUG*/

/**@file    presol_dualinfer.c
 * @ingroup PRESOLVERS
 * @brief   dual inference presolver
 * @author  Dieter Weninger
 *
 * This presolver exploits dual information for primal variable fixings.
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

#include "presol_dualinfer.h"

#define PRESOL_NAME             "dualinfer"
#define PRESOL_DESC             "exploit dual informations for variable fixings"
#define PRESOL_PRIORITY         20010000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY                TRUE     /**< should presolver be delayed, if other presolvers found reductions? */

#define MAX_LOOPS                      7     /**< maximal number of dual bound strengthening loops */


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

/** relation type */
enum Relationtype
{
   LESSEQUAL,
   GREATEREQUAL,
   EQUALITY,
   RANGED
};
typedef enum Relationtype RELATIONTYPE;



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
#ifdef SCIP_DEBUG
   const char**          rowname;            /**< name of row */
#endif
   int                   nrows;              /**< complete number of rows */
   SCIP_Real*            lhs;                /**< left hand side per row */
   SCIP_Real*            rhs;                /**< right hand side per row */
   SCIP_Bool*            islhsinfinite;      /**< is left hand side -infinity */
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

   assert(vars != NULL);
   assert(vals != NULL);

   rowidx = matrix->nrows;

   matrix->lhs[rowidx] = lhs;
   matrix->rhs[rowidx] = rhs;

   if( SCIPisInfinity(scip, -matrix->lhs[rowidx]) )
      matrix->islhsinfinite[rowidx] = TRUE;
   if( SCIPisInfinity(scip, matrix->rhs[rowidx]) )
      matrix->isrhsinfinite[rowidx] = TRUE;

#ifdef SCIP_DEBUG
   matrix->rowname[rowidx] = name;
#endif

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

         assert(matrix->ncols > col);

         assert(!SCIPisInfinity(scip, matrix->lb[col]));
         assert(!SCIPisInfinity(scip, -matrix->ub[col]));

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
         else if( val < 0.0 )
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

      if( (matrix->minactivityneginf[row] + matrix->minactivityposinf[row]) > 0 )
         matrix->minactivity[row] = -SCIPinfinity(scip);
      if( (matrix->maxactivityneginf[row] + matrix->maxactivityposinf[row]) > 0 )
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

   /* allocate memory for status of finiteness of row sides */
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &matrix->islhsinfinite, nconss) );
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
      SCIPfreeMemoryArray(scip, &((*matrix)->islhsinfinite));

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

/** return if row is redundant or not */
static
SCIP_Bool isRowRedundant(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   row                 /**< row index */
   )
{
   SCIP_Bool redundant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real minact;
   SCIP_Real maxact;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);

   lhs = matrix->lhs[row];
   rhs = matrix->rhs[row];
   minact = matrix->minactivity[row];
   maxact = matrix->maxactivity[row];

   redundant = SCIPisFeasLE(scip, lhs, minact) && SCIPisFeasLE(scip, maxact, rhs);

   return redundant;
}

/** initialize shadow prices */
static
void initShadowPrices(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   RELATIONTYPE*         rowrelation,        /**< row relation type */
   SCIP_Real*            lowershadow,        /**< lower shadow prices */
   SCIP_Real*            uppershadow         /**< upper shadow prices */
   )
{
   int r;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(rowrelation != NULL);
   assert(lowershadow != NULL);
   assert(uppershadow != NULL);

   for( r = 0; r < matrix->nrows; r++ )
   {
      if( rowrelation[r] == EQUALITY )
      {
         lowershadow[r] = -SCIPinfinity(scip);
         uppershadow[r] = SCIPinfinity(scip);
      }
      else if( rowrelation[r] == RANGED )
      {
         lowershadow[r] = -SCIPinfinity(scip);
         uppershadow[r] = SCIPinfinity(scip);
      }
      else if( rowrelation[r] == LESSEQUAL )
      {
         lowershadow[r] = 0;
         uppershadow[r] = SCIPinfinity(scip);
      }
      else if( rowrelation[r] == GREATEREQUAL )
      {
         lowershadow[r] = -SCIPinfinity(scip);
         uppershadow[r] = 0;
      }
      else
      {
         lowershadow[r] = 0;
         uppershadow[r] = 0;
      }
   }
}

/** search for singleton columns and update shadow prices */
static
SCIP_RETCODE singletonColumns(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   RELATIONTYPE*         rowrelation,        /**< row relation type */
   SCIP_Real*            lowershadow,        /**< lower shadows */
   SCIP_Real*            uppershadow,        /**< upper shadows */
   int*                  nfitsinglecols      /**< number of fitting singleton columns */
   )
{
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(rowrelation != NULL);
   assert(lowershadow != NULL);
   assert(uppershadow != NULL);
   assert(nfitsinglecols != NULL);

   for( c = 0; c < matrix->ncols; c++ )
   {
      if( matrix->colmatcnt[c] == 1 && SCIPvarGetType(matrix->vars[c]) == SCIP_VARTYPE_CONTINUOUS &&
         SCIPisInfinity(scip, SCIPvarGetUbGlobal(matrix->vars[c]))
         /*&& SCIPisEQ(scip,SCIPvarGetLbGlobal(matrix->vars[c]),0)*/ )
      {
         int row;

         row = *(matrix->colmatind + matrix->colmatbeg[c]);

         if( rowrelation[row] == LESSEQUAL || rowrelation[row] == GREATEREQUAL )
         {
            SCIP_Real val;
            SCIP_Real tmp;

            val = *(matrix->colmatval + matrix->colmatbeg[c]);
            tmp = -SCIPvarGetObj(matrix->vars[c]) / val;

            (*nfitsinglecols)++;

            if( val > 0 )
            {
               if( tmp > lowershadow[row] )
                  lowershadow[row] = tmp;
            }
            else
            {
               assert(val < 0);

               if( tmp < uppershadow[row] )
                  uppershadow[row] = tmp;
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** calculate min/max costs per column */
static
void costCalculation(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   RELATIONTYPE*         rowrelation,        /**< row relation type */
   SCIP_Real*            lowershadow,        /**< lower shadows */
   SCIP_Real*            uppershadow,        /**< upper shadows */
   SCIP_Real*            lowercosts,         /**< lower shadow costs */
   SCIP_Real*            uppercosts          /**< upper shadow costs */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   SCIP_Bool mininfinite;
   SCIP_Bool maxinfinite;
   int row;
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(rowrelation != NULL);
   assert(lowershadow != NULL);
   assert(uppershadow != NULL);
   assert(lowercosts != NULL);
   assert(uppercosts != NULL);

   for( c = 0; c < matrix->ncols; c++)
   {
      mininfinite = FALSE;
      maxinfinite = FALSE;
      lowercosts[c] = 0;
      uppercosts[c] = 0;

      colpnt = matrix->colmatind + matrix->colmatbeg[c];
      colend = colpnt + matrix->colmatcnt[c];
      valpnt = matrix->colmatval + matrix->colmatbeg[c];

      for( ; (colpnt < colend) && (!mininfinite || !maxinfinite); colpnt++, valpnt++ )
      {
         row = *colpnt;
         val = *valpnt;

         if( !isRowRedundant(scip, matrix, row) )
         {
            if( val >= 0.0 )
            {
               mininfinite = mininfinite || SCIPisInfinity(scip, -lowershadow[row]);
               maxinfinite = maxinfinite || SCIPisInfinity(scip, uppershadow[row]);

               if( !mininfinite )
               {
                  lowercosts[c] += val * lowershadow[row];
               }
               if( !maxinfinite )
               {
                  uppercosts[c] += val * uppershadow[row];
               }
            }
            else
            {
               mininfinite = mininfinite || SCIPisInfinity(scip, uppershadow[row]);
               maxinfinite = maxinfinite || SCIPisInfinity(scip, -lowershadow[row]);

               if( !mininfinite )
               {
                  lowercosts[c] += val * uppershadow[row];
               }
               if( !maxinfinite )
               {
                  uppercosts[c] += val * lowershadow[row];
               }
            }
         }
      }

      if( mininfinite )
         lowercosts[c] = -SCIPinfinity(scip);
      if( maxinfinite )
         uppercosts[c] = SCIPinfinity(scip);
   }
}

/** fix variables at their bounds out of the min/max costs */
static
void fixColumns(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   RELATIONTYPE*         rowrelation,        /**< row relation type */
   SCIP_Real*            lowercosts,         /**< lower shadow costs */
   SCIP_Real*            uppercosts,         /**< upper shadow costs */
   int*                  npossiblefixings,   /**< number of possible fixings */
   FIXINGDIRECTION*      varstofix           /**< array holding information for later upper/lower bound fixing */
   )
{
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(rowrelation != NULL);
   assert(lowercosts != NULL);
   assert(uppercosts != NULL);
   assert(npossiblefixings != NULL);
   assert(varstofix != NULL);

   for( c = 0; c < matrix->ncols; c++ )
   {
      SCIP_Real objval;

      objval = -SCIPvarGetObj(matrix->vars[c]);

      if( !SCIPisInfinity(scip, lowercosts[c]) && !SCIPisInfinity(scip, -lowercosts[c]) )
      {
         if( lowercosts[c] > objval )
         {
            if( !SCIPisInfinity(scip, -SCIPvarGetLbGlobal(matrix->vars[c])) )
            {
               varstofix[c] = FIXATLB;
               (*npossiblefixings)++;
            }
         }
      }
      if( !SCIPisInfinity(scip, uppercosts[c]) && !SCIPisInfinity(scip, -uppercosts[c]) )
      {
         if( uppercosts[c] < objval )
         {
            if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(matrix->vars[c])) )
            {
               varstofix[c] = FIXATUB;
               (*npossiblefixings)++;
            }
         }
      }
   }
}

/** dual cost fixing */
static
SCIP_RETCODE costFixing(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   RELATIONTYPE*         rowrelation,        /**< row relation type */
   FIXINGDIRECTION*      varstofix,          /**< array holding information for later upper/lower bound fixing */
   int*                  nfitsinglecols,     /**< number of continuous singleton columns */
   int*                  npossiblefixings    /**< number of possible fixings */
   )
{
   SCIP_Real* lowershadow;
   SCIP_Real* uppershadow;
   SCIP_Real* lowercosts;
   SCIP_Real* uppercosts;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(rowrelation != NULL);
   assert(varstofix != NULL);
   assert(nfitsinglecols != NULL);
   assert(npossiblefixings != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &lowershadow, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &uppershadow, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lowercosts, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &uppercosts, matrix->ncols) );

   initShadowPrices(scip, matrix, rowrelation, lowershadow, uppershadow);

   SCIP_CALL( singletonColumns(scip, matrix, rowrelation, lowershadow, uppershadow, nfitsinglecols) );

   costCalculation(scip, matrix, rowrelation, lowershadow, uppershadow, lowercosts, uppercosts);

   fixColumns(scip, matrix, rowrelation, lowercosts, uppercosts, npossiblefixings, varstofix);

   SCIPfreeBufferArray(scip, &uppercosts);
   SCIPfreeBufferArray(scip, &lowercosts);
   SCIPfreeBufferArray(scip, &uppershadow);
   SCIPfreeBufferArray(scip, &lowershadow);

   return SCIP_OKAY;
}

/** calculate maximal column activity from one continuous variable without one row */
static
SCIP_Real getMaxColActWithoutRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   withoutrow,         /**< exclude row index */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of dual variables */
   RELATIONTYPE*         rowrelation         /**< row relations */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   SCIP_Real maxcolactivity;
   int row;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);
   assert(0 <= withoutrow && withoutrow < matrix->nrows);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(rowrelation != NULL);
   assert(SCIPvarGetType(matrix->vars[col]) == SCIP_VARTYPE_CONTINUOUS);

   maxcolactivity = 0;

   colpnt = matrix->colmatind + matrix->colmatbeg[col];
   colend = colpnt + matrix->colmatcnt[col];
   valpnt = matrix->colmatval + matrix->colmatbeg[col];

   for(; (colpnt < colend); colpnt++, valpnt++ )
   {
      row = *colpnt;
      val = *valpnt;

      assert(rowrelation[row] == LESSEQUAL || rowrelation[row] == GREATEREQUAL);

      if( row == withoutrow )
         continue;

      if( rowrelation[row] == LESSEQUAL )
         val = -val;

      if( val > 0 )
      {
         assert(!SCIPisInfinity(scip, ubdual[row]));
         maxcolactivity += val * ubdual[row];
      }
      else if( val < 0.0 )
      {
         assert(!SCIPisInfinity(scip, -lbdual[row]));
         maxcolactivity += val * lbdual[row];
      }
   }
   return maxcolactivity;
}

/** calculate minimal column activity from one continuous variable without one row */
static
SCIP_Real getMinColActWithoutRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   withoutrow,         /**< exclude row index */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of dual variables */
   RELATIONTYPE*         rowrelation         /**< row relations */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   SCIP_Real mincolactivity;
   int row;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);
   assert(0 <= withoutrow && withoutrow < matrix->nrows);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(rowrelation != NULL);
   assert(SCIPvarGetType(matrix->vars[col]) == SCIP_VARTYPE_CONTINUOUS);

   mincolactivity = 0;

   colpnt = matrix->colmatind + matrix->colmatbeg[col];
   colend = colpnt + matrix->colmatcnt[col];
   valpnt = matrix->colmatval + matrix->colmatbeg[col];

   for(; (colpnt < colend); colpnt++, valpnt++ )
   {
      row = *colpnt;
      val = *valpnt;

      assert(rowrelation[row] == LESSEQUAL || rowrelation[row] == GREATEREQUAL);

      if( row == withoutrow )
         continue;

      if( rowrelation[row] == LESSEQUAL )
         val = -val;

      if( val > 0 )
      {
         assert(!SCIPisInfinity(scip, -lbdual[row]));
         mincolactivity += val * lbdual[row];
      }
      else if( val < 0.0 )
      {
         assert(!SCIPisInfinity(scip, ubdual[row]));
         mincolactivity += val * ubdual[row];
      }
   }
   return mincolactivity;
}

/** calculate minimal/maximal column residual activities */
static
void calcColActivityResiduals(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   col,                /**< column index */
   int                   row,                /**< row index */
   SCIP_Real             val,                /**< matrix coefficient */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of dual variables */
   RELATIONTYPE*         rowrelation,        /**< row relations */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactposinf,    /**< number of maximal column activity positive infinities */
   int*                  maxcolactneginf,    /**< number of maximal column activity negative infinities */
   int*                  mincolactposinf,    /**< number of minimal column activity positive infinities */
   int*                  mincolactneginf,    /**< number of minimal column activity negative infinities */
   SCIP_Real*            mincolresact,       /**< minimal column residual activity */
   SCIP_Real*            maxcolresact        /**< maximal column residual activity */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= col && col < matrix->ncols);
   assert(0 <= row && row < matrix->nrows);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(rowrelation != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactposinf != NULL);
   assert(maxcolactneginf != NULL);
   assert(mincolactposinf != NULL);
   assert(mincolactneginf != NULL);
   assert(mincolresact != NULL);
   assert(maxcolresact != NULL);
   assert(rowrelation[row] == LESSEQUAL || rowrelation[row] == GREATEREQUAL);
   assert(SCIPvarGetType(matrix->vars[col]) == SCIP_VARTYPE_CONTINUOUS);

   /* we assume always a >= relation */
   if( rowrelation[row] == LESSEQUAL )
      val = -val;

   if( val > 0.0 )
   {
      if( SCIPisInfinity(scip, ubdual[row]) )
      {
         assert(maxcolactposinf[col] >= 1);
         if( maxcolactposinf[col] == 1 && maxcolactneginf[col] == 0 )
            *maxcolresact = getMaxColActWithoutRow(scip, matrix, col, row, lbdual, ubdual, rowrelation);
         else
            *maxcolresact = SCIPinfinity(scip);
      }
      else
      {
         if( (maxcolactneginf[col]+maxcolactposinf[col]) > 0 )
            *maxcolresact = SCIPinfinity(scip);
         else
            *maxcolresact = maxcolact[col] - val * ubdual[row];
      }

      if( SCIPisInfinity(scip,-lbdual[row]) )
      {
         assert(mincolactneginf[col] >= 1);
         if( mincolactneginf[col] == 1 && mincolactposinf[col] == 0 )
            *mincolresact = getMinColActWithoutRow(scip, matrix, col, row, lbdual, ubdual, rowrelation);
         else
            *mincolresact = -SCIPinfinity(scip);
      }
      else
      {
         if( (mincolactneginf[col]+mincolactposinf[col]) > 0 )
            *mincolresact = -SCIPinfinity(scip);
         else
            *mincolresact = mincolact[col] - val * lbdual[row];
      }
   }
   else if( val < 0.0 )
   {
      if( SCIPisInfinity(scip,-lbdual[row]) )
      {
         assert(maxcolactneginf[col] >= 1);
         if( maxcolactneginf[col] == 1 && maxcolactposinf[col] == 0 )
            *maxcolresact = getMaxColActWithoutRow(scip, matrix, col, row, lbdual, ubdual, rowrelation);
         else
            *maxcolresact = SCIPinfinity(scip);
      }
      else
      {
         if( (maxcolactneginf[col]+maxcolactposinf[col]) > 0 )
            *maxcolresact = SCIPinfinity(scip);
         else
            *maxcolresact = maxcolact[col] - val * lbdual[row];
      }

      if( SCIPisInfinity(scip,ubdual[row]) )
      {
         assert(mincolactposinf[col] >= 1);
         if( mincolactposinf[col] == 1 && mincolactneginf[col] == 0 )
            *mincolresact = getMinColActWithoutRow(scip, matrix, col, row, lbdual, ubdual, rowrelation);
         else
            *mincolresact = -SCIPinfinity(scip);
      }
      else
      {
         if( (mincolactneginf[col]+mincolactposinf[col]) > 0 )
            *mincolresact = -SCIPinfinity(scip);
         else
            *mincolresact = mincolact[col] - val * ubdual[row];
      }
   }
}

/** calculate minimal/maximal column activity on continuous variables */
static
void calcColActivity(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of dual variables */
   int                   startcol,           /**< start column index */
   int                   stopcol,            /**< stop column index */
   RELATIONTYPE*         rowrelation,        /**< constraint relations */
   SCIP_Bool*            excludevar,         /**< indicating of variable is within an equation or ranged row */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactposinf,    /**< number of maximal column activity positive infinities */
   int*                  maxcolactneginf,    /**< number of maximal column activity negative infinities */
   int*                  mincolactposinf,    /**< number of minimal column activity positive infinities */
   int*                  mincolactneginf     /**< number of minimal column activity negative infinities */
   )
{
   SCIP_Real* valpnt;
   int* colpnt;
   int* colend;
   SCIP_Real val;
   int row;
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(0 <= startcol && startcol < matrix->ncols);
   assert(0 < stopcol && stopcol <= matrix->ncols);
   assert(rowrelation != NULL);
   assert(excludevar != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactposinf != NULL);
   assert(maxcolactneginf != NULL);
   assert(mincolactposinf != NULL);
   assert(mincolactneginf != NULL);

   for( c = startcol; c < stopcol; c++ )
   {
      excludevar[c] = FALSE;
      mincolact[c] = 0;
      maxcolact[c] = 0;
      maxcolactposinf[c] = 0;
      maxcolactneginf[c] = 0;
      mincolactposinf[c] = 0;
      mincolactneginf[c] = 0;

      /* no activity calculation for non-continuous variables */
      if( SCIPvarGetType(matrix->vars[c]) != SCIP_VARTYPE_CONTINUOUS )
         continue;

      colpnt = matrix->colmatind + matrix->colmatbeg[c];
      colend = colpnt + matrix->colmatcnt[c];
      valpnt = matrix->colmatval + matrix->colmatbeg[c];

      /* calculate column activities */
      for(; (colpnt < colend); colpnt++, valpnt++ )
      {
         row = *colpnt;
         val = *valpnt;

         if( rowrelation[row] == LESSEQUAL || rowrelation[row] == GREATEREQUAL )
         {
            if( rowrelation[row] == LESSEQUAL )
               val = -val;

            if( val > 0 )
            {
               if(SCIPisInfinity(scip, ubdual[row]))
                  maxcolactposinf[c]++;
               else
                  maxcolact[c] += val * ubdual[row];

               if(SCIPisInfinity(scip, -lbdual[row]))
                  mincolactneginf[c]++;
               else
                  mincolact[c] += val * lbdual[row];
            }
            else if( val < 0.0 )
            {
               if(SCIPisInfinity(scip, -lbdual[row]))
                  maxcolactneginf[c]++;
               else
                  maxcolact[c] += val * lbdual[row];

               if(SCIPisInfinity(scip, ubdual[row]))
                  mincolactposinf[c]++;
               else
                  mincolact[c] += val * ubdual[row];
            }
         }
         else if( rowrelation[row] == EQUALITY || rowrelation[row] == RANGED )
         {
            excludevar[c] = TRUE;
            break;
         }
      }

      /* consider lower bounds on variables */
      if( SCIPisPositive(scip, SCIPvarGetLbGlobal(matrix->vars[c])) )
         maxcolactposinf[c]++;

      /* consider upper bounds on variables */
      if( !SCIPisInfinity(scip, SCIPvarGetUbGlobal(matrix->vars[c])) )
         mincolactposinf[c]++;

      /* update column activities */
      if( (mincolactneginf[c] + mincolactposinf[c]) > 0 )
         mincolact[c] = -SCIPinfinity(scip);
      if( (maxcolactneginf[c] + maxcolactposinf[c]) > 0 )
         maxcolact[c] = SCIPinfinity(scip);
   }
}


/** update bounds on dual variables */
static
void updateDualBounds(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   SCIP_Real             objval,             /**< objective function value */
   SCIP_Real             val,                /**< matrix coefficient */
   int                   row,                /**< row index */
   RELATIONTYPE*         rowrelation,        /**< row relation type */
   SCIP_Real             mincolresact,       /**< minimal column residual activity */
   SCIP_Real*            lbdual,             /**< dual lower bounds */
   SCIP_Real*            ubdual,             /**< dual upper bounds */
   int*                  boundchanges,       /**< number of bound changes */
   SCIP_Bool*            updateinfcnt        /**< flag indicating to update infinity counters */
   )
{
   SCIP_Real newlbdual;
   SCIP_Real newubdual;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(rowrelation != NULL);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(boundchanges != NULL);
   assert(updateinfcnt != NULL);
   assert(rowrelation[row] == LESSEQUAL || rowrelation[row] == GREATEREQUAL);

   *updateinfcnt = FALSE;

   if( rowrelation[row] == LESSEQUAL )
      val = -val;

   if( val > 0 )
   {
      if( !SCIPisInfinity(scip, mincolresact) &&
         !SCIPisInfinity(scip, -mincolresact) )
      {
         /* upper bound calculation on dual variable */
         newubdual = (objval - mincolresact) / val;
         if( newubdual < ubdual[row] )
         {
            if( SCIPisInfinity(scip, ubdual[row]) )
               *updateinfcnt = TRUE;

            ubdual[row] = newubdual;
            (*boundchanges)++;
         }
      }
   }
   else if( val < 0 )
   {
      if( !SCIPisInfinity(scip, mincolresact) &&
         !SCIPisInfinity(scip, -mincolresact) )
      {
         /* lower bound calclulation on dual variable */
         newlbdual = (objval - mincolresact) / val;
         if( newlbdual > lbdual[row] )
         {
            if( SCIPisInfinity(scip, -lbdual[row]) )
               *updateinfcnt = TRUE;

            lbdual[row] = newlbdual;
            (*boundchanges)++;
         }
      }
   }
}

/** update minimal/maximal column activity infinity counters */
static
void infCounterUpdate(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   SCIP_Real             val,                /**< matrix coefficient */
   int                   row,                /**< row index */
   SCIP_Real*            lbdual,             /**< lower bounds of dual variables */
   SCIP_Real*            ubdual,             /**< upper bounds of dual variables */
   RELATIONTYPE*         rowrelation,        /**< row relation */
   SCIP_Bool*            excludevar,         /**< indicating if variable is within an equation or ranged row */
   SCIP_Real*            mincolact,          /**< minimal column activities */
   SCIP_Real*            maxcolact,          /**< maximal column activities */
   int*                  maxcolactposinf,    /**< number of maximal column activity positive infinities */
   int*                  maxcolactneginf,    /**< number of maximal column activity negative infinities */
   int*                  mincolactposinf,    /**< number of minimal column activity positive infinities */
   int*                  mincolactneginf     /**< number of minimal column activity negative infinities */
   )
{
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   SCIP_Real aij;
   int c;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(0 <= row && row < matrix->nrows);
   assert(lbdual != NULL);
   assert(ubdual != NULL);
   assert(rowrelation != NULL);
   assert(excludevar != NULL);
   assert(mincolact != NULL);
   assert(maxcolact != NULL);
   assert(maxcolactposinf != NULL);
   assert(maxcolactneginf != NULL);
   assert(mincolactposinf != NULL);
   assert(mincolactneginf != NULL);
   assert(rowrelation[row] == LESSEQUAL || rowrelation[row] == GREATEREQUAL);

   rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
   rowend = rowpnt + matrix->rowmatcnt[row];
   valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

   /* look at all columns entries present within row and update the corresponding infinity counters.
      if one counter gets to zero, then calculate this column activity completely new */

   if( rowrelation[row] == LESSEQUAL )
      val = -val;

   if( val > 0 )
   {
      /* we did an upper dual bound change from +infinity to a value < infinity */
      for(; (rowpnt < rowend); rowpnt++, valpnt++ )
      {
         c = *rowpnt;
         aij = *valpnt;

         if( excludevar[c] || SCIPvarGetType(matrix->vars[c]) != SCIP_VARTYPE_CONTINUOUS )
            continue;

         if( rowrelation[row] == LESSEQUAL )
            aij = -aij;

         if( aij > 0 )
         {
            assert(maxcolactposinf[c] > 0);
            maxcolactposinf[c]--;

            if( maxcolactposinf[c] == 0 )
            {
               calcColActivity(scip, matrix, lbdual, ubdual, c, c+1,
                  rowrelation, excludevar, mincolact, maxcolact,
                  maxcolactposinf, maxcolactneginf, mincolactposinf, mincolactneginf);
            }
         }
         else if( aij < 0 )
         {
            assert(mincolactposinf[c] > 0);
            mincolactposinf[c]--;

            if( mincolactposinf[c] == 0 )
            {
               calcColActivity(scip, matrix, lbdual, ubdual, c, c+1,
                  rowrelation, excludevar, mincolact, maxcolact,
                  maxcolactposinf, maxcolactneginf, mincolactposinf, mincolactneginf);
            }
         }
      }
   }
   else if( val < 0 )
   {
      /* we did a lower dual bound change from -infinity to a value > -infinity */
      for(; (rowpnt < rowend); rowpnt++, valpnt++ )
      {
         c = *rowpnt;
         aij = *valpnt;

         if( excludevar[c] || SCIPvarGetType(matrix->vars[c]) != SCIP_VARTYPE_CONTINUOUS )
            continue;

         if( rowrelation[row] == LESSEQUAL )
            aij = -aij;

         if( aij > 0 )
         {
            assert(mincolactneginf[c] > 0);
            mincolactneginf[c]--;

            if( mincolactneginf[c] == 0 )
               calcColActivity(scip, matrix, lbdual, ubdual, c, c+1,
                  rowrelation, excludevar, mincolact, maxcolact,
                  maxcolactposinf, maxcolactneginf, mincolactposinf, mincolactneginf);
         }
         else if( aij < 0 )
         {
            assert(maxcolactneginf[c] > 0);
            maxcolactneginf[c]--;

            if( maxcolactneginf[c] == 0 )
               calcColActivity(scip, matrix, lbdual, ubdual, c, c+1,
                  rowrelation, excludevar, mincolact, maxcolact,
                  maxcolactposinf, maxcolactneginf, mincolactposinf, mincolactneginf);
         }
      }
   }
}

/** dual bound strengthening on continuous variables */
static
SCIP_RETCODE dualBoundStrengthening(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   RELATIONTYPE*         rowrelation,        /**< row relation type */
   FIXINGDIRECTION*      varstofix,          /**< array holding information for later upper/lower bound fixing */
   int*                  nfitsinglecols,     /**< number of continuous singleton columns */
   int*                  npossiblefixings    /**< found number of possible fixings */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* lbdual;
   SCIP_Real* ubdual;
   SCIP_Real* mincolact;
   SCIP_Real* maxcolact;
   SCIP_Bool* excludevar;
   SCIP_Bool* varissingcol;
   int* maxcolactposinf;
   int* maxcolactneginf;
   int* mincolactposinf;
   int* mincolactneginf;
   int* colpnt;
   int* colend;
   int cnt;
   int boundchanges;
   int loops;
   int c;
   int r;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(rowrelation != NULL);
   assert(varstofix != NULL);
   assert(nfitsinglecols != NULL);
   assert(npossiblefixings != NULL);

   /* all variables must have a lower bound >= 0 */
   for( c = 0; c < matrix->ncols; c++ )
   {
      if( SCIPisNegative(scip, SCIPvarGetLbGlobal(matrix->vars[c])) )
         return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varissingcol, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &excludevar, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincolact, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolact, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolactposinf, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcolactneginf, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincolactposinf, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &mincolactneginf, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &lbdual, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ubdual, matrix->nrows) );

   for( r = 0; r < matrix->nrows; r++ )
   {
      /* initialize dual bounds */
      lbdual[r] = 0.0;
      ubdual[r] = SCIPinfinity(scip);
   }

   /* calculate bounds on dual variables from singleton continuous columns */
   cnt = 0;
   for( c = 0; c < matrix->ncols; c++ )
   {
      int row;

      row = *(matrix->colmatind + matrix->colmatbeg[c]);
      varissingcol[c] = FALSE;

      if( matrix->colmatcnt[c] == 1 && rowrelation[row] != EQUALITY && rowrelation[row] != RANGED )
      {
         varissingcol[c] = TRUE;

         if( SCIPvarGetType(matrix->vars[c]) == SCIP_VARTYPE_CONTINUOUS &&
            /*SCIPisZero(scip, SCIPvarGetLbGlobal(matrix->vars[c])) &&*/
            SCIPisInfinity(scip, SCIPvarGetUbGlobal(matrix->vars[c])) )
         {
            SCIP_Real objval;
            SCIP_Real val;
            SCIP_Real newubdual;
            SCIP_Real newlbdual;

            objval = SCIPvarGetObj(matrix->vars[c]);
            val = *(matrix->colmatval + matrix->colmatbeg[c]);

            if( rowrelation[row] == LESSEQUAL )
               val = -val;

            if( val > 0 )
            {
               newubdual = objval/val;
               if( newubdual < ubdual[row] )
               {
                  ubdual[row] = newubdual;
                  cnt++;
               }
            }
            else if( val < 0 )
            {
               newlbdual = objval/val;
               if( newlbdual > lbdual[row] )
               {
                  lbdual[row] = newlbdual;
                  cnt++;
               }
            }
         }
      }
   }
   (*nfitsinglecols) = cnt;

   /* do dual bound strengthening only if fitting singleton continuous variables are present */
   if( cnt > 0 )
   {
      loops = 0;
      boundchanges = 1;

      /* continue dual bound strengthening only if boundchanges are occurred and
       * a maximal number of processing loops is not exceeded
       */
      while( boundchanges && loops < MAX_LOOPS )
      {
         loops++;
         boundchanges = 0;

         /* calculate minimal and maximal column activity */
         calcColActivity(scip, matrix, lbdual, ubdual, 0, matrix->ncols,
            rowrelation, excludevar, mincolact, maxcolact,
            maxcolactposinf, maxcolactneginf, mincolactposinf, mincolactneginf);

         for( c = 0 ; c < matrix->ncols; c++ )
         {
            SCIP_Real objval;

            /* do no dual bound strengthing for variables which are present within
             * equalities or ranged rows, for singleton columns and non-continuous variables
             */
            if( excludevar[c] || varissingcol[c] || (SCIPvarGetType(matrix->vars[c]) != SCIP_VARTYPE_CONTINUOUS))
               continue;

            objval = SCIPvarGetObj(matrix->vars[c]);

            colpnt = matrix->colmatind + matrix->colmatbeg[c];
            colend = colpnt + matrix->colmatcnt[c];
            valpnt = matrix->colmatval + matrix->colmatbeg[c];

            for(; (colpnt < colend); colpnt++, valpnt++ )
            {
               int row;
               SCIP_Real val;
               SCIP_Real mincolresact;
               SCIP_Real maxcolresact;
               SCIP_Bool updateinfcnt;

               row = *colpnt;
               val = *valpnt;
               mincolresact = -SCIPinfinity(scip);
               maxcolresact = SCIPinfinity(scip);
               updateinfcnt = FALSE;

               /* calulate column activity residuals */
               calcColActivityResiduals(scip, matrix, c, row, val, lbdual, ubdual, rowrelation,
                  mincolact, maxcolact, maxcolactposinf, maxcolactneginf,
                  mincolactposinf, mincolactneginf, &mincolresact, &maxcolresact);

               /* update dual bounds */
               updateDualBounds(scip, matrix, objval, val, row, rowrelation, mincolresact,
                  lbdual, ubdual, &boundchanges, &updateinfcnt);

               /* update infinity counters if bound changed properly */
               if( updateinfcnt )
                  infCounterUpdate(scip, matrix, val, row, lbdual, ubdual,
                     rowrelation, excludevar, mincolact, maxcolact,
                     maxcolactposinf, maxcolactneginf, mincolactposinf, mincolactneginf);
            }
         }
      }

      if( boundchanges > 0 )
      {
         /* calculate minimal/maximal column activity the last time */
         calcColActivity(scip, matrix, lbdual, ubdual, 0, matrix->ncols,
            rowrelation, excludevar, mincolact, maxcolact,
            maxcolactposinf, maxcolactneginf, mincolactposinf, mincolactneginf);
      }

      for( c = 0; c < matrix->ncols; c++ )
      {
         SCIP_Real objval;
         objval = SCIPvarGetObj(matrix->vars[c]);

         /* exclude non-continuous variables and variables which are
            present within equalities */
         if( (SCIPvarGetType(matrix->vars[c]) != SCIP_VARTYPE_CONTINUOUS) || excludevar[c] )
            continue;

         /* apply complementary slackness: (yA-c)_i > 0 => x_i = 0 */
         if( SCIPisLT(scip, maxcolact[c], objval) && varstofix[c] == NOFIX )
         {
            varstofix[c] = FIXATLB;
            (*npossiblefixings)++;
         }
      }
   }

   SCIPfreeBufferArray(scip, &ubdual);
   SCIPfreeBufferArray(scip, &lbdual);
   SCIPfreeBufferArray(scip, &mincolactneginf);
   SCIPfreeBufferArray(scip, &mincolactposinf);
   SCIPfreeBufferArray(scip, &maxcolactneginf);
   SCIPfreeBufferArray(scip, &maxcolactposinf);
   SCIPfreeBufferArray(scip, &maxcolact);
   SCIPfreeBufferArray(scip, &mincolact);
   SCIPfreeBufferArray(scip, &excludevar);
   SCIPfreeBufferArray(scip, &varissingcol);

   return SCIP_OKAY;
}

/** initialize row relations */
static
void initRowRelation(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   RELATIONTYPE*         rowrelation         /**< row relation type */
   )
{
   int r;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(rowrelation != NULL);

   for( r = 0; r < matrix->nrows; r++ )
   {
      assert(!SCIPisInfinity(scip, -matrix->lhs[r]) || !SCIPisInfinity(scip, matrix->rhs[r]));

      if( SCIPisInfinity(scip, -matrix->lhs[r]) )
      {
         /* <= */
         rowrelation[r] = LESSEQUAL;
      }
      else if( SCIPisInfinity(scip, matrix->rhs[r]) )
      {
         /* >= */
         rowrelation[r] = GREATEREQUAL;
      }
      else if( SCIPisEQ(scip, matrix->lhs[r], matrix->rhs[r]) )
      {
         /* = */
         rowrelation[r] = EQUALITY;
      }
      else
      {
         /* ranged */
         rowrelation[r] = RANGED;
      }
   }
}


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyDualinfer)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDualinfer(scip) );

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDualinfer)
{  /*lint --e{715}*/
   CONSTRAINTMATRIX* matrix;
   SCIP_Bool initialized;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( (SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING) || SCIPinProbing(scip) )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* initialize constraint matrix */
   matrix = NULL;
   initialized = FALSE;
   SCIP_CALL( initMatrix(scip, &matrix, &initialized) );

   if( initialized )
   {
      FIXINGDIRECTION* varstofix;
      RELATIONTYPE* rowrelation;
      int nconvarsfixed;
      int nintvarsfixed;
      int nbinvarsfixed;
      int npossiblefixings;
      int nfitsinglecols;
      int i;

      nconvarsfixed = 0;
      nintvarsfixed = 0;
      nbinvarsfixed = 0;
      npossiblefixings = 0;
      nfitsinglecols = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &rowrelation, matrix->nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varstofix, matrix->ncols) );

      /* determine row relation types */
      initRowRelation(scip, matrix, rowrelation);

      /* initialize varstofix to 0 (= NOFIX) */
      BMSclearMemoryArray(varstofix, matrix->ncols);

      /* dual cost fixing */
      SCIP_CALL( costFixing(scip, matrix, rowrelation, varstofix, &nfitsinglecols, &npossiblefixings) );

      if( SCIPgetNContVars(scip) > 1 )
      {
         /* dual bound strengthening */
         SCIP_CALL( dualBoundStrengthening(scip, matrix, rowrelation, varstofix, &nfitsinglecols, &npossiblefixings) );
      }

      if( npossiblefixings > 0 )
      {
         /* look for fixable variables */
         for( i = matrix->ncols - 1; i >= 0; --i )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            SCIP_VAR* var;

            if( varstofix[i] == FIXATLB )
            {
               SCIP_Real lb;

               var = matrix->vars[i];
               lb = SCIPvarGetLbLocal(var);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPisFeasIntegral(scip, lb));

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
               *result = SCIP_SUCCESS;

               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  nconvarsfixed++;
               else if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
                  nbinvarsfixed++;
               else
                  nintvarsfixed++;
            }
            else if( varstofix[i] == FIXATUB )
            {
               SCIP_Real ub;

               var = matrix->vars[i];
               ub = SCIPvarGetUbLocal(var);
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS || SCIPisFeasIntegral(scip, ub));

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
               *result = SCIP_SUCCESS;

               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                  nconvarsfixed++;
               else if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
                  nbinvarsfixed++;
               else
                  nintvarsfixed++;
            }
         }
      }

      SCIPfreeBufferArray(scip, &varstofix);
      SCIPfreeBufferArray(scip, &rowrelation);

      if( (*result) != SCIP_CUTOFF && (nconvarsfixed + nintvarsfixed + nbinvarsfixed) > 0 )
      {
         SCIPdebugMessage("### %d vars [%d column singletons] ===>>> fixed [cont: %d, int: %d, bin: %d]\n",
            matrix->ncols, nfitsinglecols, nconvarsfixed, nintvarsfixed, nbinvarsfixed);
      }
   }

   freeMatrix(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the dual inference presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolDualinfer(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   presoldata = NULL;

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_DELAY, presolExecDualinfer, presoldata) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDualinfer) );

   return SCIP_OKAY;
}
