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
 * @todo Also run on general CIPs, if the number of locks of the investigated variables comes only from (upgraded)
 * linear constraints.
 *
 * @todo Instead of choosing variables for comparison out of one row, we should try to use 'hashvalues' for columns that
 *       indicate in which constraint type (<=, >=, or ranged row / ==) they are existing. Then sort the variables (and
 *       corresponding data) after the ranged row/equation hashvalue and only try to derive dominance on variables with
 *       the same hashvalue on ranged row/equation and fitting hashvalues for the other constraint types.
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

#include "presol_domcol.h"

#define PRESOL_NAME            "domcol"
#define PRESOL_DESC            "dominated column presolver"
#define PRESOL_PRIORITY         20000000     /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1     /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY                TRUE     /**< should presolver be delayed, if other presolvers found reductions? */

#define DEFAULT_MAXPAIRS         1000000     /**< maximal number of pair comparisons for at least one variable fixing */

/*
 * Data structures
 */

/** control parameters */
struct SCIP_PresolData
{
   int                   maxpairs;           /**< maximal number of pair comparisons for at least one variable fixing */
};

/** type of fixing direction */
enum Fixingdirection
{
   FIXATLB = -1,         /**< fix variable at lower bound */
   NOFIX   =  0,         /**< do not fix variable */
   FIXATUB =  1          /**< fix variable at upper bound */
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
#ifdef SCIP_DEBUG
   const char**          rowname;            /**< name of row */
#endif
   int                   nrows;              /**< complete number of rows */
   SCIP_Real*            lhs;                /**< left hand side per row */
   SCIP_Real*            rhs;                /**< right hand side per row */
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

   assert(vars != NULL);
   assert(vals != NULL);

   rowidx = matrix->nrows;

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
   }
   assert(!SCIPisInfinity(scip, -matrix->lhs[rowidx]));


#ifdef SCIP_DEBUG
   matrix->rowname[rowidx] = name;
#endif

   matrix->rowmatbeg[rowidx] = matrix->nnonzs;

   for( j = 0; j < nvars; j++ )
   {
      matrix->rowmatval[matrix->nnonzs] = factor * vals[j];
      probindex = SCIPvarGetProbindex(vars[j]);
      assert(matrix->vars[probindex] == vars[j]);

      assert(0 <= probindex && probindex < matrix->ncols);
      matrix->rowmatind[matrix->nnonzs] = probindex;
      (matrix->nnonzs)++;
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

#ifdef SCIP_DEBUG
/** print a row from the constraint matrix */
static
void printRow(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   row                 /**< row index for printing */
   )
{
   int* rowpnt;
   int* rowend;
   int c;
   SCIP_Real val;
   SCIP_Real* valpnt;
   char relation;

   relation='-';
   if( !matrix->isrhsinfinite &&
      SCIPisEQ(scip, matrix->lhs[row], matrix->rhs[row]))
   {
      relation='e';
   }
   else if( !matrix->isrhsinfinite &&
      !SCIPisEQ(scip, matrix->lhs[row], matrix->rhs[row]))
   {
      relation='r';
   }
   else
   {
      relation='g';
   }

   rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
   rowend = rowpnt + matrix->rowmatcnt[row];
   valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

   SCIPdebugPrintf("\n(L:%g) [%c] %g  <=", (matrix->minactivityposinf[row] + matrix->minactivityneginf[row] > 0) ? -SCIPinfinity(scip) : matrix->minactivity[row], relation, matrix->lhs[row]);
   for(; (rowpnt < rowend); rowpnt++, valpnt++)
   {
      c = *rowpnt;
      val = *valpnt;
      SCIPdebugPrintf("  %g{%s[idx:%d][bnd:%g,%g]}",
         val, SCIPvarGetName(matrix->vars[c]), c,
         SCIPvarGetLbGlobal(matrix->vars[c]), SCIPvarGetUbGlobal(matrix->vars[c]));
   }
   SCIPdebugPrintf(" <=  %g (U:%g)", (matrix->maxactivityposinf[row] + matrix->maxactivityneginf[row] > 0) ? SCIPinfinity(scip) : matrix->rhs[row], matrix->maxactivity[row]);
}

/** print all rows from the constraint matrix containing a variable */
static
SCIP_RETCODE printRowsOfCol(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int                   col                 /**< column index for printing */
   )
{
   int numrows;
   int* rows;
   int i;
   int* colpnt;
   int* colend;

   numrows = matrix->colmatcnt[col];

   SCIP_CALL( SCIPallocBufferArray(scip, &rows, numrows) );

   colpnt = matrix->colmatind + matrix->colmatbeg[col];
   colend = colpnt + matrix->colmatcnt[col];
   for( i = 0; (colpnt < colend); colpnt++, i++ )
   {
      rows[i] = *colpnt;
   }

   SCIPdebugPrintf("\n-------");
   SCIPdebugPrintf("\ncol %d number rows: %d",col,numrows);
   for( i = 0; i < numrows; i++ )
   {
      printRow(scip, matrix, rows[i]);
   }
   SCIPdebugPrintf("\n-------");

   SCIPfreeBufferArray(scip, &rows);

   return SCIP_OKAY;
}

/** print information about a dominance relation */
static
SCIP_RETCODE printDomRelInfo(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   SCIP_VAR*             dominatingvar,      /**< dominating variable */
   int                   dominatingidx,      /**< index of dominating variable */
   SCIP_VAR*             dominatedvar,       /**< dominated variable */
   int                   dominatedidx,       /**< index of dominated variable */
   SCIP_Real             dominatingub,       /**< predicted upper bound of dominating variable */
   SCIP_Real             dominatingwclb      /**< worst case lower bound of dominating variable */
   )
{
   char type;

   assert(SCIPvarGetType(dominatingvar)==SCIPvarGetType(dominatedvar));

   switch(SCIPvarGetType(dominatingvar))
   {
   case SCIP_VARTYPE_CONTINUOUS:
      type='C';
      break;
   case SCIP_VARTYPE_BINARY:
      type='B';
      break;
   case SCIP_VARTYPE_INTEGER:
   case SCIP_VARTYPE_IMPLINT:
      type='I';
      break;
   default:
      SCIPerrorMessage("unknown variable type\n");
      SCIPABORT();
      return SCIP_INVALIDDATA; /*lint !e527*/
   }

   SCIPdebugPrintf("\n\n### [%c], obj:%g->%g,\t%s[idx:%d](nrows:%d)->%s[idx:%d](nrows:%d)\twclb=%g, ub'=%g, ub=%g",
      type, SCIPvarGetObj(dominatingvar), SCIPvarGetObj(dominatedvar),
      SCIPvarGetName(dominatingvar), dominatingidx, matrix->colmatcnt[dominatingidx],
      SCIPvarGetName(dominatedvar), dominatedidx, matrix->colmatcnt[dominatedidx],
      dominatingwclb, dominatingub, SCIPvarGetUbGlobal(dominatingvar));

   SCIP_CALL( printRowsOfCol(scip, matrix, dominatingidx) );

   return SCIP_OKAY;
}
#endif

#ifndef NDEBUG
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

   for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
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
#endif

/** get minimum/maximum residual activity without the specified variable */
static
void getActivityResiduals(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   SCIP_VAR*             var,                /**< variable for activity residual calculation */
   SCIP_Real             val,                /**< coefficient of this variable in this row */
   int                   row,                /**< row index */
   int                   col,                /**< column index */
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

   lb = SCIPvarGetLbGlobal(var);
   ub = SCIPvarGetUbGlobal(var);

   assert(!SCIPisInfinity(scip, lb));
   assert(!SCIPisInfinity(scip, -ub));

   if( val >= 0.0 )
   {
      if( SCIPisInfinity(scip, ub) )
      {
         assert(matrix->maxactivityposinf[row] >= 1);
         assert(matrix->maxactivityneginf[row] >= 0);
         if( matrix->maxactivityposinf[row] == 1 && matrix->maxactivityneginf[row] == 0 )
         {
            *maxresactivity = matrix->maxactivity[row];
            assert(SCIPisFeasEQ(scip, *maxresactivity, getMaxActSingleRowWithoutCol(scip, matrix, row, col)));
         }
         else
            *maxresactivity = SCIPinfinity(scip);
      }
      else
      {
         if( (matrix->maxactivityneginf[row] + matrix->maxactivityposinf[row]) > 0 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = matrix->maxactivity[row] - val * ub;
      }

      if( SCIPisInfinity(scip, -lb) )
      {
         assert(matrix->minactivityneginf[row] >= 1);
         assert(matrix->minactivityposinf[row] >= 0);
         if( matrix->minactivityneginf[row] == 1 && matrix->minactivityposinf[row] == 0 )
         {
            *minresactivity = matrix->minactivity[row];
            assert(SCIPisFeasEQ(scip, *minresactivity, getMinActSingleRowWithoutCol(scip, matrix, row, col)));
         }
         else
            *minresactivity = -SCIPinfinity(scip);
      }
      else
      {
         if( (matrix->minactivityneginf[row] + matrix->minactivityposinf[row]) > 0 )
            *minresactivity = -SCIPinfinity(scip);
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
         {
            *maxresactivity = matrix->maxactivity[row];
            assert(SCIPisFeasEQ(scip, *maxresactivity, getMaxActSingleRowWithoutCol(scip, matrix, row, col)));
         }

         else
            *maxresactivity = SCIPinfinity(scip);
      }
      else
      {
         if( (matrix->maxactivityneginf[row] + matrix->maxactivityposinf[row]) > 0 )
            *maxresactivity = SCIPinfinity(scip);
         else
            *maxresactivity = matrix->maxactivity[row] - val * lb;
      }

      if( SCIPisInfinity(scip, ub) )
      {
         assert(matrix->minactivityposinf[row] >= 1);
         assert(matrix->minactivityneginf[row] >= 0);
         if( matrix->minactivityposinf[row] == 1 && matrix->minactivityneginf[row] == 0 )
         {
            *minresactivity = matrix->minactivity[row];
            assert(SCIPisFeasEQ(scip, *minresactivity, getMinActSingleRowWithoutCol(scip, matrix, row, col)));
         }
         else
            *minresactivity = -SCIPinfinity(scip);
      }
      else
      {
         if( (matrix->minactivityneginf[row] + matrix->minactivityposinf[row]) > 0 )
            *minresactivity = -SCIPinfinity(scip);
         else
            *minresactivity = matrix->minactivity[row] - val * ub;
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
   SCIP_Real*            calculatedwclb,     /**< predicted worst case lower bound */
   SCIP_Bool*            lbcalculated,       /**< was a lower bound calculated? */
   SCIP_Real*            calculatedlb,       /**< predicted lower bound */
   SCIP_Bool*            wcubcalculated,     /**< was a worst case upper bound calculated? */
   SCIP_Real*            calculatedwcub      /**< calculated worst case upper bound */
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
   assert(lbcalculated != NULL);
   assert(calculatedlb != NULL);
   assert(wcubcalculated != NULL);
   assert(calculatedwcub != NULL);

   assert(!SCIPisZero(scip, valdominating));
   assert(matrix->vars[coldominating] != NULL);

   vardominating = matrix->vars[coldominating];
   vardominated = matrix->vars[coldominated];

   *ubcalculated = FALSE;
   *wclbcalculated = FALSE;
   *lbcalculated = FALSE;
   *wcubcalculated = FALSE;

   /* no rowbound analysis for multiaggregated variables, which should not exist, because the matrix only consists of
    * active variables
    */
   assert(SCIPvarGetStatus(vardominating) != SCIP_VARSTATUS_MULTAGGR && SCIPvarGetStatus(vardominated) != SCIP_VARSTATUS_MULTAGGR);

   lhs = matrix->lhs[row];
   rhs = matrix->rhs[row];
   assert(!SCIPisInfinity(scip, lhs));
   assert(!SCIPisInfinity(scip, -rhs));

   /* get minimum/maximum activity of this row without the dominating variable */
   getActivityResiduals(scip, matrix, vardominating, valdominating,
      row, coldominating, &minresactivity, &maxresactivity);

   assert(!SCIPisInfinity(scip, minresactivity));
   assert(!SCIPisInfinity(scip, -maxresactivity));

   lbdominated = SCIPvarGetLbGlobal(vardominated);
   ubdominated = SCIPvarGetUbGlobal(vardominated);

   assert(!SCIPisInfinity(scip, lbdominated));
   assert(!SCIPisInfinity(scip, -ubdominated));

   *calculatedub = SCIPinfinity(scip);
   *calculatedwclb = -SCIPinfinity(scip);
   *calculatedlb = -SCIPinfinity(scip);
   *calculatedwcub = SCIPinfinity(scip);

   /* predictive rowbound analysis */

   if( valdominating > 0.0 )
   {
      /* upper bound calculation */
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* <=, = or ranged */
         if( !SCIPisInfinity(scip, -minresactivity) )
         {
            if( SCIPisGE(scip, valdominated, 0.0) )
               *calculatedub = (rhs - minresactivity)/valdominating;
            else
            {
               assert(!SCIPisInfinity(scip, ubdominated));
               if( !SCIPisInfinity(scip, -lbdominated) )
                  *calculatedub = (rhs - (minresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
               else
                  return SCIP_OKAY;
            }

            *ubcalculated = TRUE;
         }
      }

      /* wort case calculation of lower bound */
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* >=, = or ranged */
         if( !SCIPisInfinity(scip, -minresactivity) )
         {
            if( SCIPisGE(scip, valdominated, 0.0) )
               *calculatedwclb = (lhs - minresactivity)/valdominating;
            else
            {
               assert(!SCIPisInfinity(scip, ubdominated));
               if( !SCIPisInfinity(scip, -lbdominated) )
                  *calculatedwclb = (lhs - (minresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
               else
                  return SCIP_OKAY;
            }

            *wclbcalculated = TRUE;
         }
         else
         {
            /* worst case lower bound is infinity */
            *calculatedwclb = SCIPinfinity(scip);
            *wclbcalculated = TRUE;
         }
      }

      /* lower bound calculation */
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* >=, = or ranged */
         if( !SCIPisInfinity(scip, maxresactivity) )
         {
            if( SCIPisLE(scip, valdominated, 0.0) )
               *calculatedlb = (lhs - maxresactivity)/valdominating;
            else
            {
               assert(!SCIPisInfinity(scip, ubdominated));
               if( !SCIPisInfinity(scip, -lbdominated) )
                  *calculatedlb = (lhs - (maxresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
               else
                  return SCIP_OKAY;
            }

            *lbcalculated = TRUE;
         }
      }

      /* worst case calculation of upper bound */
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* <=, = or ranged */
         if( !SCIPisInfinity(scip, maxresactivity) )
         {
            if( SCIPisLE(scip, valdominated, 0.0) )
               *calculatedwcub = (rhs - maxresactivity)/valdominating;
            else
            {
               assert(!SCIPisInfinity(scip, ubdominated));
               if( !SCIPisInfinity(scip, -lbdominated) )
                  *calculatedwcub = (rhs - (maxresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
               else
                  return SCIP_OKAY;
            }

            *wcubcalculated = TRUE;
         }
         else
         {
            /* worst case upper bound is -infinity */
            *calculatedwcub = -SCIPinfinity(scip);
            *wcubcalculated = TRUE;
         }
      }
   }
   else
   {
      /* upper bound calculation */
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* >=, = or ranged */
         if( !SCIPisInfinity(scip, maxresactivity) )
         {
            if( SCIPisLE(scip, valdominated, 0.0) )
               *calculatedub = (lhs - maxresactivity)/valdominating;
            else
            {
               assert(!SCIPisInfinity(scip, ubdominated));
               if( !SCIPisInfinity(scip, -lbdominated) )
                  *calculatedub = (lhs - (maxresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
               else
                  return SCIP_OKAY;
            }

            *ubcalculated = TRUE;
         }
      }

      /* worst case calculation of lower bound */
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* <=, = or ranged */
         if( !SCIPisInfinity(scip, maxresactivity) )
         {
            if( SCIPisLE(scip, valdominated, 0.0) )
               *calculatedwclb = (rhs - maxresactivity)/valdominating;
            else
            {
               assert(!SCIPisInfinity(scip, ubdominated));
               if( !SCIPisInfinity(scip, -lbdominated) )
                  *calculatedwclb = (rhs - (maxresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
               else
                  return SCIP_OKAY;
            }

            *wclbcalculated = TRUE;
         }
         else
         {
            /* worst case lower bound is infinity */
            *calculatedwclb = SCIPinfinity(scip);
            *wclbcalculated = TRUE;
         }
      }

      /* lower bound calculation */
      if( !SCIPisInfinity(scip, rhs) )
      {
         /* >=, = or ranged */
         if( !SCIPisInfinity(scip, -minresactivity) )
         {
            if( SCIPisGE(scip, valdominated, 0.0) )
               *calculatedlb = (rhs - minresactivity)/valdominating;
            else
            {
               assert(!SCIPisInfinity(scip, ubdominated));
               if( !SCIPisInfinity(scip, -lbdominated) )
                  *calculatedlb = (rhs - (minresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
               else
                  return SCIP_OKAY;
            }

            *lbcalculated = TRUE;
         }
      }

      /* worst case calculation of upper bound */
      if( !SCIPisInfinity(scip, -lhs) )
      {
         /* >=, = or ranged */
         if( !SCIPisInfinity(scip, -minresactivity) )
         {
            if( SCIPisGE(scip, valdominated, 0.0) )
               *calculatedwcub = (lhs - minresactivity)/valdominating;
            else
            {
               assert(!SCIPisInfinity(scip, ubdominated));
               if( !SCIPisInfinity(scip, -lbdominated) )
                  *calculatedwcub = (lhs - (minresactivity - (valdominated * ubdominated) + (valdominated * lbdominated)))/valdominating;
               else
                  return SCIP_OKAY;
            }

            *wcubcalculated = TRUE;
         }
         else
         {
            /* worst case upper bound is -infinity */
            *calculatedwcub = -SCIPinfinity(scip);
            *wcubcalculated = TRUE;
         }
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
   SCIP_Real*            wclowerbound,       /**< predicted worst case lower bound */
   SCIP_Real*            lowerbound,         /**< predicted lower bound */
   SCIP_Real*            wcupperbound        /**< predicted worst case upper bound */
   )
{
   SCIP_Bool ubcalculated;
   SCIP_Bool wclbcalculated;
   SCIP_Bool lbcalculated;
   SCIP_Bool wcubcalculated;
   SCIP_Real newub;
   SCIP_Real newwclb;
   SCIP_Real newlb;
   SCIP_Real newwcub;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(row < matrix->nrows);
   assert(col1 < matrix->ncols);
   assert(col2 < matrix->ncols);
   assert(upperbound != NULL);
   assert(wclowerbound != NULL);
   assert(lowerbound != NULL);
   assert(wcupperbound != NULL);

   /* do predictive rowbound analysis */
   SCIP_CALL( calcVarBounds(scip, matrix, row, col1, val1, col2, val2,
         &ubcalculated, &newub, &wclbcalculated, &newwclb,
         &lbcalculated, &newlb, &wcubcalculated, &newwcub) );

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
   if( lbcalculated )
   {
      if( newlb > *lowerbound )
         *lowerbound = newlb;
   }
   if( wcubcalculated )
   {
      if( newwcub < *wcupperbound )
         *wcupperbound = newwcub;
   }

   return SCIP_OKAY;
}

/** detect parallel columns by using the algorithm of Bixby and Wagner
 *  see paper: "A note on Detecting Simple Redundancies in Linear Systems", June 1986
 */
static
SCIP_RETCODE detectParallelCols(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   int*                  pclass,             /**< parallel column classes */
   SCIP_Bool*            varineq             /**< indicating if variable is within an equation */
   )
{
   SCIP_Real* valpnt;
   SCIP_Real* values;
   SCIP_Real* scale;
   int* classsizes;
   int* pcset;
   int* rowpnt;
   int* rowend;
   int* colindices;
   int* pcs;
   SCIP_Real startval;
   SCIP_Real aij;
   SCIP_Bool isequality;
   int startpc;
   int startk;
   int startt;
   int pcsetfill;
   int colidx;
   int k;
   int t;
   int m;
   int i;
   int r;
   int newpclass;
   int pc;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(pclass != NULL);
   assert(varineq != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &classsizes, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scale, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcset, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &values, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colindices, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &pcs, matrix->ncols) );

   /* init */
   BMSclearMemoryArray(scale, matrix->ncols);
   BMSclearMemoryArray(pclass, matrix->ncols);
   BMSclearMemoryArray(classsizes, matrix->ncols);
   classsizes[0] = matrix->ncols;
   pcsetfill = 0;
   for( t = 1; t < matrix->ncols; ++t )
      pcset[pcsetfill++] = t;

   /* loop over all rows */
   for( r = 0; r < matrix->nrows; ++r )
   {
      isequality = !matrix->isrhsinfinite[r];

      /* we consider only equations or ranged rows */
      if( isequality )
      {
         rowpnt = matrix->rowmatind + matrix->rowmatbeg[r];
         rowend = rowpnt + matrix->rowmatcnt[r];
         valpnt = matrix->rowmatval + matrix->rowmatbeg[r];

         i = 0;
         for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
         {
            aij = *valpnt;
            colidx = *rowpnt;

#ifdef SCIP_DEBUG
            if( SCIPisZero(scip, aij) )
            {
               SCIPdebugMessage("Matrix coefficient is very small !\n");
            }
#endif
            /* remember variable was within an equation or ranged row present */
            varineq[colidx] = isequality;

            if( scale[colidx] == 0.0 )
               scale[colidx] = aij;
            assert(scale[colidx] != 0.0);

            colindices[i] = colidx;
            values[i] = aij / scale[colidx];
            pc = pclass[colidx];
            assert(pc < matrix->ncols);

            /* update class sizes and pclass set */
            assert(classsizes[pc] > 0);
            classsizes[pc]--;
            if( classsizes[pc] == 0 )
            {
               assert(pcsetfill < matrix->ncols);
               pcset[pcsetfill++] = pc;
            }
            pcs[i] = pc;

            i++;
         }

         /* sort on the pclass values */
         if( i > 1 )
            SCIPsortIntIntReal(pcs, colindices, values, i);

         k = 0;
         while( TRUE ) /*lint !e716*/
         {
            assert(k < i);
            startpc = pcs[k];
            startk = k;

            /* find pclass-sets */
            while( k < i && pcs[k] == startpc )
               k++;

            /* sort on the A values which have equal pclass values */
            if( k - startk > 1 )
               SCIPsortRealInt(&(values[startk]), &(colindices[startk]), k - startk);

            t = 0;
            while( TRUE ) /*lint !e716*/
            {
               assert(startk + t < i);
               startval = values[startk + t];
               startt = t;

               /* find A-sets */
               while( t < k - startk && SCIPisEQ(scip, startval, values[startk + t]) )
                  t++;

               /* get new pclass */
               newpclass = pcset[0];
               assert(pcsetfill > 0);
               pcset[0] = pcset[--pcsetfill];

               /* renumbering */
               for( m = startk + startt; m < startk + t; m++ )
               {
                  assert(m < i);
                  assert(colindices[m] < matrix->ncols);
                  assert(newpclass < matrix->ncols);

                  pclass[colindices[m]] = newpclass;
                  classsizes[newpclass]++;
               }

               if( t == k - startk )
                  break;
            }

            if( k == matrix->rowmatcnt[r] )
               break;
         }
      }
   }

   SCIPfreeBufferArray(scip, &pcs);
   SCIPfreeBufferArray(scip, &colindices);
   SCIPfreeBufferArray(scip, &values);
   SCIPfreeBufferArray(scip, &pcset);
   SCIPfreeBufferArray(scip, &scale);
   SCIPfreeBufferArray(scip, &classsizes);

   return SCIP_OKAY;
}


/** try to find possible variable fixings */
static
void findFixings(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix structure */
   SCIP_VAR*             dominatingvar,      /**< dominating variable */
   int                   dominatingidx,      /**< column index of the dominating variable */
   SCIP_Real             dominatingub,       /**< predicted upper bound of the dominating variable */
   SCIP_Real             dominatingwclb,     /**< predicted worst case lower bound of the dominating variable */
   SCIP_Real             dominatinglb,       /**< predicted lower bound of the dominating variable */
   SCIP_Real             dominatingwcub,     /**< predicted worst case upper bound of the dominating variable */
   SCIP_VAR*             dominatedvar,       /**< dominated variable */
   int                   dominatedidx,       /**< column index of the dominated variable */
   FIXINGDIRECTION*      varstofix,          /**< array holding fixing information */
   SCIP_Bool             onlybinvars,        /**< flag indicating only binary variables are present */
   SCIP_Bool             onlyoneone,         /**< when onlybinvars is TRUE, flag indicates if both binary variables are in clique */
   int*                  npossiblefixings    /**< counter for possible fixings */
   )
{
   /* we compare only variables from the same type */
   assert(SCIPvarGetType(dominatingvar) == SCIPvarGetType(dominatedvar) ||
      SCIPvarIsBinary(dominatingvar) == SCIPvarIsBinary(dominatedvar) ||
      (SCIPvarGetType(dominatingvar) == SCIP_VARTYPE_INTEGER && SCIPvarGetType(dominatedvar) == SCIP_VARTYPE_IMPLINT) ||
      (SCIPvarGetType(dominatedvar) == SCIP_VARTYPE_INTEGER && SCIPvarGetType(dominatingvar) == SCIP_VARTYPE_IMPLINT));

   if( matrix->colmatcnt[dominatingidx] == 1 && matrix->colmatcnt[dominatedidx] == 1 )
   {
      /* we have a x->y dominance relation and only one equality constraint
       * where the dominating variable x with a infinity upper bound and the
       * dominated variable y are present, then the dominated variable y
       * can be fixed at their lower bound.
       */
      int row;
      row = *(matrix->colmatind + matrix->colmatbeg[dominatedidx]);

      if( SCIPisEQ(scip, matrix->lhs[row], matrix->rhs[row]) &&
         SCIPisInfinity(scip, SCIPvarGetUbGlobal(dominatingvar)) )
      {
         if( varstofix[dominatedidx] == NOFIX )
         {
            varstofix[dominatedidx] = FIXATLB;
            (*npossiblefixings)++;
#ifdef SCIP_DEBUG
            SCIP_CALL_ABORT( printDomRelInfo(scip,matrix,dominatingvar,dominatingidx,
                  dominatedvar,dominatedidx,dominatingub,dominatingwclb) );
#endif
            return;
         }
      }
   }

   if( SCIPisPositive(scip, SCIPvarGetObj(dominatedvar)) )
   {
      if( !SCIPisInfinity(scip, -dominatingwclb) &&
         SCIPisLE(scip, dominatingwclb, SCIPvarGetUbGlobal(dominatingvar)) )
      {
         /* we have a x->y dominance relation with a positive obj coefficient
          * of the dominated variable y. we need to secure feasibility
          * by testing if the predicted lower worst case bound is less equal the
          * current upper bound. it is possible, that the lower worst case bound
          * is infinity and the upper bound of the dominating variable x is
          * infinity too.
          */
         if( varstofix[dominatedidx] == NOFIX )
         {
            varstofix[dominatedidx] = FIXATLB;
            (*npossiblefixings)++;
#ifdef SCIP_DEBUG
            SCIP_CALL_ABORT( printDomRelInfo(scip, matrix, dominatingvar, dominatingidx,
                  dominatedvar, dominatedidx, dominatingub, dominatingwclb) );
#endif
         }
      }
   }

   if( !SCIPisInfinity(scip, dominatingub) &&
      SCIPisLE(scip, dominatingub, SCIPvarGetUbGlobal(dominatingvar)) )
   {
      /* we have a x->y dominance relation with an arbitrary obj coefficient
       * of the dominating variable x. in all cases we have to look
       * if the predicted upper bound of the dominating variable is great enough.
       * by testing, that the predicted upper bound is not infinity we avoid problems
       * with x->y e.g.
       *    min  -x -y
       *    s.t. -x -y <= -1
       *    0<=x<=1, 0<=y<=1
       * where y is not at their lower bound.
       */
      if( varstofix[dominatedidx] == NOFIX )
      {
         varstofix[dominatedidx] = FIXATLB;
         (*npossiblefixings)++;
#ifdef SCIP_DEBUG
         SCIP_CALL_ABORT( printDomRelInfo(scip, matrix, dominatingvar, dominatingidx,
               dominatedvar, dominatedidx, dominatingub, dominatingwclb) );
#endif
      }
   }

   if( SCIPisNegative(scip, SCIPvarGetObj(dominatingvar)) )
   {
      /* we have a x->y dominance relation with a negative obj coefficient
       * of the dominating variable x. if the worst case upper bound is
       * greater equal than upper bound, we fix x at the upper bound
       */
      if( !SCIPisInfinity(scip, dominatingwcub) &&
         SCIPisGE(scip, dominatingwcub, SCIPvarGetUbGlobal(dominatingvar)) )
      {
         if( varstofix[dominatingidx] == NOFIX )
         {
            varstofix[dominatingidx] = FIXATUB;
            (*npossiblefixings)++;
         }
      }
   }

   if( !SCIPisInfinity(scip, -dominatinglb) &&
      SCIPisGE(scip, dominatinglb, SCIPvarGetUbGlobal(dominatingvar)) )
   {
       /* we have a x->y dominance relation with an arbitrary obj coefficient
        * of the dominating variable x. if the predicted lower bound is greater
        * equal than upper bound, we fix x at the upper bound.
        */
      if( varstofix[dominatingidx] == NOFIX )
      {
         varstofix[dominatingidx] = FIXATUB;
         (*npossiblefixings)++;
      }
   }

   if( onlybinvars )
   {
      if( varstofix[dominatedidx] == NOFIX && (onlyoneone || SCIPvarsHaveCommonClique(dominatingvar, TRUE, dominatedvar, TRUE, TRUE)) )
      {
         /* We have a (1->1)-clique with dominance relation (x->y) (x dominates y).
          * From this dominance relation, we know (1->0) is possible and not worse than (0->1)
          * concerning the objective function. It follows that only (1->0) or (0->0) are possible,
          * but in both cases y has the value 0 => y=0.
          */
         varstofix[dominatedidx] = FIXATLB;
         (*npossiblefixings)++;
      }

      if( varstofix[dominatingidx] == NOFIX && SCIPvarsHaveCommonClique(dominatingvar, FALSE, dominatedvar, FALSE, TRUE) )
      {
         /* We have a (0->0)-clique with dominance relation x->y (x dominates y).
          * From this dominance relation, we know (1->0) is possible and not worse than (0->1)
          * concerning the objective function. It follows that only (1->0) or (1->1) are possible,
          * but in both cases x has the value 1 => x=1
          */
         varstofix[dominatingidx] = FIXATUB;
         (*npossiblefixings)++;
      }
   }
   else
      assert(!onlyoneone);

}

/** find dominance relation between variable pairs */
static
SCIP_RETCODE findDominancePairs(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   SCIP_PRESOLDATA*      presoldata,         /**< presolver data */
   int*                  searchcols,         /**< indexes of variables for pair comparisons */
   int                   searchsize,         /**< number of variables for pair comparisons */
   SCIP_Bool             onlybinvars,        /**< flag indicating searchcols has only binary variable indexes */
   FIXINGDIRECTION*      varstofix,          /**< array holding information for later upper/lower bound fixing */
   int*                  npossiblefixings,   /**< found number of possible fixings */
   SCIP_Longint*         ndomrelations       /**< found number of dominance relations */
   )
{
   SCIP_Real* vals1;
   SCIP_Real* vals2;
   SCIP_Real tmpupperboundcol1;
   SCIP_Real tmpupperboundcol2;
   SCIP_Real tmpwclowerboundcol1;
   SCIP_Real tmpwclowerboundcol2;
   SCIP_Real tmplowerboundcol1;
   SCIP_Real tmplowerboundcol2;
   SCIP_Real tmpwcupperboundcol1;
   SCIP_Real tmpwcupperboundcol2;
   int* rows1;
   int* rows2;
   int nrows1;
   int nrows2;
   SCIP_Bool col1domcol2;
   SCIP_Bool col2domcol1;
   SCIP_Bool onlyoneone;
   int cnt1;
   int cnt2;
   int col1;
   int col2;
   int r1;
   int r2;
   int paircnt;
   int nlastpossiblefixings;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(presoldata != NULL);
   assert(searchcols != NULL);
   assert(varstofix != NULL);
   assert(npossiblefixings != NULL);
   assert(ndomrelations != NULL);

   paircnt = 0;
   nlastpossiblefixings = *npossiblefixings;

   /* pair comparisons */
   for( cnt1 = 0; cnt1 < searchsize; cnt1++ )
   {
      for( cnt2 = cnt1+1; cnt2 < searchsize; cnt2++ )
      {
         /* get indexes of this variable pair */
         col1 = searchcols[cnt1];
         col2 = searchcols[cnt2];

         onlyoneone = FALSE;

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
            if( SCIPisLT(scip, SCIPvarGetLbGlobal(matrix->vars[col1]), 0.0) ||
               SCIPisLT(scip, SCIPvarGetLbGlobal(matrix->vars[col2]), 0.0) )
            {
               col1domcol2 = FALSE;
               col2domcol1 = FALSE;
            }
         }

         if( paircnt == presoldata->maxpairs )
         {
            if( !(*npossiblefixings > nlastpossiblefixings) )
            {
               /* not enough fixings found, stop searching */
               return SCIP_OKAY;
            }
            nlastpossiblefixings = *npossiblefixings;
            paircnt = 0;
         }
         paircnt++;

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
         tmplowerboundcol1 = -SCIPinfinity(scip);
         tmplowerboundcol2 = tmplowerboundcol1;
         tmpwcupperboundcol1 = SCIPinfinity(scip);
         tmpwcupperboundcol2 = tmpwcupperboundcol1;

         /* compare rows of this column pair */
         while( (col1domcol2 || col2domcol1) && (r1 < nrows1 || r2 < nrows2))
         {
            assert((r1 >= nrows1-1) || (rows1[r1] < rows1[r1+1]));
            assert((r2 >= nrows2-1) || (rows2[r2] < rows2[r2+1]));

            /* there is a nonredundant row containing column 1 but not column 2 */
            if( r1 < nrows1 && (r2 == nrows2 || rows1[r1] < rows2[r2]) )
            {
               /* dominance depends on the type of relation */
               if( !matrix->isrhsinfinite[rows1[r1]] )
               {
                  /* no dominance relation for equations or ranged rows */
                  col2domcol1 = FALSE;
                  col1domcol2 = FALSE;
               }
               else
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals1[r1] > 0.0 )
                     col2domcol1 = FALSE;
                  else if( vals1[r1] < 0.0 )
                     col1domcol2 = FALSE;
               }

               r1++;
            }
            /* there is a nonredundant row containing column 2, but not column 1 */
            else if( r2 < nrows2 && (r1 == nrows1 || rows1[r1] > rows2[r2]) )
            {
               /* dominance depends on the type of relation */
               if( !matrix->isrhsinfinite[rows2[r2]] )
               {
                  /* no dominance relation for equations or ranged rows */
                  col2domcol1 = FALSE;
                  col1domcol2 = FALSE;
               }
               else
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals2[r2] < 0.0 )
                     col2domcol1 = FALSE;
                  else if( vals2[r2] > 0.0 )
                     col1domcol2 = FALSE;
               }

               r2++;
            }
            /* if both columns appear in a common row, compare the coefficients */
            else
            {
               assert(r1 < nrows1 && r2 < nrows2);
               assert(rows1[r1] == rows2[r2]);

               /* dominance depends on the type of inequality */
               if( !matrix->isrhsinfinite[rows1[r1]] )
               {
                  /* no dominance relation if coefficients differ for equations or ranged rows */
                  if( !SCIPisEQ(scip, vals1[r1], vals2[r2]) )
                  {
                     col2domcol1 = FALSE;
                     col1domcol2 = FALSE;
                  }

                  if( onlybinvars )
                  {
                     if( !onlyoneone && (matrix->minactivityposinf[rows1[r1]] + matrix->minactivityneginf[rows1[r1]] == 0) && SCIPisFeasGE(scip, matrix->minactivity[rows1[r1]] + MIN(vals1[r1], vals2[r2]), matrix->rhs[rows1[r1]]) )
                        onlyoneone = TRUE;
                  }
               }
               else
               {
                  /* >= relation, larger coefficients dominate smaller ones */
                  if( vals1[r1] >= vals2[r2] )
                     col2domcol1 = FALSE;
                  else
                     col1domcol2 = FALSE;
               }

               if( !onlyoneone && ((vals1[r1] < 0 && vals2[r2] < 0) || (vals1[r1] > 0 && vals2[r2] > 0)) )
               {
                  if( col1domcol2 )
                  {
                     /* update bounds for column 1 */
                     SCIP_CALL( updateBounds(scip, matrix, rows1[r1], col1, vals1[r1], col2, vals2[r2],
                           &tmpupperboundcol1, &tmpwclowerboundcol1, &tmplowerboundcol1, &tmpwcupperboundcol1) );
                  }

                  if( col2domcol1 )
                  {
                     /* update bounds for column 2 */
                     SCIP_CALL( updateBounds(scip, matrix, rows2[r2], col2, vals2[r2], col1, vals1[r1],
                           &tmpupperboundcol2, &tmpwclowerboundcol2, &tmplowerboundcol2, &tmpwcupperboundcol2) );
                  }
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
            if( !matrix->isrhsinfinite[rows1[r1]] )
            {
               col2domcol1 = FALSE;
               col1domcol2 = FALSE;
               break;
            }
            r1++;
         }
         while( r2 < nrows2 )
         {
            if( !matrix->isrhsinfinite[rows2[r2]] )
            {
               col2domcol1 = FALSE;
               col1domcol2 = FALSE;
               break;
            }
            r2++;
         }

         if( col1domcol2 || col2domcol1 )
            (*ndomrelations)++;

         if( col1domcol2 && col2domcol1 )
         {
            /* prefer the variable as dominating variable with the greater upper bound */
            if( SCIPisGE(scip, SCIPvarGetUbGlobal(matrix->vars[col1]), SCIPvarGetUbGlobal(matrix->vars[col2])) )
               col2domcol1 = FALSE;
            else
               col1domcol2 = FALSE;
         }

         /* use dominance relation and clique/bound-information to find variable fixings */
         if( col1domcol2 )
         {
            findFixings(scip, matrix, matrix->vars[col1], col1,
               tmpupperboundcol1, tmpwclowerboundcol1,
               tmplowerboundcol1, tmpwcupperboundcol1, matrix->vars[col2],
               col2, varstofix, onlybinvars, onlyoneone, npossiblefixings);
         }
         else if( col2domcol1 )
         {
            findFixings(scip, matrix, matrix->vars[col2], col2,
               tmpupperboundcol2, tmpwclowerboundcol2,
               tmplowerboundcol2, tmpwcupperboundcol2, matrix->vars[col1],
               col1, varstofix, onlybinvars, onlyoneone, npossiblefixings);
         }
      }
   }

   return SCIP_OKAY;
}

/** try to fix singleton column continuous variables */
static
SCIP_RETCODE singletonColumnStuffing(
   SCIP*                 scip,               /**< SCIP main data structure */
   CONSTRAINTMATRIX*     matrix,             /**< matrix containing the constraints */
   SCIP_Bool*            varsprocessed,      /**< array indicating that this variable has been processed */
   FIXINGDIRECTION*      varstofix,          /**< array holding fixing information */
   int*                  npossiblefixings    /**< number of possible fixings */
   )
{
   SCIP_VAR* var;
   SCIP_Real* valpnt;
   SCIP_Real* colratios;
   SCIP_Real* colcoeffs;
   SCIP_Bool* rowprocessed;
   int* rowpnt;
   int* rowend;
   int* colindices;
   int* colnozerolb;
   SCIP_Real constant;
   SCIP_Real val;
   SCIP_Real value;
   SCIP_Real boundoffset;
   SCIP_Bool tryfixing;
   int idx;
   int col;
   int row;
   int fillcnt;
   int colidx;
   int k;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(varsprocessed != NULL);
   assert(varstofix != NULL);
   assert(npossiblefixings != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &colindices, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colratios, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &colcoeffs, matrix->ncols) );

   SCIP_CALL( SCIPallocBufferArray(scip, &colnozerolb, matrix->ncols) );
   BMSclearMemoryArray(colnozerolb, matrix->ncols);

   SCIP_CALL( SCIPallocBufferArray(scip, &rowprocessed, matrix->nrows) );
   BMSclearMemoryArray(rowprocessed, matrix->nrows);

   for( col = 0; col < matrix->ncols; col++ )
   {
      /* we look only at rows with minimal one continuous singleton column */
      if( matrix->colmatcnt[col] == 1 && SCIPvarGetType(matrix->vars[col]) == SCIP_VARTYPE_CONTINUOUS )
      {
         row = *(matrix->colmatind + matrix->colmatbeg[col]);
         if( rowprocessed[row] )
            continue;

         rowprocessed[row] = TRUE;

         if( matrix->isrhsinfinite[row] )
         {
            /* singleton column pushing for >= relation */
            fillcnt = 0;
            tryfixing = TRUE;
            constant = 0.0;

            rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
            rowend = rowpnt + matrix->rowmatcnt[row];
            valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

            for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
            {
               val = *valpnt;
               colidx = *rowpnt;
               var = matrix->vars[colidx];
               assert(var != NULL);
               assert(val != 0.0);

               if( SCIPisGE(scip, SCIPvarGetLbGlobal(var), 0.0) )
               {
                  if( SCIPisNegative(scip, val) )
                  {
                     /* do we have a continuous singleton column */
                     if( matrix->colmatcnt[colidx] == 1 && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
                     {
                        if( SCIPisPositive(scip, SCIPvarGetLbGlobal(var)) )
                        {
                           constant += val * SCIPvarGetLbGlobal(var);
                           colnozerolb[fillcnt] = 1;
                        }

                        colratios[fillcnt] = SCIPvarGetObj(var) / val;
                        colindices[fillcnt] = colidx;
                        colcoeffs[fillcnt] = val;
                        fillcnt++;
                     }
                     else
                     {
                        /* discrete variables or variables in more than one row are present */
                        if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
                        {
                           tryfixing = FALSE;
                           break;
                        }
                        constant += val * SCIPvarGetUbGlobal(var);
                     }
                  }
                  else if( SCIPisPositive(scip, val) )
                  {
                     constant += val * SCIPvarGetLbGlobal(var);
                  }
               }
               else
               {
                  tryfixing = FALSE;
                  break;
               }
            }

            if( tryfixing )
            {
               SCIPsortRealRealIntInt(colratios, colcoeffs, colindices, colnozerolb, fillcnt);

               /* try to fix continuous singleton columns by their ratio */
               for( k = fillcnt - 1; k >= 0; k-- )
               {
                  boundoffset = 0.0;
                  idx = colindices[k];
                  value = colcoeffs[k] * SCIPvarGetUbGlobal(matrix->vars[idx]);

                  if( colnozerolb[k] )
                  {
                     boundoffset = colcoeffs[k] * SCIPvarGetLbGlobal(matrix->vars[idx]);
                  }

                  if( matrix->colmatcnt[idx] == 1 &&
                     SCIPvarGetType(matrix->vars[idx]) == SCIP_VARTYPE_CONTINUOUS &&
                     SCIPisNegative(scip, SCIPvarGetObj(matrix->vars[idx])) )
                  {
                     if( SCIPisGE(scip, value, matrix->lhs[row] - constant + boundoffset) )
                     {
                        constant += value;
                        varstofix[idx] = FIXATUB;
                        varsprocessed[idx] = TRUE;
                        (*npossiblefixings)++;

                        if( colnozerolb[k] )
                           constant -= boundoffset;
                     }
                     else
                        break;
                  }
               }
            }
         }


         if( matrix->isrhsinfinite[row] )
         {
            /* singleton column pulling for >= relation */
            fillcnt = 0;
            tryfixing = TRUE;
            constant = 0.0;

            rowpnt = matrix->rowmatind + matrix->rowmatbeg[row];
            rowend = rowpnt + matrix->rowmatcnt[row];
            valpnt = matrix->rowmatval + matrix->rowmatbeg[row];

            for( ; (rowpnt < rowend); rowpnt++, valpnt++ )
            {
               val = *valpnt;
               colidx = *rowpnt;
               var = matrix->vars[colidx];
               assert(var != NULL);
               assert(val != 0.0);

               if( SCIPisGE(scip, SCIPvarGetLbGlobal(var), 0.0) )
               {
                  if( SCIPisPositive(scip, val) )
                  {
                     /* do we have a continuous singleton column */
                     if( matrix->colmatcnt[colidx] == 1 && SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
                        && SCIPisPositive(scip, SCIPvarGetObj(var)) )
                     {
                        if( SCIPisPositive(scip, SCIPvarGetLbGlobal(var)) )
                        {
                           constant += val * SCIPvarGetLbGlobal(var);
                           colnozerolb[fillcnt] = 1;
                        }

                        colratios[fillcnt] = SCIPvarGetObj(var) / val;
                        colindices[fillcnt] = colidx;
                        colcoeffs[fillcnt] = val;
                        fillcnt++;
                     }
                     else
                     {
                        /* discrete variables or variables which are present within
                         * more than one row are estimated at their upper bound
                         */
                        if( SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)) )
                        {
                           tryfixing = FALSE;
                           break;
                        }
                        constant += val * SCIPvarGetUbGlobal(var);
                     }
                  }
                  else if( SCIPisNegative(scip, val) )
                  {
                     /* consider lower bound for negative coefficients */
                     constant += val * SCIPvarGetLbGlobal(var);
                  }
               }
               else
               {
                  tryfixing = FALSE;
                  break;
               }
            }

            if( tryfixing )
            {
               SCIPsortRealRealIntInt(colratios, colcoeffs, colindices, colnozerolb, fillcnt);

               /* try to fix continuous singleton columns by their ratio */
               for( k = 0; k < fillcnt; k++ )
               {
                  boundoffset = 0.0;
                  idx = colindices[k];
                  value = colcoeffs[k] * SCIPvarGetUbGlobal(matrix->vars[idx]);

                  if( colnozerolb[k] )
                  {
                     boundoffset = colcoeffs[k] * SCIPvarGetLbGlobal(matrix->vars[idx]);
                  }

                  if( matrix->colmatcnt[idx] == 1 &&
                     SCIPvarGetType(matrix->vars[idx]) == SCIP_VARTYPE_CONTINUOUS &&
                     SCIPisPositive(scip, SCIPvarGetObj(matrix->vars[idx])) )
                  {
                     if( SCIPisLE(scip, value, matrix->lhs[row] - constant + boundoffset) )
                     {
                        constant += value;
                        varstofix[idx] = FIXATUB;
                        varsprocessed[idx] = TRUE;
                        (*npossiblefixings)++;

                        if( colnozerolb[k] )
                           constant -= boundoffset;
                     }
                     else
                        break;
                  }
               }
            }
         }
      }
   }

   SCIPfreeBufferArray(scip, &rowprocessed);
   SCIPfreeBufferArray(scip, &colnozerolb);
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
SCIP_DECL_PRESOLCOPY(presolCopyDomcol)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolDomcol(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeDomcol)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecDomcol)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;
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

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* initialize constraint matrix */
   matrix = NULL;
   initialized = FALSE;
   SCIP_CALL( initMatrix(scip, &matrix, &initialized) );

   if( initialized )
   {
      int npossiblefixings;
      SCIP_Longint ndomrelations;
      int v;
      int r;
      FIXINGDIRECTION* varstofix;
      SCIP_Bool* varsprocessed;
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
#ifdef SCIP_DEBUG
      int nconvarsfixed = 0;
      int nintvarsfixed = 0;
      int nbinvarsfixed = 0;
#endif
      int* pclass;
      int* colidx;
      int pclassstart;
      int pc;
      SCIP_Bool* varineq;

      assert(SCIPgetNVars(scip) == matrix->ncols);

      npossiblefixings = 0;
      ndomrelations = 0;
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

      SCIP_CALL( SCIPallocBufferArray(scip, &pclass, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colidx, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varineq, nvars) );
      for( v = 0; v < nvars; v++ )
      {
         colidx[v] = v;
         varineq[v] = FALSE;
      }

      /* before doing dominated column presolving we stuff singleton coninuous columns.
       * this sometimes helps to do a more effective predictive bound analysis.
       */
      SCIP_CALL( singletonColumnStuffing(scip, matrix, varsprocessed, varstofix, &npossiblefixings) );

      /* 1.stage: we search for dominance relations only within parallel columns
       *          concerning equalities and ranged rows
       */

      SCIP_CALL( detectParallelCols(scip, matrix, pclass, varineq) );
      SCIPsortIntInt(pclass, colidx, nvars);

      varcount = 0;

      pc = 0;
      while( pc < nvars )
      {
         int varidx;

         varidx = 0;
         nconfill = 0;
         nintfill = 0;
         nbinfill = 0;

         pclassstart = pclass[pc];
         while( pc < nvars && pclassstart == pclass[pc] )
         {
            varidx = colidx[pc];

            /* we observe only variables which are not processed and present within
             * equalities or ranged rows
             */
            if( !varsprocessed[varidx] && varineq[varidx] )
            {
               /* we search only for dominance relations between the same variable type */
               if( SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_CONTINUOUS )
               {
                  consearchcols[nconfill++] = varidx;
               }
               else if( SCIPvarIsBinary(matrix->vars[varidx]) )
               {
                  binsearchcols[nbinfill++] = varidx;
               }
               else
               {
                  assert(SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_IMPLINT);
                  intsearchcols[nintfill++] = varidx;
               }
            }
            ++pc;
         }

         /* search for dominance relations between continuous variables */
         if( nconfill > 1 )
         {
            SCIP_CALL( findDominancePairs(scip, matrix, presoldata, consearchcols, nconfill, FALSE,
                  varstofix, &npossiblefixings, &ndomrelations) );

            for( v = 0; v < nconfill; ++v )
               varsprocessed[consearchcols[v]] = TRUE;

            varcount += nconfill;
         }
         else if( nconfill == 1 )
         {
            if( varineq[varidx] )
               varsprocessed[consearchcols[0]] = TRUE;
         }

         /* search for dominance relations between integer and impl-integer variables */
         if( nintfill > 1 )
         {
            SCIP_CALL( findDominancePairs(scip, matrix, presoldata, intsearchcols, nintfill, FALSE,
                  varstofix, &npossiblefixings, &ndomrelations) );

            for( v = 0; v < nintfill; ++v )
               varsprocessed[intsearchcols[v]] = TRUE;

            varcount += nintfill;
         }
         else if( nintfill == 1 )
         {
            if( varineq[varidx] )
               varsprocessed[intsearchcols[0]] = TRUE;
         }

         /* search for dominance relations between binary variables */
         if( nbinfill > 1 )
         {
            SCIP_CALL( findDominancePairs(scip, matrix, presoldata, binsearchcols, nbinfill, TRUE,
                  varstofix, &npossiblefixings, &ndomrelations) );

            for( v = 0; v < nbinfill; ++v )
               varsprocessed[binsearchcols[v]] = TRUE;

            varcount += nbinfill;
         }
         else if( nbinfill == 1 )
         {
            if( varineq[varidx] )
               varsprocessed[binsearchcols[0]] = TRUE;
         }

         /* break if no vars are left */
         if( varcount >= nvars )
            break;
      }

      /* 2.stage: we search for dominance relations of the left columns
       *          by row-sparsity order
       */

      /* sort rows per sparsity monotonically increasing */
      SCIPsortIntInt(rowsparsity, rowidxsorted, nrows);

      for( r = 0; r < nrows; ++r )
      {
         int rowidx;
         int* rowpnt;
         int* rowend;

         /* break if the time limit was reached; since the check is expensive,
          * we only check all 1000 constraints
          */
         if( (r % 1000 == 0) && SCIPisStopped(scip) )
            break;

         rowidx = rowidxsorted[r];
         rowpnt = matrix->rowmatind + matrix->rowmatbeg[rowidx];
         rowend = rowpnt + matrix->rowmatcnt[rowidx];

         if( matrix->rowmatcnt[rowidx] == 1 )
            continue;

         nconfill = 0;
         nintfill = 0;
         nbinfill = 0;

         for( ; rowpnt < rowend; rowpnt++ )
         {
            if( !(varsprocessed[*rowpnt]) )
            {
               int varidx;
               varidx = *rowpnt;

               /* we search only for dominance relations between the same variable type */
               if( SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_CONTINUOUS )
               {
                  consearchcols[nconfill++] = varidx;
               }
               else if( SCIPvarIsBinary(matrix->vars[varidx]) )
               {
                  binsearchcols[nbinfill++] = varidx;
               }
               else
               {
                  assert(SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(matrix->vars[varidx]) == SCIP_VARTYPE_IMPLINT);
                  intsearchcols[nintfill++] = varidx;
               }
            }
         }

         /* search for dominance relations between continuous variables */
         if( nconfill > 1 )
         {
            SCIP_CALL( findDominancePairs(scip, matrix, presoldata, consearchcols, nconfill, FALSE,
                  varstofix, &npossiblefixings, &ndomrelations) );

            for( v = 0; v < nconfill; ++v )
               varsprocessed[consearchcols[v]] = TRUE;

            varcount += nconfill;
         }

         /* search for dominance relations between integer and impl-integer variables */
         if( nintfill > 1 )
         {
            SCIP_CALL( findDominancePairs(scip, matrix, presoldata, intsearchcols, nintfill, FALSE,
                  varstofix, &npossiblefixings, &ndomrelations) );

            for( v = 0; v < nintfill; ++v )
               varsprocessed[intsearchcols[v]] = TRUE;

            varcount += nintfill;
         }

         /* search for dominance relations between binary variables */
         if( nbinfill > 1 )
         {
            SCIP_CALL( findDominancePairs(scip, matrix, presoldata, binsearchcols, nbinfill, TRUE,
                  varstofix, &npossiblefixings, &ndomrelations) );

            for( v = 0; v < nbinfill; ++v )
               varsprocessed[binsearchcols[v]] = TRUE;

            varcount += nbinfill;
         }

         /* break if no vars are left */
         if( varcount >= nvars )
            break;
      }

      if( npossiblefixings > 0 )
      {
         int oldnfixedvars = *nfixedvars;

         /* look for fixable variables */
         for( v = matrix->ncols - 1; v >= 0; --v )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;
            SCIP_VAR* var;

            if( varstofix[v] == FIXATLB )
            {
               SCIP_Real lb;

               var = matrix->vars[v];
               lb = SCIPvarGetLbGlobal(var);

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

                  break;
               }
               assert(fixed);
               (*nfixedvars)++;

#ifdef SCIP_DEBUG
               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               {
                  nconvarsfixed++;
               }
               else if( SCIPvarIsBinary(var) )
               {
                  nbinvarsfixed++;
               }
               else
               {
                  assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
                  nintvarsfixed++;
               }
#endif
            }
            else if( varstofix[v] == FIXATUB )
            {
               SCIP_Real ub;

               var = matrix->vars[v];
               ub = SCIPvarGetUbGlobal(var);

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

                  break;
               }
               assert(fixed);

               (*nfixedvars)++;

#ifdef SCIP_DEBUG
               if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
               {
                  nconvarsfixed++;
               }
               else if( SCIPvarIsBinary(var) )
               {
                  nbinvarsfixed++;
               }
               else
               {
                  assert(SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
                  nintvarsfixed++;
               }
#endif
            }
         }

         if( *result != SCIP_CUTOFF && *nfixedvars > oldnfixedvars )
            *result = SCIP_SUCCESS;
      }

      SCIPfreeBufferArray(scip, &varineq);
      SCIPfreeBufferArray(scip, &colidx);
      SCIPfreeBufferArray(scip, &pclass);
      SCIPfreeBufferArray(scip, &rowsparsity);
      SCIPfreeBufferArray(scip, &rowidxsorted);
      SCIPfreeBufferArray(scip, &binsearchcols);
      SCIPfreeBufferArray(scip, &intsearchcols);
      SCIPfreeBufferArray(scip, &consearchcols);
      SCIPfreeBufferArray(scip, &varsprocessed);
      SCIPfreeBufferArray(scip, &varstofix);

#ifdef SCIP_DEBUG
      if( (nconvarsfixed + nintvarsfixed + nbinvarsfixed) > 0 )
      {
         SCIPdebugMessage("### %d vars [%lld dom] ===>>> fixed [cont: %d, int: %d, bin: %d], %scutoff detected\n",
            matrix->ncols, ndomrelations, nconvarsfixed, nintvarsfixed, nbinvarsfixed, (*result != SCIP_CUTOFF) ? "no " : "");
      }
#endif
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
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create domcol presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );

   /* include presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_DELAY, presolExecDomcol, presoldata) );
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyDomcol) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeDomcol) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/domcol/maxpairs",
         "maximal number of pair comparisons for at least one variable fixing",
         &presoldata->maxpairs, FALSE, DEFAULT_MAXPAIRS, 10, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
