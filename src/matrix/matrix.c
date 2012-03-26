/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    matrix.c
 * @ingroup
 * @brief   constraint matrix data structure
 * @author  Dieter Weninger
 * @author  Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <assert.h>
#include <string.h>

#include "scip/cons_knapsack.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"

#include "matrix.h"

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
      assert( requiredsize <= *nvars );
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
   SCIP_Real             rhs,                /**< right hand side */
   int                   considx             /**< row index */
   )
{
   int j;
   int probindex;

   assert(vars != NULL);
   assert(vals != NULL);

   matrix->lhs[considx] = lhs;
   matrix->rhs[considx] = rhs;

   matrix->rowmatbeg[considx] = matrix->nnonzs;

   for( j = 0; j < nvars; j++ )
   {
      matrix->rowmatval[matrix->nnonzs] = vals[j];
      probindex = SCIPvarGetProbindex(vars[j]);
      assert(matrix->vars[probindex] == vars[j]);

      assert(0 <= probindex && probindex < matrix->ncols);
      matrix->rowmatind[matrix->nnonzs] = probindex;
      matrix->nnonzs = matrix->nnonzs + 1;
   }

   matrix->rowmatcnt[considx] = matrix->nnonzs - matrix->rowmatbeg[considx];

   return SCIP_OKAY;
}

/** add empty row to matrix */
static
SCIP_RETCODE addEmptyRow(
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix */
   int                   considx             /**< row/constraint index */
   )
{
   matrix->lhs[considx] = 0;
   matrix->rhs[considx] = 0;

   matrix->rowmatbeg[considx] = matrix->nnonzs;
   matrix->rowmatcnt[considx] = matrix->nnonzs - matrix->rowmatbeg[considx];

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
   SCIP_Real             rhs,                /**< right hand side */
   int                   considx             /**< constraint index */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(vars != NULL);
   assert(lhs <= rhs);

   if( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) )
      return SCIP_OKAY; /* constraint is redundant */

   activevars = NULL;
   activevals = NULL;
   nactivevars = nvars;
   activeconstant = 0;

   if( nvars > 0 )
   {
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
   }

   /* adapt left and right hand side */
   if( !SCIPisInfinity(scip, -lhs) )
   {
      lhs -= activeconstant;
   }
   if( !SCIPisInfinity(scip, rhs) )
   {
      rhs -= activeconstant;
   }
   /* add single row to matrix */
   if( nactivevars > 0 )
   {
      SCIP_CALL( addRow(matrix, activevars, activevals, nactivevars, lhs, rhs, considx) );
   }
   else
   {
      SCIP_CALL( addEmptyRow(matrix, considx) );
   }

   if( nvars > 0 )
   {
      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &activevals);
      SCIPfreeBufferArray(scip, &activevars);
   }

   return SCIP_OKAY;
}

/** get number of active variables from given variables */
static
SCIP_RETCODE getNumberActiveVars(
   SCIP*                 scip,               /**< current scip instance */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_Real*            vals,               /**< variable coefficients */
   int                   nvars,              /**< number of variables */
   int*                  nactvars            /**< number of active variables */
   )
{
   int v;
   SCIP_VAR** activevars;
   SCIP_Real* activevals;
   int nactivevars;
   SCIP_Real activeconstant;

   assert(scip != NULL);
   assert(vars != NULL);
   assert(nactvars != NULL);

   v = 0;
   activevars = NULL;
   activevals = NULL;
   nactivevars = nvars;
   activeconstant = 0;

   if( nvars > 0 )
   {
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
   }

   if( nvars > 0 )
   {
      /* free buffer arrays */
      SCIPfreeBufferArray(scip, &activevals);
      SCIPfreeBufferArray(scip, &activevars);
   }

   *nactvars = nactivevars;

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

   for( i = 0; i < matrix->ncols; i++ )
   {
      fillidx[i] = 0;
      matrix->colmatcnt[i] = 0;
   }

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
         matrix->colmatval[matrix->colmatbeg[colidx]+fillidx[colidx]] = *valpnt;
         matrix->colmatind[matrix->colmatbeg[colidx]+fillidx[colidx]] = i;
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
SCIP_RETCODE initMatrix(
   SCIP*                 scip,               /**< current scip instance */
   CONSTRAINTMATRIX*     matrix,             /**< constraint matrix object to be initialized */
   SCIP_Bool*            initialized         /**< was the initialization successful? */
   )
{
   int v;
   int c;
   int i;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   SCIP_CONSHDLR* conshdlr;
   const char* conshdlrname;
   int nnonzstmp;
   int conscounter;
   SCIP_CONSHDLR** conshdlrs;
   int nconshdlrs;
   int nconss;

   nnonzstmp = 0;
   conscounter = 0;

   assert(scip != NULL);
   assert(matrix != NULL);
   assert(initialized != NULL);

   /* set everything to zero */
   matrix->colmatval = NULL;
   matrix->colmatind = NULL;
   matrix->colmatbeg = NULL;
   matrix->colmatcnt = NULL;
   matrix->ncols = 0;
   matrix->lb = NULL;
   matrix->ub = NULL;
   matrix->vars = NULL;
   matrix->rowmatval = NULL;
   matrix->rowmatind = NULL;
   matrix->rowmatbeg = NULL;
   matrix->rowmatcnt = NULL;
   matrix->nrows = 0;
   matrix->lhs = NULL;
   matrix->rhs = NULL;
   matrix->conss = NULL;
   matrix->minactivity = NULL;
   matrix->maxactivity = NULL;
   matrix->minactivityneginf = NULL;
   matrix->minactivityposinf = NULL;
   matrix->maxactivityneginf = NULL;
   matrix->maxactivityposinf = NULL;
   matrix->nnonzs = 0;

   /* return if no variables or constraints are present */
   if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
      return SCIP_OKAY;

   /* loop over all constraint handlers and collect
      the number of enforced constraints */
   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);
   nconss = 0;
   for( i = 0; i < nconshdlrs; ++i )
   {
      nconss += SCIPconshdlrGetNEnfoConss(conshdlrs[i]);
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->conss, nconss) );

   /* copy the constraints */
   nconss = 0;
   for( i = 0; i < nconshdlrs; ++i )
   {
      SCIP_CONS** conshdlrconss;
      int nconshdlrconss;

      conshdlrconss = SCIPconshdlrGetEnfoConss(conshdlrs[i]);
      nconshdlrconss = SCIPconshdlrGetNEnfoConss(conshdlrs[i]);

      for( c = 0; c < nconshdlrconss; ++c )
      {
         matrix->conss[nconss] = conshdlrconss[c];
         nconss++;
      }
   }
   matrix->nrows = nconss;
   matrix->ncols = SCIPgetNVars(scip);

   SCIP_CALL( SCIPduplicateBufferArray(scip, &matrix->vars, SCIPgetVars(scip), matrix->ncols) );

   /* count the number of nonzeros exact */
   for( c = 0; c < matrix->nrows; c++ )
   {
      SCIP_VAR** cvars;
      int ncvars;
      SCIP_Real* cvals;
      int nactvars;

      cons = matrix->conss[c];
      assert( cons != NULL );

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( TRUE == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         cvars = SCIPgetVarsLinear(scip, cons);
         ncvars = SCIPgetNVarsLinear(scip ,cons);
         cvals = SCIPgetValsLinear(scip, cons);
         nactvars = 0;
         SCIP_CALL( getNumberActiveVars(scip,cvars,cvals,ncvars,&nactvars) );
         assert(nactvars >= 0);
         nnonzstmp += nactvars;
         conscounter++;
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         cvars = SCIPgetVarsSetppc(scip, cons);
         ncvars = SCIPgetNVarsSetppc(scip, cons);
         cvals = NULL;
         nactvars = 0;
         SCIP_CALL( getNumberActiveVars(scip,cvars,cvals,ncvars,&nactvars) );
         assert(nactvars >= 0);
         nnonzstmp += nactvars;
         conscounter++;
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         cvars = SCIPgetVarsLogicor(scip, cons);
         ncvars = SCIPgetNVarsLogicor(scip, cons);
         cvals = NULL;
         nactvars = 0;
         SCIP_CALL( getNumberActiveVars(scip,cvars,cvals,ncvars,&nactvars) );
         assert(nactvars >= 0);
         nnonzstmp += nactvars;
         conscounter++;
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;
         SCIP_Real* consvals;

         weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, SCIPgetNVarsKnapsack(scip, cons)) );
         for( v = 0; v < SCIPgetNVarsKnapsack(scip, cons); v++ )
            consvals[v] = weights[v];

         cvars = SCIPgetVarsKnapsack(scip, cons);
         ncvars = SCIPgetNVarsKnapsack(scip, cons);

         cvals = NULL;
         nactvars = 0;
         SCIP_CALL( getNumberActiveVars(scip,cvars,consvals,ncvars,&nactvars) );
         assert(nactvars >= 0);
         nnonzstmp += nactvars;
         conscounter++;

         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         SCIP_VAR** consvars;
         SCIP_Real* consvals;

         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         ncvars = 2;
         nactvars = 0;
         SCIP_CALL( getNumberActiveVars(scip,consvars,consvals,ncvars,&nactvars) );
         assert(nactvars >= 0);
         nnonzstmp += nactvars;

         SCIPfreeBufferArray(scip, &consvals);
         SCIPfreeBufferArray(scip, &consvars);
         conscounter++;
      }
   }

   /* do nothing if we have unsupported constraint types or no entries */
   if( matrix->nrows > conscounter || nnonzstmp == 0 || matrix->nrows == 0 || matrix->ncols == 0 )
   {
      SCIPfreeBufferArray(scip, &matrix->conss);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatbeg, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatcnt, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lb, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->ub, matrix->ncols) );

   for( v = 0; v < matrix->ncols; v++ )
   {
      var = matrix->vars[v];
      assert( var != NULL );

      matrix->lb[v] = SCIPvarGetLbGlobal(var);
      matrix->ub[v] = SCIPvarGetUbGlobal(var);
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatbeg, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatcnt, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lhs, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rhs, matrix->nrows) );

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivity, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivity, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivityneginf, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->minactivityposinf, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivityneginf, matrix->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->maxactivityposinf, matrix->nrows) );

   /* loop a second time over constraints and add them to the matrix */
   for( c = 0; c < matrix->nrows; c++ )
   {
      cons = matrix->conss[c];
      assert( cons != NULL );

      conshdlr = SCIPconsGetHdlr(cons);
      assert( conshdlr != NULL );

      conshdlrname = SCIPconshdlrGetName(conshdlr);
      assert( TRUE == SCIPconsIsTransformed(cons) );

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         SCIP_CALL( addConstraint(scip,
               matrix,
               SCIPgetVarsLinear(scip, cons),
               SCIPgetValsLinear(scip, cons),
               SCIPgetNVarsLinear(scip, cons),
               SCIPgetLhsLinear(scip, cons),
               SCIPgetRhsLinear(scip, cons),
               c) );
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         SCIP_Real lhs;
         SCIP_Real rhs;

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

         SCIP_CALL( addConstraint(scip,
               matrix,
               SCIPgetVarsSetppc(scip, cons),
               NULL,
               SCIPgetNVarsSetppc(scip, cons),
               lhs,
               rhs,
               c) );
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         SCIP_CALL( addConstraint(scip,
               matrix,
               SCIPgetVarsLogicor(scip, cons),
               NULL,
               SCIPgetNVarsLogicor(scip, cons),
               1.0,
               SCIPinfinity(scip),
               c) );
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         SCIP_Longint* weights;
         SCIP_Real* consvals;

         weights = SCIPgetWeightsKnapsack(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, SCIPgetNVarsKnapsack(scip, cons)) );
         for( v = 0; v < SCIPgetNVarsKnapsack(scip, cons); v++ )
            consvals[v] = weights[v];

         SCIP_CALL( addConstraint(scip,
               matrix,
               SCIPgetVarsKnapsack(scip, cons),
               consvals,
               SCIPgetNVarsKnapsack(scip, cons),
               -SCIPinfinity(scip),
               (SCIP_Real)SCIPgetCapacityKnapsack(scip, cons),
               c) );

         SCIPfreeBufferArray(scip, &consvals);
      }
      else if( strcmp(conshdlrname, "varbound") == 0 )
      {
         SCIP_VAR** consvars;
         SCIP_Real* consvals;

         SCIP_CALL( SCIPallocBufferArray(scip, &consvars, 2) );
         SCIP_CALL( SCIPallocBufferArray(scip, &consvals, 2) );

         consvars[0] = SCIPgetVarVarbound(scip, cons);
         consvars[1] = SCIPgetVbdvarVarbound(scip, cons);

         consvals[0] = 1.0;
         consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

         SCIP_CALL( addConstraint(scip,
               matrix,
               consvars,
               consvals,
               2,
               SCIPgetLhsVarbound(scip, cons),
               SCIPgetRhsVarbound(scip, cons),
               c) );

         SCIPfreeBufferArray(scip, &consvals);
         SCIPfreeBufferArray(scip, &consvars);
      }
   }

   assert(nnonzstmp == matrix->nnonzs);

   /* calculate row activity bounds */
   SCIP_CALL( calcActivityBounds(scip,matrix) );

   /* transform row major format into column major format */
   SCIP_CALL( setColumnMajorFormat(scip,matrix) );

   *initialized = TRUE;

   return SCIP_OKAY;
}


/** frees the constraint matrix */
void freeMatrix(
   SCIP*                 scip,               /**< current SCIP instance */
   CONSTRAINTMATRIX**    matrix              /**< constraint matrix object */
   )
{
   assert(scip != NULL);
   assert(matrix != NULL);

   if( (*matrix)->nnonzs > 0 )
   {
      assert((*matrix) != NULL);

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

      SCIPfreeBufferArray(scip, &((*matrix)->conss));

      (*matrix)->nrows = 0;
      (*matrix)->ncols = 0;
      (*matrix)->nnonzs = 0;
   }

   SCIPfreeBufferArrayNull(scip, &((*matrix)->vars));

   SCIPfreeBuffer(scip, matrix);
}
