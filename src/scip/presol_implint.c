/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_implint.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief  Presolver that detects implicit integer variables
 * @author Rolf van der Hulst
 */

/* TODO: explore integer to implicit integer conversion */
/* TODO: fix to work in MINLP context */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/cons_and.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_or.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_xor.h"

#include "scip/presol_implint.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_network.h"
#include "scip/pub_presol.h"
#include "scip/pub_var.h"

#include "scip/scip_cons.h"
#include "scip/scip_copy.h"
#include "scip/scip_general.h"
#include "scip/scip_message.h"
#include "scip/scip_mem.h"
#include "scip/scip_nlp.h"
#include "scip/scip_numerics.h"
#include "scip/scip_param.h"
#include "scip/scip_presol.h"
#include "scip/scip_pricer.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_timing.h"
#include "scip/scip_var.h"

#define PRESOL_NAME             "implint"
#define PRESOL_DESC             "detects implicit integer variables"
/* We want to run as late as possible, but before symmetry detection.
 * The main reason for this is that symmetry detection may add linear constraints that
 * impede the detection of implied integrality, but do not break implied integrality itself.
 * Also, symmetry methods rely on the fact that each variable in an orbit is integral,
 * as otherwise certain reductions may break. So it is currently not safe to run implied integrality detection
 * after symmetry methods are applied. */
#define PRESOL_PRIORITY         (-900000) /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS        0 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_COLUMNROWRATIO  50.0
#define DEFAULT_NUMERICSLIMIT   1e6
#define DEFAULT_CONVERTINTEGERS FALSE

/** presolver data */
struct SCIP_PresolData
{
   SCIP_Bool             computedimplints;   /**< Were implied integers already computed? */

   SCIP_Real             columnrowratio;     /**< Use the network row addition algorithm when the column to row ratio
                                               * becomes larger than this threshold. Otherwise, use column addition. */
   SCIP_Real             numericslimit;      /**< A row that contains variables with coefficients that are greater in
                                                * absolute value than this limit is not considered for
                                                * implied integrality detection. */
   SCIP_Bool             convertintegers;    /**< Controls whether implied integrality is inferred for integer variables */
};

/** constraint matrix data structure in column and row major format.
 *  Contains only the linear terms, and marks the presence of non-linear terms */
struct ImplintMatrix
{
   SCIP_Real*            colmatval;          /**< coefficients in column major format */
   int*                  colmatind;          /**< row indexes in column major format */
   int*                  colmatbeg;          /**< column storage offset */
   int*                  colmatcnt;          /**< number of row entries per column */
   int                   ncols;              /**< complete number of columns */
   SCIP_Real*            lb;                 /**< lower bound per variable */
   SCIP_Real*            ub;                 /**< upper bound per variable */
   SCIP_Bool*            colintegral;
   // implied integral? bounds integral? number of +-1 nonzeros?
   // contained in nonlinear term? ntimes operand / resultant in logical constraints?
   // nconstraints (different from nnonz because of multiple row constraints)
   // npmonenonzeros in integral equality rows
   SCIP_VAR**            vars;

   SCIP_Real*            rowmatval;          /**< coefficients in row major format */
   int*                  rowmatind;          /**< column indexed in row major format */
   int*                  rowmatbeg;          /**< row storage offset */
   int*                  rowmatcnt;          /**< number of column entries per row */
   // boundsintegral, rowequality, rowbadnumerics, rowncontinuous, rowncontinuouspmone, rownintegralpmone

   int                   nrows;              /**< complete number of rows */
   SCIP_Real*            lhs;                /**< left hand side per row */
   SCIP_Real*            rhs;                /**< right hand side per row */

   SCIP_CONS**           rowcons;            /**< Constraint where the row originated from */

   int                   nnonzs;             /**< sparsity counter */
   int                   memnonz;
};
typedef struct ImplintMatrix IMPLINT_MATRIX;

static
SCIP_Real* matrixGetColumnValues(
   IMPLINT_MATRIX* matrix,
   int column
   )
{
   assert(matrix);
   return matrix->colmatval + matrix->colmatbeg[column];
}

static
int* matrixGetColumnIndices(
   IMPLINT_MATRIX* matrix,
   int column
   )
{
   assert(matrix);
   return matrix->colmatind + matrix->colmatbeg[column];
}

static
int matrixGetColumnNNonz(
   IMPLINT_MATRIX* matrix,
   int column
   )
{
   assert(matrix);
   return matrix->colmatcnt[column];
}

static
SCIP_Real* matrixGetRowValues(
   IMPLINT_MATRIX* matrix,
   int row
)
{
   assert(matrix);
   return matrix->rowmatval + matrix->rowmatbeg[row];
}

static
int* matrixGetRowIndices(
   IMPLINT_MATRIX* matrix,
   int row
)
{
   assert(matrix);
   return matrix->rowmatind + matrix->rowmatbeg[row];
}

static
int matrixGetRowNNonz(
   IMPLINT_MATRIX* matrix,
   int row
)
{
   assert(matrix);
   return matrix->rowmatcnt[row];
}

static
int matrixGetNRows(
   IMPLINT_MATRIX* matrix
   )
{
   assert(matrix);
   return matrix->nrows;
}

static
int matrixGetNCols(
   IMPLINT_MATRIX* matrix
   )
{
   assert(matrix);
   return matrix->ncols;
}

static
SCIP_VAR* matrixGetVar(
   IMPLINT_MATRIX* matrix,
   int column
   )
{
   assert(matrix);
   assert(column >= 0 && column < matrix->ncols);
   return matrix->vars[column];
}

static
SCIP_Bool matrixColIsIntegral(
   IMPLINT_MATRIX* matrix,
   int column
   )
{
   assert(matrix);
   assert(column >= 0 && column < matrix->ncols);
   return matrix->colintegral[column];
}

static
SCIP_Real matrixGetColLb(
   IMPLINT_MATRIX* matrix,
   int column
   )
{
   assert(matrix);
   assert(column >= 0 && column < matrix->ncols);
   return matrix->lb[column];
}

static
SCIP_Real matrixGetColUb(
   IMPLINT_MATRIX* matrix,
   int column
   )
{
   assert(matrix);
   assert(column >= 0 && column < matrix->ncols);
   return matrix->ub[column];
}

static
SCIP_Real matrixGetRowLhs(
   IMPLINT_MATRIX* matrix,
   int row
   )
{
   assert(matrix);
   assert(row >= 0 && row < matrix->nrows);
   return matrix->lhs[row];
}

static
SCIP_Real matrixGetRowRhs(
   IMPLINT_MATRIX* matrix,
   int row
   )
{
   assert(matrix);
   assert(row >= 0 && row < matrix->nrows);
   return matrix->rhs[row];
}

/** transforms given variables, scalars and constant to the corresponding active variables, scalars and constant */
static
SCIP_RETCODE getActiveVariables(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_VAR***           vars,               /**< vars array to get active variables for */
   SCIP_Real**           scalars,            /**< scalars a_1, ..., a_n in linear sum a_1*x_1 + ... + a_n*x_n + c */
   int*                  nvars,              /**< pointer to number of variables and values in vars and vals array */
   SCIP_Real*            constant            /**< pointer to constant c in linear sum a_1*x_1 + ... + a_n*x_n + c */
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

   SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, *nvars, constant, &requiredsize) );

   if( requiredsize > *nvars )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, vars, requiredsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, scalars, requiredsize) );

      /* call function a second time with enough memory */
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, *vars, *scalars, nvars, requiredsize, constant, &requiredsize) );
   }
   assert(requiredsize == *nvars);

   return SCIP_OKAY;
}

/** add one row to the constraint matrix */
static
SCIP_RETCODE matrixAddRow(
   SCIP*                 scip,               /**< SCIP data structure */
   IMPLINT_MATRIX*       matrix,             /**< constraint matrix */
   SCIP_VAR**            vars,               /**< variables of this row */
   SCIP_Real*            vals,               /**< coefficients of this row */
   int                   nvars,              /**< number of variables of this row */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_CONS*            cons                /**< constraint where the row originated from */
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
      /* ignore variables with very small coefficients */
      if( SCIPisZero(scip, vals[j]) )
         continue;

      assert(matrix->nnonzs < matrix->memnonz);
      matrix->rowmatval[matrix->nnonzs] = vals[j];
      probindex = SCIPvarGetProbindex(vars[j]);
      assert(0 <= probindex && probindex < matrix->ncols);
      matrix->rowmatind[matrix->nnonzs] = probindex;
      (matrix->nnonzs)++;
   }

   matrix->rowmatcnt[rowidx] = matrix->nnonzs - matrix->rowmatbeg[rowidx];
   matrix->rowcons[rowidx] = cons;

   ++(matrix->nrows);

   return SCIP_OKAY;
}

static
SCIP_RETCODE addLinearConstraint(
   SCIP*                 scip,               /**< current scip instance */
   IMPLINT_MATRIX*       matrix,             /**< constraint matrix */
   SCIP_VAR**            vars,               /**< variables of this constraint */
   SCIP_Real*            vals,               /**< variable coefficients of this constraint.
                                              **< If set to NULL, all values are assumed to be equal to 1.0. */
   int                   nvars,              /**< number of variables */
   SCIP_Real             lhs,                /**< left hand side */
   SCIP_Real             rhs,                /**< right hand side */
   SCIP_CONS*            cons                /**< constraint belonging to the row */
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
      SCIP_CALL( matrixAddRow(scip, matrix, activevars, activevals, nactivevars, lhs, rhs, cons) );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevals);
   SCIPfreeBufferArray(scip, &activevars);

   return SCIP_OKAY;
}

static
SCIP_RETCODE addAndOrLinearization(
   SCIP*                 scip,               /**< current scip instance */
   IMPLINT_MATRIX*       matrix,             /**< constraint matrix */
   SCIP_CONS*            cons,               /**< The constraint that is linearized */
   SCIP_VAR**            operands,           /**< variables of this constraint */
   int                   noperands,          /**< number of operands */
   SCIP_VAR*             resultant,          /**< Resultant variable */
   SCIP_Bool             isAndCons           /**< Indicates if the constraint is an AND or OR linearization */
)
{
   SCIP_Real* vals;
   SCIP_VAR** vars;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int i;
   assert(noperands > 0);
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, noperands + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, noperands + 1) );
   /* First, add all the constraints of the form resultant <= operand */

   if( isAndCons )
   {
      lhs = -SCIPinfinity(scip);
      rhs = 0.0;
   }
   else
   {
      lhs = 0.0;
      rhs = SCIPinfinity(scip);
   }
   vals[0] = 1.0;
   vals[1] = -1.0;
   vars[0] = resultant;
   for( i = 0; i < noperands ; ++i )
   {
      vars[1] = operands[i];

      SCIP_CALL( addLinearConstraint(scip, matrix, vars, vals, 2, lhs, rhs, cons) );
   }
   for( i = 0; i < noperands ; ++i )
   {
      vars[i+1] = operands[i];
      vals[i+1] = -1.0;
   }

   if( isAndCons )
   {
      lhs = 1 - noperands;
      rhs = SCIPinfinity(scip);
   }
   else
   {
      lhs = -SCIPinfinity(scip);
      rhs = 0.0;
   }
   SCIP_CALL( addLinearConstraint(scip, matrix, vars, vals, noperands + 1, lhs, rhs, cons) );
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &vals);
   return SCIP_OKAY;
}

static
SCIP_RETCODE addXorLinearization(
   SCIP*                 scip,               /**< current scip instance */
   IMPLINT_MATRIX*       matrix,             /**< constraint matrix */
   SCIP_CONS*            cons                /**< The constraint that is linearized */
)
{
   SCIP_Real* vals;
   SCIP_VAR** vars;
   SCIP_VAR** operands = SCIPgetVarsXor(scip, cons);
   int noperands = SCIPgetNVarsXor(scip, cons);
   assert(noperands > 0);
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, noperands + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, noperands + 1) );
   if( noperands == 3)
   {
      assert(SCIPgetIntVarXor(scip, cons) == NULL);
      /** in the special case of 3 variables and c = 0, the following linear system is created:
       *    + x - y - z <= 0
       *    - x + y - z <= 0
       *    - x - y + z <= 0
       *    + x + y + z <= 2
       *  in the special case of 3 variables and c = 1, the following linear system is created:
       *    - x + y + z <= 1
       *    + x - y + z <= 1
       *    + x + y - z <= 1
       *    - x - y - z <= -1
       */
      SCIP_Bool rhs = SCIPgetRhsXor(scip, cons);
      SCIP_Real scale = rhs == 0 ? -1.0 : 1.0;
      SCIP_Real rhsVal = rhs;

      for( int i = 0; i < 3; ++i )
      {
         for( int j = 0; j < 3; ++j )
         {
            vals[j] = (i == j ) ? -scale : scale;
         }

         SCIP_CALL( addLinearConstraint(scip, matrix, operands, vals, 3, -SCIPinfinity(scip), rhsVal, cons) );
      }
      for( int j = 0; j < 3; ++j )
      {
         vals[j] = -scale;
      }

      rhsVal = 2.0 - 3 * rhs;
      SCIP_CALL( addLinearConstraint(scip, matrix, operands, vals, 3, -SCIPinfinity(scip), rhsVal, cons) );
   }
   else
   {
      /* First, add all the constraints of the form resultant <= operand */
      for( int i = 0; i < noperands; ++i )
      {
         vars[i] = operands[i];
         vals[i] = 1.0;
      }
      vars[noperands] = SCIPgetIntVarXor(scip, cons);
      vals[noperands] = -2.0;
      assert(vars[noperands] != NULL);

      SCIP_Bool rhs = SCIPgetRhsXor(scip, cons);

      SCIP_CALL( addLinearConstraint(scip, matrix, vars, vals, noperands + 1, rhs, rhs, cons) );
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &vals);
   return SCIP_OKAY;
}

/** transform row major format into column major format */
static
SCIP_RETCODE matrixSetColumnMajor(
   SCIP*                 scip,               /**< current scip instance */
   IMPLINT_MATRIX*       matrix              /**< constraint matrix */
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

/* @todo; use symmetry constraints to guide variable ordering for integral columns.
 * Symmetric variables should always all be either network or non-network */
/* @todo optionally skip construction of integral constraints if we do not run detection on integer variables */
static
SCIP_RETCODE matrixCreate(
   SCIP* scip,
   IMPLINT_MATRIX** pmatrix,
   SCIP_Bool* success
   )
{
   int nconshdlrs;
   SCIP_CONSHDLR** conshdlrs;
   int nmatrixrows;
   int i;
   int nconshdlrconss;
   const char* conshdlrname;
   int j;
   IMPLINT_MATRIX* matrix;
   SCIP_VAR** vars;
   int nvars;
   int nnonzstmp;

   *success = FALSE;

   /* return if no variables or constraints are present */
   if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
      return SCIP_OKAY;

   assert(SCIPgetNActivePricers(scip) == 0);

   /* loop over all constraint handlers and collect the number of checked constraints that contribute rows
    * to the matrix. */
   nconshdlrs = SCIPgetNConshdlrs(scip);
   conshdlrs = SCIPgetConshdlrs(scip);
   nmatrixrows = 0;
   nnonzstmp = 0;

   for( i = 0; i < nconshdlrs; ++i )
   {
      nconshdlrconss = SCIPconshdlrGetNCheckConss(conshdlrs[i]);

      if( nconshdlrconss > 0 )
      {
         conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);

         /* constraint handlers which can always be represented by a single row */
         if( strcmp(conshdlrname, "linear") == 0 || strcmp(conshdlrname, "knapsack") == 0
               || strcmp(conshdlrname, "setppc") == 0 || strcmp(conshdlrname, "logicor") == 0
               || strcmp(conshdlrname, "varbound") == 0 )
         {
            nmatrixrows += nconshdlrconss;
         }
         else if( strcmp(conshdlrname, "and") == 0 )
         {
            /* The linearization of AND constraints is modelled using n+1 inequalities on n variables. */
            SCIP_CONS** checked = SCIPconshdlrGetCheckConss(conshdlrs[i]);
            for( j = 0; j < nconshdlrconss; ++j )
            {
               nmatrixrows = nmatrixrows + SCIPgetNVarsAnd(scip, checked[j]) + 1;
               nnonzstmp += SCIPgetNVarsAnd(scip, checked[j]);
            }
         }
         else if( strcmp(conshdlrname, "or") == 0 )
         {
            /* The linearization of OR constraints is modelled using n+1 inequalities on n variables. */
            SCIP_CONS** checked = SCIPconshdlrGetCheckConss(conshdlrs[i]);
            for( j = 0; j < nconshdlrconss; ++j )
            {
               nmatrixrows = nmatrixrows + SCIPgetNVarsOr(scip, checked[j]) + 1;
               nnonzstmp += SCIPgetNVarsOr(scip, checked[j]);
            }
         }
         else if( strcmp(conshdlrname, "xor") == 0 )
         {
            /* The relaxation of XOR constraints is handled differently depending on their size
             * For 3 variables, the integer hull of the constraint is added as it is only 4 rows.
             * For more variables, the LP only has a single row. */
            SCIP_CONS** checked = SCIPconshdlrGetCheckConss(conshdlrs[i]);
            for( j = 0; j < nconshdlrconss; ++j )
            {
               int nxorvars = SCIPgetNVarsXor(scip, checked[j]);
               if( nxorvars == 3 )
               {
                  nmatrixrows += 4;
                  nnonzstmp += 12;
               }
               else
                  nmatrixrows += 1;
            }
         }
         else if( strcmp(conshdlrname, "orbitope_pp") == 0 || strcmp(conshdlrname, "orbitope_full") == 0
               || strcmp(conshdlrname, "orbisack") == 0 || strcmp(conshdlrname, "symresack") == 0 )
         {

         }
         else
         {
            /* TODO: support indicator, nonlinear, superindicator, sos1, sos2, bounddisjunction, linking conshdlrs */
            return SCIP_OKAY;
         }
      }
   }

   if( nmatrixrows == 0 )
      return SCIP_OKAY;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* approximate number of nonzeros by taking for each variable the number of up- and downlocks;
    * this counts nonzeros in equalities twice, but can be at most two times as high as the exact number
    */
   for( i = 0; i < nvars; ++i )
   {
      nnonzstmp += SCIPvarGetNLocksDownType(vars[i], SCIP_LOCKTYPE_MODEL);
      nnonzstmp += SCIPvarGetNLocksUpType(vars[i], SCIP_LOCKTYPE_MODEL);
   }

   if( nnonzstmp == 0 )
      return SCIP_OKAY;

   *success = TRUE;

   /* build the matrix structure */
   SCIP_CALL( SCIPallocBuffer(scip, pmatrix) );
   matrix = *pmatrix;

   SCIP_CALL( SCIPduplicateBufferArray(scip, &matrix->vars, vars, nvars) );
   matrix->ncols = nvars;
   matrix->memnonz = nnonzstmp;
   matrix->nnonzs = 0;
   matrix->nrows = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatbeg, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatcnt, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lb, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->ub, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colintegral, matrix->ncols) );

   /* init bounds */
   for( i = 0; i < matrix->ncols; ++i )
   {
      matrix->lb[i] = SCIPvarGetLbGlobal(vars[i]);
      matrix->ub[i] = SCIPvarGetUbGlobal(vars[i]);
      matrix->colintegral[i] = SCIPvarIsIntegral(vars[i]);
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatbeg, nmatrixrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatcnt, nmatrixrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lhs, nmatrixrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rhs, nmatrixrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowcons, nmatrixrows) );


   /* loop a second time over constraints handlers and add supported constraints to the matrix */
   for( i = 0; i < nconshdlrs && *success; ++i)
   {
      SCIP_CONS** conshdlrconss;
      int c;

      conshdlrname = SCIPconshdlrGetName(conshdlrs[i]);
      conshdlrconss = SCIPconshdlrGetCheckConss(conshdlrs[i]);
      nconshdlrconss = SCIPconshdlrGetNCheckConss(conshdlrs[i]);

      if( strcmp(conshdlrname, "linear") == 0 )
      {
         for( c = 0; c < nconshdlrconss; ++c )
         {
            SCIP_CONS* cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            if( SCIPconsIsModifiable(cons) )
            {
               *success = FALSE;
               break;
            }

            SCIP_CALL( addLinearConstraint(scip, matrix, SCIPgetVarsLinear(scip, cons), SCIPgetValsLinear(scip, cons),
                  SCIPgetNVarsLinear(scip, cons), SCIPgetLhsLinear(scip, cons), SCIPgetRhsLinear(scip, cons), cons) );
         }
      }
      else if( strcmp(conshdlrname, "knapsack") == 0 )
      {
         if( nconshdlrconss > 0 )
         {
            SCIP_Real* consvals;
            int nrowvars;

            SCIP_CALL( SCIPallocBufferArray(scip, &consvals, nvars) );

            for( c = 0; c < nconshdlrconss; ++c )
            {
               SCIP_Longint* weights;
               SCIP_Real rhs;
               SCIP_CONS* cons = conshdlrconss[c];
               assert(SCIPconsIsTransformed(cons));

               if( SCIPconsIsModifiable(cons) )
               {
                  *success = FALSE;
                  break;
               }
               weights = SCIPgetWeightsKnapsack(scip, cons);
               nrowvars = SCIPgetNVarsKnapsack(scip, cons);
               for( int v = 0; v < nrowvars; ++v )
                  consvals[v] = (SCIP_Real)weights[v];

               rhs = (SCIP_Real) SCIPgetCapacityKnapsack(scip, cons);
               SCIP_CALL( addLinearConstraint(scip, matrix, SCIPgetVarsKnapsack(scip, cons), consvals, nrowvars,
                     -SCIPinfinity(scip), rhs, cons) );
            }

            SCIPfreeBufferArray(scip, &consvals);
         }
      }
      else if( strcmp(conshdlrname, "setppc") == 0 )
      {
         for( c = 0; c < nconshdlrconss; ++c )
         {
            SCIP_Real lhs;
            SCIP_Real rhs;

            SCIP_CONS* cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            /* do not include constraints that can be altered due to column generation */
            if( SCIPconsIsModifiable(cons) )
            {
               *success = FALSE;
               break;
            }

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
                  SCIPABORT();
            }

            SCIP_CALL( addLinearConstraint(scip, matrix, SCIPgetVarsSetppc(scip, cons), NULL,
                                           SCIPgetNVarsSetppc(scip, cons), lhs, rhs, cons) );
         }
      }
      else if( strcmp(conshdlrname, "logicor") == 0 )
      {
         for( c = 0; c < nconshdlrconss; ++c )
         {
            SCIP_CONS* cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            if( SCIPconsIsModifiable(cons) )
            {
               *success = FALSE;
               break;
            }

            SCIP_CALL( addLinearConstraint(scip, matrix, SCIPgetVarsLogicor(scip, cons), NULL,
                  SCIPgetNVarsLogicor(scip, cons), 1.0, SCIPinfinity(scip), cons) );
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

            for( c = 0; c < nconshdlrconss; ++c )
            {
               SCIP_CONS* cons = conshdlrconss[c];
               assert(SCIPconsIsTransformed(cons));

               if( SCIPconsIsModifiable(cons) )
               {
                  *success = FALSE;
                  break;
               }

               consvars[0] = SCIPgetVarVarbound(scip, cons);
               consvars[1] = SCIPgetVbdvarVarbound(scip, cons);
               consvals[0] = 1.0;
               consvals[1] = SCIPgetVbdcoefVarbound(scip, cons);

               SCIP_CALL( addLinearConstraint(scip, matrix, consvars, consvals, 2, SCIPgetLhsVarbound(scip, cons),
                     SCIPgetRhsVarbound(scip, cons), cons) );
            }
            SCIPfreeBufferArray(scip, &consvals);
            SCIPfreeBufferArray(scip, &consvars);
         }
      }
      else if( strcmp(conshdlrname, "and") == 0 )
      {
         for( c = 0; c < nconshdlrconss; ++c )
         {
            SCIP_CONS* cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            if( SCIPconsIsModifiable(cons) )
            {
               *success = FALSE;
               break;
            }

            SCIP_CALL( addAndOrLinearization(scip, matrix, cons, SCIPgetVarsAnd(scip, cons),
                  SCIPgetNVarsAnd(scip, cons), SCIPgetResultantAnd(scip, cons), TRUE) );
         }
      }
      else if( strcmp(conshdlrname, "or") == 0 )
      {
         for( c = 0; c < nconshdlrconss; ++c )
         {
            SCIP_CONS* cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            if( SCIPconsIsModifiable(cons) )
            {
               *success = FALSE;
               break;
            }

            SCIP_CALL( addAndOrLinearization(scip, matrix, cons, SCIPgetVarsOr(scip, cons),
                                             SCIPgetNVarsOr(scip, cons), SCIPgetResultantOr(scip, cons), FALSE) );
         }
      }
      else if( strcmp(conshdlrname, "xor") == 0 )
      {
         for( c = 0; c < nconshdlrconss; ++c )
         {
            SCIP_CONS* cons = conshdlrconss[c];
            assert(SCIPconsIsTransformed(cons));

            if( SCIPconsIsModifiable(cons) )
            {
               *success = FALSE;
               break;
            }

            SCIP_CALL( addXorLinearization(scip, matrix, cons) );
         }
      }
   }

   if( *success )
   {
      SCIP_CALL( matrixSetColumnMajor(scip, matrix) );
   }
   else
   {

      SCIPfreeBufferArray(scip, &matrix->rowcons);

      SCIPfreeBufferArray(scip, &matrix->rhs);
      SCIPfreeBufferArray(scip, &matrix->lhs);
      SCIPfreeBufferArray(scip, &matrix->rowmatcnt);
      SCIPfreeBufferArray(scip, &matrix->rowmatbeg);
      SCIPfreeBufferArray(scip, &matrix->rowmatind);
      SCIPfreeBufferArray(scip, &matrix->rowmatval);

      SCIPfreeBufferArray(scip, &matrix->colintegral);
      SCIPfreeBufferArray(scip, &matrix->ub);
      SCIPfreeBufferArray(scip, &matrix->lb);
      SCIPfreeBufferArray(scip, &matrix->colmatcnt);
      SCIPfreeBufferArray(scip, &matrix->colmatbeg);
      SCIPfreeBufferArray(scip, &matrix->colmatind);
      SCIPfreeBufferArray(scip, &matrix->colmatval);
      SCIPfreeBufferArrayNull(scip, &matrix->vars);

      SCIPfreeBuffer(scip, pmatrix);
   }
   return SCIP_OKAY;
}

static
void matrixFree(
   SCIP* scip,
   IMPLINT_MATRIX** pmatrix
   )
{
   assert(scip != NULL);
   assert(pmatrix != NULL);
   IMPLINT_MATRIX* matrix = *pmatrix;
   if( matrix != NULL )
   {
      assert(matrix->colmatval != NULL);
      assert(matrix->colmatind != NULL);
      assert(matrix->colmatbeg != NULL);
      assert(matrix->colmatcnt != NULL);
      assert(matrix->lb != NULL);
      assert(matrix->ub != NULL);
      assert(matrix->colintegral != NULL);

      assert(matrix->rowmatval != NULL);
      assert(matrix->rowmatind != NULL);
      assert(matrix->rowmatbeg != NULL);
      assert(matrix->rowmatcnt != NULL);
      assert(matrix->lhs != NULL);
      assert(matrix->rhs != NULL);

      SCIPfreeBufferArray(scip, &(matrix->rowcons));

      SCIPfreeBufferArray(scip, &(matrix->rhs));
      SCIPfreeBufferArray(scip, &(matrix->lhs));
      SCIPfreeBufferArray(scip, &(matrix->rowmatcnt));
      SCIPfreeBufferArray(scip, &(matrix->rowmatbeg));
      SCIPfreeBufferArray(scip, &(matrix->rowmatind));
      SCIPfreeBufferArray(scip, &(matrix->rowmatval));

      SCIPfreeBufferArray(scip, &(matrix->colintegral));
      SCIPfreeBufferArray(scip, &(matrix->ub));
      SCIPfreeBufferArray(scip, &(matrix->lb));
      SCIPfreeBufferArray(scip, &(matrix->colmatcnt));
      SCIPfreeBufferArray(scip, &(matrix->colmatbeg));
      SCIPfreeBufferArray(scip, &(matrix->colmatind));
      SCIPfreeBufferArray(scip, &(matrix->colmatval));

      matrix->nrows = 0;
      matrix->ncols = 0;
      matrix->nnonzs = 0;

      SCIPfreeBufferArrayNull(scip, &(matrix->vars));

      SCIPfreeBuffer(scip, &matrix);
   }
}
/** Struct that contains information about the blocks/components of the submatrix given by the continuous columns.
 *
 * Note that currently, the matrix represents exactly the SCIP_MATRIX created by SCIPmatrixCreate(), but this may change
 * once MINLP problems are also accounted for.
 * @todo extend the plugin to also work for MINLP problems. This changes the computed matrix.
 */
struct MatrixComponents
{
   int nmatrixrows;                          /**< Number of rows in the matrix for the linear part of the problem */
   int nmatrixcols;                          /**< Number of columns in the matrix for the linear part of the problem */

   int* rowcomponent;                        /**< Maps a row to the index of the component it belongs to */
   int* colcomponent;                        /**< Maps a column to the index of the component it belongs to */

   int* componentrows;                       /**< Flattened array of arrays of rows that are in a given component. */
   int* componentcols;                       /**< Flattened array of arrays of columns that are in a given component. */
   int* componentrowend;                     /**< The index of componentrows where the given component ends. */
   int* componentcolend;                     /**< The index of componentcols where the given component ends. */
   int ncomponents;                          /**< The number of components. */
};
typedef struct MatrixComponents MATRIX_COMPONENTS;

/** A temporary data structure that stores some statistics/data on the rows and columns.
 *
 * This is freed again after implied integral detection is finished.
 */
struct MatrixStatistics
{
   SCIP_Bool* rowintegral;                   /**< Are all row entries of non-continuous columns and the row sides integral? */
   SCIP_Bool* rowequality;                   /**< Is the row an equality? */
   SCIP_Bool* rowbadnumerics;                /**< Does the row contain large entries that make numerics difficult? */
   int* rownnonz;                            /**< Number of nonzeros in the row */
   int* rowncontinuous;                      /**< The number of those nonzeros that are in continuous columns */
   int* rowncontinuouspmone;                 /**< The number of +-1 entries in continuous columns */
   SCIP_Bool* colintegralbounds;             /**< Does the column have integral bounds? */
};
typedef struct MatrixStatistics MATRIX_STATISTICS;

struct IntegerCandidateData{
   int column;                               /**< The candidate column to make implied integer */
   int numContPlanarEntries;                 /**< The number of nonzeros that have a row in a planar component */
   int numContNetworkEntries;                /**< The number of nonzeros that have a row in a pure network component */
   int numContTransNetworkEntries;           /**< The number of nonzeroes that have a row in a pure transposed network component */
};
typedef struct IntegerCandidateData INTEGER_CANDIDATE_DATA;

/** Creates the matrix components data structure */
static
SCIP_RETCODE createMatrixComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   IMPLINT_MATRIX *      matrix,             /**< The constraint matrix */
   MATRIX_COMPONENTS**   pmatrixcomponents   /**< Pointer to create the matrix components data structure */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, pmatrixcomponents) );
   MATRIX_COMPONENTS* comp = *pmatrixcomponents;

   int nrows = matrixGetNRows(matrix);
   int ncols = matrixGetNCols(matrix);

   comp->nmatrixrows = nrows;
   comp->nmatrixcols = ncols;

   SCIP_CALL( SCIPallocBufferArray(scip, &comp->rowcomponent, nrows) );
   for( int i = 0; i < nrows; ++i )
   {
      comp->rowcomponent[i] = -1;
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->colcomponent, ncols) );
   for( int i = 0; i < ncols; ++i )
   {
      comp->colcomponent[i] = -1;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &comp->componentrows, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->componentcols, ncols) );
   /* There will be at most ncols components */
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->componentrowend, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->componentcolend, ncols) );

   comp->ncomponents = 0;

   return SCIP_OKAY;
}

/** Frees the matrix components data structure */
static
void freeMatrixComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   MATRIX_COMPONENTS**   pmatrixcomponents   /**< Pointer to the allocated matrix components data structure */
   )
{
   MATRIX_COMPONENTS* comp = *pmatrixcomponents;

   /* Make sure to free in reverse */
   SCIPfreeBufferArray(scip, &comp->componentcolend);
   SCIPfreeBufferArray(scip, &comp->componentrowend);
   SCIPfreeBufferArray(scip, &comp->componentcols);
   SCIPfreeBufferArray(scip, &comp->componentrows);
   SCIPfreeBufferArray(scip, &comp->colcomponent);
   SCIPfreeBufferArray(scip, &comp->rowcomponent);

   SCIPfreeBuffer(scip, pmatrixcomponents);
}

/** Finds the representative of an element in the disjoint set datastructure.
 *
 * Afterwards compresses the path to speed up subsequent queries.
 */
static
int disjointSetFind(
   int*                  disjointset,        /**< The array storing the disjoint set representatives */
   int                   ind                 /**< The index to find */
   )
{
   assert(disjointset != NULL);

   int current = ind;
   int next;
   /* traverse down tree */
   while( (next = disjointset[current]) >= 0 )
   {
      current = next;
   }
   int root = current;

   /* compress indices along path */
   current = ind;
   while( (next = disjointset[current]) >= 0 )
   {
      disjointset[current] = root;
      current = next;
   }

   return root;
}

/** Merges two sets/elements into one set. Returns the index of the merged element.
 *
 * The provided elements to be merged must be representative (i.e. returned by disjointSetFind()).
 */
static
int disjointSetMerge(
   int*                  disjointset,        /**< The array storing the disjoint set representatives */
   int                   first,              /**< The first index to merge */
   int                   second              /**< The second index to merge */
   )
{
   assert(disjointset);
   assert(disjointset[first] <= -1);
   assert(disjointset[second] <= -1);
   assert(first != second); /* We cannot merge a node into itself */

   /* The rank is stored as a negative number: we decrement it making the negative number larger.
    * The rank is an upper bound on the height of the tree. We want the new root to be the one with 'largest' rank,
    * so smallest number. This way, we ensure that the tree remains shallow. If they are equal, we decrement.
    */
   int firstRank = disjointset[first];
   int secondRank = disjointset[second];
   if( firstRank > secondRank )
   {
      SCIPswapInts(&first, &second);
   }
   /* first becomes representative */
   disjointset[second] = first;
   if( firstRank == secondRank )
   {
      --disjointset[first];
   }

   return first;
}

/** Computes the continuous connected components, i.e. the connected components of the submatrix given by all
 * continuous columns.
 */
static
SCIP_RETCODE computeContinuousComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   IMPLINT_MATRIX*       matrix,             /**< The constraint matrix to compute the components for */
   MATRIX_COMPONENTS*    comp                /**< The connected components data structure to store the components in */
   )
{
   /* We let rows and columns share an index by mapping row index i to artificial column index i + nmatrixcols */
   int* disjointset = NULL;
   SCIP_CALL( SCIPallocBufferArray(scip, &disjointset, comp->nmatrixcols + comp->nmatrixrows) );
   /* First n entries belong to columns, last entries to rows */
   for( int i = 0; i < comp->nmatrixcols + comp->nmatrixrows; ++i )
   {
      disjointset[i] = -1;
   }

   for( int col = 0; col < comp->nmatrixcols; ++col )
   {
      if( matrixColIsIntegral(matrix, col) )
         continue;

      int colnnonzs = matrixGetColumnNNonz(matrix, col);
      int* colrows = matrixGetColumnIndices(matrix, col);

      int colrep = disjointSetFind(disjointset, col);
      for( int i = 0; i < colnnonzs; ++i )
      {
         int colrow = colrows[i];
         int ind = colrow + comp->nmatrixcols;
         int rowrep = disjointSetFind(disjointset, ind);
         if( colrep != rowrep )
         {
            colrep = disjointSetMerge(disjointset, colrep, rowrep);
         }
      }
   }

   /* Now, fill in the relevant data. */
   int* representativecomponent;
   SCIP_CALL( SCIPallocBufferArray(scip, &representativecomponent, comp->nmatrixcols + comp->nmatrixrows) );

   for( int i = 0; i < comp->nmatrixcols + comp->nmatrixrows; ++i )
   {
      representativecomponent[i] = -1;
   }
   comp->ncomponents = 0;
   for( int col = 0; col < comp->nmatrixcols; ++col )
   {
      if( matrixColIsIntegral(matrix, col) )
         continue;

      int colroot = disjointSetFind(disjointset, col);
      int component = representativecomponent[colroot];
      if( component < 0 )
      {
         assert(component == -1);
         /* add new component */
         component = comp->ncomponents;
         representativecomponent[colroot] = component;
         comp->componentcolend[component] = 0;
         comp->componentrowend[component] = 0;
         ++comp->ncomponents;
      }
      comp->colcomponent[col] = component;
      ++comp->componentcolend[component];
   }
   for( int row = 0; row < comp->nmatrixrows; ++row )
   {
      int rowroot = disjointSetFind(disjointset, row + comp->nmatrixcols);
      int component = representativecomponent[rowroot];
      if( component < 0 )
      {
         assert(component == -1);
         /* Any rows that have roots that we have not seen yet are rows that have no continuous columns
          * We can safely skip these for finding the continuous connected components
          */
         continue;
      }
      comp->rowcomponent[row] = component;
      ++comp->componentrowend[component];
   }
   if( comp->ncomponents >= 1 )
   {
      for( int i = 1; i < comp->ncomponents; ++i )
      {
         comp->componentrowend[i] += comp->componentrowend[i-1];
         comp->componentcolend[i] += comp->componentcolend[i-1];
      }
      int * componentnextrowindex;
      int * componentnextcolindex;
      SCIP_CALL( SCIPallocBufferArray(scip, &componentnextrowindex, comp->ncomponents) );
      SCIP_CALL( SCIPallocBufferArray(scip, &componentnextcolindex, comp->ncomponents) );
      componentnextrowindex[0] = 0;
      componentnextcolindex[0] = 0;
      for( int i = 1; i < comp->ncomponents; ++i )
      {
         componentnextcolindex[i] = comp->componentcolend[i-1];
         componentnextrowindex[i] = comp->componentrowend[i-1];
      }

      for( int col = 0; col < comp->nmatrixcols; ++col )
      {
         int component = comp->colcomponent[col];
         if( component < 0 )
         {
            assert(component == -1);
            continue;
         }
         int ind = componentnextcolindex[component];
         comp->componentcols[ind] = col;
         ++componentnextcolindex[component];
      }
      for( int row = 0; row < comp->nmatrixrows; ++row )
      {
         int component = comp->rowcomponent[row];
         if( component < 0 )
         {
            assert(component == -1);
            continue;
         }
         int ind = componentnextrowindex[component];
         comp->componentrows[ind] = row;
         ++componentnextrowindex[component];
      }

#ifndef NDEBUG
      for( int i = 0; i < comp->ncomponents; ++i )
      {
         assert(componentnextrowindex[i] == comp->componentrowend[i]);
         assert(componentnextcolindex[i] == comp->componentcolend[i]);
      }
#endif

      SCIPfreeBufferArray(scip, &componentnextcolindex);
      SCIPfreeBufferArray(scip, &componentnextrowindex);
   }

   SCIPfreeBufferArray(scip, &representativecomponent);
   SCIPfreeBufferArray(scip, &disjointset);

   return SCIP_OKAY;
}

/** Creates the matrix statistics data structure */
static
SCIP_RETCODE computeMatrixStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   IMPLINT_MATRIX*       matrix,             /**< The constraint matrix to compute the statistics for */
   MATRIX_COMPONENTS*    comp,               /**< Datastructure that contains the components of the matrix */
   MATRIX_STATISTICS**   pstats,             /**< Pointer to allocate the statistics data structure at */
   SCIP_Real             numericslimit       /**< The limit beyond which we consider integrality of coefficients
                                              *   to be unreliable */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, pstats) );
   MATRIX_STATISTICS* stats = *pstats;

   int nrows = matrixGetNRows(matrix);
   int ncols = matrixGetNCols(matrix);

   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowintegral, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowequality, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowbadnumerics, nrows) );

   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rownnonz, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowncontinuous, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &stats->rowncontinuouspmone, nrows) );

   SCIP_CALL( SCIPallocBufferArray(scip, &stats->colintegralbounds, ncols) );

   for( int i = 0; i < nrows; ++i )
   {
      SCIP_Real lhs = matrixGetRowLhs(matrix, i);
      SCIP_Real rhs = matrixGetRowRhs(matrix, i);
      int* cols = matrixGetRowIndices(matrix, i);
      SCIP_Real* vals = matrixGetRowValues(matrix, i);
      int nnonz = matrixGetRowNNonz(matrix, i);
      stats->rownnonz[i] = nnonz;
      stats->rowequality[i] = !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && SCIPisEQ(scip, lhs, rhs);

      SCIP_Bool integral = ( SCIPisInfinity(scip, -lhs) || SCIPisIntegral(scip, lhs) )
                        && ( SCIPisInfinity(scip, rhs) || SCIPisIntegral(scip, rhs) );
      SCIP_Bool badnumerics = FALSE;

      int ncontinuous = 0;
      int ncontinuouspmone = 0;
      for( int j = 0; j < nnonz; ++j )
      {
         SCIP_Bool continuous = !matrixColIsIntegral(matrix, cols[j]);
         SCIP_Real value = vals[j];
         if( continuous )
         {
            ++ncontinuous;
            if( SCIPisEQ(scip, ABS(value), 1.0) )
            {
               ++ncontinuouspmone;
            }
         }
         else
         {
            integral = integral && SCIPisIntegral(scip, value);
         }
         if( ABS(value) > numericslimit )
         {
            badnumerics = TRUE;
         }
      }

      stats->rowncontinuous[i] = ncontinuous;
      stats->rowncontinuouspmone[i] = ncontinuouspmone;
      stats->rowintegral[i] = integral;
      stats->rowbadnumerics[i] = badnumerics;
   }

   for( int i = 0; i < ncols; ++i )
   {
      SCIP_Real lb = matrixGetColLb(matrix, i);
      SCIP_Real ub = matrixGetColUb(matrix, i);
      stats->colintegralbounds[i] = ( SCIPisInfinity(scip, -lb) || SCIPisIntegral(scip, lb) )
                                 && ( SCIPisInfinity(scip, ub) || SCIPisIntegral(scip, ub) );

      /* Check that integer variables have integer bounds, as expected. */
      assert(!matrixColIsIntegral(matrix, i) || stats->colintegralbounds[i]);
   }

   return SCIP_OKAY;
}

/** Frees the matrix statistics data structure */
static
void freeMatrixStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   MATRIX_STATISTICS**   pstats              /**< Pointer to the statistics data structure to be freed */
   )
{
   MATRIX_STATISTICS* stats= *pstats;

   /* Make sure, for performance, that these frees occur in reverse */
   SCIPfreeBufferArray(scip, &stats->colintegralbounds);
   SCIPfreeBufferArray(scip, &stats->rowncontinuouspmone);
   SCIPfreeBufferArray(scip, &stats->rowncontinuous);
   SCIPfreeBufferArray(scip, &stats->rownnonz);
   SCIPfreeBufferArray(scip, &stats->rowbadnumerics);
   SCIPfreeBufferArray(scip, &stats->rowequality);
   SCIPfreeBufferArray(scip, &stats->rowintegral);

   SCIPfreeBuffer(scip, pstats);
}

/** Given the continuous components and statistics on the matrix, detect components that have implied integral variables
 * by checking if the component is a (transposed) network matrix and if all the bounds/sides/coefficients are integral.
 *
 * For every component, we detect if the associated matrix is either a network matrix or a transposed network matrix
 * (or both, in which case it represents a planar graph).
 * We choose to check if it is a (transposed) network matrix either in a row-wise or in a column-wise fashion,
 * depending on the size of the component. Finally, every variable that is in a network matrix or transposed network
 * matrix is derived to be weakly implied integral.
 */
static
SCIP_RETCODE findImpliedIntegers(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< The data belonging to the presolver */
   IMPLINT_MATRIX*       matrix,             /**< The constraint matrix to compute implied integral variables for */
   MATRIX_COMPONENTS*    comp,               /**< The continuous connected components of the matrix */
   MATRIX_STATISTICS*    stats,              /**< Statistics of the matrix */
   int*                  nchgvartypes        /**< Pointer to count the number of changed variable types */
   )
{
   /* TODO: some checks to prevent expensive memory initialization if not necessary (e.g. there must be some candidates) */
   SCIP_NETMATDEC* dec = NULL;
   SCIP_CALL( SCIPnetmatdecCreate(SCIPblkmem(scip), &dec, comp->nmatrixrows, comp->nmatrixcols) );

   SCIP_NETMATDEC* transdec = NULL;
   SCIP_CALL( SCIPnetmatdecCreate(SCIPblkmem(scip), &transdec, comp->nmatrixcols, comp->nmatrixrows) );

   /* Because the rows may also contain non-continuous columns, we need to remove these from the array that we
    * pass to the network matrix decomposition method. We use these working arrays for this purpose.
    */
   SCIP_Real* tempValArray;
   int* tempIdxArray;
   SCIP_CALL( SCIPallocBufferArray(scip, &tempValArray, comp->nmatrixcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tempIdxArray, comp->nmatrixcols) );

   SCIP_Bool * compNetworkValid;
   SCIP_Bool * compTransNetworkValid;
   SCIP_CALL(SCIPallocBufferArray(scip,&compNetworkValid,comp->ncomponents));
   SCIP_CALL(SCIPallocBufferArray(scip,&compTransNetworkValid,comp->ncomponents));

   for( int component = 0; component < comp->ncomponents; ++component )
   {
      int startrow = (component == 0) ? 0 : comp->componentrowend[component-1];
      int nrows = comp->componentrowend[component] - startrow;
      SCIP_Bool componentokay = TRUE;
      for( int i = startrow; i < startrow + nrows; ++i )
      {
         int row = comp->componentrows[i];
         if( stats->rowncontinuous[row] != stats->rowncontinuouspmone[row] || !stats->rowintegral[row] )
         {
            componentokay = FALSE;
            break;
         }
         if( stats->rowbadnumerics[row] )
         {
            componentokay = FALSE;
            break;
         }
      }
      int startcol = ( component == 0 ) ? 0 : comp->componentcolend[component-1];
      int ncols = comp->componentcolend[component] - startcol;

      for( int i = startcol; i < startcol + ncols && componentokay; ++i )
      {
         int col = comp->componentcols[i];
         if( !stats->colintegralbounds[col] )
         {
            componentokay = FALSE;
            break;
         }
      }

      if( !componentokay )
      {
         compNetworkValid[component] = FALSE;
         compTransNetworkValid[component] = FALSE;
         continue;
      }

      /* Check if the component is a network matrix */
      SCIP_Bool componentnetwork = TRUE;

      /* We use the row-wise algorithm only if the number of columns is much larger than the number of rows.
       * Generally, the column-wise algorithm will be faster, but in these extreme cases, the row algorithm is faster.
       * Only very few instances should use the row-wise algorithm.
       */
      if( nrows * presoldata->columnrowratio < ncols )
      {
         for( int i = startrow; i < startrow + nrows && componentnetwork; ++i )
         {
            int row = comp->componentrows[i];
            int nrownnoz = matrixGetRowNNonz(matrix, row);
            int* rowcols = matrixGetRowIndices(matrix, row);
            SCIP_Real* rowvals = matrixGetRowValues(matrix, row);
            int ncontnonz = 0;
            for( int j = 0; j < nrownnoz; ++j )
            {
               int col = rowcols[j];
               if( !matrixColIsIntegral(matrix, col) )
               {
                  tempIdxArray[ncontnonz] = col;
                  tempValArray[ncontnonz] = rowvals[j];
                  ++ncontnonz;
                  assert(SCIPisEQ(scip, ABS(rowvals[j]), 1.0));
               }
            }

            SCIP_CALL( SCIPnetmatdecTryAddRow(dec, row, tempIdxArray, tempValArray, ncontnonz, &componentnetwork) );
         }
      }
      else
      {
         for( int i = startcol; i < startcol + ncols && componentnetwork; ++i )
         {
            int col = comp->componentcols[i];
            int ncolnnonz = matrixGetColumnNNonz(matrix, col);
            int* colrows = matrixGetColumnIndices(matrix, col);
            SCIP_Real* colvals = matrixGetColumnValues(matrix, col);
            SCIP_CALL( SCIPnetmatdecTryAddCol(dec, col, colrows, colvals, ncolnnonz, &componentnetwork) );
         }
      }

      if( !componentnetwork )
      {
         SCIPnetmatdecRemoveComponent(dec, &comp->componentrows[startrow], nrows, &comp->componentcols[startcol], ncols);
      }

      compNetworkValid[component] = componentnetwork;

      /* Avoid redundant work; we don't need to check if component is both network and transposed network in the case
       * where do not want to extend implied integrality to the integers */
      if( !presoldata->convertintegers && componentnetwork ){
         continue;
      }

      SCIP_Bool componenttransnetwork = TRUE;

      /* For the transposed matrix, the situation is exactly reversed because the row/column algorithms are swapped */
      if( nrows <= ncols * presoldata->columnrowratio )
      {
         for( int i = startrow; i < startrow + nrows && componenttransnetwork; ++i )
         {
            int row = comp->componentrows[i];
            int nrownnoz = matrixGetRowNNonz(matrix, row);
            int* rowcols = matrixGetRowIndices(matrix, row);
            SCIP_Real* rowvals = matrixGetRowValues(matrix, row);
            int ncontnonz = 0;
            for( int j = 0; j < nrownnoz; ++j )
            {
               int col = rowcols[j];
               if( !matrixColIsIntegral(matrix, col) )
               {
                  tempIdxArray[ncontnonz] = col;
                  tempValArray[ncontnonz] = rowvals[j];
                  ++ncontnonz;
                  assert(SCIPisEQ(scip, ABS(rowvals[j]), 1.0));
               }
            }

            SCIP_CALL( SCIPnetmatdecTryAddCol(transdec, row, tempIdxArray, tempValArray, ncontnonz,
                                              &componenttransnetwork) );
         }
      }
      else
      {
         for( int i = startcol; i < startcol + ncols && componenttransnetwork; ++i )
         {
            int col = comp->componentcols[i];
            int ncolnnonz = matrixGetColumnNNonz(matrix, col);
            int* colrows = matrixGetColumnIndices(matrix, col);
            SCIP_Real* colvals = matrixGetColumnValues(matrix, col);
            SCIP_CALL( SCIPnetmatdecTryAddRow(transdec, col, colrows, colvals, ncolnnonz, &componenttransnetwork) );
         }
      }

      if( !componenttransnetwork )
      {
         SCIPnetmatdecRemoveComponent(transdec, &comp->componentcols[startcol], ncols,
                                      &comp->componentrows[startrow], nrows);
      }

      compTransNetworkValid[component] = componenttransnetwork;
   }

   if( presoldata->convertintegers )
   {

      int numCandidates = 0;
      INTEGER_CANDIDATE_DATA* candidates;
      SCIP_CALL(SCIPallocBufferArray(scip, &candidates, comp->nmatrixcols));

      /* Integer columns that are +- 1 and do not have any nonzero entries in bad rows are candidates */
      for( int col = 0; col < comp->nmatrixcols; ++col )
      {
         if( !matrixColIsIntegral(matrix, col) )
            continue;

         int ncolnnonz = matrixGetColumnNNonz(matrix, col);
         int* colrows = matrixGetColumnIndices(matrix, col);
         SCIP_Real* colvals = matrixGetColumnValues(matrix, col);
         SCIP_Bool badColumn = FALSE;
         INTEGER_CANDIDATE_DATA* data = &candidates[numCandidates];
         data->numContNetworkEntries = 0;
         data->numContPlanarEntries = 0;
         data->numContTransNetworkEntries = 0;
         data->column = col;
         for( int i = 0; i < ncolnnonz; ++i )
         {
            int entryrow = colrows[i];
            if( !stats->rowintegral[entryrow] || stats->rowbadnumerics[entryrow] ||
                !SCIPisEQ(scip, ABS(colvals[i]), 1.0))
            {
               badColumn = TRUE;
               break;
            }
            int rowcomponent = comp->rowcomponent[entryrow];
            if( rowcomponent != -1 )
            {
               SCIP_Bool networkValid = compNetworkValid[rowcomponent];
               SCIP_Bool transNetworkValid = compTransNetworkValid[rowcomponent];
               if( networkValid && transNetworkValid )
               {
                  data->numContPlanarEntries++;
               }
               else if( networkValid )
               {
                  data->numContNetworkEntries++;
               }
               else if( transNetworkValid )
               {
                  data->numContTransNetworkEntries++;
               }
               else
               {
                  badColumn = TRUE;
                  break;
               }
            }
         }
         if( badColumn )
            continue;
         ++numCandidates;
      }

      SCIP_Real* candidateScores;
      SCIP_CALL(SCIPallocBufferArray(scip, &candidateScores, numCandidates));

      for( int i = 0; i < numCandidates; ++i )
      {
         int col = candidates[i].column;
         int nnonzs = matrixGetColumnNNonz(matrix, col);

         /* Higher score; we prefer to pick this variable first. We generally prefer to detect implied integrality of
          * general integer variables over binary variables. We break ties using the number of nonzeros in the column
          * @TODO; test different scores / alternatives
          * @TODO; detect when we only have columns that  extend the network / co-network portions. In this case, we can take both at the same time.*/
         double score = 0.0;
         switch( SCIPvarGetType(matrixGetVar(matrix, col)) )
         {
            case SCIP_VARTYPE_BINARY:
               score += 10.0;
               break;
            case SCIP_VARTYPE_INTEGER:
               score += 100.0;
               break;
            case SCIP_VARTYPE_CONTINUOUS:
               break;
            default:
               SCIPerrorMessage("unknown variable type\n");
               return SCIP_INVALIDDATA;
         } /*lint !e788*/
         score += nnonzs * 0.1;
         candidateScores[i] = -score;
      }
      INTEGER_CANDIDATE_DATA** ptrArray;
      SCIP_CALL(SCIPallocBufferArray(scip, &ptrArray, numCandidates));
      for( int i = 0; i < numCandidates; ++i )
      {
         ptrArray[i] = &candidates[i];
      }
      SCIPsortDownRealPtr(candidateScores, (void**) ptrArray, numCandidates);

      int integerNetwork = 0;
      int integerTransNetwork = 0;

      for( int i = 0; i < numCandidates; ++i )
      {
         INTEGER_CANDIDATE_DATA* candidate = ptrArray[i];
         if( candidate->numContTransNetworkEntries == 0 )
         {
            int col = candidate->column;
            int* colrows = matrixGetColumnIndices(matrix, col);
            SCIP_Real* colvals = matrixGetColumnValues(matrix, col);
            int ncolnnonz = matrixGetColumnNNonz(matrix, col);
            SCIP_Bool success;
            SCIP_CALL(SCIPnetmatdecTryAddCol(dec, col, colrows, colvals, ncolnnonz, &success));
            if( success )
               integerNetwork++;
         }

         if( candidate->numContNetworkEntries == 0 )
         {
            int col = candidate->column;
            int* colrows = matrixGetColumnIndices(matrix, col);
            SCIP_Real* colvals = matrixGetColumnValues(matrix, col);
            int ncolnnonz = matrixGetColumnNNonz(matrix, col);
            SCIP_Bool success;
            SCIP_CALL(SCIPnetmatdecTryAddRow(transdec, col, colrows, colvals, ncolnnonz, &success));
            if( success )
               integerTransNetwork++;
         }
      }
      if( integerNetwork >= integerTransNetwork )
      {
         /* We add all integer columns from the network matrix */
         for( int i = 0; i < numCandidates; ++i )
         {
            int col = ptrArray[i]->column;
            if( SCIPnetmatdecContainsColumn(dec, col))
            {
               SCIP_VAR* var = matrixGetVar(matrix, col);
               SCIP_Bool infeasible = FALSE;
               SCIP_CALL( SCIPchgVarImplType(scip, var, SCIP_IMPLINTTYPE_WEAK, &infeasible) );

               ( *nchgvartypes )++;
               assert(!infeasible);
            }
         }

      }
      else
      {
         /* We add all integer columns from the transposed network matrix */
         for( int i = 0; i < numCandidates; ++i )
         {
            int col = ptrArray[i]->column;
            if( SCIPnetmatdecContainsRow(transdec, col))
            {
               SCIP_VAR* var = matrixGetVar(matrix, col);
               SCIP_Bool infeasible = FALSE;
               SCIP_CALL( SCIPchgVarImplType(scip, var, SCIP_IMPLINTTYPE_WEAK, &infeasible) );

               ( *nchgvartypes )++;
               assert(!infeasible);
            }
         }
      }

      SCIPfreeBufferArray(scip, &ptrArray);
      SCIPfreeBufferArray(scip, &candidateScores);
      SCIPfreeBufferArray(scip, &candidates);

   }
   /* Add continuous columns; here we can take both normal and transposed components at the same time. */
   for( int component = 0; component < comp->ncomponents; ++component )
   {
      if( !compNetworkValid[component] && !compTransNetworkValid[component] )
         continue;
      int startcol = ( component == 0 ) ? 0 : comp->componentcolend[component - 1];
      int ncols = comp->componentcolend[component] - startcol;
      for( int i = startcol; i < startcol + ncols; ++i )
      {
         int col = comp->componentcols[i];
         assert(SCIPnetmatdecContainsColumn(dec, col) || SCIPnetmatdecContainsRow(transdec, col));
         SCIP_VAR* var = matrixGetVar(matrix, col);
         SCIP_Bool infeasible = FALSE;
         SCIP_CALL( SCIPchgVarImplType(scip, var, SCIP_IMPLINTTYPE_WEAK, &infeasible) );
         (*nchgvartypes)++;
         assert(!infeasible);
      }
   }

   SCIPfreeBufferArray(scip, &compTransNetworkValid);
   SCIPfreeBufferArray(scip, &compNetworkValid);

   SCIPfreeBufferArray(scip, &tempIdxArray);
   SCIPfreeBufferArray(scip, &tempValArray);

   SCIPnetmatdecFree(&transdec);
   SCIPnetmatdecFree(&dec);

   return SCIP_OKAY;
}

/*
 * Callback methods of presolver
 */

/** copy method for presolver plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PRESOLCOPY(presolCopyImplint)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolImplint(scip) );

   return SCIP_OKAY;
}

/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeImplint)
{
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeBlockMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}

/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecImplint)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;

   /* TODO: re-check these conditions again */
   /* Disable implicit integer detection if we are probing or in NLP context */
   if( ( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING ) || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
   {
      return SCIP_OKAY;
   }
   /* Since implied integrality detection relies on rows being static, we disable it for branch-and-price applications*/
   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 )
   {
      return SCIP_OKAY;
   }

   /* Only run if we would otherwise terminate presolving */
   if( !SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   /* For now, we don't re-run detection for primal heuristics / subscips */
   if( SCIPgetSubscipDepth(scip) > 0 )
      return SCIP_OKAY;

   /* Exit early if we already ran or if there are no variables to upgrade */
   SCIP_PRESOLDATA* presoldata = SCIPpresolGetData(presol);
   if( presoldata->computedimplints || (!presoldata->convertintegers && SCIPgetNContVars(scip) == 0) )
   {
      return SCIP_OKAY;
   }

   presoldata->computedimplints = TRUE;

   *result = SCIP_DIDNOTFIND;

   SCIP_Real starttime = SCIPgetSolvingTime(scip);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
                   "   (%.1fs) implied integrality detection started\n", starttime);

   SCIP_Bool success = TRUE;
   IMPLINT_MATRIX* matrix = NULL;
   SCIP_CALL( matrixCreate(scip, &matrix, &success) );
   if( !success )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
                      "   (%.1fs) implied integrality detection stopped because problem is not an MILP\n",
                      SCIPgetSolvingTime(scip));
      return SCIP_OKAY;
   }

   int beforechanged = *nchgvartypes;
   MATRIX_COMPONENTS* comp = NULL;
   MATRIX_STATISTICS* stats = NULL;
   SCIP_CALL( createMatrixComponents(scip, matrix, &comp) );
   SCIP_CALL( computeMatrixStatistics(scip, matrix, comp, &stats, presoldata->numericslimit) );
   SCIP_CALL( computeContinuousComponents(scip, matrix, comp) );
   SCIP_CALL( findImpliedIntegers(scip, presoldata, matrix, comp, stats, nchgvartypes) );
   int afterchanged = *nchgvartypes;

   SCIP_Real endtime = SCIPgetSolvingTime(scip);
   if( afterchanged == beforechanged )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   (%.1fs) no implied integral variables detected (time: %.2fs)\n",
            endtime, endtime - starttime);
      *result = SCIP_DIDNOTFIND;
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   (%.1fs) %d implied integral variables detected (time: %.2fs)\n",
            endtime, afterchanged-beforechanged,endtime-starttime);

      *result = SCIP_SUCCESS;
   }
   freeMatrixStatistics(scip, &stats);
   freeMatrixComponents(scip, &comp);
   matrixFree(scip, &matrix);

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the implint presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolImplint(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create implint presolver data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &presoldata) );

   presoldata->computedimplints = FALSE;

   /* include implint presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
                                     PRESOL_TIMING, presolExecImplint, presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyImplint) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeImplint) );

   SCIP_CALL( SCIPaddRealParam(scip,
                               "presolving/implint/columnrowratio",
                               "use the network row addition algorithm when the column to row ratio becomes larger than this threshold",
                               &presoldata->columnrowratio, TRUE, DEFAULT_COLUMNROWRATIO, 0.0, 1e98, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
                               "presolving/implint/numericslimit",
                               "a row that contains variables with coefficients that are greater in absolute value than this limit is not considered for implied integrality detection",
                               &presoldata->numericslimit, TRUE, DEFAULT_NUMERICSLIMIT, 1.0, 1e98, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
                               "presolving/implint/convertintegers",
                               "Controls whether implied integrality is inferred for integer variables",
                               &presoldata->convertintegers, FALSE, DEFAULT_CONVERTINTEGERS, NULL, NULL) );
   return SCIP_OKAY;
}
