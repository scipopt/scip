/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/* TODO: support more constraint types: cons_nonlinear, cons_indicator and symmetry constraints */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_implint.h"
#include "scip/pub_cons.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_network.h"
#include "scip/pub_presol.h"
#include "scip/pub_var.h"

#include "scip/scip_cons.h"
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

#include "scip/cons_and.h"
#include "scip/cons_linear.h"
#include "scip/cons_logicor.h"
#include "scip/cons_knapsack.h"
#include "scip/cons_or.h"
#include "scip/cons_setppc.h"
#include "scip/cons_varbound.h"
#include "scip/cons_xor.h"

#define PRESOL_NAME             "implint"
#define PRESOL_DESC             "detects implicit integer variables"

/* We want to run as late as possible, but before symmetry detection.
 * The main reason for this is that symmetry detection may add linear constraints that
 * impede the detection of implied integrality, but do not break implied integrality itself.
 * Also, symmetry methods rely on the fact that each variable in an orbit is integral,
 * as otherwise certain reductions may break. So it is currently not safe to run implied integrality detection
 * after symmetry methods are applied. */
#define PRESOL_PRIORITY         -900000      /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS        0            /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_EXHAUSTIVE /* timing of the presolver (fast, medium, or exhaustive) */

#define DEFAULT_CONVERTINTEGERS FALSE        /**< should implied integrality also be detected for enforced integral variables? */
#define DEFAULT_COLUMNROWRATIO  50.0         /**< use the network row addition algorithm when the column to row ratio becomes larger than this threshold */
#define DEFAULT_NUMERICSLIMIT   1e8          /**< a row that contains variables with coefficients that are greater in absolute value than this limit is not considered for implied integrality detection */

/** presolver data */
struct SCIP_PresolData
{
   SCIP_Bool             computedimplints;   /**< were implied integers already computed? */
   SCIP_Bool             convertintegers;    /**< should implied integrality also be detected for enforced integral variables? */
   SCIP_Real             columnrowratio;     /**< use the network row addition algorithm when the column to row ratio
                                              *   becomes larger than this threshold, otherwise, use column addition */
   SCIP_Real             numericslimit;      /**< a row that contains variables with coefficients that are greater in
                                              *   absolute value than this limit is not considered for
                                              *   implied integrality detection */
};

/** constraint matrix data structure in column and row major format
 *  Contains only the linear terms, and marks the presence of non-linear terms.
 */
struct ImplintMatrix
{
   SCIP_Real*            colmatval;          /**< coefficients in column major format */
   int*                  colmatind;          /**< row indexes in column major format */
   int*                  colmatbeg;          /**< column storage offset */
   int*                  colmatcnt;          /**< number of row entries per column */
   int                   ncols;              /**< complete number of columns */
   SCIP_Real*            lb;                 /**< lower bound per variable */
   SCIP_Real*            ub;                 /**< upper bound per variable */
   SCIP_Bool*            colintegral;        /**< whether column is integral */
   SCIP_Bool*            colimplintegral;    /**< whether the column is implied integral */
   SCIP_Bool*            colinnonlinterm;    /**< is the column involved in some nonlinear term? */
   /* TODO: fields for more involved detection and scoring:
    * bounds integral? number of +-1 nonzeros?
    * ntimes operand / resultant in logical constraints?
    * nconstraints (different from nnonz because of multiple row constraints)
    * npmonenonzeros in integral equality rows */

   SCIP_VAR**            colvar;             /**< variable described by column */

   SCIP_Real*            rowmatval;          /**< coefficients in row major format */
   int*                  rowmatind;          /**< column indexed in row major format */
   int*                  rowmatbeg;          /**< row storage offset */
   int*                  rowmatcnt;          /**< number of column entries per row */

   int                   nrows;              /**< complete number of rows */
   SCIP_Real*            lhs;                /**< left hand side per row */
   SCIP_Real*            rhs;                /**< right hand side per row */

   SCIP_CONS**           rowcons;            /**< constraint described by row */

   int                   nnonzs;             /**< sparsity counter */
   int                   nnonzssize;         /**< size of the nonzero arrays */
};
typedef struct ImplintMatrix IMPLINT_MATRIX;

/** struct that contains information about the blocks/components of the submatrix given by the continuous columns */
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

/** a temporary data structure that stores some statistics/data on the rows and columns */
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

/** struct that contains some information for each integer variable that is a candidate for implied integrality detection */
struct IntegerCandidateData
{
   int column;                               /**< The candidate column to make implied integer */
   int numContPlanarEntries;                 /**< The number of nonzeros that have a row in a planar component */
   int numContNetworkEntries;                /**< The number of nonzeros that have a row in a pure network component */
   int numContTransNetworkEntries;           /**< The number of nonzeroes that have a row in a pure transposed network component */
};
typedef struct IntegerCandidateData INTEGER_CANDIDATE_DATA;

/** gets a pointer to the array of nonzero values for the nonzeros in the given column */
static
SCIP_Real* matrixGetColumnVals(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->colmatval + matrix->colmatbeg[column];
}

/** gets a pointer to the array of row indices for the nonzeros in the given column */
static
int* matrixGetColumnInds(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->colmatind + matrix->colmatbeg[column];
}

/** gets the number of nonzeros in the given column */
static
int matrixGetColumnNNonzs(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->colmatcnt[column];
}

/** gets a pointer to the array of nonzero values for the nonzeros in the given row */
static
SCIP_Real* matrixGetRowVals(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   row                 /**< the row */
   )
{
   assert(matrix != NULL);
   assert(row >= 0);
   assert(row < matrix->nrows);

   return matrix->rowmatval + matrix->rowmatbeg[row];
}

/** gets a pointer to the array of column indices for the nonzeros in the given row */
static
int* matrixGetRowInds(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   row                 /**< the row */
   )
{
   assert(matrix != NULL);
   assert(row >= 0);
   assert(row < matrix->nrows);

   return matrix->rowmatind + matrix->rowmatbeg[row];
}

/** gets the number of nonzeros in the given row */
static
int matrixGetRowNNonzs(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   row                 /**< the row */
   )
{
   assert(matrix != NULL);
   assert(row >= 0);
   assert(row < matrix->nrows);

   return matrix->rowmatcnt[row];
}

/** returns the number of rows in the matrix */
static
int matrixGetNRows(
   IMPLINT_MATRIX*       matrix              /**< the matrix data structure */
   )
{
   assert(matrix != NULL);

   return matrix->nrows;
}

/** returns the number of columns in the matrix */
static
int matrixGetNCols(
   IMPLINT_MATRIX*       matrix              /**< the matrix data structure */
   )
{
   assert(matrix != NULL);

   return matrix->ncols;
}

/** returns the variable associated with the column */
static
SCIP_VAR* matrixGetVar(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->colvar[column];
}

/** returns TRUE if the given column originates from an integral variable */
static
SCIP_Bool matrixColIsIntegral(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->colintegral[column];
}

/** returns TRUE if the given column originates from an implied integral variable */
static
SCIP_Bool matrixColIsImpliedIntegral(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->colimplintegral[column];
}

/** returns TRUE if the given column occurs in a nonlinear expression in some constraint */
static
SCIP_Bool matrixColInNonlinearTerm(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->colinnonlinterm[column];
}

/** returns the lower bound of the given column */
static
SCIP_Real matrixGetColLb(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->lb[column];
}

/** returns the upper bound of the given column */
static
SCIP_Real matrixGetColUb(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   column              /**< the column */
   )
{
   assert(matrix != NULL);
   assert(column >= 0);
   assert(column < matrix->ncols);

   return matrix->ub[column];
}

/** returns the left hand side of the given row */
static
SCIP_Real matrixGetRowLhs(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   row                 /**< the row */
   )
{
   assert(matrix != NULL);
   assert(row >= 0);
   assert(row < matrix->nrows);

   return matrix->lhs[row];
}

/** returns the right hand side of the given row */
static
SCIP_Real matrixGetRowRhs(
   IMPLINT_MATRIX*       matrix,             /**< the matrix data structure */
   int                   row                 /**< the row */
   )
{
   assert(matrix != NULL);
   assert(row >= 0);
   assert(row < matrix->nrows);

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
   int probindex;
   int rowidx;
   int j;

   assert(vars != NULL);
   assert(vals != NULL);

   rowidx = matrix->nrows;

   matrix->lhs[rowidx] = lhs;
   matrix->rhs[rowidx] = rhs;
   matrix->rowmatbeg[rowidx] = matrix->nnonzs;

   for( j = 0; j < nvars; ++j )
   {
      /* ignore variables with very small coefficients */
      if( SCIPisZero(scip, vals[j]) )
         continue;

      assert(matrix->nnonzs < matrix->nnonzssize);
      matrix->rowmatval[matrix->nnonzs] = vals[j];
      probindex = SCIPvarGetProbindex(vars[j]);
      assert(0 <= probindex && probindex < matrix->ncols);
      matrix->rowmatind[matrix->nnonzs] = probindex;
      ++matrix->nnonzs;
   }

   matrix->rowmatcnt[rowidx] = matrix->nnonzs - matrix->rowmatbeg[rowidx];
   matrix->rowcons[rowidx] = cons;

   ++matrix->nrows;

   return SCIP_OKAY;
}

/** transforms the weighted sum to active variables and then adds the given linear constraint to the matrix */
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
   assert(nvars >= 1 || (!SCIPisPositive(scip, lhs) && !SCIPisNegative(scip, rhs)));

   /* constraint is redundant */
   if( nvars == 0 || ( SCIPisInfinity(scip, -lhs) && SCIPisInfinity(scip, rhs) ) )
      return SCIP_OKAY;

   activevars = NULL;
   activevals = NULL;
   nactivevars = nvars;
   activeconstant = 0.0;

   /* duplicate variable and value array */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &activevars, vars, nactivevars) );
   if( vals != NULL )
   {
      SCIP_CALL( SCIPduplicateBufferArray(scip, &activevals, vals, nactivevars) );
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &activevals, nactivevars) );

      for( v = 0; v < nactivevars; ++v )
         activevals[v] = 1.0;
   }

   /* retransform given variables to active variables */
   SCIP_CALL( getActiveVariables(scip, &activevars, &activevals, &nactivevars, &activeconstant) );

   /* adapt left and right hand side */
   if( !SCIPisInfinity(scip, -lhs) )
      lhs -= activeconstant;
   if( !SCIPisInfinity(scip, rhs) )
      rhs -= activeconstant;

   assert(nactivevars >= 1 || (!SCIPisPositive(scip, lhs) && !SCIPisNegative(scip, rhs)));

   /* add single row to matrix */
   if( nactivevars >= 1 )
   {
      /**@todo normalize by greatest common divisor of coefficients for integral columns */
      SCIP_CALL( matrixAddRow(scip, matrix, activevars, activevals, nactivevars, lhs, rhs, cons) );
   }

   /* free buffer arrays */
   SCIPfreeBufferArray(scip, &activevals);
   SCIPfreeBufferArray(scip, &activevars);

   return SCIP_OKAY;
}

/** adds the linearization of a given AND constraint or OR constraint to the constraint matrix */
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

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, noperands + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, noperands + 1) );

   /* add all the constraints of the form resultant <= operand */
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

   vars[0] = resultant;
   vals[0] = 1.0;
   vals[1] = -1.0;

   for( i = 0; i < noperands; ++i )
   {
      vars[1] = operands[i];

      SCIP_CALL( addLinearConstraint(scip, matrix, vars, vals, 2, lhs, rhs, cons) );
   }

   /* add the constraint of the form noperands - 1 + resultant >= sum operands */
   if( isAndCons )
   {
      lhs = 1.0 - noperands;
      rhs = SCIPinfinity(scip);
   }
   else
   {
      lhs = -SCIPinfinity(scip);
      rhs = 0.0;
   }

   for( i = 0; i < noperands; ++i )
   {
      vars[i + 1] = operands[i];
      vals[i + 1] = -1.0;
   }

   SCIP_CALL( addLinearConstraint(scip, matrix, vars, vals, noperands + 1, lhs, rhs, cons) );

   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** adds the linearization of a given XOR constraint to the constraint matrix */
static
SCIP_RETCODE addXorLinearization(
   SCIP*                 scip,               /**< current scip instance */
   IMPLINT_MATRIX*       matrix,             /**< constraint matrix */
   SCIP_CONS*            cons,               /**< The constraint that is linearized */
   SCIP_VAR**            operands,           /**< variables of this constraint */
   int                   noperands,          /**< number of operands */
   SCIP_VAR*             intvar,             /**< the intvar of the xor constraint */
   SCIP_Real             rhs                 /**< the right hand side of the xor constraint */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* vals;
   int i;
   int j;
   int k;

   SCIP_CALL( SCIPallocBufferArray(scip, &vals, noperands + 1) );

   if( intvar != NULL )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, noperands + 1) );

      /* add intvar constraint */
      for( j = 0; j < noperands; ++j )
      {
         vars[j] = operands[j];
         vals[j] = 1.0;
      }
      vars[noperands] = intvar;
      vals[noperands] = -2.0;

      SCIP_CALL( addLinearConstraint(scip, matrix, vars, vals, noperands + 1, rhs, rhs, cons) );

      SCIPfreeBufferArray(scip, &vars);
   }
   else if( noperands == 3 )
   {
      /*  in the special case of 3 variables and c = 0, the following linear system is created:
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
      SCIP_Real scale = rhs == 0.0 ? 1.0 : -1.0; /*lint !e777*/

      for( i = 0; i < noperands; ++i )
      {
         for( j = 0; j < noperands; ++j )
            vals[j] = (i == j) ? scale : -scale;

         SCIP_CALL( addLinearConstraint(scip, matrix, operands, vals, noperands, -SCIPinfinity(scip), rhs, cons) );
      }

      for( j = 0; j < noperands; ++j )
         vals[j] = scale;

      SCIP_CALL( addLinearConstraint(scip, matrix, operands, vals, noperands, -SCIPinfinity(scip), 2.0 - rhs * noperands, cons) );
   }
   else if( noperands < 3 )
   {
      for( j = 0; j < noperands; ++j )
         vals[j] = (j <= rhs) ? 1.0 : -1.0;

      SCIP_CALL( addLinearConstraint(scip, matrix, operands, vals, noperands, rhs, rhs, cons) );
   }
   else
   {
      /* long XOR constraints are represented nonlinearly, so the relevant active variables are marked non-linear */
      SCIP_VAR** aggrvars;
      SCIP_Real* scalars;
      SCIP_Real constant;
      int naggrvars;
      int col;

      SCIP_CALL( SCIPallocBufferArray(scip, &scalars, matrix->ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &aggrvars, matrix->ncols) );

      /* we can transform all variables together in the sum because the constraint can be interpreted as a nonlinear
       * constraint of the form sum(operands) % 2 = rhs. Thus, it is okay if we have cancellations in the sum.
       */
      for( j = 0; j < noperands; ++j )
      {
         scalars[j] = 1.0;
         aggrvars[j] = operands[j];
      }

      naggrvars = noperands;
      SCIP_CALL( getActiveVariables(scip, &aggrvars, &scalars, &naggrvars, &constant) );

      for( k = 0; k < naggrvars; ++k )
      {
         /* if the variable has an even coefficient, it does not contribute to the modulo constraint */
         if( !SCIPisIntegral(scip, 0.5 * scalars[k]) )
         {
            col = SCIPvarGetProbindex(aggrvars[k]);
            assert(col >= 0);
            assert(col < matrix->ncols);
            matrix->colinnonlinterm[col] = TRUE;
         }
      }

      SCIPfreeBufferArray(scip, &aggrvars);
      SCIPfreeBufferArray(scip, &scalars);
   }

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
   SCIP_Real* valpnt;
   int* rowpnt;
   int* rowend;
   int* fillidx;
   int colidx;
   int i;

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

   for( i = 0; i < matrix->nrows; ++i )
   {
      rowpnt = matrix->rowmatind + matrix->rowmatbeg[i];
      rowend = rowpnt + matrix->rowmatcnt[i];
      for( ; rowpnt < rowend; ++rowpnt )
      {
         colidx = *rowpnt;
         ++matrix->colmatcnt[colidx];
      }
   }

   matrix->colmatbeg[0] = 0;
   for( i = 0; i < matrix->ncols - 1; ++i )
      matrix->colmatbeg[i+1] = matrix->colmatbeg[i] + matrix->colmatcnt[i];

   for( i = 0; i < matrix->nrows; ++i )
   {
      rowpnt = matrix->rowmatind + matrix->rowmatbeg[i];
      rowend = rowpnt + matrix->rowmatcnt[i];
      valpnt = matrix->rowmatval + matrix->rowmatbeg[i];

      for( ; rowpnt != rowend; ++rowpnt, ++valpnt )
      {
         colidx = *rowpnt;
         assert(colidx < matrix->ncols);
         matrix->colmatind[matrix->colmatbeg[colidx] + fillidx[colidx]] = i;
         matrix->colmatval[matrix->colmatbeg[colidx] + fillidx[colidx]] = *valpnt;
         ++fillidx[colidx];
      }
   }

   SCIPfreeBufferArray(scip, &fillidx);

   return SCIP_OKAY;
}

/* @todo: skip construction of integral constraints if we do not run detection on integer variables */
/* @todo: use symmetry constraints to guide variable ordering for integral columns because
 *        symmetric variables should always all be either network or non-network
 */
/** create the matrix from the current transformed problem */
static
SCIP_RETCODE matrixCreate(
   SCIP*                 scip,               /**< the scip data structure */
   IMPLINT_MATRIX**      pmatrix             /**< pointer to create the matrix at */
   )
{
   SCIP_CONSHDLR** conshdlrs;
   SCIP_VAR** vars;
   IMPLINT_MATRIX* matrix;
   SCIP_Bool success;
   const char* conshdlrname;
   int nconshdlrs;
   int nmatrixrows;
   int nconshdlrconss;
   int nvars;
   int nnonzstmp;
   int i;
   int j;

   *pmatrix = NULL;

   /* return if no variables or constraints are present */
   if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
      return SCIP_OKAY;

   assert(SCIPgetNActivePricers(scip) == 0);

   /* loop over all constraint handlers and collect the number of checked constraints that contribute rows
    * to the matrix */
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
            nmatrixrows += nconshdlrconss;
         else if( strcmp(conshdlrname, "and") == 0 )
         {
            /* the linearization of AND constraints is modelled using n + 1 inequalities on n + 1 variables */
            SCIP_CONS** checked = SCIPconshdlrGetCheckConss(conshdlrs[i]);
            for( j = 0; j < nconshdlrconss; ++j )
            {
               int nandvars = SCIPgetNVarsAnd(scip, checked[j]);
               nmatrixrows += nandvars + 1;
               if( nandvars > 1 )
                  nnonzstmp += nandvars - 1;
            }
         }
         else if( strcmp(conshdlrname, "or") == 0 )
         {
            /* the linearization of OR constraints is modelled using n + 1 inequalities on n + 1 variables */
            SCIP_CONS** checked = SCIPconshdlrGetCheckConss(conshdlrs[i]);
            for( j = 0; j < nconshdlrconss; ++j )
            {
               int norvars = SCIPgetNVarsOr(scip, checked[j]);
               nmatrixrows += norvars + 1;
               if( norvars > 1 )
                  nnonzstmp += norvars - 1;
            }
         }
         else if( strcmp(conshdlrname, "xor") == 0 )
         {
            /* the relaxation of XOR constraints is handled differently depending on the integer variable and size:
             * with integer variable or less than three variables, the representation is a single row;
             * without integer variable and three variables, the convex hull of the constraint is added with four rows;
             * otherwise, the constraint is considered nonlinear because the convex hull representation is exponential
             */
            SCIP_CONS** checked = SCIPconshdlrGetCheckConss(conshdlrs[i]);
            for( j = 0; j < nconshdlrconss; ++j )
            {
               int nxorvars = SCIPgetNVarsXor(scip, checked[j]);
               if( SCIPgetIntVarXor(scip, checked[j]) != NULL || nxorvars < 3 )
                  nmatrixrows += 1;
               else if( nxorvars == 3 )
               {
                  nmatrixrows += 4;
                  nnonzstmp += 6;
               }
            }
         }
         else
         {
            /* @todo: support symmetry, linking, sos1, sos2, bounddisjunction, nonlinear, indicator, superindicator conshdlrs */
            return SCIP_OKAY;
         }
      }
   }

   if( nmatrixrows == 0 )
      return SCIP_OKAY;

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* approximate number of nonzeros by taking for each variable the number of down- and uplocks;
    * this counts nonzeros in equalities twice, but can be at most two times as high as the exact number
    */
   for( i = 0; i < nvars; ++i )
      nnonzstmp += SCIPvarGetNLocksDownType(vars[i], SCIP_LOCKTYPE_MODEL) + SCIPvarGetNLocksUpType(vars[i], SCIP_LOCKTYPE_MODEL);

   if( nnonzstmp == 0 )
      return SCIP_OKAY;

   success = TRUE;

   /* build the matrix structure */
   SCIP_CALL( SCIPallocBuffer(scip, pmatrix) );
   matrix = *pmatrix;

   SCIP_CALL( SCIPduplicateBufferArray(scip, &matrix->colvar, vars, nvars) );

   matrix->ncols = nvars;
   matrix->nnonzssize = nnonzstmp;
   matrix->nnonzs = 0;
   matrix->nrows = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatbeg, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colmatcnt, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lb, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->ub, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colintegral, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colimplintegral, matrix->ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->colinnonlinterm, matrix->ncols) );

   /* init bounds */
   for( i = 0; i < matrix->ncols; ++i )
   {
      matrix->lb[i] = SCIPvarGetLbGlobal(vars[i]);
      matrix->ub[i] = SCIPvarGetUbGlobal(vars[i]);
      matrix->colintegral[i] = SCIPvarIsIntegral(vars[i]);
      matrix->colimplintegral[i] = SCIPvarIsImpliedIntegral(vars[i]);
      matrix->colinnonlinterm[i] = FALSE;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatval, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatind, nnonzstmp) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatbeg, nmatrixrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowmatcnt, nmatrixrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->lhs, nmatrixrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rhs, nmatrixrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &matrix->rowcons, nmatrixrows) );

   /* loop a second time over constraints handlers and add supported constraints to the matrix */
   for( i = 0; i < nconshdlrs; ++i )
   {
      SCIP_CONS** conshdlrconss;
      int c;
      int v;

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
               success = FALSE;
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
                  success = FALSE;
                  break;
               }
               weights = SCIPgetWeightsKnapsack(scip, cons);
               nrowvars = SCIPgetNVarsKnapsack(scip, cons);
               for( v = 0; v < nrowvars; ++v )
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
               success = FALSE;
               break;
            }

            switch( SCIPgetTypeSetppc(scip, cons) )
            {
            case SCIP_SETPPCTYPE_PARTITIONING:
               lhs = 1.0;
               rhs = 1.0;
               break;
            case SCIP_SETPPCTYPE_PACKING:
               lhs = -SCIPinfinity(scip);
               rhs = 1.0;
               break;
            case SCIP_SETPPCTYPE_COVERING:
               lhs = 1.0;
               rhs = SCIPinfinity(scip);
               break;
            default:
               SCIPABORT();
               return SCIP_ERROR;
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
               success = FALSE;
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
                  success = FALSE;
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
               success = FALSE;
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
               success = FALSE;
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
               success = FALSE;
               break;
            }

            SCIP_CALL( addXorLinearization(scip, matrix, cons, SCIPgetVarsXor(scip, cons), SCIPgetNVarsXor(scip, cons),
                                           SCIPgetIntVarXor(scip, cons), (SCIP_Real) SCIPgetRhsXor(scip, cons)) );
         }
      }
   }

   if( success )
   {
      SCIP_CALL( matrixSetColumnMajor(scip, matrix) );
      /**@todo scale continuous columns by least common multiple of coefficients */
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
      SCIPfreeBufferArray(scip, &matrix->colinnonlinterm);
      SCIPfreeBufferArray(scip, &matrix->colimplintegral);
      SCIPfreeBufferArray(scip, &matrix->colintegral);
      SCIPfreeBufferArray(scip, &matrix->ub);
      SCIPfreeBufferArray(scip, &matrix->lb);
      SCIPfreeBufferArray(scip, &matrix->colmatcnt);
      SCIPfreeBufferArray(scip, &matrix->colmatbeg);
      SCIPfreeBufferArray(scip, &matrix->colmatind);
      SCIPfreeBufferArray(scip, &matrix->colmatval);
      SCIPfreeBufferArrayNull(scip, &matrix->colvar);
      SCIPfreeBuffer(scip, pmatrix);
   }

   return SCIP_OKAY;
}

/** frees the matrix from memory */
static
void matrixFree(
   SCIP*                 scip,               /**< the scip data structure */
   IMPLINT_MATRIX**      pmatrix             /**< pointer to the allocated matrix */
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
      assert(matrix->colimplintegral != NULL);
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
      SCIPfreeBufferArray(scip, &(matrix->colinnonlinterm));
      SCIPfreeBufferArray(scip, &(matrix->colimplintegral));
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

      SCIPfreeBufferArrayNull(scip, &(matrix->colvar));
      SCIPfreeBuffer(scip, &matrix);
   }
}

/** creates the matrix components data structure */
static
SCIP_RETCODE createMatrixComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   IMPLINT_MATRIX *      matrix,             /**< The constraint matrix */
   MATRIX_COMPONENTS**   pmatrixcomponents   /**< Pointer to create the matrix components data structure */
   )
{
   int i;

   SCIP_CALL( SCIPallocBuffer(scip, pmatrixcomponents) );
   MATRIX_COMPONENTS* comp = *pmatrixcomponents;

   int nrows = matrixGetNRows(matrix);
   int ncols = matrixGetNCols(matrix);

   comp->nmatrixrows = nrows;
   comp->nmatrixcols = ncols;

   SCIP_CALL( SCIPallocBufferArray(scip, &comp->rowcomponent, nrows) );
   for( i = 0; i < nrows; ++i )
   {
      comp->rowcomponent[i] = -1;
   }
   SCIP_CALL( SCIPallocBufferArray(scip, &comp->colcomponent, ncols) );
   for( i = 0; i < ncols; ++i )
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

/** frees the matrix components data structure */
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

/** finds the representative of an element in the disjoint set datastructure
 *  Afterwards compresses the path to speed up subsequent queries.
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

/** merges two sets/elements into one set. Returns the index of the merged element
 *  The provided elements to be merged must be representative (i.e. returned by disjointSetFind()).
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

/** computes the connected components of the submatrix given by all continuous columns */
static
SCIP_RETCODE computeContinuousComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   IMPLINT_MATRIX*       matrix,             /**< the constraint matrix to compute the components for */
   MATRIX_COMPONENTS*    comp,               /**< the connected components data structure to store the components in */
   SCIP_Bool             includeimplints     /**< should implied integral variables be treated continuous? */
   )
{
   int* disjointset;
   int* representativecomponent;
   int* componentnextrowindex;
   int* componentnextcolindex;
   int col;
   int row;
   int i;

   /* let columns and rows share an index by mapping row index i to artificial column index i + nmatrixcols */
   SCIP_CALL( SCIPallocBufferArray(scip, &disjointset, comp->nmatrixcols + comp->nmatrixrows) );
   for( i = 0; i < comp->nmatrixcols + comp->nmatrixrows; ++i )
      disjointset[i] = -1;

   for( col = 0; col < comp->nmatrixcols; ++col )
   {
      if( matrixColIsIntegral(matrix, col) && ( !includeimplints || !matrixColIsImpliedIntegral(matrix, col) ) )
         continue;

      int* colrows = matrixGetColumnInds(matrix, col);
      int colnnonzs = matrixGetColumnNNonzs(matrix, col);
      int colrep = disjointSetFind(disjointset, col);

      for( i = 0; i < colnnonzs; ++i )
      {
         int colrow = colrows[i];
         int ind = colrow + comp->nmatrixcols;
         int rowrep = disjointSetFind(disjointset, ind);

         if( colrep != rowrep )
            colrep = disjointSetMerge(disjointset, colrep, rowrep);
      }
   }

   /* fill in the relevant data */
   SCIP_CALL( SCIPallocBufferArray(scip, &representativecomponent, comp->nmatrixcols + comp->nmatrixrows) );
   for( i = 0; i < comp->nmatrixcols + comp->nmatrixrows; ++i )
      representativecomponent[i] = -1;

   comp->ncomponents = 0;

   for( col = 0; col < comp->nmatrixcols; ++col )
   {
      if( matrixColIsIntegral(matrix,col) && ( !includeimplints || !matrixColIsImpliedIntegral(matrix, col) ) )
         continue;

      int colroot = disjointSetFind(disjointset, col);
      int component = representativecomponent[colroot];

      /* add new component */
      if( component < 0 )
      {
         assert(component == -1);
         component = comp->ncomponents;
         representativecomponent[colroot] = component;
         comp->componentcolend[component] = 0;
         comp->componentrowend[component] = 0;
         ++comp->ncomponents;
      }

      comp->colcomponent[col] = component;
      ++comp->componentcolend[component];
   }

   for( row = 0; row < comp->nmatrixrows; ++row )
   {
      int rowroot = disjointSetFind(disjointset, row + comp->nmatrixcols);
      int component = representativecomponent[rowroot];

      /* any unseen row can be skipped because it has no continuous column */
      if( component < 0 )
      {
         assert(component == -1);
         continue;
      }

      comp->rowcomponent[row] = component;
      ++comp->componentrowend[component];
   }

   if( comp->ncomponents >= 1 )
   {
      for( i = 1; i < comp->ncomponents; ++i )
      {
         comp->componentcolend[i] += comp->componentcolend[i - 1];
         comp->componentrowend[i] += comp->componentrowend[i - 1];
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &componentnextcolindex, comp->ncomponents) );
      SCIP_CALL( SCIPallocBufferArray(scip, &componentnextrowindex, comp->ncomponents) );

      componentnextcolindex[0] = 0;
      componentnextrowindex[0] = 0;

      for( i = 1; i < comp->ncomponents; ++i )
      {
         componentnextcolindex[i] = comp->componentcolend[i - 1];
         componentnextrowindex[i] = comp->componentrowend[i - 1];
      }

      for( col = 0; col < comp->nmatrixcols; ++col )
      {
         int component = comp->colcomponent[col];

         if( component < 0 )
         {
            assert(component == -1);
            continue;
         }

         comp->componentcols[componentnextcolindex[component]] = col;
         ++componentnextcolindex[component];
      }

      for( row = 0; row < comp->nmatrixrows; ++row )
      {
         int component = comp->rowcomponent[row];

         if( component < 0 )
         {
            assert(component == -1);
            continue;
         }

         comp->componentrows[componentnextrowindex[component]] = row;
         ++componentnextrowindex[component];
      }

#ifndef NDEBUG
      for( i = 0; i < comp->ncomponents; ++i )
      {
         assert(componentnextcolindex[i] == comp->componentcolend[i]);
         assert(componentnextrowindex[i] == comp->componentrowend[i]);
      }
#endif

      SCIPfreeBufferArray(scip, &componentnextrowindex);
      SCIPfreeBufferArray(scip, &componentnextcolindex);
   }

   SCIPfreeBufferArray(scip, &representativecomponent);
   SCIPfreeBufferArray(scip, &disjointset);

   return SCIP_OKAY;
}

/** creates the matrix statistics data structure */
static
SCIP_RETCODE computeMatrixStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   IMPLINT_MATRIX*       matrix,             /**< The constraint matrix to compute the statistics for */
   MATRIX_STATISTICS**   pstats,             /**< Pointer to allocate the statistics data structure at */
   SCIP_Real             numericslimit       /**< The limit beyond which we consider integrality of coefficients
                                              *   to be unreliable */
   )
{
   int i;
   int j;

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

   for( i = 0; i < nrows; ++i )
   {
      SCIP_Real lhs = matrixGetRowLhs(matrix, i);
      SCIP_Real rhs = matrixGetRowRhs(matrix, i);
      int* cols = matrixGetRowInds(matrix, i);
      SCIP_Real* vals = matrixGetRowVals(matrix, i);
      int nnonz = matrixGetRowNNonzs(matrix, i);
      stats->rownnonz[i] = nnonz;
      stats->rowequality[i] = !SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs) && SCIPisEQ(scip, lhs, rhs);

      SCIP_Bool integral = ( SCIPisInfinity(scip, -lhs) || SCIPisIntegral(scip, lhs) )
                        && ( SCIPisInfinity(scip, rhs) || SCIPisIntegral(scip, rhs) );
      SCIP_Bool badnumerics = FALSE;

      int ncontinuous = 0;
      int ncontinuouspmone = 0;
      for( j = 0; j < nnonz; ++j )
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
            /* @todo for exact version of plugin, adjust to tight check */
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

   for( i = 0; i < ncols; ++i )
   {
      SCIP_Real lb = matrixGetColLb(matrix, i);
      SCIP_Real ub = matrixGetColUb(matrix, i);
      /* @todo for exact version of plugin, adjust to tight check */
      stats->colintegralbounds[i] = ( SCIPisInfinity(scip, -lb) || SCIPisIntegral(scip, lb) )
                                 && ( SCIPisInfinity(scip, ub) || SCIPisIntegral(scip, ub) );

      /* Check that integer variables have integer bounds, as expected. */
      assert(!matrixColIsIntegral(matrix, i) || stats->colintegralbounds[i]);
   }

   return SCIP_OKAY;
}

/** frees the matrix statistics data structure */
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

/** detects components of implied integral variables
 *  Given the continuous components and statistics on the matrix, each component is checked if the associated matrix
 *  describes either a network or a transposed network (or both, in which case it is represented by a planar graph) and
 *  whether bounds/sides/coefficients are integral.
 *  We choose to check if it is a (transposed) network matrix either in a row-wise or in a column-wise fashion,
 *  depending on the size of the component. Finally, every variable that is in a network matrix or transposed network
 *  matrix is derived to be weakly implied integral.
 */
static
SCIP_RETCODE findImpliedIntegers(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PRESOLDATA*      presoldata,         /**< data belonging to the presolver */
   IMPLINT_MATRIX*       matrix,             /**< constraint matrix to compute implied integral variables for */
   MATRIX_COMPONENTS*    comp,               /**< continuous connected components of the matrix */
   MATRIX_STATISTICS*    stats,              /**< statistics of the matrix */
   int*                  nchgvartypes        /**< pointer to count the number of changed variable types */
   )
{
   assert(presoldata != NULL);

   SCIP_NETMATDEC* dec = NULL;
   SCIP_NETMATDEC* transdec = NULL;
   SCIP_Real* tempValArray;
   SCIP_Bool* compNetworkValid;
   SCIP_Bool* compTransNetworkValid;
   SCIP_Bool runintdetection = presoldata->convertintegers && SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) >= 1;
   int* tempIdxArray;
   int component;
   int col;
   int row;
   int i;
   int j;

   /* TODO: some checks to prevent expensive memory initialization if not necessary.
    * For example, make sure that there exist some +-1 candidate columns exist before performing these allocations.
    */
   SCIP_CALL( SCIPnetmatdecCreate(SCIPblkmem(scip), &dec, comp->nmatrixrows, comp->nmatrixcols) );
   SCIP_CALL( SCIPnetmatdecCreate(SCIPblkmem(scip), &transdec, comp->nmatrixcols, comp->nmatrixrows) );

   /* Because the rows may also contain non-continuous columns, we need to remove these from the array that we
    * pass to the network matrix decomposition method. We use these working arrays for this purpose.
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &tempValArray, comp->nmatrixcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tempIdxArray, comp->nmatrixcols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &compNetworkValid, comp->ncomponents) );
   SCIP_CALL( SCIPallocBufferArray(scip, &compTransNetworkValid, comp->ncomponents) );

   for( component = 0; component < comp->ncomponents; ++component )
   {
      int startrow = (component == 0) ? 0 : comp->componentrowend[component - 1];
      int nrows = comp->componentrowend[component] - startrow;
      SCIP_Bool componentokay = TRUE;

      for( i = startrow; i < startrow + nrows; ++i )
      {
         row = comp->componentrows[i];

         if( stats->rowncontinuous[row] != stats->rowncontinuouspmone[row] || !stats->rowintegral[row] || stats->rowbadnumerics[row] )
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

      int startcol = (component == 0) ? 0 : comp->componentcolend[component - 1];
      int ncols = comp->componentcolend[component] - startcol;

      for( i = startcol; i < startcol + ncols; ++i )
      {
         col = comp->componentcols[i];

         if( !stats->colintegralbounds[col] || matrixColInNonlinearTerm(matrix, col) )
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

      /* check if the component is a network matrix */
      SCIP_Bool componentnetwork = TRUE;

      /* We use the row-wise algorithm only if the number of columns is much larger than the number of rows.
       * Generally, the column-wise algorithm will be faster, but in these extreme cases, the row algorithm is faster.
       * Only very few instances should use the row-wise algorithm.
       */
      if( nrows * presoldata->columnrowratio < ncols )
      {
         for( i = startrow; i < startrow + nrows && componentnetwork; ++i )
         {
            row = comp->componentrows[i];
            SCIP_Real* rowvals = matrixGetRowVals(matrix, row);
            int* rowcols = matrixGetRowInds(matrix, row);
            int rownnonzs = matrixGetRowNNonzs(matrix, row);
            int contnnonzs = 0;

            for( j = 0; j < rownnonzs; ++j )
            {
               int rowcol = rowcols[j];

               if( !matrixColIsIntegral(matrix, rowcol) )
               {
                  tempIdxArray[contnnonzs] = rowcol;
                  tempValArray[contnnonzs] = rowvals[j];
                  ++contnnonzs;
                  assert(SCIPisEQ(scip, ABS(rowvals[j]), 1.0));
               }
            }

            SCIP_CALL( SCIPnetmatdecTryAddRow(dec, row, tempIdxArray, tempValArray, contnnonzs, &componentnetwork) );
         }
      }
      else
      {
         for( i = startcol; i < startcol + ncols && componentnetwork; ++i )
         {
            col = comp->componentcols[i];
            SCIP_Real* colvals = matrixGetColumnVals(matrix, col);
            int* colrows = matrixGetColumnInds(matrix, col);
            int colnnonzs = matrixGetColumnNNonzs(matrix, col);

            SCIP_CALL( SCIPnetmatdecTryAddCol(dec, col, colrows, colvals, colnnonzs, &componentnetwork) );
         }
      }

      if( !componentnetwork )
         SCIPnetmatdecRemoveComponent(dec, &comp->componentrows[startrow], nrows, &comp->componentcols[startcol], ncols);

      compNetworkValid[component] = componentnetwork;

      /* we don't need to check if component is both network and transposed network in case we do not want to extend
       * implied integrality to the enforced integeral variables
       */
      if( componentnetwork && !runintdetection )
      {
         compTransNetworkValid[component] = FALSE;
         continue;
      }

      SCIP_Bool componenttransnetwork = TRUE;

      /* for the transposed matrix, the situation is exactly reversed because the row/column algorithms are swapped */
      if( nrows <= ncols * presoldata->columnrowratio )
      {
         for( i = startrow; i < startrow + nrows && componenttransnetwork; ++i )
         {
            row = comp->componentrows[i];
            SCIP_Real* rowvals = matrixGetRowVals(matrix, row);
            int* rowcols = matrixGetRowInds(matrix, row);
            int rownnonzs = matrixGetRowNNonzs(matrix, row);
            int contnnonzs = 0;

            for( j = 0; j < rownnonzs; ++j )
            {
               int rowcol = rowcols[j];

               if( !matrixColIsIntegral(matrix, rowcol) )
               {
                  tempIdxArray[contnnonzs] = rowcol;
                  tempValArray[contnnonzs] = rowvals[j];
                  ++contnnonzs;
                  assert(SCIPisEQ(scip, ABS(rowvals[j]), 1.0));
               }
            }

            SCIP_CALL( SCIPnetmatdecTryAddCol(transdec, row, tempIdxArray, tempValArray, contnnonzs, &componenttransnetwork) );
         }
      }
      else
      {
         for( i = startcol; i < startcol + ncols && componenttransnetwork; ++i )
         {
            col = comp->componentcols[i];
            SCIP_Real* colvals = matrixGetColumnVals(matrix, col);
            int* colrows = matrixGetColumnInds(matrix, col);
            int colnnonzs = matrixGetColumnNNonzs(matrix, col);

            SCIP_CALL( SCIPnetmatdecTryAddRow(transdec, col, colrows, colvals, colnnonzs, &componenttransnetwork) );
         }
      }

      if( !componenttransnetwork )
         SCIPnetmatdecRemoveComponent(transdec, &comp->componentcols[startcol], ncols, &comp->componentrows[startrow], nrows);

      compTransNetworkValid[component] = componenttransnetwork;
   }

   /* add continuous columns; here we can take both normal or transposed components */
   for( component = 0; component < comp->ncomponents; ++component )
   {
      if( !compNetworkValid[component] && !compTransNetworkValid[component] )
         continue;

      int startcol = (component == 0) ? 0 : comp->componentcolend[component - 1];
      int endcol = comp->componentcolend[component];

      for( i = startcol; i < endcol; ++i )
      {
         col = comp->componentcols[i];
         assert(SCIPnetmatdecContainsColumn(dec, col) || SCIPnetmatdecContainsRow(transdec, col));
         SCIP_VAR* var = matrixGetVar(matrix, col);
         assert(!SCIPvarIsIntegral(var));
         SCIP_Bool infeasible;

         SCIP_CALL( SCIPchgVarImplType(scip, var, SCIP_IMPLINTTYPE_WEAK, &infeasible) );
         assert(!infeasible);
         ++(*nchgvartypes);
      }
   }

   /* detect implied integrality for integer columns; first, we compute valid columns that have only +-1 entries in
    * rows that are integral; then, we sort these and greedily attempt to add them to the (transposed) network matrix
    */
   if( runintdetection )
   {
      MATRIX_COMPONENTS* implintcomp;
      SCIP_Bool* implCompNetworkValid;
      SCIP_Bool* implCompTransNetworkValid;

      /**@todo avoid work when there is no implied integer variables by taking the original components instead */
      SCIP_CALL( createMatrixComponents(scip, matrix, &implintcomp) );
      SCIP_CALL( computeContinuousComponents(scip, matrix, implintcomp, TRUE) );

      SCIP_CALL( SCIPallocBufferArray(scip, &implCompNetworkValid, implintcomp->ncomponents) );
      SCIP_CALL( SCIPallocBufferArray(scip, &implCompTransNetworkValid, implintcomp->ncomponents) );

      /* extend network and transposed network decomposition by missing implied integral columns */
      for( component = 0; component < implintcomp->ncomponents; ++component )
      {
         SCIP_Bool componentnetwork = TRUE;
         SCIP_Bool componenttransnetwork = TRUE;

         int startrow = (component == 0) ? 0 : implintcomp->componentrowend[component - 1];
         int endrow = implintcomp->componentrowend[component];

         for( i = startrow; i < endrow; ++i )
         {
            row = implintcomp->componentrows[i];
            int contcomponent = comp->rowcomponent[row];

            /* integrality and numerics of rows in continuous components is already checked */
            if( contcomponent != -1 )
            {
               componentnetwork = componentnetwork && compNetworkValid[contcomponent];
               componenttransnetwork = componenttransnetwork && compTransNetworkValid[contcomponent];
            }
            else if( !stats->rowintegral[row] || stats->rowbadnumerics[row] )
            {
               componentnetwork = FALSE;
               componenttransnetwork = FALSE;
               break;
            }
         }

         if( !componentnetwork && !componenttransnetwork )
         {
            implCompNetworkValid[component] = FALSE;
            implCompTransNetworkValid[component] = FALSE;
            continue;
         }

         int startcol = (component == 0) ? 0 : implintcomp->componentcolend[component - 1];
         int endcol = implintcomp->componentcolend[component];

         for( i = startcol; i < endcol; ++i )
         {
            col = implintcomp->componentcols[i];
            int contcomponent = comp->colcomponent[col];

            /* unity and linearity of columns in continuous components is already checked */
            if( contcomponent != -1 )
            {
               componentnetwork = componentnetwork && compNetworkValid[contcomponent];
               componenttransnetwork = componenttransnetwork && compTransNetworkValid[contcomponent];
            }
            else
            {
               assert(stats->colintegralbounds[col]);

               SCIP_Real* colvals = matrixGetColumnVals(matrix, col);
               int colnnonz = matrixGetColumnNNonzs(matrix, col);
               SCIP_Bool implpmone = !matrixColInNonlinearTerm(matrix, col);

               for( j = 0; j < colnnonz && implpmone; ++j )
               {
                  if( !SCIPisEQ(scip, ABS(colvals[j]), 1.0) )
                     implpmone = FALSE;
               }

               if( !implpmone )
               {
                  componentnetwork = FALSE;
                  componenttransnetwork = FALSE;
                  break;
               }
            }
         }

         if( !componentnetwork && !componenttransnetwork )
         {
            implCompNetworkValid[component] = FALSE;
            implCompTransNetworkValid[component] = FALSE;
            continue;
         }

         /* try extending the network and transposed network by the implied integral columns of the component */
         for( i = startcol; i < endcol; ++i )
         {
            col = implintcomp->componentcols[i];
            int contcomponent = comp->colcomponent[col];

            if( contcomponent != -1 )
            {
               assert(!matrixColIsIntegral(matrix, col));
               continue;
            }
            assert(matrixColIsImpliedIntegral(matrix, col));

            SCIP_Real* colvals = matrixGetColumnVals(matrix, col);
            int* colrows = matrixGetColumnInds(matrix, col);
            int colnnonz = matrixGetColumnNNonzs(matrix, col);

            /* If a column can not be added, this does not invalidate implied integrality but means that the
             * integrality constraints of adjacent columns may be required for a differnt reason. Thus, we do not need
             * to remove components here altogether, like we did before.
             */
            if( componentnetwork )
            {
               assert(!SCIPnetmatdecContainsColumn(dec, col));
               SCIP_CALL( SCIPnetmatdecTryAddCol(dec, col, colrows, colvals, colnnonz, &componentnetwork) );
            }

            if( componenttransnetwork )
            {
               assert(!SCIPnetmatdecContainsRow(transdec, col));
               SCIP_CALL( SCIPnetmatdecTryAddRow(transdec, col, colrows, colvals, colnnonz, &componenttransnetwork) );
            }

            if( !componentnetwork && !componenttransnetwork )
               break;
         }

         implCompNetworkValid[component] = componentnetwork;
         implCompTransNetworkValid[component] = componenttransnetwork;
      }

      INTEGER_CANDIDATE_DATA* candidates;
      int numCandidates = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &candidates, comp->nmatrixcols) );

      /* candidates are non-implied integral columns with +- 1 entries without any nonzeros in bad rows */
      for( col = 0; col < comp->nmatrixcols; ++col )
      {
         if( !SCIPvarIsNonimpliedIntegral(matrixGetVar(matrix, col)) || matrixColInNonlinearTerm(matrix, col) )
            continue;
         assert(matrixColIsIntegral(matrix, col));

         SCIP_Real* colvals = matrixGetColumnVals(matrix, col);
         int* colrows = matrixGetColumnInds(matrix, col);
         int colnnonz = matrixGetColumnNNonzs(matrix, col);
         INTEGER_CANDIDATE_DATA* data = candidates + numCandidates;
         SCIP_Bool badColumn = FALSE;

         data->column = col;
         data->numContNetworkEntries = 0;
         data->numContPlanarEntries = 0;
         data->numContTransNetworkEntries = 0;

         for( i = 0; i < colnnonz; ++i )
         {
            int colrow = colrows[i];

            if( !stats->rowintegral[colrow] || stats->rowbadnumerics[colrow]
               || !SCIPisEQ(scip, ABS(colvals[i]), 1.0) )
            {
               badColumn = TRUE;
               break;
            }

            int rowcomponent = implintcomp->rowcomponent[colrow];

            if( rowcomponent != -1 )
            {
               SCIP_Bool networkValid = implCompNetworkValid[rowcomponent];
               SCIP_Bool transNetworkValid = implCompTransNetworkValid[rowcomponent];

               if( networkValid && transNetworkValid )
                  ++data->numContPlanarEntries;
               else if( networkValid )
                  ++data->numContNetworkEntries;
               else if( transNetworkValid )
                  ++data->numContTransNetworkEntries;
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

      SCIP_CALL( SCIPallocBufferArray(scip, &candidateScores, numCandidates) );

      /* higher score: pick this variable first */
      for( i = 0; i < numCandidates; ++i )
      {
         col = candidates[i].column;
         int nnonzs = matrixGetColumnNNonzs(matrix, col);

         /* @TODO test different scores / alternatives */
         /* we generally prefer to detect implied integrality of general integer variables over binary variables */
         if( SCIPvarGetType(matrixGetVar(matrix, col)) == SCIP_VARTYPE_BINARY )
            candidateScores[i] = 10.0;
         else
         {
            assert(SCIPvarGetType(matrixGetVar(matrix, col)) == SCIP_VARTYPE_INTEGER);
            candidateScores[i] = 100.0;
         }

         /* we break ties using the number of nonzeros in the column */
         candidateScores[i] -= 0.001 * nnonzs;

         /* @TODO detect when all columns only extend the network / transposed components, then we can take both */
      }

      int* indArray;
      int integerNetwork = 0;
      int integerTransNetwork = 0;

      SCIP_CALL( SCIPallocBufferArray(scip, &indArray, numCandidates) );

      for( i = 0; i < numCandidates; ++i )
         indArray[i] = i;

      SCIPsortDownRealInt(candidateScores, indArray, numCandidates);

      for( i = 0; i < numCandidates; ++i )
      {
         INTEGER_CANDIDATE_DATA* candidate = candidates + indArray[i];

         if( candidate->numContTransNetworkEntries == 0 )
         {
            col = candidate->column;
            SCIP_Real* colvals = matrixGetColumnVals(matrix, col);
            int* colrows = matrixGetColumnInds(matrix, col);
            int colnnonz = matrixGetColumnNNonzs(matrix, col);
            SCIP_Bool success;

            SCIP_CALL( SCIPnetmatdecTryAddCol(dec, col, colrows, colvals, colnnonz, &success) );

            if( success )
               ++integerNetwork;
         }

         if( candidate->numContNetworkEntries == 0 )
         {
            col = candidate->column;
            SCIP_Real* colvals = matrixGetColumnVals(matrix, col);
            int* colrows = matrixGetColumnInds(matrix, col);
            int colnnonz = matrixGetColumnNNonzs(matrix, col);
            SCIP_Bool success;

            SCIP_CALL( SCIPnetmatdecTryAddRow(transdec, col, colrows, colvals, colnnonz, &success) );

            if( success )
               ++integerTransNetwork;
         }
      }

      /* we add all enforced integral columns from the network matrix */
      if( integerNetwork >= integerTransNetwork )
      {
         for( i = 0; i < numCandidates; ++i )
         {
            col = candidates[indArray[i]].column;

            if( SCIPnetmatdecContainsColumn(dec, col) )
            {
               SCIP_VAR* var = matrixGetVar(matrix, col);
               assert(SCIPvarIsNonimpliedIntegral(var));
               SCIP_Bool infeasible;

               SCIP_CALL( SCIPchgVarImplType(scip, var, SCIP_IMPLINTTYPE_WEAK, &infeasible) );
               assert(!infeasible);
               ++(*nchgvartypes);
            }
         }
      }
      /* we add all enforced integral columns from the transposed network matrix */
      else
      {
         for( i = 0; i < numCandidates; ++i )
         {
            col = candidates[indArray[i]].column;

            if( SCIPnetmatdecContainsRow(transdec, col) )
            {
               SCIP_VAR* var = matrixGetVar(matrix, col);
               assert(SCIPvarIsNonimpliedIntegral(var));
               SCIP_Bool infeasible = FALSE;

               SCIP_CALL( SCIPchgVarImplType(scip, var, SCIP_IMPLINTTYPE_WEAK, &infeasible) );
               assert(!infeasible);
               ++(*nchgvartypes);
            }
         }
      }

      SCIPfreeBufferArray(scip, &indArray);
      SCIPfreeBufferArray(scip, &candidateScores);
      SCIPfreeBufferArray(scip, &candidates);

      SCIPfreeBufferArray(scip, &implCompTransNetworkValid);
      SCIPfreeBufferArray(scip, &implCompNetworkValid);
      freeMatrixComponents(scip, &implintcomp);
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
   SCIP_PRESOLDATA* sourcepresoldata;
   SCIP_PRESOLDATA* targetpresoldata;

   assert(scip != NULL);
   assert(presol != NULL);
   assert(strcmp(SCIPpresolGetName(presol), PRESOL_NAME) == 0);

   /* call inclusion method of presolver */
   SCIP_CALL( SCIPincludePresolImplint(scip) );

   /* copy computedimplints flag */
   sourcepresoldata = SCIPpresolGetData(presol);
   assert(sourcepresoldata != NULL);
   targetpresoldata = SCIPpresolGetData(SCIPfindPresol(scip, PRESOL_NAME));
   assert(targetpresoldata != NULL);
   targetpresoldata->computedimplints = sourcepresoldata->computedimplints;

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

   /* TODO: check these conditions again */
   /* disable implied integrality detection if we are probing or in NLP context */
   if( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPinProbing(scip) || SCIPisNLPEnabled(scip) )
      return SCIP_OKAY;

   /* skip implied integrality detection in branch-and-price, since it relies on rows being static */
   if( SCIPisStopped(scip) || SCIPgetNActivePricers(scip) > 0 || !SCIPallowStrongDualReds(scip) )
      return SCIP_OKAY;

   /* only run if we would otherwise terminate presolving */
   if( !SCIPisPresolveFinished(scip) )
      return SCIP_OKAY;

   SCIP_PRESOLDATA* presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   /* terminate if it already ran */
   if( presoldata->computedimplints )
      return SCIP_OKAY;

   presoldata->computedimplints = TRUE;

   *result = SCIP_DIDNOTFIND;

   /* exit early if there are no variables to upgrade */
   if( SCIPgetNContVars(scip) == 0
      && ( !presoldata->convertintegers || SCIPgetNBinVars(scip) + SCIPgetNIntVars(scip) == 0 ) )
      return SCIP_OKAY;

   SCIP_Real starttime = SCIPgetSolvingTime(scip);
   SCIP_Real endtime;
   IMPLINT_MATRIX* matrix = NULL;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "   (%.1fs) implied integrality detection started\n", starttime);

   SCIP_CALL( matrixCreate(scip, &matrix) );

   if( matrix == NULL )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL,
            "   (%.1fs) implied integrality detection stopped because problem contains unsuitable constraints\n",
            SCIPgetSolvingTime(scip));
      return SCIP_OKAY;
   }

   MATRIX_COMPONENTS* comp = NULL;
   MATRIX_STATISTICS* stats = NULL;
   int beforechanged = *nchgvartypes;
   int afterchanged;

   /* run implied integrality detection algorithm */
   SCIP_CALL( createMatrixComponents(scip, matrix, &comp) );
   SCIP_CALL( computeMatrixStatistics(scip, matrix, &stats, presoldata->numericslimit) );
   SCIP_CALL( computeContinuousComponents(scip, matrix, comp, FALSE) );
   SCIP_CALL( findImpliedIntegers(scip, presoldata, matrix, comp, stats, nchgvartypes) );

   afterchanged = *nchgvartypes;
   endtime = SCIPgetSolvingTime(scip);

   if( afterchanged == beforechanged )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   (%.1fs) no implied integral variables detected (time: %.2fs)\n",
            endtime, endtime - starttime);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "   (%.1fs) %d implied integral variables detected (time: %.2fs)\n",
            endtime, afterchanged - beforechanged, endtime - starttime);

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

   /* include implint presolver */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
         PRESOL_TIMING, presolExecImplint, presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyImplint) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeImplint) );

   presoldata->computedimplints = FALSE;

   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/implint/convertintegers",
         "should implied integrality also be detected for enforced integral variables?",
         &presoldata->convertintegers, FALSE, DEFAULT_CONVERTINTEGERS, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/implint/columnrowratio",
         "use the network row addition algorithm when the column to row ratio becomes larger than this threshold",
         &presoldata->columnrowratio, TRUE, DEFAULT_COLUMNROWRATIO, 0.0, 1e98, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/implint/numericslimit",
         "a row that contains variables with coefficients that are greater in absolute value than this limit is not considered for implied integrality detection",
         &presoldata->numericslimit, TRUE, DEFAULT_NUMERICSLIMIT, 1.0, 1e98, NULL, NULL) );

   return SCIP_OKAY;
}
