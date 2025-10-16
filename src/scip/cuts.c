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

/**@file   cuts.c
 * @ingroup OTHER_CFILES
 * @brief  methods for aggregation of rows
 * @author Jakob Witzig
 * @author Leona Gottwald
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "blockmemshell/memory.h"
#include "scip/cuts.h"
#include "scip/certificate.h"
#include "scip/dbldblarith.h"
#include "scip/intervalarith.h"
#include "scip/lp.h"
#include "scip/misc.h"
#include "scip/pub_lp.h"
#include "scip/pub_lpexact.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_misc_select.h"
#include "scip/pub_misc_sort.h"
#include "scip/pub_var.h"
#include "scip/scip_certificate.h"
#include "scip/scip_cut.h"
#include "scip/scip_exact.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_var.h"
#include "scip/struct_lp.h"
#include "scip/struct_lpexact.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_certificate.h"
#include "scip/type_certificate.h"
#include "scip/rational.h"

/* =========================================== general static functions =========================================== */
#ifdef SCIP_DEBUG
static
void printCutQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            cutcoefs,           /**< non-zero coefficients of cut */
   QUAD(SCIP_Real        cutrhs),            /**< right hand side of the MIR row */
   int*                  cutinds,            /**< indices of problem variables for non-zero coefficients */
   int                   cutnnz,             /**< number of non-zeros in cut */
   SCIP_Bool             ignoresol,
   SCIP_Bool             islocal
   )
{
   SCIP_Real QUAD(activity);
   SCIP_VAR** vars;

   assert(scip != NULL);
   vars = SCIPgetVars(scip);

   SCIPdebugMsg(scip, "CUT:");
   QUAD_ASSIGN(activity, 0.0);

   /**! [SnippetCodeStyleInLoopDeclaration] */
   for( int i = 0; i < cutnnz; ++i )
   {
      SCIP_Real QUAD(coef);

      QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]);

      if( SCIPvarGetType(vars[cutinds[i]]) == SCIP_VARTYPE_BINARY )
         SCIPdebugMsgPrint(scip, " %+g<%s>[B]", QUAD_TO_DBL(coef), SCIPvarGetName(vars[cutinds[i]]));
      else if( SCIPvarGetType(vars[cutinds[i]]) == SCIP_VARTYPE_INTEGER )
         SCIPdebugMsgPrint(scip, " %+g<%s>[I]", QUAD_TO_DBL(coef), SCIPvarGetName(vars[cutinds[i]]));
      else
         SCIPdebugMsgPrint(scip, " %+g<%s>[C]", QUAD_TO_DBL(coef), SCIPvarGetName(vars[cutinds[i]]));

      if( ! ignoresol )
      {
         SCIPquadprecProdQD(coef, coef, (sol == NULL ? SCIPvarGetLPSol(vars[cutinds[i]]) : SCIPgetSolVal(scip, sol, vars[cutinds[i]])));
      }
      else
      {
         if( cutcoefs[i] > 0.0 )
         {
            SCIPquadprecProdQD(coef, coef, (islocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]])));
         }
         else
         {
            SCIPquadprecProdQD(coef, coef, (islocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]])));
         }
      }

      SCIPquadprecSumQQ(activity, activity, coef);
   }
   /**! [SnippetCodeStyleInLoopDeclaration] */
   SCIPdebugMsgPrint(scip, " <= %.6f (activity: %g)\n", QUAD_TO_DBL(cutrhs), QUAD_TO_DBL(activity));
}

static
void printCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            cutcoefs,           /**< non-zero coefficients of cut */
   SCIP_Real             cutrhs,             /**< right hand side of the MIR row */
   int*                  cutinds,            /**< indices of problem variables for non-zero coefficients */
   int                   cutnnz,             /**< number of non-zeros in cut */
   SCIP_Bool             ignoresol,
   SCIP_Bool             islocal
   )
{
   SCIP_Real activity;
   SCIP_VAR** vars;
   int i;

   assert(scip != NULL);
   vars = SCIPgetVars(scip);

   SCIPdebugMsg(scip, "CUT:");
   activity = 0.0;
   for( i = 0; i < cutnnz; ++i )
   {
      SCIP_Real coef;

      coef = cutcoefs[cutinds[i]];

      if( SCIPvarGetType(vars[cutinds[i]]) == SCIP_VARTYPE_BINARY )
         SCIPdebugMsgPrint(scip, " %+g<%s>[B]", coef, SCIPvarGetName(vars[cutinds[i]]));
      else if( SCIPvarGetType(vars[cutinds[i]]) == SCIP_VARTYPE_INTEGER )
         SCIPdebugMsgPrint(scip, " %+g<%s>[I]", coef, SCIPvarGetName(vars[cutinds[i]]));
      else
         SCIPdebugMsgPrint(scip, " %+g<%s>[C]", coef, SCIPvarGetName(vars[cutinds[i]]));

      if( ! ignoresol )
      {
         coef = coef * (sol == NULL ? SCIPvarGetLPSol(vars[cutinds[i]]) : SCIPgetSolVal(scip, sol, vars[cutinds[i]]));
      }
      else
      {
         if( cutcoefs[i] > 0.0 )
         {
            coef = coef * (islocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]));
         }
         else
         {
            coef = coef * (islocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]));
         }
      }

      activity += coef;
   }
   SCIPdebugMsgPrint(scip, " <= %.6f (activity: %g)\n", cutrhs, activity);
}
#endif

/** macro to make sure a value is not equal to zero, i.e. NONZERO(x) != 0.0
 *  will be TRUE for every x including 0.0
 *
 *  To avoid branches it will add 1e-100 with the same sign as x to x which will
 *  be rounded away for any sane non-zero value but will make sure the value is
 *  never exactly 0.0.
 */
#define NONZERO(x)   (COPYSIGN(1e-100, (x)) + (x))

/** add a scaled row to a dense vector indexed over the problem variables and keep the
 *  index of non-zeros up-to-date
 */
static
SCIP_RETCODE varVecAddScaledRowCoefs(
   int*RESTRICT          inds,               /**< pointer to array with variable problem indices of non-zeros in variable vector */
   SCIP_Real*RESTRICT    vals,               /**< array with values of variable vector */
   int*RESTRICT          nnz,                /**< number of non-zeros coefficients of variable vector */
   SCIP_ROW*             row,                /**< row coefficients to add to variable vector */
   SCIP_Real             scale               /**< scale for adding given row to variable vector */
   )
{
   int i;

   assert(inds != NULL);
   assert(vals != NULL);
   assert(nnz != NULL);
   assert(row != NULL);

   /* add the non-zeros to the aggregation row and keep non-zero index up to date */
   for( i = 0 ; i < row->len; ++i )
   {
      SCIP_Real val;
      int probindex;

      probindex = row->cols[i]->var_probindex;
      val = vals[probindex];

      if( val == 0.0 )
         inds[(*nnz)++] = probindex;

      val += row->vals[i] * scale;

      /* the value must not be exactly zero due to sparsity pattern */
      val = NONZERO(val);

      assert(val != 0.0);
      vals[probindex] = val;
   }

   return SCIP_OKAY;
}

/** add a scaled row to a dense vector indexed over the problem variables and keep the
 *  index of non-zeros up-to-date
 *
 *  This is the quad precision version of varVecAddScaledRowCoefs().
 */
static
SCIP_RETCODE varVecAddScaledRowCoefsQuad(
   int*RESTRICT          inds,               /**< pointer to array with variable problem indices of non-zeros in variable vector */
   SCIP_Real*RESTRICT    vals,               /**< array with values of variable vector */
   int*RESTRICT          nnz,                /**< number of non-zeros coefficients of variable vector */
   SCIP_ROW*             row,                /**< row coefficients to add to variable vector */
   SCIP_Real             scale               /**< scale for adding given row to variable vector */
   )
{
   int i;

   assert(inds != NULL);
   assert(vals != NULL);
   assert(nnz != NULL);
   assert(row != NULL);

   /* add the non-zeros to the aggregation row and keep non-zero index up to date */
   for( i = 0 ; i < row->len; ++i )
   {
      SCIP_Real QUAD(scaledrowval);
      SCIP_Real QUAD(val);
      int probindex;

      probindex = row->cols[i]->var_probindex;
      QUAD_ARRAY_LOAD(val, vals, probindex);

      if( QUAD_HI(val) == 0.0 )
         inds[(*nnz)++] = probindex;

      SCIPquadprecProdDD(scaledrowval, row->vals[i], scale);
      SCIPquadprecSumQQ(val, val, scaledrowval);

      /* the value must not be exactly zero due to sparsity pattern */
      QUAD_HI(val) = NONZERO(QUAD_HI(val));
      assert(QUAD_HI(val) != 0.0);

      QUAD_ARRAY_STORE(vals, probindex, val);
   }

   return SCIP_OKAY;
}

/** add a scaled row to a dense vector indexed over the problem variables and keep the
 *  index of non-zeros up-to-date
 *
 *  This is the quad precision version of varVecAddScaledRowCoefs() with a quad precision scaling factor.
 */
static
SCIP_RETCODE varVecAddScaledRowCoefsQuadScale(
   int*RESTRICT          inds,               /**< pointer to array with variable problem indices of non-zeros in variable vector */
   SCIP_Real*RESTRICT    vals,               /**< array with values of variable vector */
   int*RESTRICT          nnz,                /**< number of non-zeros coefficients of variable vector */
   SCIP_ROW*             row,                /**< row coefficients to add to variable vector */
   QUAD(SCIP_Real        scale)              /**< scale for adding given row to variable vector */
   )
{
   int i;

   assert(inds != NULL);
   assert(vals != NULL);
   assert(nnz != NULL);
   assert(row != NULL);

   /* add the non-zeros to the aggregation row and keep non-zero index up to date */
   for( i = 0 ; i < row->len; ++i )
   {
      SCIP_Real QUAD(val);
      SCIP_Real QUAD(rowval);
      int probindex;

      probindex = row->cols[i]->var_probindex;
      QUAD_ARRAY_LOAD(val, vals, probindex);

      if( QUAD_HI(val) == 0.0 )
      {
         inds[(*nnz)++] = probindex;
         SCIPquadprecProdQD(val, scale, row->vals[i]);
      }
      else
      {
         SCIPquadprecProdQD(rowval, scale, row->vals[i]);
         SCIPquadprecSumQQ(val, val, rowval);
      }

      /* the value must not be exactly zero due to sparsity pattern */
      QUAD_HI(val) = NONZERO(QUAD_HI(val));
      assert(QUAD_HI(val) != 0.0);

      QUAD_ARRAY_STORE(vals, probindex, val);
   }

   return SCIP_OKAY;
}

/** add a scaled row to a dense vector indexed over the problem variables and keep the index of non-zeros up-to-date
 *
 *  In the safe variant, we need to transform all variables (implicitly) to nonnegative variables using their
 *  upper/lower bounds.  When adding \f$\lambda * (c^Tx \le d)\f$ to \f$a^Tx \le b\f$, this results in:
 *
 *    \f{align*}{
 *      m_i & =a_i+\lambda c_i \\
 *      U \cap L & = \emptyset \\
 *      U & = \{ i : x_i \le u_i\} \\
 *      L & = \{ i : x_i \ge l_i\} \\
 *      \sum_{i \in U} \overline{m_i}x_i + \sum_{i \in L}\underline{m_i}x_i
 *                    & \le b+ \lambda d + \sum_{i \in U, u_i > 0}(\overline{m_i}-\underline{m_i})u_i + \sum_{i \in L, l_i < 0}(\underline{m_i}-\overline{m_i})l_i
 *    \f}
 *
 *  This methods sums up the left hand side, and stores the change of the rhs due to the variable bounds in rhschange.
 *
 *  @note this method is safe for usage in exact solving mode
 */
static
SCIP_RETCODE varVecAddScaledRowCoefsSafely(
   SCIP*                 scip,               /**< scip data structure */
   int*                  inds,               /**< pointer to array with variable problem indices of non-zeros in variable vector */
   SCIP_Real*            vals,               /**< array with values of variable vector */
   int*                  nnz,                /**< number of non-zeros coefficients of variable vector */
   SCIP_ROW*             row,                /**< row coefficients to add to variable vector */
   SCIP_Real             scale,              /**< scale for adding given row to variable vector */
   SCIP_Real*            rhschange,          /**< change in rhs due to variable conjugation */
   SCIP_Bool*            success             /**< was the addition successful? */
   )
{
   int i;
   SCIP_ROUNDMODE previousroundmode;
   SCIP_VAR* var;
   SCIP_ROWEXACT* rowexact;

   assert(SCIPisExact(scip));
   assert(inds != NULL);
   assert(vals != NULL);
   assert(nnz != NULL);
   assert(row != NULL);
   assert(rhschange != NULL);
   assert(success != NULL);
   assert(*success);

   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();

   *rhschange = 0;
   rowexact = SCIProwGetRowExact(row);

   /* add the non-zeros to the aggregation row and keep non-zero index up to date */
   for( i = 0 ; i < row->len; ++i )
   {
      SCIP_Real val;
      SCIP_INTERVAL valinterval;
      int probindex;

      probindex = row->cols[i]->var_probindex;
      var = row->cols[i]->var;
      val = vals[probindex];

      if( val == 0.0 )
      {
         assert(*nnz < SCIPgetNVars(scip));
         inds[(*nnz)++] = probindex;
      }

      if( val == SCIP_INVALID ) /*lint !e777*/
         val = 0.0;

      SCIPintervalSetBounds(&valinterval, rowexact->valsinterval[i].inf, rowexact->valsinterval[i].sup);
      SCIPintervalMulScalar(SCIPinfinity(scip), &valinterval, valinterval, scale);
      SCIPintervalAddScalar(SCIPinfinity(scip), &valinterval, valinterval, val);

      if( SCIPisInfinity(scip, REALABS(valinterval.inf)) || SCIPisInfinity(scip, REALABS(valinterval.sup)) )
      {
         *success = FALSE;
         SCIPintervalSetRoundingMode(previousroundmode);
         return SCIP_OKAY;
      }

      if( SCIPvarGetLbGlobal(var) > -SCIPinfinity(scip) && SCIPvarGetLbGlobal(var) >= 0 )
         val = valinterval.inf; 
      else if(SCIPvarGetUbGlobal(var) < SCIPinfinity(scip) && SCIPvarGetUbGlobal(var) <= 0 )
         val = valinterval.sup;
      else if( SCIPvarGetLbGlobal(var) > -SCIPinfinity(scip) )
      {
         val = valinterval.inf;
         SCIPintervalSetRoundingModeUpwards();
         *rhschange += (valinterval.sup - valinterval.inf) * (-SCIPvarGetLbGlobal(var));
      }
      else if( SCIPvarGetUbGlobal(var) < SCIPinfinity(scip) )
      {
         val = valinterval.sup;
         SCIPintervalSetRoundingModeUpwards();
         *rhschange += (valinterval.sup - valinterval.inf) * (SCIPvarGetUbGlobal(var));
      }
      else
      {
         *success = FALSE;
         SCIPintervalSetRoundingMode(previousroundmode);
         return SCIP_OKAY;
      }

      /* we can't set the value to 0 or the sparsity pattern does not work. We can't perturb it slightly because we are solving
       * exactly; this is taken care of in removeZerosSafely */
      if( val == 0.0 )
         val = SCIP_INVALID;

      vals[probindex] = val;
   }

   SCIPintervalSetRoundingMode(previousroundmode);

   return SCIP_OKAY;
}

/** calculates the cut efficacy for the given solution */
static
SCIP_Real calcEfficacy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to calculate the efficacy for (NULL for LP solution) */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   SCIP_Real             cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int                   cutnnz              /**< the number of non-zeros in the cut */
   )
{
   SCIP_VAR** vars;
   SCIP_Real norm = 0.0;
   SCIP_Real activity = 0.0;
   int i;

   assert(scip != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);

   vars = SCIPgetVars(scip);

   switch( scip->set->sepa_efficacynorm )
   {
   case 'e':
      for( i = 0; i < cutnnz; ++i )
      {
         activity += cutcoefs[i] * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         norm += SQR(cutcoefs[i]);
      }
      norm = sqrt(norm);
      break;
   case 'm':
      for( i = 0; i < cutnnz; ++i )
      {
         SCIP_Real absval;

         activity += cutcoefs[i] * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         absval = REALABS(cutcoefs[i]);
         norm = MAX(norm, absval);
      }
      break;
   case 's':
      for( i = 0; i < cutnnz; ++i )
      {
         activity += cutcoefs[i] * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         norm += REALABS(cutcoefs[i]);
      }
      break;
   case 'd':
      for( i = 0; i < cutnnz; ++i )
      {
         activity += cutcoefs[i] * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         if( !SCIPisZero(scip, cutcoefs[i]) )
            norm = 1.0;
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c'\n", scip->set->sepa_efficacynorm);
      assert(FALSE); /*lint !e506*/
   }

   return (activity - cutrhs) / MAX(1e-6, norm);
}

/** calculates the efficacy norm of the given aggregation row, which depends on the "separating/efficacynorm" parameter */
static
SCIP_Real calcEfficacyNormQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            vals,               /**< array of the non-zero coefficients in the vector; this is a quad precision array! */
   int*                  inds,               /**< array of the problem indices of variables with a non-zero coefficient in the vector */
   int                   nnz                 /**< the number of non-zeros in the vector */
   )
{
   SCIP_Real norm = 0.0;
   SCIP_Real QUAD(coef);
   int i;

   assert(scip != NULL);
   assert(scip->set != NULL);

   switch( scip->set->sepa_efficacynorm )
   {
   case 'e':
      for( i = 0; i < nnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, vals, inds[i]);
         norm += SQR(QUAD_TO_DBL(coef));
      }
      norm = sqrt(norm);
      break;
   case 'm':
      for( i = 0; i < nnz; ++i )
      {
         SCIP_Real absval;
         QUAD_ARRAY_LOAD(coef, vals, inds[i]);

         absval = REALABS(QUAD_TO_DBL(coef));
         norm = MAX(norm, absval);
      }
      break;
   case 's':
      for( i = 0; i < nnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, vals, inds[i]);
         norm += REALABS(QUAD_TO_DBL(coef));
      }
      break;
   case 'd':
      for( i = 0; i < nnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, vals, inds[i]);
         if( !SCIPisZero(scip, QUAD_TO_DBL(coef)) )
         {
            norm = 1.0;
            break;
         }
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c.'\n", scip->set->sepa_efficacynorm);
      assert(FALSE); /*lint !e506*/
   }

   return norm;
}

/** calculates the cut efficacy for the given solution; the cut coefs are stored densely */
static
SCIP_Real calcEfficacyDenseStorage(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to calculate the efficacy for (NULL for LP solution) */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut; this is a quad precision array! */
   SCIP_Real             cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int                   cutnnz              /**< the number of non-zeros in the cut */
   )
{
   SCIP_VAR** vars;
   SCIP_Real norm = 0.0;
   SCIP_Real activity = 0.0;
   SCIP_Real coef;
   int i;

   assert(scip != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
   assert(scip->set != NULL);

   vars = SCIPgetVars(scip);

   switch( scip->set->sepa_efficacynorm )
   {
   case 'e':
      for( i = 0; i < cutnnz; ++i )
      {
         coef = cutcoefs[cutinds[i]];
         activity += coef * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         norm += SQR(coef);
      }
      norm = SQR(norm);
      break;
   case 'm':
      for( i = 0; i < cutnnz; ++i )
      {
         SCIP_Real absval;

         coef = cutcoefs[cutinds[i]];
         activity += coef * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         absval = REALABS(coef);
         norm = MAX(norm, absval);
      }
      break;
   case 's':
      for( i = 0; i < cutnnz; ++i )
      {
         coef = cutcoefs[cutinds[i]];
         activity += coef * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         norm += REALABS(coef);
      }
      break;
   case 'd':
      for( i = 0; i < cutnnz; ++i )
      {
         coef = cutcoefs[cutinds[i]];
         activity += coef * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         if( !SCIPisZero(scip, coef) )
            norm = 1.0;
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c.'\n", scip->set->sepa_efficacynorm);
      assert(FALSE); /*lint !e506*/
   }

   return (activity - cutrhs) / MAX(1e-6, norm);
}

/** calculates the cut efficacy for the given solution; the cut coefs are stored densely and in quad precision */
static
SCIP_Real calcEfficacyDenseStorageQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to calculate the efficacy for (NULL for LP solution) */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut; this is a quad precision array! */
   SCIP_Real             cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int                   cutnnz              /**< the number of non-zeros in the cut */
   )
{
   SCIP_VAR** vars;
   SCIP_Real norm = 0.0;
   SCIP_Real activity = 0.0;
   SCIP_Real QUAD(coef);
   int i;

   assert(scip != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
   assert(scip->set != NULL);

   vars = SCIPgetVars(scip);

   switch( scip->set->sepa_efficacynorm )
   {
   case 'e':
      for( i = 0; i < cutnnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]);
         activity += QUAD_TO_DBL(coef) * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         norm += SQR(QUAD_TO_DBL(coef));
      }
      norm = sqrt(norm);
      break;
   case 'm':
      for( i = 0; i < cutnnz; ++i )
      {
         SCIP_Real absval;

         QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]);
         activity += QUAD_TO_DBL(coef) * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         absval = REALABS(QUAD_TO_DBL(coef));
         norm = MAX(norm, absval);
      }
      break;
   case 's':
      for( i = 0; i < cutnnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]);
         activity += QUAD_TO_DBL(coef) * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         norm += REALABS(QUAD_TO_DBL(coef));
      }
      break;
   case 'd':
      for( i = 0; i < cutnnz; ++i )
      {
         QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]);
         activity += QUAD_TO_DBL(coef) * SCIPgetSolVal(scip, sol, vars[cutinds[i]]);
         if( !SCIPisZero(scip, QUAD_TO_DBL(coef)) )
            norm = 1.0;
      }
      break;
   default:
      SCIPerrorMessage("invalid efficacy norm parameter '%c.'\n", scip->set->sepa_efficacynorm);
      assert(FALSE); /*lint !e506*/
   }

   return (activity - cutrhs) / MAX(1e-6, norm);
}

/** safely (in the exact solving mode sense) remove all items with |a_i| or |u_i - l_i)| below the given value
 *
 *  Returns TRUE if the cut became redundant.
 *  If it is a local cut, use local bounds, otherwise, use global bounds.
 *
 * * @note this method is safe for usage in exact solving mode
 */
static
SCIP_Bool removeZerosSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             minval,             /**< minimal absolute value of coefficients that should not be removed */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz              /**< the number of non-zeros in the cut */
   )
{
   int i;
   SCIP_VAR** vars;
   SCIP_ROUNDMODE previousroundmode;

   assert(SCIPisExact(scip));

   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();

   vars = SCIPgetVars(scip);

   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real val;
      SCIP_Real lb;
      SCIP_Real ub;
      int v;
      SCIP_Bool isfixed;

      v = cutinds[i];
      val = cutcoefs[v];

      if( val == SCIP_INVALID ) /*lint !e777*/
         val = 0.0;

      /* for now we always use global bounds in exact solving mode (could be improved for local cuts in the future) */
      lb = SCIPvarGetLbGlobal(vars[v]);
      ub = SCIPvarGetUbGlobal(vars[v]);

      if( !(SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub)) && SCIPisEQ(scip, ub, lb) )
         isfixed = TRUE;
      else
         isfixed = FALSE;

      if( EPSZ(val, minval) || isfixed )
      {
         /* adjust left and right hand sides with max contribution */
         if( val < 0.0 )
         {
            if( SCIPisInfinity(scip, ub) )
            {
               SCIPintervalSetRoundingMode(previousroundmode);
               return TRUE;
            }
            else
               *cutrhs += (-val) * ub;
         }
         else
         {
            if( SCIPisInfinity(scip, -lb) )
            {
               SCIPintervalSetRoundingMode(previousroundmode);
               return TRUE;
            }
            else
               *cutrhs += (-val) * lb;
         }

         val = 0.0;
         cutcoefs[v] = val;

         /* remove non-zero entry */
         --(*cutnnz);
         cutinds[i] = cutinds[*cutnnz];
      }
      else
         ++i;
   }

   /* relax rhs to 0, if it's very close to 0 */
   if( *cutrhs < 0.0 && *cutrhs >= -SCIPepsilon(scip) )
      *cutrhs = 0.0;

   SCIPintervalSetRoundingMode(previousroundmode);

   return FALSE;
}

/** safely remove all items with |a_i| or |u_i - l_i)| below the given value
 *
 *  Returns TRUE if the cut became redundant.
 *  If it is a local cut, use local bounds, otherwise, use global bounds.
 *
 *  @note this method is safe for usage in exact solving mode
 */
static
SCIP_Bool removeZerosQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             minval,             /**< minimal absolute value of coefficients that should not be removed */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   QUAD(SCIP_Real*       cutrhs),            /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz              /**< the number of non-zeros in the cut */
   )
{
   int i;
   SCIP_VAR** vars;

   vars = SCIPgetVars(scip);

   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real QUAD(val);
      SCIP_Real lb;
      SCIP_Real ub;
      int v;
      SCIP_Bool isfixed;

      v = cutinds[i];
      QUAD_ARRAY_LOAD(val, cutcoefs, v);

      if( cutislocal )
      {
         lb = SCIPvarGetLbLocal(vars[v]);
         ub = SCIPvarGetUbLocal(vars[v]);
      }
      else
      {
         lb = SCIPvarGetLbGlobal(vars[v]);
         ub = SCIPvarGetUbGlobal(vars[v]);
      }

      if( !(SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub)) && SCIPisEQ(scip, ub, lb) )
         isfixed = TRUE;
      else
         isfixed = FALSE;

      if( isfixed || EPSZ(QUAD_TO_DBL(val), minval) )
      {
         if( REALABS(QUAD_TO_DBL(val)) > QUAD_EPSILON )
         {
            /* adjust right hand side with max contribution */
            if( QUAD_TO_DBL(val) < 0.0 )
            {
               if( SCIPisInfinity(scip, ub) )
                  return TRUE;
               else
               {
                  SCIPquadprecProdQD(val, val, ub);
                  SCIPquadprecSumQQ(*cutrhs, *cutrhs, -val);
               }
            }
            else
            {
               if( SCIPisInfinity(scip, -lb) )
                  return TRUE;
               else
               {
                  SCIPquadprecProdQD(val, val, lb);
                  SCIPquadprecSumQQ(*cutrhs, *cutrhs, -val);
               }
            }
         }

         QUAD_ASSIGN(val, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, val);

         /* remove non-zero entry */
         --(*cutnnz);
         cutinds[i] = cutinds[*cutnnz];
      }
      else
         ++i;
   }

   /* relax rhs to 0, if it's very close to 0 */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   return FALSE;
}

/** safely remove all items with |a_i| or |u_i - l_i| below the given value
 *
 *  Returns TRUE if the cut became redundant.
 *  If it is a local cut, use local bounds, otherwise, use global bounds.
 */
static
SCIP_Bool removeZeros(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             minval,             /**< minimal absolute value of coefficients that should not be removed */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   QUAD(SCIP_Real*       cutrhs),            /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz              /**< the number of non-zeros in the cut */
   )
{
   int i;
   SCIP_VAR** vars;

   vars = SCIPgetVars(scip);

   /* loop over non-zeros and remove values below minval; values above QUAD_EPSILON are cancelled with their bound
    * to avoid numerical rounding errors
    */
   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real val;
      SCIP_Real lb;
      SCIP_Real ub;
      int v;
      SCIP_Bool isfixed;
      SCIP_Real QUAD(quadprod);

      v = cutinds[i];
      val = cutcoefs[v];

      if( cutislocal )
      {
         lb = SCIPvarGetLbLocal(vars[v]);
         ub = SCIPvarGetUbLocal(vars[v]);
      }
      else
      {
         lb = SCIPvarGetLbGlobal(vars[v]);
         ub = SCIPvarGetUbGlobal(vars[v]);
      }

      if( !(SCIPisInfinity(scip, -lb) || SCIPisInfinity(scip, ub)) && SCIPisEQ(scip, ub, lb) )
         isfixed = TRUE;
      else
         isfixed = FALSE;

      if( EPSZ(val, minval) || isfixed )
      {
         if( REALABS(val) > QUAD_EPSILON )
         {
            /* adjust left and right hand sides with max contribution */
            if( val < 0.0 )
            {
               if( SCIPisInfinity(scip, ub) )
                  return TRUE;
               else
               {
                  SCIPquadprecProdDD(quadprod, -val, ub);
                  SCIPquadprecSumQQ(*cutrhs, *cutrhs, quadprod);
               }
            }
            else
            {
               if( SCIPisInfinity(scip, -lb) )
                  return TRUE;
               else
               {
                  SCIPquadprecProdDD(quadprod, -val, lb);
                  SCIPquadprecSumQQ(*cutrhs, *cutrhs, quadprod);
               }
            }
         }

         cutcoefs[v] = 0.0;

         /* remove non-zero entry */
         --(*cutnnz);
         cutinds[i] = cutinds[*cutnnz];
      }
      else
         ++i;
   }

   /* relax rhs to 0, if it's very close to 0 */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   return FALSE;
}

/** compare absolute values of coefficients in quad precision */
static
SCIP_DECL_SORTINDCOMP(compareAbsCoefsQuad)
{
   SCIP_Real abscoef1;
   SCIP_Real abscoef2;
   SCIP_Real QUAD(coef1);
   SCIP_Real QUAD(coef2);
   SCIP_Real* coefs = (SCIP_Real*) dataptr;

   QUAD_ARRAY_LOAD(coef1, coefs, ind1);
   QUAD_ARRAY_LOAD(coef2, coefs, ind2);

   abscoef1 = REALABS(QUAD_TO_DBL(coef1));
   abscoef2 = REALABS(QUAD_TO_DBL(coef2));

   if( abscoef1 < abscoef2 )
      return -1;
   if( abscoef2 < abscoef1 )
      return 1;

   return 0;
}

/** compare absolute values of coefficients */
static
SCIP_DECL_SORTINDCOMP(compareAbsCoefs)
{
   SCIP_Real abscoef1;
   SCIP_Real abscoef2;
   SCIP_Real* coefs = (SCIP_Real*) dataptr;

   abscoef1 = REALABS(coefs[ind1]);
   abscoef2 = REALABS(coefs[ind2]);

   if( abscoef1 < abscoef2 )
      return -1;
   if( abscoef2 < abscoef1 )
      return 1;

   return 0;
}

/** change given coefficient to new given value, adjust right hand side using the variables bound;
 *  returns TRUE if the right hand side would need to be changed to infinity and FALSE otherwise
 */
static
SCIP_Bool chgCoeffWithBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable the coefficient belongs to */
   SCIP_Real             oldcoeff,           /**< old coefficient value */
   SCIP_Real             newcoeff,           /**< new coefficient value */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   QUAD(SCIP_Real*       cutrhs)             /**< pointer to adjust right hand side of cut */
   )
{
   SCIP_Real QUAD(delta);

   SCIPquadprecSumDD(delta, newcoeff, -oldcoeff);

   if( QUAD_TO_DBL(delta) > QUAD_EPSILON )
   {
      SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);

      if( SCIPisInfinity(scip, ub) )
         return TRUE;
      else
      {
         SCIPquadprecProdQD(delta, delta, ub);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, delta);
      }
   }
   else if( QUAD_TO_DBL(delta) < -QUAD_EPSILON )
   {
      SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);

      if( SCIPisInfinity(scip, -lb) )
         return TRUE;
      else
      {
         SCIPquadprecProdQD(delta, delta, lb);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, delta);
      }
   }

   return FALSE;
}

/** change given (quad) coefficient to new given value, adjust right hand side using the variables bound;
 *  returns TRUE if the right hand side would need to be changed to infinity and FALSE otherwise
 */
static
SCIP_Bool chgQuadCoeffWithBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable the coefficient belongs to */
   QUAD(SCIP_Real        oldcoeff),          /**< old coefficient value */
   SCIP_Real             newcoeff,           /**< new coefficient value */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   QUAD(SCIP_Real*       cutrhs)             /**< pointer to adjust right hand side of cut */
   )
{
   SCIP_Real QUAD(delta);

   SCIPquadprecSumQD(delta, -oldcoeff, newcoeff);

   if( QUAD_TO_DBL(delta) > QUAD_EPSILON )
   {
      SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);

      if( SCIPisInfinity(scip, ub) )
         return TRUE;
      else
      {
         SCIPquadprecProdQD(delta, delta, ub);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, delta);
      }
   }
   else if( QUAD_TO_DBL(delta) < -QUAD_EPSILON )
   {
      SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);

      if( SCIPisInfinity(scip, -lb) )
         return TRUE;
      else
      {
         SCIPquadprecProdQD(delta, delta, lb);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, delta);
      }
   }

   return FALSE;
}

/** change given coefficient to new given value, adjust right hand side using the variables bound;
 *  returns TRUE if the right hand side would need to be changed to infinity and FALSE otherwise
 */
static
SCIP_Bool chgCoeffWithBoundSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable the coefficient belongs to */
   SCIP_Real             oldcoeff,           /**< old coefficient value */
   SCIP_Real             newcoeff,           /**< new coefficient value */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutrhs              /**< pointer to adjust right hand side of cut */
   )
{
   SCIP_INTERVAL delta;
   SCIP_ROUNDMODE previousroundmode;

   assert(SCIPisExact(scip));

   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();

   SCIPintervalSet(&delta, newcoeff);
   SCIPintervalSubScalar(SCIPinfinity(scip), &delta, delta, oldcoeff);

   if( SCIPintervalGetSup(delta) > 0 )
   {
      SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);

      if( SCIPisInfinity(scip, ub) )
      {
         SCIPintervalSetRoundingMode(previousroundmode);
         return TRUE;
      }
      else
      {
         SCIPintervalMulScalar(SCIPinfinity(scip), &delta, delta, ub);
         *cutrhs += SCIPintervalGetSup(delta);
      }
   }
   else if( SCIPintervalGetInf(delta) < 0 )
   {
      SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);

      if( SCIPisInfinity(scip, -lb) )
      {
         SCIPintervalSetRoundingMode(previousroundmode);
         return TRUE;
      }
      else
      {
         SCIPintervalMulScalar(SCIPinfinity(scip), &delta, delta, lb);
         *cutrhs += SCIPintervalGetSup(delta);
      }
   }
   else
   {
      return TRUE;
   }

   return FALSE;
}


/** scales the cut and then tightens the coefficients of the given cut based on the maximal activity;
 *  see cons_linear.c consdataTightenCoefs() for details; the cut is given in a semi-sparse quad precision array;
 *
 *  This is the quad precision version of cutTightenCoefs() below.
 */
static
SCIP_RETCODE cutTightenCoefsQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   QUAD(SCIP_Real*       cutrhs),            /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< the number of non-zeros in the cut */
   SCIP_Bool*            redundant           /**< whether the cut was detected to be redundant */
   )
{
   int i;
   int nintegralvars;
   SCIP_Bool isintegral = TRUE;
   SCIP_VAR** vars;
   SCIP_Real QUAD(maxacttmp);
   SCIP_Real maxact;
   SCIP_Real maxabsintval = 0.0;
   SCIP_Real maxabscontval = 0.0;

   QUAD_ASSIGN(maxacttmp, 0.0);

   vars = SCIPgetVars(scip);
   nintegralvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

   assert(redundant != NULL);
   *redundant = FALSE;

   /* compute maximal activity and maximal absolute coefficient values for all and for integral variables in the cut */
   for( i = 0; i < *cutnnz; ++i )
   {
      SCIP_Real QUAD(val);

      assert(cutinds[i] >= 0);
      assert(vars[cutinds[i]] != NULL);

      QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);

      if( QUAD_TO_DBL(val) < 0.0 )
      {
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, -lb) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, -QUAD_TO_DBL(val));
         else
         {
            maxabscontval = MAX(maxabscontval, -QUAD_TO_DBL(val));
            isintegral = FALSE;
         }

         SCIPquadprecProdQD(val, val, lb);
         SCIPquadprecSumQQ(maxacttmp, maxacttmp, val);
      }
      else
      {
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, ub) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, QUAD_TO_DBL(val));
         else
         {
            maxabscontval = MAX(maxabscontval, QUAD_TO_DBL(val));
            isintegral = FALSE;
         }

         SCIPquadprecProdQD(val, val, ub);
         SCIPquadprecSumQQ(maxacttmp, maxacttmp, val);
      }
   }

   maxact = QUAD_TO_DBL(maxacttmp);

   /* cut is redundant in activity bounds */
   if( SCIPisFeasLE(scip, maxact, QUAD_TO_DBL(*cutrhs)) )
   {
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   /* cut is only on integral variables, try to scale to integral coefficients */
   if( isintegral )
   {
      SCIP_Real equiscale;
      SCIP_Real intscalar;
      SCIP_Bool success;
      SCIP_Real* intcoeffs;

      SCIP_CALL( SCIPallocBufferArray(scip, &intcoeffs, *cutnnz) );

      equiscale = 1.0 / MIN((maxact - QUAD_TO_DBL(*cutrhs)), maxabsintval);

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(val);

         QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);
         SCIPquadprecProdQD(val, val, equiscale);

         intcoeffs[i] = QUAD_TO_DBL(val);
      }

      SCIP_CALL( SCIPcalcIntegralScalar(intcoeffs, *cutnnz, -SCIPsumepsilon(scip), SCIPsumepsilon(scip),
            (SCIP_Longint)scip->set->sepa_maxcoefratio, scip->set->sepa_maxcoefratio, &intscalar, &success) );

      SCIPfreeBufferArray(scip, &intcoeffs);

      if( success )
      {
         /* if successful, apply the scaling */
         intscalar *= equiscale;
         SCIPquadprecProdQD(*cutrhs, *cutrhs, intscalar);

         for( i = 0; i < *cutnnz; )
         {
            SCIP_Real QUAD(val);
            SCIP_Real intval;

            QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);
            SCIPquadprecProdQD(val, val, intscalar);

            intval = SCIPround(scip, QUAD_TO_DBL(val));

            if( chgQuadCoeffWithBound(scip, vars[cutinds[i]], QUAD(val), intval, cutislocal, QUAD(cutrhs)) )
            {
               /* TODO maybe change the coefficient to the other value instead of discarding the cut? */
               *redundant = TRUE;
               return SCIP_OKAY;
            }

            if( intval != 0.0 )
            {
               QUAD_ASSIGN(val, intval);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], val);
               ++i;
            }
            else
            {
               /* this must not be -0.0, otherwise the clean buffer memory is not cleared properly */
               QUAD_ASSIGN(val, 0.0);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], val);
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
            }
         }

         SCIPquadprecEpsFloorQ(*cutrhs, *cutrhs, SCIPfeastol(scip)); /*lint !e666*/

         /* recompute the maximal activity after scaling to integral values */
         QUAD_ASSIGN(maxacttmp, 0.0);
         maxabsintval = 0.0;

         for( i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real QUAD(val);

            assert(cutinds[i] >= 0);
            assert(vars[cutinds[i]] != NULL);

            QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);

            if( QUAD_TO_DBL(val) < 0.0 )
            {
               SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, -QUAD_TO_DBL(val));

               SCIPquadprecProdQD(val, val, lb);

               SCIPquadprecSumQQ(maxacttmp, maxacttmp, val);
            }
            else
            {
               SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, QUAD_TO_DBL(val));

               SCIPquadprecProdQD(val, val, ub);

               SCIPquadprecSumQQ(maxacttmp, maxacttmp, val);
            }
         }

         assert(EPSISINT(QUAD_TO_DBL(maxacttmp), 1e-4));
         SCIPquadprecSumQD(maxacttmp, maxacttmp, 0.5);
         SCIPquadprecFloorQ(maxacttmp, maxacttmp);
      }
      else
      {
         /* otherwise, apply the equilibrium scaling */
         isintegral = FALSE;

         /* perform the scaling */
         SCIPquadprecProdQD(maxacttmp, maxacttmp, equiscale);
         SCIPquadprecProdQD(*cutrhs, *cutrhs, equiscale);
         maxabsintval *= equiscale;

         for( i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real QUAD(val);

            QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);
            SCIPquadprecProdQD(val, val, equiscale);
            QUAD_ARRAY_STORE(cutcoefs, cutinds[i], val);
         }
      }
   }
   else
   {
      /* cut has integer and continuous variables, so scale it to equilibrium */
      SCIP_Real scale;
      SCIP_Real maxabsval;

      maxabsval = maxact - QUAD_TO_DBL(*cutrhs);
      maxabsval = MIN(maxabsval, maxabsintval);
      maxabsval = MAX(maxabsval, maxabscontval);

      scale = 1.0 / maxabsval; /*lint !e795*/

      /* perform the scaling */
      SCIPquadprecProdQD(maxacttmp, maxacttmp, scale);
      SCIPquadprecProdQD(*cutrhs, *cutrhs, scale);
      maxabsintval *= scale;

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(val);

         QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);
         SCIPquadprecProdQD(val, val, scale);
         QUAD_ARRAY_STORE(cutcoefs, cutinds[i], val);
      }
   }

   maxact = QUAD_TO_DBL(maxacttmp);

   /* check again for redundancy after scaling */
   if( SCIPisFeasLE(scip, maxact, QUAD_TO_DBL(*cutrhs)) )
   {
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   /* no coefficient tightening can be performed since the precondition doesn't hold for any of the variables */
   if( SCIPisGT(scip, maxact - maxabsintval, QUAD_TO_DBL(*cutrhs)) )
      return SCIP_OKAY;

   /* first sort indices, so that in the following sort, the order for coefficients with same absolute value does not depend on how cutinds was initially ordered */
   SCIPsortInt(cutinds, *cutnnz);
   SCIPsortDownInd(cutinds, compareAbsCoefsQuad, (void*) cutcoefs, *cutnnz);

   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real QUAD(val);

      if( cutinds[i] >= nintegralvars )
      {
         ++i;
         continue;
      }

      QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);

      assert(SCIPvarIsIntegral(vars[cutinds[i]]));

      if( QUAD_TO_DBL(val) < 0.0 && SCIPisLE(scip, maxact + QUAD_TO_DBL(val), QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, *cutrhs, -maxacttmp);

         if( isintegral )
         {
            /* if cut is integral, the true coefficient must also be integral; thus round it to exact integral value */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) > QUAD_TO_DBL(val) )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQQ(delta, -val, coef);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);
            SCIPdebugMsg(scip, "tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
               QUAD_TO_DBL(val), QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp), lb,
               cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]));

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(!SCIPisPositive(scip, QUAD_TO_DBL(coef)));

            if( SCIPisNegative(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
            }
            else
            {
               QUAD_ASSIGN(coef, 0.0);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else if( QUAD_TO_DBL(val) > 0.0 && SCIPisLE(scip, maxact - QUAD_TO_DBL(val), QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, maxacttmp, -*cutrhs);

         if( isintegral )
         {
            /* if cut is integral, the true coefficient must also be integral; thus round it to exact integral value */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) < QUAD_TO_DBL(val) )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQQ(delta, -val, coef);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);
            SCIPdebugMsg(scip, "tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
               QUAD_TO_DBL(val), QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp),
               cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]), ub);

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(SCIPisGE(scip, QUAD_TO_DBL(coef), 0.0));

            if( SCIPisPositive(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
            }
            else
            {
               QUAD_ASSIGN(coef, 0.0);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else /* due to sorting we can stop completely if the precondition was not fulfilled for this variable */
         break;

      ++i;
   }

   return SCIP_OKAY;
}

/** multiplies a parameter for a variable in a row safely (using variable bounds and increasing the rhs)
 *
 *  @return the scaled value
 */
static
SCIP_Real scaleValSafely(
   SCIP*                 scip,               /**< SCIP structure */
   SCIP_Real             val,                /**< the value that should be scaled */
   SCIP_Real             scale,              /**< scaling factor */
   SCIP_Bool             cutislocal,         /**< should local or global bounds be used */
   SCIP_VAR*             var,                /**< the variable that is relevant */
   SCIP_Real*            rhschange,          /**< resulting change in rhs of row */
   SCIP_Bool*            success             /**< was the operation successful? (false if no bounds) */
   )
{
   SCIP_ROUNDMODE previousroundmode;
   SCIP_INTERVAL valinterval;
   SCIP_Real ub;
   SCIP_Real lb;
   SCIP_Real newval = 0.0;

   assert(SCIPisExact(scip));

   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();

   *rhschange = 0;

   SCIPintervalSet(&valinterval, val);
   SCIPintervalMulScalar(SCIPinfinity(scip), &valinterval, valinterval, scale);

   lb = cutislocal ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
   ub = cutislocal ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);

   *success = TRUE;

   if( lb > -SCIPinfinity(scip) && lb >= 0 )
   {
      SCIPdebugMessage("Lb positive, no change in rhs needed \n");
      newval = SCIPintervalGetInf(valinterval);
   }
   else if(ub < SCIPinfinity(scip) && ub <= 0 )
   {
      SCIPdebugMessage("Ub negative, no change in rhs needed \n");
      newval = SCIPintervalGetSup(valinterval);
   }
   else if( lb > -SCIPinfinity(scip) )
   {
      newval = SCIPintervalGetInf(valinterval);
      SCIPintervalSetRoundingModeUpwards();
      *rhschange = (SCIPintervalGetSup(valinterval) - SCIPintervalGetInf(valinterval)) * (-lb);
      SCIPdebugMessage("Using lb %.17g corrected by %.17g. Change to rhs: %.17g \n", SCIPvarGetLbGlobal(var), -SCIPintervalGetSup(valinterval) + SCIPintervalGetInf(valinterval), *rhschange);
   }
   else if( ub < SCIPinfinity(scip) )
   {
      newval = SCIPintervalGetSup(valinterval);
      SCIPintervalSetRoundingModeUpwards();
      *rhschange = (SCIPintervalGetSup(valinterval) - SCIPintervalGetInf(valinterval)) * (ub);
      SCIPdebugMessage("Using ub %.17g corrected by %.17g. Change to rhs: %.17g \n", SCIPvarGetUbGlobal(var), SCIPintervalGetSup(valinterval) - SCIPintervalGetInf(valinterval), *rhschange);
   }
   else
   {
      *success = FALSE;
      SCIPintervalSetRoundingMode(previousroundmode);
      return newval;
   }

   SCIPintervalSetRoundingMode(previousroundmode);

   return newval;
}

/** scales the cut and then tightens the coefficients of the given cut based on the maximal activity;
 *  see cons_linear.c consdataTightenCoefs() for details;
 *
 *  This is the safe version of cutTightenCoefs() below.
 */
static
SCIP_RETCODE cutTightenCoefsSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< the number of non-zeros in the cut */
   SCIP_Bool*            redundant           /**< whether the cut was detected to be redundant */
   )
{ /*lint --e{644}*/
   int i;
   int nintegralvars;
   SCIP_Bool isintegral = TRUE;
   SCIP_VAR** vars;
   SCIP_INTERVAL maxact;
   SCIP_INTERVAL tmp;
   SCIP_Real maxabsintval = 0.0;
   SCIP_Real maxabscontval = 0.0;
   SCIP_ROUNDMODE previousroundmode;
   SCIP_MIRINFO* mirinfo = NULL;

   assert(SCIPisExact(scip));

   if( SCIPisCertified(scip)   )
   {
      mirinfo = SCIPgetCertificate(scip)->mirinfo[SCIPgetCertificate(scip)->nmirinfos - 1];
      assert(mirinfo != NULL);
   }

   SCIPintervalSet(&maxact, 0.0);

   vars = SCIPgetVars(scip);
   nintegralvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

   assert(redundant != NULL);
   *redundant = FALSE;

   /* compute maximal activity and maximal absolute coefficient values for all and for integral variables in the cut */
   for( i = 0; i < *cutnnz; ++i )
   {
      SCIP_Real val;

      assert(cutinds[i] >= 0);
      assert(vars[cutinds[i]] != NULL);

      val = cutcoefs[cutinds[i]];

      if( val < 0.0 )
      {
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, -lb) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, -val);
         else
         {
            maxabscontval = MAX(maxabscontval, -val);
            isintegral = FALSE;
         }

         SCIPintervalSet(&tmp, val);
         SCIPintervalMulScalar(SCIPinfinity(scip), &tmp, tmp, lb);
         SCIPintervalAdd(SCIPinfinity(scip), &maxact, maxact, tmp);
      }
      else
      {
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, ub) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, val);
         else
         {
            maxabscontval = MAX(maxabscontval, val);
            isintegral = FALSE;
         }

         SCIPintervalSet(&tmp, val);
         SCIPintervalMulScalar(SCIPinfinity(scip), &tmp, tmp, ub);
         SCIPintervalAdd(SCIPinfinity(scip), &maxact, maxact, tmp);
      }
   }

   /* cut is redundant in activity bounds */
   if( SCIPintervalGetSup(maxact) <= *cutrhs )
   {
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   previousroundmode = SCIPintervalGetRoundingMode();

   /* cut is only on integral variables, try to scale to integral coefficients */
   if( isintegral )
   {
      SCIP_Real equiscale;
      SCIP_Real intscalar;
      SCIP_Bool success;
      SCIP_Real* intcoeffs;
      SCIP_Real rhschange;

      SCIP_CALL( SCIPallocBufferArray(scip, &intcoeffs, *cutnnz) );

      equiscale = 1.0 / MAX(maxabscontval, maxabsintval);

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real val;

         val = cutcoefs[cutinds[i]];
         val *= equiscale;

         intcoeffs[i] = val;
      }

      SCIP_CALL( SCIPcalcIntegralScalar(intcoeffs, *cutnnz, -SCIPsumepsilon(scip), SCIPepsilon(scip),
            (SCIP_Longint)scip->set->sepa_maxcoefratio, scip->set->sepa_maxcoefratio, &intscalar, &success) );

      SCIPfreeBufferArray(scip, &intcoeffs);

      if( success )
      {
         /* if successful, apply the scaling */
         SCIPintervalSetRoundingModeUpwards();
         intscalar *= equiscale;

         *cutrhs *= intscalar;

         if( SCIPisCertified(scip) )
         {
            assert(mirinfo != NULL);
            mirinfo->scale = intscalar;
         }

         for( i = 0; i < *cutnnz; )
         {
            SCIP_Real val;
            SCIP_Real intval;

            val = cutcoefs[cutinds[i]];
            val = scaleValSafely(scip, val, intscalar, cutislocal, vars[cutinds[i]], &rhschange, &success);

            *cutrhs += rhschange;

            if( !success )
            {
               /* TODO maybe change the coefficient to the other value instead of discarding the cut? */
               *redundant = TRUE;
               SCIPintervalSetRoundingMode(previousroundmode);
               return SCIP_OKAY;
            }
            intval = SCIPround(scip, val);

            if( chgCoeffWithBoundSafely(scip, vars[cutinds[i]], val, intval, cutislocal, cutrhs) )
            {
               /* TODO maybe change the coefficient to the other value instead of discarding the cut? */
               SCIPintervalSetRoundingMode(previousroundmode);
               *redundant = TRUE;
               return SCIP_OKAY;
            }

            if( intval != 0.0 )
            {
               val = intval;
               cutcoefs[cutinds[i]] = val;
               ++i;
            }
            else
            {
               /* this must not be -0.0, otherwise the clean buffer memory is not cleared properly */
               val = 0.0;
               cutcoefs[cutinds[i]] = val;
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
            }
         }

         if( SCIPisCertified(scip) )
         {
            assert(mirinfo != NULL);
            mirinfo->unroundedrhs = *cutrhs;
         }
         *cutrhs = floor(*cutrhs); /*lint !e835*/

         /* recompute the maximal activity after scaling to integral values */
         SCIPintervalSet(&maxact, 0.0);
         maxabsintval = 0.0;

         for( i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real val;

            assert(cutinds[i] >= 0);
            assert(vars[cutinds[i]] != NULL);

            val = cutcoefs[cutinds[i]];

            if( val < 0.0 )
            {
               SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, -val);

               SCIPintervalSet(&tmp, val);
               SCIPintervalMulScalar(SCIPinfinity(scip), &tmp, tmp, lb);
               SCIPintervalAdd(SCIPinfinity(scip), &maxact, maxact, tmp);
            }
            else
            {
               SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, val);

               SCIPintervalSet(&tmp, val);
               SCIPintervalMulScalar(SCIPinfinity(scip), &tmp, tmp, ub);
               SCIPintervalAdd(SCIPinfinity(scip), &maxact, maxact, tmp);
            }
         }

         assert(EPSISINT(SCIPintervalGetSup(maxact), 1e-4));/*lint !e666*/
         /* check again for redundancy */
         if( SCIPisFeasLE(scip, SCIPintervalGetSup(maxact), *cutrhs) )
         {
            *redundant = TRUE;
            SCIPintervalSetRoundingMode(previousroundmode);
            return SCIP_OKAY;
         }
      }
      else
      {
         /* otherwise, apply the equilibrium scaling */

         /* perform the scaling */
         SCIPintervalMulScalar(SCIPinfinity(scip), &maxact, maxact, equiscale);
         SCIPintervalSetRoundingModeUpwards();

         *cutrhs *= equiscale;
         maxabsintval *= equiscale;
         if( SCIPisCertified(scip) )
         {
            assert(mirinfo != NULL);
            mirinfo->scale = equiscale;
         }

         for( i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real val;

            val = cutcoefs[cutinds[i]];
            val = scaleValSafely(scip, val, equiscale, cutislocal, vars[cutinds[i]], &rhschange, &success);

            if( !success )
            {
               *redundant = TRUE;
               SCIPintervalSetRoundingMode(previousroundmode);
               return SCIP_OKAY;
            }
            cutcoefs[cutinds[i]] = val;
            *cutrhs += rhschange;
         }
      }
   }
   else
   {
      /* cut has integer and continuous variables, so scale it to equilibrium */
      SCIP_Real scale;
      SCIP_Real maxabsval;
      SCIP_Bool success;
      SCIP_Real rhschange;

      maxabsval = SCIPintervalGetSup(maxact) - *cutrhs;
      maxabsval = MIN(maxabsval, maxabsintval);
      maxabsval = MAX(maxabsval, maxabscontval);

      scale = 1.0 / maxabsval; /*lint !e795*/

      /* perform the scaling */
      SCIPintervalSet(&maxact, scale);
      SCIPintervalSetRoundingModeUpwards();
      *cutrhs *= scale;
      maxabsintval *= scale;
      if( SCIPisCertified(scip) )
      {
         assert(mirinfo != NULL);
         mirinfo->scale = scale;
      }

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real val;

         val = cutcoefs[cutinds[i]];
         val = scaleValSafely(scip, val, scale, cutislocal, vars[cutinds[i]], &rhschange, &success);

         *cutrhs += rhschange;

         if( !success )
         {
            *redundant = TRUE;
            SCIPintervalSetRoundingMode(previousroundmode);
            return SCIP_OKAY;
         }
         cutcoefs[cutinds[i]] = val;
      }
   }

   /* no coefficient tightening can be performed since the precondition doesn't hold for any of the variables */
   if( SCIPisGT(scip, SCIPintervalGetSup(maxact) - maxabsintval, *cutrhs) )
   {
      SCIPintervalSetRoundingMode(previousroundmode);
      return SCIP_OKAY;
   }

   SCIPsortDownInd(cutinds, compareAbsCoefs, (void*) cutcoefs, *cutnnz);

#ifdef SCIP_DISABLED_CODE
   /** @todo implement and certify coefficient tightening for cuts in exact solving mode */
   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *cutnnz && FALSE; )
   {
      SCIP_Real QUAD(val);

      if( cutinds[i] >= nintegralvars )
      {
         ++i;
         continue;
      }

      QUAD_ARRAY_LOAD(val, cutcoefs, cutinds[i]);

      assert(SCIPvarIsIntegral(vars[cutinds[i]]));

      if( QUAD_TO_DBL(val) < 0.0 && SCIPisLE(scip, maxact + QUAD_TO_DBL(val), QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, *cutrhs, -maxacttmp);

         if( isintegral )
         {
            /* if cut is integral, the true coefficient must also be integral; thus round it to exact integral value */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) > QUAD_TO_DBL(val) )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQQ(delta, -val, coef);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
               QUAD_TO_DBL(val), QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp), lb,
               cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]));

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(!SCIPisPositive(scip, QUAD_TO_DBL(coef)));

            if( SCIPisNegative(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
            }
            else
            {
               QUAD_ASSIGN(coef, 0.0);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else if( QUAD_TO_DBL(val) > 0.0 && SCIPisLE(scip, maxact - QUAD_TO_DBL(val), QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, maxacttmp, -*cutrhs);

         if( isintegral )
         {
            /* if cut is integral, the true coefficient must also be integral; thus round it to exact integral value */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) < QUAD_TO_DBL(val) )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQQ(delta, -val, coef);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);
            SCIPdebugPrintf("tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
               QUAD_TO_DBL(val), QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp),
               cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]), ub);

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(SCIPisGE(scip, QUAD_TO_DBL(coef), 0.0));

            if( SCIPisPositive(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
            }
            else
            {
               QUAD_ASSIGN(coef, 0.0);
               QUAD_ARRAY_STORE(cutcoefs, cutinds[i], coef);
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else /* due to sorting we can stop completely if the precondition was not fulfilled for this variable */
         break;

      ++i;
   }
#endif /*lint --e{438}*/

   return SCIP_OKAY;
}


/** scales the cut and then tightens the coefficients of the given cut based on the maximal activity;
 *  see cons_linear.c consdataTightenCoefs() for details; the cut is given in a semi-sparse array;
 */
static
SCIP_RETCODE cutTightenCoefs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   QUAD(SCIP_Real*       cutrhs),            /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< the number of non-zeros in the cut */
   SCIP_Bool*            redundant           /**< pointer to return whtether the cut was detected to be redundant */
   )
{
   int i;
   int nintegralvars;
   SCIP_Bool isintegral = TRUE;
   SCIP_VAR** vars;
   SCIP_Real QUAD(maxacttmp);
   SCIP_Real maxact;
   SCIP_Real maxabsintval = 0.0;
   SCIP_Real maxabscontval = 0.0;

   assert(!SCIPisExact(scip));

   QUAD_ASSIGN(maxacttmp, 0.0);

   vars = SCIPgetVars(scip);
   nintegralvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);

   assert(redundant != NULL);
   *redundant = FALSE;

   /* compute maximal activity and maximal absolute coefficient values for all and for integral variables in the cut */
   for( i = 0; i < *cutnnz; ++i )
   {
      SCIP_Real val;
      SCIP_Real QUAD(quadprod);

      assert(cutinds[i] >= 0);
      assert(vars[cutinds[i]] != NULL);

      val = cutcoefs[cutinds[i]];

      if( val < 0.0 )
      {
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, -lb) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, -val);
         else
         {
            maxabscontval = MAX(maxabscontval, -val);
            isintegral = FALSE;
         }

         SCIPquadprecProdDD(quadprod, val, lb);
         SCIPquadprecSumQQ(maxacttmp, maxacttmp, quadprod);
      }
      else
      {
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, ub) )
            return SCIP_OKAY;

         if( cutinds[i] < nintegralvars )
            maxabsintval = MAX(maxabsintval, val);
         else
         {
            maxabscontval = MAX(maxabscontval, val);
            isintegral = FALSE;
         }

         SCIPquadprecProdDD(quadprod, val, ub);
         SCIPquadprecSumQQ(maxacttmp, maxacttmp, quadprod);
      }
   }

   maxact = QUAD_TO_DBL(maxacttmp);

   /* cut is redundant in activity bounds */
   if( SCIPisFeasLE(scip, maxact, QUAD_TO_DBL(*cutrhs)) )
   {
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   /* cut is only on integral variables, try to scale to integral coefficients */
   if( isintegral )
   {
      SCIP_Real equiscale;
      SCIP_Real intscalar;
      SCIP_Bool success;
      SCIP_Real* intcoeffs;

      SCIP_CALL( SCIPallocBufferArray(scip, &intcoeffs, *cutnnz) );

      equiscale = 1.0 / MIN((maxact - QUAD_TO_DBL(*cutrhs)), maxabsintval);

      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real val;

         val = equiscale * cutcoefs[cutinds[i]];

         intcoeffs[i] = val;
      }

      SCIP_CALL( SCIPcalcIntegralScalar(intcoeffs, *cutnnz, -SCIPsumepsilon(scip), SCIPsumepsilon(scip),
            (SCIP_Longint)scip->set->sepa_maxcoefratio, scip->set->sepa_maxcoefratio, &intscalar, &success) );

      SCIPfreeBufferArray(scip, &intcoeffs);

      if( success )
      {
         /* if successful, apply the scaling */
         intscalar *= equiscale;
         SCIPquadprecProdQD(*cutrhs, *cutrhs, intscalar);

         for( i = 0; i < *cutnnz; )
         {
            SCIP_Real val;
            SCIP_Real intval;

            val = cutcoefs[cutinds[i]];
            val *= intscalar;

            intval = SCIPround(scip, val);

            if( chgCoeffWithBound(scip, vars[cutinds[i]], val, intval, cutislocal, QUAD(cutrhs)) )
            {
               /* TODO maybe change the coefficient to the other value instead of discarding the cut? */
               *redundant = TRUE;
               return SCIP_OKAY;
            }

            if( intval != 0.0 )
            {
               cutcoefs[cutinds[i]] = intval;
               ++i;
            }
            else
            {
               /* this must not be -0.0, otherwise the clean buffer memory is not cleared properly */
               cutcoefs[cutinds[i]] = 0.0;
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
            }
         }

         SCIPquadprecEpsFloorQ(*cutrhs, *cutrhs, SCIPfeastol(scip)); /*lint !e666*/

         /* recompute the maximal activity after scaling to integral values */
         QUAD_ASSIGN(maxacttmp, 0.0);
         maxabsintval = 0.0;

         for( i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real val;
            SCIP_Real QUAD(quadprod);

            assert(cutinds[i] >= 0);
            assert(vars[cutinds[i]] != NULL);

            val = cutcoefs[cutinds[i]];

            if( val < 0.0 )
            {
               SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, -val);

               SCIPquadprecProdDD(quadprod, val, lb);
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, quadprod);
            }
            else
            {
               SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

               maxabsintval = MAX(maxabsintval, val);

               SCIPquadprecProdDD(quadprod, val, ub);
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, quadprod);
            }
         }

         assert(EPSISINT(QUAD_TO_DBL(maxacttmp), 1e-4));
         SCIPquadprecSumQD(maxacttmp, maxacttmp, 0.5);
         SCIPquadprecFloorQ(maxacttmp, maxacttmp);
      }
      else
      {
         /* otherwise, apply the equilibrium scaling */
         isintegral = FALSE;

         /* perform the scaling */
         SCIPquadprecProdQD(maxacttmp, maxacttmp, equiscale);
         SCIPquadprecProdQD(*cutrhs, *cutrhs, equiscale);
         maxabsintval *= equiscale;

         for( i = 0; i < *cutnnz; ++i )
            cutcoefs[cutinds[i]] *= equiscale;
      }
   }
   else
   {
      /* cut has integer and continuous variables, so scale it to equilibrium */
      SCIP_Real scale;
      SCIP_Real maxabsval;

      maxabsval = maxact - QUAD_TO_DBL(*cutrhs);
      maxabsval = MIN(maxabsval, maxabsintval);
      maxabsval = MAX(maxabsval, maxabscontval);

      scale = 1.0 / maxabsval; /*lint !e795*/

      /* perform the scaling */
      SCIPquadprecProdQD(maxacttmp, maxacttmp, scale);
      SCIPquadprecProdQD(*cutrhs, *cutrhs, scale);
      maxabsintval *= scale;

      for( i = 0; i < *cutnnz; ++i )
         cutcoefs[cutinds[i]] *= scale;
   }

   maxact = QUAD_TO_DBL(maxacttmp);

   /* check again for redundancy after scaling */
   if( SCIPisFeasLE(scip, maxact, QUAD_TO_DBL(*cutrhs)) )
   {
      *redundant = TRUE;
      return SCIP_OKAY;
   }

   /* no coefficient tightening can be performed since the precondition doesn't hold for any of the variables */
   if( SCIPisGT(scip, maxact - maxabsintval, QUAD_TO_DBL(*cutrhs)) )
      return SCIP_OKAY;

   /* first sort indices, so that in the following sort, the order for coefficients with same absolute value does not depend on how cutinds was initially ordered */
   SCIPsortInt(cutinds, *cutnnz);
   SCIPsortDownInd(cutinds, compareAbsCoefs, (void*) cutcoefs, *cutnnz);

   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *cutnnz; )
   {
      SCIP_Real val;

      if( cutinds[i] >= nintegralvars )
      {
         ++i;
         continue;
      }

      val = cutcoefs[cutinds[i]];

      assert(SCIPvarIsIntegral(vars[cutinds[i]]));

      if( val < 0.0 && SCIPisLE(scip, maxact + val, QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, -maxacttmp, *cutrhs);

         if( isintegral )
         {
            /* if cut is integral, the true coefficient must also be integral; thus round it to exact integral value */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) > val )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQD(delta, coef, -val);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);
            SCIPdebugMsg(scip, "tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
               val, QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp), lb,
               cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]));

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(!SCIPisPositive(scip, QUAD_TO_DBL(coef)));

            if( SCIPisNegative(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               cutcoefs[cutinds[i]] = QUAD_TO_DBL(coef);
            }
            else
            {
               cutcoefs[cutinds[i]] = 0.0;
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else if( val > 0.0 && SCIPisLE(scip, maxact - val, QUAD_TO_DBL(*cutrhs)) )
      {
         SCIP_Real QUAD(coef);
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         SCIPquadprecSumQQ(coef, maxacttmp, -*cutrhs);

         if( isintegral )
         {
            /* if cut is integral, the true coefficient must also be integral; thus round it to exact integral value */
            assert(SCIPisFeasIntegral(scip, QUAD_TO_DBL(coef)));
            QUAD_ASSIGN(coef, SCIPround(scip, QUAD_TO_DBL(coef)));
         }

         if( QUAD_TO_DBL(coef) < val )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumQD(delta, coef, -val);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQQ(tmp, delta, *cutrhs);
            SCIPdebugMsg(scip, "tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
               val, QUAD_TO_DBL(coef), QUAD_TO_DBL(*cutrhs), QUAD_TO_DBL(tmp),
               cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]), ub);

            QUAD_ASSIGN_Q(*cutrhs, tmp);

            assert(! SCIPisNegative(scip, QUAD_TO_DBL(coef)));

            if( SCIPisPositive(scip, QUAD_TO_DBL(coef)) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               cutcoefs[cutinds[i]] = QUAD_TO_DBL(coef);
            }
            else
            {
               cutcoefs[cutinds[i]] = 0.0;
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               continue;
            }
         }
      }
      else /* due to sorting we can stop completely if the precondition was not fulfilled for this variable */
         break;

      ++i;
   }

   return SCIP_OKAY;
}

/** perform activity based coefficient tightening on the given cut; returns TRUE if the cut was detected
 *  to be redundant due to activity bounds
 *
 *  See also cons_linear.c:consdataTightenCoefs().
 */
SCIP_Bool SCIPcutsTightenCoefficients(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut local? */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< the right hand side of the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< the number of non-zeros in the cut */
   int*                  nchgcoefs           /**< number of changed coefficients */
   )
{
   int i;
   int nintegralvars;
   SCIP_VAR** vars;
   SCIP_Real* absvals;
   SCIP_Real QUAD(maxacttmp);
   SCIP_Real maxact;
   SCIP_Real maxabsval = 0.0;
   SCIP_Bool redundant = FALSE;

   assert(nchgcoefs != NULL);

   QUAD_ASSIGN(maxacttmp, 0.0);

   vars = SCIPgetVars(scip);
   nintegralvars = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &absvals, *cutnnz) );

   assert(nchgcoefs != NULL);
   *nchgcoefs = 0;

   for( i = 0; i < *cutnnz; ++i )
   {
      SCIP_Real QUAD(quadprod);

      assert(cutinds[i] >= 0);
      assert(vars[cutinds[i]] != NULL);

      if( cutcoefs[i] < 0.0 )
      {
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, -lb) )
            goto TERMINATE;

         if( cutinds[i] < nintegralvars )
         {
            maxabsval = MAX(maxabsval, -cutcoefs[i]);
            absvals[i] = -cutcoefs[i];
         }
         else
            absvals[i] = 0.0;

         SCIPquadprecProdDD(quadprod, lb, cutcoefs[i]);
         SCIPquadprecSumQQ(maxacttmp, maxacttmp, quadprod);
      }
      else
      {
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         if( SCIPisInfinity(scip, ub) )
            goto TERMINATE;

         if( cutinds[i] < nintegralvars )
         {
            maxabsval = MAX(maxabsval, cutcoefs[i]);
            absvals[i] = cutcoefs[i];
         }
         else
            absvals[i] = 0.0;

         SCIPquadprecProdDD(quadprod, ub, cutcoefs[i]);
         SCIPquadprecSumQQ(maxacttmp, maxacttmp, quadprod);
      }
   }

   maxact = QUAD_TO_DBL(maxacttmp);

   /* cut is redundant in activity bounds */
   if( SCIPisFeasLE(scip, maxact, *cutrhs) )
   {
      redundant = TRUE;
      goto TERMINATE;
   }

   /* terminate, because coefficient tightening cannot be performed; also excludes the case in which no integral variable is present */
   if( SCIPisGT(scip, maxact - maxabsval, *cutrhs) )
      goto TERMINATE;

   SCIPsortDownRealRealInt(absvals, cutcoefs, cutinds, *cutnnz);
   SCIPfreeBufferArray(scip, &absvals);

   /* loop over the integral variables and try to tighten the coefficients; see cons_linear for more details */
   for( i = 0; i < *cutnnz; ++i )
   {
      /* due to sorting, we can exit if we reached a continuous variable: all further integral variables have 0 coefficents anyway */
      if( cutinds[i] >= nintegralvars )
         break;

      assert(SCIPvarIsIntegral(vars[cutinds[i]]));

      if( cutcoefs[i] < 0.0 && SCIPisLE(scip, maxact + cutcoefs[i], *cutrhs) )
      {
         SCIP_Real coef = (*cutrhs) - maxact;
         SCIP_Real lb = cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]);

         coef = SCIPfloor(scip, coef);

         if( coef > cutcoefs[i] )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumDD(delta, coef, -cutcoefs[i]);
            SCIPquadprecProdQD(delta, delta, lb);

            SCIPquadprecSumQD(tmp, delta, *cutrhs);
            SCIPdebugMsg(scip, "tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
               cutcoefs[i], coef, (*cutrhs), QUAD_TO_DBL(tmp), lb,
               cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]));

            *cutrhs = QUAD_TO_DBL(tmp);

            assert(!SCIPisPositive(scip, coef));

            ++(*nchgcoefs);

            if( SCIPisNegative(scip, coef) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               cutcoefs[i] = coef;
            }
            else
            {
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               cutcoefs[i] = cutcoefs[*cutnnz];
               continue;
            }
         }
      }
      else if( cutcoefs[i] > 0.0 && SCIPisLE(scip, maxact - cutcoefs[i], *cutrhs) )
      {
         SCIP_Real coef = maxact - (*cutrhs);
         SCIP_Real ub = cutislocal ? SCIPvarGetUbLocal(vars[cutinds[i]]) : SCIPvarGetUbGlobal(vars[cutinds[i]]);

         coef = SCIPceil(scip, coef);

         if( coef < cutcoefs[i] )
         {
            SCIP_Real QUAD(delta);
            SCIP_Real QUAD(tmp);

            SCIPquadprecSumDD(delta, coef, -cutcoefs[i]);
            SCIPquadprecProdQD(delta, delta, ub);

            SCIPquadprecSumQD(tmp, delta, *cutrhs);
            SCIPdebugMsg(scip, "tightened coefficient from %g to %g; rhs changed from %g to %g; the bounds are [%g,%g]\n",
               cutcoefs[i], coef, (*cutrhs), QUAD_TO_DBL(tmp),
               cutislocal ? SCIPvarGetLbLocal(vars[cutinds[i]]) : SCIPvarGetLbGlobal(vars[cutinds[i]]), ub);

            *cutrhs = QUAD_TO_DBL(tmp);

            assert(!SCIPisNegative(scip, coef));

            ++(*nchgcoefs);

            if( SCIPisPositive(scip, coef) )
            {
               SCIPquadprecSumQQ(maxacttmp, maxacttmp, delta);
               maxact = QUAD_TO_DBL(maxacttmp);
               cutcoefs[i] = coef;
            }
            else
            {
               --(*cutnnz);
               cutinds[i] = cutinds[*cutnnz];
               cutcoefs[i] = cutcoefs[*cutnnz];
               continue;
            }
         }
      }
      else /* due to sorting we can stop completely if the precondition was not fulfilled for this variable */
         break;
   }

  TERMINATE:
   SCIPfreeBufferArrayNull(scip, &absvals);

   return redundant;
}

/* =========================================== aggregation row =========================================== */


/** create an empty aggregation row
 *
 *  @note By default, this data structure uses quad precision via double-double arithmetic, i.e., it allocates a
 *        SCIP_Real array of length two times SCIPgetNVars() for storing the coefficients.  In exact solving mode, we
 *        cannot use quad precision because we need to control the ronding mode, hence only the first SCIPgetNVars()
 *        entries are used.
 */
SCIP_RETCODE SCIPaggrRowCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow             /**< pointer to return aggregation row */
   )
{
   int nvars;
   assert(scip != NULL);
   assert(aggrrow != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, aggrrow) );

   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*aggrrow)->vals, QUAD_ARRAY_SIZE(nvars)) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*aggrrow)->inds, nvars) );

   BMSclearMemoryArray((*aggrrow)->vals, QUAD_ARRAY_SIZE(nvars));

   (*aggrrow)->local = FALSE;
   (*aggrrow)->nnz = 0;
   (*aggrrow)->rank = 0;
   QUAD_ASSIGN((*aggrrow)->rhs, 0.0);
   (*aggrrow)->rowsinds = NULL;
   (*aggrrow)->slacksign = NULL;
   (*aggrrow)->rowweights = NULL;
   (*aggrrow)->nrows = 0;
   (*aggrrow)->rowssize = 0;

   return SCIP_OKAY;
}

/** free a aggregation row */
void SCIPaggrRowFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow             /**< pointer to aggregation row that should be freed */
   )
{
   int nvars;

   assert(scip != NULL);
   assert(aggrrow != NULL);

   nvars = SCIPgetNVars(scip);

   SCIPfreeBlockMemoryArray(scip, &(*aggrrow)->inds, nvars);
   SCIPfreeBlockMemoryArray(scip, &(*aggrrow)->vals, QUAD_ARRAY_SIZE(nvars)); /*lint !e647*/
   SCIPfreeBlockMemoryArrayNull(scip, &(*aggrrow)->rowsinds, (*aggrrow)->rowssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*aggrrow)->slacksign, (*aggrrow)->rowssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(*aggrrow)->rowweights, (*aggrrow)->rowssize);
   SCIPfreeBlockMemory(scip, aggrrow);
}

/** output aggregation row to file stream */
void SCIPaggrRowPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< pointer to return aggregation row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP_VAR** vars;
   SCIP_MESSAGEHDLR* messagehdlr;
   int i;

   assert(scip != NULL);
   assert(aggrrow != NULL);

   vars = SCIPgetVars(scip);
   assert(vars != NULL);

   messagehdlr = SCIPgetMessagehdlr(scip);
   assert(messagehdlr);

   /* print coefficients */
   if( aggrrow->nnz == 0 )
      SCIPmessageFPrintInfo(messagehdlr, file, "0 ");

   for( i = 0; i < aggrrow->nnz; ++i )
   {
      SCIP_Real QUAD(val);

      QUAD_ARRAY_LOAD(val, aggrrow->vals, aggrrow->inds[i]);
      assert(SCIPvarGetProbindex(vars[aggrrow->inds[i]]) == aggrrow->inds[i]);
      SCIPmessageFPrintInfo(messagehdlr, file, "%+.15g<%s> ", QUAD_TO_DBL(val), SCIPvarGetName(vars[aggrrow->inds[i]]));
   }

   /* print right hand side */
   SCIPmessageFPrintInfo(messagehdlr, file, "<= %.15g\n", QUAD_TO_DBL(aggrrow->rhs));
}

/** copy a aggregation row */
SCIP_RETCODE SCIPaggrRowCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow,            /**< pointer to return aggregation row */
   SCIP_AGGRROW*         source              /**< source aggregation row */
   )
{
   int nvars;

   assert(scip != NULL);
   assert(aggrrow != NULL);
   assert(source != NULL);

   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBlockMemory(scip, aggrrow) );

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->vals, source->vals, QUAD_ARRAY_SIZE(nvars)) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->inds, source->inds, nvars) );
   (*aggrrow)->nnz = source->nnz;
   QUAD_ASSIGN_Q((*aggrrow)->rhs, source->rhs);

   if( source->nrows > 0 )
   {
      assert(source->rowsinds != NULL);
      assert(source->slacksign != NULL);
      assert(source->rowweights != NULL);

      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->rowsinds, source->rowsinds, source->nrows) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->slacksign, source->slacksign, source->nrows) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*aggrrow)->rowweights, source->rowweights, source->nrows) );
   }
   else
   {
      (*aggrrow)->rowsinds = NULL;
      (*aggrrow)->slacksign = NULL;
      (*aggrrow)->rowweights = NULL;
   }

   (*aggrrow)->nrows = source->nrows;
   (*aggrrow)->rowssize = source->nrows;
   (*aggrrow)->rank = source->rank;
   (*aggrrow)->local = source->local;

   return SCIP_OKAY;
}

/** add weighted row to aggregation row */
SCIP_RETCODE SCIPaggrRowAddRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_ROW*             row,                /**< row to add to aggregation row */
   SCIP_Real             weight,             /**< scale for adding given row to aggregation row */
   int                   sidetype            /**< specify row side type (-1 = lhs, 0 = automatic, 1 = rhs) */
   )
{
   SCIP_Real QUAD(quadprod);
   SCIP_Real sideval;
   SCIP_Bool uselhs;
   int i;

   assert(row->lppos >= 0);

   /* update local flag */
   aggrrow->local = aggrrow->local || row->local;

   /* update rank */
   aggrrow->rank = MAX(row->rank, aggrrow->rank);

   i = aggrrow->nrows++;

   if( aggrrow->nrows > aggrrow->rowssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, aggrrow->nrows);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowsinds, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->slacksign, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowweights, aggrrow->rowssize, newsize) );
      aggrrow->rowssize = newsize;
   }
   aggrrow->rowsinds[i] = SCIProwGetLPPos(row);
   aggrrow->rowweights[i] = weight;

   if( sidetype == -1 )
   {
      assert( ! SCIPisInfinity(scip, -row->lhs) );
      uselhs = TRUE;
   }
   else if( sidetype == 1 )
   {
      assert( ! SCIPisInfinity(scip, row->rhs) );
      uselhs = FALSE;
   }
   else
   {
      /* Automatically decide, whether we want to use the left or the right hand side of the row in the summation.
       * If possible, use the side that leads to a positive slack value in the summation.
       */
      if( SCIPisInfinity(scip, row->rhs) || (!SCIPisInfinity(scip, -row->lhs) && weight < 0.0) )
         uselhs = TRUE;
      else
         uselhs = FALSE;
   }

   if( uselhs )
   {
      aggrrow->slacksign[i] = -1;
      sideval = row->lhs - row->constant;
      if( row->integral )
         sideval = SCIPceil(scip, sideval); /* row is integral: round left hand side up */
   }
   else
   {
      aggrrow->slacksign[i] = +1;
      sideval = row->rhs - row->constant;
      if( row->integral )
         sideval = SCIPfloor(scip, sideval); /* row is integral: round right hand side up */
   }

   SCIPquadprecProdDD(quadprod, weight, sideval);
   SCIPquadprecSumQQ(aggrrow->rhs, aggrrow->rhs, quadprod);

   /* add up coefficients */
   SCIP_CALL( varVecAddScaledRowCoefsQuad(aggrrow->inds, aggrrow->vals, &aggrrow->nnz, row, weight) );

   return SCIP_OKAY;
}

/** add weighted row to aggregation row
 *
 *  @note this method is the variant of SCIPaggrRowAddRow that is safe to use in exact solving mode
 */
SCIP_RETCODE SCIPaggrRowAddRowSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_ROW*             row,                /**< row to add to aggregation row */
   SCIP_Real             weight,             /**< scale for adding given row to aggregation row */
   int                   sidetype,           /**< specify row side type (-1 = lhs, 0 = automatic, 1 = rhs) */
   SCIP_Bool*            success             /**< was the row added successfully */
   )
{
   SCIP_Real sideval;
   SCIP_Real sidevalchg;
   SCIP_Bool uselhs;
   SCIP_ROW* userow;
   SCIP_ROWEXACT* rowexact;
   SCIP_ROUNDMODE previousroundmode;
   int i;

   assert(SCIPisExact(scip));
   assert(success != NULL);

   /* update local flag */
   aggrrow->local = aggrrow->local || row->local;

   /* update rank */
   aggrrow->rank = MAX(row->rank, aggrrow->rank);

   i = aggrrow->nrows++;

   if( aggrrow->nrows > aggrrow->rowssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, aggrrow->nrows);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowsinds, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->slacksign, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowweights, aggrrow->rowssize, newsize) );
      aggrrow->rowssize = newsize;
   }
   aggrrow->rowsinds[i] = SCIProwGetLPPos(row);
   aggrrow->rowweights[i] = weight;

   if( sidetype == -1 )
   {
      assert(!SCIPisInfinity(scip, -row->lhs));
      uselhs = TRUE;
   }
   else if( sidetype == 1 )
   {
      assert(!SCIPisInfinity(scip, row->rhs));
      uselhs = FALSE;
   }
   else
   {
      /* Automatically decide, whether we want to use the left or the right hand side of the row in the summation.
       * If possible, use the side that leads to a positive slack value in the summation.
       */
      if( SCIPisInfinity(scip, row->rhs) || (!SCIPisInfinity(scip, -row->lhs) && weight < 0.0) )
         uselhs = TRUE;
      else
         uselhs = FALSE;
   }
   rowexact = SCIProwGetRowExact(row);
   if( !SCIProwExactHasFpRelax(rowexact) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   else if( SCIProwExactGetRowRhs(rowexact) != NULL && weight >= 0.0 )
      userow = SCIProwExactGetRowRhs(rowexact);
   else
      userow = row;

   aggrrow->slacksign[i] = uselhs ? -1 : 1;

   previousroundmode = SCIPintervalGetRoundingMode();

   if( uselhs )
   {
      SCIPintervalSetRoundingModeDownwards();
      sideval = userow->lhs - userow->constant;
#ifdef SCIP_DISABLED_CODE
      /* this is disabled because we can't certify it yet in exact solving mode; if enabled change also in cutsSubstituteMIRSafely() */
      /* row is integral? round left hand side up */
      if( userow->integral )
         sideval = ceil(sideval)
#endif
   }
   else
   {
      SCIPintervalSetRoundingModeUpwards();
      sideval = userow->rhs - userow->constant;
#ifdef SCIP_DISABLED_CODE
      /* this is disabled because we can't certify it yet in exact solving mode; if enabled change also in cutsSubstituteMIRSafely() */
      /* row is integral? round right hand side down */
      if( userow->integral )
         sideval = floor(sideval);
#endif
   }

   SCIPintervalSetRoundingModeUpwards();
   sidevalchg = QUAD_TO_DBL(aggrrow->rhs);
   sidevalchg += sideval * weight;
   QUAD_ASSIGN(aggrrow->rhs, sidevalchg);

   /* add up coefficients */
   *success = TRUE;
   SCIP_CALL( varVecAddScaledRowCoefsSafely(scip, aggrrow->inds, aggrrow->vals, &aggrrow->nnz, userow, weight, &sidevalchg, success) );

   sidevalchg += QUAD_TO_DBL(aggrrow->rhs);
   QUAD_ASSIGN(aggrrow->rhs, sidevalchg);

   SCIPintervalSetRoundingMode(previousroundmode);

   return SCIP_OKAY;
}

/** Removes a given variable @p var from position @p pos the aggregation row and updates the right-hand side according
 *  to sign of the coefficient, i.e., rhs -= coef * bound, where bound = lb if coef >= 0 and bound = ub, otherwise.
 *
 *  @note: The choice of global or local bounds depend on the validity (global or local) of the aggregation row.
 *
 *  @note: The list of non-zero indices will be updated by swapping the last non-zero index to @p pos.
 */
void SCIPaggrRowCancelVarWithBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_VAR*             var,                /**< variable that should be removed */
   int                   pos,                /**< position of the variable in the aggregation row */
   SCIP_Bool*            valid               /**< pointer to return whether the aggregation row is still valid */
   )
{
   SCIP_Real QUAD(val);
   int v;

   assert(valid != NULL);
   assert(pos >= 0);

   v = aggrrow->inds[pos];
   assert(v == SCIPvarGetProbindex(var));

   QUAD_ARRAY_LOAD(val, aggrrow->vals, v);

   *valid = TRUE;

   /* adjust left and right hand sides with max contribution */
   if( QUAD_TO_DBL(val) < 0.0 )
   {
      SCIP_Real ub = aggrrow->local ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);

      if( SCIPisInfinity(scip, ub) )
         QUAD_ASSIGN(aggrrow->rhs, SCIPinfinity(scip));
      else
      {
         SCIPquadprecProdQD(val, val, ub);
         SCIPquadprecSumQQ(aggrrow->rhs, aggrrow->rhs, -val);
      }
   }
   else
   {
      SCIP_Real lb = aggrrow->local ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);

      if( SCIPisInfinity(scip, -lb) )
         QUAD_ASSIGN(aggrrow->rhs, SCIPinfinity(scip));
      else
      {
         SCIPquadprecProdQD(val, val, lb);
         SCIPquadprecSumQQ(aggrrow->rhs, aggrrow->rhs, -val);
      }
   }

   QUAD_ASSIGN(val, 0.0);
   QUAD_ARRAY_STORE(aggrrow->vals, v, val);

   /* remove non-zero entry */
   --(aggrrow->nnz);
   aggrrow->inds[pos] = aggrrow->inds[aggrrow->nnz];

   if( SCIPisInfinity(scip, QUAD_HI(aggrrow->rhs)) )
      *valid = FALSE;
}

/** add the objective function with right-hand side @p rhs and scaled by @p scale to the aggregation row */
SCIP_RETCODE SCIPaggrRowAddObjectiveFunction(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_Real             rhs,                /**< right-hand side of the artificial row */
   SCIP_Real             scale               /**< scalar */
   )
{
   SCIP_VAR** vars;
   SCIP_Real QUAD(val);
   int nvars;

   assert(scip != NULL);
   assert(aggrrow != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* add all variables straight forward if the aggregation row is empty */
   if( aggrrow->nnz == 0 )
   {
      int i;
      for( i = 0; i < nvars; ++i )
      {
         assert(SCIPvarGetProbindex(vars[i]) == i);

         /* skip all variables with zero objective coefficient */
         if( SCIPisZero(scip, scale * SCIPvarGetObj(vars[i])) )
            continue;

         QUAD_ASSIGN(val, scale * SCIPvarGetObj(vars[i])); /*lint !e665*/
         QUAD_ARRAY_STORE(aggrrow->vals, i, val);
         aggrrow->inds[aggrrow->nnz++] = i;
      }

      /* add right-hand side value */
      QUAD_ASSIGN(aggrrow->rhs, scale * rhs); /*lint !e665*/
   }
   else
   {
      int i;
      SCIP_Real QUAD(quadprod);
      /* add the non-zeros to the aggregation row and keep non-zero index up to date */
      for( i = 0 ; i < nvars; ++i )
      {
         SCIP_Real varobj;
         assert(SCIPvarGetProbindex(vars[i]) == i);

         /* skip all variables with zero objective coefficient */
         if( SCIPisZero(scip, scale * SCIPvarGetObj(vars[i])) )
            continue;

         QUAD_ARRAY_LOAD(val, aggrrow->vals, i); /* val = aggrrow->vals[i] */

         if( QUAD_HI(val) == 0.0 )
            aggrrow->inds[aggrrow->nnz++] = i;

         varobj = SCIPvarGetObj(vars[i]);
         SCIPquadprecProdDD(quadprod, scale, varobj);
         SCIPquadprecSumQQ(val, val, quadprod);

         /* the value must not be exactly zero due to sparsity pattern */
         QUAD_HI(val) = NONZERO(QUAD_HI(val));
         assert(QUAD_HI(val) != 0.0);

         QUAD_ARRAY_STORE(aggrrow->vals, i, val);
      }

      /* add right-hand side value */
      SCIPquadprecProdDD(quadprod, scale, rhs);
      SCIPquadprecSumQQ(aggrrow->rhs, aggrrow->rhs, quadprod);
   }

   return SCIP_OKAY;
}

/** add weighted constraint to the aggregation row */
SCIP_RETCODE SCIPaggrRowAddCustomCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   int*                  inds,               /**< variable problem indices in constraint to add to the aggregation row */
   SCIP_Real*            vals,               /**< values of constraint to add to the aggregation row */
   int                   len,                /**< length of constraint to add to the aggregation row */
   SCIP_Real             rhs,                /**< right hand side of constraint to add to the aggregation row */
   SCIP_Real             weight,             /**< (positive) scale for adding given constraint to the aggregation row */
   int                   rank,               /**< rank to use for given constraint */
   SCIP_Bool             local               /**< is constraint only valid locally */
   )
{
   SCIP_Real QUAD(quadprod);
   int i;

   assert(weight >= 0.0);
   assert(!SCIPisInfinity(scip, REALABS(weight * rhs)));

   /* update local flag */
   aggrrow->local = aggrrow->local || local;

   /* update rank */
   aggrrow->rank = MAX(rank, aggrrow->rank);

   /* add right hand side value */
   SCIPquadprecProdDD(quadprod, weight, rhs);
   SCIPquadprecSumQQ(aggrrow->rhs, aggrrow->rhs, quadprod);

   /* add the non-zeros to the aggregation row and keep non-zero index up to date */
   for( i = 0 ; i < len; ++i )
   {
      SCIP_Real QUAD(val);
      int probindex = inds[i];

      QUAD_ARRAY_LOAD(val, aggrrow->vals, probindex); /* val = aggrrow->vals[probindex] */

      if( QUAD_HI(val) == 0.0 )
         aggrrow->inds[aggrrow->nnz++] = probindex;

      SCIPquadprecProdDD(quadprod, vals[i], weight);
      SCIPquadprecSumQQ(val, val, quadprod);

      /* the value must not be exactly zero due to sparsity pattern */
      QUAD_HI(val) = NONZERO(QUAD_HI(val));
      assert(QUAD_HI(val) != 0.0);

      QUAD_ARRAY_STORE(aggrrow->vals, probindex, val);
   }

   return SCIP_OKAY;
}

/** version for use in exact solving mode of SCIPaggrRowClear() */
void SCIPaggrRowClearSafely(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   int i;

   /* in exact solving mode, we do not use quad precision, because we need to control the rounding mode; hence, we only
    * use and clear the first SCIPgetNVars() entries
    */
   for( i = 0; i < aggrrow->nnz; ++i )
   {
      aggrrow->vals[aggrrow->inds[i]] = 0.0;
   }

   aggrrow->nnz = 0;
   aggrrow->nrows = 0;
   aggrrow->rank = 0;
   QUAD_ASSIGN(aggrrow->rhs, 0.0);
   aggrrow->local = FALSE;
}

/** clear all entries int the aggregation row but don't free memory */
void SCIPaggrRowClear(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   int i;
   SCIP_Real QUAD(tmp);

   QUAD_ASSIGN(tmp, 0.0);

   for( i = 0; i < aggrrow->nnz; ++i )
   {
      QUAD_ARRAY_STORE(aggrrow->vals, aggrrow->inds[i], tmp);
   }

   aggrrow->nnz = 0;
   aggrrow->nrows = 0;
   aggrrow->rank = 0;
   QUAD_ASSIGN(aggrrow->rhs, 0.0);
   aggrrow->local = FALSE;
}

/** calculates the efficacy norm of the given aggregation row, which depends on the "separating/efficacynorm" parameter
 *
 *  @return the efficacy norm of the given aggregation row, which depends on the "separating/efficacynorm" parameter
 */
SCIP_Real SCIPaggrRowCalcEfficacyNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   return calcEfficacyNormQuad(scip, aggrrow->vals, aggrrow->inds, aggrrow->nnz);
}

/** adds one row to the aggregation row
 *
 *  @note this method differs from SCIPaggrRowAddRow() by providing some additional parameters required for
 *        SCIPaggrRowSumRows()
 */
static
SCIP_RETCODE addOneRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_ROW*             row,                /**< the row to add */
   SCIP_Real             weight,             /**< weight of row to add */
   SCIP_Bool             sidetypebasis,      /**< choose sidetypes of row (lhs/rhs) based on basis information? */
   SCIP_Bool             allowlocal,         /**< should local rows allowed to be used? */
   int                   negslack,           /**< should negative slack variables allowed to be used? (0: no, 1: only for integral rows, 2: yes) */
   int                   maxaggrlen,         /**< maximal length of aggregation row */
   SCIP_Bool*            rowtoolong          /**< is the aggregated row too long */
   )
{
   SCIP_Real QUAD(quadprod);
   SCIP_Real sideval;
   SCIP_Bool uselhs;
   int i;

   assert( rowtoolong != NULL );
   *rowtoolong = FALSE;

   if( SCIPisFeasZero(scip, weight) || SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !allowlocal) )
   {
      return SCIP_OKAY;
   }

   if( sidetypebasis && !SCIPisEQ(scip, SCIProwGetLhs(row), SCIProwGetRhs(row)) )
   {
      SCIP_BASESTAT stat = SCIProwGetBasisStatus(row);

      if( stat == SCIP_BASESTAT_LOWER )
      {
         assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );
         uselhs = TRUE;
      }
      else if( stat == SCIP_BASESTAT_UPPER )
      {
         assert( ! SCIPisInfinity(scip, SCIProwGetRhs(row)) );
         uselhs = FALSE;
      }
      else if( SCIPisInfinity(scip, SCIProwGetRhs(row)) || (weight < 0.0 && ! SCIPisInfinity(scip, -SCIProwGetLhs(row))) )
         uselhs = TRUE;
      else
         uselhs = FALSE;
   }
   else if( (weight < 0.0 && !SCIPisInfinity(scip, -row->lhs)) || SCIPisInfinity(scip, row->rhs) )
      uselhs = TRUE;
   else
      uselhs = FALSE;

   if( uselhs )
   {
      assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );

      if( weight > 0.0 && ((negslack == 0) || (negslack == 1 && !row->integral)) )
         return SCIP_OKAY;

      sideval = row->lhs - row->constant;
      /* row is integral? round left hand side up */
      if( row->integral )
         sideval = SCIPceil(scip, sideval);
   }
   else
   {
      assert( ! SCIPisInfinity(scip, SCIProwGetRhs(row)) );

      if( weight < 0.0 && ((negslack == 0) || (negslack == 1 && !row->integral)) )
         return SCIP_OKAY;

      sideval = row->rhs - row->constant;
      /* row is integral? round right hand side down */
      if( row->integral )
         sideval = SCIPfloor(scip, sideval);
   }

   /* add right hand side, update rank and local flag */
   SCIPquadprecProdDD(quadprod, sideval, weight);
   SCIPquadprecSumQQ(aggrrow->rhs, aggrrow->rhs, quadprod);
   aggrrow->rank = MAX(aggrrow->rank, row->rank);
   aggrrow->local = aggrrow->local || row->local;

   /* ensure the array for storing the row information is large enough */
   i = aggrrow->nrows++;
   if( aggrrow->nrows > aggrrow->rowssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, aggrrow->nrows);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowsinds, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->slacksign, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowweights, aggrrow->rowssize, newsize) );
      aggrrow->rowssize = newsize;
   }

   /* add information of addditional row */
   aggrrow->rowsinds[i] = row->lppos;
   aggrrow->rowweights[i] = weight;
   aggrrow->slacksign[i] = uselhs ? -1 : 1;

   /* add up coefficients */
   SCIP_CALL( varVecAddScaledRowCoefsQuad(aggrrow->inds, aggrrow->vals, &aggrrow->nnz, row, weight) );

   /* check if row is too long now */
   if( aggrrow->nnz > maxaggrlen )
      *rowtoolong = TRUE;

   return SCIP_OKAY;
}

/** adds one row to the aggregation row
 *
 *  @note this method differs from SCIPaggrRowAddRowSafely() by providing some additional parameters required for
 *        SCIPaggrRowSumRows()
 *
 *  @note this method is the variant of addOneRow() that is safe to use in exact solving mode
 */
static
SCIP_RETCODE addOneRowSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_ROW*             row,                /**< the row to add */
   SCIP_Real             weight,             /**< weight of row to add */
   SCIP_Bool             sidetypebasis,      /**< choose sidetypes of row (lhs/rhs) based on basis information? */
   SCIP_Bool             allowlocal,         /**< should local rows allowed to be used? */
   int                   negslack,           /**< should negative slack variables allowed to be used? (0: no, 1: only for integral rows, 2: yes) */
   int                   maxaggrlen,         /**< maximal length of aggregation row */
   SCIP_Bool*            rowtoolong,         /**< is the aggregated row too long */
   SCIP_Bool*            rowused,            /**< was the row really added? */
   SCIP_Bool*            success,            /**< was the row added successfully? */
   SCIP_Bool*            lhsused             /**< was the lhs or the rhs of the row used? */
   )
{
   SCIP_Real sideval;
   SCIP_Real sidevalchg;
   SCIP_Bool uselhs;
   SCIP_ROW* userow;
   SCIP_ROWEXACT* rowexact;
   SCIP_ROUNDMODE previousroundmode;
   int i;

   assert(SCIPisExact(scip));
   assert(rowtoolong != NULL);
   *rowtoolong = FALSE;
   *rowused = FALSE;

   if( SCIPisFeasZero(scip, weight) || SCIProwIsModifiable(row) || (SCIProwIsLocal(row) && !allowlocal) )
   {
      return SCIP_OKAY;
   }

   if( sidetypebasis && !SCIPisEQ(scip, SCIProwGetLhs(row), SCIProwGetRhs(row)) )
   {
      SCIP_BASESTAT stat = SCIProwGetBasisStatus(row);

      if( stat == SCIP_BASESTAT_LOWER )
      {
         assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );
         uselhs = TRUE;
      }
      else if( stat == SCIP_BASESTAT_UPPER )
      {
         assert( ! SCIPisInfinity(scip, SCIProwGetRhs(row)) );
         uselhs = FALSE;
      }
      else if( SCIPisInfinity(scip, SCIProwGetRhs(row)) || (weight < 0.0 && ! SCIPisInfinity(scip, -SCIProwGetLhs(row))) )
         uselhs = TRUE;
      else
         uselhs = FALSE;
   }
   else if( (weight < 0.0 && !SCIPisInfinity(scip, -row->lhs)) || SCIPisInfinity(scip, row->rhs) )
      uselhs = TRUE;
   else
      uselhs = FALSE;

   rowexact = SCIProwGetRowExact(row);
   if( !SCIProwExactHasFpRelax(rowexact) )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }
   else if( SCIProwExactGetRowRhs(rowexact) != NULL && weight >= 0.0 )
      userow = SCIProwExactGetRowRhs(rowexact);
   else
      userow = row;

   previousroundmode = SCIPintervalGetRoundingMode();

   if( uselhs )
   {
      *lhsused = TRUE;
      assert( ! SCIPisInfinity(scip, -SCIProwGetLhs(row)) );

      if( weight > 0.0 && ((negslack == 0) || (negslack == 1 && !row->integral)) )
         return SCIP_OKAY;

      SCIPintervalSetRoundingModeDownwards();

      sideval = userow->lhs - userow->constant;
#ifdef SCIP_DISABLED_CODE
      /* this is disabled because we can't certify it yet in exact solving mode; if enabled change also in cutsSubstituteMIRSafely() */
      /* row is integral? round left hand side up */
      if( userow->integral )
         sideval = ceil(sideval);
#endif
   }
   else
   {
      *lhsused = FALSE;
      assert( ! SCIPisInfinity(scip, SCIProwGetRhs(row)) );

      if( weight < 0.0 && ((negslack == 0) || (negslack == 1 && !row->integral)) )
         return SCIP_OKAY;

      SCIPintervalSetRoundingModeUpwards();

      sideval = userow->rhs - userow->constant;
#ifdef SCIP_DISABLED_CODE
      /* this is disabled because we can't certify it yet in exact solving mode; if enabled change also in cutsSubstituteMIRSafely() */
      /* row is integral? round right hand side down */
      if( userow->integral )
         sideval = floor(sideval);
#endif
   }

   SCIPintervalSetRoundingModeUpwards();

   sidevalchg = QUAD_TO_DBL(aggrrow->rhs);
   sidevalchg += sideval * weight;
   QUAD_ASSIGN(aggrrow->rhs, sidevalchg);

   aggrrow->rank = MAX(aggrrow->rank, userow->rank);
   aggrrow->local = aggrrow->local || userow->local;

   /* ensure the array for storing the row information is large enough */
   i = aggrrow->nrows++;
   *rowused = TRUE;
   if( aggrrow->nrows > aggrrow->rowssize )
   {
      int newsize = SCIPcalcMemGrowSize(scip, aggrrow->nrows);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowsinds, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->slacksign, aggrrow->rowssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &aggrrow->rowweights, aggrrow->rowssize, newsize) );
      aggrrow->rowssize = newsize;
   }

   /* add information of addditional row */
   aggrrow->rowsinds[i] = row->lppos;
   aggrrow->rowweights[i] = weight;
   aggrrow->slacksign[i] = uselhs ? -1 : 1;

   /* add up coefficients */
   SCIP_CALL( varVecAddScaledRowCoefsSafely(scip, aggrrow->inds, aggrrow->vals, &aggrrow->nnz, userow, weight, &sidevalchg, success) );

   sidevalchg += QUAD_TO_DBL(aggrrow->rhs);
   QUAD_ASSIGN(aggrrow->rhs, sidevalchg);

   /* check if row is too long now */
   if( aggrrow->nnz > maxaggrlen )
      *rowtoolong = TRUE;

   SCIPintervalSetRoundingMode(previousroundmode);

   return SCIP_OKAY;
}

/** aggregate rows using the given weights; the current content of the aggregation row, \p aggrrow, is overwritten
 *
 *  @note this method is safe for usage in exact solving mode
 */
SCIP_RETCODE SCIPaggrRowSumRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  rowinds,            /**< array to store indices of non-zero entries of the weights array, or NULL */
   int                   nrowinds,           /**< number of non-zero entries in weights array, -1 if rowinds is NULL */
   SCIP_Bool             sidetypebasis,      /**< choose sidetypes of row (lhs/rhs) based on basis information? */
   SCIP_Bool             allowlocal,         /**< should local rows allowed to be used? */
   int                   negslack,           /**< should negative slack variables allowed to be used? (0: no, 1: only for integral rows, 2: yes) */
   int                   maxaggrlen,         /**< maximal length of aggregation row */
   SCIP_Bool*            valid               /**< is the aggregation valid */
   )
{
   SCIP_AGGRROW* certificaterow = NULL;
   SCIP_ROW** rows;
   SCIP_ROW** usedrows = NULL;
   SCIP_ROW** negslackrows = NULL;
   SCIP_VAR** vars;
   SCIP_Real* usedweights = NULL;
   SCIP_Real* negslackweights = NULL;
   int nrows;
   int nvars;
   int k;
   int nusedrows;
   int nnegslackrows;
   SCIP_Bool rowtoolong;
   SCIP_Bool rowused, rowusedcert, lhsused;

   assert( scip != NULL );
   assert( aggrrow != NULL );
   assert( valid != NULL );

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   if( SCIPisExact(scip) )
      SCIPaggrRowClearSafely(aggrrow);
   else
      SCIPaggrRowClear(aggrrow);
   *valid = TRUE;
   lhsused = FALSE;
   nusedrows = 0;
   nnegslackrows = 0;

   SCIPdebugMessage("Summing up %d rows in aggrrow \n", nrowinds);

   if( SCIPisCertified(scip) )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &usedrows, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &usedweights, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &negslackrows, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &negslackweights, nrows) );
      SCIP_CALL( SCIPaggrRowCreate(scip, &certificaterow) );
   }

   if( rowinds != NULL && nrowinds > -1 )
   {
      for( k = 0; k < nrowinds; ++k )
      {
         if( !SCIPisExact(scip) )
         {
            SCIP_CALL( addOneRow(scip, aggrrow, rows[rowinds[k]], weights[rowinds[k]], sidetypebasis, allowlocal, negslack, maxaggrlen, &rowtoolong) );

            if( rowtoolong )
               *valid = FALSE;
         }
         else  /*lint --e{644}*/
         {
            SCIPdebugMessage("Adding %g times row: ", weights[rowinds[k]]);
            SCIPdebug(SCIPprintRow(scip, rows[rowinds[k]], NULL));
            SCIP_CALL( addOneRowSafely(scip, aggrrow, rows[rowinds[k]], weights[rowinds[k]], sidetypebasis, allowlocal,
                  negslack, maxaggrlen, &rowtoolong, &rowused, valid, &lhsused) );

            if( rowtoolong )
               *valid = FALSE;

            if( !(*valid) )
               break;

            if( certificaterow != NULL )
            {
               assert(usedrows != NULL);
               assert(usedweights != NULL);
               assert(negslackrows != NULL);
               assert(negslackweights != NULL);

               SCIP_ROW* row = rows[rowinds[k]];
               SCIP_Bool integral = FALSE;

               /* just exclude the negative continuous slacks for the certificate rows */
               if( row->integral &&
                  ((!lhsused && SCIPrealIsExactlyIntegral(row->rhs) &&  SCIPrealIsExactlyIntegral(row->constant)) ||
                  (lhsused && SCIPrealIsExactlyIntegral(row->lhs) &&  SCIPrealIsExactlyIntegral(row->constant))) )
               {
                  SCIPdebugMessage("row has integral slack\n");
                  rowusedcert = FALSE;
                  integral = TRUE;
               }
               else
               {
                  /* the certificate row may exceed the limit maxaggrlen */
                  SCIP_CALL( addOneRowSafely(scip, certificaterow, rows[rowinds[k]], weights[rowinds[k]], sidetypebasis,
                        allowlocal, 0, nvars, &rowtoolong, &rowusedcert, valid, &lhsused) );
                  assert(!rowtoolong);
               }
               if( rowusedcert )
               {
                  usedrows[nusedrows] = rows[rowinds[k]];
                  usedweights[nusedrows] = weights[rowinds[k]];
                  nusedrows++;
               }
               if( rowused && !rowusedcert && !integral )
               {
                  SCIPdebugMessage("row has negative continous slack\n");
                  assert( (lhsused && weights[rowinds[k]] >= 0) || ((!lhsused) && weights[rowinds[k]] <= 0) || row->integral );
                  negslackrows[nnegslackrows] = rows[rowinds[k]];
                  negslackweights[nnegslackrows] = -weights[rowinds[k]];
                  nnegslackrows++;
               }
            }
         }

         if( !(*valid) )
            break;
      }
   }
   else
   {
      for( k = 0; k < nrows; ++k )
      {
         if( weights[k] != 0.0 )
         {
            if( !SCIPisExact(scip) )
            {
               SCIP_CALL( addOneRow(scip, aggrrow, rows[k], weights[k], sidetypebasis, allowlocal, negslack, maxaggrlen, &rowtoolong) );

               if( rowtoolong )
                  *valid = FALSE;
            }
            else
            {
               SCIPdebugMessage("Adding %g times row: ", weights[k]);
               SCIPdebug(SCIPprintRow(scip, rows[k], NULL));
               SCIP_CALL( addOneRowSafely(scip, aggrrow, rows[k], weights[k], sidetypebasis, allowlocal, negslack,
                     maxaggrlen, &rowtoolong, &rowused, valid, &lhsused) );

               if( rowtoolong )
                  *valid = FALSE;

               if( !(*valid) )
                  break;

               if( certificaterow != NULL )
               {
                  assert(usedrows != NULL);
                  assert(usedweights != NULL);
                  assert(negslackrows != NULL);
                  assert(negslackweights != NULL);

                  SCIP_ROW* row = rows[k];
                  SCIP_Bool integral = FALSE;

                  /* just exclude the negative continuous slacks for the certificate rows */
                  if( row->integral &&
                     ((!lhsused && SCIPrealIsExactlyIntegral(row->rhs) &&  SCIPrealIsExactlyIntegral(row->constant)) ||
                     (lhsused && SCIPrealIsExactlyIntegral(row->lhs) &&  SCIPrealIsExactlyIntegral(row->constant))) )
                  {
                     rowusedcert = FALSE;
                     SCIPdebugMessage("row has integral slack\n");
                     integral = TRUE;
                  }
                  else
                  {
                     /* the certificate row may exceed the limit maxaggrlen */
                     SCIP_CALL( addOneRowSafely(scip, certificaterow, rows[k], weights[k], sidetypebasis, allowlocal, 0,
                           nvars, &rowtoolong, &rowusedcert, valid, &lhsused) );
                     assert(!rowtoolong);
                  }
                  if( rowusedcert )
                  {
                     usedrows[nusedrows] = rows[k];
                     usedweights[nusedrows] = weights[k];
                     nusedrows++;
                  }
                  if( rowused && !rowusedcert && !integral )
                  {
                     SCIPdebugMessage("row has negative continous slack\n");
                     assert( (lhsused && weights[k] >= 0) || ((!lhsused) && weights[k] <= 0) || row->integral );
                     negslackrows[nnegslackrows] = rows[k];
                     negslackweights[nnegslackrows] = -weights[k];
                     nnegslackrows++;
                  }
               }
            }

            if( !(*valid) )
               break;
         }
      }
   }

   if( *valid )
   {
      SCIPaggrRowRemoveZeros(scip, aggrrow, FALSE, valid);

      if( certificaterow != NULL )
      {
         SCIPaggrRowRemoveZeros(scip, certificaterow, FALSE, valid);
         SCIP_CALL( SCIPaddCertificateAggrInfo(scip, certificaterow, usedrows, usedweights, certificaterow->nrows,
               negslackrows, negslackweights, nnegslackrows) );
      }
   }

   if( certificaterow != NULL )
   {
      SCIPaggrRowFree(scip, &certificaterow);
      SCIPfreeBufferArray(scip, &negslackweights);
      SCIPfreeBufferArray(scip, &negslackrows);
      SCIPfreeBufferArray(scip, &usedweights);
      SCIPfreeBufferArray(scip, &usedrows);
   }

   return SCIP_OKAY;
}

/** checks for cut redundancy and performs activity based coefficient tightening;
 *  removes coefficients that are zero with QUAD_EPSILON tolerance and uses variable bounds
 *  to remove small coefficients (relative to the maximum absolute coefficient)
 */
static
SCIP_RETCODE postprocessCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut a local cut */
   int*                  cutinds,            /**< variable problem indices of non-zeros in cut */
   SCIP_Real*            cutcoefs,           /**< non-zeros coefficients of cut */
   int*                  nnz,                /**< number non-zeros coefficients of cut */
   SCIP_Real*            cutrhs,             /**< right hand side of cut */
   SCIP_Bool*            success             /**< pointer to return whether post-processing was successful or cut is redundant */
   )
{
   int i;
   SCIP_Bool redundant;
   SCIP_Real maxcoef;
   SCIP_Real minallowedcoef;
   SCIP_Real QUAD(rhs);

   assert(scip != NULL);
   assert(cutinds != NULL);
   assert(cutcoefs != NULL);
   assert(cutrhs != NULL);
   assert(success != NULL);

   *success = FALSE;

   QUAD_ASSIGN(rhs, *cutrhs);

   if( removeZeros(scip, SCIPfeastol(scip), cutislocal, cutcoefs, QUAD(&rhs), cutinds, nnz) )
   {
      /* right hand side was changed to infinity -> cut is redundant */
      return SCIP_OKAY;
   }

   if( *nnz == 0 )
      return SCIP_OKAY;

   SCIP_CALL( cutTightenCoefs(scip, cutislocal, cutcoefs, QUAD(&rhs), cutinds, nnz, &redundant) );

   if( redundant )
   {
      /* cut is redundant */
      return SCIP_OKAY;
   }

   maxcoef = 0.0;
   for( i = 0; i < *nnz; ++i )
   {
      SCIP_Real absval = REALABS(cutcoefs[cutinds[i]]);
      maxcoef = MAX(absval, maxcoef);
   }

   maxcoef /= scip->set->sepa_maxcoefratio;
   minallowedcoef = SCIPsumepsilon(scip);
   minallowedcoef = MAX(minallowedcoef, maxcoef);

   *success = ! removeZeros(scip, minallowedcoef, cutislocal, cutcoefs, QUAD(&rhs), cutinds, nnz);
   *cutrhs = QUAD_TO_DBL(rhs);

   return SCIP_OKAY;
}


/** checks for cut redundancy and performs activity based coefficient tightening;
 *  removes coefficients that are zero with QUAD_EPSILON tolerance and uses variable bounds
 *  to remove small coefficients (relative to the maximum absolute coefficient).
 *  The cutcoefs must be a quad precision array, i.e. allocated with size
 *  QUAD_ARRAY_SIZE(nvars) and accessed with QUAD_ARRAY_LOAD and QUAD_ARRAY_STORE
 *  macros.
 */
static
SCIP_RETCODE postprocessCutQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut a local cut */
   int*                  cutinds,            /**< variable problem indices of non-zeros in cut */
   SCIP_Real*            cutcoefs,           /**< non-zeros coefficients of cut */
   int*                  nnz,                /**< number non-zeros coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< right hand side of cut */
   SCIP_Bool*            success             /**< pointer to return whether the cleanup was successful or if it is useless */
   )
{
   int i;
   SCIP_Bool redundant;
   SCIP_Real maxcoef;
   SCIP_Real minallowedcoef;

   assert(scip != NULL);
   assert(cutinds != NULL);
   assert(cutcoefs != NULL);
   assert(QUAD_HI(cutrhs) != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( removeZerosQuad(scip, SCIPfeastol(scip), cutislocal, cutcoefs, QUAD(cutrhs), cutinds, nnz) )
   {
      /* right hand side was changed to infinity -> cut is redundant */
      return SCIP_OKAY;
   }

   if( *nnz == 0 )
      return SCIP_OKAY;

   SCIP_CALL( cutTightenCoefsQuad(scip, cutislocal, cutcoefs, QUAD(cutrhs), cutinds, nnz, &redundant) );
   if( redundant )
   {
      /* cut is redundant */
      return SCIP_OKAY;
   }

   maxcoef = 0.0;
   for( i = 0; i < *nnz; ++i )
   {
      SCIP_Real abscoef;
      SCIP_Real QUAD(coef);
      QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[i]); /* coef = cutcoefs[cutinds[i]] */
      abscoef = REALABS(QUAD_TO_DBL(coef));
      maxcoef = MAX(abscoef, maxcoef);
   }

   maxcoef /= scip->set->sepa_maxcoefratio;
   minallowedcoef = SCIPsumepsilon(scip);
   minallowedcoef = MAX(minallowedcoef, maxcoef);

   *success = ! removeZerosQuad(scip, minallowedcoef, cutislocal, cutcoefs, QUAD(cutrhs), cutinds, nnz);

   return SCIP_OKAY;
}

/** checks for cut redundancy and performs activity based coefficient tightening;
 *  removes coefficients that are zero with QUAD_EPSILON tolerance and uses variable bounds
 *  to remove small coefficients (relative to the maximum absolute coefficient).
 *  The cutcoefs must be a quad precision array, i.e. allocated with size
 *  QUAD_ARRAY_SIZE(nvars) and accessed with QUAD_ARRAY_LOAD and QUAD_ARRAY_STORE
 *  macros.
 */
static
SCIP_RETCODE postprocessCutSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             cutislocal,         /**< is the cut a local cut */
   int*                  cutinds,            /**< variable problem indices of non-zeros in cut */
   SCIP_Real*            cutcoefs,           /**< non-zeros coefficients of cut */
   int*                  nnz,                /**< number non-zeros coefficients of cut */
   SCIP_Real*            cutrhs,             /**< right hand side of cut */
   SCIP_Bool*            success             /**< pointer to return whether the cleanup was successful or if it is useless */
   )
{
   int i;
   SCIP_Bool redundant;
   SCIP_Real maxcoef;
   SCIP_Real minallowedcoef;

   assert(SCIPisExact(scip));

   assert(scip != NULL);
   assert(cutinds != NULL);
   assert(cutcoefs != NULL);
   assert(cutrhs != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( removeZerosSafely(scip, SCIPfeastol(scip), cutcoefs, cutrhs, cutinds, nnz) )
   {
      /* right hand side was changed to infinity -> cut is redundant */
      return SCIP_OKAY;
   }

   if( *nnz == 0 )
      return SCIP_OKAY;

   SCIP_CALL( cutTightenCoefsSafely(scip, cutislocal, cutcoefs, cutrhs, cutinds, nnz, &redundant) );
   if( redundant )
   {
      /* cut is redundant */
      return SCIP_OKAY;
   }

   maxcoef = 0.0;
   for( i = 0; i < *nnz; ++i )
   {
      SCIP_Real abscoef;
      SCIP_Real coef;
      coef = cutcoefs[cutinds[i]];
      abscoef = REALABS(coef);
      maxcoef = MAX(abscoef, maxcoef);
   }

   maxcoef /= scip->set->sepa_maxcoefratio;
   minallowedcoef = SCIPsumepsilon(scip);
   minallowedcoef = MAX(minallowedcoef, maxcoef);

   *success = ! removeZerosSafely(scip, minallowedcoef, cutcoefs, cutrhs, cutinds, nnz);

   return SCIP_OKAY;
}

/** removes almost zero entries from the aggregation row. */
void SCIPaggrRowRemoveZeros(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_Bool             useglbbounds,       /**< consider global bound although the cut is local? */
   SCIP_Bool*            valid               /**< pointer to return whether the aggregation row is still valid */
   )
{
   assert(aggrrow != NULL);
   assert(valid != NULL);

   if( SCIPisExact(scip) )
   {
      SCIP_Real rhs;
      rhs = QUAD_TO_DBL(aggrrow->rhs);
      *valid = !removeZerosSafely(scip, SCIPsumepsilon(scip), aggrrow->vals, &rhs, aggrrow->inds, &aggrrow->nnz);
      QUAD_ASSIGN(aggrrow->rhs, rhs);
      return;
   }

   *valid = ! removeZerosQuad(scip, SCIPsumepsilon(scip), useglbbounds ? FALSE : aggrrow->local, aggrrow->vals,
      QUAD(&aggrrow->rhs), aggrrow->inds, &aggrrow->nnz);
}

/** get number of aggregated rows */
int SCIPaggrRowGetNRows(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->nrows;
}

/** get array with lp positions of rows used in aggregation */
int* SCIPaggrRowGetRowInds(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   assert(aggrrow != NULL);
   assert(aggrrow->rowsinds != NULL || aggrrow->nrows == 0);

   return aggrrow->rowsinds;
}

/** get array with weights of aggregated rows */
SCIP_Real* SCIPaggrRowGetRowWeights(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   )
{
   assert(aggrrow != NULL);
   assert(aggrrow->rowweights != NULL || aggrrow->nrows == 0);

   return aggrrow->rowweights;
}

/** checks whether a given row has been added to the aggregation row */
SCIP_Bool SCIPaggrRowHasRowBeenAdded(
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_ROW*             row                 /**< row for which it is checked whether it has been added to the aggregation */
   )
{
   int i;
   int rowind;

   assert(aggrrow != NULL);
   assert(row != NULL);

   rowind = SCIProwGetLPPos(row);

   for( i = 0; i < aggrrow->nrows; ++i )
   {
      if( aggrrow->rowsinds[i] == rowind )
         return TRUE;
   }

   return FALSE;
}

/** gets the array of corresponding variable problem indices for each non-zero in the aggregation row */
int* SCIPaggrRowGetInds(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->inds;
}

/** gets the number of non-zeros in the aggregation row */
int SCIPaggrRowGetNNz(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->nnz;
}

/** gets the rank of the aggregation row */
int SCIPaggrRowGetRank(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->rank;
}

/** checks if the aggregation row is only valid locally */
SCIP_Bool SCIPaggrRowIsLocal(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return aggrrow->local;
}

/** gets the right hand side of the aggregation row */
SCIP_Real SCIPaggrRowGetRhs(
   SCIP_AGGRROW*         aggrrow             /**< aggregation row */
   )
{
   assert(aggrrow != NULL);

   return QUAD_TO_DBL(aggrrow->rhs);
}

/* =========================================== c-MIR =========================================== */

#define MAXCMIRSCALE               1e+6 /**< maximal scaling (scale/(1-f0)) allowed in c-MIR calculations */

/* In order to derive cuts, we partition the variable array up in (not necessarily contiguous) sections.
 * The only requirement we place on these sections is that section i can only have variable bounds variables whose section
 * is strictly greater than i. This way, we can process the variable array in a 'linear' manner. */

/* @todo maintain a DAG for used varbounds and use topological ordering instead, this would also allow
 * variable bounds on variables of the same section to be used */

#define NSECTIONS 6

typedef struct MIR_Data
{
   int                   totalnnz; /* The total number of nonzeros in all of the sections */
   int*                  secindices[NSECTIONS]; /* The indices of the variables belonging to the section */
   int                   secnnz[NSECTIONS];  /* The number of nonzero indices in the section */

   SCIP_Bool             isenfint[NSECTIONS];/**< Does the section have an integrality constraint? */
   SCIP_Bool             isimplint[NSECTIONS];/**< Is the section implied integer variables? */

   /* Settings for cut derivation, per section */
   int                   usevbds[NSECTIONS]; /**< Should variable bound substitution be done for this section? */

   /* Problem data that we reuse often */
   SCIP_VAR**            vars;               /**< pointer to SCIPs variable array */
   int                   nvars;              /**< total number of variables */
   int                   nbinvars;           /**< total number of non-implint binary variables */
   int                   nintvars;           /**< total number of non-implint integer variables */
   int                   nbinimplvars;       /**< total number of implint binary variables */
   int                   nintimplvars;       /**< total number of implint integer variables */
   int                   ncontimplvars;      /**< total number of implint continuous variables */
   int                   ncontvars;          /**< total number of non-implied continuous variables */

   SCIP_Real*            cutcoefs;           /**< working cut indices value array */
   SCIP_Real             QUAD(cutrhs);       /**< the working right hand side of the cut*/

   int*                  cutinds;            /**< working cut variable problem index array */
   int                   ncutinds;           /**< number of values in the working cut variable problem index array */
} MIR_DATA;

/** Returns the section of a variable.
 *
 *  For now, this is equal to the variable type section of the variable in the problem.
 */
static
int varSection(
   MIR_DATA*             data,               /**< The MIR separation data */
   int                   probindex           /**< Problem index of a variable */
   )
{
   int limit;

   assert(data != NULL);

   limit = data->nvars - data->ncontvars;
   if( probindex >= limit )
      return 0;

   limit -= data->ncontimplvars;
   if( probindex >= limit )
      return 1;

   limit -= data->nintimplvars;
   if( probindex >= limit )
      return 2;

   limit -= data->nbinimplvars;
   if( probindex >= limit )
      return 3;

   limit -= data->nintvars;
   if( probindex >= limit )
      return 4;

   assert(limit == data->nbinvars);

   return 5;
}

/** finds the best lower bound of the variable to use for MIR transformation.
 *
 *  Currently, we use a slightly different function for the exact MIR cuts than for the normal MIR cuts due to differences
 *  in how the codes can handle variable bound substitution. This function can only be used with the safe MIR code. */
/*  @todo make behavior identical to the unsafe MIR cut computation */
static
SCIP_RETCODE findBestLbSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   int                   usevbds,            /**< should variable bounds be used in bound transformation? (0: no, 1: only binary, 2: all) */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestlb,             /**< pointer to store best bound value */
   SCIP_Real*            simplebound,        /**< pointer to store simple bound value */
   int*                  bestlbtype          /**< pointer to store best bound type (-2: local bound, -1: global bound, >= 0 variable bound index) */
   )
{
   assert(bestlb != NULL);
   assert(bestlbtype != NULL);
   assert(usevbds >= 0 && usevbds <= 2);

   *bestlb = SCIPvarGetLbGlobal(var);
   *bestlbtype = -1;

   if( allowlocal )
   {
      SCIP_Real loclb;

      loclb = SCIPvarGetLbLocal(var);
      if( SCIPisGT(scip, loclb, *bestlb) )
      {
         *bestlb = loclb;
         *bestlbtype = -2;
      }
   }

   *simplebound = *bestlb;

   if( usevbds && !SCIPvarIsIntegral(var) )
   {
      SCIP_Real bestvlb;
      int bestvlbidx;

      SCIP_CALL( SCIPgetVarClosestVlb(scip, var, sol, &bestvlb, &bestvlbidx) );
      if( bestvlbidx >= 0 && (bestvlb > *bestlb || (*bestlbtype < 0 && SCIPisGE(scip, bestvlb, *bestlb))) )
      {
         SCIP_VAR** vlbvars;
         SCIP_VAR* vlbvar;

         /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
         /**@todo this check is not needed for continuous variables; but allowing all but binary variables
          *       to be replaced by variable bounds seems to be buggy (wrong result on gesa2)
          */
         vlbvars = SCIPvarGetVlbVars(var);
         assert(vlbvars != NULL);
         vlbvar = vlbvars[bestvlbidx];
         assert(vlbvar != NULL);
         if( ( usevbds == 2 || ( SCIPvarGetType(vlbvar) == SCIP_VARTYPE_BINARY
            && !SCIPvarIsImpliedIntegral(vlbvar) ) )
            && SCIPvarGetProbindex(vlbvar) < SCIPvarGetProbindex(var) )
         {
            *bestlb = bestvlb;
            *bestlbtype = bestvlbidx;
         }
      }
   }

   return SCIP_OKAY;
}

/** finds the best upper bound of the variable to use for MIR transformation.
 * currently, we use a slightly different function for the exact MIR cuts than for the normal MIR cuts due to differences
 * in how the codes can handle variable bound substitution. This function can only be used with the safe MIR code. */
/* @todo make behavior identical to the unsafe MIR cut computation */
static
SCIP_RETCODE findBestUbSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   int                   usevbds,            /**< should variable bounds be used in bound transformation? (0: no, 1: only binary, 2: all) */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestub,             /**< pointer to store best bound value */
   SCIP_Real*            simplebound,        /**< pointer to store simple bound */
   int*                  bestubtype          /**< pointer to store best bound type (-2: local bound, -1: global bound, >= 0 variable bound index) */
   )
{
   assert(bestub != NULL);
   assert(bestubtype != NULL);

   *bestub = SCIPvarGetUbGlobal(var);
   *bestubtype = -1;

   if( allowlocal )
   {
      SCIP_Real locub;

      locub = SCIPvarGetUbLocal(var);
      if( SCIPisLT(scip, locub, *bestub) )
      {
         *bestub = locub;
         *bestubtype = -2;
      }
   }

   *simplebound = *bestub;

   if( usevbds && !SCIPvarIsIntegral(var) )
   {
      SCIP_Real bestvub;
      int bestvubidx;

      SCIP_CALL( SCIPgetVarClosestVub(scip, var, sol, &bestvub, &bestvubidx) );
      if( bestvubidx >= 0 && (bestvub < *bestub || (*bestubtype < 0 && SCIPisLE(scip, bestvub, *bestub))) )
      {
         SCIP_VAR** vubvars;
         SCIP_VAR* vubvar;

         /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
         /**@todo this check is not needed for continuous variables; but allowing all but binary variables
          *       to be replaced by variable bounds seems to be buggy (wrong result on gesa2)
          */
         vubvars = SCIPvarGetVubVars(var);
         assert(vubvars != NULL);
         vubvar = vubvars[bestvubidx];
         assert(vubvar != NULL);
         if( ( usevbds == 2 || ( SCIPvarGetType(vubvar) == SCIP_VARTYPE_BINARY
            && !SCIPvarIsImpliedIntegral(vubvar) ) )
            && SCIPvarGetProbindex(vubvar) < SCIPvarGetProbindex(var) )
         {
            *bestub = bestvub;
            *bestubtype = bestvubidx;
         }
      }
   }

   return SCIP_OKAY;
}

/** determine the best bounds with respect to the given solution for complementing the given variable */
/* @todo make behavior identical to the unsafe MIR cut computation */
static
SCIP_RETCODE determineBestBoundsSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to determine best bound for */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   int                   usevbds,            /**< should variable bounds be used in bound transformation? (0: no, 1: only binary, 2: all) */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   SCIP_Bool             ignoresol,          /**< should the LP solution be ignored? (eg, apply MIR to dualray) */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real*            bestlb,             /**< pointer to store best lower bound of variable */
   SCIP_Real*            bestub,             /**< pointer to store best upper bound of variable */
   int*                  bestlbtype,         /**< pointer to store type of the best lower bound of variable (-2: local bound, -1: global bound, >= 0 variable bound index) */
   int*                  bestubtype,         /**< pointer to store type of best upper bound of variable (-2: local bound, -1: global bound, >= 0 variable bound index) */
   SCIP_BOUNDTYPE*       selectedbound,      /**< pointer to store whether the lower bound or the upper bound should be preferred */
   SCIP_Bool*            freevariable        /**< pointer to store if this is a free variable */
   )
{
   SCIP_Real simplelb;
   SCIP_Real simpleub;
   int v;

   v = SCIPvarGetProbindex(var);

   /* check if the user specified a bound to be used */
   if( boundsfortrans != NULL && boundsfortrans[v] > -3 )
   {
      assert(!SCIPvarIsIntegral(var) || boundsfortrans[v] == -2 || boundsfortrans[v] == -1);
      assert(boundtypesfortrans != NULL);

      /* user has explicitly specified a bound to be used */
      if( boundtypesfortrans[v] == SCIP_BOUNDTYPE_LOWER )
      {
         /* user wants to use lower bound */
         *bestlbtype = boundsfortrans[v];
         if( *bestlbtype == -1 )
            *bestlb = SCIPvarGetLbGlobal(var); /* use global standard lower bound */
         else if( *bestlbtype == -2 )
            *bestlb = SCIPvarGetLbLocal(var);  /* use local standard lower bound */
         else
         {
            SCIP_VAR** vlbvars;
            SCIP_Real* vlbcoefs;
            SCIP_Real* vlbconsts;
            int k;

            assert(!ignoresol);

            /* use the given variable lower bound */
            vlbvars = SCIPvarGetVlbVars(var);
            vlbcoefs = SCIPvarGetVlbCoefs(var);
            vlbconsts = SCIPvarGetVlbConstants(var);
            k = boundsfortrans[v];
            assert(k >= 0 && k < SCIPvarGetNVlbs(var));
            assert(vlbvars != NULL);
            assert(vlbcoefs != NULL);
            assert(vlbconsts != NULL);

            *bestlb = vlbcoefs[k] * (sol == NULL ? SCIPvarGetLPSol(vlbvars[k]) : SCIPgetSolVal(scip, sol, vlbvars[k])) + vlbconsts[k];
         }

         assert(!SCIPisInfinity(scip, - *bestlb));
         *selectedbound = SCIP_BOUNDTYPE_LOWER;

         /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
         SCIP_CALL( findBestUbSafely(scip, var, sol, fixintegralrhs ? usevbds : 0, allowlocal && fixintegralrhs, bestub, &simpleub, bestubtype) );
      }
      else
      {
         assert(boundtypesfortrans[v] == SCIP_BOUNDTYPE_UPPER);

         /* user wants to use upper bound */
         *bestubtype = boundsfortrans[v];
         if( *bestubtype == -1 )
            *bestub = SCIPvarGetUbGlobal(var); /* use global standard upper bound */
         else if( *bestubtype == -2 )
            *bestub = SCIPvarGetUbLocal(var);  /* use local standard upper bound */
         else
         {
            SCIP_VAR** vubvars;
            SCIP_Real* vubcoefs;
            SCIP_Real* vubconsts;
            int k;

            assert(!ignoresol);

            /* use the given variable upper bound */
            vubvars = SCIPvarGetVubVars(var);
            vubcoefs = SCIPvarGetVubCoefs(var);
            vubconsts = SCIPvarGetVubConstants(var);
            k = boundsfortrans[v];
            assert(k >= 0 && k < SCIPvarGetNVubs(var));
            assert(vubvars != NULL);
            assert(vubcoefs != NULL);
            assert(vubconsts != NULL);

            /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
            *bestub = vubcoefs[k] * (sol == NULL ? SCIPvarGetLPSol(vubvars[k]) : SCIPgetSolVal(scip, sol, vubvars[k])) + vubconsts[k];
         }

         assert(!SCIPisInfinity(scip, *bestub));
         *selectedbound = SCIP_BOUNDTYPE_UPPER;

         /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
         SCIP_CALL( findBestLbSafely(scip, var, sol, fixintegralrhs ? usevbds : 0, allowlocal && fixintegralrhs, bestlb, &simplelb, bestlbtype) );
      }
   }
   else
   {
      SCIP_Real varsol;

      /* bound selection should be done automatically */

      /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
      SCIP_CALL( findBestLbSafely(scip, var, sol, usevbds, allowlocal, bestlb, &simplelb, bestlbtype) );

      /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
      SCIP_CALL( findBestUbSafely(scip, var, sol, usevbds, allowlocal, bestub, &simpleub, bestubtype) );

      /* check, if variable is free variable */
      if( SCIPisInfinity(scip, - *bestlb) && SCIPisInfinity(scip, *bestub) )
      {
         /* we found a free variable in the row with non-zero coefficient
            *  -> MIR row can't be transformed in standard form
            */
         *freevariable = TRUE;
         return SCIP_OKAY;
      }

      if( !ignoresol )
      {
         /* select transformation bound */
         varsol = (sol == NULL ? SCIPvarGetLPSol(var) : SCIPgetSolVal(scip, sol, var));

         if( SCIPisInfinity(scip, *bestub) ) /* if there is no ub, use lb */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( SCIPisInfinity(scip, - *bestlb) ) /* if there is no lb, use ub */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( SCIPisLT(scip, varsol, (1.0 - boundswitch) * (*bestlb) + boundswitch * (*bestub)) )
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( SCIPisGT(scip, varsol, (1.0 - boundswitch) * (*bestlb) + boundswitch * (*bestub)) )
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( *bestlbtype == -1 )  /* prefer global standard bounds */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( *bestubtype == -1 )  /* prefer global standard bounds */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( ((*bestlbtype) >= 0 || (*bestubtype) >= 0) && !SCIPisEQ(scip, *bestlb - simplelb, simpleub - *bestub) )
         {
            if( *bestlb - simplelb > simpleub - *bestub )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else
               *selectedbound = SCIP_BOUNDTYPE_UPPER;
         }
         else if( *bestlbtype >= 0 )   /* prefer variable bounds over local bounds */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( *bestubtype >= 0 )   /* prefer variable bounds over local bounds */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else                         /* no decision yet? just use lower bound */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
      }
      else
      {
         SCIP_Real glbub = SCIPvarGetUbGlobal(var);
         SCIP_Real glblb = SCIPvarGetLbGlobal(var);
         SCIP_Real distlb = REALABS(glblb - *bestlb);
         SCIP_Real distub = REALABS(glbub - *bestub);

         assert(!SCIPisInfinity(scip, - *bestlb) || !SCIPisInfinity(scip, *bestub));

         if( SCIPisInfinity(scip, - *bestlb) )
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( !SCIPisNegative(scip, *bestlb) )
         {
            if( SCIPisInfinity(scip, *bestub) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else if( SCIPisZero(scip, glblb) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else if( SCIPisLE(scip, distlb, distub) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else
               *selectedbound = SCIP_BOUNDTYPE_UPPER;
         }
         else
         {
            assert(!SCIPisInfinity(scip, - *bestlb));
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         }
      }
   }

   return SCIP_OKAY; /*lint !e438*/
}

/** finds the best lower bound of the variable to use for MIR transformation */
static
SCIP_RETCODE findBestLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   int                   usevbds,            /**< should variable bounds be used in bound transformation? (0: no, 1: only binary, 2: all) */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestlb,             /**< pointer to store best bound value */
   int*                  bestlbtype          /**< pointer to store best bound type (-2: local bound, -1: global bound, >= 0 variable bound index) */
   )
{
   assert(bestlb != NULL);
   assert(bestlbtype != NULL);
   assert(usevbds >= 0 && usevbds <= 2);

   *bestlb = SCIPvarGetLbGlobal(var);
   *bestlbtype = -1;

   if( allowlocal )
   {
      SCIP_Real loclb;

      loclb = SCIPvarGetLbLocal(var);
      if( SCIPisGT(scip, loclb, *bestlb) )
      {
         *bestlb = loclb;
         *bestlbtype = -2;
      }
   }

   if( usevbds > 0 )
   {
      SCIP_Real bestvlb;
      int bestvlbidx;

      SCIP_CALL( SCIPgetVarClosestVlb(scip, var, sol, &bestvlb, &bestvlbidx) );
      if( bestvlbidx >= 0 && (bestvlb > *bestlb || (*bestlbtype < 0 && SCIPisGE(scip, bestvlb, *bestlb))) )
      {
         SCIP_VAR** vlbvars;
         SCIP_VAR* vlbvar;

         /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
         /**@todo this check is not needed for continuous variables; but allowing all but binary variables
          *       to be replaced by variable bounds seems to be buggy (wrong result on gesa2)
          */
         vlbvars = SCIPvarGetVlbVars(var);
         assert(vlbvars != NULL);
         vlbvar = vlbvars[bestvlbidx];
         assert(vlbvar != NULL);
         if( ( usevbds == 2 || ( SCIPvarGetType(vlbvar) == SCIP_VARTYPE_BINARY
            && !SCIPvarIsImpliedIntegral(vlbvar) ) ) )
         {
            assert(SCIPvarGetProbindex(vlbvar) < SCIPvarGetProbindex(var));
            *bestlb = bestvlb;
            *bestlbtype = bestvlbidx;
         }
      }
   }

   return SCIP_OKAY;
}

/** finds the best upper bound of the variable to use for MIR transformation */
static
SCIP_RETCODE findBestUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   int                   usevbds,            /**< should variable bounds be used in bound transformation? (0: no, 1: only binary, 2: all) */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestub,             /**< pointer to store best bound value */
   int*                  bestubtype          /**< pointer to store best bound type (-2: local bound, -1: global bound, >= 0 variable bound index) */
   )
{
   assert(bestub != NULL);
   assert(bestubtype != NULL);
   assert(usevbds >= 0 && usevbds <= 2);

   *bestub = SCIPvarGetUbGlobal(var);
   *bestubtype = -1;

   if( allowlocal )
   {
      SCIP_Real locub;

      locub = SCIPvarGetUbLocal(var);
      if( SCIPisLT(scip, locub, *bestub) )
      {
         *bestub = locub;
         *bestubtype = -2;
      }
   }

   if( usevbds > 0 )
   {
      SCIP_Real bestvub;
      int bestvubidx;

      SCIP_CALL( SCIPgetVarClosestVub(scip, var, sol, &bestvub, &bestvubidx) );
      if( bestvubidx >= 0 && (bestvub < *bestub || (*bestubtype < 0 && SCIPisLE(scip, bestvub, *bestub))) )
      {
         SCIP_VAR** vubvars;
         SCIP_VAR* vubvar;

         /* we have to avoid cyclic variable bound usage, so we enforce to use only variable bounds variables of smaller index */
         /**@todo this check is not needed for continuous variables; but allowing all but binary variables
          *       to be replaced by variable bounds seems to be buggy (wrong result on gesa2)
          */
         vubvars = SCIPvarGetVubVars(var);
         assert(vubvars != NULL);
         vubvar = vubvars[bestvubidx];
         assert(vubvar != NULL);
         if( ( usevbds == 2 || ( SCIPvarGetType(vubvar) == SCIP_VARTYPE_BINARY
            && !SCIPvarIsImpliedIntegral(vubvar) ) ) )
         {
            assert( SCIPvarGetProbindex(vubvar) < SCIPvarGetProbindex(var) );
            *bestub = bestvub;
            *bestubtype = bestvubidx;
         }
      }
   }

   return SCIP_OKAY;
}


/** finds the best lower bound of the variable to use for MIR transformation.
 *  Differs from findBestLB() in that it allows more variable bound substitutions based on the variable sections. */
static
SCIP_RETCODE findMIRBestLb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   MIR_DATA*             data,               /**< the MIR data that specifies the variable sections */
   int                   usevbds,            /**< should variable bounds be used in bound transformation? (0: no, 1: only binary, 2: all) */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestlb,             /**< pointer to store best bound value */
   SCIP_Real*            simplebound,        /**< pointer to store simple bound value */
   int*                  bestlbtype          /**< pointer to store best bound type (-2: local bound, -1: global bound, >= 0 variable bound index) */
   )
{
   assert(bestlb != NULL);
   assert(bestlbtype != NULL);
   assert(usevbds >= 0 && usevbds <= 2);

   *bestlb = SCIPvarGetLbGlobal(var);
   *bestlbtype = -1;

   if( allowlocal )
   {
      SCIP_Real loclb;

      loclb = SCIPvarGetLbLocal(var);
      if( SCIPisGT(scip, loclb, *bestlb) )
      {
         *bestlb = loclb;
         *bestlbtype = -2;
      }
   }

   *simplebound = *bestlb;

   if( usevbds > 0 )
   {
      int nvlbs = SCIPvarGetNVlbs(var);

      if( nvlbs > 0 )
      {
         SCIP_Real bestvlb = SCIP_REAL_MIN;
         int bestvlbtype = -1;
         int boundedsection = varSection(data, SCIPvarGetProbindex(var));

         SCIP_VAR** vlbvars;
         SCIP_Real* vlbcoefs;
         SCIP_Real* vlbconsts;
         int i;

         vlbvars = SCIPvarGetVlbVars(var);
         vlbcoefs = SCIPvarGetVlbCoefs(var);
         vlbconsts = SCIPvarGetVlbConstants(var);

         /* search best VLB */
         for( i = 0; i < nvlbs; i++ )
         {
            /* For now, we only allow variable bounds from sections that are strictly greater to prevent cyclic usage.*/
            /** @todo: We don't use the caching mechanism of SCIPvarGetClosestVLB() because the cached variable bound
             * may be illegal. Building a local cache here may be worth it. */
            if( SCIPvarIsActive(vlbvars[i]) && boundedsection < varSection(data, SCIPvarGetProbindex(vlbvars[i])) &&
                  (usevbds == 2 || SCIPvarIsBinary(vlbvars[i])) )
            {
               SCIP_Real vlbsol;
               SCIP_Real vlbbnd;

               vlbsol = SCIPgetSolVal(scip, sol, vlbvars[i]);
               vlbbnd = vlbcoefs[i] * vlbsol + vlbconsts[i];

               if( vlbbnd > bestvlb )
               {
                  bestvlb = vlbbnd;
                  bestvlbtype = i;
               }
            }
         }

         if( bestvlbtype >= 0 && SCIPisGE(scip, bestvlb, *bestlb) )
         {
            *bestlb = bestvlb;
            *bestlbtype = bestvlbtype;
         }
      }
   }

   return SCIP_OKAY;
}

/** finds the best upper bound of the variable to use for MIR transformation.
 *  Differs from findBestUB() in that it allows more variable bound substitutions based on the variable sections. */
static
SCIP_RETCODE findMIRBestUb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   MIR_DATA*             data,               /**< the MIR data that specifies the variable sections */
   int                   usevbds,            /**< should variable bounds be used in bound transformation? (0: no, 1: only binary, 2: all) */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            bestub,             /**< pointer to store best bound value */
   SCIP_Real*            simplebound,        /**< pointer to store simple bound */
   int*                  bestubtype          /**< pointer to store best bound type (-2: local bound, -1: global bound, >= 0 variable bound index) */
   )
{
   assert(bestub != NULL);
   assert(bestubtype != NULL);
   assert(usevbds >= 0 && usevbds <= 2);

   *bestub = SCIPvarGetUbGlobal(var);
   *bestubtype = -1;

   if( allowlocal )
   {
      SCIP_Real locub;

      locub = SCIPvarGetUbLocal(var);
      if( SCIPisLT(scip, locub, *bestub) )
      {
         *bestub = locub;
         *bestubtype = -2;
      }
   }

   *simplebound = *bestub;

   if( usevbds > 0 )
   {
      int nvubs = SCIPvarGetNVubs(var);

      if( nvubs > 0 )
      {
         SCIP_Real bestvub = SCIP_REAL_MAX;
         int bestvubtype = -1;
         int boundedsection = varSection(data, SCIPvarGetProbindex(var));

         SCIP_VAR** vubvars;
         SCIP_Real* vubcoefs;
         SCIP_Real* vubconsts;
         int i;

         vubvars = SCIPvarGetVubVars(var);
         vubcoefs = SCIPvarGetVubCoefs(var);
         vubconsts = SCIPvarGetVubConstants(var);

         /* search best VUB */
         for( i = 0; i < nvubs; i++ )
         {
            /* For now, we only allow variable bounds from sections that are strictly greater to prevent cyclic usage.*/
            /** @todo: We don't use the caching mechanism of SCIPvarGetClosestVLB() because the cached variable bound
             * may be illegal. Building a local cache here may be worth it. */
            if( SCIPvarIsActive(vubvars[i]) && boundedsection < varSection(data, SCIPvarGetProbindex(vubvars[i])) &&
               (usevbds == 2 || SCIPvarIsBinary(vubvars[i])) )
            {
               SCIP_Real vubsol;
               SCIP_Real vubbnd;

               vubsol = SCIPgetSolVal(scip, sol, vubvars[i]);
               vubbnd = vubcoefs[i] * vubsol + vubconsts[i];

               if( vubbnd < bestvub )
               {
                  bestvub = vubbnd;
                  bestvubtype = i;
               }
            }
         }

         if( bestvubtype >= 0 && SCIPisLE(scip, bestvub, *bestub) )
         {
            *bestub = bestvub;
            *bestubtype = bestvubtype;
         }
      }
   }

   return SCIP_OKAY;
}

/** determine the best bounds with respect to the given solution for complementing the given variable */
static
SCIP_RETCODE determineBestBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to determine best bound for */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   MIR_DATA*             data,               /**< the MIR data that specifies the variable sections */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   int                   usevbds,            /**< should variable bounds be used in bound transformation? (0: no, 1: only binary, 2: all) */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   SCIP_Bool             ignoresol,          /**< should the LP solution be ignored? (eg, apply MIR to dualray) */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real*            bestlb,             /**< pointer to store best lower bound of variable */
   SCIP_Real*            bestub,             /**< pointer to store best upper bound of variable */
   int*                  bestlbtype,         /**< pointer to store type of best lower bound of variable (-2: local bound, -1: global bound, >= 0 variable bound index) */
   int*                  bestubtype,         /**< pointer to store type of best upper bound of variable (-2: local bound, -1: global bound, >= 0 variable bound index) */
   SCIP_BOUNDTYPE*       selectedbound,      /**< pointer to store whether the lower bound or the upper bound should be preferred */
   SCIP_Bool*            freevariable        /**< pointer to store if this is a free variable */
   )
{
   SCIP_Real simplelb;
   SCIP_Real simpleub;
   int v;

   assert(usevbds >= 0 && usevbds <= 2);

   v = SCIPvarGetProbindex(var);

   /* check if the user specified a bound to be used */
   if( boundsfortrans != NULL && boundsfortrans[v] > -3 )
   {
      assert(!SCIPvarIsIntegral(var) || boundsfortrans[v] == -2 || boundsfortrans[v] == -1);
      assert(boundtypesfortrans != NULL);

      /* user has explicitly specified a bound to be used */
      if( boundtypesfortrans[v] == SCIP_BOUNDTYPE_LOWER )
      {
         /* user wants to use lower bound */
         *bestlbtype = boundsfortrans[v];
         if( *bestlbtype == -1 )
            *bestlb = SCIPvarGetLbGlobal(var); /* use global standard lower bound */
         else if( *bestlbtype == -2 )
            *bestlb = SCIPvarGetLbLocal(var);  /* use local standard lower bound */
         else
         {
            SCIP_Real vlbsol;
            SCIP_VAR** vlbvars;
            SCIP_Real* vlbcoefs;
            SCIP_Real* vlbconsts;
            int k;

            assert(!ignoresol);

            /* use the given variable lower bound */
            vlbvars = SCIPvarGetVlbVars(var);
            vlbcoefs = SCIPvarGetVlbCoefs(var);
            vlbconsts = SCIPvarGetVlbConstants(var);
            k = boundsfortrans[v];
            assert(k >= 0 && k < SCIPvarGetNVlbs(var));
            assert(vlbvars != NULL);
            assert(vlbcoefs != NULL);
            assert(vlbconsts != NULL);

            vlbsol = SCIPgetSolVal(scip, sol, vlbvars[k]);

            *bestlb = vlbcoefs[k] * vlbsol + vlbconsts[k];
         }

         assert(!SCIPisInfinity(scip, - *bestlb));
         *selectedbound = SCIP_BOUNDTYPE_LOWER;

         /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
         SCIP_CALL( findMIRBestUb(scip, var, sol, data, fixintegralrhs ? usevbds : 0, allowlocal && fixintegralrhs, bestub, &simpleub, bestubtype) );
      }
      else
      {
         assert(boundtypesfortrans[v] == SCIP_BOUNDTYPE_UPPER);

         /* user wants to use upper bound */
         *bestubtype = boundsfortrans[v];
         if( *bestubtype == -1 )
            *bestub = SCIPvarGetUbGlobal(var); /* use global standard upper bound */
         else if( *bestubtype == -2 )
            *bestub = SCIPvarGetUbLocal(var);  /* use local standard upper bound */
         else
         {
            SCIP_Real vubsol;
            SCIP_VAR** vubvars;
            SCIP_Real* vubcoefs;
            SCIP_Real* vubconsts;
            int k;

            assert(!ignoresol);

            /* use the given variable upper bound */
            vubvars = SCIPvarGetVubVars(var);
            vubcoefs = SCIPvarGetVubCoefs(var);
            vubconsts = SCIPvarGetVubConstants(var);
            k = boundsfortrans[v];
            assert(k >= 0 && k < SCIPvarGetNVubs(var));
            assert(vubvars != NULL);
            assert(vubcoefs != NULL);
            assert(vubconsts != NULL);

            vubsol = SCIPgetSolVal(scip, sol, vubvars[k]);
            *bestub = vubcoefs[k] * vubsol + vubconsts[k];
         }

         assert(!SCIPisInfinity(scip, *bestub));
         *selectedbound = SCIP_BOUNDTYPE_UPPER;

         /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
         SCIP_CALL( findMIRBestLb(scip, var, sol, data, fixintegralrhs ? usevbds : 0, allowlocal && fixintegralrhs, bestlb, &simplelb, bestlbtype) );
      }
   }
   else
   {
      SCIP_Real varsol;

      /* bound selection should be done automatically */

      /* find closest lower bound in standard lower bound (and variable lower bounds for continuous variables) */
      SCIP_CALL( findMIRBestLb(scip, var, sol, data, usevbds, allowlocal, bestlb, &simplelb, bestlbtype) );

      /* find closest upper bound in standard upper bound (and variable upper bounds for continuous variables) */
      SCIP_CALL( findMIRBestUb(scip, var, sol, data, usevbds, allowlocal, bestub, &simpleub, bestubtype) );

      /* check, if variable is free variable */
      if( SCIPisInfinity(scip, - *bestlb) && SCIPisInfinity(scip, *bestub) )
      {
         /* we found a free variable in the row with non-zero coefficient
            *  -> MIR row can't be transformed in standard form
            */
         *freevariable = TRUE;
         return SCIP_OKAY;
      }

      if( !ignoresol )
      {
         /* select transformation bound */
         varsol = SCIPgetSolVal(scip, sol, var);

         if( SCIPisInfinity(scip, *bestub) ) /* if there is no ub, use lb */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( SCIPisInfinity(scip, - *bestlb) ) /* if there is no lb, use ub */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( SCIPisLT(scip, varsol, (1.0 - boundswitch) * (*bestlb) + boundswitch * (*bestub)) )
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( SCIPisGT(scip, varsol, (1.0 - boundswitch) * (*bestlb) + boundswitch * (*bestub)) )
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( *bestlbtype == -1 )  /* prefer global standard bounds */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( *bestubtype == -1 )  /* prefer global standard bounds */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( ((*bestlbtype) >= 0 || (*bestubtype) >= 0) && !SCIPisEQ(scip, *bestlb - simplelb, simpleub - *bestub) )
         {
            if( *bestlb - simplelb > simpleub - *bestub )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else
               *selectedbound = SCIP_BOUNDTYPE_UPPER;
         }
         else if( *bestlbtype >= 0 )   /* prefer variable bounds over local bounds */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         else if( *bestubtype >= 0 )   /* prefer variable bounds over local bounds */
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else                         /* no decision yet? just use lower bound */
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
      }
      else
      {
         SCIP_Real glbub = SCIPvarGetUbGlobal(var);
         SCIP_Real glblb = SCIPvarGetLbGlobal(var);
         SCIP_Real distlb = REALABS(glblb - *bestlb);
         SCIP_Real distub = REALABS(glbub - *bestub);

         assert(!SCIPisInfinity(scip, - *bestlb) || !SCIPisInfinity(scip, *bestub));

         if( SCIPisInfinity(scip, - *bestlb) )
            *selectedbound = SCIP_BOUNDTYPE_UPPER;
         else if( !SCIPisNegative(scip, *bestlb) )
         {
            if( SCIPisInfinity(scip, *bestub) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else if( SCIPisZero(scip, glblb) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else if( SCIPisLE(scip, distlb, distub) )
               *selectedbound = SCIP_BOUNDTYPE_LOWER;
            else
               *selectedbound = SCIP_BOUNDTYPE_UPPER;
         }
         else
         {
            assert(!SCIPisInfinity(scip, - *bestlb));
            *selectedbound = SCIP_BOUNDTYPE_LOWER;
         }
      }
   }

   return SCIP_OKAY; /*lint !e438*/
}

/** Performs bound substitution for a MIR cut */
static
void doMIRBoundSubstitution(
   SCIP*                 scip,               /**< SCIP datastructure */
   MIR_DATA*             data,               /**< the MIR data structure for this cut */
   int                   varsign,            /**< stores the sign of the transformed variable in summation */
   int                   boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Real             boundval,           /**< array of best bound to be used for the substitution for each nonzero index */
   int                   probindex,          /**< problem index of variable to perform the substitution step for */
   SCIP_Bool*            localbdsused        /**< pointer to updated whether a local bound was used for substitution */
   )
{
   SCIP_Real QUAD(coef);
   SCIP_Real QUAD(tmp);

   assert(!SCIPisInfinity(scip, -varsign * boundval));

   QUAD_ARRAY_LOAD(coef, data->cutcoefs, probindex);

   /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
   if( boundtype < 0 )
   {
      SCIPquadprecProdQD(tmp, coef, boundval);
      SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, -tmp);
      *localbdsused = *localbdsused || ( boundtype == -2 );
   }
   else
   {
      SCIP_VAR** vbdvars;
      SCIP_Real* vbdcoefs;
      SCIP_Real* vbdconsts;
      SCIP_Real QUAD(zcoef);
      int zidx;
      SCIP_VAR* var = SCIPgetVars(scip)[probindex];

      if( varsign == +1 )
      {
         vbdvars = SCIPvarGetVlbVars(var);
         vbdcoefs = SCIPvarGetVlbCoefs(var);
         vbdconsts = SCIPvarGetVlbConstants(var);
         assert(0 <= boundtype && boundtype < SCIPvarGetNVlbs(var));
      }
      else
      {
         vbdvars = SCIPvarGetVubVars(var);
         vbdcoefs = SCIPvarGetVubCoefs(var);
         vbdconsts = SCIPvarGetVubConstants(var);
         assert(0 <= boundtype && boundtype < SCIPvarGetNVubs(var));
      }

      assert(vbdvars != NULL);
      assert(vbdcoefs != NULL);
      assert(vbdconsts != NULL);
      assert(SCIPvarIsActive(vbdvars[boundtype]));

      zidx = SCIPvarGetProbindex(vbdvars[boundtype]);

      SCIPquadprecProdQD(tmp, coef, vbdconsts[boundtype]);
      SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, -tmp);

      /* check if integral variable already exists in the row */
      QUAD_ARRAY_LOAD(zcoef, data->cutcoefs, zidx);

      /* If it is new, add it to the indices */
      if( QUAD_HI(zcoef) == 0.0 )
      {
         int section = varSection(data, zidx);
         assert(section > varSection(data, probindex));

         data->secindices[section][data->secnnz[section]] = zidx;
         ++data->secnnz[section];
         ++data->totalnnz;
      }

      SCIPquadprecProdQD(tmp, coef, vbdcoefs[boundtype]);
      SCIPquadprecSumQQ(zcoef, zcoef, tmp);

      QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
      assert(QUAD_HI(zcoef) != 0.0);

      QUAD_ARRAY_STORE(data->cutcoefs, zidx, zcoef);
   }
}

/** performs the bound substitution step with the given variable or simple bounds for the variable with the given problem index
 *
 *  @note this method is safe for usage in exact solving mode
 *
 *  @todo make behavior identical to the unsafe MIR cut computation
 */
static
void performBoundSubstitutionSafely(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*            cutcoefs,           /**< array of cut coefficients */
   SCIP_Real*            cutrhs,             /**< pointer to right hand side of the cut */
   int                   varsign,            /**< stores the sign of the transformed variable in summation */
   int                   boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Real             boundval,           /**< array of best bound to be used for the substitution for each nonzero index */
   int                   probindex,          /**< problem index of variable to perform the substitution step for */
   SCIP_Bool*            localbdsused        /**< pointer to updated whether a local bound was used for substitution */
   )
{
   SCIP_Real coef;
   SCIP_ROUNDMODE previousroundmode;

   assert(!SCIPisInfinity(scip, -varsign * boundval));
   assert(SCIPisExact(scip));

   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();

   coef = cutcoefs[probindex];

   /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
   if( boundtype < 0 )
   {
      *cutrhs += coef * (-boundval);
      *localbdsused = *localbdsused || (boundtype == -2);
   }
   else
   {
      /* we don't support vlbs in exact mode yet */
      assert(!SCIPisExact(scip));
      SCIPerrorMessage("variable lower bounds not implemented in exact solving mode yet \n");
      SCIPABORT();
   }

   SCIPintervalSetRoundingMode(previousroundmode); /*lint !e644*/
}

/** performs the bound substitution step with the simple bound for the variable with the given problem index
 *
 *  @note this method is safe for usage in exact solving mode
 *
 *  @todo make behavior identical to the unsafe MIR cut computation
 */
static
void performBoundSubstitutionSimpleSafely(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*            cutcoefs,           /**< array of cut coefficients */
   SCIP_Real*            cutrhs,             /**< pointer to right hand side of the cut */
   int                   boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Real             boundval,           /**< array of best bound to be used for the substitution for each nonzero index */
   int                   probindex,          /**< problem index of variable to perform the substitution step for */
   SCIP_Bool*            localbdsused        /**< pointer to updated whether a local bound was used for substitution */
   )
{
   SCIP_Real coef;
   SCIP_ROUNDMODE previousroundmode;

   assert(!SCIPisInfinity(scip, ABS(boundval)));
   assert(SCIPisExact(scip));

   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();

   coef = cutcoefs[probindex];

   /* must be a standard bound */
   assert( boundtype < 0 );

   *cutrhs += coef * (-boundval);

   *localbdsused = *localbdsused || (boundtype == -2);

   SCIPintervalSetRoundingMode(previousroundmode); /*lint !e644*/
}

/** performs the bound substitution step with the given variable or simple bounds for the variable with the given problem index */
static
void performBoundSubstitution(
   SCIP*                 scip,               /**< SCIP datastructure */
   int*                  cutinds,            /**< index array of nonzeros in the cut */
   SCIP_Real*            cutcoefs,           /**< array of cut coefficients */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of the cut */
   int*                  nnz,                /**< pointer to number of nonzeros of the cut */
   int                   varsign,            /**< stores the sign of the transformed variable in summation */
   int                   boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Real             boundval,           /**< array of best bound to be used for the substitution for each nonzero index */
   int                   probindex,          /**< problem index of variable to perform the substitution step for */
   SCIP_Bool*            localbdsused        /**< pointer to updated whether a local bound was used for substitution */
   )
{
   SCIP_Real QUAD(coef);
   SCIP_Real QUAD(tmp);

   assert(!SCIPisInfinity(scip, -varsign * boundval));
   assert(!SCIPisExact(scip));

   QUAD_ARRAY_LOAD(coef, cutcoefs, probindex);

   /* standard (bestlbtype < 0) or variable (bestlbtype >= 0) lower bound? */
   if( boundtype < 0 )
   {
      SCIPquadprecProdQD(tmp, coef, boundval);
      SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
      *localbdsused = *localbdsused || (boundtype == -2);
   }
   else
   {
      SCIP_VAR** vbdvars;
      SCIP_Real* vbdcoefs;
      SCIP_Real* vbdconsts;
      SCIP_Real QUAD(zcoef);
      int zidx;
      SCIP_VAR* var = SCIPgetVars(scip)[probindex];

      if( varsign == +1 )
      {
         vbdvars = SCIPvarGetVlbVars(var);
         vbdcoefs = SCIPvarGetVlbCoefs(var);
         vbdconsts = SCIPvarGetVlbConstants(var);
         assert(0 <= boundtype && boundtype < SCIPvarGetNVlbs(var));
      }
      else
      {
         vbdvars = SCIPvarGetVubVars(var);
         vbdcoefs = SCIPvarGetVubCoefs(var);
         vbdconsts = SCIPvarGetVubConstants(var);
         assert(0 <= boundtype && boundtype < SCIPvarGetNVubs(var));
      }

      assert(vbdvars != NULL);
      assert(vbdcoefs != NULL);
      assert(vbdconsts != NULL);
      assert(SCIPvarIsActive(vbdvars[boundtype]));

      zidx = SCIPvarGetProbindex(vbdvars[boundtype]);

      SCIPquadprecProdQD(tmp, coef, vbdconsts[boundtype]);
      SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);

      /* check if integral variable already exists in the row */
      QUAD_ARRAY_LOAD(zcoef, cutcoefs, zidx);

      if( QUAD_HI(zcoef) == 0.0 )
         cutinds[(*nnz)++] = zidx;

      SCIPquadprecProdQD(tmp, coef, vbdcoefs[boundtype]);
      SCIPquadprecSumQQ(zcoef, zcoef, tmp);

      QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
      assert(QUAD_HI(zcoef) != 0.0);

      QUAD_ARRAY_STORE(cutcoefs, zidx, zcoef);
   }
}

/** performs the bound substitution step with the simple bound for the variable with the given problem index */
static
void performBoundSubstitutionSimple(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*            cutcoefs,           /**< array of cut coefficients */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of the cut */
   int                   boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Real             boundval,           /**< array of best bound to be used for the substitution for each nonzero index */
   int                   probindex,          /**< problem index of variable to perform the substitution step for */
   SCIP_Bool*            localbdsused        /**< pointer to updated whether a local bound was used for substitution */
   )
{
   SCIP_Real QUAD(coef);
   SCIP_Real QUAD(tmp);

   assert(!SCIPisInfinity(scip, ABS(boundval)));
   assert(!SCIPisExact(scip));

   QUAD_ARRAY_LOAD(coef, cutcoefs, probindex);

   /* must be a standard bound */
   assert( boundtype < 0 );

   SCIPquadprecProdQD(tmp, coef, boundval);
   SCIPquadprecSumQQ(*cutrhs, *cutrhs, -tmp);
   *localbdsused = *localbdsused || (boundtype == -2);
}

/** Transform equation \f$ a \cdot x = b; lb \leq x \leq ub \f$ into standard form
 *    \f$ a^\prime \cdot x^\prime = b,\; 0 \leq x^\prime \leq ub' \f$.
 *
 *  Transform variables (lb or ub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \mbox{if lb is used in transformation},\\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if ub is used in transformation},
 *  \end{array}
 *  \f]
 *  and move the constant terms \f$ a_j\, lb_j \f$ or \f$ a_j\, ub_j \f$ to the rhs.
 *
 *  Transform variables (vlb or vub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - (bl_j\, zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \mbox{if vlb is used in transf.} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if vub is used in transf.}
 *  \end{array}
 *  \f]
 *  move the constant terms \f$ a_j\, dl_j \f$ or \f$ a_j\, du_j \f$ to the rhs, and update the coefficient of the VLB variable:
 *  \f[
 *  \begin{array}{ll}
 *    a_{zl_j} := a_{zl_j} + a_j\, bl_j,& \mbox{or} \\
 *    a_{zu_j} := a_{zu_j} + a_j\, bu_j &
 *  \end{array}
 *  \f]
 *
 *  @note this method is safe for usage in exact solving mode
 *
 *  @todo make behavior identical to the unsafe MIR cut computation
 */
static
SCIP_RETCODE cutsTransformMIRSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   SCIP_Bool             ignoresol,          /**< should the LP solution be ignored? (eg, apply MIR to dualray) */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   SCIP_Real*            cutrhs,             /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in MIR row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{ /*lint --e{644}*/
   SCIP_Real* bestlbs;
   SCIP_Real* bestubs;
   int* bestlbtypes;
   int* bestubtypes;
   SCIP_BOUNDTYPE* selectedbounds;
   int i;
   int aggrrowintstart;
   int nvars;
   int firstcontvar;
   SCIP_VAR** vars;
   SCIP_MIRINFO* mirinfo = NULL;

   SCIP_ROUNDMODE previousroundmode;

   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);
   assert(SCIPisExact(scip));

   if( SCIPisCertified(scip)   )
      mirinfo = SCIPgetCertificate(scip)->mirinfo[SCIPgetCertificate(scip)->nmirinfos - 1];
   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();

   *freevariable = FALSE;
   *localbdsused = FALSE;

   /* allocate temporary memory to store best bounds and bound types */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbs, 2*(*nnz)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubs, 2*(*nnz)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbtypes, 2*(*nnz)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubtypes, 2*(*nnz)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selectedbounds, 2*(*nnz)) );

   /* start with continuous variables, because using variable bounds can affect the untransformed integral
    * variables, and these changes have to be incorporated in the transformation of the integral variables
    * (continuous variables have largest problem indices!)
    */
   SCIPsortDownInt(cutinds, *nnz);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   firstcontvar = nvars - SCIPgetNContVars(scip);

   /* determine the best bounds for the continuous variables */
   for( i = 0; i < *nnz && cutinds[i] >= firstcontvar; ++i )
   {
      SCIP_CALL( determineBestBoundsSafely(scip, vars[cutinds[i]], sol, boundswitch, usevbds ? 2 : 0, allowlocal, fixintegralrhs,
            ignoresol, boundsfortrans, boundtypesfortrans,
            bestlbs + i, bestubs + i, bestlbtypes + i, bestubtypes + i, selectedbounds + i, freevariable) );

      if( *freevariable )
         goto TERMINATE;
   }

   /* remember start of integer variables in the aggrrow */
   aggrrowintstart = i;

   /* perform bound substitution for continuous variables */
   for( i = 0; i < aggrrowintstart; ++i )
   {
      int v = cutinds[i];

      if( selectedbounds[i] == SCIP_BOUNDTYPE_LOWER )
      {
         assert(!SCIPisInfinity(scip, -bestlbs[i]));

         /* use lower bound as transformation bound: x'_j := x_j - lb_j */
         boundtype[i] = bestlbtypes[i];
         varsign[i] = +1;

         performBoundSubstitutionSafely(scip, cutcoefs, cutrhs, varsign[i], boundtype[i], bestlbs[i], v, localbdsused);
      }
      else
      {
         assert(!SCIPisInfinity(scip, bestubs[i]));

         /* use upper bound as transformation bound: x'_j := ub_j - x_j */
         boundtype[i] = bestubtypes[i];
         varsign[i] = -1;

         performBoundSubstitutionSafely(scip, cutcoefs, cutrhs, varsign[i], boundtype[i], bestubs[i], v, localbdsused);
      }

      if( SCIPisCertified(scip) )
      {
         assert(mirinfo != NULL);
         if( boundtype[i] == -2 )
         {
            mirinfo->localbdused[v] = TRUE;
            mirinfo->nlocalvars++;
         }
         mirinfo->upperused[v] = (varsign[i] == -1);
      }
   }

   /* remove integral variables that now have a zero coefficient due to variable bound usage of continuous variables
    * and determine the bound to use for the integer variables that are left
    */
   while( i < *nnz )
   {
      int v = cutinds[i];
      assert(cutinds[i] < firstcontvar);

      /* determine the best bounds for the integral variable, usevbd can be set to 0 here as vbds are only used for continuous variables */
      SCIP_CALL( determineBestBoundsSafely(scip, vars[v], sol, boundswitch, 0, allowlocal, fixintegralrhs,
            ignoresol, boundsfortrans, boundtypesfortrans,
            bestlbs + i, bestubs + i, bestlbtypes + i, bestubtypes + i, selectedbounds + i, freevariable) );

      /* increase i */
      ++i;

      if( *freevariable )
         goto TERMINATE;
   }

   /* now perform the bound substitution on the remaining integral variables which only uses standard bounds */
   for( i = aggrrowintstart; i < *nnz; ++i )
   {
      int v = cutinds[i];

      /* perform bound substitution */
      if( selectedbounds[i] == SCIP_BOUNDTYPE_LOWER )
      {
         assert(!SCIPisInfinity(scip, - bestlbs[i]));
         assert(bestlbtypes[i] < 0);

         /* use lower bound as transformation bound: x'_j := x_j - lb_j */
         boundtype[i] = bestlbtypes[i];
         varsign[i] = +1;

         performBoundSubstitutionSimpleSafely(scip, cutcoefs, cutrhs, boundtype[i], bestlbs[i], v, localbdsused);
      }
      else
      {
         assert(!SCIPisInfinity(scip, bestubs[i]));
         assert(bestubtypes[i] < 0);

         /* use upper bound as transformation bound: x'_j := ub_j - x_j */
         boundtype[i] = bestubtypes[i];
         varsign[i] = -1;

         performBoundSubstitutionSimpleSafely(scip, cutcoefs, cutrhs, boundtype[i], bestubs[i], v, localbdsused);
      }

      if( SCIPisCertified(scip) )
      {
         assert(mirinfo != NULL);
         if( boundtype[i] == -2 )
         {
            mirinfo->localbdused[v] = TRUE;
            mirinfo->nlocalvars++;
         }
         mirinfo->upperused[v] = (varsign[i] == -1);
      }
   }

  TERMINATE:
   SCIPintervalSetRoundingMode(previousroundmode); /*lint !e644*/

   /*free temporary memory */
   SCIPfreeBufferArray(scip, &selectedbounds);
   SCIPfreeBufferArray(scip, &bestubtypes);
   SCIPfreeBufferArray(scip, &bestlbtypes);
   SCIPfreeBufferArray(scip, &bestubs);
   SCIPfreeBufferArray(scip, &bestlbs);

   return SCIP_OKAY;
}

/** Transform equation \f$ a \cdot x = b; lb \leq x \leq ub \f$ into standard form
 *    \f$ a^\prime \cdot x^\prime = b,\; 0 \leq x^\prime \leq ub' \f$.
 *
 *  Transform variables (lb or ub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \mbox{if lb is used in transformation},\\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if ub is used in transformation},
 *  \end{array}
 *  \f]
 *  and move the constant terms \f$ a_j\, lb_j \f$ or \f$ a_j\, ub_j \f$ to the rhs.
 *
 *  Transform variables (vlb or vub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - (bl_j\, zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \mbox{if vlb is used in transf.} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if vub is used in transf.}
 *  \end{array}
 *  \f]
 *  move the constant terms \f$ a_j\, dl_j \f$ or \f$ a_j\, du_j \f$ to the rhs, and update the coefficient of the VLB variable:
 *  \f[
 *  \begin{array}{ll}
 *    a_{zl_j} := a_{zl_j} + a_j\, bl_j,& \mbox{or} \\
 *    a_{zu_j} := a_{zu_j} + a_j\, bu_j &
 *  \end{array}
 *  \f]
 */
static
SCIP_RETCODE cutsTransformMIR(
   SCIP*                 scip,               /**< SCIP datastructure */
   MIR_DATA*             data,               /**< the MIR data structure for this cut */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   SCIP_Bool             ignoresol,          /**< should the LP solution be ignored? (eg, apply MIR to dualray) */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in MIR row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_Real* bestlbs;
   SCIP_Real* bestubs;
   int* bestlbtypes;
   int* bestubtypes;
   SCIP_BOUNDTYPE* selectedbounds;
   SCIP_Real QUAD(coef);
   int totalnnz;
   int s;
   int i;

   assert(data != NULL);
   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);

   totalnnz = data->totalnnz;

   *freevariable = FALSE;
   *localbdsused = FALSE;

   int allocsize = MIN(NSECTIONS * totalnnz, data->nvars);
   /* allocate temporary memory to store best bounds and bound types */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbs, allocsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubs, allocsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbtypes, allocsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubtypes, allocsize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selectedbounds, allocsize) );

   /* transform the cut, one variable section at a time */
   for( s = 0; s < NSECTIONS; ++s )
   {
      int* indices = data->secindices[s];
      int cutindsstart = data->ncutinds;
      int usevbds = data->usevbds[s];

      i = 0;
      while( i < data->secnnz[s] )
      {
         int cutindex;
         int v = indices[i];

         /* due to variable bound usage, cancellation may have occurred */
         QUAD_ARRAY_LOAD(coef, data->cutcoefs, v);
         if( EPSZ(QUAD_TO_DBL(coef), QUAD_EPSILON) )
         {
            QUAD_ASSIGN(coef, 0.0);
            QUAD_ARRAY_STORE(data->cutcoefs, v, coef);
            --data->secnnz[s];
            --data->totalnnz;
            indices[i] = indices[data->secnnz[s]];
            /* do not increase the index */
            continue;
         }

         cutindex = data->ncutinds;
         assert(cutindex < allocsize);
         SCIP_CALL( determineBestBounds(scip, data->vars[v], sol, data, boundswitch, usevbds, allowlocal, fixintegralrhs,
               ignoresol, boundsfortrans, boundtypesfortrans,
               bestlbs + cutindex, bestubs + cutindex, bestlbtypes + cutindex, bestubtypes + cutindex,
               selectedbounds + cutindex, freevariable) );

         data->cutinds[cutindex] = v;
         ++data->ncutinds;

         ++i;

         /* if there is a free variable, we terminate because we cannot derive a MIR cut */
         if( *freevariable )
         {
            int j;
            int k;

            data->ncutinds = 0;

            /* if we terminate early, we need to make sure all the zeros in the cut coefficient array are cancelled */
            for( j = 0; j < NSECTIONS; ++j )
            {
               int* indexlist = data->secindices[j];
               for( k = 0; k < data->secnnz[j]; ++k )
               {
                  data->cutinds[data->ncutinds] = indexlist[k];
                  ++data->ncutinds;
               }
            }
            goto TERMINATE;
         }
      }

      /* perform bound substitution for added variables */
      for( i = cutindsstart; i < data->ncutinds; ++i )
      {
         SCIP_Real bestbnd;
         int v = data->cutinds[i];

         if( selectedbounds[i] == SCIP_BOUNDTYPE_LOWER )
         {
            assert(!SCIPisInfinity(scip, -bestlbs[i]));

            /* use lower bound as transformation bound: x'_j := x_j - lb_j */
            boundtype[i] = bestlbtypes[i];
            varsign[i] = +1;
            bestbnd = bestlbs[i];
         }
         else
         {
            assert(!SCIPisInfinity(scip, bestubs[i]));

            /* use upper bound as transformation bound: x'_j := ub_j - x_j */
            boundtype[i] = bestubtypes[i];
            varsign[i] = -1;
            bestbnd = bestubs[i];
         }
         doMIRBoundSubstitution(scip, data, varsign[i], boundtype[i], bestbnd, v, localbdsused);
      }
   }

   if( fixintegralrhs )
   {
      SCIP_Real f0;

      /* check if rhs is fractional */
      f0 = EPSFRAC(QUAD_TO_DBL(data->cutrhs), SCIPsumepsilon(scip));
      if( f0 < minfrac || f0 > maxfrac )
      {
         SCIP_Real bestviolgain;
         SCIP_Real bestnewf0;
         int besti;

         /* choose complementation of one variable differently such that f0 is in correct range */
         besti = -1;
         bestviolgain = -1e+100;
         bestnewf0 = 1.0;
         for( i = 0; i < data->ncutinds; i++ )
         {
            int v;

            v = data->cutinds[i];
            assert(0 <= v && v < data->nvars);

            QUAD_ARRAY_LOAD(coef, data->cutcoefs, v);
            assert(!EPSZ(QUAD_TO_DBL(coef), QUAD_EPSILON));

            if( boundtype[i] < 0
               && ((varsign[i] == +1 && !SCIPisInfinity(scip, bestubs[i]) && bestubtypes[i] < 0)
                  || (varsign[i] == -1 && !SCIPisInfinity(scip, -bestlbs[i]) && bestlbtypes[i] < 0)) )
            {
               SCIP_Real fj;
               SCIP_Real newfj;
               SCIP_Real newrhs;
               SCIP_Real newf0;
               SCIP_Real solval;
               SCIP_Real viol;
               SCIP_Real newviol;
               SCIP_Real violgain;

               /* currently:              a'_j =  varsign * a_j  ->  f'_j =  a'_j - floor(a'_j)
                * after complementation: a''_j = -varsign * a_j  -> f''_j = a''_j - floor(a''_j) = 1 - f'_j
                *                        rhs'' = rhs' + varsign * a_j * (lb_j - ub_j)
                * cut violation from f0 and fj:   f'_0 -  f'_j *  x'_j
                * after complementation:         f''_0 - f''_j * x''_j
                *
                * for continuous variables, we just set f'_j = f''_j = |a'_j|
                */
               newrhs = QUAD_TO_DBL(data->cutrhs) + varsign[i] * QUAD_TO_DBL(coef) * (bestlbs[i] - bestubs[i]);
               newf0 = EPSFRAC(newrhs, SCIPsumepsilon(scip));

               if( newf0 < minfrac || newf0 > maxfrac )
                  continue;
               if( v >= data->nvars - data->ncontvars )
               {
                  fj = REALABS(QUAD_TO_DBL(coef));
                  newfj = fj;
               }
               else
               {
                  fj = SCIPfrac(scip, varsign[i] * QUAD_TO_DBL(coef));
                  newfj = SCIPfrac(scip, -varsign[i] * QUAD_TO_DBL(coef));
               }

               if( !ignoresol )
               {
                  solval = SCIPgetSolVal(scip, sol, data->vars[v]);
                  viol = f0 - fj * (varsign[i] == +1 ? solval - bestlbs[i] : bestubs[i] - solval);
                  newviol = newf0 - newfj * (varsign[i] == -1 ? solval - bestlbs[i] : bestubs[i] - solval);
                  violgain = newviol - viol;
               }
               else
               {
                  /* todo: this should be done, this can improve the dualray significantly */
                  SCIPerrorMessage("Cannot handle closest bounds with ignoring the LP solution.\n");
                  return SCIP_INVALIDCALL;
               }

               /* prefer larger violations; for equal violations, prefer smaller f0 values since then the possibility that
                * we f_j > f_0 is larger and we may improve some coefficients in rounding
                */
               if( SCIPisGT(scip, violgain, bestviolgain) || (SCIPisGE(scip, violgain, bestviolgain) && newf0 < bestnewf0) )
               {
                  besti = i;
                  bestviolgain = violgain;
                  bestnewf0 = newf0;
               }
            }
         }

         if( besti >= 0 )
         {
            SCIP_Real QUAD(tmp);

            assert(besti < data->ncutinds);
            assert(boundtype[besti] < 0);
            assert(!SCIPisInfinity(scip, -bestlbs[besti]));
            assert(!SCIPisInfinity(scip, bestubs[besti]));

            QUAD_ARRAY_LOAD(coef, data->cutcoefs, data->cutinds[besti]);
            QUAD_SCALE(coef, varsign[besti]);

            /* switch the complementation of this variable */
            SCIPquadprecSumDD(tmp, bestlbs[besti], - bestubs[besti]);
            SCIPquadprecProdQQ(tmp, tmp, coef);
            SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);

            if( varsign[besti] == +1 )
            {
               /* switch to upper bound */
               assert(bestubtypes[besti] < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
               boundtype[besti] = bestubtypes[besti];
               varsign[besti] = -1;
            }
            else
            {
               /* switch to lower bound */
               assert(bestlbtypes[besti] < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
               boundtype[besti] = bestlbtypes[besti];
               varsign[besti] = +1;
            }
            *localbdsused = *localbdsused || (boundtype[besti] == -2);
         }
      }
   }

   TERMINATE:

   /*free temporary memory */
   SCIPfreeBufferArray(scip, &selectedbounds);
   SCIPfreeBufferArray(scip, &bestubtypes);
   SCIPfreeBufferArray(scip, &bestlbtypes);
   SCIPfreeBufferArray(scip, &bestubs);
   SCIPfreeBufferArray(scip, &bestlbs);

   return SCIP_OKAY;
}

/** Calculate fractionalities \f$ f_0 := b - down(b), f_j := a^\prime_j - down(a^\prime_j) \f$, and derive MIR cut \f$ \tilde{a} \cdot x' \leq down(b) \f$
 * \f[
 * \begin{array}{rll}
 *  integers :&  \tilde{a}_j = down(a^\prime_j),                        & if \qquad f_j \leq f_0 \\
 *            &  \tilde{a}_j = down(a^\prime_j) + (f_j - f_0)/(1 - f_0),& if \qquad f_j >  f_0 \\
 *  continuous:& \tilde{a}_j = 0,                                       & if \qquad a^\prime_j \geq 0 \\
 *             & \tilde{a}_j = a^\prime_j/(1 - f_0),                    & if \qquad a^\prime_j <  0
 * \end{array}
 * \f]
 *
 *  Transform inequality back to \f$ \hat{a} \cdot x \leq rhs \f$:
 *
 *  (lb or ub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{if lb was used in transformation}, \\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{if ub was used in transformation},
 * \end{array}
 * \f]
 *  and move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j \cdot lb_j = -\hat{a}_j \cdot lb_j,& \mbox{or} \\
 *     \tilde{a}_j \cdot ub_j = -\hat{a}_j \cdot ub_j &
 * \end{array}
 * \f]
 *  to the rhs.
 *
 *  (vlb or vub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - (bl_j \cdot zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{(vlb)} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{(vub)}
 * \end{array}
 * \f]
 *  move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j\, dl_j = -\hat{a}_j\, dl_j,& \mbox{or} \\
 *     \tilde{a}_j\, du_j = -\hat{a}_j\, du_j &
 * \end{array}
 * \f]
 *  to the rhs, and update the VB variable coefficients:
 * \f[
 * \begin{array}{ll}
 *    \hat{a}_{zl_j} := \hat{a}_{zl_j} - \tilde{a}_j\, bl_j = \hat{a}_{zl_j} - \hat{a}_j\, bl_j,& \mbox{or} \\
 *    \hat{a}_{zu_j} := \hat{a}_{zu_j} + \tilde{a}_j\, bu_j = \hat{a}_{zu_j} - \hat{a}_j\, bu_j &
 * \end{array}
 * \f]
 *
 *  @note this method is safe for usage in exact solving mode
 */
static
SCIP_RETCODE cutsRoundMIRSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*RESTRICT    cutcoefs,           /**< array of coefficients of cut */
   SCIP_Real*RESTRICT    cutrhs,             /**< pointer to right hand side of cut */
   int*RESTRICT          cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*RESTRICT          nnz,                /**< number of non-zeros in cut */
   int*RESTRICT          varsign,            /**< stores the sign of the transformed variable in summation */
   int*RESTRICT          boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub) */
   SCIP_INTERVAL         f0                  /**< fractional value of rhs */
   )
{
   SCIP_INTERVAL onedivoneminusf0;
   int i;
   int firstcontvar;
   SCIP_VAR** vars;
   int ndelcontvars;
   SCIP_ROUNDMODE previousroundmode;
   SCIP_MIRINFO* mirinfo = NULL;
   SCIP_INTERVAL tmpinterval;

   assert(cutrhs != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(boundtype != NULL);
   assert(varsign != NULL);
   assert(0.0 < SCIPintervalGetInf(f0) && SCIPintervalGetSup(f0) < 1.0);
   assert(SCIPisExact(scip));

   /* round up at first, since we are dividing and divisor should be as large as possible,
    * then switch to down since we are working on lhs */
   /* we need to careate the split-data for certification here, since part of the f_j > f_0 variables goes into the continuous part of the split */
   if( SCIPisCertified(scip) )
      mirinfo = SCIPgetCertificate(scip)->mirinfo[SCIPgetCertificate(scip)->nmirinfos - 1];

   previousroundmode = SCIPintervalGetRoundingMode();
   tmpinterval = f0;
   SCIPintervalSetBounds(&tmpinterval, -SCIPintervalGetSup(f0), -SCIPintervalGetInf(f0));
   SCIPintervalAddScalar(SCIPinfinity(scip), &tmpinterval, tmpinterval, 1.0);
   SCIPintervalSet(&onedivoneminusf0, 1.0);
   SCIPintervalDiv(SCIPinfinity(scip), &onedivoneminusf0, onedivoneminusf0, tmpinterval);
   SCIPintervalSetRoundingModeDownwards();

   /* Loop backwards to process integral variables first and be able to delete coefficients of integral variables
    * without destroying the ordering of the aggrrow's non-zeros.
    * (due to sorting in cutsTransformMIR the ordering is continuous before integral)
    */

   firstcontvar = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   vars = SCIPgetVars(scip);
#ifndef NDEBUG
   /*in debug mode check that all continuous variables of the aggrrow come before the integral variables */
   i = 0;
   while( i < *nnz && cutinds[i] >= firstcontvar )
      ++i;

   while( i < *nnz )
   {
      assert(cutinds[i] < firstcontvar);
      ++i;
   }
#endif

   /* round down everything on lhs (excepts for the denominator part above) */
   SCIPintervalSetRoundingModeDownwards();

   for( i = *nnz - 1; i >= 0 && cutinds[i] < firstcontvar; --i )
   {
      SCIP_VAR* var;
      SCIP_INTERVAL cutaj;

      int v;

      v = cutinds[i];
      assert(0 <= v && v < SCIPgetNVars(scip));

      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[i] == +1 || varsign[i] == -1);

      /* work on lhs -> round down */
      SCIPintervalSetRoundingModeDownwards();

      /* calculate the coefficient in the retransformed cut */
      {
         SCIP_Real aj;
         SCIP_Real downaj;
         SCIP_Real fj;

         aj = cutcoefs[v] * varsign[i];

         downaj = floor(aj);
         fj = aj - downaj;
         assert(fj >= -SCIPepsilon(scip) && fj <= 1.0);

         if( SCIPisLE(scip, fj, SCIPintervalGetInf(f0)) )
         {
            SCIPintervalSet(&cutaj, downaj);

            if( SCIPisCertified(scip) && mirinfo != NULL )
            {
               SCIP_RATIONAL* boundval;

               mirinfo->splitcoefficients[v] = SCIPintervalGetInf(cutaj); /*lint !e644*/
               assert(!SCIPisInfinity(scip, fabs(cutaj.inf)));
               if( mirinfo->upperused[v] )
               {
                  mirinfo->splitcoefficients[v] *= -1;
                  boundval = mirinfo->localbdused[v] ? SCIPvarGetUbLocalExact(var) : SCIPvarGetUbGlobalExact(var);
               }
               else
               {
                  boundval = mirinfo->localbdused[v] ? SCIPvarGetLbLocalExact(var) : SCIPvarGetLbGlobalExact(var);
               }
               SCIPrationalAddProdReal(mirinfo->rhs, boundval, mirinfo->splitcoefficients[v]);
            }
         }
         else
         {
            SCIPintervalSet(&tmpinterval, aj);
            SCIPintervalSubScalar(SCIPinfinity(scip), &tmpinterval, tmpinterval, downaj);
            SCIPintervalSub(SCIPinfinity(scip), &tmpinterval, tmpinterval, f0);
            SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, onedivoneminusf0);
            SCIPintervalAddScalar(SCIPinfinity(scip), &cutaj, tmpinterval, downaj);

            if( SCIPisCertified(scip) && mirinfo != NULL )
            {
               SCIP_RATIONAL* boundval;
               mirinfo->splitcoefficients[v] = downaj;
               mirinfo->splitcoefficients[v] += 1.0;
               if( mirinfo->upperused[v] )
               {
                  mirinfo->splitcoefficients[v] *= -1;
                  boundval = mirinfo->localbdused[v] ? SCIPvarGetUbLocalExact(var) : SCIPvarGetUbGlobalExact(var);
               }
               else
               {
                  boundval = mirinfo->localbdused[v] ? SCIPvarGetLbLocalExact(var) : SCIPvarGetLbGlobalExact(var);
               }
               SCIPrationalAddProdReal(mirinfo->rhs, boundval, mirinfo->splitcoefficients[v]);
            }
         }

         SCIPintervalMulScalar(SCIPinfinity(scip), &cutaj, cutaj, (double) varsign[i]);
      }

      /* integral var uses standard bound */
      assert(boundtype[i] < 0);

      if( cutaj.inf != 0.0 || cutaj.sup != 0 )
      {
         /* we have to use the inf of the cutaj-interval both times! */
         SCIPintervalSetRoundingModeUpwards();

         /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
         if( varsign[i] == +1 )
         {
            /* lower bound was used */
            if( boundtype[i] == -1 )
            {
               assert(SCIPrationalIsEQReal(SCIPvarGetLbGlobalExact(var), SCIPvarGetLbGlobal(var)));
               assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
               SCIPintervalSetRational(&tmpinterval, SCIPvarGetLbGlobalExact(var));
               SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, cutaj);
               *cutrhs += SCIPintervalGetSup(tmpinterval);
            }
            else
            {
               assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
               SCIPintervalSetRational(&tmpinterval, SCIPvarGetLbLocalExact(var));
               SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, cutaj);
               *cutrhs += SCIPintervalGetSup(tmpinterval);
            }
         }
         else
         {
            /* upper bound was used */
            if( boundtype[i] == -1 )
            {
               assert(SCIPrationalIsEQReal(SCIPvarGetUbGlobalExact(var), SCIPvarGetUbGlobal(var)));
               assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
               SCIPintervalSetRational(&tmpinterval, SCIPvarGetUbGlobalExact(var));
               SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, cutaj);
               *cutrhs += SCIPintervalGetSup(tmpinterval);
            }
            else
            {
               assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
               SCIPintervalSetRational(&tmpinterval, SCIPvarGetUbLocalExact(var));
               SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, cutaj);
               *cutrhs += SCIPintervalGetSup(tmpinterval);
            }
         }
      }

      /* remove zero cut coefficients from cut, only remove exactly 0 in exact solving mode
       * we can only do this here, since the sup might be positive and impact the rhs of the cut */
      if( cutaj.inf == 0.0 )
      {
         cutcoefs[v] = 0.0;
         --*nnz;
         cutinds[i] = cutinds[*nnz];
         continue;
      }

      cutcoefs[v] = SCIPintervalGetInf(cutaj);
   }

   /* adapt lhs -> round down */
   SCIPintervalSetRoundingModeDownwards();

   /* now process the continuous variables; postpone deletetion of zeros till all continuous variables have been processed */
   ndelcontvars = 0;
   while( i >= ndelcontvars )
   {
      SCIP_VAR* var;
      SCIP_INTERVAL cutaj;
      SCIP_Real aj;
      int v;

      v = cutinds[i];
      assert(0 <= v && v < SCIPgetNVars(scip));

      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[i] == +1 || varsign[i] == -1);
      assert( v >= firstcontvar );

      /* adapt lhs -> round down */
      SCIPintervalSetRoundingModeDownwards();

      /* calculate the coefficient in the retransformed cut */
      aj = cutcoefs[v];

      if( aj * varsign[i] >= 0.0 )
         SCIPintervalSet(&cutaj, 0.0);
      else
      {
         SCIPintervalMulScalar(SCIPinfinity(scip), &cutaj, onedivoneminusf0, aj); /* cutaj = varsign[i] * aj * onedivoneminusf0; // a^_j */
      }

      /* remove zero cut coefficients from cut; move a continuous var from the beginning
       * to the current position, so that all integral variables stay behind the continuous
       * variables
       */
      if( EPSZ(SCIPintervalGetInf(cutaj), QUAD_EPSILON) && (SCIPintervalGetInf(cutaj) >= 0.0) )
      {
         SCIPintervalSet(&cutaj, 0.0);
         cutcoefs[v] = 0.0;
         cutinds[i] = cutinds[ndelcontvars];
         varsign[i] = varsign[ndelcontvars];
         boundtype[i] = boundtype[ndelcontvars];
         ++ndelcontvars;
         continue;
      }

      cutcoefs[v] = SCIPintervalGetInf(cutaj);

      SCIPintervalSetRoundingModeUpwards();

      /* check for variable bound use */
      if( boundtype[i] < 0 )
      {
         /* standard bound */

         /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
         if( varsign[i] == +1 )
         {
            /* lower bound was used */
            if( boundtype[i] == -1 )
            {
               assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
               SCIPintervalSetRational(&tmpinterval, SCIPvarGetLbGlobalExact(var));
               SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, cutaj);
               *cutrhs += SCIPintervalGetSup(tmpinterval);
            }
            else
            {
               assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
               SCIPintervalSetRational(&tmpinterval, SCIPvarGetLbLocalExact(var));
               SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, cutaj);
               *cutrhs += SCIPintervalGetSup(tmpinterval);
            }
         }
         else
         {
            /* upper bound was used */
            if( boundtype[i] == -1 )
            {
               assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
               SCIPintervalSetRational(&tmpinterval, SCIPvarGetUbGlobalExact(var));
               SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, cutaj);
               *cutrhs += SCIPintervalGetSup(tmpinterval);
            }
            else
            {
               assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
               SCIPintervalSetRational(&tmpinterval, SCIPvarGetUbLocalExact(var));
               SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, tmpinterval, cutaj);
               *cutrhs += SCIPintervalGetSup(tmpinterval);
            }
         }
      }
      else
      {
         SCIPerrorMessage("varbounds not yet implemented in exact SCIP \n");
         return SCIP_ERROR;
      }

      /* advance to next variable */
      --i;
   }

   /* fill the empty position due to deleted continuous variables */
   if( ndelcontvars > 0 )
   {
      assert(ndelcontvars <= *nnz);
      *nnz -= ndelcontvars;
      if( *nnz < ndelcontvars )
      {
         BMScopyMemoryArray(cutinds, cutinds + ndelcontvars, *nnz);
      }
      else
      {
         BMScopyMemoryArray(cutinds, cutinds + *nnz, ndelcontvars);
      }
   }

   /* reset rounding mode, also set the rhs->data in the mirinfo */
   SCIPintervalSetRoundingMode(previousroundmode);

   return SCIP_OKAY;
}

#ifdef SCIP_DISABLED_CODE
/** Calculate fractionalities \f$ f_0 := b - down(b), f_j := a^\prime_j - down(a^\prime_j) \f$, and derive MIR cut \f$ \tilde{a} \cdot x' \leq down(b) \f$
 * \f[
 * \begin{array}{rll}
 *  integers :&  \tilde{a}_j = down(a^\prime_j),                        & if \qquad f_j \leq f_0 \\
 *            &  \tilde{a}_j = down(a^\prime_j) + (f_j - f_0)/(1 - f_0),& if \qquad f_j >  f_0 \\
 *  continuous:& \tilde{a}_j = 0,                                       & if \qquad a^\prime_j \geq 0 \\
 *             & \tilde{a}_j = a^\prime_j/(1 - f_0),                    & if \qquad a^\prime_j <  0
 * \end{array}
 * \f]
 *
 *  Transform inequality back to \f$ \hat{a} \cdot x \leq rhs \f$:
 *
 *  (lb or ub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{if lb was used in transformation} \\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{if ub was used in transformation}
 * \end{array}
 * \f]
 *  and move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j \cdot lb_j = -\hat{a}_j \cdot lb_j,& \mbox{or} \\
 *     \tilde{a}_j \cdot ub_j = -\hat{a}_j \cdot ub_j &
 * \end{array}
 * \f]
 *  to the rhs.
 *
 *  (vlb or vub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - (bl_j \cdot zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{(vlb)} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{(vub)}
 * \end{array}
 * \f]
 *  move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j\, dl_j = -\hat{a}_j\, dl_j,& \mbox{or} \\
 *     \tilde{a}_j\, du_j = -\hat{a}_j\, du_j &
 * \end{array}
 * \f]
 *  to the rhs, and update the VB variable coefficients:
 * \f[
 * \begin{array}{ll}
 *    \hat{a}_{zl_j} := \hat{a}_{zl_j} - \tilde{a}_j\, bl_j = \hat{a}_{zl_j} - \hat{a}_j\, bl_j,& \mbox{or} \\
 *    \hat{a}_{zu_j} := \hat{a}_{zu_j} + \tilde{a}_j\, bu_j = \hat{a}_{zu_j} - \hat{a}_j\, bu_j &
 * \end{array}
 * \f]
 *
 *  @note this method is safe for usage in exact solving mode
 */
static
SCIP_RETCODE cutsRoundMIRRational(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*RESTRICT    cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*RESTRICT cutrhs),          /**< pointer to right hand side of cut */
   int*RESTRICT          cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*RESTRICT          nnz,                /**< number of non-zeros in cut */
   int*RESTRICT          varsign,            /**< stores the sign of the transformed variable in summation */
   int*RESTRICT          boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub) */
   SCIP_RATIONAL*        f0                  /**< fractional value of rhs */
   )
{
   SCIP_RATIONAL* tmp;
   SCIP_RATIONAL* onedivoneminusf0;
   int i;
   int firstcontvar;
   SCIP_VAR** vars;
   int ndelcontvars;
   SCIP_ROUNDMODE previousroundmode;
   SCIP_MIRINFO* mirinfo;

   assert(QUAD_HI(cutrhs) != NULL);
   assert(cutcoefs != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(boundtype != NULL);
   assert(varsign != NULL);
   assert(SCIPrationalIsPositive(f0) && SCIPrationalIsLTReal(f0, 1.0));
   assert(SCIPisExact(scip));

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &onedivoneminusf0) );
   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmp) );

   /* round up at first, since we are dividing and divisor should be as large as possible,
    * then switch to down since we are working on lhs */
   /* we need to careate the split-data for certification here, since part of the f_j > f_0 variables goes into the continuous part of the split */
   if( SCIPisCertified(scip) )
      mirinfo = SCIPgetCertificate(scip)->mirinfo[SCIPgetCertificate(scip)->nmirinfos - 1];

   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPrationalSetReal(tmp, 1.0);
   SCIPrationalDiff(tmp, tmp, f0);
   SCIPrationalSetReal(onedivoneminusf0, 1.0);
   SCIPrationalDiv(onedivoneminusf0, onedivoneminusf0, tmp);

   /* Loop backwards to process integral variables first and be able to delete coefficients of integral variables
    * without destroying the ordering of the aggrrow's non-zeros.
    * (due to sorting in cutsTransformMIR the ordering is continuous before integral)
    */

   firstcontvar = SCIPgetNVars(scip) - SCIPgetNContVars(scip);
   vars = SCIPgetVars(scip);
#ifndef NDEBUG
   /*in debug mode check that all continuous variables of the aggrrow come before the integral variables */
   i = 0;
   while( i < *nnz && cutinds[i] >= firstcontvar )
      ++i;

   while( i < *nnz )
   {
      assert(cutinds[i] < firstcontvar);
      ++i;
   }
#endif

   for( i = *nnz - 1; i >= 0 && cutinds[i] < firstcontvar; --i )
   {
      SCIP_VAR* var;
      SCIP_RATIONAL* cutaj;
      SCIP_Real QUAD(cutajquad);
      int v;

      v = cutinds[i];
      assert(0 <= v && v < SCIPgetNVars(scip));

      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &cutaj) );

      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[i] == +1 || varsign[i] == -1);

      /* calculate the coefficient in the retransformed cut */
      {
         SCIP_Real QUAD(aj);
         SCIP_Real downaj;
         SCIP_RATIONAL* fj;

         SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &fj) );

         QUAD_ARRAY_LOAD(aj, cutcoefs, v);
         QUAD_SCALE(aj, varsign[i]);
         SCIPrationalSetReal(tmp, aj);

         downaj = floor(QUAD_TO_DBL(aj));
         SCIPrationalDiffReal(fj, tmp, downaj);

         if( SCIPrationalIsLE(fj, f0) )
         {
            SCIPrationalSetReal(cutaj, downaj);

            if( SCIPisCertified(scip) )
            {
               SCIP_RATIONAL* boundval;

               mirinfo->splitcoefficients[v] = downaj;
               if( mirinfo->upperused[v] )
               {
                  mirinfo->splitcoefficients[v] *= -1;
                  boundval = mirinfo->localbdused[v] ? SCIPvarGetUbLocalExact(var) : SCIPvarGetUbGlobalExact(var);
               }
               else
               {
                  boundval = mirinfo->localbdused[v] ? SCIPvarGetLbLocalExact(var) : SCIPvarGetLbGlobalExact(var);
               }
               SCIPrationalAddProdReal(mirinfo->rhs, boundval, mirinfo->splitcoefficients[v]);
            }
         }
         else
         {
            SCIPrationalSetReal(tmp, QUAD_TO_DBL(aj));
            SCIPrationalDiffReal(tmp, tmp, downaj);
            SCIPrationalDiff(tmp, tmp, f0);
            SCIPrationalMult(tmp, tmp, onedivoneminusf0);
            SCIPrationalAddReal(cutaj, tmp, downaj);

            if( SCIPisCertified(scip) )
            {
               SCIP_RATIONAL* boundval;

               mirinfo->splitcoefficients[v] = QUAD_TO_DBL(downaj);
               mirinfo->splitcoefficients[v] += 1.0;
               if( mirinfo->upperused[v] )
               {
                  mirinfo->splitcoefficients[v] *= -1;
                  boundval = mirinfo->localbdused[v] ? SCIPvarGetUbLocalExact(var) : SCIPvarGetUbGlobalExact(var);
               }
               else
               {
                  boundval = mirinfo->localbdused[v] ? SCIPvarGetLbLocalExact(var) : SCIPvarGetLbGlobalExact(var);
               }
               SCIPrationalAddProdReal(mirinfo->rhs, boundval, mirinfo->splitcoefficients[v]);
            }
         }

         SCIPrationalMultReal(cutaj, cutaj, varsign[i]);

         SCIPrationalFreeBuffer(SCIPbuffer(scip), &fj);
      }

      /* remove zero cut coefficients from cut, only remove positive coefficients in exact solving mode */
      if( SCIPrationalIsZero(cutaj) )
      {
         QUAD_ASSIGN(cutajquad, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, cutajquad);
         --*nnz;
         cutinds[i] = cutinds[*nnz];
         SCIPrationalFreeBuffer(SCIPbuffer(scip), &cutaj);
         continue;
      }

      QUAD_ASSIGN(cutajquad, SCIPrationalRoundReal(cutaj, SCIP_R_ROUND_DOWNWARDS));

      QUAD_ARRAY_STORE(cutcoefs, v, cutajquad);

      /* integral var uses standard bound */
      assert(boundtype[i] < 0);

      SCIPintervalSetRoundingModeUpwards();

      /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
      if( varsign[i] == +1 )
      {
         /* lower bound was used */
         if( boundtype[i] == -1 )
         {
            assert(SCIPrationalRoundReal(SCIPvarGetLbGlobalExact(var), SCIP_R_ROUND_DOWNWARDS) ==  SCIPvarGetLbGlobal(var));
            assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
            SCIPrationalMult(tmp, cutaj, SCIPvarGetLbGlobalExact(var));
            SCIPquadprecSumQD(*cutrhs, *cutrhs, SCIPrationalRoundReal(tmp, SCIP_R_ROUND_UPWARDS)); /* rhs += cutaj * SCIPvarGetLbGlobal(var) */
         }
         else
         {
            assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
            SCIPrationalMult(tmp, cutaj, SCIPvarGetLbLocalExact(var));
            SCIPquadprecSumQD(*cutrhs, *cutrhs, SCIPrationalRoundReal(tmp, SCIP_R_ROUND_UPWARDS)); /* rhs += cutaj * SCIPvarGetLbLocal(var) */
         }
      }
      else
      {
         /* upper bound was used */
         if( boundtype[i] == -1 )
         {
            assert(SCIPrationalRoundReal(SCIPvarGetUbGlobalExact(var), SCIP_R_ROUND_UPWARDS) ==  SCIPvarGetUbGlobal(var));
            assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
            SCIPrationalMult(tmp, cutaj, SCIPvarGetUbGlobalExact(var));
            SCIPquadprecSumQD(*cutrhs, *cutrhs, SCIPrationalRoundReal(tmp, SCIP_R_ROUND_UPWARDS)); /* rhs += cutaj * SCIPvarGetUbGlobal(var) */
         }
         else
         {
            assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
            SCIPrationalMult(tmp, cutaj, SCIPvarGetUbLocalExact(var));
            SCIPquadprecSumQD(*cutrhs, *cutrhs, SCIPrationalRoundReal(tmp, SCIP_R_ROUND_UPWARDS)); /* rhs += cutaj * SCIPvarGetUbLocal(var) */
         }
      }
      SCIPrationalFreeBuffer(SCIPbuffer(scip), &cutaj);
   }

   /* now process the continuous variables; postpone deletetion of zeros till all continuous variables have been processed */
   ndelcontvars = 0;
   while( i >= ndelcontvars )
   {
      SCIP_VAR* var;
      SCIP_RATIONAL* cutaj;
      SCIP_RATIONAL* tmprational;
      SCIP_Real QUAD(cutajquad);
      SCIP_Real QUAD(aj);
      int v;

      /* adapt lhs -> round down */
      SCIPintervalSetRoundingModeDownwards();

      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &cutaj) );
      SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmprational) );

      v = cutinds[i];
      assert(0 <= v && v < SCIPgetNVars(scip));

      var = vars[v];
      assert(var != NULL);
      assert(SCIPvarGetProbindex(var) == v);
      assert(varsign[i] == +1 || varsign[i] == -1);
      assert( v >= firstcontvar );

      /* calculate the coefficient in the retransformed cut */
      QUAD_ARRAY_LOAD(aj, cutcoefs, v);

      if( QUAD_TO_DBL(aj) * varsign[i] >= 0.0 )
         SCIPrationalSetReal(cutaj, 0.0);
      else
      {
         SCIPrationalSetRational(cutaj, onedivoneminusf0);
         SCIPrationalMultReal(cutaj, cutaj, QUAD_TO_DBL(aj)); /* cutaj = varsign[i] * aj * onedivoneminusf0; // a^_j */
      }

      /* remove zero cut coefficients from cut; move a continuous var from the beginning
       * to the current position, so that all integral variables stay behind the continuous
       * variables
       */
      if( EPSZ(SCIPrationalGetReal(cutaj), QUAD_EPSILON) && SCIPrationalIsGEReal(cutaj, 0.0) )
      {
         assert(SCIPrationalIsZero(cutaj));
         SCIPrationalSetReal(cutaj, 0.0);
         QUAD_ASSIGN_Q(cutajquad, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, cutajquad);
         cutinds[i] = cutinds[ndelcontvars];
         varsign[i] = varsign[ndelcontvars];
         boundtype[i] = boundtype[ndelcontvars];
         ++ndelcontvars;

         SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmprational);
         SCIPrationalFreeBuffer(SCIPbuffer(scip), &cutaj);

         continue;
      }

      QUAD_ASSIGN(cutajquad, SCIPrationalRoundReal(cutaj, SCIP_R_ROUND_DOWNWARDS));
      QUAD_ARRAY_STORE(cutcoefs, v, cutajquad);

      SCIPintervalSetRoundingModeUpwards();

      /* check for variable bound use */
      if( boundtype[i] < 0 )
      {
         /* standard bound */

         /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
         if( varsign[i] == +1 )
         {
            /* lower bound was used */
            if( boundtype[i] == -1 )
            {
               assert(SCIPrationalRoundReal(SCIPvarGetLbGlobalExact(var), SCIP_R_ROUND_DOWNWARDS) ==  SCIPvarGetLbGlobal(var));
               assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
               SCIPrationalMult(tmprational, cutaj, SCIPvarGetLbGlobalExact(var));
               SCIPquadprecSumQD(*cutrhs, *cutrhs, SCIPrationalRoundReal(tmprational, SCIP_R_ROUND_UPWARDS));
            }
            else
            {
               assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
               SCIPrationalMult(tmprational, cutaj, SCIPvarGetLbLocalExact(var));
               SCIPquadprecSumQD(*cutrhs, *cutrhs, SCIPrationalRoundReal(tmprational, SCIP_R_ROUND_UPWARDS));
            }
         }
         else
         {
            /* upper bound was used */
            if( boundtype[i] == -1 )
            {
               assert(SCIPrationalRoundReal(SCIPvarGetUbGlobalExact(var), SCIP_R_ROUND_UPWARDS) ==  SCIPvarGetUbGlobal(var));
               assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
               SCIPrationalMult(tmprational, cutaj, SCIPvarGetUbGlobalExact(var));
               SCIPquadprecSumQD(*cutrhs, *cutrhs, SCIPrationalRoundReal(tmprational, SCIP_R_ROUND_UPWARDS));
            }
            else
            {
               assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
               SCIPrationalMult(tmprational, cutaj, SCIPvarGetUbLocalExact(var));
               SCIPquadprecSumQD(*cutrhs, *cutrhs, SCIPrationalRoundReal(tmprational, SCIP_R_ROUND_UPWARDS));
            }
         }

         SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmprational);
         SCIPrationalFreeBuffer(SCIPbuffer(scip), &cutaj);
      }
      else
      {
#ifdef SCIP_DISABLED_CODE
         SCIP_VAR** vbz;
         SCIP_Real* vbb;
         SCIP_Real* vbd;
         SCIP_Real QUAD(zcoef);
         int vbidx;
         int zidx;

         assert(!SCIPisExact(scip));

         /* variable bound */
         vbidx = boundtype[i];

         /* change mirrhs and cutaj of integer variable z_j of variable bound */
         if( varsign[i] == +1 )
         {
            /* variable lower bound was used */
            assert(0 <= vbidx && vbidx < SCIPvarGetNVlbs(var));
            vbz = SCIPvarGetVlbVars(var);
            vbb = SCIPvarGetVlbCoefs(var);
            vbd = SCIPvarGetVlbConstants(var);
         }
         else
         {
            /* variable upper bound was used */
            assert(0 <= vbidx && vbidx < SCIPvarGetNVubs(var));
            vbz = SCIPvarGetVubVars(var);
            vbb = SCIPvarGetVubCoefs(var);
            vbd = SCIPvarGetVubConstants(var);
         }
         assert(SCIPvarIsActive(vbz[vbidx]));
         zidx = SCIPvarGetProbindex(vbz[vbidx]);
         assert(0 <= zidx && zidx < firstcontvar);

         SCIPquadprecProdQD(tmp, cutaj, vbd[vbidx]);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);

         SCIPquadprecProdQD(tmp, cutaj, vbb[vbidx]);
         QUAD_ARRAY_LOAD(zcoef, cutcoefs, zidx);

         /* update sparsity pattern */
         if( QUAD_HI(zcoef) == 0.0 )
            cutinds[(*nnz)++] = zidx;

         SCIPquadprecSumQQ(zcoef, zcoef, -tmp);
         QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
         QUAD_ARRAY_STORE(cutcoefs, zidx, zcoef);
         assert(QUAD_HI(zcoef) != 0.0);
#endif
      }

      /* advance to next variable */
      --i;
   }

   /* fill the empty position due to deleted continuous variables */
   if( ndelcontvars > 0 )
   {
      assert(ndelcontvars <= *nnz);
      *nnz -= ndelcontvars;
      if( *nnz < ndelcontvars )
      {
         BMScopyMemoryArray(cutinds, cutinds + ndelcontvars, *nnz);
      }
      else
      {
         BMScopyMemoryArray(cutinds, cutinds + *nnz, ndelcontvars);
      }
   }

   /* reset rounding mode, also set the rhs->data in the mirinfo */
   SCIPintervalSetRoundingMode(previousroundmode);

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmp);
   SCIPrationalFreeBuffer(SCIPbuffer(scip), &onedivoneminusf0);

   return SCIP_OKAY;
}
#endif


/** Calculate fractionalities \f$ f_0 := b - down(b), f_j := a^\prime_j - down(a^\prime_j) \f$, and derive MIR cut \f$ \tilde{a} \cdot x' \leq down(b) \f$
 * \f[
 * \begin{array}{rll}
 *  integers :&  \tilde{a}_j = down(a^\prime_j),                        & if \qquad f_j \leq f_0 \\
 *            &  \tilde{a}_j = down(a^\prime_j) + (f_j - f_0)/(1 - f_0),& if \qquad f_j >  f_0 \\
 *  continuous:& \tilde{a}_j = 0,                                       & if \qquad a^\prime_j \geq 0 \\
 *             & \tilde{a}_j = a^\prime_j/(1 - f_0),                    & if \qquad a^\prime_j <  0
 * \end{array}
 * \f]
 *
 *  Transform inequality back to \f$ \hat{a} \cdot x \leq rhs \f$:
 *
 *  (lb or ub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{if lb was used in transformation} \\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{if ub was used in transformation}
 * \end{array}
 * \f]
 *  and move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j \cdot lb_j = -\hat{a}_j \cdot lb_j,& \mbox{or} \\
 *     \tilde{a}_j \cdot ub_j = -\hat{a}_j \cdot ub_j &
 * \end{array}
 * \f]
 *  to the rhs.
 *
 *  (vlb or vub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - (bl_j \cdot zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{(vlb)} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{(vub)}
 * \end{array}
 * \f]
 *  move the constant terms
 * \f[
 * \begin{array}{cl}
 *    -\tilde{a}_j\, dl_j = -\hat{a}_j\, dl_j,& \mbox{or} \\
 *     \tilde{a}_j\, du_j = -\hat{a}_j\, du_j &
 * \end{array}
 * \f]
 *  to the rhs, and update the VB variable coefficients:
 * \f[
 * \begin{array}{ll}
 *    \hat{a}_{zl_j} := \hat{a}_{zl_j} - \tilde{a}_j\, bl_j = \hat{a}_{zl_j} - \hat{a}_j\, bl_j,& \mbox{or} \\
 *    \hat{a}_{zu_j} := \hat{a}_{zu_j} + \tilde{a}_j\, bu_j = \hat{a}_{zu_j} - \hat{a}_j\, bu_j &
 * \end{array}
 * \f]
 */
static
SCIP_RETCODE cutsRoundMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   MIR_DATA*             data,               /**< the MIR data structure for this cut */
   int*RESTRICT          varsign,            /**< stores the sign of the transformed variable in summation */
   int*RESTRICT          boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub) */
   QUAD(SCIP_Real        f0)                 /**< fractional value of rhs */
   )
{
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(onedivoneminusf0);
   int s;
   int cutindex;
   int i;

   assert(data != NULL);
   assert(boundtype != NULL);
   assert(varsign != NULL);
   assert(0.0 < QUAD_TO_DBL(f0) && QUAD_TO_DBL(f0) < 1.0);

   SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0);
   SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0);

   /* Loop backwards through the sections, so that the reversing of varbound substitutions does not prematurely effect
    * the coefficients of variables in other sections, because the section index of a variable bound must always be
    * higher than that of the bounded variable. */
   cutindex = data->ncutinds - 1;
   for( s = NSECTIONS - 1; s >= 0; --s )
   {
      int* indices = data->secindices[s];
      int nnz = data->secnnz[s];

      SCIP_Bool enfintegral = data->isenfint[s];
      SCIP_Bool implintegral = data->isimplint[s];

      /* iterate backwards over indices in section, so we can easily shrink the section if we find zeros */
      for( i = nnz - 1; i >= 0 ; --i )
      {
         int v;
         int sign;
         int type;
         SCIP_Real QUAD(cutaj);
         SCIP_Real QUAD(aj);
         SCIP_VAR* var;

         v = indices[i];
         assert(0 <= v && v < data->nvars);
         assert(data->cutinds[cutindex] == v);

         sign = varsign[cutindex];
         assert(sign == +1 || sign == -1);
         type = boundtype[cutindex];

         --cutindex;

         var = data->vars[v];
         assert(var != NULL);
         assert(SCIPvarGetProbindex(var) == v );

         QUAD_ARRAY_LOAD(aj, data->cutcoefs, v);

         if( enfintegral || implintegral )
         {
            /* variable is integral */
            SCIP_Real QUAD(downaj);
            SCIP_Real QUAD(fj);

            QUAD_SCALE(aj, sign);

            SCIPquadprecEpsFloorQ(downaj, aj, SCIPepsilon(scip)); /*lint !e666*/
            SCIPquadprecSumQQ(fj, aj, -downaj);
            assert(QUAD_TO_DBL(fj) >= -SCIPepsilon(scip) && QUAD_TO_DBL(fj) < 1.0);

            if( SCIPisLE(scip, QUAD_TO_DBL(fj), QUAD_TO_DBL(f0)) )
            {
               QUAD_ASSIGN_Q(cutaj, downaj); /* a^_j */
            }
            else
            {
               SCIPquadprecSumQQ(tmp, fj, -f0);
               SCIPquadprecProdQQ(tmp, tmp, onedivoneminusf0);
               SCIPquadprecSumQQ(cutaj, tmp, downaj);
            }
            QUAD_SCALE(cutaj, sign);
         }
         else
         {
            /* variable is continuous */
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);

            if( QUAD_TO_DBL(aj) * sign >= 0.0 )
               QUAD_ASSIGN(cutaj, 0.0);
            else
               SCIPquadprecProdQQ(cutaj, onedivoneminusf0, aj); /* cutaj = aj * onedivoneminusf0 */
         }

         /* remove coefficient from cut if it becomes zero */
         if( EPSZ(QUAD_TO_DBL(cutaj), QUAD_EPSILON) )
         {
            QUAD_ASSIGN(cutaj, 0.0);
            QUAD_ARRAY_STORE(data->cutcoefs, v, cutaj);
            --data->totalnnz;
            --data->secnnz[s];
            indices[i] = indices[data->secnnz[s]];
            continue;
         }

         /* store the updated coefficient */
         QUAD_ARRAY_STORE(data->cutcoefs, v, cutaj);

         /* undo bound transformations */
         if( type < 0 )
         {
            /* standard bound */
            /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
            if( sign == +1 )
            {
               /* lower bound was used */
               if( type == -1 )
               {
                  assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
                  SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbGlobal(var));
                  SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
                  SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbLocal(var));
                  SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);
               }
            }
            else
            {
               /* upper bound was used */
               if( type == -1 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
                  SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbGlobal(var));
                  SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbLocal(var));
                  SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);
               }
            }
         }
         else
         {
            /* variable bound */
            SCIP_VAR** vbz;
            SCIP_Real* vbb;
            SCIP_Real* vbd;
            SCIP_Real QUAD(zcoef);
            int vbidx;
            int zidx;

            /* variable bound */
            vbidx = type;

            /* change mirrhs and cutaj of integer variable z_j of variable bound */
            if( sign == +1 )
            {
               /* variable lower bound was used */
               assert(0 <= vbidx && vbidx < SCIPvarGetNVlbs(var));
               vbz = SCIPvarGetVlbVars(var);
               vbb = SCIPvarGetVlbCoefs(var);
               vbd = SCIPvarGetVlbConstants(var);
            }
            else
            {
               /* variable upper bound was used */
               assert(0 <= vbidx && vbidx < SCIPvarGetNVubs(var));
               vbz = SCIPvarGetVubVars(var);
               vbb = SCIPvarGetVubCoefs(var);
               vbd = SCIPvarGetVubConstants(var);
            }
            assert(SCIPvarIsActive(vbz[vbidx]));
            zidx = SCIPvarGetProbindex(vbz[vbidx]);
            assert(varSection(data, zidx) > s);

            SCIPquadprecProdQD(tmp, cutaj, vbd[vbidx]);
            SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);

            SCIPquadprecProdQD(tmp, cutaj, vbb[vbidx]);
            QUAD_ARRAY_LOAD(zcoef, data->cutcoefs, zidx);

            /* update sparsity pattern */
            if( QUAD_HI(zcoef) == 0.0 )
            {
               int zsection = varSection(data, zidx);
               data->secindices[zsection][data->secnnz[zsection]] = zidx;
               ++data->secnnz[zsection];
               ++data->totalnnz;
            }

            SCIPquadprecSumQQ(zcoef, zcoef, -tmp);
            QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
            QUAD_ARRAY_STORE(data->cutcoefs, zidx, zcoef);
            assert(QUAD_HI(zcoef) != 0.0);
         }
      }
   }

   /* Finally, store the relevant data in cutinds which is the array used by the other functions */
   data->ncutinds = 0;
   for( s = 0; s < NSECTIONS; ++s )
   {
      int* indices = data->secindices[s];
      int nnz = data->secnnz[s];
      for( i = 0; i < nnz; ++i )
      {
         data->cutinds[data->ncutinds] = indices[i];
         ++data->ncutinds;
      }
   }
   return SCIP_OKAY;
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row: \f$ a^\prime_r = scale \cdot weight[r] \cdot slacksign[r]. \f$
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 *  \f[
 *  \begin{array}{rll}
 *    integers : & \hat{a}_r = \tilde{a}_r = down(a^\prime_r),                        & \mbox{if}\qquad f_r \leq f_0 \\
 *               & \hat{a}_r = \tilde{a}_r = down(a^\prime_r) + (f_r - f_0)/(1 - f_0),& \mbox{if}\qquad f_r >  f_0 \\
 *    continuous:& \hat{a}_r = \tilde{a}_r = 0,                                       & \mbox{if}\qquad a^\prime_r \geq 0 \\
 *               & \hat{a}_r = \tilde{a}_r = a^\prime_r/(1 - f_0),                    & \mbox{if}\qquad a^\prime_r <  0
 *  \end{array}
 *  \f]
 *
 *  Substitute \f$ \hat{a}_r \cdot s_r \f$ by adding \f$ \hat{a}_r \f$ times the slack's definition to the cut.
 *
 *  @note this method is safe for usage in exact solving mode
 *
 *  @todo certify and use integrality of row in exact solving mode
 *
 *  @todo make behavior identical to the unsafe MIR cut computation
 */
static
SCIP_RETCODE cutsSubstituteMIRSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   int*                  rowinds,            /**< sparsity pattern of used rows */
   int                   nrowinds,           /**< number of used rows */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   SCIP_Real*            cutrhs,             /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   SCIP_INTERVAL         f0                  /**< fractional value of rhs */
   )
{  /*lint --e{715}*/
   SCIP_ROW** rows;
   SCIP_ROW* userow;
   SCIP_ROWEXACT* rowexact;
   SCIP_INTERVAL onedivoneminusf0, tmpinterval;
   SCIP_ROUNDMODE previousroundmode;
   SCIP_AGGREGATIONINFO* aggrinfo = NULL;
   SCIP_MIRINFO* mirinfo = NULL;
   SCIP_Real mult;
   SCIP_Real splitcoef;
   SCIP_Real slackweight;
   SCIP_Bool slackroundeddown;
   int i;
   int currentnegslackrow;

   assert(scip != NULL);
   assert(weights != NULL || nrowinds == 0);
   assert(slacksign != NULL || nrowinds == 0);
   assert(rowinds != NULL || nrowinds == 0);
   assert(scale > 0.0);
   assert(cutcoefs != NULL);
   assert(cutrhs != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(0.0 < SCIPintervalGetInf(f0) && SCIPintervalGetSup(f0) < 1.0);

   assert(SCIPisExact(scip));

   /* compute 1/(1-f0) in interval arithmetic */
   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetBounds(&tmpinterval, -SCIPintervalGetSup(f0), -SCIPintervalGetInf(f0));
   SCIPintervalAddScalar(SCIPinfinity(scip), &tmpinterval, tmpinterval, 1.0);
   SCIPintervalSet(&onedivoneminusf0, 1.0);
   SCIPintervalDiv(SCIPinfinity(scip), &onedivoneminusf0, onedivoneminusf0, tmpinterval);

   if( SCIPisCertified(scip)   )
   {
      aggrinfo = SCIPgetCertificate(scip)->aggrinfo[SCIPgetCertificate(scip)->naggrinfos -1];
      mirinfo = SCIPgetCertificate(scip)->mirinfo[SCIPgetCertificate(scip)->nmirinfos - 1];
   }

   rows = SCIPgetLPRows(scip);
   currentnegslackrow = 0;
   for( i = 0; i < nrowinds; i++ )
   {
      SCIP_ROW* row;
      SCIP_INTERVAL ar;
      SCIP_INTERVAL cutar;
      int r;
      SCIP_Bool integralslack = FALSE;

      r = rowinds[i]; /*lint !e613*/
      assert(0 <= r && r < SCIPgetNLPRows(scip));
      assert(slacksign[i] == -1 || slacksign[i] == +1); /*lint !e613*/
      assert(!SCIPisZero(scip, weights[i])); /*lint !e613*/

      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      if( slacksign[i] == 1 )
         SCIPintervalSetRoundingModeDownwards();
      else
         SCIPintervalSetRoundingModeUpwards();

      /* get the slack's coefficient a'_r = weights[i] * scale in the aggregated row */
      SCIPintervalSet(&ar, weights[i]);
      SCIPintervalMulScalar(SCIPinfinity(scip), &ar, ar, scale);
      SCIPintervalMulScalar(SCIPinfinity(scip), &ar, ar, (double) slacksign[i]);

      /* calculate slack variable's coefficient a^_r in the cut */
      if( row->integral &&
            ((slacksign[i] == +1 && SCIPrealIsExactlyIntegral(row->rhs) &&  SCIPrealIsExactlyIntegral(row->constant))
            || (slacksign[i] == -1 && SCIPrealIsExactlyIntegral(row->lhs) &&  SCIPrealIsExactlyIntegral(row->constant))) ) /*lint !e613*/
      {
         /* slack variable is always integral:
          *    a^_r = a~_r = down(a'_r)                      , if f_r <= f0
          *    a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
          */
         SCIP_Real downar;
         SCIP_INTERVAL fr;

         SCIPdebugMessage("resubstituting integer slack for row %s\n", row->name);
         downar = floor(ar.inf);
         SCIPintervalSubScalar(SCIPinfinity(scip), &fr, ar, downar);

         integralslack = TRUE;

         if( SCIPisLE(scip, fr.inf, f0.inf) )
         {
            SCIPintervalSet(&cutar, downar);
            splitcoef = downar;
            slackweight = weights[i];
            slackroundeddown = TRUE;
            SCIPintervalMul(SCIPinfinity(scip), &fr, fr, onedivoneminusf0);
            SCIPdebugMessage("fractionality %g, f0 %g -> round down to %g\n", fr.inf, f0.inf, splitcoef);
         }
         else
         {
            SCIPintervalSetBounds(&cutar, ar.inf, ar.sup);
            SCIPintervalSubScalar(SCIPinfinity(scip), &cutar, cutar, downar);
            SCIPintervalSub(SCIPinfinity(scip), &cutar, cutar, f0);
            SCIPintervalMul(SCIPinfinity(scip), &cutar, cutar, onedivoneminusf0);
            SCIPintervalAddScalar(SCIPinfinity(scip), &cutar, cutar, downar);
            splitcoef = downar + 1;
            slackweight = weights[i];
            slackroundeddown = FALSE;
            SCIPdebugMessage("fractionality %g, f0 %g -> round up! splitcoef %g sub-coefficient %g", fr.inf, f0.inf, splitcoef, cutar.inf);
         }
      }
      else
      {
         /* slack variable is continuous:
          *    a^_r = a~_r = 0                               , if a'_r >= 0
          *    a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
          */
         if( SCIPintervalGetInf(ar) >= 0.0 )
            continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
         else
         {
            SCIPintervalMul(SCIPinfinity(scip), &cutar, onedivoneminusf0, ar); /* cutaj = varsign[i] * aj * onedivoneminusf0; // a^_j */
            SCIPdebugMessage("resubstituting negative continuous slack for row %s with coef %g\n", row->name, cutar.inf);
         }
      }

      rowexact = SCIProwGetRowExact(row);
      assert(SCIProwExactHasFpRelax(rowexact));
      if( SCIProwExactGetRowRhs(rowexact) != NULL && slacksign[i] == 1.0 )
         userow = SCIProwExactGetRowRhs(rowexact);
      else
         userow = row;

      SCIPintervalMulScalar(SCIPinfinity(scip), &cutar, cutar, (double) -slacksign[i]);

      if( slacksign[i] == -1 )
         mult = cutar.inf;
      else
         mult = cutar.sup;

      if( SCIPisCertified(scip) && integralslack) /*lint --e{644}*/
      {
         assert(mirinfo != NULL);
         /* save the value for the split disjunction for the integer slack and the continous part (for rounded up we
          * subtract 1-f); multiply by -slacksign (same as above) since slack = side - row
          */
         mirinfo->slackrows[mirinfo->nslacks] = userow;
         SCIP_CALL( SCIPcaptureRow(scip, userow) );
         mirinfo->slackcoefficients[mirinfo->nslacks] = splitcoef * (-slacksign[i]);
         mirinfo->slacksign[mirinfo->nslacks] = slacksign[i];
         assert(SCIPrealIsExactlyIntegral(splitcoef));
         mirinfo->slackweight[mirinfo->nslacks] = slackweight;
         mirinfo->slackscale[mirinfo->nslacks] = scale;
         mirinfo->slackusedcoef[mirinfo->nslacks] = mult;

         /* save the value that goes into the certificate aggregation row (either downar or ar) */
         mirinfo->slackroundeddown[mirinfo->nslacks] = slackroundeddown;
         if( slackroundeddown )
            mirinfo->nrounddownslacks++;
         mirinfo->nslacks++;
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( EPSZ(SCIPintervalGetInf(cutar), QUAD_EPSILON) && (SCIPintervalGetInf(cutar) >= 0.0) )
         continue;

      /* depending on the slack's sign, we have
       * - sign = 1: s = rhs - a^Tx >= 0
       * - sign = -1: s = lhs - a^Tx <= 0
       */
      {
         SCIP_Bool success = TRUE;
         SCIP_Real sidevalchg;

         if( SCIPisCertified(scip) && !integralslack )
         {
            assert(aggrinfo != NULL);
            assert(aggrinfo->negslackweights[currentnegslackrow] == -weights[i]); /*lint !e777*/
            aggrinfo->substfactor[currentnegslackrow] = mult;
            currentnegslackrow++;
         }

         SCIP_CALL( varVecAddScaledRowCoefsSafely(scip, cutinds, cutcoefs, nnz, userow, mult, &sidevalchg, &success) );
         assert(success);

         /* move to rhs -> need to round up */
         SCIPintervalSetRoundingModeUpwards();
         *cutrhs += sidevalchg;
      }

      /* move slack's constant to the right hand side */
      if( slacksign[i] == +1 ) /*lint !e613*/
      {
         SCIP_INTERVAL rowrhs;

         SCIPintervalSetRoundingModeUpwards();
         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
         assert(!SCIPisInfinity(scip, userow->rhs));
         SCIPintervalSet(&rowrhs, userow->rhs);
         SCIPintervalSubScalar(SCIPinfinity(scip), &rowrhs, rowrhs, userow->constant);
#ifdef SCIP_DISABLED_CODE
         /* this is disabled because we can't certify it yet in exact solving mode; if enabled change also in addOneRowSafely() */
         if( row->integral )
         {
            /* the right hand side was implicitly rounded down in row aggregation */
            QUAD_ASSIGN(rowrhs, floor(QUAD_TO_DBL(rowrhs)));
         }
#endif
         SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, cutar, rowrhs);
         *cutrhs += SCIPintervalGetSup(tmpinterval);
      }
      else
      {
         SCIP_INTERVAL rowlhs;

         SCIPintervalSetRoundingModeUpwards();
         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
         assert(!SCIPisInfinity(scip, -userow->lhs));
         SCIPintervalSet(&rowlhs, userow->lhs);
         SCIPintervalSubScalar(SCIPinfinity(scip), &rowlhs, rowlhs, userow->constant);
#ifdef SCIP_DISABLED_CODE
         /* this is disabled because we can't certify it yet in exact solving mode; if enabled change also in addOneRowSafely() */
         if( row->integral )
         {
            /* the left hand side was implicitly rounded up in row aggregation */
            QUAD_ASSIGN(rowlhs, floor(QUAD_TO_DBL(rowlhs)));
         }
#endif
         SCIPintervalMul(SCIPinfinity(scip), &tmpinterval, cutar, rowlhs);
         *cutrhs += SCIPintervalGetSup(tmpinterval);
      }
   }

   /* relax rhs to zero, if it's very close to 0 */
   if( *cutrhs < 0.0 && *cutrhs >= SCIPepsilon(scip) )
      *cutrhs = 0.0;

   SCIPintervalSetRoundingMode(previousroundmode);

   return SCIP_OKAY;
}

#ifdef SCIP_DISABLED_CODE
/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row: \f$ a^\prime_r = scale * weight[r] * slacksign[r]. \f$
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 *  \f[
 *  \begin{array}{rll}
 *    integers : & \hat{a}_r = \tilde{a}_r = down(a^\prime_r),                      & \mbox{if}\qquad f_r <= f0 \\
 *               & \hat{a}_r = \tilde{a}_r = down(a^\prime_r) + (f_r - f0)/(1 - f0),& \mbox{if}\qquad f_r >  f0 \\
 *    continuous:& \hat{a}_r = \tilde{a}_r = 0,                                     & \mbox{if}\qquad a^\prime_r >= 0 \\
 *               & \hat{a}_r = \tilde{a}_r = a^\prime_r/(1 - f0),                   & \mbox{if}\qquad a^\prime_r <  0
 *  \end{array}
 *  \f]
 *
 *  Substitute \f$ \hat{a}_r \cdot s_r \f$ by adding \f$ \hat{a}_r \f$ times the slack's definition to the cut.
 *
 * @note this method is safe for usage in exact solving mode
 */
static
SCIP_RETCODE cutsSubstituteMIRRational(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   int*                  rowinds,            /**< sparsity pattern of used rows */
   int                   nrowinds,           /**< number of used rows */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   SCIP_RATIONAL*        f0                  /**< fractional value of rhs */
   )
{  /*lint --e{715}*/
   SCIP_ROW** rows;
   SCIP_ROW* userow;
   SCIP_ROWEXACT* rowexact;
   SCIP_RATIONAL* onedivoneminusf0;
   SCIP_RATIONAL* tmprational;
   SCIP_ROUNDMODE previousroundmode;
   SCIP_AGGREGATIONINFO* aggrinfo;
   SCIP_Real mult;
   int i;
   int currentnegslackrow;
   SCIP_RATIONAL* ar;
   SCIP_RATIONAL* cutar;
   int r;

   assert(scip != NULL);
   assert(weights != NULL || nrowinds == 0);
   assert(slacksign != NULL || nrowinds == 0);
   assert(rowinds != NULL || nrowinds == 0);
   assert(scale > 0.0);
   assert(cutcoefs != NULL);
   assert(QUAD_HI(cutrhs) != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(SCIPisExact(scip));

   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &tmprational) );
   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &onedivoneminusf0) );
   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &ar) );
   SCIP_CALL( SCIPrationalCreateBuffer(SCIPbuffer(scip), &cutar) );

   /* compute 1/(1-f0) in interval arithmetic */
   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPrationalMultReal(onedivoneminusf0, f0, -1);
   SCIPrationalAddReal(onedivoneminusf0, onedivoneminusf0, 1.0);
   SCIPrationalInvert(onedivoneminusf0, onedivoneminusf0);

   if( SCIPisCertified(scip)   )
   {
      aggrinfo = SCIPgetCertificate(scip)->aggrinfo[SCIPgetCertificate(scip)->naggrinfos -1];
   }

   rows = SCIPgetLPRows(scip);
   currentnegslackrow = 0;
   for( i = 0; i < nrowinds; i++ )
   {
      SCIP_ROW* row;

      r = rowinds[i]; /*lint !e613*/
      assert(0 <= r && r < SCIPgetNLPRows(scip));
      assert(slacksign[i] == -1 || slacksign[i] == +1); /*lint !e613*/
      assert(!SCIPisZero(scip, weights[i])); /*lint !e613*/

      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      if( slacksign[i] == 1 )
         SCIPintervalSetRoundingModeDownwards();
      else
         SCIPintervalSetRoundingModeUpwards();

      /* get the slack's coefficient a'_r = weights[i] * scale in the aggregated row */
      SCIPrationalSetReal(ar, weights[i]);
      SCIPrationalMultReal(ar, ar, scale);
      SCIPrationalMultReal(ar, ar, slacksign[i]);

      /* calculate slack variable's coefficient a^_r in the cut */
#ifdef SCIP_DISABLED_CODE
      if( row->integral && !SCIPisExact(scip)
         && ((slacksign[i] == +1 && SCIPisFeasIntegral(scip, row->rhs - row->constant))
            || (slacksign[i] == -1 && SCIPisFeasIntegral(scip, row->lhs - row->constant))) ) /*lint !e613*/
      {
         /* slack variable is always integral:
          *    a^_r = a~_r = down(a'_r)                      , if f_r <= f0
          *    a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
          */
         if( !SCIPisExact(scip) )
            downar = EPSFLOOR(ar, QUAD_EPSILON);
         else
            downar = floor(ar);

         SCIPquadprecSumDD(fr, ar, -downar);
         if( SCIPisLE(scip, QUAD_TO_DBL(fr), QUAD_TO_DBL(f0)) && (!SCIPisExact(scip) || QUAD_TO_DBL(fr) <= QUAD_TO_DBL(f0)) )
         {
            QUAD_ASSIGN(cutar, downar);
         }
         else
         {
            SCIPquadprecSumQQ(cutar, fr, -f0);
            SCIPquadprecProdQQ(cutar, cutar, onedivoneminusf0);
            SCIPquadprecSumQD(cutar, cutar, downar);
         }
      }
      else
#endif
      {
         /* slack variable is continuous:
          *    a^_r = a~_r = 0                               , if a'_r >= 0
          *    a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
          */
         if( !SCIPrationalIsNegative(ar) )
            continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
         else
         {
            SCIPrationalMult(cutar, ar, onedivoneminusf0); /* cutaj = varsign[i] * aj * onedivoneminusf0; // a^_j */
         }
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( SCIPrationalIsZero(cutar) )
         continue;

      /* depending on the slack's sign, we have
       *   sign = 1: s = rhs - a^Tx >= 0
           sign = -1: s = lhs - a^Tx <= 0
       */

      rowexact = SCIProwGetRowExact(row);
      assert(SCIProwExactHasFpRelax(rowexact));
      if( SCIProwExactGetRowRhs(rowexact) != NULL && slacksign[i] == 1.0 )
         userow = SCIProwExactGetRowRhs(rowexact);
      else
         userow = row;

      {
         SCIP_Bool success = TRUE;
         SCIP_Real sidevalchg;

         SCIPrationalMultReal(cutar, cutar, -slacksign[i]);
         if( slacksign[i] == -1 )
            mult = SCIPrationalRoundReal(cutar, SCIP_R_ROUND_DOWNWARDS);
         else
            mult = SCIPrationalRoundReal(cutar, SCIP_R_ROUND_UPWARDS);

         if( SCIPisCertified(scip) )
         {
            assert(aggrinfo->negslackweights[currentnegslackrow] == -weights[i]);
            aggrinfo->substfactor[currentnegslackrow] = mult;
            currentnegslackrow++;
         }

         SCIP_CALL( varVecAddScaledRowCoefsSafely(scip, cutinds, cutcoefs, nnz, userow, mult, &sidevalchg, &success) );
         assert(success);

         /* move to rhs -> need to round up */
         SCIPintervalSetRoundingModeUpwards();
         *cutrhs += sidevalchg;
      }

      /* move slack's constant to the right hand side */
      if( slacksign[i] == +1 ) /*lint !e613*/
      {
         SCIP_INTERVAL valinterval;
         SCIP_INTERVAL cutarinterval;

         SCIPintervalSetRoundingModeUpwards();
         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
         assert(!SCIPisInfinity(scip, userow->rhs));
#ifdef SCIP_DISABLED_CODE
         /* this is disabled because we can't certify it yet in exact solving mode; if enabled change also in addOneRowSafely() */
         if( row->integral )
         {
            /* the right hand side was implicitly rounded down in row aggregation */
            QUAD_ASSIGN(rowrhs, floor(QUAD_TO_DBL(rowrhs)));
         }
#endif
         SCIPintervalSet(&valinterval, userow->rhs);
         SCIPintervalSubScalar(SCIPinfinity(scip), &valinterval, valinterval, userow->constant);
         SCIPintervalSetRational(&cutarinterval, cutar);
         SCIPintervalMul(SCIPinfinity(scip), &valinterval, valinterval, cutarinterval);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, SCIPintervalGetSup(valinterval));
      }
      else
      {
         SCIP_INTERVAL valinterval;
         SCIP_INTERVAL cutarinterval;

         SCIPintervalSetRoundingModeUpwards();
         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
         assert(!SCIPisInfinity(scip, -userow->lhs));
#ifdef SCIP_DISABLED_CODE
          /* this is disabled because we can't certify it yet in exact solving mode; if enabled change also in addOneRowSafely() */
         if( row->integral )
         {
            /* the left hand side was implicitly rounded up in row aggregation */
            QUAD_ASSIGN(rowlhs, floor(QUAD_TO_DBL(rowlhs)));
         }
#endif
         SCIPintervalSet(&valinterval, userow->lhs);
         SCIPintervalSubScalar(SCIPinfinity(scip), &valinterval, valinterval, userow->constant);
         SCIPintervalSetRational(&cutarinterval, cutar);
         SCIPintervalMul(SCIPinfinity(scip), &valinterval, valinterval, cutarinterval);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, SCIPintervalGetSup(valinterval));
      }
   }

   /* relax rhs to zero, if it's very close to 0 */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   if( SCIPisExact(scip) )
      SCIPintervalSetRoundingMode(previousroundmode);

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &cutar);
   SCIPrationalFreeBuffer(SCIPbuffer(scip), &ar);

   SCIPrationalFreeBuffer(SCIPbuffer(scip), &onedivoneminusf0);
   SCIPrationalFreeBuffer(SCIPbuffer(scip), &tmprational);

   return SCIP_OKAY;
}
#endif

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row: \f$ a^\prime_r = scale * weight[r] * slacksign[r]. \f$
 *
 *  Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
 *  \f[
 *  \begin{array}{rll}
 *    integers : & \hat{a}_r = \tilde{a}_r = down(a^\prime_r),                      & \mbox{if}\qquad f_r <= f0 \\
 *               & \hat{a}_r = \tilde{a}_r = down(a^\prime_r) + (f_r - f0)/(1 - f0),& \mbox{if}\qquad f_r >  f0 \\
 *    continuous:& \hat{a}_r = \tilde{a}_r = 0,                                     & \mbox{if}\qquad a^\prime_r >= 0 \\
 *               & \hat{a}_r = \tilde{a}_r = a^\prime_r/(1 - f0),                   & \mbox{if}\qquad a^\prime_r <  0
 *  \end{array}
 *  \f]
 *
 *  Substitute \f$ \hat{a}_r \cdot s_r \f$ by adding \f$ \hat{a}_r \f$ times the slack's definition to the cut.
 */
static
SCIP_RETCODE cutsSubstituteMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   int*                  rowinds,            /**< sparsity pattern of used rows */
   int                   nrowinds,           /**< number of used rows */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   QUAD(SCIP_Real        f0)                 /**< fractional value of rhs */
   )
{  /*lint --e{715}*/
   SCIP_ROW** rows;
   SCIP_Real QUAD(onedivoneminusf0);
   int i;

   assert(scip != NULL);
   assert(weights != NULL || nrowinds == 0);
   assert(slacksign != NULL || nrowinds == 0);
   assert(rowinds != NULL || nrowinds == 0);
   assert(scale > 0.0);
   assert(cutcoefs != NULL);
   assert(QUAD_HI(cutrhs) != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(0.0 < QUAD_TO_DBL(f0) && QUAD_TO_DBL(f0) < 1.0);
   assert(!SCIPisExact(scip));

   SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0);
   SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0);

   rows = SCIPgetLPRows(scip);
   for( i = 0; i < nrowinds; i++ )
   {
      SCIP_ROW* row;
      SCIP_Real QUAD(ar);
      SCIP_Real QUAD(downar);
      SCIP_Real QUAD(cutar);
      SCIP_Real QUAD(fr);
      SCIP_Real QUAD(tmp);
      SCIP_Real QUAD(myprod);
      int r;

      r = rowinds[i]; /*lint !e613*/
      assert(0 <= r && r < SCIPgetNLPRows(scip));
      assert(slacksign[i] == -1 || slacksign[i] == +1); /*lint !e613*/
      assert(!SCIPisZero(scip, weights[i])); /*lint !e613*/

      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* get the slack's coefficient a'_r in the aggregated row */
      SCIPquadprecProdDD(ar, slacksign[i] * scale, weights[i]);

      /* calculate slack variable's coefficient a^_r in the cut */
      if( row->integral )
      {
         /* slack variable is always integral:
          *    a^_r = a~_r = down(a'_r)                      , if f_r <= f0
          *    a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
          */
         SCIPquadprecEpsFloorQ(downar, ar, SCIPepsilon(scip)); /*lint !e666*/
         SCIPquadprecSumQQ(fr, ar, -downar);
         assert(QUAD_TO_DBL(fr) >= -SCIPepsilon(scip) && QUAD_TO_DBL(fr) < 1.0);

         if( SCIPisLE(scip, QUAD_TO_DBL(fr), QUAD_TO_DBL(f0)) )
            QUAD_ASSIGN_Q(cutar, downar); /* a^_r */
         else
         {
            SCIPquadprecSumQQ(cutar, fr, -f0);
            SCIPquadprecProdQQ(cutar, cutar, onedivoneminusf0);
            SCIPquadprecSumQQ(cutar, cutar, downar);
         }
      }
      else
      {
         /* slack variable is continuous:
          *    a^_r = a~_r = 0                               , if a'_r >= 0
          *    a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
          */
         if( QUAD_TO_DBL(ar) >= 0.0 )
            continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
         else
            SCIPquadprecProdQQ(cutar, onedivoneminusf0, ar);
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( EPSZ(QUAD_TO_DBL(cutar), QUAD_EPSILON) )
         continue;

      /* depending on the slack's sign, we have
       *   a*x + c + s == rhs  =>  s == - a*x - c + rhs,  or  a*x + c - s == lhs  =>  s == a*x + c - lhs
       * substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
       */
      SCIPquadprecProdQD(myprod, cutar, -slacksign[i]);

      /* add the slack's definition multiplied with a^_j to the cut */
      SCIP_CALL( varVecAddScaledRowCoefsQuadScale(cutinds, cutcoefs, nnz, row, QUAD(myprod)) );

      /* move slack's constant to the right hand side */
      if( slacksign[i] == +1 ) /*lint !e613*/
      {
         SCIP_Real QUAD(rowrhs);

         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
         assert(!SCIPisInfinity(scip, row->rhs));
         QUAD_ASSIGN(rowrhs, row->rhs - row->constant);
         if( row->integral )
         {
            /* the right hand side was implicitly rounded down in row aggregation */
            SCIPquadprecEpsFloorQ(rowrhs, rowrhs, SCIPepsilon(scip)); /*lint !e666*/
         }
         SCIPquadprecProdQQ(tmp, myprod, rowrhs);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
      }
      else
      {
         SCIP_Real QUAD(rowlhs);

         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
         assert(!SCIPisInfinity(scip, -row->lhs));
         QUAD_ASSIGN(rowlhs, row->lhs - row->constant);
         if( row->integral )
         {
            /* the left hand side was implicitly rounded up in row aggregation */
            SCIPquadprecEpsCeilQ(rowlhs, rowlhs, SCIPepsilon(scip)); /*lint !e666*/
         }
         SCIPquadprecProdQQ(tmp, myprod, rowlhs);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, tmp);
      }
   }

   /* relax rhs to zero, if it's very close to 0 */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   return SCIP_OKAY;
}

/** calculates an MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because
 *  these rows cannot participate in an MIR cut.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note this method is safe for usage in exact solving mode
 *
 *  @todo make behavior identical to the unsafe MIR cut computation
 */
static
SCIP_RETCODE calcMIRSafely(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to the aggrrow; must be positive */
   SCIP_AGGRROW*         aggrrow,            /**< aggrrow to compute MIR cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut if its efficacy improves cutefficacy */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut if its efficacy improves cutefficacy */
   int*                  cutinds,            /**< array to store the indices of non-zero coefficients in the cut if its efficacy improves cutefficacy */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut if its efficacy improves cutefficacy */
   SCIP_Real*            cutefficacy,        /**< pointer to store efficacy of cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut or NULL if it improves cutefficacy */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally if it improves cutefficacy */
   SCIP_Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut and it improves cutefficacy */
   )
{
   int i;
   int nvars;
   int tmpnnz;
   int* varsign;
   int* boundtype;
   int* tmpinds;
   SCIP_Real* tmpcoefs;

   SCIP_Real rhs;
   SCIP_Real downrhs;
   SCIP_Bool freevariable;
   SCIP_Bool localbdsused;
   SCIP_Bool tmpislocal;

   SCIP_ROUNDMODE previousroundmode;
   SCIP_INTERVAL f0interval;

   assert(aggrrow != NULL);
   assert(SCIPisPositive(scip, scale));
   assert(SCIPisExact(scip));
   assert(success != NULL);

   SCIPdebugMsg(scip, "calculating MIR cut (scale: %g)\n", scale);

   *success = FALSE;
   *cutislocal = FALSE;

   /* allocate temporary memory */
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtype, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpinds, nvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &tmpcoefs, nvars) );

   /* initialize cut with aggregation */
   tmpnnz = aggrrow->nnz;
   tmpislocal = aggrrow->local;

   previousroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeUpwards();

   if( SCIPisCertified(scip) )
   {
      SCIP_CALL( SCIPaddCertificateMirInfo(scip) );
   }

   rhs = QUAD_TO_DBL(aggrrow->rhs) * scale;

   if( tmpnnz > 0 )
   {
      BMScopyMemoryArray(tmpinds, aggrrow->inds, tmpnnz);

      for( i = 0; i < tmpnnz; ++i )
      {
         SCIP_Real coef;
         int k = aggrrow->inds[i];

         coef = aggrrow->vals[k];
         coef *= scale;
         tmpcoefs[k] = coef;

         assert(coef != 0.0);
      }

      SCIPdebugMsg(scip, "Initial row:\n");
      SCIPdebug(printCut(scip, sol, tmpcoefs, rhs, tmpinds, tmpnnz, FALSE, FALSE));

      /* Transform equation  a*x == b, lb <= x <= ub  into standard form
       *   a'*x' == b, 0 <= x' <= ub'.
       *
       * Transform variables (lb or ub):
       *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
       *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
       * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
       *
       * Transform variables (vlb or vub):
       *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
       *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
       * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
       *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
       *   a_{zu_j} := a_{zu_j} + a_j * bu_j
       */
      SCIP_CALL( cutsTransformMIRSafely(scip, sol, boundswitch, usevbds, allowlocal, fixintegralrhs, FALSE,
            boundsfortrans, boundtypesfortrans, tmpcoefs, &rhs, tmpinds, &tmpnnz, varsign, boundtype, &freevariable, &localbdsused) );
      assert(allowlocal || !localbdsused);
      tmpislocal = tmpislocal || localbdsused;

      if( freevariable )
         goto TERMINATE;

      SCIPdebugMsg(scip, "Aggregated and transformed:\n");
      SCIPdebug(printCut(scip, sol, tmpcoefs, rhs, tmpinds, tmpnnz, FALSE, FALSE));
   }

   /* Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers :  a~_j = down(a'_j)                      , if f_j <= f_0
    *             a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
    * continuous: a~_j = 0                               , if a'_j >= 0
    *             a~_j = a'_j/(1 - f0)                   , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */

   downrhs = floor(rhs);

   if( SCIPisCertified(scip)   )
   {
      SCIP_MIRINFO* mirinfo = SCIPgetCertificate(scip)->mirinfo[SCIPgetCertificate(scip)->nmirinfos - 1];
      SCIPrationalSetReal(mirinfo->rhs, downrhs);
      SCIPrationalSetReal(mirinfo->frac, rhs);
      SCIPrationalDiffReal(mirinfo->frac, mirinfo->frac, downrhs);
   }

   SCIPintervalSet(&f0interval, rhs);
   SCIPintervalSubScalar(SCIPinfinity(scip), &f0interval, f0interval, downrhs);

   if( f0interval.inf < minfrac || f0interval.sup > maxfrac )
      goto TERMINATE;

   /* We multiply the coefficients of the base inequality roughly by scale/(1-f0).
    * If this gives a scalar that is very big, we better do not generate this cut.
    */
   if( REALABS(scale)/(1.0 - f0interval.inf) > MAXCMIRSCALE )
      goto TERMINATE;

   /* renormalize f0 value */
   rhs = downrhs;

   if( tmpnnz > 0 )
   {
      SCIP_CALL( cutsRoundMIRSafely(scip, tmpcoefs, &rhs, tmpinds, &tmpnnz, varsign, boundtype, f0interval) ); /*lint !e644*/

      SCIPdebugMsg(scip, "After MIR rounding:\n");
      SCIPdebug(printCut(scip, sol, tmpcoefs, rhs, tmpinds, tmpnnz, FALSE, FALSE));
   }

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a^_r = a~_r = down(a'_r)                      , if f_r <= f0
    *               a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
    *   continuous: a^_r = a~_r = 0                               , if a'_r >= 0
    *               a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
    *
    * Substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
    */

   SCIP_CALL( cutsSubstituteMIRSafely(scip, aggrrow->rowweights, aggrrow->slacksign, aggrrow->rowsinds,
      aggrrow->nrows, scale, tmpcoefs, &rhs, tmpinds, &tmpnnz, f0interval) );

   SCIPdebugMsg(scip, "After slack substitution:\n");
   SCIPdebug( printCut(scip, sol, tmpcoefs, rhs, tmpinds, tmpnnz, FALSE, FALSE) );

   /* we work on rhs -> round up */
   SCIPintervalSetRoundingModeUpwards();

   if( postprocess )
   {
      SCIP_CALL( postprocessCutSafely(scip, tmpislocal, tmpinds, tmpcoefs, &tmpnnz, &rhs, success) );
   }
   else
   {
      *success = !removeZerosSafely(scip, SCIPsumepsilon(scip), tmpcoefs, &rhs, tmpinds, &tmpnnz);
   }

   SCIPdebugMsg(scip, "After post processing:\n");
   SCIPdebug( printCut(scip, sol, tmpcoefs, rhs, tmpinds, tmpnnz, FALSE, FALSE) );

   if( *success )
   {
      SCIP_Real mirefficacy = calcEfficacyDenseStorage(scip, sol, tmpcoefs, rhs, tmpinds, tmpnnz);

      if( SCIPisEfficacious(scip, mirefficacy) && (cutefficacy == NULL || mirefficacy > *cutefficacy) )
      {
         BMScopyMemoryArray(cutinds, tmpinds, tmpnnz);
         *cutnnz = tmpnnz;
         *cutrhs = rhs;
         *cutislocal = tmpislocal;

         /* clean tmpcoefs and go back to double precision */
         for( i = 0; i < *cutnnz; ++i )
         {
            int j = cutinds[i];

            cutcoefs[i] = tmpcoefs[j];
            tmpcoefs[j] = 0.0;
         }

         if( cutefficacy != NULL )
            *cutefficacy = mirefficacy;

         if( cutrank != NULL )
            *cutrank = aggrrow->rank + 1;
      }
      else
      {
         *success = FALSE;
      }
   }

  TERMINATE:

   /* reset the rounding mode in exact mode */
   SCIPintervalSetRoundingMode(previousroundmode); /*lint !e644*/

   if( !(*success) )
   {
      for( i = 0; i < tmpnnz; ++i )
      {
         tmpcoefs[tmpinds[i]] = 0.0;
      }
   }

   /* free temporary memory */
   SCIPfreeCleanBufferArray(scip, &tmpcoefs);
   SCIPfreeBufferArray(scip, &tmpinds);
   SCIPfreeBufferArray(scip, &boundtype);
   SCIPfreeBufferArray(scip, &varsign);

   return SCIP_OKAY;
}


/** calculates an MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because
 *  these rows cannot participate in an MIR cut.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 *
 *  @note this method is safe for usage in exact solving mode
 */
SCIP_RETCODE SCIPcalcMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   int                   vartypeusevbds,     /**< for all variable types with index smaller than this number, variable
                                              *   type substitution is allowed. The indices are: 0: continuous,
                                              *   1: continuous implint., 2: integer implint, 3: binary implint,
                                              *   4: integer, 5: binary */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to the aggrrow; must be positive */
   SCIP_AGGRROW*         aggrrow,            /**< aggrrow to compute MIR cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut if its efficacy improves cutefficacy */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut if its efficacy improves cutefficacy */
   int*                  cutinds,            /**< array to store the indices of non-zero coefficients in the cut if its efficacy improves cutefficacy */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut if its efficacy improves cutefficacy */
   SCIP_Real*            cutefficacy,        /**< pointer to store efficacy of cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut or NULL if it improves cutefficacy */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally if it improves cutefficacy */
   SCIP_Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut and it improves cutefficacy */
   )
{
   MIR_DATA* data;
   int* varsign;
   int* boundtype;
   SCIP_Real QUAD(downrhs);
   SCIP_Real QUAD(f0);
   SCIP_Bool freevariable;
   SCIP_Bool localbdsused;
   SCIP_Bool tmpislocal;

   assert(aggrrow != NULL);
   assert(SCIPisPositive(scip, scale));
   assert(success != NULL);

   if( SCIPisExact(scip) )
   {
      /* TODO: update exactSCIP cuts to behave identically with respect to implied integrality */
      return calcMIRSafely(scip, sol, postprocess, boundswitch, vartypeusevbds > 0 ? TRUE : FALSE, allowlocal, fixintegralrhs,
         boundsfortrans, boundtypesfortrans, minfrac, maxfrac, scale, aggrrow, cutcoefs, cutrhs,
         cutinds, cutnnz, cutefficacy, cutrank, cutislocal, success);
   }
   SCIPdebugMsg(scip, "calculating MIR cut (scale: %g)\n", scale);

   *success = FALSE;

   /* Setup data to track cut and initialize the cut with aggregation */
   int l;
   int nnz;

   assert(vartypeusevbds >= 0 && vartypeusevbds < NSECTIONS);

   SCIP_CALL(SCIPallocBuffer(scip, &data));

   nnz = aggrrow->nnz;
   data->totalnnz = nnz;

   /* initialize sections */
   for( l = 0; l < NSECTIONS; ++l )
   {
      SCIP_CALL(SCIPallocBufferArray(scip, &data->secindices[l], nnz));
      data->secnnz[l] = 0;
      /* Cont. | cont impl. | int impl. | bin impl. | int | bin */
      assert(NSECTIONS == 6); /* If the section definition is changed, the below lines should also be adjusted to match */
      data->isenfint[l] = l >= 2 ? TRUE : FALSE;
      data->isimplint[l] = l >= 1 && l <= 3 ? TRUE : FALSE;
      /* Use variable bounds for the sections specified by the user */
      data->usevbds[l] = l < vartypeusevbds ? 2 : 0;
   }

   /* Problem data needs to be initialized before cut data as it is used to partition the variables into the sections */
   data->vars = SCIPgetVars(scip);
   data->nvars = SCIPgetNVars(scip);
   data->nbinvars = SCIPgetNBinVars(scip);
   data->nintvars = SCIPgetNIntVars(scip);
   data->nbinimplvars = SCIPgetNBinImplVars(scip);
   data->nintimplvars = SCIPgetNIntImplVars(scip);
   data->ncontimplvars = SCIPgetNContImplVars(scip);
   data->ncontvars = SCIPgetNContVars(scip);

   SCIP_CALL(SCIPallocCleanBufferArray(scip, &( data->cutcoefs ), QUAD_ARRAY_SIZE(data->nvars)));
   SCIP_CALL(SCIPallocBufferArray(scip, &data->cutinds, data->nvars));

   SCIPquadprecProdQD(data->cutrhs, aggrrow->rhs, scale);

   if( nnz > 0 )
   {
      /* Initalize cut with the aggregation */
      BMScopyMemoryArray(data->cutinds, aggrrow->inds, nnz);

      for( l = 0; l < nnz; ++l )
      {
         SCIP_Real QUAD(coef);
         int m = aggrrow->inds[l];

         QUAD_ARRAY_LOAD(coef, aggrrow->vals, m);

         SCIPquadprecProdQD(coef, coef, scale);

         QUAD_ARRAY_STORE(data->cutcoefs, m, coef);

         assert(QUAD_HI(coef) != 0.0);
      }

      /* Sort the array by problem index and add the variables to their sections */
      SCIPsortDownInt(data->cutinds, nnz);
      for( l = 0; l < nnz; ++l )
      {
         int section = varSection(data, data->cutinds[l]);
         data->secindices[section][data->secnnz[section]] = data->cutinds[l];
         ++data->secnnz[section];
      }
   }

   SCIPdebugMsg(scip, "Initial row:\n");
   SCIPdebug( printCutQuad(scip, sol, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, nnz, FALSE, FALSE) );

   data->ncutinds = 0;
   tmpislocal = aggrrow->local;

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, data->nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtype, data->nvars) );

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a'*x' == b, 0 <= x' <= ub'.
    *
    * Transform variables (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
    * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
    *
    * Transform variables (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
    * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
    *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
    *   a_{zu_j} := a_{zu_j} + a_j * bu_j
    */
   if( data->totalnnz > 0 )
   {
      SCIP_CALL( cutsTransformMIR(scip, data, sol, boundswitch, allowlocal, fixintegralrhs, FALSE,
            boundsfortrans, boundtypesfortrans, minfrac, maxfrac,
            varsign, boundtype, &freevariable, &localbdsused) );
      assert(allowlocal || !localbdsused);
      tmpislocal = tmpislocal || localbdsused;

      if( freevariable )
         goto TERMINATE;

      SCIPdebugMsg(scip, "Aggregated and transformed:\n");
      SCIPdebug(printCutQuad(scip, sol, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, data->ncutinds, FALSE, FALSE));
   }

   /* Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers :  a~_j = down(a'_j)                      , if f_j <= f_0
    *             a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
    * continuous: a~_j = 0                               , if a'_j >= 0
    *             a~_j = a'_j/(1 - f0)                   , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */

   SCIPquadprecEpsFloorQ(downrhs, data->cutrhs, SCIPepsilon(scip)); /*lint !e666*/
   SCIPquadprecSumQQ(f0, data->cutrhs, -downrhs);
   assert(QUAD_TO_DBL(f0) >= -SCIPepsilon(scip) && QUAD_TO_DBL(f0) < 1.0);

   if( QUAD_TO_DBL(f0) < minfrac || QUAD_TO_DBL(f0) > maxfrac )
      goto TERMINATE;

   /* We multiply the coefficients of the base inequality roughly by scale/(1-f0).
    * If this gives a scalar that is very big, we better do not generate this cut.
    */
   if( REALABS(scale)/(1.0 - QUAD_TO_DBL(f0)) > MAXCMIRSCALE )
      goto TERMINATE;

   /* renormalize f0 value */
   SCIPquadprecSumDD(f0, QUAD_HI(f0), QUAD_LO(f0));

   QUAD_ASSIGN_Q(data->cutrhs, downrhs);

   if( data->totalnnz > 0 )
   {
      SCIP_CALL( cutsRoundMIR(scip, data, varsign, boundtype, QUAD(f0)) );

      SCIPdebugMsg(scip, "After MIR rounding:\n");
      SCIPdebug(printCutQuad(scip, sol, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, data->ncutinds, FALSE, FALSE));
   }

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a^_r = a~_r = down(a'_r)                      , if f_r <= f0
    *               a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
    *   continuous: a^_r = a~_r = 0                               , if a'_r >= 0
    *               a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
    *
    * Substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
    */
   SCIP_CALL( cutsSubstituteMIR(scip, aggrrow->rowweights, aggrrow->slacksign, aggrrow->rowsinds,
         aggrrow->nrows, scale,
         data->cutcoefs, QUAD(&data->cutrhs), data->cutinds, &data->ncutinds, QUAD(f0)) );

   SCIPdebugMsg(scip, "After slack substitution:\n");
   SCIPdebug(printCutQuad(scip, sol, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, data->ncutinds, FALSE, FALSE));

   if( postprocess )
   {
      /* remove all nearly-zero coefficients from MIR row and relax the right hand side correspondingly in order to
       * prevent numerical rounding errors
       */
      SCIP_CALL( postprocessCutQuad(scip, tmpislocal, data->cutinds, data->cutcoefs, &data->ncutinds, QUAD(&data->cutrhs), success) );
   }
   else
   {
      *success = ! removeZerosQuad(scip, SCIPsumepsilon(scip), tmpislocal, data->cutcoefs, QUAD(&data->cutrhs), data->cutinds, &data->ncutinds);
   }

   SCIPdebugMsg(scip, "After post processing:\n");
   SCIPdebug( printCutQuad(scip, sol, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, data->ncutinds, FALSE, FALSE) );

   if( *success )
   {
      SCIP_Real mirefficacy = calcEfficacyDenseStorageQuad(scip, sol, data->cutcoefs, QUAD_TO_DBL(data->cutrhs), data->cutinds, data->ncutinds);

      if( SCIPisEfficacious(scip, mirefficacy) && (cutefficacy == NULL || mirefficacy > *cutefficacy) )
      {
         BMScopyMemoryArray(cutinds, data->cutinds, data->ncutinds);
         *cutnnz = data->ncutinds;
         *cutrhs = QUAD_TO_DBL(data->cutrhs);
         *cutislocal = tmpislocal;

         /* clean tmpcoefs and go back to double precision */
         for(int i = 0; i < *cutnnz; ++i )
         {
            SCIP_Real QUAD(coef);
            int j = cutinds[i];

            QUAD_ARRAY_LOAD(coef, data->cutcoefs, j);

            cutcoefs[i] = QUAD_TO_DBL(coef);
            QUAD_ASSIGN(coef, 0.0);
            QUAD_ARRAY_STORE(data->cutcoefs, j, coef);
         }

         if( cutefficacy != NULL )
            *cutefficacy = mirefficacy;

         if( cutrank != NULL )
            *cutrank = aggrrow->rank + 1;
      }
      else
      {
         *success = FALSE;
      }
   }

  TERMINATE:
   if( !(*success) )
   {
      SCIP_Real QUAD(tmp);

      QUAD_ASSIGN(tmp, 0.0);
      for(int i = 0; i < data->ncutinds; ++i )
      {
         QUAD_ARRAY_STORE(data->cutcoefs, data->cutinds[i], tmp);
      }
   }

#ifndef NDEBUG
   for( int i = 0; i < QUAD_ARRAY_SIZE(data->nvars); ++i )
   {
      if(data->cutcoefs[i] != 0.0)
      {
         SCIPdebugMsg(scip, "coefs have not been reset\n");
         SCIPABORT();
      }
   }
#endif

   SCIPfreeBufferArray(scip, &boundtype);
   SCIPfreeBufferArray(scip, &varsign);

   if( data->cutinds != NULL )
      SCIPfreeBufferArray(scip, &data->cutinds);

   if( data->cutcoefs != NULL )
      SCIPfreeCleanBufferArray(scip, &data->cutcoefs);

   for( int s = NSECTIONS - 1; s >= 0; --s )
   {
      SCIPfreeBufferArray(scip, &data->secindices[s]);
   }

   SCIPfreeBuffer(scip, &data);

   return SCIP_OKAY;
}

/** compute the efficacy of the MIR cut for the given values without computing the cut.
 *  This is used for the CMIR cut generation heuristic.
 */
static
SCIP_Real computeMIREfficacy(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*RESTRICT    coefs,              /**< array with coefficients in row */
   SCIP_Real*RESTRICT    solvals,            /**< solution values of variables in the row */
   SCIP_Real             rhs,                /**< right hand side of MIR cut */
   SCIP_Real             contactivity,       /**< aggregated activity of continuous variables in the row */
   SCIP_Real             contsqrnorm,        /**< squared norm of continuous variables */
   SCIP_Real             delta,              /**< delta value to compute the violation for */
   int                   nvars,              /**< number of variables in the row, i.e. the size of coefs and solvals arrays */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac             /**< maximal fractionality of rhs to produce MIR cut for */
   )
{
   int i;
   SCIP_Real f0pluseps;
   SCIP_Real f0;
   SCIP_Real onedivoneminusf0;
   SCIP_Real scale;
   SCIP_Real downrhs;
   SCIP_Real norm;
   SCIP_Real contscale;

   scale = 1.0 / delta;
   rhs *= scale;
   downrhs = SCIPfloor(scip, rhs);
   f0 = rhs - downrhs;

   if( f0 < minfrac || f0 > maxfrac )
      return 0.0;

   onedivoneminusf0 = 1.0 / (1.0 - f0);

   contscale = scale * onedivoneminusf0;

   /* We multiply the coefficients of the base inequality roughly by scale/(1-f0).
    * If this gives a scalar that is very big, we better do not generate this cut.
    */
   if( contscale > MAXCMIRSCALE )
      return 0.0;

   rhs = downrhs;
   rhs -= contscale * contactivity;
   norm = SQR(contscale) * contsqrnorm;

   assert(!SCIPisFeasZero(scip, f0));
   assert(!SCIPisFeasZero(scip, 1.0 - f0));

   f0pluseps = f0 + SCIPepsilon(scip);

   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real floorai = SCIPfloor(scip, scale * coefs[i]);
      SCIP_Real fi = (scale * coefs[i]) - floorai;

      if( fi > f0pluseps )
         floorai += (fi - f0) * onedivoneminusf0;

      rhs -= solvals[i] * floorai;
      norm += SQR(floorai);
   }

   norm = sqrt(norm);

   return - rhs / MAX(norm, 1e-6);
}

/** calculates an MIR cut out of an aggregation of LP rows
 *
 *  Given the aggregation, it is transformed to a mixed knapsack set via complementation (using bounds or variable bounds)
 *  Then, different scalings of the mkset are used to generate a MIR and the best is chosen.
 *  One of the steps of the MIR is to round the coefficients of the integer variables down,
 *  so one would prefer to have integer coefficients for integer variables which are far away from their bounds in the
 *  mkset.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcutGenerationHeuristicCMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   int                   vartypeusevbds,     /**< for all variable types with index smaller than this number, variable
                                              *   type substitution is allowed. The indices are: 0: continuous,
                                              *   1: continuous implint., 2: integer implint, 3: binary implint,
                                              *   4: integer, 5: binary */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int                   maxtestdelta,       /**< maximum number of deltas to test */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_AGGRROW*         aggrrow,            /**< aggrrow to compute MIR cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store efficacy of best cut; only cuts that are strictly better than the value of
                                              *   this efficacy on input to this function are returned */
   int*                  cutrank,            /**< pointer to return rank of generated cut (or NULL) */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid and efficacious cut was returned */
   )
{
   int i;
   int firstcontvar;
   int nvars;
   int intstart;
   int ntmpcoefs;
   int* varsign;
   int* boundtype;
   int* mksetinds;
   SCIP_Real* mksetcoefs;
   SCIP_Real QUAD(mksetrhs);
   int mksetnnz;
   SCIP_Real* bounddist;
   int* bounddistpos;
   int nbounddist;
   SCIP_Real* tmpcoefs;
   SCIP_Real* tmpvalues;
   SCIP_Real* deltacands;
   int ndeltacands;
   SCIP_Real bestdelta;
   SCIP_Real bestefficacy;
   SCIP_Real maxabsmksetcoef;
   SCIP_VAR** vars;
   SCIP_Bool freevariable;
   SCIP_Bool localbdsused;
   SCIP_Real contactivity;
   SCIP_Real contsqrnorm;
   MIR_DATA* data;

   assert(aggrrow != NULL);
   assert(aggrrow->nrows + aggrrow->nnz >= 1);
   assert(success != NULL);

   *success = FALSE;
   nvars = SCIPgetNVars(scip);
   firstcontvar = nvars - SCIPgetNContVars(scip) - SCIPgetNContImplVars(scip);
   vars = SCIPgetVars(scip);

   /* allocate temporary memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtype, nvars) );

   SCIP_CALL( SCIPallocBufferArray(scip, &tmpcoefs, nvars + aggrrow->nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvalues, nvars + aggrrow->nrows) );
   /* The +4 comes from a few rules that create extra delta candidates, see usages of ndeltacands. */
   SCIP_CALL( SCIPallocBufferArray(scip, &deltacands, nvars + 4) );

   /* we only compute bound distance for integer variables; by variable bound substitution, the number of integer variables
    * can grow significantly. Hence, these allocations are length nvars */
   SCIP_CALL( SCIPallocBufferArray(scip, &bounddist, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounddistpos, nvars) );

   /* initialize mkset with the unscaled aggregation */
   {
      int l;
      int nnz;

      assert(vartypeusevbds >= 0 && vartypeusevbds < NSECTIONS);

      SCIP_CALL( SCIPallocBuffer(scip, &data) );

      nnz = aggrrow->nnz;
      data->totalnnz = nnz;

      /* initialize sections */
      for( l = 0; l < NSECTIONS; ++l )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &data->secindices[l], nnz) );
         data->secnnz[l] = 0;
         /* Cont. | cont impl. | int impl. | bin impl. | int | bin */
         assert(NSECTIONS == 6); /* If the section definition is changed, the below lines should also be adjusted to match */
         data->isenfint[l] = l >= 2 ? TRUE : FALSE;
         data->isimplint[l] = l >= 1 && l <= 3 ? TRUE : FALSE;
         /* Use variable bounds for the sections specified by the user */
         data->usevbds[l] = l < vartypeusevbds ? 2 : 0;
      }

      /* Problem data needs to be initialized before cut data as it is used to partition the variables into the sections */
      data->vars = SCIPgetVars(scip);
      data->nvars = SCIPgetNVars(scip);
      data->nbinvars = SCIPgetNBinVars(scip);
      data->nintvars = SCIPgetNIntVars(scip);
      data->nbinimplvars = SCIPgetNBinImplVars(scip);
      data->nintimplvars = SCIPgetNIntImplVars(scip);
      data->ncontimplvars = SCIPgetNContImplVars(scip);
      data->ncontvars = SCIPgetNContVars(scip);

      SCIP_CALL( SCIPallocCleanBufferArray(scip, &(data->cutcoefs), QUAD_ARRAY_SIZE(data->nvars)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &data->cutinds, data->nvars) );

      SCIPquadprecProdQD(data->cutrhs, aggrrow->rhs, 1.0);

      if( nnz > 0 )
      {
         /* Initalize cut with the aggregation */
         BMScopyMemoryArray(data->cutinds, aggrrow->inds, nnz);

         for( l = 0; l < nnz; ++l )
         {
            SCIP_Real QUAD(coef);
            int m = aggrrow->inds[l];

            QUAD_ARRAY_LOAD(coef, aggrrow->vals, m);

            SCIPquadprecProdQD(coef, coef, 1.0);

            QUAD_ARRAY_STORE(data->cutcoefs, m, coef);

            assert(QUAD_HI(coef) != 0.0);
         }

         /* Sort the array by problem index and add the variables to their sections */
         SCIPsortDownInt(data->cutinds, nnz);
         for( l = 0; l < nnz; ++l )
         {
            int section = varSection(data, data->cutinds[l]);
            data->secindices[section][data->secnnz[section]] = data->cutinds[l];
            ++data->secnnz[section];
         }
      }

      data->ncutinds = 0;
   }
   *cutislocal = aggrrow->local;

   /* Transform equation  a*x == b, lb <= x <= ub  into standard form
    *   a'*x' == b, 0 <= x' <= ub'.
    *
    * Transform variables (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
    * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
    *
    * Transform variables (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
    * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
    *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
    *   a_{zu_j} := a_{zu_j} + a_j * bu_j
    */
   SCIP_CALL( cutsTransformMIR(scip, data, sol, boundswitch, allowlocal, FALSE, FALSE,
         boundsfortrans, boundtypesfortrans, minfrac, maxfrac, varsign, boundtype, &freevariable, &localbdsused) );
   assert(allowlocal || !localbdsused);

   /* Use aliases to stay more consistent with the old code. mksetrhs needs to synchronize its values data->cutrhs
    * again before calling SCIProundMIR()! */
   mksetinds = data->cutinds;
   mksetcoefs = data->cutcoefs;
   mksetnnz = data->ncutinds;
   QUAD_ASSIGN_Q(mksetrhs, data->cutrhs);

   if( freevariable )
      goto TERMINATE;

   SCIPdebugMsg(scip, "transformed aggrrow row:\n");
   SCIPdebug( printCutQuad(scip, sol, mksetcoefs, QUAD(mksetrhs), mksetinds, mksetnnz, FALSE, FALSE) );

   /* found positions of integral variables that are strictly between their bounds */
   maxabsmksetcoef = -1.0;
   nbounddist = 0;

   assert(mksetnnz <= nvars);
   for( i = mksetnnz - 1; i >= 0 && mksetinds[i] < firstcontvar; --i )
   {
      SCIP_VAR* var = vars[mksetinds[i]];
      SCIP_Real primsol = SCIPgetSolVal(scip, sol, var);
      SCIP_Real lb = SCIPvarGetLbLocal(var);
      SCIP_Real ub = SCIPvarGetUbLocal(var);
      SCIP_Real QUAD(coef);

      QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[i]);

      if( SCIPisEQ(scip, primsol, lb) || SCIPisEQ(scip, primsol, ub) )
         continue;

      bounddist[nbounddist] = MIN(ub - primsol, primsol - lb);
      bounddistpos[nbounddist] = i;
      deltacands[nbounddist] = QUAD_TO_DBL(coef);
      ++nbounddist;
   }
   assert(nbounddist <= nvars);

   /* no fractional variable; so abort here */
   if( nbounddist == 0 )
      goto TERMINATE;

   intstart = i + 1;

   /* Check that the continuous and implied integer variables and integer variables are partitioned */
#ifndef NDEBUG
   for( int j = 0; j < intstart; ++j )
   {
      assert(SCIPvarGetType(vars[mksetinds[j]]) == SCIP_VARTYPE_CONTINUOUS);
   }
   for( int j = intstart; j < data->ncutinds; ++j )
   {
      assert(SCIPvarIsIntegral(vars[mksetinds[j]]));
   }
#endif

   ndeltacands = nbounddist;
   assert(ndeltacands <= nvars);
   SCIPsortDownRealRealInt(bounddist, deltacands, bounddistpos, nbounddist);

   {
      SCIP_Real intscale;
      SCIP_Bool intscalesuccess;

      SCIP_CALL( SCIPcalcIntegralScalar(deltacands, nbounddist, -SCIPepsilon(scip), SCIPsumepsilon(scip), (SCIP_Longint)10000, 10000.0, &intscale, &intscalesuccess) );

      if( intscalesuccess )
      {
         SCIP_Real intf0;
         SCIP_Real intscalerhs;
         SCIP_Real delta;

         intscalerhs = QUAD_TO_DBL(data->cutrhs) * intscale;
         delta = 1.0 / intscale;
         intf0 = intscalerhs - floor(intscalerhs);

         if( ! SCIPisFeasIntegral(scip, intf0) )
         {
            if( intf0 < minfrac || intf0 > maxfrac )
            {
               intscale *= SCIPceil(scip, MAX(minfrac, (1.0 - maxfrac)) / MIN(intf0, (1.0 - intf0)));
               intscalerhs = QUAD_TO_DBL(data->cutrhs) * intscale;
               delta = 1.0 / intscale;
               intf0 = intscalerhs - floor(intscalerhs);
            }

            if( intf0 >= minfrac && intf0 <= maxfrac )
            {
               if( ! SCIPisEQ(scip, delta, 1.0) )
                  deltacands[ndeltacands++] = delta;

               if( intf0 < maxfrac )
               {
                  SCIP_Real delta2;

                  delta2 = 1.0 / (intscale * SCIPfloor(scip, maxfrac / intf0));

                  if( ! SCIPisEQ(scip, delta, delta2) && ! SCIPisEQ(scip, delta2, 1.0) )
                     deltacands[ndeltacands++] = delta2;
               }
            }
         }
      }
   }

   for( i = 0; i < nbounddist; ++i )
   {
      SCIP_Real absmksetcoef;

      absmksetcoef = REALABS(deltacands[i]);
      maxabsmksetcoef = MAX(absmksetcoef, maxabsmksetcoef);

      deltacands[i] = absmksetcoef;
   }

   /* also test 1.0 and maxabsmksetcoef + 1.0 as last delta values */
   if( maxabsmksetcoef != -1.0 )
      deltacands[ndeltacands++] = maxabsmksetcoef + 1.0;

   deltacands[ndeltacands++] = 1.0;

   maxtestdelta = MIN(ndeltacands, maxtestdelta);

   /* For each delta
    * Calculate fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j) , and derive MIR cut
    *   a~*x' <= down(b)
    * integers :  a~_j = down(a'_j)                      , if f_j <= f_0
    *             a~_j = down(a'_j) + (f_j - f0)/(1 - f0), if f_j >  f_0
    * continuous: a~_j = 0                               , if a'_j >= 0
    *             a~_j = a'_j/(1 - f0)                   , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */

   ntmpcoefs = 0;
   assert(mksetnnz <= nvars);
   for( i = intstart; i < mksetnnz; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real solval;
      SCIP_Real QUAD(coef);

      var = vars[mksetinds[i]];

      /* get the soltion value of the continuous variable */
      solval = SCIPgetSolVal(scip, sol, var);

      /* now compute the solution value in the transform space considering complementation */
      if( boundtype[i] == -1 )
      {
         /* variable was complemented with global (simple) bound */
         if( varsign[i] == -1 )
            solval = SCIPvarGetUbGlobal(var) - solval;
         else
            solval = solval - SCIPvarGetLbGlobal(var);
      }
      else
      {
         assert(boundtype[i] == -2);

         /* variable was complemented with local (simple) bound */
         if( varsign[i] == -1 )
            solval = SCIPvarGetUbLocal(var) - solval;
         else
            solval = solval - SCIPvarGetLbLocal(var);
      }

      tmpvalues[ntmpcoefs] = solval;
      QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[i]);
      tmpcoefs[ntmpcoefs] = varsign[i] * QUAD_TO_DBL(coef);
      ++ntmpcoefs;
   }

   assert(ntmpcoefs == mksetnnz - intstart);

   contactivity = 0.0;
   contsqrnorm = 0.0;
   for( i = 0; i < intstart; ++i )
   {
      SCIP_Real solval;
      SCIP_Real QUAD(mksetcoef);

      QUAD_ARRAY_LOAD(mksetcoef, mksetcoefs, mksetinds[i]);

      if( varsign[i] * QUAD_TO_DBL(mksetcoef) >= 0.0 )
         continue;

      /* get the soltion value of the continuous variable */
      solval = SCIPgetSolVal(scip, sol, vars[mksetinds[i]]);

      /* now compute the solution value in the transform space considering complementation */
      switch( boundtype[i] )
      {
      case -1:
         /* variable was complemented with global (simple) bound */
         if( varsign[i] == -1 )
            solval = SCIPvarGetUbGlobal(vars[mksetinds[i]]) - solval;
         else
            solval = solval - SCIPvarGetLbGlobal(vars[mksetinds[i]]);
         break;
      case -2:
         /* variable was complemented with local (simple) bound */
         if( varsign[i] == -1 )
            solval = SCIPvarGetUbLocal(vars[mksetinds[i]]) - solval;
         else
            solval = solval - SCIPvarGetLbLocal(vars[mksetinds[i]]);
         break;
      default:
         /* variable was complemented with a variable bound */
         if( varsign[i] == -1 )
         {
            SCIP_Real coef;
            SCIP_Real constant;
            SCIP_Real vbdsolval;

            coef = SCIPvarGetVubCoefs(vars[mksetinds[i]])[boundtype[i]];
            constant = SCIPvarGetVubConstants(vars[mksetinds[i]])[boundtype[i]];
            vbdsolval = SCIPgetSolVal(scip, sol, SCIPvarGetVubVars(vars[mksetinds[i]])[boundtype[i]]);

            solval = (coef * vbdsolval + constant) - solval;
         }
         else
         {
            SCIP_Real coef;
            SCIP_Real constant;
            SCIP_Real vbdsolval;

            coef = SCIPvarGetVlbCoefs(vars[mksetinds[i]])[boundtype[i]];
            constant = SCIPvarGetVlbConstants(vars[mksetinds[i]])[boundtype[i]];
            vbdsolval = SCIPgetSolVal(scip, sol, SCIPvarGetVlbVars(vars[mksetinds[i]])[boundtype[i]]);

            solval = solval - (coef * vbdsolval + constant);
         }
      }

      contactivity += solval * (QUAD_TO_DBL(mksetcoef) * varsign[i]);
      contsqrnorm += QUAD_TO_DBL(mksetcoef) * QUAD_TO_DBL(mksetcoef);
   }

   {
      SCIP_ROW** rows;

      rows = SCIPgetLPRows(scip);
      assert(ntmpcoefs <= nvars);
      for( i = 0; i < aggrrow->nrows; ++i )
      {
         SCIP_ROW* row;
         SCIP_Real slackval;

         row = rows[aggrrow->rowsinds[i]];

         if( (aggrrow->rowweights[i] * aggrrow->slacksign[i]) >= 0.0 && !row->integral )
            continue;

         /* compute solution value of slack variable */
         slackval = SCIPgetRowSolActivity(scip, row, sol);

         if( aggrrow->slacksign[i] == +1 )
         {
            /* right hand side */
            assert(!SCIPisInfinity(scip, row->rhs));

            slackval = row->rhs - slackval;
         }
         else
         {
            /* left hand side */
            assert(aggrrow->slacksign[i] == -1);
            assert(!SCIPisInfinity(scip, -row->lhs));

            slackval = slackval - row->lhs;
         }

         if( row->integral )
         {
            /* if row is integral add variable to tmp arrays */
            tmpvalues[ntmpcoefs] = slackval;
            tmpcoefs[ntmpcoefs] = aggrrow->rowweights[i] * aggrrow->slacksign[i];
            ++ntmpcoefs;
         }
         else
         {
            SCIP_Real slackcoeff = (aggrrow->rowweights[i] * aggrrow->slacksign[i]);

            /* otherwise add it to continuous activity */
            contactivity += slackval * slackcoeff;
            contsqrnorm += SQR(slackcoeff);
         }
      }
   }

   /* try all candidates for delta and remember best */
   bestdelta = SCIP_INVALID;
   bestefficacy = -SCIPinfinity(scip);

   for( i = 0; i < maxtestdelta; ++i )
   {
      int j;
      SCIP_Real efficacy;

      /* check if we have seen this value of delta before */
      SCIP_Bool deltaseenbefore = FALSE;
      for( j = 0; j < i; ++j )
      {
         if( SCIPisEQ(scip, deltacands[i], deltacands[j]) )
         {
            deltaseenbefore = TRUE;
            break;
         }
      }

      /* skip this delta value and allow one more delta value if available */
      if( deltaseenbefore )
      {
         maxtestdelta = MIN(maxtestdelta + 1, ndeltacands);
         continue;
      }

      efficacy = computeMIREfficacy(scip, tmpcoefs, tmpvalues, QUAD_TO_DBL(mksetrhs), contactivity, contsqrnorm, deltacands[i], ntmpcoefs, minfrac, maxfrac);

      if( efficacy > bestefficacy )
      {
         bestefficacy = efficacy;
         bestdelta = deltacands[i];
      }
   }

   /* no delta was found that yielded any cut */
   if( bestdelta == SCIP_INVALID ) /*lint !e777*/
      goto TERMINATE;

   /* try bestdelta divided by 2, 4 and 8 */
   {
      SCIP_Real basedelta = bestdelta;
      for( i = 2; i <= 8 ; i *= 2 )
      {
         SCIP_Real efficacy;
         SCIP_Real delta;

         delta = basedelta / i;

         efficacy = computeMIREfficacy(scip, tmpcoefs, tmpvalues, QUAD_TO_DBL(mksetrhs), contactivity, contsqrnorm, delta, ntmpcoefs, minfrac, maxfrac);

         if( efficacy > bestefficacy )
         {
            bestefficacy = efficacy;
            bestdelta = delta;
         }
      }
   }

   /* try to improve efficacy by switching complementation of integral variables that are not at their bounds
    * in order of non-increasing bound distance
    */
   for( i = 0; i < nbounddist; ++i )
   {
      int k;
      SCIP_Real newefficacy;
      SCIP_Real QUAD(newrhs);
      SCIP_Real QUAD(quadprod);
      SCIP_Real bestlb;
      SCIP_Real bestub;
      SCIP_Real oldsolval;
      int bestlbtype;
      int bestubtype;

      k = bounddistpos[i];

      SCIP_CALL( findBestLb(scip, vars[mksetinds[k]], sol, 0, allowlocal, &bestlb, &bestlbtype) );

      if( SCIPisInfinity(scip, -bestlb) )
         continue;

      SCIP_CALL( findBestUb(scip, vars[mksetinds[k]], sol, 0, allowlocal, &bestub, &bestubtype) );

      if( SCIPisInfinity(scip, bestub) )
         continue;

      /* switch the complementation of this variable */
#ifndef NDEBUG
      {
         SCIP_Real QUAD(coef);
         QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[k]);
         assert(SCIPisEQ(scip, tmpcoefs[k - intstart], varsign[k] * QUAD_TO_DBL(coef)));
      }
#endif

      /* compute this: newrhs = mksetrhs + tmpcoefs[k - intstart] * (bestlb - bestub); */
      SCIPquadprecProdDD(quadprod, tmpcoefs[k - intstart], bestlb - bestub);
      SCIPquadprecSumQQ(newrhs, mksetrhs, quadprod);
      tmpcoefs[k - intstart] = -tmpcoefs[k - intstart];

      oldsolval = tmpvalues[k - intstart];
      tmpvalues[k - intstart] = varsign[k] == +1 ? bestub - SCIPgetSolVal(scip, sol, vars[mksetinds[k]]) : SCIPgetSolVal(scip, sol, vars[mksetinds[k]]) - bestlb;

      /* compute new violation */
      newefficacy = computeMIREfficacy(scip, tmpcoefs, tmpvalues, QUAD_TO_DBL(newrhs), contactivity, contsqrnorm, bestdelta, ntmpcoefs, minfrac, maxfrac);

      /* check if violation was increased */
      if( newefficacy > bestefficacy )
      {
         /* keep change of complementation */
         bestefficacy = newefficacy;
         QUAD_ASSIGN_Q(mksetrhs, newrhs);

         if( varsign[k] == +1 )
         {
            /* switch to upper bound */
            assert(bestubtype < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
            boundtype[k] = bestubtype;
            varsign[k] = -1;
         }
         else
         {
            /* switch to lower bound */
            assert(bestlbtype < 0); /* cannot switch to a variable bound (would lead to further coef updates) */
            boundtype[k] = bestlbtype;
            varsign[k] = +1;
         }

         localbdsused = localbdsused || (boundtype[k] == -2);
      }
      else
      {
         /* undo the change of the complementation */
         tmpcoefs[k - intstart] = -tmpcoefs[k - intstart];
         tmpvalues[k - intstart] = oldsolval;
      }
   } /*lint !e438*/

   if( bestefficacy > 0.0 )
   {
      SCIP_Real mirefficacy;
      SCIP_Real QUAD(downrhs);
      SCIP_Real QUAD(f0);
      SCIP_Real scale;

      scale = 1.0 / bestdelta;
      SCIPquadprecProdQD(mksetrhs, mksetrhs, scale);
      SCIPquadprecEpsFloorQ(downrhs, mksetrhs, SCIPepsilon(scip)); /*lint !e666*/
      SCIPquadprecSumQQ(f0, mksetrhs, -downrhs);
      assert(QUAD_TO_DBL(f0) >= -SCIPepsilon(scip) && QUAD_TO_DBL(f0) < 1.0);

      /* renormalize f0 value */
      SCIPquadprecSumDD(f0, QUAD_HI(f0), QUAD_LO(f0));

      /* scale entries by the chosen scale factor */
      for( i = 0; i < mksetnnz; ++i )
      {
         SCIP_Real QUAD(coef);

         QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[i]);
         SCIPquadprecProdQD(coef, coef, scale);
         QUAD_ARRAY_STORE(mksetcoefs, mksetinds[i], coef);
      }
      SCIPdebugMsg(scip, "applied best scale (=%.13g):\n", scale);
      SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(mksetrhs), mksetinds, mksetnnz, FALSE, FALSE));

      QUAD_ASSIGN_Q(mksetrhs, downrhs);

      QUAD_ASSIGN_Q(data->cutrhs, mksetrhs);

      SCIP_CALL( cutsRoundMIR(scip, data, varsign, boundtype, QUAD(f0)) );

      SCIPdebugMsg(scip, "rounded MIR cut:\n");
      SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(data->cutrhs), mksetinds, data->ncutinds, FALSE, FALSE));

      /* substitute aggregated slack variables:
       *
       * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
       * variable only appears in its own row:
       *    a'_r = scale * weight[r] * slacksign[r].
       *
       * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
       *   integers :  a^_r = a~_r = down(a'_r)                      , if f_r <= f0
       *               a^_r = a~_r = down(a'_r) + (f_r - f0)/(1 - f0), if f_r >  f0
       *   continuous: a^_r = a~_r = 0                               , if a'_r >= 0
       *               a^_r = a~_r = a'_r/(1 - f0)                   , if a'_r <  0
       *
       * Substitute a^_r * s_r by adding a^_r times the slack's definition to the cut.
       */
      SCIP_CALL( cutsSubstituteMIR(scip, aggrrow->rowweights, aggrrow->slacksign, aggrrow->rowsinds,
            aggrrow->nrows, scale, mksetcoefs, QUAD(&data->cutrhs), mksetinds, &data->ncutinds, QUAD(f0)) );

      SCIPdebugMsg(scip, "substituted slacks in MIR cut:\n");
      SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(data->cutrhs), mksetinds, data->ncutinds, FALSE, FALSE));

#ifndef NDEBUG
      {
         SCIP_Real efficacy = -QUAD_TO_DBL(data->cutrhs);
         for( i = 0; i < data->ncutinds; ++i )
         {
            SCIP_Real QUAD(coef);
            QUAD_ARRAY_LOAD(coef, mksetcoefs, mksetinds[i]);
            efficacy += QUAD_TO_DBL(coef) * SCIPgetSolVal(scip, sol, vars[mksetinds[i]]);
         }

         if( !EPSZ(SCIPrelDiff(efficacy, bestefficacy), 1e-4) )
         {
            SCIPdebugMsg(scip, "efficacy of cmir cut is different than expected efficacy: %f != %f\n", efficacy, bestefficacy);
         }
      }
#endif

      *cutislocal = *cutislocal || localbdsused;

      /* remove all nearly-zero coefficients from MIR row and relax the right hand side correspondingly in order to
       * prevent numerical rounding errors
       */
      if( postprocess )
      {
         SCIP_CALL( postprocessCutQuad(scip, *cutislocal, mksetinds, mksetcoefs, &data->ncutinds, QUAD(&data->cutrhs), success) );
      }
      else
      {
         *success = ! removeZerosQuad(scip, SCIPsumepsilon(scip), *cutislocal, mksetcoefs, QUAD(&data->cutrhs), mksetinds, &data->ncutinds);
      }

      SCIPdebugMsg(scip, "post-processed cut (success = %s):\n", *success ? "TRUE" : "FALSE");
      SCIPdebug(printCutQuad(scip, sol, mksetcoefs, QUAD(data->cutrhs), mksetinds, data->ncutinds, FALSE, FALSE));

      if( *success )
      {
         mirefficacy = calcEfficacyDenseStorageQuad(scip, sol, mksetcoefs, QUAD_TO_DBL(data->cutrhs), mksetinds, data->ncutinds);

         if( SCIPisEfficacious(scip, mirefficacy) && SCIPisGT(scip, mirefficacy, *cutefficacy) )
         {
            BMScopyMemoryArray(cutinds, mksetinds, data->ncutinds);
            for( i = 0; i < data->ncutinds; ++i )
            {
               SCIP_Real QUAD(coef);
               int j = cutinds[i];

               QUAD_ARRAY_LOAD(coef, mksetcoefs, j);

               cutcoefs[i] = QUAD_TO_DBL(coef);
               QUAD_ASSIGN(coef, 0.0);
               QUAD_ARRAY_STORE(mksetcoefs, j, coef);
            }
            *cutnnz = data->ncutinds;
            *cutrhs = QUAD_TO_DBL(data->cutrhs);
            *cutefficacy = mirefficacy;
            if( cutrank != NULL )
               *cutrank = aggrrow->rank + 1;
            *cutislocal = *cutislocal || localbdsused;
         }
         else
            *success = FALSE;
      }
   }

  TERMINATE:
   /* if we aborted early we need to clean the mksetcoefs */
   if( !(*success) )
   {
      SCIP_Real QUAD(tmp);
      QUAD_ASSIGN(tmp, 0.0);

      for( i = 0; i < data->ncutinds; ++i )
      {
         QUAD_ARRAY_STORE(mksetcoefs, mksetinds[i], tmp);
      }
   }

#ifndef NDEBUG
   for( i = 0; i < QUAD_ARRAY_SIZE(nvars); ++i )
   {
      if(mksetcoefs[i] != 0.0)
      {
         SCIPdebugMsg(scip, "mksetcoefs have not been reset\n");
         SCIPABORT();
      }
   }
#endif

   if( data->cutinds != NULL )
      SCIPfreeBufferArray(scip, &data->cutinds);

   if( data->cutcoefs != NULL )
      SCIPfreeCleanBufferArray(scip, &data->cutcoefs);

   for( int s = NSECTIONS - 1; s >= 0; --s )
   {
      SCIPfreeBufferArray(scip, &data->secindices[s]);
   }

   SCIPfreeBuffer(scip, &data);
   /* free temporary memory */
   SCIPfreeBufferArray(scip, &bounddistpos);
   SCIPfreeBufferArray(scip, &bounddist);
   SCIPfreeBufferArray(scip, &deltacands);
   SCIPfreeBufferArray(scip, &tmpvalues);
   SCIPfreeBufferArray(scip, &tmpcoefs);
   SCIPfreeBufferArray(scip, &boundtype);
   SCIPfreeBufferArray(scip, &varsign);

   return SCIP_OKAY;
}

/* =========================================== flow cover =========================================== */

#define NO_EXACT_KNAPSACK

#ifndef NO_EXACT_KNAPSACK
#define MAXDNOM                  1000LL
#define MINDELTA                  1e-03
#define MAXDELTA                  1e-09
#define MAXSCALE                 1000.0
#define MAXDYNPROGSPACE         1000000
#endif

#define MAXABSVBCOEF               1e+5 /**< maximal absolute coefficient in variable bounds used for snf relaxation */
#define MAXBOUND                  1e+10   /**< maximal value of normal bounds used for snf relaxation */

/** structure that contains all data required to perform the sequence independent lifting
 */
typedef
struct LiftingData
{
   SCIP_Real*            M;                  /**< \f$ M_0 := 0.0 \f$ and \f$ M_i := M_i-1 + m_i \f$ */
   SCIP_Real*            m;                  /**< non-increasing array of variable upper bound coefficients
                                              *   for all variables in \f$ C^{++} \f$  and \f$ L^- \f$,
                                              *   where \f$ C = C^+ \cup C^- \f$ is the flowcover and
                                              *   \f$ C^{++} := \{ j \in C^+ \mid u_j > \lambda \} \f$
                                              *   \f$ L^- := \{ j \in (N^- \setminus C^-) \mid u_j > \lambda \} \f$
                                              */
   int                   r;                  /**< size of array m */
   int                   t;                  /**< index of smallest value in m that comes from a variable in \f$ C^{++} \f$ */
   SCIP_Real             d1;                 /**< right hand side of single-node-flow set plus the sum of all \f$ u_j \f$ for \f$ j \in C^- \f$ */
   SCIP_Real             d2;                 /**< right hand side of single-node-flow set plus the sum of all \f$ u_j \f$ for \f$ j \in N^- \f$ */
   SCIP_Real             lambda;             /**< excess of the flowcover */
   SCIP_Real             mp;                 /**< smallest variable bound coefficient of variable in \f$ C^{++} (min_{j \in C++} u_j) \f$ */
   SCIP_Real             ml;                 /**< \f$ ml := min(\lambda, \sum_{j \in C^+ \setminus C^{++}} u_j) \f$ */
} LIFTINGDATA;

/** structure that contains all the data that defines the single-node-flow relaxation of an aggregation row */
typedef
struct SNF_Relaxation
{
   int*                  transvarcoefs;      /**< coefficients of all vars in relaxed set */
   SCIP_Real*            transbinvarsolvals; /**< sol val of bin var in vub of all vars in relaxed set */
   SCIP_Real*            transcontvarsolvals;/**< sol val of all real vars in relaxed set */
   SCIP_Real*            transvarvubcoefs;   /**< coefficient in vub of all vars in relaxed set */
   int                   ntransvars;         /**< number of vars in relaxed set */
   SCIP_Real             transrhs;           /**< rhs in relaxed set */
   int*                  origbinvars;        /**< associated original binary var for all vars in relaxed set */
   int*                  origcontvars;       /**< associated original continuous var for all vars in relaxed set */
   SCIP_Real*            aggrcoefsbin;       /**< aggregation coefficient of the original binary var used to define the
                                              *   continuous variable in the relaxed set */
   SCIP_Real*            aggrcoefscont;      /**< aggregation coefficient of the original continuous var used to define the
                                              *   continuous variable in the relaxed set */
   SCIP_Real*            aggrconstants;      /**< aggregation constant used to define the continuous variable in the relaxed set */
} SNF_RELAXATION;

/** get solution value and index of variable lower bound (with binary variable) which is closest to the current LP
 *  solution value of a given variable; candidates have to meet certain criteria in order to ensure the nonnegativity
 *  of the variable upper bound imposed on the real variable in the 0-1 single node flow relaxation associated with the
 *  given variable
 */
static
SCIP_RETCODE getClosestVlb(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< given active problem variable */
   SCIP_SOL*             sol,                /**< solution to use for variable bound; NULL for LP solution */
   SCIP_Real*            rowcoefs,           /**< (dense) array of coefficients of row */
   int8_t*               binvarused,         /**< array that stores if a binary variable was already used (+1)
                                              *   was not used (0) or was not used but is contained in the row (-1) */
   SCIP_Real             bestsub,            /**< closest simple upper bound of given variable */
   SCIP_Real             rowcoef,            /**< coefficient of given variable in current row */
   SCIP_Real*            closestvlb,         /**< pointer to store the LP sol value of the closest variable lower bound */
   int*                  closestvlbidx       /**< pointer to store the index of the closest vlb; -1 if no vlb was found */
   )
{
   int nvlbs;
   int nbinvars;

   assert(scip != NULL);
   assert(var != NULL);
   assert(bestsub == SCIPvarGetUbGlobal(var) || bestsub == SCIPvarGetUbLocal(var)); /*lint !e777*/
   assert(!SCIPisInfinity(scip, bestsub));
   assert(!EPSZ(rowcoef, QUAD_EPSILON));
   assert(rowcoefs != NULL);
   assert(binvarused != NULL);
   assert(closestvlb != NULL);
   assert(closestvlbidx != NULL);

   nvlbs = SCIPvarGetNVlbs(var);
   nbinvars = SCIPgetNBinVars(scip);

   *closestvlbidx = -1;
   *closestvlb = -SCIPinfinity(scip);
   if( nvlbs > 0 )
   {
      SCIP_VAR** vlbvars;
      SCIP_Real* vlbcoefs;
      SCIP_Real* vlbconsts;
      int i;

      vlbvars = SCIPvarGetVlbVars(var);
      vlbcoefs = SCIPvarGetVlbCoefs(var);
      vlbconsts = SCIPvarGetVlbConstants(var);

      for( i = 0; i < nvlbs; i++ )
      {
         SCIP_Real rowcoefbinvar;
         SCIP_Real val1;
         SCIP_Real val2;
         SCIP_Real vlbsol;
         SCIP_Real rowcoefsign;
         int probidxbinvar;

         if( bestsub > vlbconsts[i] )
            continue;

         /* for numerical reasons, ignore variable bounds with large absolute coefficient and
          * those which lead to an infinite variable bound coefficient (val2) in snf relaxation
          */
         if( REALABS(vlbcoefs[i]) > MAXABSVBCOEF  )
            continue;

         /* use only variable lower bounds l~_i * x_i + d_i with x_i binary which are active */
         probidxbinvar = SCIPvarGetProbindex(vlbvars[i]);

         /* if the variable is not active the problem index is -1, so we cast to unsigned int before the comparison which
          * ensures that the problem index is between 0 and nbinvars - 1
          */
         if( (unsigned int)probidxbinvar >= (unsigned int)nbinvars )
            continue;

         assert(SCIPvarIsBinary(vlbvars[i]));

         /* check if current variable lower bound l~_i * x_i + d_i imposed on y_j meets the following criteria:
          * (let a_j  = coefficient of y_j in current row,
          *      u_j  = closest simple upper bound imposed on y_j,
          *      c_i  = coefficient of x_i in current row)
          *   0. no other non-binary variable y_k has used a variable bound with x_i to get transformed variable y'_k yet
          * if a_j > 0:
          *   1. u_j <= d_i
          *   2. a_j ( u_j - d_i ) + c_i <= 0
          *   3. a_j l~_i + c_i <= 0
          * if a_j < 0:
          *   1. u_j <= d_i
          *   2. a_j ( u_j - d_i ) + c_i >= 0
          *   3. a_j l~_i + c_i >= 0
          */

         /* has already been used in the SNF relaxation */
         if( binvarused[probidxbinvar] == 1 )
            continue;

         /* get the row coefficient */
         {
            SCIP_Real QUAD(tmp);
            QUAD_ARRAY_LOAD(tmp, rowcoefs, probidxbinvar);
            rowcoefbinvar = QUAD_TO_DBL(tmp);
         }
         rowcoefsign = COPYSIGN(1.0, rowcoef);

         val2 = rowcoefsign * ((rowcoef * vlbcoefs[i]) + rowcoefbinvar);

         /* variable lower bound does not meet criteria */
         if( val2 > 0.0 || SCIPisInfinity(scip, -val2) )
            continue;

         val1 = rowcoefsign * ((rowcoef * (bestsub - vlbconsts[i])) + rowcoefbinvar);

         /* variable lower bound does not meet criteria */
         if( val1 > 0.0 )
            continue;

         vlbsol = vlbcoefs[i] * SCIPgetSolVal(scip, sol, vlbvars[i]) + vlbconsts[i];
         if( vlbsol > *closestvlb )
         {
            *closestvlb = vlbsol;
            *closestvlbidx = i;
         }
         assert(*closestvlbidx >= 0);
      }
   }

   return SCIP_OKAY;
}

/** get LP solution value and index of variable upper bound (with binary variable) which is closest to the current LP
 *  solution value of a given variable; candidates have to meet certain criteria in order to ensure the nonnegativity
 *  of the variable upper bound imposed on the real variable in the 0-1 single node flow relaxation associated with the
 *  given variable
 */
static
SCIP_RETCODE getClosestVub(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< given active problem variable */
   SCIP_SOL*             sol,                /**< solution to use for variable bound; NULL for LP solution */
   SCIP_Real*            rowcoefs,           /**< (dense) array of coefficients of row */
   int8_t*               binvarused,         /**< array that stores if a binary variable was already used (+1)
                                              *   was not used (0) or was not used but is contained in the row (-1)
                                              */
   SCIP_Real             bestslb,            /**< closest simple lower bound of given variable */
   SCIP_Real             rowcoef,            /**< coefficient of given variable in current row */
   SCIP_Real*            closestvub,         /**< pointer to store the LP sol value of the closest variable upper bound */
   int*                  closestvubidx       /**< pointer to store the index of the closest vub; -1 if no vub was found */
   )
{
   int nvubs;
   int nbinvars;

   assert(scip != NULL);
   assert(var != NULL);
   assert(bestslb == SCIPvarGetLbGlobal(var) || bestslb == SCIPvarGetLbLocal(var)); /*lint !e777*/
   assert(!SCIPisInfinity(scip, - bestslb));
   assert(!EPSZ(rowcoef, QUAD_EPSILON));
   assert(rowcoefs != NULL);
   assert(binvarused != NULL);
   assert(closestvub != NULL);
   assert(closestvubidx != NULL);

   nvubs = SCIPvarGetNVubs(var);
   nbinvars = SCIPgetNBinVars(scip);

   *closestvubidx = -1;
   *closestvub = SCIPinfinity(scip);
   if( nvubs > 0 )
   {
      SCIP_VAR** vubvars;
      SCIP_Real* vubcoefs;
      SCIP_Real* vubconsts;
      int i;

      vubvars = SCIPvarGetVubVars(var);
      vubcoefs = SCIPvarGetVubCoefs(var);
      vubconsts = SCIPvarGetVubConstants(var);

      for( i = 0; i < nvubs; i++ )
      {
         SCIP_Real rowcoefbinvar;
         SCIP_Real val1;
         SCIP_Real val2;
         SCIP_Real vubsol;
         SCIP_Real rowcoefsign;
         int probidxbinvar;

         if( bestslb < vubconsts[i] )
            continue;

         /* for numerical reasons, ignore variable bounds with large absolute coefficient and
          * those which lead to an infinite variable bound coefficient (val2) in snf relaxation
          */
         if( REALABS(vubcoefs[i]) > MAXABSVBCOEF  )
            continue;

         /* use only variable upper bound u~_i * x_i + d_i with x_i binary and which are active */
         probidxbinvar = SCIPvarGetProbindex(vubvars[i]);

         /* if the variable is not active the problem index is -1, so we cast to unsigned int before the comparison which
          * ensures that the problem index is between 0 and nbinvars - 1
          */
         if( (unsigned int)probidxbinvar >= (unsigned int)nbinvars )
            continue;

         assert(SCIPvarIsBinary(vubvars[i]));

         /* checks if current variable upper bound u~_i * x_i + d_i meets the following criteria
          * (let a_j  = coefficient of y_j in current row,
          *      l_j  = closest simple lower bound imposed on y_j,
          *      c_i  = coefficient of x_i in current row)
          *   0. no other non-binary variable y_k has used a variable bound with x_i to get transformed variable y'_k
          * if a > 0:
          *   1. l_j >= d_i
          *   2. a_j ( l_i - d_i ) + c_i >= 0
          *   3. a_j u~_i + c_i >= 0
          * if a < 0:
          *   1. l_j >= d_i
          *   2. a_j ( l_j - d_i ) + c_i <= 0
          *   3. a_j u~_i + c_i <= 0
          */

         /* has already been used in the SNF relaxation */
         if( binvarused[probidxbinvar] == 1 )
            continue;

         /* get the row coefficient */
         {
            SCIP_Real QUAD(tmp);
            QUAD_ARRAY_LOAD(tmp, rowcoefs, probidxbinvar);
            rowcoefbinvar = QUAD_TO_DBL(tmp);
         }
         rowcoefsign = COPYSIGN(1.0, rowcoef);

         val2 = rowcoefsign * ((rowcoef * vubcoefs[i]) + rowcoefbinvar);

         /* variable upper bound does not meet criteria */
         if( val2 < 0.0 || SCIPisInfinity(scip, val2) )
            continue;

         val1 = rowcoefsign * ((rowcoef * (bestslb - vubconsts[i])) + rowcoefbinvar);

         /* variable upper bound does not meet criteria */
         if( val1 < 0.0 )
            continue;

         vubsol = vubcoefs[i] * SCIPgetSolVal(scip, sol, vubvars[i]) + vubconsts[i];
         if( vubsol < *closestvub )
         {
            *closestvub = vubsol;
            *closestvubidx = i;
         }
         assert(*closestvubidx >= 0);
      }
   }

   return SCIP_OKAY;
}

/** determines the bounds to use for constructing the single-node-flow relaxation of a variable in
 *  the given row.
 */
static
SCIP_RETCODE determineBoundForSNF(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to use for variable bound; NULL for LP solution */
   SCIP_VAR**            vars,               /**< array of problem variables */
   SCIP_Real*            rowcoefs,           /**< (dense) array of variable coefficients in the row */
   int*                  rowinds,            /**< array with positions of non-zero values in the rowcoefs array */
   int                   varposinrow,        /**< position of variable in the rowinds array for which the bounds should be determined */
   int8_t*               binvarused,         /**< array that stores if a binary variable was already used (+1)
                                              *   was not used (0) or was not used but is contained in the row (-1)
                                              */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Real*            bestlb,             /**< pointer to store best lower bound for transformation */
   SCIP_Real*            bestub,             /**< pointer to store best upper bound for transformation */
   SCIP_Real*            bestslb,            /**< pointer to store best simple lower bound for transformation */
   SCIP_Real*            bestsub,            /**< pointer to store best simple upper bound for transformation */
   int*                  bestlbtype,         /**< pointer to store type of best lower bound (-2: local bound, -1: global bound, >= 0 variable bound index) */
   int*                  bestubtype,         /**< pointer to store type of best upper bound (-2: local bound, -1: global bound, >= 0 variable bound index) */
   int*                  bestslbtype,        /**< pointer to store type of best simple lower bound */
   int*                  bestsubtype,        /**< pointer to store type of best simple upper bound */
   SCIP_BOUNDTYPE*       selectedbounds,     /**< pointer to store the preferred bound for the transformation */
   SCIP_Bool*            freevariable        /**< pointer to store if variable is a free variable */
   )
{
   SCIP_VAR* var;

   SCIP_Real rowcoef;
   SCIP_Real solval;

   int probidx;

   bestlb[varposinrow] = -SCIPinfinity(scip);
   bestub[varposinrow] = SCIPinfinity(scip);
   bestlbtype[varposinrow] = -3;
   bestubtype[varposinrow] = -3;

   probidx = rowinds[varposinrow];
   var = vars[probidx];
   {
      SCIP_Real QUAD(tmp);
      QUAD_ARRAY_LOAD(tmp, rowcoefs, probidx);
      rowcoef = QUAD_TO_DBL(tmp);
   }

   assert(!EPSZ(rowcoef, QUAD_EPSILON));

   /* get closest simple lower bound and closest simple upper bound */
   SCIP_CALL( findBestLb(scip, var, sol, 0, allowlocal, &bestslb[varposinrow], &bestslbtype[varposinrow]) );
   SCIP_CALL( findBestUb(scip, var, sol, 0, allowlocal, &bestsub[varposinrow], &bestsubtype[varposinrow]) );

   /* do not use too large bounds */
   if( bestslb[varposinrow] <= -MAXBOUND )
      bestslb[varposinrow] = -SCIPinfinity(scip);

   if( bestsub[varposinrow] >= MAXBOUND )
      bestsub[varposinrow] = SCIPinfinity(scip);

   solval = SCIPgetSolVal(scip, sol, var);

   SCIPdebugMsg(scip, "  %d: %g <%s, idx=%d, lp=%g, [%g(%d),%g(%d)]>:\n", varposinrow, rowcoef, SCIPvarGetName(var), probidx,
      solval, bestslb[varposinrow], bestslbtype[varposinrow], bestsub[varposinrow], bestsubtype[varposinrow]);

   /* mixed integer set cannot be relaxed to 0-1 single node flow set because both simple bounds are -infinity
    * and infinity, respectively
    */
   if( SCIPisInfinity(scip, -bestslb[varposinrow]) && SCIPisInfinity(scip, bestsub[varposinrow]) )
   {
      *freevariable = TRUE;
      return SCIP_OKAY;
   }

   /* get closest lower bound that can be used to define the real variable y'_j in the 0-1 single node flow
    * relaxation
    */
   if( !SCIPisInfinity(scip, bestsub[varposinrow]) )
   {
      bestlb[varposinrow] = bestslb[varposinrow];
      bestlbtype[varposinrow] = bestslbtype[varposinrow];

      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvlb;
         int bestvlbidx;

         SCIP_CALL( getClosestVlb(scip, var, sol, rowcoefs, binvarused, bestsub[varposinrow], rowcoef, &bestvlb, &bestvlbidx) );
         if( SCIPisGT(scip, bestvlb, bestlb[varposinrow]) )
         {
            bestlb[varposinrow] = bestvlb;
            bestlbtype[varposinrow] = bestvlbidx;
         }
      }
   }

   /* get closest upper bound that can be used to define the real variable y'_j in the 0-1 single node flow
    * relaxation
    */
   if( !SCIPisInfinity(scip, -bestslb[varposinrow]) )
   {
      bestub[varposinrow] = bestsub[varposinrow];
      bestubtype[varposinrow] = bestsubtype[varposinrow];

      if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      {
         SCIP_Real bestvub;
         int bestvubidx;

         SCIP_CALL( getClosestVub(scip, var, sol, rowcoefs, binvarused, bestslb[varposinrow], rowcoef, &bestvub, &bestvubidx) );
         if( SCIPisLT(scip, bestvub, bestub[varposinrow]) )
         {
            bestub[varposinrow] = bestvub;
            bestubtype[varposinrow] = bestvubidx;
         }
      }
   }
   SCIPdebugMsg(scip, "        bestlb=%g(%d), bestub=%g(%d)\n", bestlb[varposinrow], bestlbtype[varposinrow], bestub[varposinrow], bestubtype[varposinrow]);

   /* mixed integer set cannot be relaxed to 0-1 single node flow set because there are no suitable bounds
    * to define the transformed variable y'_j
    */
   if( SCIPisInfinity(scip, -bestlb[varposinrow]) && SCIPisInfinity(scip, bestub[varposinrow]) )
   {
      *freevariable = TRUE;
      return SCIP_OKAY;
   }

   *freevariable = FALSE;

   /* select best upper bound if it is closer to the LP value of y_j and best lower bound otherwise and use this bound
    * to define the real variable y'_j with 0 <= y'_j <= u'_j x_j in the 0-1 single node flow relaxation;
    * prefer variable bounds
    */
   if( SCIPisEQ(scip, solval, (1.0 - boundswitch) * bestlb[varposinrow] + boundswitch * bestub[varposinrow]) && bestlbtype[varposinrow] >= 0 )
   {
      selectedbounds[varposinrow] = SCIP_BOUNDTYPE_LOWER;
   }
   else if( SCIPisEQ(scip, solval, (1.0 - boundswitch) * bestlb[varposinrow] + boundswitch * bestub[varposinrow])
      && bestubtype[varposinrow] >= 0 )
   {
      selectedbounds[varposinrow] = SCIP_BOUNDTYPE_UPPER;
   }
   else if( SCIPisLE(scip, solval, (1.0 - boundswitch) * bestlb[varposinrow] + boundswitch * bestub[varposinrow]) )
   {
      selectedbounds[varposinrow] = SCIP_BOUNDTYPE_LOWER;
   }
   else
   {
      assert(SCIPisGT(scip, solval, (1.0 - boundswitch) * bestlb[varposinrow] + boundswitch * bestub[varposinrow]));
      selectedbounds[varposinrow] = SCIP_BOUNDTYPE_UPPER;
   }

   if( selectedbounds[varposinrow] == SCIP_BOUNDTYPE_LOWER && bestlbtype[varposinrow] >= 0 )
   {
      int vlbvarprobidx;
      SCIP_VAR** vlbvars = SCIPvarGetVlbVars(var);

       /* mark binary variable of vlb so that it is not used for other continuous variables
       * by setting it's position in the aggrrow to a negative value
       */
      vlbvarprobidx = SCIPvarGetProbindex(vlbvars[bestlbtype[varposinrow]]);
      binvarused[vlbvarprobidx] = 1;
   }
   else if( selectedbounds[varposinrow] == SCIP_BOUNDTYPE_UPPER && bestubtype[varposinrow] >= 0 )
   {
      int vubvarprobidx;
      SCIP_VAR** vubvars = SCIPvarGetVubVars(var);

       /* mark binary variable of vub so that it is not used for other continuous variables
       * by setting it's position in the aggrrow to a negative value
       */
      vubvarprobidx = SCIPvarGetProbindex(vubvars[bestubtype[varposinrow]]);
      binvarused[vubvarprobidx] = 1;
   }

   return SCIP_OKAY; /*lint !e438*/
}

/** construct a 0-1 single node flow relaxation (with some additional simple constraints) of a mixed integer set
 *  corresponding to the given aggrrow a * x <= rhs
 */
static
SCIP_RETCODE constructSNFRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            rowcoefs,           /**< array of coefficients of row */
   QUAD(SCIP_Real        rowrhs),            /**< pointer to right hand side of row */
   int*                  rowinds,            /**< array of variables problem indices for non-zero coefficients in row */
   int                   nnz,                /**< number of non-zeros in row */
   SNF_RELAXATION*       snf,                /**< stores the sign of the transformed variable in summation */
   SCIP_Bool*            success,            /**< stores whether the transformation was valid */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_VAR** vars;
   int i;
   int nnonbinvarsrow;
   int8_t* binvarused;
   int nbinvars;
   SCIP_Real QUAD(transrhs);

   /* arrays to store the selected bound for each non-binary variable in the row */
   SCIP_Real* bestlb;
   SCIP_Real* bestub;
   SCIP_Real* bestslb;
   SCIP_Real* bestsub;
   int* bestlbtype;
   int* bestubtype;
   int* bestslbtype;
   int* bestsubtype;
   SCIP_BOUNDTYPE* selectedbounds;

   *success = FALSE;

   SCIPdebugMsg(scip, "--------------------- construction of SNF relaxation ------------------------------------\n");

   nbinvars = SCIPgetNBinVars(scip);
   vars = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &bestlb, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestub, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestslb, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestsub, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbtype, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubtype, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestslbtype, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestsubtype, nnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selectedbounds, nnz) );

   /* sort descending to have continuous variables first */
   SCIPsortDownInt(rowinds, nnz);

   /* array to store whether a binary variable is in the row (-1) or has been used (1) due to variable bound usage */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &binvarused, nbinvars) );

   for( i = nnz - 1; i >= 0 && rowinds[i] < nbinvars; --i )
      binvarused[rowinds[i]] = -1;

   nnonbinvarsrow = i + 1;
   /* determine the bounds to use for transforming the non-binary variables */
   for( i = 0; i < nnonbinvarsrow; ++i )
   {
      SCIP_Bool freevariable;

      assert(rowinds[i] >= nbinvars);

      SCIP_CALL( determineBoundForSNF(scip, sol, vars, rowcoefs, rowinds, i, binvarused, allowlocal, boundswitch,
            bestlb, bestub, bestslb, bestsub, bestlbtype, bestubtype, bestslbtype, bestsubtype, selectedbounds, &freevariable) );

      if( freevariable )
      {
         int j;

         /* clear binvarused at indices of binary variables of row */
         for( j = nnz - 1; j >= nnonbinvarsrow; --j )
            binvarused[rowinds[j]] = 0;

         /* clear binvarused at indices of selected variable bounds */
         for( j = 0; j < i; ++j )
         {
            if( selectedbounds[j] == SCIP_BOUNDTYPE_LOWER && bestlbtype[j] >= 0 )
            {
               SCIP_VAR** vlbvars = SCIPvarGetVlbVars(vars[rowinds[j]]);
               binvarused[SCIPvarGetProbindex(vlbvars[bestlbtype[j]])] = 0;
            }
            else if( selectedbounds[j] == SCIP_BOUNDTYPE_UPPER && bestubtype[j] >= 0 )
            {
               SCIP_VAR** vubvars = SCIPvarGetVubVars(vars[rowinds[j]]);
               binvarused[SCIPvarGetProbindex(vubvars[bestubtype[j]])] = 0;
            }
         }

         /* terminate */
         goto TERMINATE;
      }
   }

   *localbdsused = FALSE;
   QUAD_ASSIGN_Q(transrhs, rowrhs);
   snf->ntransvars = 0;

   assert(snf->transvarcoefs != NULL); /* for lint */
   assert(snf->transvarvubcoefs != NULL);
   assert(snf->transbinvarsolvals != NULL);
   assert(snf->transcontvarsolvals != NULL);
   assert(snf->aggrconstants != NULL);
   assert(snf->aggrcoefscont != NULL);
   assert(snf->origcontvars != NULL);
   assert(snf->origbinvars != NULL);
   assert(snf->aggrcoefsbin != NULL);

   /* transform non-binary variables */
   for( i = 0; i < nnonbinvarsrow; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real QUAD(rowcoef);
      SCIP_Real solval;
      int probidx;

      probidx = rowinds[i];
      var = vars[probidx];
      QUAD_ARRAY_LOAD(rowcoef, rowcoefs, probidx);
      assert(!EPSZ(QUAD_TO_DBL(rowcoef), QUAD_EPSILON));
      solval = SCIPgetSolVal(scip, sol, var);

      assert(probidx >= nbinvars);

      if( selectedbounds[i] == SCIP_BOUNDTYPE_LOWER )
      {
         /* use bestlb to define y'_j */

         assert(!SCIPisInfinity(scip, bestsub[i]));
         assert(!SCIPisInfinity(scip, - bestlb[i]));
         assert(bestsubtype[i] == -1 || bestsubtype[i] == -2);
         assert(bestlbtype[i] > -3 && bestlbtype[i] < SCIPvarGetNVlbs(var));

         /* store for y_j that bestlb is the bound used to define y'_j and that y'_j is the associated real variable
          * in the relaxed set
          */
         snf->origcontvars[snf->ntransvars] = probidx;

         if( bestlbtype[i] < 0 )
         {
            SCIP_Real QUAD(val);
            SCIP_Real QUAD(contsolval);
            SCIP_Real QUAD(rowcoeftimesbestsub);

            /* use simple lower bound in bestlb = l_j <= y_j <= u_j = bestsub to define
             *   y'_j = - a_j ( y_j - u_j ) with 0 <= y'_j <=   a_j ( u_j - l_j ) x_j and x_j = 1    if a_j > 0
             *   y'_j =   a_j ( y_j - u_j ) with 0 <= y'_j <= - a_j ( u_j - l_j ) x_j and x_j = 1    if a_j < 0,
             * put j into the set
             *   N2   if a_j > 0
             *   N1   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j u_j
             */
            SCIPquadprecSumDD(val, bestsub[i], -bestlb[i]);
            SCIPquadprecProdQQ(val, val, rowcoef);
            SCIPquadprecSumDD(contsolval, solval, -bestsub[i]);
            SCIPquadprecProdQQ(contsolval, contsolval, rowcoef);

            if( bestlbtype[i] == -2 || bestsubtype[i] == -2 )
               *localbdsused = TRUE;

            SCIPquadprecProdQD(rowcoeftimesbestsub, rowcoef, bestsub[i]);

            /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
            snf->origbinvars[snf->ntransvars] = -1;
            snf->aggrcoefsbin[snf->ntransvars] = 0.0;

            if( QUAD_TO_DBL(rowcoef) >= 0.0 )
            {
               snf->transvarcoefs[snf->ntransvars] = - 1;
               snf->transvarvubcoefs[snf->ntransvars] = QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = 1.0;
               snf->transcontvarsolvals[snf->ntransvars] = - QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrconstants[snf->ntransvars] = QUAD_TO_DBL(rowcoeftimesbestsub);
               snf->aggrcoefscont[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
            }
            else
            {
               snf->transvarcoefs[snf->ntransvars] = 1;
               snf->transvarvubcoefs[snf->ntransvars] = - QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = 1.0;
               snf->transcontvarsolvals[snf->ntransvars] = QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrconstants[snf->ntransvars] = - QUAD_TO_DBL(rowcoeftimesbestsub);
               snf->aggrcoefscont[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
            }
            SCIPquadprecSumQQ(transrhs, transrhs, -rowcoeftimesbestsub);

            SCIPdebugMsg(scip, "    --> bestlb used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=1), rhs=%g-(%g*%g)=%g\n",
               snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars, snf->transvarvubcoefs[snf->ntransvars],
               snf->ntransvars, QUAD_TO_DBL(transrhs) + QUAD_TO_DBL(rowcoeftimesbestsub), QUAD_TO_DBL(rowcoef), bestsub[i], QUAD_TO_DBL(transrhs));
         }
         else
         {
            SCIP_Real QUAD(rowcoefbinary);
            SCIP_Real varsolvalbinary;
            SCIP_Real QUAD(val);
            SCIP_Real QUAD(contsolval);
            SCIP_Real QUAD(rowcoeftimesvlbconst);
            int vlbvarprobidx;

            SCIP_VAR** vlbvars = SCIPvarGetVlbVars(var);
            SCIP_Real* vlbconsts = SCIPvarGetVlbConstants(var);
            SCIP_Real* vlbcoefs = SCIPvarGetVlbCoefs(var);

            /* use variable lower bound in bestlb = l~_j x_j + d_j <= y_j <= u_j = bestsub to define
             *   y'_j = - ( a_j ( y_j - d_j ) + c_j x_j ) with 0 <= y'_j <= - ( a_j l~_j + c_j ) x_j    if a_j > 0
             *   y'_j =     a_j ( y_j - d_j ) + c_j x_j   with 0 <= y'_j <=   ( a_j l~_j + c_j ) x_j    if a_j < 0,
             * where c_j is the coefficient of x_j in the row, put j into the set
             *   N2   if a_j > 0
             *   N1   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j d_j
             */

            vlbvarprobidx = SCIPvarGetProbindex(vlbvars[bestlbtype[i]]);
            assert(binvarused[vlbvarprobidx] == 1);
            assert(vlbvarprobidx < nbinvars);

            QUAD_ARRAY_LOAD(rowcoefbinary, rowcoefs, vlbvarprobidx);
            varsolvalbinary = SCIPgetSolVal(scip, sol, vlbvars[bestlbtype[i]]);

            SCIPquadprecProdQD(val, rowcoef, vlbcoefs[bestlbtype[i]]);
            SCIPquadprecSumQQ(val, val, rowcoefbinary);
            {
               SCIP_Real QUAD(tmp);

               SCIPquadprecProdQD(tmp, rowcoefbinary, varsolvalbinary);
               SCIPquadprecSumDD(contsolval, solval, - vlbconsts[bestlbtype[i]]);
               SCIPquadprecProdQQ(contsolval, contsolval, rowcoef);
               SCIPquadprecSumQQ(contsolval, contsolval, tmp);
            }

            SCIPquadprecProdQD(rowcoeftimesvlbconst, rowcoef, vlbconsts[bestlbtype[i]]);

            /* clear the binvarpos array, since the variable has been processed */
            binvarused[vlbvarprobidx] = 0;

            /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
            snf->origbinvars[snf->ntransvars] = vlbvarprobidx;

            if( QUAD_TO_DBL(rowcoef) >= 0.0 )
            {
               snf->transvarcoefs[snf->ntransvars] = - 1;
               snf->transvarvubcoefs[snf->ntransvars] = - QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = varsolvalbinary;
               snf->transcontvarsolvals[snf->ntransvars] = - QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefsbin[snf->ntransvars] = - QUAD_TO_DBL(rowcoefbinary);
               snf->aggrcoefscont[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = QUAD_TO_DBL(rowcoeftimesvlbconst);
            }
            else
            {
               snf->transvarcoefs[snf->ntransvars] = 1;
               snf->transvarvubcoefs[snf->ntransvars] = QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = varsolvalbinary;
               snf->transcontvarsolvals[snf->ntransvars] = QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefsbin[snf->ntransvars] = QUAD_TO_DBL(rowcoefbinary);
               snf->aggrcoefscont[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = - QUAD_TO_DBL(rowcoeftimesvlbconst);
            }
            SCIPquadprecSumQQ(transrhs, transrhs, -rowcoeftimesvlbconst);

            SCIPdebugMsg(scip, "    --> bestlb used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s), rhs=%g-(%g*%g)=%g\n",
               snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars, snf->transvarvubcoefs[snf->ntransvars],
               snf->ntransvars, SCIPvarGetName(vlbvars[bestlbtype[i]]), QUAD_TO_DBL(transrhs) + QUAD_TO_DBL(rowcoeftimesvlbconst), QUAD_TO_DBL(rowcoef),
               vlbconsts[bestlbtype[i]], QUAD_TO_DBL(transrhs) );
         }
      }
      else
      {
         /* use bestub to define y'_j */

         assert(!SCIPisInfinity(scip, bestub[i]));
         assert(!SCIPisInfinity(scip, - bestslb[i]));
         assert(bestslbtype[i] == -1 || bestslbtype[i] == -2);
         assert(bestubtype[i] > -3 && bestubtype[i] < SCIPvarGetNVubs(var));

         /* store for y_j that y'_j is the associated real variable
          * in the relaxed set
          */
         snf->origcontvars[snf->ntransvars] = probidx;

         if( bestubtype[i] < 0 )
         {
            SCIP_Real QUAD(val);
            SCIP_Real QUAD(contsolval);
            SCIP_Real QUAD(rowcoeftimesbestslb);

            /* use simple upper bound in bestslb = l_j <= y_j <= u_j = bestub to define
             *   y'_j =   a_j ( y_j - l_j ) with 0 <= y'_j <=   a_j ( u_j - l_j ) x_j and x_j = 1    if a_j > 0
             *   y'_j = - a_j ( y_j - l_j ) with 0 <= y'_j <= - a_j ( u_j - l_j ) x_j and x_j = 1    if a_j < 0,
             * put j into the set
             *   N1   if a_j > 0
             *   N2   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j l_j
             */
            SCIPquadprecSumDD(val, bestub[i], - bestslb[i]);
            SCIPquadprecProdQQ(val, val, rowcoef);
            SCIPquadprecSumDD(contsolval, solval, - bestslb[i]);
            SCIPquadprecProdQQ(contsolval, contsolval, rowcoef);

            if( bestubtype[i] == -2 || bestslbtype[i] == -2 )
               *localbdsused = TRUE;

            SCIPquadprecProdQD(rowcoeftimesbestslb, rowcoef, bestslb[i]);

            /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
            snf->origbinvars[snf->ntransvars] = -1;
            snf->aggrcoefsbin[snf->ntransvars] = 0.0;

            if( QUAD_TO_DBL(rowcoef) >= 0.0 )
            {
               snf->transvarcoefs[snf->ntransvars] = 1;
               snf->transvarvubcoefs[snf->ntransvars] = QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = 1.0;
               snf->transcontvarsolvals[snf->ntransvars] = QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefscont[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = - QUAD_TO_DBL(rowcoeftimesbestslb);
            }
            else
            {
               snf->transvarcoefs[snf->ntransvars] = - 1;
               snf->transvarvubcoefs[snf->ntransvars] = - QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = 1.0;
               snf->transcontvarsolvals[snf->ntransvars] = - QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefscont[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = QUAD_TO_DBL(rowcoeftimesbestslb);
            }
            SCIPquadprecSumQQ(transrhs, transrhs, -rowcoeftimesbestslb);

            SCIPdebugMsg(scip, "    --> bestub used for trans: ... %s y'_%d + ..., Y'_%d <= %g x_%d (=1), rhs=%g-(%g*%g)=%g\n",
               snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars, snf->transvarvubcoefs[snf->ntransvars],
               snf->ntransvars, QUAD_TO_DBL(transrhs) + QUAD_TO_DBL(rowcoeftimesbestslb), QUAD_TO_DBL(rowcoef), bestslb[i], QUAD_TO_DBL(transrhs));
         }
         else
         {
            SCIP_Real QUAD(rowcoefbinary);
            SCIP_Real varsolvalbinary;
            SCIP_Real QUAD(val);
            SCIP_Real QUAD(contsolval);
            SCIP_Real QUAD(rowcoeftimesvubconst);
            int vubvarprobidx;

            SCIP_VAR** vubvars = SCIPvarGetVubVars(var);
            SCIP_Real* vubconsts = SCIPvarGetVubConstants(var);
            SCIP_Real* vubcoefs = SCIPvarGetVubCoefs(var);

            /* use variable upper bound in bestslb = l_j <= y_j <= u~_j x_j + d_j = bestub to define
             *   y'_j =     a_j ( y_j - d_j ) + c_j x_j   with 0 <= y'_j <=   ( a_j u~_j + c_j ) x_j    if a_j > 0
             *   y'_j = - ( a_j ( y_j - d_j ) + c_j x_j ) with 0 <= y'_j <= - ( a_j u~_j + c_j ) x_j    if a_j < 0,
             * where c_j is the coefficient of x_j in the row, put j into the set
             *   N1   if a_j > 0
             *   N2   if a_j < 0
             * and update the right hand side of the constraint in the relaxation
             *   rhs = rhs - a_j d_j
             */

            vubvarprobidx = SCIPvarGetProbindex(vubvars[bestubtype[i]]);
            assert(binvarused[vubvarprobidx] == 1);
            assert(vubvarprobidx < nbinvars);

            QUAD_ARRAY_LOAD(rowcoefbinary, rowcoefs, vubvarprobidx);
            varsolvalbinary = SCIPgetSolVal(scip, sol, vubvars[bestubtype[i]]);

            /* clear the binvarpos array, since the variable has been processed */
            binvarused[vubvarprobidx] = 0;

            SCIPquadprecProdQD(val, rowcoef, vubcoefs[bestubtype[i]]);
            SCIPquadprecSumQQ(val, val, rowcoefbinary);
            {
               SCIP_Real QUAD(tmp);
               SCIPquadprecProdQD(tmp, rowcoefbinary, varsolvalbinary);
               SCIPquadprecSumDD(contsolval, solval, - vubconsts[bestubtype[i]]);
               SCIPquadprecProdQQ(contsolval, contsolval, rowcoef);
               SCIPquadprecSumQQ(contsolval, contsolval, tmp);
            }

            SCIPquadprecProdQD(rowcoeftimesvubconst, rowcoef, vubconsts[bestubtype[i]]);
            /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
            snf->origbinvars[snf->ntransvars] = vubvarprobidx;

            if( QUAD_TO_DBL(rowcoef) >= 0.0 )
            {
               snf->transvarcoefs[snf->ntransvars] = 1;
               snf->transvarvubcoefs[snf->ntransvars] = QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = varsolvalbinary;
               snf->transcontvarsolvals[snf->ntransvars] = QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefsbin[snf->ntransvars] = QUAD_TO_DBL(rowcoefbinary);
               snf->aggrcoefscont[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = - QUAD_TO_DBL(rowcoeftimesvubconst);
            }
            else
            {
               snf->transvarcoefs[snf->ntransvars] = - 1;
               snf->transvarvubcoefs[snf->ntransvars] = - QUAD_TO_DBL(val);
               snf->transbinvarsolvals[snf->ntransvars] = varsolvalbinary;
               snf->transcontvarsolvals[snf->ntransvars] = - QUAD_TO_DBL(contsolval);

               /* aggregation information for y'_j */
               snf->aggrcoefsbin[snf->ntransvars] = - QUAD_TO_DBL(rowcoefbinary);
               snf->aggrcoefscont[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
               snf->aggrconstants[snf->ntransvars] = QUAD_TO_DBL(rowcoeftimesvubconst);
            }
            SCIPquadprecSumQQ(transrhs, transrhs, -rowcoeftimesvubconst);

            /* store for x_j that y'_j is the associated real variable in the 0-1 single node flow relaxation */

            SCIPdebugMsg(scip, "    --> bestub used for trans: ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s), rhs=%g-(%g*%g)=%g\n",
               snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars, snf->transvarvubcoefs[snf->ntransvars],
               snf->ntransvars, SCIPvarGetName(vubvars[bestubtype[i]]), QUAD_TO_DBL(transrhs) + QUAD_TO_DBL(rowcoeftimesvubconst), QUAD_TO_DBL(rowcoef),
               vubconsts[bestubtype[i]], QUAD_TO_DBL(transrhs));
         }
      }

      /* make sure the coefficient is not negative due to small numerical rounding errors */
      assert(snf->transvarvubcoefs[snf->ntransvars] > -QUAD_EPSILON);
      snf->transvarvubcoefs[snf->ntransvars] = MAX(snf->transvarvubcoefs[snf->ntransvars], 0.0);

      ++snf->ntransvars;
   }

   snf->transrhs = QUAD_TO_DBL(transrhs);

   /* transform remaining binary variables of row */
   for( i = nnonbinvarsrow; i < nnz; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real QUAD(rowcoef);
      int probidx;
      SCIP_Real val;
      SCIP_Real contsolval;
      SCIP_Real varsolval;

      probidx = rowinds[i];
      /* variable should be binary */
      assert(probidx >= 0);
      assert(probidx < nbinvars);

      /* binary variable was processed together with a non-binary variable */
      if( binvarused[probidx] == 0 )
         continue;

      /* binary variable was not processed yet, so the binvarused value sould be -1 */
      assert(binvarused[probidx] == -1);

      /* set binvarused to zero since it has been processed */
      binvarused[probidx] = 0;

      var = vars[probidx];
      QUAD_ARRAY_LOAD(rowcoef, rowcoefs, probidx);
      assert(!EPSZ(QUAD_TO_DBL(rowcoef), QUAD_EPSILON));

      varsolval = SCIPgetSolVal(scip, sol, var);
      SCIPdebugMsg(scip, "  %d: %g <%s, idx=%d, lp=%g, [%g, %g]>:\n", i, QUAD_TO_DBL(rowcoef), SCIPvarGetName(var), probidx, varsolval,
         SCIPvarGetLbGlobal(var), SCIPvarGetUbGlobal(var));

      /* define
       *    y'_j =   c_j x_j with 0 <= y'_j <=   c_j x_j    if c_j > 0
       *    y'_j = - c_j x_j with 0 <= y'_j <= - c_j x_j    if c_j < 0,
       * where c_j is the coefficient of x_j in the row and put j into the set
       *    N1   if c_j > 0
       *    N2   if c_j < 0.
       */
      val = QUAD_TO_DBL(rowcoef);
      contsolval = QUAD_TO_DBL(rowcoef) * varsolval;

      /* store aggregation information for y'_j for transforming cuts for the SNF relaxation back to the problem variables later */
      snf->origbinvars[snf->ntransvars] = probidx;
      snf->origcontvars[snf->ntransvars] = -1;
      snf->aggrcoefscont[snf->ntransvars] = 0.0;
      snf->aggrconstants[snf->ntransvars] = 0.0;

      if( QUAD_TO_DBL(rowcoef) >= 0.0 )
      {
         snf->transvarcoefs[snf->ntransvars] = 1;
         snf->transvarvubcoefs[snf->ntransvars] = val;
         snf->transbinvarsolvals[snf->ntransvars] = varsolval;
         snf->transcontvarsolvals[snf->ntransvars] = contsolval;

         /* aggregation information for y'_j */
         snf->aggrcoefsbin[snf->ntransvars] = QUAD_TO_DBL(rowcoef);
      }
      else
      {
         snf->transvarcoefs[snf->ntransvars] = - 1;
         snf->transvarvubcoefs[snf->ntransvars] = - val;
         snf->transbinvarsolvals[snf->ntransvars] = varsolval;
         snf->transcontvarsolvals[snf->ntransvars] = - contsolval;

         /* aggregation information for y'_j */
         snf->aggrcoefsbin[snf->ntransvars] = - QUAD_TO_DBL(rowcoef);
      }

      assert(snf->transvarcoefs[snf->ntransvars] == 1 || snf->transvarcoefs[snf->ntransvars] == - 1 );
      assert(SCIPisFeasGE(scip, snf->transbinvarsolvals[snf->ntransvars], 0.0)
         && SCIPisFeasLE(scip, snf->transbinvarsolvals[snf->ntransvars], 1.0));
      assert(SCIPisFeasGE(scip, snf->transvarvubcoefs[snf->ntransvars], 0.0)
         && !SCIPisInfinity(scip, snf->transvarvubcoefs[snf->ntransvars]));

      SCIPdebugMsg(scip, "   --> ... %s y'_%d + ..., y'_%d <= %g x_%d (=%s))\n", snf->transvarcoefs[snf->ntransvars] == 1 ? "+" : "-", snf->ntransvars, snf->ntransvars,
         snf->transvarvubcoefs[snf->ntransvars], snf->ntransvars, SCIPvarGetName(var) );

      /* updates number of variables in transformed problem */
      snf->ntransvars++;
   }

   /* construction was successful */
   *success = TRUE;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "constraint in constructed 0-1 single node flow relaxation: ");
   for( i = 0; i < snf->ntransvars; i++ )
   {
      SCIPdebugMsgPrint(scip, "%s y'_%d ", snf->transvarcoefs[i] == 1 ? "+" : "-", i);
   }
   SCIPdebugMsgPrint(scip, "<= %g\n", snf->transrhs);
#endif

  TERMINATE:

   SCIPfreeCleanBufferArray(scip, &binvarused);
   SCIPfreeBufferArray(scip, &selectedbounds);
   SCIPfreeBufferArray(scip, &bestsubtype);
   SCIPfreeBufferArray(scip, &bestslbtype);
   SCIPfreeBufferArray(scip, &bestubtype);
   SCIPfreeBufferArray(scip, &bestlbtype);
   SCIPfreeBufferArray(scip, &bestsub);
   SCIPfreeBufferArray(scip, &bestslb);
   SCIPfreeBufferArray(scip, &bestub);
   SCIPfreeBufferArray(scip, &bestlb);

   return SCIP_OKAY;
}

/** allocate buffer arrays for storing the single-node-flow relaxation */
static
SCIP_RETCODE allocSNFRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< pointer to snf relaxation to be destroyed */
   int                   nvars               /**< number of active problem variables */
   )
{
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->transvarcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->transbinvarsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->transcontvarsolvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->transvarvubcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->origbinvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->origcontvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->aggrcoefsbin, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->aggrcoefscont, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &snf->aggrconstants, nvars) );

   return SCIP_OKAY;
}

/** free buffer arrays for storing the single-node-flow relaxation */
static
void destroySNFRelaxation(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf                 /**< pointer to snf relaxation to be destroyed */
   )
{
   SCIPfreeBufferArray(scip, &snf->aggrconstants);
   SCIPfreeBufferArray(scip, &snf->aggrcoefscont);
   SCIPfreeBufferArray(scip, &snf->aggrcoefsbin);
   SCIPfreeBufferArray(scip, &snf->origcontvars);
   SCIPfreeBufferArray(scip, &snf->origbinvars);
   SCIPfreeBufferArray(scip, &snf->transvarvubcoefs);
   SCIPfreeBufferArray(scip, &snf->transcontvarsolvals);
   SCIPfreeBufferArray(scip, &snf->transbinvarsolvals);
   SCIPfreeBufferArray(scip, &snf->transvarcoefs);
}

/** solve knapsack problem in maximization form with "<" constraint approximately by greedy; if needed, one can provide
 *  arrays to store all selected items and all not selected items
 */
static
SCIP_RETCODE SCIPsolveKnapsackApproximatelyLT(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nitems,             /**< number of available items */
   SCIP_Real*            weights,            /**< item weights */
   SCIP_Real*            profits,            /**< item profits */
   SCIP_Real             capacity,           /**< capacity of knapsack */
   int*                  items,              /**< item numbers */
   int*                  solitems,           /**< array to store items in solution, or NULL */
   int*                  nonsolitems,        /**< array to store items not in solution, or NULL */
   int*                  nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*                  nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   SCIP_Real*            solval              /**< pointer to store optimal solution value, or NULL */
   )
{
   SCIP_Real* tempsort;
   SCIP_Real solitemsweight;
   SCIP_Real mediancapacity;
   int j;
   int i;
   int criticalitem;

   assert(weights != NULL);
   assert(profits != NULL);
   assert(SCIPisFeasGE(scip, capacity, 0.0));
   assert(!SCIPisInfinity(scip, capacity));
   assert(items != NULL);
   assert(nitems >= 0);

   if( solitems != NULL )
   {
      *nsolitems = 0;
      *nnonsolitems = 0;
   }
   if( solval != NULL )
      *solval = 0.0;

   /* allocate memory for temporary array used for sorting; array should contain profits divided by corresponding weights (p_1 / w_1 ... p_n / w_n )*/
   SCIP_CALL( SCIPallocBufferArray(scip, &tempsort, nitems) );

   /* initialize temporary array */
   for( i = nitems - 1; i >= 0; --i )
      tempsort[i] = profits[i] / weights[i];

   /* decrease capacity slightly to make it tighter than the original capacity */
   mediancapacity = capacity * (1 - SCIPfeastol(scip));

   /* rearrange items around  */
   SCIPselectWeightedDownRealRealInt(tempsort, profits, items, weights, mediancapacity, nitems, &criticalitem);

   /* free temporary array */
   SCIPfreeBufferArray(scip, &tempsort);

   /* select items as long as they fit into the knapsack */
   solitemsweight = 0.0;
   for( j = 0; j < nitems && SCIPisFeasLT(scip, solitemsweight + weights[j], capacity); j++ )
   {
      if( solitems != NULL )
      {
         solitems[*nsolitems] = items[j];
         (*nsolitems)++;
      }
      if( solval != NULL )
         (*solval) += profits[j];
      solitemsweight += weights[j];
   }

   /* continue to put items into the knapsack if they entirely fit */
   for( ; j < nitems; j++ )
   {
      if( SCIPisFeasLT(scip, solitemsweight + weights[j], capacity) )
      {
         if( solitems != NULL )
         {
            solitems[*nsolitems] = items[j];
            (*nsolitems)++;
         }
         if( solval != NULL )
            (*solval) += profits[j];
         solitemsweight += weights[j];
      }
      else if( solitems != NULL )
      {
         nonsolitems[*nnonsolitems] = items[j];
         (*nnonsolitems)++;
      }
   }

   return SCIP_OKAY;
}


/** build the flow cover which corresponds to the given exact or approximate solution of KP^SNF; given unfinished
 *  flow cover contains variables which have been fixed in advance
 */
static
void buildFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  coefs,              /**< coefficient of all real variables in N1&N2 */
   SCIP_Real*            vubcoefs,           /**< coefficient in vub of all real variables in N1&N2 */
   SCIP_Real             rhs,                /**< right hand side of 0-1 single node flow constraint */
   int*                  solitems,           /**< items in knapsack */
   int*                  nonsolitems,        /**< items not in knapsack */
   int                   nsolitems,          /**< number of items in knapsack */
   int                   nnonsolitems,       /**< number of items not in knapsack */
   int*                  nflowcovervars,     /**< pointer to store number of variables in flow cover */
   int*                  nnonflowcovervars,  /**< pointer to store number of variables not in flow cover */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */
   QUAD(SCIP_Real*       flowcoverweight),   /**< pointer to store weight of flow cover */
   SCIP_Real*            lambda              /**< pointer to store lambda */
   )
{
   int j;
   SCIP_Real QUAD(tmp);

   assert(scip != NULL);
   assert(coefs != NULL);
   assert(vubcoefs != NULL);
   assert(solitems != NULL);
   assert(nonsolitems != NULL);
   assert(nsolitems >= 0);
   assert(nnonsolitems >= 0);
   assert(nflowcovervars != NULL && *nflowcovervars >= 0);
   assert(nnonflowcovervars != NULL && *nnonflowcovervars >= 0);
   assert(flowcoverstatus != NULL);
   assert(QUAD_HI(flowcoverweight) != NULL);
   assert(lambda != NULL);

   /* get flowcover status for each item */
   for( j = 0; j < nsolitems; j++ )
   {
      /* j in N1 with z°_j = 1 => j in N1\C1 */
      if( coefs[solitems[j]] == 1 )
      {
         flowcoverstatus[solitems[j]] = -1;
         (*nnonflowcovervars)++;
      }
      /* j in N2 with z_j = 1 => j in C2 */
      else
      {
         assert(coefs[solitems[j]] == -1);
         flowcoverstatus[solitems[j]] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(*flowcoverweight, *flowcoverweight, -vubcoefs[solitems[j]]);
      }
   }
   for( j = 0; j < nnonsolitems; j++ )
   {
      /* j in N1 with z°_j = 0 => j in C1 */
      if( coefs[nonsolitems[j]] == 1 )
      {
         flowcoverstatus[nonsolitems[j]] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(*flowcoverweight, *flowcoverweight, vubcoefs[nonsolitems[j]]);
      }
      /* j in N2 with z_j = 0 => j in N2\C2 */
      else
      {
         assert(coefs[nonsolitems[j]] == -1);
         flowcoverstatus[nonsolitems[j]] = -1;
         (*nnonflowcovervars)++;
      }
   }

   /* get lambda = sum_{j in C1} u_j - sum_{j in C2} u_j - rhs */
   SCIPquadprecSumQD(tmp, *flowcoverweight, -rhs);
   *lambda = QUAD_TO_DBL(tmp);
}

#ifndef NO_EXACT_KNAPSACK

/** checks whether the given scalar scales the given value to an integral number with error in the given bounds */
static
SCIP_Bool isIntegralScalar(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real downval;
   SCIP_Real upval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   downval = floor(sval);
   upval = ceil(sval);

   return (SCIPrelDiff(sval, downval) <= maxdelta || SCIPrelDiff(sval, upval) >= mindelta);
}

/** get integral number with error in the bounds which corresponds to given value scaled by a given scalar;
 *  should be used in connection with isIntegralScalar()
 */
static
SCIP_Longint getIntegralVal(
   SCIP_Real             val,                /**< value that should be scaled to an integral value */
   SCIP_Real             scalar,             /**< scalar that should be tried */
   SCIP_Real             mindelta,           /**< minimal relative allowed difference of scaled coefficient s*c and integral i */
   SCIP_Real             maxdelta            /**< maximal relative allowed difference of scaled coefficient s*c and integral i */
   )
{
   SCIP_Real sval;
   SCIP_Real upval;
   SCIP_Longint intval;

   assert(mindelta <= 0.0);
   assert(maxdelta >= 0.0);

   sval = val * scalar;
   upval = ceil(sval);

   if( SCIPrelDiff(sval, upval) >= mindelta )
      intval = (SCIP_Longint) upval;
   else
      intval = (SCIP_Longint) (floor(sval));

   return intval;
}

/** get a flow cover (C1, C2) for a given 0-1 single node flow set
 *    {(x,y) in {0,1}^n x R^n : sum_{j in N1} y_j - sum_{j in N2} y_j <= b, 0 <= y_j <= u_j x_j},
 *  i.e., get sets C1 subset N1 and C2 subset N2 with sum_{j in C1} u_j - sum_{j in C2} u_j = b + lambda and lambda > 0
 */
static
SCIP_RETCODE getFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< the single node flow relaxation */
   int*                  nflowcovervars,     /**< pointer to store number of variables in flow cover */
   int*                  nnonflowcovervars,  /**< pointer to store number of variables not in flow cover */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */
   SCIP_Real*            lambda,             /**< pointer to store lambda */
   SCIP_Bool*            found               /**< pointer to store whether a cover was found */
   )
{
   SCIP_Real* transprofitsint;
   SCIP_Real* transprofitsreal;
   SCIP_Real* transweightsreal;
   SCIP_Longint* transweightsint;
   int* items;
   int* itemsint;
   int* nonsolitems;
   int* solitems;
   SCIP_Real QUAD(flowcoverweight);
   SCIP_Real QUAD(flowcoverweightafterfix);
   SCIP_Real n1itemsweight;
   SCIP_Real n2itemsminweight;
   SCIP_Real scalar;
   SCIP_Real transcapacityreal;
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
   SCIP_Bool kpexact;
#endif
   SCIP_Bool scalesuccess;
   SCIP_Bool transweightsrealintegral;
   SCIP_Longint transcapacityint;
   int nflowcovervarsafterfix;
   int nitems;
   int nn1items;
   int nnonflowcovervarsafterfix;
   int nnonsolitems;
   int nsolitems;
   int j;

   assert(scip != NULL);
   assert(snf->transvarcoefs != NULL);
   assert(snf->transbinvarsolvals != NULL);
   assert(snf->transvarvubcoefs != NULL);
   assert(snf->ntransvars > 0);
   assert(nflowcovervars != NULL);
   assert(nnonflowcovervars != NULL);
   assert(flowcoverstatus != NULL);
   assert(lambda != NULL);
   assert(found != NULL);

   SCIPdebugMsg(scip, "--------------------- get flow cover ----------------------------------------------------\n");

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &items, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &itemsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofitsreal, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofitsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsreal, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, snf->ntransvars) );

   BMSclearMemoryArray(flowcoverstatus, snf->ntransvars);
   *found = FALSE;
   *nflowcovervars = 0;
   *nnonflowcovervars = 0;

   QUAD_ASSIGN(flowcoverweight, 0.0);
   nflowcovervarsafterfix = 0;
   nnonflowcovervarsafterfix = 0;
   QUAD_ASSIGN(flowcoverweightafterfix, 0.0);
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
   kpexact = FALSE;
#endif

   /* fix some variables in advance according to the following fixing strategy
    *   put j into N1\C1,          if j in N1 and x*_j = 0,
    *   put j into C1,             if j in N1 and x*_j = 1,
    *   put j into C2,             if j in N2 and x*_j = 1,
    *   put j into N2\C2,          if j in N2 and x*_j = 0
    * and get the set of the remaining variables
    */
   SCIPdebugMsg(scip, "0. Fix some variables in advance:\n");
   nitems = 0;
   nn1items = 0;
   n1itemsweight = 0.0;
   n2itemsminweight = SCIP_REAL_MAX;
   for( j = 0; j < snf->ntransvars; j++ )
   {
      assert(snf->transvarcoefs[j] == 1 || snf->transvarcoefs[j] == -1);
      assert(SCIPisFeasGE(scip, snf->transbinvarsolvals[j], 0.0) && SCIPisFeasLE(scip,  snf->transbinvarsolvals[j], 1.0));
      assert(SCIPisFeasGE(scip, snf->transvarvubcoefs[j], 0.0));

      /* if u_j = 0, put j into N1\C1 and N2\C2, respectively */
      if( SCIPisFeasZero(scip, snf->transvarvubcoefs[j]) )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         continue;
      }

      /* x*_j is fractional */
      if( !SCIPisFeasIntegral(scip,  snf->transbinvarsolvals[j]) )
      {
         items[nitems] = j;
         nitems++;
         if( snf->transvarcoefs[j] == 1 )
         {
            n1itemsweight += snf->transvarvubcoefs[j];
            nn1items++;
         }
         else
            n2itemsminweight = MIN(n2itemsminweight, snf->transvarvubcoefs[j]);
      }
      /* j is in N1 and x*_j = 0 */
      else if( snf->transvarcoefs[j] == 1 &&  snf->transbinvarsolvals[j] < 0.5 )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMsg(scip, "     <%d>: in N1-C1\n", j);
      }
      /* j is in N1 and x*_j = 1 */
      else if( snf->transvarcoefs[j] == 1 &&  snf->transbinvarsolvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(flowcoverweight, flowcoverweight, snf->transvarvubcoefs[j]);
         SCIPdebugMsg(scip, "     <%d>: in C1\n", j);
      }
      /* j is in N2 and x*_j = 1 */
      else if( snf->transvarcoefs[j] == -1 &&  snf->transbinvarsolvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(flowcoverweight, flowcoverweight, -snf->transvarvubcoefs[j]);
         SCIPdebugMsg(scip, "     <%d>: in C2\n", j);
      }
      /* j is in N2 and x*_j = 0 */
      else
      {
         assert(snf->transvarcoefs[j] == -1 &&  snf->transbinvarsolvals[j] < 0.5);
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMsg(scip, "     <%d>: in N2-C2\n", j);
      }
   }
   assert((*nflowcovervars) + (*nnonflowcovervars) + nitems == snf->ntransvars);
   assert(nn1items >= 0);

   /* to find a flow cover, transform the following knapsack problem
    *
    * (KP^SNF)      max sum_{j in N1} ( x*_j - 1 ) z_j + sum_{j in N2} x*_j z_j
    *                   sum_{j in N1}          u_j z_j - sum_{j in N2} u_j  z_j > b
    *                                         z_j in {0,1} for all j in N1 & N2
    *
    * 1. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have
    *    positive weights and the constraint is a "<" constraint, by complementing all variables in N1
    *
    *    (KP^SNF_rat)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}          u_j z°_j + sum_{j in N2} u_j  z_j < - b + sum_{j in N1} u_j
    *                                                 z°_j in {0,1} for all j in N1
    *                                                  z_j in {0,1} for all j in N2,
    *    and solve it approximately under consideration of the fixing,
    * or
    * 2. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have
    *    positive integer weights and the constraint is a "<=" constraint, by complementing all variables in N1
    *    and multiplying the constraint by a suitable scalar C
    *
    *    (KP^SNF_int)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}        C u_j z°_j + sum_{j in N2} C u_j  z_j <= c
    *                                                   z°_j in {0,1} for all j in N1
    *                                                    z_j in {0,1} for all j in N2,
    *    where
    *      c = floor[ C (- b + sum_{j in N1} u_j ) ]      if frac[ C (- b + sum_{j in N1} u_j ) ] > 0
    *      c =        C (- b + sum_{j in N1} u_j )   - 1  if frac[ C (- b + sum_{j in N1} u_j ) ] = 0
    *    and solve it exactly under consideration of the fixing.
    */
   SCIPdebugMsg(scip, "1. Transform KP^SNF to KP^SNF_rat:\n");

   /* get weight and profit of variables in KP^SNF_rat and check, whether all weights are already integral */
   transweightsrealintegral = TRUE;
   for( j = 0; j < nitems; j++ )
   {
      transweightsreal[j] = snf->transvarvubcoefs[items[j]];

      if( !isIntegralScalar(transweightsreal[j], 1.0, -MINDELTA, MAXDELTA) )
         transweightsrealintegral = FALSE;

      if( snf->transvarcoefs[items[j]] == 1 )
      {
         transprofitsreal[j] = 1.0 -  snf->transbinvarsolvals[items[j]];
         SCIPdebugMsg(scip, "     <%d>: j in N1:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j],
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
      else
      {
         transprofitsreal[j] =  snf->transbinvarsolvals[items[j]];
         SCIPdebugMsg(scip, "     <%d>: j in N2:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j],
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
   }
   /* get capacity of knapsack constraint in KP^SNF_rat */
   transcapacityreal = - snf->transrhs + QUAD_TO_DBL(flowcoverweight) + n1itemsweight;
   SCIPdebugMsg(scip, "     transcapacity = -rhs(%g) + flowcoverweight(%g) + n1itemsweight(%g) = %g\n",
      snf->transrhs, QUAD_TO_DBL(flowcoverweight), n1itemsweight, transcapacityreal);

   /* there exists no flow cover if the capacity of knapsack constraint in KP^SNF_rat after fixing
    * is less than or equal to zero
    */
   if( SCIPisFeasLE(scip, transcapacityreal/10, 0.0) )
   {
      assert(!(*found));
      goto TERMINATE;
   }

   /* KP^SNF_rat has been solved by fixing some variables in advance */
   assert(nitems >= 0);
   if( nitems == 0)
   {
      /* get lambda = sum_{j in C1} u_j - sum_{j in C2} u_j - rhs */
      SCIPquadprecSumQD(flowcoverweight, flowcoverweight, -snf->transrhs);
      *lambda = QUAD_TO_DBL(flowcoverweight);
      *found = TRUE;
      goto TERMINATE;
   }

   /* Use the following strategy
    *   solve KP^SNF_int exactly,          if a suitable factor C is found and (nitems*capacity) <= MAXDYNPROGSPACE,
    *   solve KP^SNF_rat approximately,    otherwise
    */

   /* find a scaling factor C */
   if( transweightsrealintegral )
   {
      /* weights are already integral */
      scalar = 1.0;
      scalesuccess = TRUE;
   }
   else
   {
      scalesuccess = FALSE;
      SCIP_CALL( SCIPcalcIntegralScalar(transweightsreal, nitems, -MINDELTA, MAXDELTA, MAXDNOM, MAXSCALE, &scalar,
            &scalesuccess) );
   }

   /* initialize number of (non-)solution items, should be changed to a nonnegative number in all possible paths below */
   nsolitems = -1;
   nnonsolitems = -1;

   /* suitable factor C was found*/
   if( scalesuccess )
   {
      SCIP_Real tmp1;
      SCIP_Real tmp2;

      /* transform KP^SNF to KP^SNF_int */
      for( j = 0; j < nitems; ++j )
      {
         transweightsint[j] = getIntegralVal(transweightsreal[j], scalar, -MINDELTA, MAXDELTA);
         transprofitsint[j] = transprofitsreal[j];
         itemsint[j] = items[j];
      }
      if( isIntegralScalar(transcapacityreal, scalar, -MINDELTA, MAXDELTA) )
      {
         transcapacityint = getIntegralVal(transcapacityreal, scalar, -MINDELTA, MAXDELTA);
         transcapacityint -= 1;
      }
      else
         transcapacityint = (SCIP_Longint) (transcapacityreal * scalar);
      nflowcovervarsafterfix = *nflowcovervars;
      nnonflowcovervarsafterfix = *nnonflowcovervars;
      QUAD_ASSIGN_Q(flowcoverweightafterfix, flowcoverweight);

      tmp1 = (SCIP_Real) (nitems + 1);
      tmp2 = (SCIP_Real) ((transcapacityint) + 1);
      if( transcapacityint * nitems <= MAXDYNPROGSPACE && tmp1 * tmp2 <= INT_MAX / 8.0)
      {
         SCIP_Bool success;

         /* solve KP^SNF_int by dynamic programming */
         SCIP_CALL( SCIPsolveKnapsackExactly(scip, nitems, transweightsint, transprofitsint, transcapacityint,
               itemsint, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL, &success) );

         if( !success )
         {
            /* solve KP^SNF_rat approximately */
            SCIP_CALL( SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal,
                  transcapacityreal, items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL) );
         }
#if !defined(NDEBUG) || defined(SCIP_DEBUG)
         else
            kpexact = TRUE;
#endif
      }
      else
      {
         /* solve KP^SNF_rat approximately */
         SCIP_CALL( SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
               items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL) );
         assert(!kpexact);
      }
   }
   else
   {
      /* solve KP^SNF_rat approximately */
      SCIP_CALL( SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
            items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL) );
      assert(!kpexact);
   }

   assert(nsolitems != -1);
   assert(nnonsolitems != -1);

   /* build the flow cover from the solution of KP^SNF_rat and KP^SNF_int, respectively and the fixing */
   assert(*nflowcovervars + *nnonflowcovervars + nsolitems + nnonsolitems == snf->ntransvars);
   buildFlowCover(scip, snf->transvarcoefs, snf->transvarvubcoefs, snf->transrhs, solitems, nonsolitems, nsolitems, nnonsolitems, nflowcovervars,
      nnonflowcovervars, flowcoverstatus, QUAD(&flowcoverweight), lambda);
   assert(*nflowcovervars + *nnonflowcovervars == snf->ntransvars);

   /* if the found structure is not a flow cover, because of scaling, solve KP^SNF_rat approximately */
   if( SCIPisFeasLE(scip, *lambda, 0.0) )
   {
      assert(kpexact);

      /* solve KP^SNF_rat approximately */
      SCIP_CALL( SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
            items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL) );
#ifdef SCIP_DEBUG /* this time only for SCIP_DEBUG, because only then, the variable is used again  */
      kpexact = FALSE;
#endif

      /* build the flow cover from the solution of KP^SNF_rat and the fixing */
      *nflowcovervars = nflowcovervarsafterfix;
      *nnonflowcovervars = nnonflowcovervarsafterfix;
      QUAD_ASSIGN_Q(flowcoverweight, flowcoverweightafterfix);

      assert(*nflowcovervars + *nnonflowcovervars + nsolitems + nnonsolitems == snf->ntransvars);
      buildFlowCover(scip, snf->transvarcoefs, snf->transvarvubcoefs, snf->transrhs, solitems, nonsolitems, nsolitems, nnonsolitems, nflowcovervars,
         nnonflowcovervars, flowcoverstatus, QUAD(&flowcoverweight), lambda);
      assert(*nflowcovervars + *nnonflowcovervars == snf->ntransvars);
   }
   *found = SCIPisFeasGT(scip, *lambda, 0.0);

  TERMINATE:
   assert((!*found) || SCIPisFeasGT(scip, *lambda, 0.0));
#ifdef SCIP_DEBUG
   if( *found )
   {
      SCIPdebugMsg(scip, "2. %s solution:\n", kpexact ? "exact" : "approximate");
      for( j = 0; j < snf->ntransvars; j++ )
      {
         if( snf->transvarcoefs[j] == 1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMsg(scip, "     C1: + y_%d [u_%d = %g]\n", j, j, snf->transvarvubcoefs[j]);
         }
         else if( snf->transvarcoefs[j] == -1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMsg(scip, "     C2: - y_%d [u_%d = %g]\n", j, j, snf->transvarvubcoefs[j]);
         }
      }
      SCIPdebugMsg(scip, "     flowcoverweight(%g) = rhs(%g) + lambda(%g)\n", QUAD_TO_DBL(flowcoverweight), snf->transrhs, *lambda);
   }
#endif

   /* free data structures */
   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &transweightsint);
   SCIPfreeBufferArray(scip, &transweightsreal);
   SCIPfreeBufferArray(scip, &transprofitsint);
   SCIPfreeBufferArray(scip, &transprofitsreal);
   SCIPfreeBufferArray(scip, &itemsint);
   SCIPfreeBufferArray(scip, &items);

   return SCIP_OKAY;
}

#else

/** get a flow cover \f$(C1, C2)\f$ for a given 0-1 single node flow set
 *    \f${(x,y) in {0,1}^n x R^n : sum_{j in N1} y_j - sum_{j in N2} y_j <= b, 0 <= y_j <= u_j x_j}\f$,
 *  i.e., get sets \f$ C1 \subset N1 \f$ and \f$ C2 \subset N2 \f$ with
 *  \f$ \sum_{j in C1} u_j - sum_{j in C2} u_j = b + lambda \f$ and \f$ lambda > 0 \f$
 */
static
SCIP_RETCODE getFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< the 0-1 single node flow relaxation */
   int*                  nflowcovervars,     /**< pointer to store number of variables in flow cover */
   int*                  nnonflowcovervars,  /**< pointer to store number of variables not in flow cover */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */
   SCIP_Real*            lambda,             /**< pointer to store lambda */
   SCIP_Bool*            found               /**< pointer to store whether a cover was found */
   )
{
   SCIP_Real* transprofitsreal;
   SCIP_Real* transweightsreal;
   SCIP_Longint* transweightsint;
   int* items;
   int* itemsint;
   int* nonsolitems;
   int* solitems;
   SCIP_Real QUAD(flowcoverweight);
   SCIP_Real n1itemsweight;
   SCIP_Real n2itemsminweight;
   SCIP_Real transcapacityreal;
   int nitems;
#ifndef NDEBUG
   int nn1items = 0;
#endif
   int nnonsolitems;
   int nsolitems;
   int j;

   assert(scip != NULL);
   assert(snf->transvarcoefs != NULL);
   assert(snf->transbinvarsolvals != NULL);
   assert(snf->transvarvubcoefs != NULL);
   assert(snf->ntransvars > 0);
   assert(nflowcovervars != NULL);
   assert(nnonflowcovervars != NULL);
   assert(flowcoverstatus != NULL);
   assert(lambda != NULL);
   assert(found != NULL);

   SCIPdebugMsg(scip, "--------------------- get flow cover ----------------------------------------------------\n");

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &items, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &itemsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transprofitsreal, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsreal, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transweightsint, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &solitems, snf->ntransvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nonsolitems, snf->ntransvars) );

   BMSclearMemoryArray(flowcoverstatus, snf->ntransvars);
   *found = FALSE;
   *nflowcovervars = 0;
   *nnonflowcovervars = 0;

   QUAD_ASSIGN(flowcoverweight, 0.0);

   /* fix some variables in advance according to the following fixing strategy
    *   put j into N1\C1,          if j in N1 and x*_j = 0,
    *   put j into C1,             if j in N1 and x*_j = 1,
    *   put j into C2,             if j in N2 and x*_j = 1,
    *   put j into N2\C2,          if j in N2 and x*_j = 0
    * and get the set of the remaining variables
    */
   SCIPdebugMsg(scip, "0. Fix some variables in advance:\n");
   nitems = 0;
   n1itemsweight = 0.0;
   n2itemsminweight = SCIP_REAL_MAX;
   for( j = 0; j < snf->ntransvars; j++ )
   {
      assert(snf->transvarcoefs[j] == 1 || snf->transvarcoefs[j] == -1);
      assert(SCIPisFeasGE(scip, snf->transbinvarsolvals[j], 0.0) && SCIPisFeasLE(scip,  snf->transbinvarsolvals[j], 1.0));
      assert(SCIPisFeasGE(scip, snf->transvarvubcoefs[j], 0.0));

      /* if u_j = 0, put j into N1\C1 and N2\C2, respectively */
      if( SCIPisFeasZero(scip, snf->transvarvubcoefs[j]) )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         continue;
      }

      /* x*_j is fractional */
      if( !SCIPisFeasIntegral(scip,  snf->transbinvarsolvals[j]) )
      {
         items[nitems] = j;
         nitems++;
         if( snf->transvarcoefs[j] == 1 )
         {
            n1itemsweight += snf->transvarvubcoefs[j];
#ifndef NDEBUG
            nn1items++;
#endif
         }
         else
            n2itemsminweight = MIN(n2itemsminweight, snf->transvarvubcoefs[j]);
      }
      /* j is in N1 and x*_j = 0 */
      else if( snf->transvarcoefs[j] == 1 &&  snf->transbinvarsolvals[j] < 0.5 )
      {
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMsg(scip, "     <%d>: in N1-C1\n", j);
      }
      /* j is in N1 and x*_j = 1 */
      else if( snf->transvarcoefs[j] == 1 &&  snf->transbinvarsolvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(flowcoverweight, flowcoverweight, snf->transvarvubcoefs[j]);
         SCIPdebugMsg(scip, "     <%d>: in C1\n", j);
      }
      /* j is in N2 and x*_j = 1 */
      else if( snf->transvarcoefs[j] == -1 &&  snf->transbinvarsolvals[j] > 0.5 )
      {
         flowcoverstatus[j] = 1;
         (*nflowcovervars)++;
         SCIPquadprecSumQD(flowcoverweight, flowcoverweight, -snf->transvarvubcoefs[j]);
         SCIPdebugMsg(scip, "     <%d>: in C2\n", j);
      }
      /* j is in N2 and x*_j = 0 */
      else
      {
         assert(snf->transvarcoefs[j] == -1 &&  snf->transbinvarsolvals[j] < 0.5);
         flowcoverstatus[j] = -1;
         (*nnonflowcovervars)++;
         SCIPdebugMsg(scip, "     <%d>: in N2-C2\n", j);
      }
   }
   assert((*nflowcovervars) + (*nnonflowcovervars) + nitems == snf->ntransvars);
   assert(nn1items >= 0);

   /* to find a flow cover, transform the following knapsack problem
    *
    * (KP^SNF)      max sum_{j in N1} ( x*_j - 1 ) z_j + sum_{j in N2} x*_j z_j
    *                   sum_{j in N1}          u_j z_j - sum_{j in N2} u_j  z_j > b
    *                                         z_j in {0,1} for all j in N1 & N2
    *
    * 1. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have
    *    positive weights and the constraint is a "<" constraint, by complementing all variables in N1
    *
    *    (KP^SNF_rat)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}          u_j z°_j + sum_{j in N2} u_j  z_j < - b + sum_{j in N1} u_j
    *                                                 z°_j in {0,1} for all j in N1
    *                                                  z_j in {0,1} for all j in N2,
    *    and solve it approximately under consideration of the fixing,
    * or
    * 2. to a knapsack problem in maximization form, such that all variables in the knapsack constraint have
    *    positive integer weights and the constraint is a "<=" constraint, by complementing all variables in N1
    *    and multiplying the constraint by a suitable scalar C
    *
    *    (KP^SNF_int)  max sum_{j in N1} ( 1 - x*_j ) z°_j + sum_{j in N2} x*_j z_j
    *                      sum_{j in N1}        C u_j z°_j + sum_{j in N2} C u_j  z_j <= c
    *                                                   z°_j in {0,1} for all j in N1
    *                                                    z_j in {0,1} for all j in N2,
    *    where
    *      c = floor[ C (- b + sum_{j in N1} u_j ) ]      if frac[ C (- b + sum_{j in N1} u_j ) ] > 0
    *      c =        C (- b + sum_{j in N1} u_j )   - 1  if frac[ C (- b + sum_{j in N1} u_j ) ] = 0
    *    and solve it exactly under consideration of the fixing.
    */
   SCIPdebugMsg(scip, "1. Transform KP^SNF to KP^SNF_rat:\n");

   /* get weight and profit of variables in KP^SNF_rat and check, whether all weights are already integral */
   for( j = 0; j < nitems; j++ )
   {
      transweightsreal[j] = snf->transvarvubcoefs[items[j]];

      if( snf->transvarcoefs[items[j]] == 1 )
      {
         transprofitsreal[j] = 1.0 -  snf->transbinvarsolvals[items[j]];
         SCIPdebugMsg(scip, "     <%d>: j in N1:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j],
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
      else
      {
         transprofitsreal[j] =  snf->transbinvarsolvals[items[j]];
         SCIPdebugMsg(scip, "     <%d>: j in N2:   w_%d = %g, p_%d = %g %s\n", items[j], items[j], transweightsreal[j],
            items[j], transprofitsreal[j], SCIPisIntegral(scip, transweightsreal[j]) ? "" : "  ----> NOT integral");
      }
   }
   /* get capacity of knapsack constraint in KP^SNF_rat */
   transcapacityreal = - snf->transrhs + QUAD_TO_DBL(flowcoverweight) + n1itemsweight; /*lint !e644*/
   SCIPdebugMsg(scip, "     transcapacity = -rhs(%g) + flowcoverweight(%g) + n1itemsweight(%g) = %g\n",
      snf->transrhs, QUAD_TO_DBL(flowcoverweight), n1itemsweight, transcapacityreal);

   /* there exists no flow cover if the capacity of knapsack constraint in KP^SNF_rat after fixing
    * is less than or equal to zero
    */
   if( SCIPisFeasLE(scip, transcapacityreal/10, 0.0) )
   {
      assert(!(*found));
      goto TERMINATE;
   }

   /* KP^SNF_rat has been solved by fixing some variables in advance */
   assert(nitems >= 0);
   if( nitems == 0 )
   {
      /* get lambda = sum_{j in C1} u_j - sum_{j in C2} u_j - rhs */
      SCIPquadprecSumQD(flowcoverweight, flowcoverweight, -snf->transrhs);
      *lambda = QUAD_TO_DBL(flowcoverweight);
      *found = TRUE;
      goto TERMINATE;
   }

   /* Solve the KP^SNF_rat approximately */

   /* initialize number of (non-)solution items, should be changed to a nonnegative number in all possible paths below */
   nsolitems = -1;
   nnonsolitems = -1;

   /* suitable factor C was found*/
   /* solve KP^SNF_rat approximately */
   SCIP_CALL( SCIPsolveKnapsackApproximatelyLT(scip, nitems, transweightsreal, transprofitsreal, transcapacityreal,
         items, solitems, nonsolitems, &nsolitems, &nnonsolitems, NULL) );

   assert(nsolitems != -1);
   assert(nnonsolitems != -1);

   /* build the flow cover from the solution of KP^SNF_rat and KP^SNF_int, respectively and the fixing */
   assert(*nflowcovervars + *nnonflowcovervars + nsolitems + nnonsolitems == snf->ntransvars);
   buildFlowCover(scip, snf->transvarcoefs, snf->transvarvubcoefs, snf->transrhs, solitems, nonsolitems, nsolitems, nnonsolitems, nflowcovervars,
      nnonflowcovervars, flowcoverstatus, QUAD(&flowcoverweight), lambda);
   assert(*nflowcovervars + *nnonflowcovervars == snf->ntransvars);

   *found = SCIPisFeasGT(scip, *lambda, 0.0);

  TERMINATE:
   assert((!*found) || SCIPisFeasGT(scip, *lambda, 0.0));
#ifdef SCIP_DEBUG
   if( *found )
   {
      SCIPdebugMsg(scip, "2. approximate solution:\n");
      for( j = 0; j < snf->ntransvars; j++ )
      {
         if( snf->transvarcoefs[j] == 1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMsg(scip, "     C1: + y_%d [u_%d = %g]\n", j, j, snf->transvarvubcoefs[j]);
         }
         else if( snf->transvarcoefs[j] == -1 && flowcoverstatus[j] == 1 )
         {
            SCIPdebugMsg(scip, "     C2: - y_%d [u_%d = %g]\n", j, j, snf->transvarvubcoefs[j]);
         }
      }
      SCIPdebugMsg(scip, "     flowcoverweight(%g) = rhs(%g) + lambda(%g)\n", QUAD_TO_DBL(flowcoverweight), snf->transrhs, *lambda);
   }
#endif

   /* free data structures */
   SCIPfreeBufferArray(scip, &nonsolitems);
   SCIPfreeBufferArray(scip, &solitems);
   SCIPfreeBufferArray(scip, &transweightsint);
   SCIPfreeBufferArray(scip, &transweightsreal);
   SCIPfreeBufferArray(scip, &transprofitsreal);
   SCIPfreeBufferArray(scip, &itemsint);
   SCIPfreeBufferArray(scip, &items);

   return SCIP_OKAY;
}

#endif

/** evaluate the super-additive lifting function for the lifted simple generalized flowcover inequalities
 *  for a given value \f$ x \in \{ u_j \mid j \in C- \} \f$.
 */
static
SCIP_Real evaluateLiftingFunction(
   SCIP*                 scip,               /**< SCIP data structure */
   LIFTINGDATA*          liftingdata,        /**< lifting data to use */
   SCIP_Real             x                   /**< value where to evaluate lifting function */
   )
{
   SCIP_Real QUAD(tmp);
   SCIP_Real xpluslambda;
   int i;

   assert( liftingdata != NULL );

   xpluslambda = x + liftingdata->lambda;

   i = 0;
   while( i < liftingdata->r && SCIPisGT(scip, xpluslambda, liftingdata->M[i+1]) )
      ++i;

   if( i < liftingdata->t )
   {
      if( SCIPisLE(scip, liftingdata->M[i], x) )
      {
         assert(SCIPisLE(scip, xpluslambda, liftingdata->M[i+1]));
         return i * liftingdata->lambda;
      }

      assert(i > 0 && SCIPisLE(scip, liftingdata->M[i], xpluslambda) && x <= liftingdata->M[i]);

      /* return x - liftingdata->M[i] + i * liftingdata->lambda */
      SCIPquadprecProdDD(tmp, i, liftingdata->lambda);
      SCIPquadprecSumQD(tmp, tmp, x);
      SCIPquadprecSumQD(tmp, tmp, -liftingdata->M[i]);
      return QUAD_TO_DBL(tmp);
   }

   if( i < liftingdata->r )
   {
      assert(!SCIPisInfinity(scip, liftingdata->mp));

      /* p = liftingdata->m[i] - (liftingdata->mp - liftingdata->lambda) - liftingdata->ml; */
      SCIPquadprecSumDD(tmp, liftingdata->m[i], -liftingdata->mp);
      SCIPquadprecSumQD(tmp, tmp, -liftingdata->ml);
      SCIPquadprecSumQD(tmp, tmp, liftingdata->lambda);

      /* p = MAX(0.0, p); */
      if( QUAD_HI(tmp) < 0.0 )
      {
         QUAD_ASSIGN(tmp, 0.0);
      }

      SCIPquadprecSumQD(tmp, tmp, liftingdata->M[i]);
      SCIPquadprecSumQD(tmp, tmp, liftingdata->ml);

      if( SCIPisLT(scip, QUAD_TO_DBL(tmp), xpluslambda) )
         return i * liftingdata->lambda;

      assert(SCIPisFeasLE(scip, liftingdata->M[i], xpluslambda) &&
             SCIPisFeasLE(scip, xpluslambda, liftingdata->M[i] + liftingdata->ml +
             MAX(0.0, liftingdata->m[i] - (liftingdata->mp - liftingdata->lambda) - liftingdata->ml)));

      SCIPquadprecProdDD(tmp, i, liftingdata->lambda);
      SCIPquadprecSumQD(tmp, tmp, x);
      SCIPquadprecSumQD(tmp, tmp, - liftingdata->M[i]);
      return QUAD_TO_DBL(tmp);
   }

   assert(i == liftingdata->r && SCIPisLE(scip, liftingdata->M[liftingdata->r], xpluslambda));

   SCIPquadprecProdDD(tmp, liftingdata->r, liftingdata->lambda);
   SCIPquadprecSumQD(tmp, tmp, x);
   SCIPquadprecSumQD(tmp, tmp, - liftingdata->M[liftingdata->r]);
   return QUAD_TO_DBL(tmp);
}

/** computes
 * \f[
 * (\alpha_j, \beta_j) =
 *    \begin{cases}
 *       (0, 0) &\quad\text{if} M_i \leq u_j \leq M_{i+1} - \lambda \\
 *       (1, M_i - i \lambda) &\quad\text{if} M_i − \lambda < u_j < M_i \\
 *    \end{cases}
 * \f]
 */
static
void getAlphaAndBeta(
   SCIP*                 scip,               /**< SCIP data structure */
   LIFTINGDATA*          liftingdata,        /**< pointer to lifting function struct */
   SCIP_Real             vubcoef,            /**< vub coefficient to get alpha and beta for */
   int*                  alpha,              /**< get alpha coefficient for lifting */
   SCIP_Real*            beta                /**< get beta coefficient for lifting */
   )
{
   SCIP_Real vubcoefpluslambda;
   int i;

   vubcoefpluslambda = vubcoef + liftingdata->lambda;

   i = 0;
   while( i < liftingdata->r && SCIPisGT(scip, vubcoefpluslambda, liftingdata->M[i+1]) )
      ++i;

   if( SCIPisLT(scip, vubcoef, liftingdata->M[i]) )
   {
      SCIP_Real QUAD(tmp);
      assert(liftingdata->M[i] < vubcoefpluslambda);
      *alpha = 1;
      SCIPquadprecProdDD(tmp, -i, liftingdata->lambda);
      SCIPquadprecSumQD(tmp, tmp, liftingdata->M[i]);
      *beta = QUAD_TO_DBL(tmp);
   }
   else
   {
      assert(SCIPisSumLE(scip, liftingdata->M[i], vubcoef));
      assert(i == liftingdata->r || SCIPisLE(scip, vubcoefpluslambda, liftingdata->M[i+1]));
      *alpha = 0;
      *beta = 0.0;
   }
}

/** compute relevant data for performing the sequence independent lifting */
static
SCIP_RETCODE computeLiftingData(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< pointer to SNF relaxation */
   int*                  transvarflowcoverstatus, /**< pointer to store whether non-binary var is in L2 (2) or not (-1 or 1) */
   SCIP_Real             lambda,             /**< lambda */
   LIFTINGDATA*          liftingdata,        /**< pointer to lifting function struct */
   SCIP_Bool*            valid               /**< is the lifting data valid */
   )
{
   int i;
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(sumN2mC2LE);
   SCIP_Real QUAD(sumN2mC2GT);
   SCIP_Real QUAD(sumC1LE);
   SCIP_Real QUAD(sumC2);

#ifndef NDEBUG
   /* for debugging */
   liftingdata->m = NULL;
   liftingdata->M = NULL;
   liftingdata->lambda = SCIP_INVALID;
   liftingdata->t = 0;
   liftingdata->mp = SCIP_INVALID;
#endif

   SCIP_CALL( SCIPallocBufferArray(scip, &liftingdata->m, snf->ntransvars) );

   liftingdata->r = 0;
   QUAD_ASSIGN(sumN2mC2LE, 0.0);
   QUAD_ASSIGN(sumC1LE, 0.0);
   QUAD_ASSIGN(sumN2mC2GT, 0.0);
   QUAD_ASSIGN(sumC2, 0.0);

   liftingdata->mp = SCIPinfinity(scip);

   *valid = FALSE;

   for( i = 0; i < snf->ntransvars; ++i )
   {
      int s = (snf->transvarcoefs[i] + 1) + (transvarflowcoverstatus[i] + 1)/2;

      switch(s)
      {
      case 0: /* var is in N2 \ C2 */
         assert(snf->transvarvubcoefs[i] >= 0.0);
         assert(snf->transvarcoefs[i] == -1 && transvarflowcoverstatus[i] == -1);

         if( SCIPisGT(scip, snf->transvarvubcoefs[i], lambda) )
         {
            SCIPquadprecSumQD(sumN2mC2GT, sumN2mC2GT, snf->transvarvubcoefs[i]);
            liftingdata->m[liftingdata->r++] = snf->transvarvubcoefs[i];
         }
         else
         {
            SCIPquadprecSumQD(sumN2mC2LE, sumN2mC2LE, snf->transvarvubcoefs[i]);
         }
         break;
      case 1: /* var is in C2 */
         assert(snf->transvarvubcoefs[i] > 0.0);
         assert(snf->transvarcoefs[i] == -1 && transvarflowcoverstatus[i] == 1);

         SCIPquadprecSumQD(sumC2, sumC2, snf->transvarvubcoefs[i]);
         break;
      case 3: /* var is in C1 */
         assert(snf->transvarcoefs[i] == 1 && transvarflowcoverstatus[i] == 1);
         assert(snf->transvarvubcoefs[i] > 0.0);

         if( SCIPisGT(scip, snf->transvarvubcoefs[i], lambda) )
         {
            liftingdata->m[liftingdata->r++] = snf->transvarvubcoefs[i];
            liftingdata->mp = MIN(liftingdata->mp, snf->transvarvubcoefs[i]);
         }
         else
         {
            SCIPquadprecSumQD(sumC1LE, sumC1LE, snf->transvarvubcoefs[i]);
         }
         break;
      default:
         assert(s == 2);
         continue;
      }
   }

   if( SCIPisInfinity(scip, liftingdata->mp) )
   {
      SCIPfreeBufferArray(scip, &liftingdata->m);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &liftingdata->M, liftingdata->r + 1) );

   *valid = TRUE;

   SCIPquadprecSumQQ(tmp, sumC1LE, sumN2mC2LE);
   liftingdata->ml = MIN(lambda, QUAD_TO_DBL(tmp));
   SCIPquadprecSumQD(tmp, sumC2, snf->transrhs);
   liftingdata->d1 = QUAD_TO_DBL(tmp);
   SCIPquadprecSumQQ(tmp, tmp, sumN2mC2GT);
   SCIPquadprecSumQQ(tmp, tmp, sumN2mC2LE);
   liftingdata->d2 = QUAD_TO_DBL(tmp);

   SCIPsortDownReal(liftingdata->m, liftingdata->r);

   /* compute M[i] = sum_{i \in [1,r]} m[i] where m[*] is sorted decreasingly and M[0] = 0 */
   QUAD_ASSIGN(tmp, 0.0);
   for( i = 0; i < liftingdata->r; ++i)
   {
      liftingdata->M[i] = QUAD_TO_DBL(tmp);
      SCIPquadprecSumQD(tmp, tmp, liftingdata->m[i]);
   }

   liftingdata->M[liftingdata->r] = QUAD_TO_DBL(tmp);

   SCIP_UNUSED( SCIPsortedvecFindDownReal(liftingdata->m, liftingdata->mp, liftingdata->r, &liftingdata->t) );
   assert(liftingdata->m[liftingdata->t] == liftingdata->mp || SCIPisInfinity(scip, liftingdata->mp)); /*lint !e777*/

   /* compute t largest index sucht that m_t = mp
    * note that liftingdata->m[t-1] == mp due to zero based indexing of liftingdata->m
    */
   ++liftingdata->t;
   while( liftingdata->t < liftingdata->r && liftingdata->m[liftingdata->t] == liftingdata->mp ) /*lint !e777*/
      ++liftingdata->t;

   liftingdata->lambda = lambda;

   return SCIP_OKAY;
}

/** destroy data used for the sequence independent lifting */
static
void destroyLiftingData(
   SCIP*                 scip,               /**< SCIP data structure */
   LIFTINGDATA*          liftingdata         /**< pointer to lifting function struct */
   )
{
   SCIPfreeBufferArray(scip, &liftingdata->M);
   SCIPfreeBufferArray(scip, &liftingdata->m);
}

/** store the simple lifted flowcover cut defined by the given data in the given arrays
 *  the array for storing the cut coefficients must be all zeros
 */
static
SCIP_RETCODE generateLiftedFlowCoverCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SNF_RELAXATION*       snf,                /**< pointer to SNF relaxation */
   SCIP_AGGRROW*         aggrrow,            /**< aggrrow used to construct SNF relaxation */
   int*                  flowcoverstatus,    /**< pointer to store whether variable is in flow cover (+1) or not (-1) */
   SCIP_Real             lambda,             /**< lambda */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   SCIP_Real*            cutrhs,             /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   SCIP_Bool*            success             /**< was the cut successfully generated */
   )
{
   SCIP_Real QUAD(rhs);
   LIFTINGDATA liftingdata;
   int i;

   SCIP_CALL( computeLiftingData(scip, snf, flowcoverstatus, lambda, &liftingdata, success) );
   if( ! *success )
      return SCIP_OKAY;
   assert( liftingdata.m != NULL );
   assert( liftingdata.M != NULL );
   assert( liftingdata.lambda != SCIP_INVALID ); /*lint !e777*/
   assert( liftingdata.r >= 0 );
   assert( liftingdata.t >= 0 );
   assert( liftingdata.mp != SCIP_INVALID ); /*lint !e777*/

   QUAD_ASSIGN(rhs, liftingdata.d1);

   *nnz = 0;

   for( i = 0; i < snf->ntransvars; ++i )
   {
      int s = (snf->transvarcoefs[i] + 1) + (flowcoverstatus[i] + 1)/2;

      switch(s)
      {
      case 0: /* var is in N2 \ C2 */
         if( SCIPisGT(scip, snf->transvarvubcoefs[i], lambda) )
         {
            /* var is in L- */
            if( snf->origbinvars[i] != -1 )
            {
               assert(cutcoefs[snf->origbinvars[i]] == 0.0);
               cutinds[*nnz] = snf->origbinvars[i];
               cutcoefs[snf->origbinvars[i]] = -lambda;
               ++(*nnz);
            }
            else
            {
               SCIPquadprecSumQD(rhs, rhs, lambda);
            }
         }
         else
         {
            /* var is in L-- */
            if( snf->origcontvars[i] != -1 && snf->aggrcoefscont[i] != 0.0 )
            {
               assert(cutcoefs[snf->origcontvars[i]] == 0.0);
               cutinds[*nnz] = snf->origcontvars[i];
               cutcoefs[snf->origcontvars[i]] = -snf->aggrcoefscont[i];
               ++(*nnz);
            }

            if( snf->origbinvars[i] != -1 && snf->aggrcoefsbin[i] != 0.0 )
            {
               assert(cutcoefs[snf->origbinvars[i]] == 0.0);
               cutinds[*nnz] = snf->origbinvars[i];
               cutcoefs[snf->origbinvars[i]] = -snf->aggrcoefsbin[i];
               ++(*nnz);
            }

            SCIPquadprecSumQD(rhs, rhs, snf->aggrconstants[i]);
         }
         break;
      case 1: /* var is in C2 */
      {
         assert(snf->transvarvubcoefs[i] > 0.0);
         assert(snf->transvarcoefs[i] == -1 && flowcoverstatus[i] == 1);

         if( snf->origbinvars[i] != -1 )
         {
            SCIP_Real liftedbincoef = evaluateLiftingFunction(scip, &liftingdata, snf->transvarvubcoefs[i]);
            assert(cutcoefs[snf->origbinvars[i]] == 0.0);
            if( liftedbincoef != 0.0 )
            {
               cutinds[*nnz] = snf->origbinvars[i];
               cutcoefs[snf->origbinvars[i]] = -liftedbincoef;
               ++(*nnz);
               SCIPquadprecSumQD(rhs, rhs, -liftedbincoef);
            }
         }
         break;
      }
      case 2: /* var is in N1 \ C1 */
      {
         int alpha;
         SCIP_Real beta;

         assert(snf->transvarcoefs[i] == 1 && flowcoverstatus[i] == -1);

         getAlphaAndBeta(scip, &liftingdata, snf->transvarvubcoefs[i], &alpha, &beta);
         assert(alpha == 0 || alpha == 1);

         if( alpha == 1 )
         {
            SCIP_Real QUAD(binvarcoef);
            assert(beta > 0.0);

            if( snf->origcontvars[i] != -1 && snf->aggrcoefscont[i] != 0.0 )
            {
               assert(cutcoefs[snf->origcontvars[i]] == 0.0);
               cutinds[*nnz] = snf->origcontvars[i];
               cutcoefs[snf->origcontvars[i]] = snf->aggrcoefscont[i];
               ++(*nnz);
            }

            SCIPquadprecSumDD(binvarcoef, snf->aggrcoefsbin[i], -beta);
            if( snf->origbinvars[i] != -1 )
            {
               SCIP_Real tmp;

               assert(cutcoefs[snf->origbinvars[i]] == 0.0);

               tmp = QUAD_TO_DBL(binvarcoef);
               if( tmp != 0.0 )
               {
                  cutinds[*nnz] = snf->origbinvars[i];
                  cutcoefs[snf->origbinvars[i]] = tmp;
                  ++(*nnz);
               }
            }
            else
            {
               SCIPquadprecSumQQ(rhs, rhs, -binvarcoef);
            }

            SCIPquadprecSumQD(rhs, rhs, -snf->aggrconstants[i]);
         }
         break;
      }
      case 3: /* var is in C1 */
      {
         SCIP_Real bincoef = snf->aggrcoefsbin[i];
         SCIP_Real constant = snf->aggrconstants[i];

         if( snf->origbinvars[i] != -1 && SCIPisGT(scip, snf->transvarvubcoefs[i], lambda) )
         {
            /* var is in C++ */
            SCIP_Real QUAD(tmp);
            SCIP_Real QUAD(tmp2);

            SCIPquadprecSumDD(tmp, snf->transvarvubcoefs[i], -lambda);

            SCIPquadprecSumQD(tmp2, tmp, constant);
            constant = QUAD_TO_DBL(tmp2);

            SCIPquadprecSumQD(tmp2, tmp, -bincoef);
            bincoef = -QUAD_TO_DBL(tmp2);
         }

         if( snf->origbinvars[i] != -1 && bincoef != 0.0 )
         {
            assert(cutcoefs[snf->origbinvars[i]] == 0.0);
            cutinds[*nnz] = snf->origbinvars[i];
            cutcoefs[snf->origbinvars[i]] = bincoef;
            ++(*nnz);
         }

         if( snf->origcontvars[i] != -1 && snf->aggrcoefscont[i] != 0.0 )
         {
            assert(cutcoefs[snf->origcontvars[i]] == 0.0);
            cutinds[*nnz] = snf->origcontvars[i];
            cutcoefs[snf->origcontvars[i]] = snf->aggrcoefscont[i];
            ++(*nnz);
         }

         SCIPquadprecSumQD(rhs, rhs, -constant);
         break;
      }
      default:
         SCIPABORT();
      }
   }

   destroyLiftingData(scip, &liftingdata);

   {
      SCIP_ROW** rows = SCIPgetLPRows(scip);
      for( i = 0; i < aggrrow->nrows; ++i )
      {
         SCIP_ROW* row;
         SCIP_Real rowlhs;
         SCIP_Real rowrhs;
         SCIP_Real slackub;
         SCIP_Real slackcoef;

         slackcoef = aggrrow->rowweights[i] * aggrrow->slacksign[i];
         assert(slackcoef != 0.0);

         /* positive slack was implicitly handled in flow cover separation */
         if( slackcoef > 0.0 )
            continue;

         row = rows[aggrrow->rowsinds[i]];

         /* add the slack's definition multiplied with its coefficient to the cut */
         SCIP_CALL( varVecAddScaledRowCoefs(cutinds, cutcoefs, nnz, row, -aggrrow->rowweights[i]) );

         /* retrieve sides of row */
         rowlhs = row->lhs - row->constant;
         rowrhs = row->rhs - row->constant;

         if( row->integral )
         {
            rowrhs = SCIPfloor(scip, rowrhs);
            rowlhs = SCIPceil(scip, rowlhs);
         }

         slackub = rowrhs - rowlhs;

         /* move slack's constant to the right hand side, and add lambda to the right hand side if the
          * upper bound of the slack is larger than lambda, since then an artifical binary variable
          * for the slack would get coefficient -lambda
          */
         if( aggrrow->slacksign[i] == +1 )
         {
            SCIP_Real rhsslack;
            /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a^_r * (rhs - c) to the right hand side */
            assert(!SCIPisInfinity(scip, row->rhs));

            rhsslack = rowrhs - SCIPgetRowMinActivity(scip, row);
            slackub = -aggrrow->rowweights[i] * MIN(rhsslack, slackub);

            if( SCIPisGE(scip, slackub, lambda) )
               SCIPquadprecSumQD(rhs, rhs, lambda);

            SCIPquadprecSumQD(rhs, rhs, -aggrrow->rowweights[i] * rowrhs);
         }
         else
         {
            SCIP_Real lhsslack;
            /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a^_r * (c - lhs) to the right hand side */
            assert(!SCIPisInfinity(scip, -row->lhs));

            lhsslack = SCIPgetRowMaxActivity(scip, row) - rowlhs;
            slackub = aggrrow->rowweights[i] * MIN(lhsslack, slackub);

            if( SCIPisGE(scip, slackub, lambda) )
               SCIPquadprecSumQD(rhs, rhs, lambda);

            SCIPquadprecSumQD(rhs, rhs, -aggrrow->rowweights[i] * rowlhs);
         }
      }
   }

   *cutrhs = QUAD_TO_DBL(rhs);

   /* relax rhs to zero, if it's very close to 0 */
   if( *cutrhs < 0.0 && *cutrhs >= -SCIPepsilon(scip) )
      *cutrhs = 0.0;

   return SCIP_OKAY;
}

/** calculates a lifted simple generalized flow cover cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in the cut.
 *  For further details we refer to:
 *
 *  Gu, Z., Nemhauser, G. L., & Savelsbergh, M. W. (1999). Lifted flow cover inequalities for mixed 0-1 integer programs.
 *  Mathematical Programming, 85(3), 439-467.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcalcFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute flow cover cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store the efficacy of the cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid cut was returned */
   )
{
   int i;
   int nvars;
   SCIP_Bool localbdsused;
   SNF_RELAXATION snf;
   SCIP_Real lambda;
   SCIP_Real* tmpcoefs;
   int *transvarflowcoverstatus;
   int nflowcovervars;
   int nnonflowcovervars;

   nvars = SCIPgetNVars(scip);

   *success = FALSE;

   /* get data structures */
   SCIP_CALL( SCIPallocBufferArray(scip, &transvarflowcoverstatus, nvars) );
   SCIP_CALL( allocSNFRelaxation(scip,  &snf, nvars) );

   SCIPdebug( printCutQuad(scip, sol, aggrrow->vals, QUAD(aggrrow->rhs), aggrrow->inds, aggrrow->nnz, FALSE, aggrrow->local) );

   SCIP_CALL( constructSNFRelaxation(scip, sol, boundswitch, allowlocal, aggrrow->vals, QUAD(aggrrow->rhs), aggrrow->inds, aggrrow->nnz, &snf, success, &localbdsused) );

   if( ! *success )
   {
      goto TERMINATE;
   }

   *cutislocal = aggrrow->local || localbdsused;

   /* initialize lambda because gcc issues a stupid warning */
   lambda = 0.0;
   SCIP_CALL( getFlowCover(scip, &snf, &nflowcovervars, &nnonflowcovervars, transvarflowcoverstatus, &lambda, success) );

   if( ! *success )
   {
      goto TERMINATE;
   }

   SCIP_CALL( SCIPallocCleanBufferArray(scip, &tmpcoefs, nvars) );

   SCIP_CALL( generateLiftedFlowCoverCut(scip, &snf, aggrrow, transvarflowcoverstatus, lambda, tmpcoefs, cutrhs, cutinds, cutnnz, success) );
   SCIPdebugMsg(scip, "computed flowcover_%lli_%i:\n", SCIPgetNLPs(scip), SCIPgetNCuts(scip));

   /* if success is FALSE generateLiftedFlowCoverCut wont have touched the tmpcoefs array so we dont need to clean it then */
   if( *success )
   {
      if( postprocess )
      {
         SCIP_CALL( postprocessCut(scip, *cutislocal, cutinds, tmpcoefs, cutnnz, cutrhs, success) );
      }
      else
      {
         SCIP_Real QUAD(rhs);

         QUAD_ASSIGN(rhs, *cutrhs);
         *success = ! removeZeros(scip, SCIPsumepsilon(scip), *cutislocal, tmpcoefs, QUAD(&rhs), cutinds, cutnnz);
         *cutrhs = QUAD_TO_DBL(rhs);
      }

      if( *success )
      {
         /* store cut sparse and calculate efficacy */
         for( i = 0; i < *cutnnz; ++i )
         {
            int j = cutinds[i];
            assert(tmpcoefs[j] != 0.0);
            cutcoefs[i] = tmpcoefs[j];
            tmpcoefs[j] = 0.0;
         }

         if( cutefficacy != NULL )
            *cutefficacy = calcEfficacy(scip, sol, cutcoefs, *cutrhs, cutinds, *cutnnz);

         if( cutrank != NULL )
            *cutrank = aggrrow->rank + 1;
      }
      else
      {
         /* clean buffer array */
         for( i = 0; i < *cutnnz; ++i )
         {
            int j = cutinds[i];
            assert(tmpcoefs[j] != 0.0);
            tmpcoefs[j] = 0.0;
         }
      }
   }

   SCIPfreeCleanBufferArray(scip, &tmpcoefs);

  TERMINATE:
   destroySNFRelaxation(scip, &snf);
   SCIPfreeBufferArray(scip, &transvarflowcoverstatus);

   return SCIP_OKAY;
}

/* =========================================== knapsack cover =========================================== */

/** Relax the row to a possibly fractional knapsack row containing no integer or continuous variables
 *  and only having positive coefficients for binary variables. General integer and continuous variables
 *  are complemented with variable or simple bounds such that their coefficient becomes positive and then
 *  it is relaxed to zero.
 *  All remaining binary variables are complemented with simple upper or lower bounds such that their
 *  coefficient becomes positive.
 */
static
SCIP_RETCODE cutsTransformKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            localbdsused,       /**< pointer to store whether local bounds were used in transformation */
   SCIP_Bool*            success             /**< stores whether the row could successfully be transformed into a knapsack constraint.
                                              *   Returns FALSE in case a continuous or general integer variable is unbounded in the
                                              *   required direction. */
   )
{
   SCIP_Real* bestbds;
   int i;
   int aggrrowbinstart;
   int firstnonbinvar;
   SCIP_VAR** vars;

   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(success != NULL);
   assert(localbdsused != NULL);

   *success = FALSE;

   /* allocate temporary memory to store best bounds and bound types */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestbds, 2*(*nnz)) );

   /* start with continuous variables, because using variable bounds can affect the untransformed binary
    * variables, and these changes have to be incorporated in the transformation of the binary variables
    * (binary variables have the smallest problem indices!)
    */
   SCIPsortDownInt(cutinds, *nnz);

   vars = SCIPgetVars(scip);
   firstnonbinvar = SCIPgetNBinVars(scip);

   /* determine best bounds for the continuous and general integer variables such that they will have
    * a positive coefficient in the transformation */
   for( i = 0; i < *nnz && cutinds[i] >= firstnonbinvar; ++i )
   {
      SCIP_Real QUAD(coef);
      int v = cutinds[i];

      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      if( QUAD_TO_DBL(coef) > 0.0 )
      {
         /* find closest lower bound in standard lower bound or variable lower bound for continuous variable
          * so that it will have a positive coefficient */
         SCIP_CALL( findBestLb(scip, vars[v], sol, SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS ? 1 : 0, allowlocal, bestbds + i, boundtype + i) );

         /* cannot transform into knapsack */
         if( SCIPisInfinity(scip, -bestbds[i]) )
            goto TERMINATE;

         varsign[i] = +1;
      }
      else if( QUAD_TO_DBL(coef) < 0.0 )
      {
         /* find closest upper bound in standard upper bound or variable upper bound for continuous variable
          * so that it will have a positive coefficient */
         SCIP_CALL( findBestUb(scip, vars[v], sol, SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS ? 1 : 0, allowlocal, bestbds + i, boundtype + i) );

          /* cannot transform into knapsack */
         if( SCIPisInfinity(scip, bestbds[i]) )
            goto TERMINATE;

         varsign[i] = -1;
      }
   }

   /* remember start of integer variables in the aggrrow */
   aggrrowbinstart = i;

   /* perform bound substitution for continuous variables */
   for( i = 0; i < aggrrowbinstart; ++i )
   {
      SCIP_Real QUAD(coef);
      int v = cutinds[i];

      performBoundSubstitution(scip, cutinds, cutcoefs, QUAD(cutrhs), nnz, varsign[i], boundtype[i], bestbds[i], v, localbdsused);

      /* relax non-binary coefficient to zero after bound substitution */
      QUAD_ASSIGN(coef, 0.0);
      QUAD_ARRAY_STORE(cutcoefs, v, coef);
   }

   assert(i == aggrrowbinstart);

   /* remove non-binary variables because their coefficients have been set to zero after bound substitution */
   if( aggrrowbinstart != 0 )
   {
      *nnz -= aggrrowbinstart;
      BMSmoveMemoryArray(cutinds, cutinds + aggrrowbinstart, *nnz);
   }
   i = 0;

   /* after doing bound substitution of non-binary vars, some coefficients of binary vars might have changed, so here we
    * remove the ones that became 0 if any; also, we need that all remaining binary vars have positive coefficients,
    * thus we perform bound substitution with simple bounds (i.e. complementing) to achieve this.
    */
   while( i < *nnz )
   {
      SCIP_Real QUAD(coef);
      SCIP_Real bestlb;
      SCIP_Real bestub;
      SCIP_Bool setzero;
      int v = cutinds[i];

      assert(SCIPvarGetType(vars[v]) == SCIP_VARTYPE_BINARY && !SCIPvarIsImpliedIntegral(vars[v]));

      assert(v < firstnonbinvar);
      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      /* due to variable bound usage for bound substitution of continuous variables cancellation may have occurred */
      if( EPSZ(QUAD_TO_DBL(coef), QUAD_EPSILON) )
      {
         /* do not increase i, since last element is copied to the i-th position */
         setzero = TRUE;
      }
      else
      {
         /* perform bound substitution */
         if( QUAD_TO_DBL(coef) < 0.0 )
         {
            SCIP_CALL( findBestUb(scip, vars[v], sol, 0, allowlocal, &bestub, boundtype + i) );

            if( SCIPisZero(scip, bestub) )
            {
               /* binary variable is fixed to zero */
               setzero = TRUE;
               *localbdsused = *localbdsused || (boundtype[i] == -2);
            }
            else
            {
               varsign[i] = -1;

               performBoundSubstitutionSimple(scip, cutcoefs, QUAD(cutrhs), boundtype[i], bestub, v, localbdsused);
               QUAD_ARRAY_STORE(cutcoefs, v, -coef);
               setzero = FALSE;
            }
         }
         else
         {
            SCIP_CALL( findBestLb(scip, vars[v], sol, 0, allowlocal, &bestlb, boundtype + i) );

            if( !SCIPisZero(scip, bestlb) )
            {
               /* binary variable is fixed to one */
               performBoundSubstitutionSimple(scip, cutcoefs, QUAD(cutrhs), boundtype[i], bestlb, v, localbdsused);
               setzero = TRUE;
            }
            else
            {
               varsign[i] = +1;
               setzero = FALSE;
            }
         }

         assert(boundtype[i] == -1 || boundtype[i] == -2);
      }

      /* increase i or remove zero coefficient (i.e. var with 0 coef) by shifting last nonzero to current position */
      if( setzero )
      {
         QUAD_ASSIGN(coef, 0.0);
         QUAD_ARRAY_STORE(cutcoefs, v, coef);
         --(*nnz);
         cutinds[i] = cutinds[*nnz];
      }
      else
         ++i;
   }

   /* relax rhs to zero if it is close to but slightly below zero */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   *success = TRUE;
  TERMINATE:
   /*free temporary memory */
   SCIPfreeBufferArray(scip, &bestbds);

   return SCIP_OKAY;
}

/** determines the initial cover for the given (fractional) knapsack row */
static
SCIP_Bool computeInitialKnapsackCover(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   SCIP_Real             cutrhs,             /**< pointer to the right hand side of the cut */
   int                   cutnnz,             /**< pointer to the number of non-zeros in the cut */
   int*                  varsign,            /**< sign of coefficients for each nonzero in the row be transformation */
   int*                  coverstatus,        /**< array to return the coverstatus for each variable in the  knapsack row */
   int*                  coverpos,           /**< position of nonzero in the knapsack row for each variable in the cover */
   SCIP_Real*            covervals,          /**< coefficient value of each variable in the cover */
   int*                  coversize,          /**< pointer to return number of variables in the cover;
                                              *   matches the length of the associated arrays */
   QUAD(SCIP_Real*       coverweight)        /**< pointer to return the weight of the cover;
                                              *   the weight is the sum of the coefficient values of variables in the cover */
   )
{
   SCIP_VAR** vars;
   int k;
   int j;
   QUAD_ASSIGN(*coverweight, 0);
   *coversize = 0;
   j = cutnnz-1;
   vars = SCIPgetVars(scip);

   for( k = 0; k < cutnnz; ++k )
   {
      SCIP_Real solval;
      int v = cutinds[k];
      SCIP_Real QUAD(coef);
      QUAD_ARRAY_LOAD(coef, cutcoefs, v);

      solval = SCIPgetSolVal(scip, sol, vars[v]);
      if( varsign[k] == -1 )
         solval = 1 - solval;

      if( SCIPisFeasEQ(scip, solval, 1.0) )
      {
         /* every variable with solution value 1 is forced into the cover */
         coverpos[*coversize] = k;
         covervals[*coversize] = QUAD_TO_DBL(coef);
         coverstatus[k] = 1;
         *coversize += 1;
         SCIPquadprecSumQQ(*coverweight, *coverweight, coef);
      }
      else
      {
         coverpos[j] = k;
         covervals[j] = solval * QUAD_TO_DBL(coef);
         coverstatus[k] = 0;
         j -= 1;
      }
   }

   /* Use these two arrays to sort the variables by decreasing contribution
    * and pick them greedily in the while loop below until they are a cover.
    * Since the cover does not need to be minimal we do not need to remove any of the
    * variables with a high activity contribution even if they are not necessary after
    * picking the last variable.
    */
   SCIPsortDownRealInt(covervals + (*coversize), coverpos + (*coversize), cutnnz - (*coversize));

   /* overwrite covervals with the coefficients of the variables in the cover
    * as we need to sort decreasingly by those again for the lifting
    */
   while( *coversize < cutnnz &&
          SCIPisFeasLE(scip, QUAD_TO_DBL(*coverweight), cutrhs) )
   {
      int v;
      SCIP_Real QUAD(coef);
      k = coverpos[*coversize];
      v = cutinds[k];
      coverstatus[k] = 1;
      QUAD_ARRAY_LOAD(coef, cutcoefs, v);
      covervals[*coversize] = QUAD_TO_DBL(coef);
      SCIPquadprecSumQQ(*coverweight, *coverweight, coef);
      *coversize += 1;
   }

   /* there is no cover */
   if( SCIPisFeasLE(scip, QUAD_TO_DBL(*coverweight), cutrhs) || *coversize == 0 )
      return FALSE;

   SCIPdebugMsg(scip, "coverweight is %g and right hand side is %g\n", QUAD_TO_DBL(*coverweight), cutrhs);
   assert(*coversize > 0);

   return TRUE;
}

/** prepares the data needed to evaluate the lifting function */
static
void prepareLiftingData(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*            cutcoefs,           /**< array of the non-zero coefficients in the cut */
   int*                  cutinds,            /**< array of the problem indices of variables with a non-zero coefficient in the cut */
   QUAD(SCIP_Real        cutrhs),            /**< pointer to the right hand side of the cut */
   int*                  coverpos,           /**< position of nonzero in the knapsack row for each variable in the cover */
   int                   coversize,          /**< number of variables in the cover */
   QUAD(SCIP_Real        coverweight),       /**< weight of cover */
   SCIP_Real*            covervals,          /**< coefficient value of each variable in the cover;
                                              *   on output stores the running sum of S^-(*) values */
   int*                  coverstatus,        /**< coverstatus for each variable in the cover. After calling this function
                                              *   variables in C^- will have the value -1, variables in C^+ the value 1,
                                              *   and all variables outside the cover keep the value 0. */
   QUAD(SCIP_Real*       abar),              /**< pointer to store the reciprocal value of \bar{a} */
   int*                  cplussize           /**< pointer to store the size of C^+ */
   )
{
   int k;
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(sigma);

   /* Now compute \bar{a}, the unique rational number such that for the cover C it holds that
    * b = \sum_{a_i \in C} \min(\bar{a}, a_i).
    * For that we need to sort by decreasing coefficients of the variables in the cover.
    * After the sorting the covervals array is free to be reused.
    */
   SCIPsortDownRealInt(covervals, coverpos, coversize);

   /* Now follows Algorithm 1 in the paper to compute \bar{a} */

   /* set \bar{a} = l_1 */
   QUAD_ARRAY_LOAD(*abar, cutcoefs, cutinds[coverpos[0]]);
   SCIPquadprecSumQQ(sigma, coverweight, -cutrhs);

   for( k = 1; k < coversize; ++k )
   {
      SCIP_Real QUAD(lkplus1);
      SCIP_Real QUAD(kdelta);

      /* load next coefficient l_{k+1} in sorted order of cover */
      QUAD_ARRAY_LOAD(lkplus1, cutcoefs, cutinds[coverpos[k]]);

      /* Let \delta = \bar{a} - l_{k+1} and compute k * \delta */
      SCIPquadprecSumQQ(kdelta, *abar, -lkplus1);
      SCIPquadprecProdQD(kdelta, kdelta, k);

      /* Set tmp = k * \delta - \sigma to check condition k * \delta < \sigma by tmp < 0 */
      SCIPquadprecSumQQ(tmp, kdelta, -sigma);
      if( QUAD_TO_DBL(tmp) < 0.0 )
      {
         /* Set \bar{a} = l_{k+1} and \sigma = \sigma - k*\delta */
         QUAD_ASSIGN_Q(*abar, lkplus1);
         SCIPquadprecSumQQ(sigma, sigma, -kdelta);
      }
      else
      {
         /* Set \bar{a} = \bar{a} - \sigma / k and \sigma = 0; break; */
         SCIP_Real minusoneoverk = -1.0 / k;
         SCIPquadprecProdQD(sigma, sigma, minusoneoverk);
         SCIPquadprecSumQQ(*abar, *abar, sigma);
         QUAD_ASSIGN(sigma, 0.0);
         break;
      }
   }

   if( QUAD_TO_DBL(sigma) > 0.0 )
   {
      SCIP_Real oneoverc = 1.0 / coversize;
      SCIPquadprecProdQD(*abar, cutrhs, oneoverc);
   }

   /* now we partition C into C^+ and C^-, where C^+ are all the elements of C whose weight is strictly larger than
    * \bar{a} and C^- the rest.  If a_i are the weights of the elements in C, let a_i^- = min(a_i, \bar{a}) We also
    * compute S^-(h) = sum of the h largest a_i^- and store S^-(h+1) in in covervals[h], for k = 0, ..., coversize - 1
    * (S^-(0) = 0 so it doesn't need to be stored; we use S to compute the lifted cut, see below)
    * we remember which elements of C^- in coverstatus, so that element in C^+ have coverstatus 1 and
    * elements in C^- have coverstatus -1 (elements not in C have coverstatus 0)
    */
   QUAD_ASSIGN(tmp, 0.0);
   *cplussize = 0;
   for( k = 0; k < coversize; ++k )
   {
      SCIP_Real QUAD(coef);
      SCIP_Real QUAD(coefminusabar);

      QUAD_ARRAY_LOAD(coef, cutcoefs, cutinds[coverpos[k]]);
      SCIPquadprecSumQQ(coefminusabar, coef, -*abar);
      if( QUAD_TO_DBL(coefminusabar) > 0.0 )
      {
         /* coefficient is in C^+ because it is greater than \bar{a} and contributes only \bar{a} to the sum */
         SCIPquadprecSumQQ(tmp, tmp, *abar);

         /* rather be on the safe side in numerical corner cases and relax the coefficient to exactly \bar{a}.
          * In that case the coefficient is not treated as in C^+ but as being <= \bar{a} and therefore in C^-.
          */
         if( QUAD_TO_DBL(coefminusabar) > SCIPfeastol(scip) )
            ++(*cplussize);
         else
            coverstatus[coverpos[k]] = -1;
      }
      else
      {
         /* coefficient is in C^- because it is smaller or equal to \bar{a} */
         coverstatus[coverpos[k]] = -1;
         SCIPquadprecSumQQ(tmp, tmp, coef);
      }
      covervals[k] = QUAD_TO_DBL(tmp);
      SCIPdebugMsg(scip, "S^-(%d) = %g\n", k + 1, covervals[k]);
   }

   /* set abar to its reciprocal for faster computation of the lifting coefficients */
   SCIPquadprecDivDQ(*abar, 1, *abar);
}

/** evaluate the lifting function based on the given values */
static
SCIP_Real evaluateLiftingFunctionKnapsack(
   SCIP*                 scip,               /**< SCIP datastructure */
   QUAD(SCIP_Real        x),                 /**< value to evaluate the lifting function at */
   QUAD(SCIP_Real        abar),              /**< the reciprocal value of \bar{a} */
   SCIP_Real*            covervals,          /**< the running sum of S^-(*) values */
   int                   coversize,          /**< the size of the cover */
   int                   cplussize,          /**< the size of C^+ */
   SCIP_Real*            scale               /**< pointer to update the scale to integrality when a fractional value is returned */
   )
{
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(hfrac);
   SCIP_Real cutcoef;
   SCIP_Real hreal;
   int h;

   /* the lifted value is at least the coeficient (a_k) divided by \bar{a} because the largest value
    * contributed to the running sum stored in C is \bar{a}
    * therefore we start the search for the correct h at floor(a_k / \bar{a})
    */

   SCIPdebugMsg(scip, "coef is %g, coversize is %d\n", QUAD_TO_DBL(x), coversize );

   SCIPquadprecProdQQ(hfrac, x, abar);

   /* if the coefficient is below \bar{a}, i.e. a / \bar{a} < 1 then g(a_k) = 0, otherwise g(a_k) > 0 */
   if( QUAD_TO_DBL(hfrac) < 1 )
      return 0.0;

   /* we perform h = MIN(h, coversize) in floating-point first because on some instances h was seen to exceed the range
    * of int */
   hreal = SCIPfloor(scip, QUAD_TO_DBL(hfrac));
   if( hreal > (SCIP_Real)coversize )
      h = coversize;
   else
      h = (int)hreal;

   SCIPquadprecSumQD(hfrac, hfrac, -h);

   assert(h > 0);
   if( h < cplussize && ABS(QUAD_TO_DBL(hfrac)) <= QUAD_EPSILON )
   {
      /* cutcoef can be increased by 0.5 because it is a multiple of \bar{a}
       * (This is the first non-dominated lifting function presented in the paper)
       */
      cutcoef = 0.5;
      *scale = 2.0;
   }
   else
      cutcoef = 0.0;

   /* decrease by one to make sure rounding errors or coefficients that are larger than the right hand side by themselves
    * did not push h too far */
   h--;

   /* now increase coefficient to its lifted value based on its size relative to the S^- values.
    * The coefficient a_i is lifted to the unique integer h such that S^-(h) < a_i <= S^-(h+1).
    * (todo: variables that have a coefficient above the right hand side can get an arbitrarily large coefficient but can
    *  also be trivially fixed using the base row. Currently they get the coefficient |C| which is 1 above the right hand
    *  side in the cover cut so that they can still be trivially fixed by propagating the cover cut.
    *  We do not want to apply fixings here though because the LP should stay flushed during separation.
    *  Possibly add a parameter to return additional fixings to the caller of the SCIPcalc*() functions in here
    *  and the caller can add them as cuts to the sepastore or we add them to the sepastore here?)
    */
   while( h < coversize )
   {
      SCIPquadprecSumQD(tmp, x, -covervals[h]); /* recall: covervals[h] = S^-(h+1) */
      /* compare with standard epsilon tolerance since computation involves abar, which is computed like an activity */
      if( !SCIPisPositive(scip, QUAD_TO_DBL(tmp)) )
         break;

      ++h;
   }

   cutcoef += h;

   SCIPdebugMsg(scip, "x is %g, coversize is %d, h is %d\n", QUAD_TO_DBL(x), coversize, h );
   /* the lifted coefficient is h increased possibly by 0.5 for the case checked above */
   SCIPdebugMsg(scip, "lifted coef %g < %g <= %g to %g\n", h == 0 ? 0 : covervals[h-1], QUAD_TO_DBL(x),
         covervals[h], cutcoef);

   return cutcoef;
}

/** calculates a lifted knapsack cover cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in the cut.
 *  For further details we refer to:
 *
 *  Letchford, A. N., & Souli, G. (2019). On lifted cover inequalities: A new lifting procedure with unusual properties.
 *  Operations Research Letters, 47(2), 83-87.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcalcKnapsackCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute flow cover cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store the efficacy of the cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid cut was returned */
   )
{
   int* varsign;
   int* boundtype;
   int* coverstatus;
   int* coverpos;
   int* tmpinds;
   SCIP_Real* tmpcoefs;
   SCIP_Real* covervals;
   SCIP_Real QUAD(rhs);
   SCIP_Real QUAD(coverweight);
   SCIP_Real QUAD(abar);
   SCIP_Bool transformed;
   SCIP_Bool local;
   SCIP_Real efficacy;
   SCIP_Real scale;
   int k;
   int nvars;
   int coversize;
   int cplussize;
   int nnz;

   assert(scip != NULL);
   assert(aggrrow != NULL);
   assert(cutcoefs != NULL);
   assert(cutrhs != NULL);
   assert(cutinds != NULL);
   assert(cutnnz != NULL);
   assert(cutefficacy != NULL);
   assert(cutislocal != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( aggrrow->nnz == 0 )
      return SCIP_OKAY;

   for( k = 0; k < aggrrow->nrows; ++k )
   {
      /* cannot handle negative slack variables */
      if( aggrrow->rowweights[k] * aggrrow->slacksign[k] < 0 )
         return SCIP_OKAY;
   }

   /* allocate temporary memory */
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtype, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coverstatus, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &covervals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coverpos, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpinds, nvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &tmpcoefs, QUAD_ARRAY_SIZE(nvars)) );

   /* initialize cut with aggregation */
   nnz = aggrrow->nnz;
   QUAD_ASSIGN_Q(rhs, aggrrow->rhs);

   BMScopyMemoryArray(tmpinds, aggrrow->inds, nnz);

   for( k = 0; k < nnz; ++k )
   {
      SCIP_Real QUAD(coef);
      int j = tmpinds[k];

      QUAD_ARRAY_LOAD(coef, aggrrow->vals, j);

      QUAD_HI(coef) = NONZERO(QUAD_HI(coef));
      assert(QUAD_HI(coef) != 0.0);

      QUAD_ARRAY_STORE(tmpcoefs, j, coef);
   }
   SCIPdebugMsg(scip, "Computing lifted knapsack cover for ");
   SCIPdebug(printCutQuad(scip, NULL, tmpcoefs, QUAD(rhs), tmpinds, nnz, FALSE, FALSE));

   /* Transform aggregated row into a (fractional, i.e. with possibly fractional weights) knapsack constraint.
    * Uses simple or variable lower or upper bounds to relax out continuous and general integers
    * so that only binary variables remain and complements those such that they have a positive coefficient.
    */
   local = aggrrow->local;
   SCIP_CALL( cutsTransformKnapsackCover(scip, sol, allowlocal,
         tmpcoefs, QUAD(&rhs), tmpinds, &nnz, varsign, boundtype, &local, &transformed) );

   assert(allowlocal || !local);

   if( !transformed )
      goto TERMINATE;

   SCIPdebugMsg(scip, "Transformed knapsack relaxation ");
   SCIPdebug(printCutQuad(scip, NULL, tmpcoefs, QUAD(rhs), tmpinds, nnz, FALSE, FALSE));

   if( !computeInitialKnapsackCover(scip, sol, tmpcoefs, tmpinds, QUAD_TO_DBL(rhs), nnz, varsign, coverstatus,
            coverpos, covervals, &coversize, QUAD(&coverweight)) )
      goto TERMINATE;

   SCIPdebugMsg(scip, "coverweight is %g and right hand side is %g\n", QUAD_TO_DBL(coverweight), QUAD_TO_DBL(rhs));
   assert(coversize > 0);

   /* by default do not scale the cut */
   scale = 1.0;

   if( coversize == 1 )
   {
      SCIP_Real QUAD(tmp);
      /* cover is trivial, return the fixing as cut */
      QUAD_ASSIGN(tmp, 0.0);
      for( k = 0; k < nnz; ++k )
      {
         if( coverstatus[k] == 0 )
         {
            QUAD_ARRAY_STORE(tmpcoefs, tmpinds[k], tmp);
         }
         else
         {
            tmpinds[0] = tmpinds[k];
            varsign[0] = varsign[k];
         }
      }

      nnz = 1;
      if( varsign[0] == -1 )
      {
         QUAD_ASSIGN(rhs, -1.0);
         QUAD_ASSIGN(tmp, -1.0);
      }
      else
      {
         QUAD_ASSIGN(rhs, 0.0);
         QUAD_ASSIGN(tmp, 1.0);
      }

      QUAD_ARRAY_STORE(tmpcoefs, tmpinds[0], tmp);
   }
   else
   {
      SCIP_Real QUAD(tmp);

      /* compute lifted cover inequality:
       * sum_{i \in C^-) x_i + sum_{i \in N \ C^-) g(a_i) x_i <= c - 1
       * where g(z) is equal to
       *   - 0 if z is 0 (irrelevant as there shouldn't be element with weight 0 in the knapsack)
       *   - h + 1/2 if z = k * \bar{a} for some integer k \in [1, |C^+| - 1] and S^-(h) < z <= S^-(h+1) for some h = 0, ..., coversize -1
       *   - h if S^-(h) < z <= S^-(h+1) for some h = 0, ..., coversize -1
       * the function S^- is defined above. Note that S^-(0) = 0
       * we store the cut coefficients in tmpcoef
       */

      SCIPdebugMsg(scip, "call prepareLiftingData: \n");
      /* prepare data required to evaluate lifting function */
      prepareLiftingData(scip, tmpcoefs, tmpinds, QUAD(rhs), coverpos, coversize,
            QUAD(coverweight), covervals, coverstatus, QUAD(&abar), &cplussize);

      /* compute lifted cover inequality */
      QUAD_ASSIGN(rhs, (coversize - 1));
      for( k = 0; k < nnz; )
      {
         SCIP_Real cutcoef;
         if( coverstatus[k] == -1 )
         { /* variables in C^- get the coefficients 1 */
            cutcoef = 1.0;
         }
         else
         { /* variables is either in C^+ or not in the cover and its coefficient value is computed with the lifing function */
            SCIP_Real QUAD(coef);

            SCIPdebugMsg(scip, "load QUAD(coef) from tmpcoefs[tmpinds[k] = %d]\n",tmpinds[k]);
            QUAD_ARRAY_LOAD(coef, tmpcoefs, tmpinds[k]);

            SCIPdebugMsg(scip, "coef is QUAD_HI=%g, QUAD_LO=%g, QUAD_TO_DBL = %g\n",QUAD_HI(coef), QUAD_LO(coef), QUAD_TO_DBL(coef));

            SCIPdebugMsg(scip, "call evaluateLiftingFunctionKnapsack:\n");
            cutcoef = evaluateLiftingFunctionKnapsack(scip, QUAD(coef), QUAD(abar), covervals, coversize, cplussize, &scale);

            /* if the coefficient value is zero then remove the nonzero entry and continue */
            if( cutcoef == 0.0 )
            {
               QUAD_ASSIGN(tmp, 0.0);
               QUAD_ARRAY_STORE(tmpcoefs, tmpinds[k], tmp);
               --nnz;
               coverstatus[k] = coverstatus[nnz];
               tmpinds[k] = tmpinds[nnz];
               varsign[k] = varsign[nnz];
               continue;
            }
         }

         /* directly undo the complementation before storing back the coefficient */
         if( varsign[k] == -1 )
         {
            /* variable was complemented so we have cutcoef * (1-x) = cutcoef - cutcoef * x.Thus we need to adjust the rhs
             * to rhs - cutcoef and flip the sign of cutcoef */
            cutcoef = -cutcoef;
            SCIPquadprecSumQD(rhs, rhs, cutcoef);
         }

         QUAD_ASSIGN(tmp, cutcoef);
         QUAD_ARRAY_STORE(tmpcoefs, tmpinds[k], tmp);

         ++k;
      }
   }

   /* calculate the efficacy of the computed cut and store the success flag if the efficacy exceeds the
    * one stored in the cutefficacy variable by the caller
    */
   efficacy = calcEfficacyDenseStorageQuad(scip, sol, tmpcoefs, QUAD_TO_DBL(rhs), tmpinds, nnz);
   *success = SCIPisGT(scip, efficacy, *cutefficacy);

   SCIPdebugMsg(scip, "FINAL LCI:");
   SCIPdebug(printCutQuad(scip, sol, tmpcoefs, QUAD(rhs), tmpinds, nnz, FALSE, FALSE));

   if( *success )
   {
      /* return the cut into the given arrays/pointers */
      *cutislocal = local;
      *cutrhs = scale * QUAD_TO_DBL(rhs);
      *cutnnz = nnz;

      /* store cut in given array in sparse representation and clean buffer array */
      for( k = 0; k < nnz; ++k )
      {
         SCIP_Real QUAD(coef);
         int j = tmpinds[k];

         QUAD_ARRAY_LOAD(coef, tmpcoefs, j);
         assert(QUAD_HI(coef) != 0.0);

         cutcoefs[k] = scale * QUAD_TO_DBL(coef);
         cutinds[k] = j;
         QUAD_ASSIGN(coef, 0.0);
         QUAD_ARRAY_STORE(tmpcoefs, j, coef);
      }

      assert( cutefficacy != NULL );
      /* calculate efficacy again to make sure it matches the coefficients after they where rounded to double values
       * and after the cleanup and postprocessing step was applied. */
      *cutefficacy = calcEfficacy(scip, sol, cutcoefs, *cutrhs, cutinds, nnz);

      if( cutrank != NULL )
         *cutrank = aggrrow->rank + 1;
   }

  TERMINATE:

   /* if we aborted early the tmpcoefs array needs to be cleaned */
   if( !(*success) )
   {
      SCIP_Real QUAD(tmp);
      QUAD_ASSIGN(tmp, 0.0);

      for( k = 0; k < nnz; ++k )
      {
         QUAD_ARRAY_STORE(tmpcoefs, tmpinds[k], tmp);
      }
   }
#ifndef NDEBUG
   for( k = 0; k < QUAD_ARRAY_SIZE(nvars); ++k )
   {
      if(tmpcoefs[k] != 0.0)
      {
         SCIPdebugMsg(scip, "tmpcoefs have not been reset\n");
         SCIPABORT();
      }
   }
#endif

   /* free temporary memory */
   SCIPfreeCleanBufferArray(scip, &tmpcoefs);
   SCIPfreeBufferArray(scip, &tmpinds);
   SCIPfreeBufferArray(scip, &coverpos);
   SCIPfreeBufferArray(scip, &covervals);
   SCIPfreeBufferArray(scip, &coverstatus);
   SCIPfreeBufferArray(scip, &boundtype);
   SCIPfreeBufferArray(scip, &varsign);

   return SCIP_OKAY;
}


/* =========================================== strongcg =========================================== */

/** Transform equation \f$ a \cdot x = b; lb \leq x \leq ub \f$ into standard form
 *    \f$ a^\prime \cdot x^\prime = b,\; 0 \leq x^\prime \leq ub' \f$.
 *
 *  Differs from cutsTransformMIR for continuous variables for which the lower bound must be used
 *  when in case their coefficient is positive and the upper bound in case their coefficient is
 *  negative. This forces all continuous variable to have a positive coefficient in the transformed
 *  row.
 *
 *  Transform variables (lb or ub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - lb_j,&   x_j = x^\prime_j + lb_j,&   a^\prime_j =  a_j,&   \mbox{if lb is used in transformation}\\
 *    x^\prime_j := ub_j - x_j,&   x_j = ub_j - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if ub is used in transformation}
 *  \end{array}
 *  \f]
 *  and move the constant terms \f$ a_j\, lb_j \f$ or \f$ a_j\, ub_j \f$ to the rhs.
 *
 *  Transform variables (vlb or vub):
 *  \f[
 *  \begin{array}{llll}
 *    x^\prime_j := x_j - (bl_j\, zl_j + dl_j),&   x_j = x^\prime_j + (bl_j\, zl_j + dl_j),&   a^\prime_j =  a_j,&   \mbox{if vlb is used in transf.} \\
 *    x^\prime_j := (bu_j\, zu_j + du_j) - x_j,&   x_j = (bu_j\, zu_j + du_j) - x^\prime_j,&   a^\prime_j = -a_j,&   \mbox{if vub is used in transf.}
 *  \end{array}
 *  \f]
 *  move the constant terms \f$ a_j\, dl_j \f$ or \f$ a_j\, du_j \f$ to the rhs, and update the coefficient of the VLB variable:
 *  \f[
 *  \begin{array}{ll}
 *    a_{zl_j} := a_{zl_j} + a_j\, bl_j,& \mbox{or} \\
 *    a_{zu_j} := a_{zu_j} + a_j\, bu_j &
 *  \end{array}
 *  \f]
 */
static
SCIP_RETCODE cutsTransformStrongCG(
   SCIP*                 scip,               /**< SCIP datastructure */
   MIR_DATA*             data,               /**< the MIR data structure for this cut */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable:
                                              *   vlb/vub_idx, or -1 for global lb/ub, or -2 for local lb/ub */
   SCIP_Bool*            freevariable,       /**< stores whether a free variable was found in MIR row -> invalid summation */
   SCIP_Bool*            localbdsused        /**< pointer to store whether local bounds were used in transformation */
   )
{
   SCIP_Real* bestlbs;
   SCIP_Real* bestubs;
   int* bestlbtypes;
   int* bestubtypes;
   SCIP_BOUNDTYPE* selectedbounds;
   int totalnnz;
   int s;
   int i;

   assert(data != NULL);
   assert(varsign != NULL);
   assert(boundtype != NULL);
   assert(freevariable != NULL);
   assert(localbdsused != NULL);

   totalnnz = data->totalnnz;

   *freevariable = FALSE;
   *localbdsused = FALSE;

   /* allocate temporary memory to store best bounds and bound types */
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbs, 2*totalnnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubs, 2*totalnnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestlbtypes, 2*totalnnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bestubtypes, 2*totalnnz) );
   SCIP_CALL( SCIPallocBufferArray(scip, &selectedbounds, 2*totalnnz) );

   /* transform the cut, one variable section at a time */
   for( s = 0; s < NSECTIONS; ++s )
   {
      int* indices = data->secindices[s];
      int cutindsstart = data->ncutinds;
      int usevbds = data->usevbds[s];

      i = 0;
      /* Iterate over all nonzeros in the section */
      while( i < data->secnnz[s] )
      {
         SCIP_Real QUAD(coef);
         int v = indices[i];

         /* due to variable bound usage, cancellation may have occurred */
         QUAD_ARRAY_LOAD(coef, data->cutcoefs, v);
         if( EPSZ(QUAD_TO_DBL(coef), QUAD_EPSILON) )
         {
            QUAD_ASSIGN(coef, 0.0);
            QUAD_ARRAY_STORE(data->cutcoefs, v, coef);
            --data->secnnz[s];
            --data->totalnnz;
            indices[i] = indices[data->secnnz[s]];
            /* do not increase the index */
            continue;
         }

         int cutindex = data->ncutinds;
         /* For continuous variables, we must choose the bound substitution so that they become positive in the cut */
         if( !data->isenfint[s] && !data->isimplint[s] )
         {
            if( QUAD_TO_DBL(coef) > 0.0 )
            {
               SCIP_Real simplelb;

               /* find closest lower bound in standard lower bound or variable lower bound for continuous variable so that it will have a positive coefficient */
               SCIP_CALL( findMIRBestLb(scip, data->vars[v], sol, data, usevbds, allowlocal,
                     bestlbs + cutindex, &simplelb, bestlbtypes + cutindex) );

               /* cannot create transformation for strongcg cut */
               if( SCIPisInfinity(scip, -bestlbs[cutindex]) )
               {
                  *freevariable = TRUE;
                  goto TERMINATE;
               }

               varsign[cutindex] = +1;
               selectedbounds[cutindex] = SCIP_BOUNDTYPE_LOWER;
            }
            else
            {
               SCIP_Real simpleub;

               assert(QUAD_TO_DBL(coef) < 0.0);

               /* find closest upper bound in standard upper bound or variable upper bound for continuous variable so that it will have a positive coefficient */
               SCIP_CALL( findMIRBestUb(scip, data->vars[v], sol, data, usevbds, allowlocal,
                     bestubs + cutindex, &simpleub, bestubtypes + cutindex) );

               /* cannot create transformation for strongcg cut */
               if( SCIPisInfinity(scip, bestubs[cutindex]) )
               {
                  *freevariable = TRUE;
                  goto TERMINATE;
               }

               varsign[cutindex] = -1;
               selectedbounds[cutindex] = SCIP_BOUNDTYPE_UPPER;
            }
         }
         else if( data->isimplint[s] )
         {
            /* For implied integers, we still prefer to choose the bound substitution that makes them positive, but
             * if we cannot manage to do so it is not an error, because we can still treat them as integer variables */
            SCIP_Real simplelb;
            SCIP_Real simpleub;
            SCIP_Bool lowerinf;
            SCIP_Bool upperinf;
            SCIP_Bool positive;

            /* find closest lower bound in standard lower bound or variable lower bound for continuous variable so that it will have a positive coefficient */
            SCIP_CALL( findMIRBestLb(scip, data->vars[v], sol, data, usevbds, allowlocal,
                  bestlbs + cutindex, &simplelb, bestlbtypes + cutindex) );

            /* find closest upper bound in standard upper bound or variable upper bound for continuous variable so that it will have a positive coefficient */
            SCIP_CALL( findMIRBestUb(scip, data->vars[v], sol, data, usevbds, allowlocal,
                  bestubs + cutindex, &simpleub, bestubtypes + cutindex) );

            lowerinf = SCIPisInfinity(scip, -bestlbs[cutindex]);
            upperinf = SCIPisInfinity(scip, bestubs[cutindex]);
            positive = QUAD_TO_DBL(coef) > 0.0;

            if( lowerinf && upperinf )
            {
               /* we found a free variable in the row with non-zero coefficient
                *  -> MIR row can't be transformed in standard form
                */
               *freevariable = TRUE;
               goto TERMINATE;
            }

            /* preferably, choose bound that makes value positive */
            if( (positive && lowerinf) || (!positive && !upperinf) )
            {
               varsign[cutindex] = -1;
               selectedbounds[cutindex] = SCIP_BOUNDTYPE_UPPER;
            }
            else
            {
               varsign[cutindex] = +1;
               selectedbounds[cutindex] = SCIP_BOUNDTYPE_LOWER;
            }
         }
         else
         {
            /* For explicit integers, we have no restrictions. */
            SCIP_CALL( determineBestBounds(scip, data->vars[v], sol, data, boundswitch, usevbds, allowlocal, FALSE, FALSE,
                  NULL, NULL, bestlbs + cutindex, bestubs + cutindex,
                  bestlbtypes + cutindex, bestubtypes + cutindex, selectedbounds + cutindex, freevariable) );

            if( *freevariable)
               goto TERMINATE;
         }

         data->cutinds[cutindex] = v;
         ++data->ncutinds;

         ++i;
      }

      /* perform bound substitution for all nonzeros in the section */
      for( i = cutindsstart; i < data->ncutinds; ++i )
      {
         SCIP_Real bestbnd;
         int v = data->cutinds[i];

         if( selectedbounds[i] == SCIP_BOUNDTYPE_LOWER )
         {
            assert(!SCIPisInfinity(scip, -bestlbs[i]));

            /* use lower bound as transformation bound: x'_j := x_j - lb_j */
            boundtype[i] = bestlbtypes[i];
            varsign[i] = +1;
            bestbnd = bestlbs[i];
         }
         else
         {
            assert(!SCIPisInfinity(scip, bestubs[i]));

            /* use upper bound as transformation bound: x'_j := ub_j - x_j */
            boundtype[i] = bestubtypes[i];
            varsign[i] = -1;
            bestbnd = bestubs[i];
         }

         doMIRBoundSubstitution(scip, data, varsign[i], boundtype[i], bestbnd, v, localbdsused);
      }
   }

   /* relax rhs to zero if it is close to */
   if( QUAD_TO_DBL(data->cutrhs) < 0.0 && QUAD_TO_DBL(data->cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(data->cutrhs, 0.0);

   TERMINATE:

   /* If we terminate early, we need to make sure all the zeros in the cut coefficient array are cancelled */
   if( *freevariable )
   {
      int j;
      int k;

      data->ncutinds = 0;
      for( j = 0; j < NSECTIONS; ++j )
      {
         int* indexlist = data->secindices[j];
         for( k = 0; k < data->secnnz[j]; ++k )
         {
            data->cutinds[data->ncutinds] = indexlist[k];
            ++data->ncutinds;
         }
      }
   }

   /*free temporary memory */
   SCIPfreeBufferArray(scip, &selectedbounds);
   SCIPfreeBufferArray(scip, &bestubtypes);
   SCIPfreeBufferArray(scip, &bestlbtypes);
   SCIPfreeBufferArray(scip, &bestubs);
   SCIPfreeBufferArray(scip, &bestlbs);

   return SCIP_OKAY;
}

/** Calculate fractionalities \f$ f_0 := b - down(b) \f$, \f$ f_j := a^\prime_j - down(a^\prime_j) \f$,
 *   integer \f$ k \geq 1 \f$ with \f$ 1/(k + 1) \leq f_0 < 1/k \f$  \f$ (\Rightarrow k = up(1/f_0) - 1) \f$ and
 *   integer \f$ 1 \leq p_j \leq k \f$ with \f$ f_0 + ((p_j - 1) \cdot (1 - f_0)/k) < f_j \leq f_0 + (p_j (1 - f_0)/k)\f$ \f$ (\Rightarrow p_j = up( k\,(f_j - f_0)/(1 - f_0) )) \f$
 * and derive strong CG cut \f$ \tilde{a} x^\prime \leq down(b) \f$
 * \f[
 * \begin{array}{rll}
 * integers : &  \tilde{a}_j = down(a^\prime_j)                &, if \qquad f_j \leq f_0 \\
 *            &  \tilde{a}_j = down(a^\prime_j) + p_j/(k + 1)  &, if \qquad f_j >  f_0 \\
 * continuous:&  \tilde{a}_j = 0                               &, if \qquad a^\prime_j \geq 0 \\
 *            &  \mbox{no strong CG cut found}                 &, if \qquad a^\prime_j <  0
 * \end{array}
 * \f]
 *
 * Transform inequality back to \f$ \hat{a}*x <= rhs \f$:
 *
 *  (lb or ub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - lb_j,&   x_j == x^\prime_j + lb_j,&   a^\prime_j ==  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{if lb was used in transformation} \\
 *    x^\prime_j := ub_j - x_j,&   x_j == ub_j - x^\prime_j,&   a^\prime_j == -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{if ub was used in transformation}
 * \end{array}
 * \f]
 * \f[
 *  and move the constant terms
 * \begin{array}{rl}
 *    -\tilde{a}_j * lb_j == -\hat{a}_j * lb_j, & \mbox{or} \\
 *     \tilde{a}_j * ub_j == -\hat{a}_j * ub_j &
 * \end{array}
 * \f]
 *  to the rhs.
 *
 *  (vlb or vub):
 * \f[
 * \begin{array}{lllll}
 *    x^\prime_j := x_j - (bl_j * zl_j + dl_j),&   x_j == x^\prime_j + (bl_j * zl_j + dl_j),&   a^\prime_j ==  a_j,&   \hat{a}_j :=  \tilde{a}_j,&   \mbox{(vlb)} \\
 *    x^\prime_j := (bu_j * zu_j + du_j) - x_j,&   x_j == (bu_j * zu_j + du_j) - x^\prime_j,&   a^\prime_j == -a_j,&   \hat{a}_j := -\tilde{a}_j,&   \mbox{(vub)}
 * \end{array}
 * \f]
 *  move the constant terms
 * \f[
 * \begin{array}{rl}
 *    -\tilde{a}_j * dl_j == -\hat{a}_j * dl_j,& \mbox{or} \\
 *     \tilde{a}_j * du_j == -\hat{a}_j * du_j &
 * \end{array}
 * \f]
 * to the rhs, and update the VB variable coefficients:
 * \f[
 * \begin{array}{ll}
 *    \hat{a}_{zl_j} := \hat{a}_{zl_j} - \tilde{a}_j * bl_j == \hat{a}_{zl_j} - \hat{a}_j * bl_j,& \mbox{or} \\
 *    \hat{a}_{zu_j} := \hat{a}_{zu_j} + \tilde{a}_j * bu_j == \hat{a}_{zu_j} - \hat{a}_j * bu_j &
 * \end{array}
 * \f]
 */
static
SCIP_RETCODE cutsRoundStrongCG(
   SCIP*                 scip,               /**< SCIP datastructure */
   MIR_DATA*             data,               /**< the MIR data structure for this cut */
   int*                  varsign,            /**< stores the sign of the transformed variable in summation */
   int*                  boundtype,          /**< stores the bound used for transformed variable (vlb/vub_idx or -1 for lb/ub)*/
   QUAD(SCIP_Real        f0),                /**< fractional value of rhs */
   SCIP_Real             k                   /**< factor to strengthen strongcg cut */
   )
{
   SCIP_Real QUAD(tmp);
   SCIP_Real QUAD(onedivoneminusf0);
   int cutindex;
   int s;
   int i;

   assert(data != NULL);
   assert(boundtype != NULL);
   assert(varsign != NULL);
   assert(0.0 < QUAD_TO_DBL(f0) && QUAD_TO_DBL(f0) < 1.0);

   SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0);
   SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0);

   /* Loop backwards through the sections, so that the reversing of varbound substitutions does not prematurely effect
    * the coefficients of variables in other sections, because the section index of a variable bound must always be
    * higher than that of the bounded variable. */
   cutindex = data->ncutinds - 1;
   for( s = NSECTIONS - 1; s >= 0; --s )
   {
      int* indices = data->secindices[s];
      int nnz = data->secnnz[s];
      SCIP_Bool enfintegral = data->isenfint[s];
      SCIP_Bool implintegral = data->isimplint[s];

      /* iterate backwards over indices in section, so we can easily shrink the section if we find zeros */
      for( i = nnz - 1; i >= 0 ; --i )
      {
         SCIP_Real QUAD(cutaj);
         SCIP_Real QUAD(aj);
         SCIP_VAR* var;
         int v = indices[i];
         int sign;
         int type;

         assert(0 <= v && v < data->nvars);
         assert(data->cutinds[cutindex] == v);
         sign = varsign[cutindex];
         assert(sign == +1 || sign == -1);
         type = boundtype[cutindex];

         --cutindex;

         var = data->vars[v];
         assert(var != NULL);
         assert(SCIPvarGetProbindex(var) == v );

         QUAD_ARRAY_LOAD(aj, data->cutcoefs, v);

         if( enfintegral || implintegral )
         {
            /* Variable is integral */
            SCIP_Real QUAD(downaj);
            SCIP_Real QUAD(fj);

            /* calculate the coefficient in the retransformed cut */
            QUAD_ARRAY_LOAD(aj, data->cutcoefs, v);
            QUAD_SCALE(aj, sign);
            SCIPquadprecEpsFloorQ(downaj, aj, SCIPepsilon(scip)); /*lint !e666*/
            SCIPquadprecSumQQ(fj, aj, -downaj);
            assert(QUAD_TO_DBL(fj) >= -SCIPepsilon(scip) && QUAD_TO_DBL(fj) < 1.0);

            if( SCIPisLE(scip, QUAD_TO_DBL(fj), QUAD_TO_DBL(f0)) )
               QUAD_ASSIGN_Q(cutaj, downaj); /* a_j */
            else
            {
               SCIP_Real pj;

               SCIPquadprecSumQQ(cutaj, fj, -f0);
               SCIPquadprecProdQD(cutaj, cutaj, k);
               SCIPquadprecProdQQ(cutaj, cutaj, onedivoneminusf0);
               pj = SCIPceil(scip, QUAD_TO_DBL(cutaj));
               assert(pj >= 0); /* should be >= 1, but due to rounding bias can be 0 if fj is almost equal to f0 */
               assert(pj <= k);
               SCIPquadprecDivDD(cutaj, pj, k + 1.0);
               SCIPquadprecSumQQ(cutaj, cutaj, downaj);
            }

            QUAD_SCALE(cutaj, sign);
         }
         else
         {
            /* Variable is continuous; must always be positive in strongcg cut. It will be automatically deleted. */
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS);
            assert(QUAD_TO_DBL(aj) * sign >= 0.0);
            QUAD_ASSIGN(cutaj, 0.0);
         }

         /* remove coefficient from cut if it becomes zero */
         if( EPSZ(QUAD_TO_DBL(cutaj), QUAD_EPSILON) )
         {
            QUAD_ASSIGN(cutaj, 0.0);
            QUAD_ARRAY_STORE(data->cutcoefs, v, cutaj);
            --data->totalnnz;
            --data->secnnz[s];
            indices[i] = indices[data->secnnz[s]];
            continue;
         }

         /* store the updated coefficient */
         QUAD_ARRAY_STORE(data->cutcoefs, v, cutaj);

         /* undo bound transformations. */
         if( type < 0 )
         {
            /* standard bound */
            /* move the constant term  -a~_j * lb_j == -a^_j * lb_j , or  a~_j * ub_j == -a^_j * ub_j  to the rhs */
            if( sign == +1 )
            {
               /* lower bound was used */
               if( type == -1 )
               {
                  assert(!SCIPisInfinity(scip, -SCIPvarGetLbGlobal(var)));
                  SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbGlobal(var));
                  SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, -SCIPvarGetLbLocal(var)));
                  SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetLbLocal(var));
                  SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);
               }
            }
            else
            {
               /* upper bound was used */
               if( type == -1 )
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbGlobal(var)));
                  SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbGlobal(var));
                  SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);
               }
               else
               {
                  assert(!SCIPisInfinity(scip, SCIPvarGetUbLocal(var)));
                  SCIPquadprecProdQD(tmp, cutaj, SCIPvarGetUbLocal(var));
                  SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);
               }
            }
         }
         else
         {
            /* variable bound */
            SCIP_VAR** vbz;
            SCIP_Real* vbb;
            SCIP_Real* vbd;
            SCIP_Real QUAD(zcoef);
            int vbidx;
            int zidx;

            /* variable bound */
            vbidx = type;

            /* change mirrhs and cutaj of integer variable z_j of variable bound */
            if( sign == +1 )
            {
               /* variable lower bound was used */
               assert(0 <= vbidx && vbidx < SCIPvarGetNVlbs(var));
               vbz = SCIPvarGetVlbVars(var);
               vbb = SCIPvarGetVlbCoefs(var);
               vbd = SCIPvarGetVlbConstants(var);
            }
            else
            {
               /* variable upper bound was used */
               assert(0 <= vbidx && vbidx < SCIPvarGetNVubs(var));
               vbz = SCIPvarGetVubVars(var);
               vbb = SCIPvarGetVubCoefs(var);
               vbd = SCIPvarGetVubConstants(var);
            }
            assert(SCIPvarIsActive(vbz[vbidx]));
            zidx = SCIPvarGetProbindex(vbz[vbidx]);
            assert(varSection(data, zidx) > s);

            SCIPquadprecProdQD(tmp, cutaj, vbd[vbidx]);
            SCIPquadprecSumQQ(data->cutrhs, data->cutrhs, tmp);

            SCIPquadprecProdQD(tmp, cutaj, vbb[vbidx]);
            QUAD_ARRAY_LOAD(zcoef, data->cutcoefs, zidx);

            /* update sparsity pattern */
            if( QUAD_HI(zcoef) == 0.0 )
            {
               int zsection = varSection(data, zidx);
               data->secindices[zsection][data->secnnz[zsection]] = zidx;
               ++data->secnnz[zsection];
               ++data->totalnnz;
            }

            SCIPquadprecSumQQ(zcoef, zcoef, -tmp);
            QUAD_HI(zcoef) = NONZERO(QUAD_HI(zcoef));
            QUAD_ARRAY_STORE(data->cutcoefs, zidx, zcoef);
            assert(QUAD_HI(zcoef) != 0.0);
         }
      }
   }

   /* Finally, store the relevant data in cutinds which is the array used by the other functions */
   data->ncutinds = 0;
   for( s = 0; s < NSECTIONS; ++s )
   {
      int* indices = data->secindices[s];
      int nnz = data->secnnz[s];
      for( i = 0; i < nnz; ++i )
      {
         data->cutinds[data->ncutinds] = indices[i];
         ++data->ncutinds;
      }
   }

   return SCIP_OKAY;
}

/** substitute aggregated slack variables:
 *
 *  The coefficient of the slack variable \f$s_r\f$ is equal to the row's weight times the slack's sign, because the slack
 *  variable only appears in its own row: \f$ a^\prime_r = scale \cdot weight[r] \cdot slacksign[r] \f$.
 *
 *  Depending on the slack's type (integral or continuous), its coefficient in the cut calculates as follows:
 *  \f[
 *  \begin{array}{rll}
 *    integers:  & \hat{a}_r = \tilde{a}_r = down(a^\prime_r),                  & if \qquad f_r \leq f_0 \\
 *               & \hat{a}_r = \tilde{a}_r = down(a^\prime_r) + p_r/(k + 1),    & if \qquad f_r >  f_0 \\
 *    continuous:& \hat{a}_r = \tilde{a}_r = 0,                                 & if \qquad a^\prime_r \geq 0 \\
 *               & \mbox{no strong CG cut found},                               & if \qquad a^\prime_r <  0
 *  \end{array}
 *  \f]
 *
 *  Substitute \f$ \hat{a}_r \cdot s_r \f$ by adding \f$ \hat{a}_r \f$ times the slack's definition to the cut.
 */
static
SCIP_RETCODE cutsSubstituteStrongCG(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  slacksign,          /**< stores the sign of the row's slack variable in summation */
   int*                  rowinds,            /**< sparsity pattern of used rows */
   int                   nrowinds,           /**< number of used rows */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            cutcoefs,           /**< array of coefficients of cut */
   QUAD(SCIP_Real*       cutrhs),            /**< pointer to right hand side of cut */
   int*                  cutinds,            /**< array of variables problem indices for non-zero coefficients in cut */
   int*                  nnz,                /**< number of non-zeros in cut */
   QUAD(SCIP_Real        f0),                /**< fractional value of rhs */
   SCIP_Real             k                   /**< factor to strengthen strongcg cut */
   )
{  /*lint --e{715}*/
   SCIP_ROW** rows;
   SCIP_Real QUAD(onedivoneminusf0);
   int i;

   assert(scip != NULL);
   assert(weights != NULL);
   assert(slacksign != NULL);
   assert(rowinds != NULL);
   assert(SCIPisPositive(scip, scale));
   assert(cutcoefs != NULL);
   assert(QUAD_HI(cutrhs) != NULL);
   assert(cutinds != NULL);
   assert(nnz != NULL);
   assert(0.0 < QUAD_TO_DBL(f0) && QUAD_TO_DBL(f0) < 1.0);

   SCIPquadprecSumQD(onedivoneminusf0, -f0, 1.0);
   SCIPquadprecDivDQ(onedivoneminusf0, 1.0, onedivoneminusf0);

   rows = SCIPgetLPRows(scip);
   for( i = 0; i < nrowinds; i++ )
   {
      SCIP_ROW* row;
      SCIP_Real QUAD(ar);
      SCIP_Real QUAD(downar);
      SCIP_Real QUAD(cutar);
      SCIP_Real QUAD(fr);
      SCIP_Real mul;
      int r;

      r = rowinds[i];
      assert(0 <= r && r < SCIPgetNLPRows(scip));
      assert(slacksign[i] == -1 || slacksign[i] == +1);
      assert(!SCIPisZero(scip, weights[i]));

      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->cols_index != NULL);
      assert(row->len == 0 || row->vals != NULL);

      /* get the slack's coefficient a'_r in the aggregated row */
      SCIPquadprecProdDD(ar, slacksign[i] * scale, weights[i]);

      /* calculate slack variable's coefficient a_r in the cut */
      if( row->integral )
      {
         /* slack variable is always integral */
         SCIPquadprecEpsFloorQ(downar, ar, SCIPepsilon(scip)); /*lint !e666*/
         SCIPquadprecSumQQ(fr, ar, -downar);
         assert(QUAD_TO_DBL(fr) >= -SCIPepsilon(scip) && QUAD_TO_DBL(fr) < 1.0);

         if( SCIPisLE(scip, QUAD_TO_DBL(fr), QUAD_TO_DBL(f0)) )
            QUAD_ASSIGN_Q(cutar, downar); /* a_r */
         else
         {
            SCIP_Real pr;

            SCIPquadprecSumQQ(cutar, fr, -f0);
            SCIPquadprecProdQD(cutar, cutar, k);
            SCIPquadprecProdQQ(cutar, cutar, onedivoneminusf0);
            pr = SCIPceil(scip, QUAD_TO_DBL(cutar));
            assert(pr >= 0); /* should be >= 1, but due to rounding bias can be 0 if fr is almost equal to f0 */
            assert(pr <= k);
            SCIPquadprecDivDD(cutar, pr, k + 1.0);
            SCIPquadprecSumQQ(cutar, cutar, downar);
         }
      }
      else
      {
         /* slack variable is continuous: */
         assert(QUAD_TO_DBL(ar) >= 0.0);
         continue; /* slack can be ignored, because its coefficient is reduced to 0.0 */
      }

      /* if the coefficient was reduced to zero, ignore the slack variable */
      if( EPSZ(QUAD_TO_DBL(cutar), QUAD_EPSILON) )
         continue;

      /* depending on the slack's sign, we have
       *   a*x + c + s == rhs  =>  s == - a*x - c + rhs,  or  a*x + c - s == lhs  =>  s == a*x + c - lhs
       * substitute a_r * s_r by adding a_r times the slack's definition to the cut.
       */
      mul = -slacksign[i] * QUAD_TO_DBL(cutar);

      /* add the slack's definition multiplied with a_j to the cut */
      SCIP_CALL( varVecAddScaledRowCoefsQuad(cutinds, cutcoefs, nnz, row, mul) );

      /* move slack's constant to the right hand side */
      if( slacksign[i] == +1 )
      {
         SCIP_Real rhs;

         /* a*x + c + s == rhs  =>  s == - a*x - c + rhs: move a_r * (rhs - c) to the right hand side */
         assert(!SCIPisInfinity(scip, row->rhs));
         rhs = row->rhs - row->constant;
         if( row->integral )
         {
            /* the right hand side was implicitly rounded down in row aggregation */
            rhs = SCIPfloor(scip, rhs);
         }

         SCIPquadprecProdQD(cutar, cutar, rhs);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, -cutar);
      }
      else
      {
         SCIP_Real lhs;

         /* a*x + c - s == lhs  =>  s == a*x + c - lhs: move a_r * (c - lhs) to the right hand side */
         assert(!SCIPisInfinity(scip, -row->lhs));
         lhs = row->lhs - row->constant;
         if( row->integral )
         {
            /* the left hand side was implicitly rounded up in row aggregation */
            lhs = SCIPceil(scip, lhs);
         }

         SCIPquadprecProdQD(cutar, cutar, lhs);
         SCIPquadprecSumQQ(*cutrhs, *cutrhs, cutar);
      }
   }

   /* relax rhs to zero, if it's very close to 0 */
   if( QUAD_TO_DBL(*cutrhs) < 0.0 && QUAD_TO_DBL(*cutrhs) >= -SCIPepsilon(scip) )
      QUAD_ASSIGN(*cutrhs, 0.0);

   return SCIP_OKAY;
}


/** calculates a strong CG cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in a strongcg cut
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcalcStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Bool             postprocess,        /**< apply a post-processing step to the resulting cut? */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   int                   vartypeusevbds,     /**< for all variable types with index smaller than this number, variable
                                              *   type substitution is allowed. The indices are: 0: continuous,
                                              *   1: continuous implint., 2: integer implint, 3: binary implint,
                                              *   4: integer, 5: binary */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute a strong CG cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store the efficacy of the cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid cut was returned */
   )
{
   int i;
   int nvars;
   int* varsign;
   int* boundtype;
   SCIP_Real QUAD(downrhs);
   SCIP_Real QUAD(f0);
   SCIP_Real QUAD(tmp);
   SCIP_Real large;
   SCIP_Real k;
   SCIP_Bool freevariable;
   SCIP_Bool localbdsused;
   MIR_DATA* data;

   assert(scip != NULL);
   assert(aggrrow != NULL);
   assert(SCIPisPositive(scip, scale));
   assert(cutcoefs != NULL);
   assert(cutrhs != NULL);
   assert(cutinds != NULL);
   assert(success != NULL);
   assert(cutislocal != NULL);

   SCIPdebugMsg(scip, "calculating strong CG cut (scale: %g)\n", scale);

   *success = FALSE;

   /* determine value from which fractionalities are no longer reliable within tolerance */
   large = SCIPgetHugeValue(scip) * SCIPepsilon(scip);

   /* terminate if an integral slack fractionality is unreliable or a negative continuous slack variable is present */
   for( i = 0; i < aggrrow->nrows; ++i )
   {
      if( ( scip->lp->rows[aggrrow->rowsinds[i]]->integral && ABS(aggrrow->rowweights[i] * scale) > large )
         || ( !scip->lp->rows[aggrrow->rowsinds[i]]->integral && aggrrow->rowweights[i] * aggrrow->slacksign[i] < 0.0 ) )
         return SCIP_OKAY;
   }

   /* allocate temporary memory */
   nvars = SCIPgetNVars(scip);
   SCIP_CALL( SCIPallocBufferArray(scip, &varsign, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &boundtype, nvars) );

   /* Initialize cut data */
   int l;
   int nnz;

   assert(vartypeusevbds >= 0 && vartypeusevbds < NSECTIONS);

   SCIP_CALL(SCIPallocBuffer(scip, &data));

   nnz = aggrrow->nnz;
   data->totalnnz = nnz;

   /* initialize sections */
   for( l = 0; l < NSECTIONS; ++l )
   {
      SCIP_CALL(SCIPallocBufferArray(scip, &data->secindices[l], nnz));
      data->secnnz[l] = 0;
      /* Cont. | cont impl. | int impl. | bin impl. | int | bin */
      assert(NSECTIONS == 6); /* If the section definition is changed, the below lines should also be adjusted to match */
      data->isenfint[l] = l >= 2 ? TRUE : FALSE;
      data->isimplint[l] = l >= 1 && l <= 3 ? TRUE : FALSE;
      /* Use variable bounds for the sections specified by the user */
      data->usevbds[l] = l < vartypeusevbds ? 2 : 0;
   }

   /* Problem data needs to be initialized before cut data as it is used to partition the variables into the sections */
   data->vars = SCIPgetVars(scip);
   data->nvars = SCIPgetNVars(scip);
   data->nbinvars = SCIPgetNBinVars(scip);
   data->nintvars = SCIPgetNIntVars(scip);
   data->nbinimplvars = SCIPgetNBinImplVars(scip);
   data->nintimplvars = SCIPgetNIntImplVars(scip);
   data->ncontimplvars = SCIPgetNContImplVars(scip);
   data->ncontvars = SCIPgetNContVars(scip);

   SCIP_CALL(SCIPallocCleanBufferArray(scip, &( data->cutcoefs ), QUAD_ARRAY_SIZE(data->nvars)));
   SCIP_CALL(SCIPallocBufferArray(scip, &data->cutinds, data->nvars));

   SCIPquadprecProdQD(data->cutrhs, aggrrow->rhs, scale);

   if( nnz > 0 )
   {
      /* Initalize cut with the aggregation */
      BMScopyMemoryArray(data->cutinds, aggrrow->inds, nnz);

      for( l = 0; l < nnz; ++l )
      {
         SCIP_Real QUAD(coef);
         int m = aggrrow->inds[l];

         QUAD_ARRAY_LOAD(coef, aggrrow->vals, m);

         SCIPquadprecProdQD(coef, coef, scale);

         QUAD_ARRAY_STORE(data->cutcoefs, m, coef);

         assert(QUAD_HI(coef) != 0.0);
      }

      /* Sort the array by problem index and add the variables to their sections */
      SCIPsortDownInt(data->cutinds, nnz);
      for( l = 0; l < nnz; ++l )
      {
         int section = varSection(data, data->cutinds[l]);
         data->secindices[section][data->secnnz[section]] = data->cutinds[l];
         ++data->secnnz[section];
      }
   }

   data->ncutinds = 0;
   *cutislocal = aggrrow->local;

   if( data->totalnnz > 0 )
   {
      int firstcontvar;

      /* Transform equation  a*x == b, lb <= x <= ub  into standard form
       *   a'*x' == b, 0 <= x' <= ub'.
       *
       * Transform variables (lb or ub):
       *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   if lb is used in transformation
       *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   if ub is used in transformation
       * and move the constant terms "a_j * lb_j" or "a_j * ub_j" to the rhs.
       *
       * Transform variables (vlb or vub):
       *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   if vlb is used in transf.
       *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   if vub is used in transf.
       * move the constant terms "a_j * dl_j" or "a_j * du_j" to the rhs, and update the coefficient of the VLB variable:
       *   a_{zl_j} := a_{zl_j} + a_j * bl_j, or
       *   a_{zu_j} := a_{zu_j} + a_j * bu_j
       */
      SCIP_CALL( cutsTransformStrongCG(scip, data, sol, boundswitch, allowlocal, varsign, boundtype, &freevariable, &localbdsused) );

      if( freevariable )
         goto TERMINATE;

      assert(allowlocal || !localbdsused);
      *cutislocal = *cutislocal || localbdsused;

      firstcontvar = nvars - SCIPgetNContVars(scip);

      /* terminate if an integral coefficient fractionality is unreliable */
      for( i = data->ncutinds - 1; i >= 0 && data->cutinds[i] < firstcontvar; --i )
      {
         SCIP_Real QUAD(coef);

         QUAD_ARRAY_LOAD(coef, data->cutcoefs, data->cutinds[i]);

         if( ABS(QUAD_TO_DBL(coef)) > large )
            goto TERMINATE;
      }

      SCIPdebug(printCutQuad(scip, NULL, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, data->ncutinds, FALSE, FALSE));
   }

   /* terminate if the side fractionality is unreliable */
   if( ABS(QUAD_TO_DBL(data->cutrhs)) > large )
      goto TERMINATE;

   /* Calculate
    *  - fractionalities  f_0 := b - down(b), f_j := a'_j - down(a'_j)
    *  - integer k >= 1 with 1/(k + 1) <= f_0 < 1/k
    *    (=> k = up(1/f_0) - 1)
    *  - integer 1 <= p_j <= k with f_0 + ((p_j - 1) * (1 - f_0)/k) < f_j <= f_0 + (p_j * (1 - f_0)/k)
    *    (=> p_j = up( (f_j - f_0)/((1 - f_0)/k) ))
    * and derive strong CG cut
    *   a~*x' <= (k+1) * down(b)
    * integers :  a~_j = down(a'_j)                , if f_j <= f_0
    *             a~_j = down(a'_j) + p_j/(k + 1)  , if f_j >  f_0
    * continuous: a~_j = 0                         , if a'_j >= 0
    *             no strong CG cut found           , if a'_j <  0
    *
    * Transform inequality back to a^*x <= rhs:
    *
    * (lb or ub):
    *   x'_j := x_j - lb_j,   x_j == x'_j + lb_j,   a'_j ==  a_j,   a^_j :=  a~_j,   if lb was used in transformation
    *   x'_j := ub_j - x_j,   x_j == ub_j - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   if ub was used in transformation
    * and move the constant terms
    *   -a~_j * lb_j == -a^_j * lb_j, or
    *    a~_j * ub_j == -a^_j * ub_j
    * to the rhs.
    *
    * (vlb or vub):
    *   x'_j := x_j - (bl_j * zl_j + dl_j),   x_j == x'_j + (bl_j * zl_j + dl_j),   a'_j ==  a_j,   a^_j :=  a~_j,   (vlb)
    *   x'_j := (bu_j * zu_j + du_j) - x_j,   x_j == (bu_j * zu_j + du_j) - x'_j,   a'_j == -a_j,   a^_j := -a~_j,   (vub)
    * move the constant terms
    *   -a~_j * dl_j == -a^_j * dl_j, or
    *    a~_j * du_j == -a^_j * du_j
    * to the rhs, and update the VB variable coefficients:
    *   a^_{zl_j} := a^_{zl_j} - a~_j * bl_j == a^_{zl_j} - a^_j * bl_j, or
    *   a^_{zu_j} := a^_{zu_j} + a~_j * bu_j == a^_{zu_j} - a^_j * bu_j
    */
   SCIPquadprecEpsFloorQ(downrhs, data->cutrhs, SCIPepsilon(scip)); /*lint !e666*/
   SCIPquadprecSumQQ(f0, data->cutrhs, -downrhs);
   assert(QUAD_TO_DBL(f0) >= -SCIPepsilon(scip) && QUAD_TO_DBL(f0) < 1.0);

   if( QUAD_TO_DBL(f0) < minfrac || QUAD_TO_DBL(f0) > maxfrac )
      goto TERMINATE;

   /* renormalize the f0 value */
   SCIPquadprecSumDD(f0, QUAD_HI(f0), QUAD_LO(f0));

   SCIPquadprecDivDQ(tmp, 1.0, f0);
   SCIPquadprecSumQD(tmp, tmp, -1.0);
   k = SCIPceil(scip, QUAD_TO_DBL(tmp));
   QUAD_ASSIGN_Q(data->cutrhs, downrhs);

   if( data->totalnnz > 0 )
   {
      SCIP_CALL( cutsRoundStrongCG(scip, data, varsign, boundtype, QUAD(f0), k) );
      SCIPdebug(printCutQuad(scip, sol, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, data->ncutinds, FALSE, FALSE));
   }

   /* substitute aggregated slack variables:
    *
    * The coefficient of the slack variable s_r is equal to the row's weight times the slack's sign, because the slack
    * variable only appears in its own row:
    *    a'_r = scale * weight[r] * slacksign[r].
    *
    * Depending on the slacks type (integral or continuous), its coefficient in the cut calculates as follows:
    *   integers :  a_r = a~_r = (k + 1) * down(a'_r)        , if f_r <= f0
    *               a_r = a~_r = (k + 1) * down(a'_r) + p_r  , if f_r >  f0
    *   continuous: a_r = a~_r = 0                           , if a'_r >= 0
    *               a_r = a~_r = a'_r/(1 - f0)               , if a'_r <  0
    *
    * Substitute a_r * s_r by adding a_r times the slack's definition to the cut.
    */
   SCIP_CALL( cutsSubstituteStrongCG(scip, aggrrow->rowweights, aggrrow->slacksign, aggrrow->rowsinds,
         aggrrow->nrows, scale, data->cutcoefs, QUAD(&data->cutrhs), data->cutinds, &data->ncutinds, QUAD(f0), k) );
   SCIPdebug(printCutQuad(scip, sol, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, data->ncutinds, FALSE, FALSE));

   /* remove all nearly-zero coefficients from strong CG row and relax the right hand side correspondingly in order to
    * prevent numerical rounding errors
    */
   if( postprocess )
   {
      SCIP_CALL( postprocessCutQuad(scip, *cutislocal, data->cutinds, data->cutcoefs, &data->ncutinds, QUAD(&data->cutrhs), success) );
   }
   else
   {
      *success = ! removeZerosQuad(scip, SCIPsumepsilon(scip), *cutislocal, data->cutcoefs, QUAD(&data->cutrhs), data->cutinds, &data->ncutinds);
   }
   SCIPdebug(printCutQuad(scip, sol, data->cutcoefs, QUAD(data->cutrhs), data->cutinds, data->ncutinds, FALSE, FALSE));

   if( *success )
   {
      *cutrhs = QUAD_TO_DBL(data->cutrhs);
      *cutnnz = data->ncutinds;

      /* store cut in given array in sparse representation and clean buffer array */
      for( i = 0; i < *cutnnz; ++i )
      {
         SCIP_Real QUAD(coef);
         int j = data->cutinds[i];

         QUAD_ARRAY_LOAD(coef, data->cutcoefs, j);
         assert(QUAD_HI(coef) != 0.0);

         cutcoefs[i] = QUAD_TO_DBL(coef);
         cutinds[i] = j;
         QUAD_ASSIGN(coef, 0.0);
         QUAD_ARRAY_STORE(data->cutcoefs, j, coef);
      }

      if( cutefficacy != NULL )
         *cutefficacy = calcEfficacy(scip, sol, cutcoefs, *cutrhs, cutinds, *cutnnz);

      if( cutrank != NULL )
         *cutrank = aggrrow->rank + 1;
   }

  TERMINATE:

   /* if we aborted early the temporary coefficients need to be cleaned */
   if( !(*success) )
   {
      QUAD_ASSIGN(tmp, 0.0);

      for( i = 0; i < data->ncutinds; ++i )
      {
         QUAD_ARRAY_STORE(data->cutcoefs, data->cutinds[i], tmp);
      }
   }

   if( data->cutinds != NULL )
      SCIPfreeBufferArray(scip, &data->cutinds);

   if( data->cutcoefs != NULL )
      SCIPfreeCleanBufferArray(scip, &data->cutcoefs);

   for( int s = NSECTIONS - 1; s >= 0; --s )
   {
      SCIPfreeBufferArray(scip, &data->secindices[s]);
   }

   SCIPfreeBuffer(scip, &data);

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &boundtype);
   SCIPfreeBufferArray(scip, &varsign);

   return SCIP_OKAY;
}
