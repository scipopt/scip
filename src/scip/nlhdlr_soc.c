/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   nlhdlr_soc.c
 * @ingroup DEFPLUGINS_NLHDLR
 * @brief  nonlinear handler for second order cone constraints

 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Fabian Wegscheider
 *
 * This is a nonlinear handler for second order cone constraints of the form
 *
 * \f[\sqrt{\sum_{i=1}^{n} (v_i^T x + \beta_i)^2} \leq v_{n+1}^T x + \beta_{n+1}.\f]
 *
 * Note that \f$v_i\f$, for \f$i \leq n\f$, could be 0, thus allowing a positive constant term inside the root.
 *
 * @todo test if it makes sense to only disaggregate when nterms > some parameter
 *
 */

#include <string.h>

#include "scip/nlhdlr_soc.h"
#include "scip/cons_nonlinear.h"
#include "scip/expr_pow.h"
#include "scip/expr_sum.h"
#include "scip/expr_var.h"
#include "scip/debug.h"
#include "scip/pub_nlhdlr.h"
#include "scip/lapack_calls.h"


/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "soc"
#define NLHDLR_DESC               "nonlinear handler for second-order cone structures"
#define NLHDLR_DETECTPRIORITY       100 /**< priority of the nonlinear handler for detection */
#define NLHDLR_ENFOPRIORITY         100 /**< priority of the nonlinear handler for enforcement */
#define DEFAULT_MINCUTEFFICACY     1e-5 /**< default value for parameter mincutefficacy */
#define DEFAULT_COMPEIGENVALUES    TRUE /**< default value for parameter compeigenvalues */

/*
 * Data structures
 */

/** nonlinear handler expression data. The data is structured in the following way:
 *
 *  A 'term' is one of the arguments of the quadratic terms, i.e. \f$v_i^T x + beta_i\f$.
 *  The last term is always the one on the right-hand side. This means that nterms is
 *  equal to n+1 in the above description.
 *
 *  - vars contains a list of all expressions which are treated as variables (no duplicates)
 *  - offsets contains the constants beta_i of each term
 *  - transcoefs contains the non-zero values of the transformation vectors v_i of each term
 *  - transcoefsidx contains for each entry of transcoefs the position of the respective variable in vars
 *  - termbegins contains the index at which the transcoefs of each term start, with a sentinel value
 *  - nterms is the total number of terms appearing on both sides
 *  - nvars is the total number of unique variables appearing (length of vars)
 *
 *  Note that the numbers of nonzeroes in v_i is termbegins[i+1] - termbegins[i] and that
 *  the total number of entries in transcoefs and transcoefsidx is termbegins[nterms]
 *
 *  The disaggregation is implicitly stored in the variables disvars and disrow. An SOC as
 *  described above is replaced by n smaller SOCs
 *
 *              (v_i^T x + beta_i)^2 <= disvar_i     * (v_{n+1}^T x + beta_{n+1})
 *
 *  and the row       sum_i disvar_i <= v_{n+1}^T x + beta_{n+1}.
 *
 *  The disaggregation only happens if we have more than 3 terms.
 *
 *  Example: The constraint sqrt(5 + (3x - 4y + 2)^2 + y^2 + 7z^2) <= 5x - y - 1
 *           results in the following nlhdlrexprdata:
 *
 *           vars = {x, y, z}
 *           offsets = {2, 0, 0, sqrt(5), -1}
 *           transcoefs = {3, -4, 1, sqrt(7), 5, -1}
 *           transcoefsidx = {0, 1, 1, 2, 0, 1}
 *           termbegins = {0, 2, 3, 4, 4, 6}
 *           nvars = 3
 *           nterms = 5
 *
 * @note: due to the current implementation, the constant term is the second to last term, except when the SOC was a rotated
 * SOC, e.g., 1 + x^2 - y*z, i.e., when detected by detectSocQuadraticSimple. In that case, the constant is third to
 * last term.
 */
struct SCIP_NlhdlrExprData
{
   SCIP_EXPR**           vars;               /**< expressions which (aux)variables appear on both sides (x) */
   SCIP_Real*            offsets;            /**< offsets of both sides (beta_i) */
   SCIP_Real*            transcoefs;         /**< non-zeros of linear transformation vectors (v_i) */
   int*                  transcoefsidx;      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  termbegins;         /**< starting indices of transcoefs for each term */
   int                   nvars;              /**< total number of variables appearing */
   int                   nterms;             /**< number of summands in the SQRT +1 for RHS (n+1) */

   /* variables for cone disaggregation */
   SCIP_VAR**            disvars;            /**< disaggregation variables for each term in lhs */
   SCIP_ROW*             disrow;             /**< disaggregation row */

   /* separation data */
   SCIP_Real*            varvals;            /**< current values for vars */
   SCIP_Real*            disvarvals;         /**< current values for disvars */
};

struct SCIP_NlhdlrData
{
   SCIP_Real             mincutefficacy;     /**< minimum efficacy a cut need to be added */
   SCIP_Bool             compeigenvalues;    /**< whether Eigenvalue computations should be done to detect complex cases */
};

/*
 * Local methods
 */

#ifdef SCIP_DEBUG
/** prints the nlhdlr expression data */
static
void printNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< pointer to store nonlinear handler expression data */
   )
{
   int nterms;
   int i;
   int j;

   nterms = nlhdlrexprdata->nterms;

   SCIPinfoMessage(scip, NULL, "SQRT( ");

   for( i = 0; i < nterms - 1; ++i )
   {
      int startidx;

      startidx = nlhdlrexprdata->termbegins[i];

      /* v_i is 0 */
      if( startidx == nlhdlrexprdata->termbegins[i + 1] )
      {
         assert(nlhdlrexprdata->offsets[i] != 0.0);

         SCIPinfoMessage(scip, NULL, "%g", SQR(nlhdlrexprdata->offsets[i]));
         continue;
      }

      /* v_i is not 0 */
      SCIPinfoMessage(scip, NULL, "(");

      for( j = startidx; j < nlhdlrexprdata->termbegins[i + 1]; ++j )
      {
         if( nlhdlrexprdata->transcoefs[j] != 1.0 )
            SCIPinfoMessage(scip, NULL, " %+g*", nlhdlrexprdata->transcoefs[j]);
         else
            SCIPinfoMessage(scip, NULL, " +");
         if( SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]) != NULL )
         {
            SCIPinfoMessage(scip, NULL, "%s", SCIPvarGetName(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]])));
            SCIPinfoMessage(scip, NULL, "(%p)", (void*)nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]);
         }
         else
            SCIPinfoMessage(scip, NULL, "%p", (void*)nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]);
      }
      if( nlhdlrexprdata->offsets[i] != 0.0 )
         SCIPinfoMessage(scip, NULL, " %+g", nlhdlrexprdata->offsets[i]);

      SCIPinfoMessage(scip, NULL, ")^2");

      if( i < nterms - 2 )
         SCIPinfoMessage(scip, NULL, " + ");
   }

   SCIPinfoMessage(scip, NULL, " ) <=");

   for( j = nlhdlrexprdata->termbegins[nterms-1]; j < nlhdlrexprdata->termbegins[nterms]; ++j )
   {
      if( nlhdlrexprdata->transcoefs[j] != 1.0 )
         SCIPinfoMessage(scip, NULL, " %+g*", nlhdlrexprdata->transcoefs[j]);
      else
         SCIPinfoMessage(scip, NULL, " +");
      if( SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]) != NULL )
         SCIPinfoMessage(scip, NULL, "%s", SCIPvarGetName(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]])));
      else
         SCIPinfoMessage(scip, NULL, "%p", (void*)nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]);
   }
   if( nlhdlrexprdata->offsets[nterms-1] != 0.0 )
      SCIPinfoMessage(scip, NULL, " %+g", nlhdlrexprdata->offsets[nterms-1]);

   SCIPinfoMessage(scip, NULL, "\n");
}
#endif

/** helper method to create variables for the cone disaggregation */
static
SCIP_RETCODE createDisaggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< nonlinear handler expression data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int ndisvars;
   int i;

   assert(nlhdlrexprdata != NULL);

   ndisvars = nlhdlrexprdata->nterms - 1;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexprdata->disvars, ndisvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexprdata->disvarvals, ndisvars) );

   /* create disaggregation variables representing the epigraph of (v_i^T x + beta_i)^2 / (v_{n+1}^T x + beta_{n+1}) */
   for( i = 0; i < ndisvars; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_%d", (void*) expr, i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &nlhdlrexprdata->disvars[i], name, 0.0, SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS) );
      SCIPvarMarkRelaxationOnly(nlhdlrexprdata->disvars[i]);

      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[i]) );
      SCIP_CALL( SCIPaddVarLocksType(scip, nlhdlrexprdata->disvars[i], SCIP_LOCKTYPE_MODEL, 1, 1) );
   }

   return SCIP_OKAY;
}

/** helper method to free variables for the cone disaggregation */
static
SCIP_RETCODE freeDisaggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< nonlinear handler expression data */
   )
{
   int ndisvars;
   int i;

   assert(nlhdlrexprdata != NULL);

   if( nlhdlrexprdata->disvars == NULL )
      return SCIP_OKAY;

   ndisvars = nlhdlrexprdata->nterms - 1;

   /* release variables */
   for( i = 0; i < ndisvars; ++i )
   {
      SCIP_CALL( SCIPaddVarLocksType(scip, nlhdlrexprdata->disvars[i], SCIP_LOCKTYPE_MODEL, -1, -1) );
      SCIP_CALL( SCIPreleaseVar(scip, &nlhdlrexprdata->disvars[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->disvars, ndisvars);
   SCIPfreeBlockMemoryArrayNull(scip, &nlhdlrexprdata->disvarvals, ndisvars);

   return SCIP_OKAY;
}

/** helper method to create the disaggregation row \f$\text{disvars}_i \leq v_{n+1}^T x + \beta_{n+1}\f$ */
static
SCIP_RETCODE createDisaggrRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraint handler */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata      /**< nonlinear handler expression data */
   )
{
   SCIP_Real beta;
   char name[SCIP_MAXSTRLEN];
   int ndisvars;
   int nterms;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->disrow == NULL);

   nterms = nlhdlrexprdata->nterms;
   beta = nlhdlrexprdata->offsets[nterms - 1];

   ndisvars = nterms - 1;

   /* create row 0 <= beta_{n+1} */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_row", (void*) expr);
   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &nlhdlrexprdata->disrow, conshdlr, name,
         -SCIPinfinity(scip), beta, FALSE, FALSE, TRUE) );

   /* add disvars to row */
   for( i = 0; i < ndisvars; ++i )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, nlhdlrexprdata->disrow, nlhdlrexprdata->disvars[i], 1.0) );
   }

   /* add rhs vars to row */
   for( i = nlhdlrexprdata->termbegins[nterms - 1]; i < nlhdlrexprdata->termbegins[nterms]; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[i]]);
      assert(var != NULL);

      coef = -nlhdlrexprdata->transcoefs[i];

      SCIP_CALL( SCIPaddVarToRow(scip, nlhdlrexprdata->disrow, var, coef) );
   }

   return SCIP_OKAY;
}

/** helper method to create nonlinear handler expression data */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           vars,               /**< expressions which variables appear on both sides (\f$x\f$) */
   SCIP_Real*            offsets,            /**< offsets of bot sides (\f$beta_i\f$) */
   SCIP_Real*            transcoefs,         /**< non-zeroes of linear transformation vectors (\f$v_i\f$) */
   int*                  transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  termbegins,         /**< starting indices of transcoefs for each term */
   int                   nvars,              /**< total number of variables appearing */
   int                   nterms,             /**< number of summands in the SQRT, +1 for RHS */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata      /**< pointer to store nonlinear handler expression data */
   )
{
   int ntranscoefs;

   assert(vars != NULL);
   assert(offsets != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(termbegins != NULL);
   assert(nlhdlrexprdata != NULL);

   ntranscoefs = termbegins[nterms];

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexprdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, vars, nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->offsets, offsets, nterms) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefs, transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefsidx, transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->termbegins, termbegins, nterms + 1) );
   (*nlhdlrexprdata)->nvars = nvars;
   (*nlhdlrexprdata)->nterms = nterms;

   (*nlhdlrexprdata)->disrow = NULL;
   (*nlhdlrexprdata)->disvars = NULL;

   (*nlhdlrexprdata)->varvals = NULL;
   (*nlhdlrexprdata)->disvarvals = NULL;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "created nlhdlr data for the following soc expression:\n");
   printNlhdlrExprData(scip, *nlhdlrexprdata);
   /* SCIPdebugMsg(scip, "x is %p\n", (void *)vars[0]); */
#endif

   return SCIP_OKAY;
}

/** helper method to free nonlinear handler expression data */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata      /**< pointer to free nonlinear handler expression data */
   )
{
   int ntranscoefs;

   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   /* free variables and row for cone disaggregation */
   SCIP_CALL( freeDisaggrVars(scip, *nlhdlrexprdata) );

   ntranscoefs = (*nlhdlrexprdata)->termbegins[(*nlhdlrexprdata)->nterms];

   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexprdata)->varvals, (*nlhdlrexprdata)->nvars);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->termbegins, (*nlhdlrexprdata)->nterms + 1);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefsidx, ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefs, ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->offsets, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, (*nlhdlrexprdata)->nvars);
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** set varvalrs in nlhdlrexprdata to values from given SCIP solution */
static
void updateVarVals(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_SOL*             sol,                /**< SCIP solution */
   SCIP_Bool             roundtinyfrac       /**< whether values close to integers should be rounded */
   )
{
   int i;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->varvals != NULL);

   /* update varvals */
   for( i = 0; i < nlhdlrexprdata->nvars; ++i )
   {
      nlhdlrexprdata->varvals[i] = SCIPgetSolVal(scip, sol, SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[i]));
      if( roundtinyfrac && SCIPisIntegral(scip, nlhdlrexprdata->varvals[i]) )
         nlhdlrexprdata->varvals[i] = SCIPround(scip, nlhdlrexprdata->varvals[i]);
   }

   /* update disvarvals (in unittests, this may be NULL even though nterms > 1 */
   if( nlhdlrexprdata->disvarvals != NULL )
      for( i = 0; i < nlhdlrexprdata->nterms - 1; ++i )
      {
         nlhdlrexprdata->disvarvals[i] = SCIPgetSolVal(scip, sol, nlhdlrexprdata->disvars[i]);
         if( roundtinyfrac && SCIPisIntegral(scip, nlhdlrexprdata->disvarvals[i]) )
            nlhdlrexprdata->disvarvals[i] = SCIPround(scip, nlhdlrexprdata->disvarvals[i]);
      }
}

/** evaluate a single term of the form \f$v_i^T x + \beta_i\f$ */
static
SCIP_Real evalSingleTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   int                   k                   /**< term to be evaluated */
   )
{
   SCIP_Real result;
   int i;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(0 <= k && k < nlhdlrexprdata->nterms);

   result = nlhdlrexprdata->offsets[k];

   for( i = nlhdlrexprdata->termbegins[k]; i < nlhdlrexprdata->termbegins[k + 1]; ++i )
      result += nlhdlrexprdata->transcoefs[i] * nlhdlrexprdata->varvals[nlhdlrexprdata->transcoefsidx[i]];

   return result;
}

/** computes gradient cut for a 2D or 3D SOC
 *
 *  A 3D SOC looks like
 *  \f[
 *    \sqrt{ (v_1^T x + \beta_1)^2 + (v_2^T x + \beta_2)^2 } \leq v_3^T x + \beta_3
 *  \f]
 *
 *  Let \f$f(x)\f$ be the left-hand-side. The partial derivatives of \f$f\f$ are given by
 *  \f[
 *    \frac{\delta f}{\delta x_j} = \frac{(v_1)_j(v_1^T x + \beta_1) + (v_2)_j (v_2^T x + \beta_2)}{f(x)}
 *  \f]
 *
 *  and the gradient cut is then \f$f(x^*) + \nabla f(x^*)(x - x^*) \leq v_3^T x + \beta_3\f$.
 *
 *  If \f$\beta_1 = \beta_2 = 0\f$, then the constant on the left-hand-side of the cut becomes zero:
 *  \f[
 *    f(x^*) - (\frac{(v_1)_j v_1^T x^* + (v_2)_j v_2^T x^*}{f(x^*)})_j^T x^*
 *    = f(x^*) - \frac{1}{f(x^*)} \sum_j ((v_1)_j x_j^* v_1^T x^* + (v_2)_j x_j^* v_2^T x^*)
 *    = f(x^*) - \frac{1}{f(x^*)} ((v_1^T x^*)^2 + (v_2^T x^*)^2)
 *    = f(x^*) - \frac{1}{f(x^*)} f(x^*)^2 = 0
 *  \f]
 *
 *  A 2D SOC is
 *  \f[
 *    |v_1^T x + \beta_1| \leq v_2^T x + \beta_2
 *  \f]
 *  but we build the cut using the same procedure as for 3D.
 */
static
SCIP_RETCODE generateCutSolSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep,            /**< buffer to store rowprep with cut data */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_CONS*            cons,               /**< the constraint that expr is part of */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   SCIP_Real             mincutviolation,    /**< minimal required cut violation */
   SCIP_Real             rhsval              /**< value of last term at sol */
   )
{
   SCIP_Real* transcoefs;
   SCIP_Real cutcoef;
   SCIP_Real fvalue;
   SCIP_Real valterms[2] = {0.0, 0.0}; /* for lint */
   SCIP_Real cutrhs;
   SCIP_EXPR** vars;
   SCIP_VAR* cutvar;
   SCIP_Bool offsetzero;
   int* transcoefsidx;
   int* termbegins;
   int nterms;
   int i;
   int j;

   assert(rowprep != NULL);
   assert(expr != NULL);
   assert(cons != NULL);
   assert(nlhdlrexprdata != NULL);

   vars = nlhdlrexprdata->vars;
   transcoefs = nlhdlrexprdata->transcoefs;
   transcoefsidx = nlhdlrexprdata->transcoefsidx;
   termbegins = nlhdlrexprdata->termbegins;
   nterms = nlhdlrexprdata->nterms;

   *rowprep = NULL;

   /* evaluate lhs terms and compute f(x*), check whether both beta_1 and beta_2 are zero */
   fvalue = 0.0;
   offsetzero = TRUE;
   for( i = 0; i < nterms - 1; ++i )
   {
      valterms[i] = evalSingleTerm(scip, nlhdlrexprdata, i);
      fvalue += SQR( valterms[i] );
      if( nlhdlrexprdata->offsets[i] != 0.0 )
         offsetzero = FALSE;
   }
   fvalue = sqrt(fvalue);

   /* don't generate cut if we are not violated @todo: remove this once core detects better when a nlhdlr's cons is
    * violated
    */
   if( fvalue - rhsval <= mincutviolation )
   {
      SCIPdebugMsg(scip, "do not generate cut: rhsval %g, fvalue %g violation is %g\n", rhsval, fvalue, fvalue - rhsval);
      return SCIP_OKAY;
   }

   /* if f(x*) = 0 then we are at top of cone, where we cannot generate cut */
   if( SCIPisZero(scip, fvalue) )
   {
      SCIPdebugMsg(scip, "do not generate cut for lhs=%g, cannot linearize at top of cone\n", fvalue);
      return SCIP_OKAY;
   }

   /* create cut */
   SCIP_CALL( SCIPcreateRowprep(scip, rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( SCIPensureRowprepSize(scip, *rowprep, termbegins[nterms]) );

   /* cut is f(x*) + \nabla f(x*)^T (x - x*) \leq v_n^T x + \beta_n, i.e.,
    * \nabla f(x*)^T x - v_n^T x \leq \beta_n + \nabla f(x*)^T x* - f(x*)
    * thus cutrhs is \beta_n - f(x*) + \nabla f(x*)^T x*
    * if offsetzero, then we make sure that cutrhs is exactly \beta_n
    */
   cutrhs = nlhdlrexprdata->offsets[nterms - 1];
   if( !offsetzero )
      cutrhs -= fvalue;

   /* add cut coefficients from lhs terms and compute cut's rhs */
   for( j = 0; j < nterms - 1; ++j )
   {
      for( i = termbegins[j]; i < termbegins[j + 1]; ++i )
      {
         cutvar = SCIPgetExprAuxVarNonlinear(vars[transcoefsidx[i]]);

         /* cutcoef is (the first part of) the partial derivative w.r.t cutvar */
         cutcoef = transcoefs[i] * valterms[j] / fvalue;

         SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, cutvar, cutcoef) );

         if( !offsetzero )
            cutrhs += cutcoef * nlhdlrexprdata->varvals[transcoefsidx[i]];
      }
   }

   /* add terms for v_n */
   for( i = termbegins[nterms - 1]; i < termbegins[nterms]; ++i )
   {
      cutvar = SCIPgetExprAuxVarNonlinear(vars[transcoefsidx[i]]);
      SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, cutvar, -transcoefs[i]) );
   }

   /* add side */
   SCIProwprepAddSide(*rowprep, cutrhs);

   /* set name */
   (void) SCIPsnprintf(SCIProwprepGetName(*rowprep), SCIP_MAXSTRLEN, "soc%d_%p_%" SCIP_LONGINT_FORMAT, nterms, (void*) expr, SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

/** helper method to compute and add a gradient cut for the k-th cone disaggregation
 *
 *  After the SOC constraint \f$\sqrt{\sum_{i = 0}^{n-1} (v_i^T x + \beta_i)^2} \leq v_n^T x + \beta_n\f$
 *  has been disaggregated into the row \f$\sum_{i = 0}^{n-1} y_i \leq v_n^T x + \beta_n\f$ and the smaller SOC constraints
 *  \f[
 *    (v_i^T x + \beta_i)^2 \leq (v_n^T x + \beta_n) y_i \text{ for } i \in \{0, \ldots, n -1\},
 *  \f]
 *  we want to separate one of the small rotated cones.
 *  We first transform it into standard form:
 *  \f[
 *    \sqrt{4(v_i^T x + \beta_i)^2 + (v_n^T x + \beta_n - y_i)^2} - v_n^T x - \beta_n - y_i \leq 0.
 *  \f]
 *  Let \f$f(x,y)\f$ be the left-hand-side. We now compute the gradient by
 *  \f{align*}{
 *    \frac{\delta f}{\delta x_j} &= \frac{(v_i)_j(4v_i^T x + 4\beta_i) + (v_n)_j(v_n^T x + \beta_n - y_i)}{\sqrt{4(v_i^T x + \beta_i)^2 + (v_n^T x + \beta_n - y_i)^2}} - (v_n)_j \\
 *    \frac{\delta f}{\delta y_i} &= \frac{y_i - v_n^T x -\beta_n}{\sqrt{4(v_i^T x + \beta_i)^2 + (v_n^T x + \beta_n - y_i)^2}} - 1
 *  \f}
 *  and the gradient cut is then \f$f(x^*, y^*) + \nabla f(x^*,y^*)((x,y) - (x^*, y^*)) \leq 0\f$.
 *
 *  As in \ref generateCutSolSOC(), the cut constant is zero if \f$\beta_i = \beta_n = 0\f$.
 */
static
SCIP_RETCODE generateCutSolDisagg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWPREP**        rowprep,            /**< buffer to store rowprep with cut data */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_CONS*            cons,               /**< the constraint that expr is part of */
   SCIP_NLHDLREXPRDATA*  nlhdlrexprdata,     /**< nonlinear handler expression data */
   int                   disaggidx,          /**< index of disaggregation to separate */
   SCIP_Real             mincutviolation,    /**< minimal required cut violation */
   SCIP_Real             rhsval              /**< value of the rhs term */
   )
{
   SCIP_EXPR** vars;
   SCIP_VAR** disvars;
   SCIP_Real* transcoefs;
   int* transcoefsidx;
   int* termbegins;
   SCIP_VAR* cutvar;
   SCIP_Real cutcoef;
   SCIP_Real fvalue;
   SCIP_Real disvarval;
   SCIP_Real lhsval;
   SCIP_Real constant;
   SCIP_Real denominator;
   SCIP_Bool offsetzero;
   int ncutvars;
   int nterms;
   int i;

   assert(rowprep != NULL);
   assert(expr != NULL);
   assert(cons != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(disaggidx < nlhdlrexprdata->nterms-1);

   vars = nlhdlrexprdata->vars;
   disvars = nlhdlrexprdata->disvars;
   transcoefs = nlhdlrexprdata->transcoefs;
   transcoefsidx = nlhdlrexprdata->transcoefsidx;
   termbegins = nlhdlrexprdata->termbegins;
   nterms = nlhdlrexprdata->nterms;

   /* nterms is equal to n in the description and disaggidx is in {0, ..., n - 1} */

   *rowprep = NULL;

   disvarval = nlhdlrexprdata->disvarvals[disaggidx];

   lhsval = evalSingleTerm(scip, nlhdlrexprdata, disaggidx);

   denominator = sqrt(4.0 * SQR(lhsval) + SQR(rhsval - disvarval));

   /* compute value of function to be separated (f(x*,y*)) */
   fvalue = denominator - rhsval - disvarval;

   /* if the disagg soc is not violated don't compute cut */
   if( fvalue <= mincutviolation )
   {
      SCIPdebugMsg(scip, "skip cut on disaggregation index %d as violation=%g below minviolation %g\n", disaggidx,
            fvalue, mincutviolation);
      return SCIP_OKAY;
   }

   /* if the denominator is 0 -> the constraint can't be violated, and the gradient is infinite */
   if( SCIPisZero(scip, denominator) )
   {
      SCIPdebugMsg(scip, "skip cut on disaggregation index %d as we are on top of cone (denom=%g)\n", disaggidx, denominator);
      return SCIP_OKAY;
   }

   /* compute upper bound on the number of variables in cut: vars in rhs + vars in term + disagg var */
   ncutvars = (termbegins[nterms] - termbegins[nterms-1]) + (termbegins[disaggidx + 1] - termbegins[disaggidx]) + 1;

   /* create cut */
   SCIP_CALL( SCIPcreateRowprep(scip, rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( SCIPensureRowprepSize(scip, *rowprep, ncutvars) );

   /* check whether offsets (beta) are zero, so we can know cut constant will be zero */
   offsetzero = nlhdlrexprdata->offsets[disaggidx] == 0.0 && nlhdlrexprdata->offsets[nterms-1] == 0.0;

   /* constant will be grad_f(x*,y*)^T  (x*, y*) */
   constant = 0.0;

   /* a variable could appear on the lhs and rhs, but we add the coefficients separately  */

   /* add terms for v_disaggidx */
   for( i = termbegins[disaggidx]; i < termbegins[disaggidx + 1]; ++i )
   {
      cutvar = SCIPgetExprAuxVarNonlinear(vars[transcoefsidx[i]]);
      assert(cutvar != NULL);

      /* cutcoef is (the first part of) the partial derivative w.r.t cutvar */
      cutcoef = 4.0 * lhsval * transcoefs[i] / denominator;

      SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, cutvar, cutcoef) );

      if( !offsetzero )
         constant += cutcoef * nlhdlrexprdata->varvals[transcoefsidx[i]];
   }

   /* add terms for v_n */
   for( i = termbegins[nterms - 1]; i < termbegins[nterms]; ++i )
   {
      cutvar = SCIPgetExprAuxVarNonlinear(vars[transcoefsidx[i]]);
      assert(cutvar != NULL);

      /* cutcoef is the (second part of) the partial derivative w.r.t cutvar */
      cutcoef = (rhsval - disvarval) * transcoefs[i] / denominator - transcoefs[i];

      SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, cutvar, cutcoef) );

      if( !offsetzero )
         constant += cutcoef * nlhdlrexprdata->varvals[transcoefsidx[i]];
   }

   /* add term for disvar: cutcoef is the the partial derivative w.r.t. the disaggregation variable */
   cutcoef = (disvarval - rhsval) / denominator - 1.0;
   cutvar = disvars[disaggidx];

   SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, cutvar, cutcoef) );

   if( !offsetzero )
   {
      constant += cutcoef * nlhdlrexprdata->disvarvals[disaggidx];

      /* add side */
      SCIProwprepAddSide(*rowprep, constant - fvalue);
   }

   /* set name */
   (void) SCIPsnprintf(SCIProwprepGetName(*rowprep), SCIP_MAXSTRLEN, "soc_%p_%d_%" SCIP_LONGINT_FORMAT, (void*) expr, disaggidx, SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

/** given a rowprep, does a number of cleanup and checks and, if successful, generate a cut to be added to the sepastorage */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_ROWPREP*         rowprep,            /**< rowprep from which to generate row and add as cut */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_CONS*            cons,               /**< constraint for which cut is generated, or NULL */
   SCIP_Bool             allowweakcuts,      /**< whether weak cuts are allowed */
   SCIP_RESULT*          result              /**< result pointer to update (set to SCIP_CUTOFF or SCIP_SEPARATED if cut is added) */
   )
{
   SCIP_ROW* cut;
   SCIP_Real cutefficacy;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(nlhdlrdata != NULL);
   assert(rowprep != NULL);
   assert(result != NULL);

   SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIPgetHugeValue(scip), &success) );

   if( !success )
   {
      SCIPdebugMsg(scip, "rowprep cleanup failed, skip cut\n");
      return SCIP_OKAY;
   }

   if( SCIPgetRowprepViolation(scip, rowprep, sol, NULL) <= SCIPgetLPFeastol(scip) )
   {
      SCIPdebugMsg(scip, "rowprep violation %g below LP feastol %g, skip cut\n",
         SCIPgetRowprepViolation(scip, rowprep, sol, NULL), SCIPgetLPFeastol(scip));
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, rowprep, cons) );

   cutefficacy = SCIPgetCutEfficacy(scip, sol, cut);

   SCIPdebugMsg(scip, "generated row for SOC, efficacy=%g, minefficacy=%g, allowweakcuts=%u\n",
      cutefficacy, nlhdlrdata->mincutefficacy, allowweakcuts);

   /* check whether cut is applicable */
   if( SCIPisCutApplicable(scip, cut) && (allowweakcuts || cutefficacy >= nlhdlrdata->mincutefficacy) )
   {
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPaddRow(scip, cut, FALSE, &infeasible) );

#ifdef SCIP_CONSNONLINEAR_ROWNOTREMOVABLE
      /* mark row as not removable from LP for current node, if in enforcement (==addbranchscores) (this can prevent some cycling) */
      if( addbranchscores )
         SCIPmarkRowNotRemovableLocal(scip, row);
#endif

      if( infeasible )
         *result = SCIP_CUTOFF;
      else
         *result = SCIP_SEPARATED;
   }

   /* release row */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** given a rowprep, does a number of cleanup and checks and, if successful, generate a cut to be added to the cutpool */
static
SCIP_RETCODE addCutPool(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_ROWPREP*         rowprep,            /**< rowprep from which to generate row and add as cut */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_CONS*            cons                /**< constraint for which cut is generated, or NULL */
   )
{
   SCIP_ROW* cut;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(nlhdlrdata != NULL);
   assert(rowprep != NULL);

   assert(!SCIProwprepIsLocal(rowprep));

   SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIPgetHugeValue(scip), &success) );
   /* if failed or cut is only locally valid now, then skip */
   if( !success || SCIProwprepIsLocal(rowprep) )
      return SCIP_OKAY;

   /* if row after cleanup is just a boundchange, then skip */
   if( SCIProwprepGetNVars(rowprep) <= 1 )
      return SCIP_OKAY;

   /* generate row and add to cutpool */
   SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, rowprep, cons) );

   SCIP_CALL( SCIPaddPoolCut(scip, cut) );

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** checks if an expression is quadratic and collects all occurring expressions
 *
 * @pre `expr2idx` and `occurringexprs` need to be initialized with capacity 2 * nchildren
 *
 * @note We assume that a linear term always appears before its corresponding
 * quadratic term in quadexpr; this should be ensured by canonicalize
 */
static
SCIP_RETCODE checkAndCollectQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            quadexpr,           /**< candidate for a quadratic expression */
   SCIP_HASHMAP*         expr2idx,           /**< hashmap to store expressions */
   SCIP_EXPR**           occurringexprs,     /**< array to store expressions */
   int*                  nexprs,             /**< buffer to store number of expressions */
   SCIP_Bool*            success             /**< buffer to store whether the check was successful */
   )
{
   SCIP_EXPR** children;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(quadexpr != NULL);
   assert(expr2idx != NULL);
   assert(occurringexprs != NULL);
   assert(nexprs != NULL);
   assert(success != NULL);

   *nexprs = 0;
   *success = FALSE;
   children = SCIPexprGetChildren(quadexpr);
   nchildren = SCIPexprGetNChildren(quadexpr);

   /* iterate in reverse order to ensure that quadratic terms are found before linear terms */
   for( i = nchildren - 1; i >= 0; --i )
   {
      SCIP_EXPR* child;

      child = children[i];
      if( SCIPisExprPower(scip, child) )
      {
         SCIP_EXPR* childarg;

         if( SCIPgetExponentExprPow(child) != 2.0 )
            return SCIP_OKAY;

         childarg = SCIPexprGetChildren(child)[0];

         if( !SCIPhashmapExists(expr2idx, (void*) childarg) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) childarg, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = childarg;

            ++(*nexprs);
         }
      }
      else if( SCIPisExprVar(scip, child) && SCIPvarIsBinary(SCIPgetVarExprVar(child)) )
      {
         if( !SCIPhashmapExists(expr2idx, (void*) child) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) child, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = child;

            ++(*nexprs);
         }
      }
      else if( SCIPisExprProduct(scip, child) )
      {
         SCIP_EXPR* childarg1;
         SCIP_EXPR* childarg2;

         if( SCIPexprGetNChildren(child) != 2 )
            return SCIP_OKAY;

         childarg1 = SCIPexprGetChildren(child)[0];
         childarg2 = SCIPexprGetChildren(child)[1];

         if( !SCIPhashmapExists(expr2idx, (void*) childarg1) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) childarg1, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = childarg1;

            ++(*nexprs);
         }

         if( !SCIPhashmapExists(expr2idx, (void*) childarg2) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) childarg2, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = childarg2;

            ++(*nexprs);
         }
      }
      else
      {
         /* if there is a linear term without corresponding quadratic term, it is not a SOC */
         if( !SCIPhashmapExists(expr2idx, (void*) child) )
            return SCIP_OKAY;
      }
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/* builds the constraint defining matrix and vector of a quadratic expression
 *
 * @pre `quadmatrix` and `linvector` need to be initialized with size `nexprs`^2 and `nexprs`, resp.
 */
static
void buildQuadExprMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            quadexpr,           /**< the quadratic expression */
   SCIP_HASHMAP*         expr2idx,           /**< hashmap mapping the occurring expressions to their index */
   int                   nexprs,             /**< number of occurring expressions */
   SCIP_Real*            quadmatrix,         /**< pointer to store (the lower-left triangle of) the quadratic matrix */
   SCIP_Real*            linvector           /**< pointer to store the linear vector */
   )
{
   SCIP_EXPR** children;
   SCIP_Real* childcoefs;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(quadexpr != NULL);
   assert(expr2idx != NULL);
   assert(quadmatrix != NULL);
   assert(linvector != NULL);

   children = SCIPexprGetChildren(quadexpr);
   nchildren = SCIPexprGetNChildren(quadexpr);
   childcoefs = SCIPgetCoefsExprSum(quadexpr);

   /* iterate over children to build the constraint defining matrix and vector */
   for( i = 0; i < nchildren; ++i )
   {
      int varpos;

      if( SCIPisExprPower(scip, children[i]) )
      {
         assert(SCIPgetExponentExprPow(children[i]) == 2.0);
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPexprGetChildren(children[i])[0]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPexprGetChildren(children[i])[0]);
         assert(0 <= varpos && varpos < nexprs);

         quadmatrix[varpos * nexprs + varpos] = childcoefs[i];
      }
      else if( SCIPisExprVar(scip, children[i]) && SCIPvarIsBinary(SCIPgetVarExprVar(children[i])) )
      {
         assert(SCIPhashmapExists(expr2idx, (void*) children[i]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);
         assert(0 <= varpos && varpos < nexprs);

         quadmatrix[varpos * nexprs + varpos] = childcoefs[i];
      }
      else if( SCIPisExprProduct(scip, children[i]) )
      {
         int varpos2;

         assert(SCIPexprGetNChildren(children[i]) == 2);
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPexprGetChildren(children[i])[0]));
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPexprGetChildren(children[i])[1]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPexprGetChildren(children[i])[0]);
         assert(0 <= varpos && varpos < nexprs);

         varpos2 = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPexprGetChildren(children[i])[1]);
         assert(0 <= varpos2 && varpos2 < nexprs);
         assert(varpos != varpos2);

         /* Lapack uses only the lower left triangle of the symmetric matrix */
         quadmatrix[MIN(varpos, varpos2) * nexprs + MAX(varpos, varpos2)] = childcoefs[i] / 2.0;
      }
      else
      {
         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);
         assert(0 <= varpos && varpos < nexprs);

         linvector[varpos] = childcoefs[i];
      }
   }
}

/** tries to fill the nlhdlrexprdata for a potential quadratic SOC expression
 *
 * We say "try" because the expression might still turn out not to be a SOC at this point.
 */
static
SCIP_RETCODE tryFillNlhdlrExprDataQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           occurringexprs,     /**< array of all occurring expressions (nvars many) */
   SCIP_Real*            eigvecmatrix,       /**< array containing the Eigenvectors */
   SCIP_Real*            eigvals,            /**< array containing the Eigenvalues */
   SCIP_Real*            bp,                 /**< product of linear vector b * P (eigvecmatrix^t) */
   int                   nvars,              /**< number of variables */
   int*                  termbegins,         /**< pointer to store the termbegins */
   SCIP_Real*            transcoefs,         /**< pointer to store the transcoefs */
   int*                  transcoefsidx,      /**< pointer to store the transcoefsidx */
   SCIP_Real*            offsets,            /**< pointer to store the offsets */
   SCIP_Real*            lhsconstant,        /**< pointer to store the lhsconstant */
   int*                  nterms,             /**< pointer to store the total number of terms */
   SCIP_Bool*            success             /**< whether the expression is indeed a SOC */
   )
{
   SCIP_Real sqrteigval;
   int nextterm = 0;
   int nexttranscoef = 0;
   int specialtermidx;
   int i;
   int j;

   assert(scip != NULL);
   assert(occurringexprs != NULL);
   assert(eigvecmatrix != NULL);
   assert(eigvals != NULL);
   assert(bp != NULL);
   assert(termbegins != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(offsets != NULL);
   assert(lhsconstant != NULL);
   assert(success != NULL);

   *success = FALSE;
   *nterms = 0;

   /* we have lhsconstant + x^t A x + b x <= 0 and A has a single negative eigenvalue; try to build soc;
    * we now store all the v_i^T x + beta_i on the lhs, and compute the constant
    */
   specialtermidx = -1;
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPisZero(scip, eigvals[i]) )
         continue;

      if( eigvals[i] < 0.0 )
      {
         assert(specialtermidx == -1); /* there should only be one negative eigenvalue */

         specialtermidx = i;

         *lhsconstant -= bp[i] * bp[i] / (4.0 * eigvals[i]);

         continue;
      }

      assert(eigvals[i] > 0.0);
      sqrteigval = sqrt(eigvals[i]);

      termbegins[nextterm] = nexttranscoef;
      offsets[nextterm] = bp[i] / (2.0 * sqrteigval);
      *lhsconstant -= bp[i] * bp[i] / (4.0 * eigvals[i]);

      /* set transcoefs */
      for( j = 0; j < nvars; ++j )
      {
         if( !SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
         {
            transcoefs[nexttranscoef] = sqrteigval * eigvecmatrix[i * nvars + j];
            transcoefsidx[nexttranscoef] = j;

            ++nexttranscoef;
         }
      }
      ++nextterm;
   }
   assert(specialtermidx > -1);

   /* process constant; if constant is negative -> no soc */
   if( SCIPisNegative(scip, *lhsconstant) )
      return SCIP_OKAY;

   /* we need lhsconstant to be >= 0 */
   if( *lhsconstant < 0.0 )
      *lhsconstant = 0.0;

   /* store constant term */
   if( *lhsconstant > 0.0 )
   {
      termbegins[nextterm] = nexttranscoef;
      offsets[nextterm] = sqrt(*lhsconstant);
      ++nextterm;
   }

   /* now process rhs */
   {
      SCIP_Real rhstermlb;
      SCIP_Real rhstermub;
      SCIP_Real signfactor;

      assert(-eigvals[specialtermidx] > 0.0);
      sqrteigval = sqrt(-eigvals[specialtermidx]);

      termbegins[nextterm] = nexttranscoef;
      offsets[nextterm] = -bp[specialtermidx] / (2.0 * sqrteigval);

      /* the expression can only be an soc if the resulting rhs term does not change sign;
       * the rhs term is a linear combination of variables, so estimate its bounds
       */
      rhstermlb = offsets[nextterm];
      for( j = 0; j < nvars; ++j )
      {
         SCIP_INTERVAL activity;
         SCIP_Real aux;

         if( SCIPisZero(scip, eigvecmatrix[specialtermidx * nvars + j]) )
            continue;

         SCIP_CALL( SCIPevalExprActivity(scip, occurringexprs[j]) );
         activity = SCIPexprGetActivity(occurringexprs[j]);

         if( eigvecmatrix[specialtermidx * nvars + j] > 0.0 )
         {
            aux = activity.inf;
            assert(!SCIPisInfinity(scip, aux));
         }
         else
         {
            aux = activity.sup;
            assert(!SCIPisInfinity(scip, -aux));
         }

         if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
         {
            rhstermlb = -SCIPinfinity(scip);
            break;
         }
         else
            rhstermlb += sqrteigval * eigvecmatrix[specialtermidx * nvars + j] * aux;
      }

      rhstermub = offsets[nextterm];
      for( j = 0; j < nvars; ++j )
      {
         SCIP_INTERVAL activity;
         SCIP_Real aux;

         if( SCIPisZero(scip, eigvecmatrix[specialtermidx * nvars + j]) )
            continue;

         SCIP_CALL( SCIPevalExprActivity(scip, occurringexprs[j]) );
         activity = SCIPexprGetActivity(occurringexprs[j]);

         if( eigvecmatrix[specialtermidx * nvars + j] > 0.0 )
         {
            aux = activity.sup;
            assert(!SCIPisInfinity(scip, -aux));
         }
         else
         {
            aux = activity.inf;
            assert(!SCIPisInfinity(scip, aux));
         }

         if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
         {
            rhstermub = SCIPinfinity(scip);
            break;
         }
         else
            rhstermub += sqrteigval * eigvecmatrix[specialtermidx * nvars + j] * aux;
      }

      /* since we are just interested in obtaining an interval that contains the real bounds
       * and is tight enough so that we can identify that the rhsvar does not change sign,
       * we swap the bounds in case of numerical troubles
       */
      if( rhstermub < rhstermlb )
      {
         assert(SCIPisEQ(scip, rhstermub, rhstermlb));
         SCIPswapReals(&rhstermub, &rhstermlb);
      }

      /* if rhs changes sign -> not a SOC */
      if( SCIPisLT(scip, rhstermlb, 0.0) && SCIPisGT(scip, rhstermub, 0.0) )
         return SCIP_OKAY;

      signfactor = SCIPisLE(scip, rhstermub, 0.0) ? -1.0 : 1.0;

      offsets[nextterm] *= signfactor;

      /* set transcoefs for rhs term */
      for( j = 0; j < nvars; ++j )
      {
         if( SCIPisZero(scip, eigvecmatrix[specialtermidx * nvars + j]) )
            continue;

         transcoefs[nexttranscoef] = signfactor * sqrteigval * eigvecmatrix[specialtermidx * nvars + j];
         transcoefsidx[nexttranscoef] = j;

         ++nexttranscoef;
      }

      /* if rhs is a constant this method shouldn't have been called */
      assert(nexttranscoef > termbegins[nextterm]);

      /* finish processing term */
      ++nextterm;
   }

   *nterms = nextterm;

   /* sentinel value */
   termbegins[nextterm] = nexttranscoef;

   *success = TRUE;

   return SCIP_OKAY;
}

/** detects if expr &le; auxvar is of the form sqrt(sum_i coef_i (expr_i + shift_i)^2 + const) &le; auxvar
 *
 * @note if a user inputs the above expression with `const` = -epsilon, then `const` is going to be set to 0.
 */
static
SCIP_RETCODE detectSocNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_EXPR** children;
   SCIP_EXPR* child;
   SCIP_EXPR** vars;
   SCIP_HASHMAP* expr2idx;
   SCIP_HASHSET* linexprs;
   SCIP_Real* childcoefs;
   SCIP_Real* offsets;
   SCIP_Real* transcoefs;
   SCIP_Real constant;
   SCIP_Bool issoc;
   int* transcoefsidx;
   int* termbegins;
   int nchildren;
   int nterms;
   int nvars;
   int nextentry;
   int i;

   assert(expr != NULL);
   assert(success != NULL);

   *success = FALSE;
   issoc = TRUE;

   /* relation is not "<=" -> skip */
   if( SCIPgetExprNLocksPosNonlinear(expr) == 0 )
      return SCIP_OKAY;

   assert(SCIPexprGetNChildren(expr) > 0);

   child = SCIPexprGetChildren(expr)[0];
   assert(child != NULL);

   /* check whether expression is a sqrt and has a sum as child with at least 2 children and a non-negative constant */
   if( ! SCIPisExprPower(scip, expr)
      || SCIPgetExponentExprPow(expr) != 0.5
      || !SCIPisExprSum(scip, child)
      || SCIPexprGetNChildren(child) < 2
      || SCIPgetConstantExprSum(child) < 0.0)
   {
      return SCIP_OKAY;
   }

   /* assert(SCIPvarGetLbLocal(auxvar) >= 0.0); */

   /* get children of the sum */
   children = SCIPexprGetChildren(child);
   nchildren = SCIPexprGetNChildren(child);
   childcoefs = SCIPgetCoefsExprSum(child);

   /* TODO: should we initialize the hashmap with size SCIPgetNVars() so that it never has to be resized? */
   SCIP_CALL( SCIPhashmapCreate(&expr2idx, SCIPblkmem(scip), nchildren) );
   SCIP_CALL( SCIPhashsetCreate(&linexprs, SCIPblkmem(scip), nchildren) );

   /* we create coefs array here already, since we have to fill it in first loop in case of success
    * +1 for auxvar
    */
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, nchildren+1) );

   nterms = 0;

   /* check if all children are squares or linear terms with matching square term:
    * if the i-th child is (pow, expr, 2) we store the association <|expr -> i|> in expr2idx and if expr was in
    * linexprs, we remove it from there.
    * if the i-th child is expr' (different from (pow, expr, 2)) and expr' is not a key of expr2idx, we add it
    * to linexprs.
    * if at the end there is any expr in linexpr -> we do not have a separable quadratic function.
    */
   for( i = 0; i < nchildren; ++i )
   {
      /* handle quadratic expressions children */
      if( SCIPisExprPower(scip, children[i]) && SCIPgetExponentExprPow(children[i]) == 2.0 )
      {
         SCIP_EXPR* squarearg = SCIPexprGetChildren(children[i])[0];

         if( !SCIPhashmapExists(expr2idx, (void*) squarearg) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void *) squarearg, nterms) );
         }

         if( childcoefs[i] < 0.0 )
         {
            issoc = FALSE;
            break;
         }
         transcoefs[nterms] = sqrt(childcoefs[i]);

         SCIP_CALL( SCIPhashsetRemove(linexprs, (void*) squarearg) );
         ++nterms;
      }
      /* handle binary variable children */
      else if( SCIPisExprVar(scip, children[i]) && SCIPvarIsBinary(SCIPgetVarExprVar(children[i])) )
      {
         assert(!SCIPhashmapExists(expr2idx, (void*) children[i]));
         assert(!SCIPhashsetExists(linexprs, (void*) children[i]));

         SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void *) children[i], nterms) );

         if( childcoefs[i] < 0.0 )
         {
            issoc = FALSE;
            break;
         }
         transcoefs[nterms] = sqrt(childcoefs[i]);

         ++nterms;
      }
      else
      {
         if( !SCIPhashmapExists(expr2idx, (void*) children[i]) )
         {
            SCIP_CALL( SCIPhashsetInsert(linexprs, SCIPblkmem(scip), (void*) children[i]) );
         }
      }
   }

   /* there are linear terms without corresponding quadratic terms or it was detected not to be soc */
   if( SCIPhashsetGetNElements(linexprs) > 0 || ! issoc )
   {
      SCIPfreeBufferArray(scip, &transcoefs);
      SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   /* add one to terms counter for auxvar */
   ++nterms;

   constant = SCIPgetConstantExprSum(child);

   /* compute constant of possible soc expression to check its sign */
   for( i = 0; i < nchildren; ++i )
   {
      if( ! SCIPisExprPower(scip, children[i]) || SCIPgetExponentExprPow(children[i]) != 2.0 )
      {
         int auxvarpos;

         assert(SCIPhashmapExists(expr2idx, (void*) children[i]) );
         auxvarpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);

         constant -= SQR(0.5 * childcoefs[i] / transcoefs[auxvarpos]);
      }
   }

   /* if the constant is negative -> no SOC */
   if( SCIPisNegative(scip, constant) )
   {
      SCIPfreeBufferArray(scip, &transcoefs);
      SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }
   else if( SCIPisZero(scip, constant) )
      constant = 0.0;
   assert(constant >= 0.0);

   /* at this point, we have found an SOC structure */
   *success = TRUE;

   nvars = nterms;

   /* add one to terms counter for constant term */
   if( constant > 0.0 )
      ++nterms;

   /* allocate temporary memory to collect data */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, nterms + 1) );

   /* fill in data for non constant terms of lhs; initialize their offsets */
   for( i = 0; i < nvars - 1; ++i )
   {
      transcoefsidx[i] = i;
      termbegins[i] = i;
      offsets[i] = 0.0;
   }

   /* add constant term and rhs */
   vars[nvars - 1] = expr;
   if( constant > 0.0 )
   {
      /* constant term */
      termbegins[nterms - 2] = nterms - 2;
      offsets[nterms - 2] = sqrt(constant);

      /* rhs */
      termbegins[nterms - 1] = nterms - 2;
      offsets[nterms - 1] = 0.0;
      transcoefsidx[nterms - 2] = nvars - 1;
      transcoefs[nterms - 2] = 1.0;

      /* sentinel value */
      termbegins[nterms] = nterms - 1;
   }
   else
   {
      /* rhs */
      termbegins[nterms - 1] = nterms - 1;
      offsets[nterms - 1] = 0.0;
      transcoefsidx[nterms - 1] = nvars - 1;
      transcoefs[nterms - 1] = 1.0;

      /* sentinel value */
      termbegins[nterms] = nterms;
   }

   /* request required auxiliary variables and fill vars and offsets array */
   nextentry = 0;
   for( i = 0; i < nchildren; ++i )
   {
      if( SCIPisExprPower(scip, children[i]) && SCIPgetExponentExprPow(children[i]) == 2.0 )
      {
         SCIP_EXPR* squarearg;

         squarearg = SCIPexprGetChildren(children[i])[0];
         assert(SCIPhashmapGetImageInt(expr2idx, (void*) squarearg) == nextentry);

         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, squarearg, TRUE, FALSE, FALSE, FALSE) );

         vars[nextentry] = squarearg;
         ++nextentry;
      }
      else if( SCIPisExprVar(scip, children[i]) && SCIPvarIsBinary(SCIPgetVarExprVar(children[i])) )
      {
         /* handle binary variable children: no need to request auxvar */
         assert(SCIPhashmapGetImageInt(expr2idx, (void*) children[i]) == nextentry);
         vars[nextentry] = children[i];
         ++nextentry;
      }
      else
      {
         int auxvarpos;

         assert(SCIPhashmapExists(expr2idx, (void*) children[i]));
         auxvarpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);

         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, children[i], TRUE, FALSE, FALSE, FALSE) );

         offsets[auxvarpos] = 0.5 * childcoefs[i] / transcoefs[auxvarpos];
      }
   }
   assert(nextentry == nvars - 1);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n", (void*)expr);
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, " <= auxvar\n");
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, offsets, transcoefs, transcoefsidx, termbegins,
            nvars, nterms, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

   /* free memory */
   SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
   SCIPhashmapFree(&expr2idx);
   SCIPfreeBufferArray(scip, &termbegins);
   SCIPfreeBufferArray(scip, &transcoefsidx);
   SCIPfreeBufferArray(scip, &offsets);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &transcoefs);

   return SCIP_OKAY;
}

/** helper method to detect c + sum_i coef_i expr_i^2 - coef_k expr_k^2 &le; 0
 *  and c + sum_i coef_i expr_i^2 - coef_k expr_k expr_l &le; 0
 *
 *  binary linear variables are interpreted as quadratic terms
 *
 *  @todo: extend this function to detect  c + sum_i coef_i (expr_i + const_i)^2 - ...
 *  this would probably share a lot of code with detectSocNorm
 */
static
SCIP_RETCODE detectSocQuadraticSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            enforcebelow,       /**< pointer to store whether we enforce <= (TRUE) or >= (FALSE); only valid when success is TRUE */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_EXPR** children;
   SCIP_EXPR** vars = NULL;
   SCIP_Real* childcoefs;
   SCIP_Real* offsets = NULL;
   SCIP_Real* transcoefs = NULL;
   int* transcoefsidx = NULL;
   int* termbegins = NULL;
   SCIP_Real constant;
   SCIP_Real lhsconstant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real rhssign;
   SCIP_INTERVAL expractivity;
   int ntranscoefs;
   int nposquadterms;
   int nnegquadterms;
   int nposbilinterms;
   int nnegbilinterms;
   int rhsidx;
   int lhsidx;
   int specialtermidx;
   int nchildren;
   int nnzinterms;
   int nterms;
   int nvars;
   int nextentry;
   int i;
   SCIP_Bool ishyperbolic;

   assert(expr != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is a sum with at least 2 children */
   if( ! SCIPisExprSum(scip, expr) || SCIPexprGetNChildren(expr) < 2 )
      return SCIP_OKAY;

   /* get children of the sum */
   children = SCIPexprGetChildren(expr);
   nchildren = SCIPexprGetNChildren(expr);
   constant = SCIPgetConstantExprSum(expr);

   /* we duplicate the child coefficients since we have to manipulate them */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &childcoefs, SCIPgetCoefsExprSum(expr), nchildren) ); /*lint !e666*/

   /* initialize data */
   lhsidx = -1;
   rhsidx = -1;
   nposquadterms = 0;
   nnegquadterms = 0;
   nposbilinterms = 0;
   nnegbilinterms = 0;

   /* check if all children are quadratic or binary linear and count number of positive and negative terms */
   for( i = 0; i < nchildren; ++i )
   {
      if( SCIPisExprPower(scip, children[i]) && SCIPgetExponentExprPow(children[i]) == 2.0 )
      {
         if( childcoefs[i] > 0.0 )
         {
            ++nposquadterms;
            lhsidx = i;
         }
         else
         {
            ++nnegquadterms;
            rhsidx = i;
         }
      }
      else if( SCIPisExprVar(scip, children[i]) && SCIPvarIsBinary(SCIPgetVarExprVar(children[i])) )
      {
         if( childcoefs[i] > 0.0 )
         {
            ++nposquadterms;
            lhsidx = i;
         }
         else
         {
            ++nnegquadterms;
            rhsidx = i;
         }
      }
      else if( SCIPisExprProduct(scip, children[i]) && SCIPexprGetNChildren(children[i]) == 2 )
      {
         if( childcoefs[i] > 0.0 )
         {
            ++nposbilinterms;
            lhsidx = i;
         }
         else
         {
            ++nnegbilinterms;
            rhsidx = i;
         }
      }
      else
      {
         goto CLEANUP;
      }

      /* more than one positive eigenvalue and more than one negative eigenvalue -> can't be convex */
      if( nposquadterms > 1 && nnegquadterms > 1 )
         goto CLEANUP;

      /* more than one bilinear term -> can't be handled by this method */
      if( nposbilinterms + nnegbilinterms > 1 )
         goto CLEANUP;

      /* one positive bilinear term and also at least one positive quadratic term -> not a simple SOC */
      if( nposbilinterms > 0 && nposquadterms > 0 )
         goto CLEANUP;

      /* one negative bilinear term and also at least one negative quadratic term -> not a simple SOC */
      if( nnegbilinterms > 0 && nnegquadterms > 0 )
         goto CLEANUP;
   }

   if( nposquadterms == nchildren || nnegquadterms == nchildren )
      goto CLEANUP;

   assert(nposquadterms <= 1 || nnegquadterms <= 1);
   assert(nposbilinterms + nnegbilinterms <= 1);
   assert(nposbilinterms == 0 || nposquadterms == 0);
   assert(nnegbilinterms == 0 || nnegquadterms == 0);

   /* if a bilinear term is involved, it is a hyperbolic expression */
   ishyperbolic = (nposbilinterms + nnegbilinterms > 0);

   if( conslhs == SCIP_INVALID || consrhs == SCIP_INVALID )  /*lint !e777*/
   {
      SCIP_CALL( SCIPevalExprActivity(scip, expr) );
      expractivity = SCIPexprGetActivity(expr);

      lhs = (conslhs == SCIP_INVALID ? expractivity.inf : conslhs); /*lint !e777*/
      rhs = (consrhs == SCIP_INVALID ? expractivity.sup : consrhs); /*lint !e777*/
   }
   else
   {
      lhs = conslhs;
      rhs = consrhs;
   }

   /* detect case and store lhs/rhs information */
   if( (ishyperbolic && nnegbilinterms > 0) || (!ishyperbolic && nnegquadterms < 2) )
   {
      /* we have -x*y + z^2 ... -> we want to write  z^2 ... <= x*y;
       * or we have -x^2 + y^2  ... -> we want to write y^2 ... <= x^2;
       * in any case, we need a finite rhs
       */
      assert(nnegbilinterms == 1 || nnegquadterms == 1);
      assert(rhsidx != -1);

      /* if rhs is infinity, it can't be soc
       * TODO: if it can't be soc, then we should enforce the caller so that we do not try the more complex quadratic
       * method
       */
      if( SCIPisInfinity(scip, rhs) )
         goto CLEANUP;

      specialtermidx = rhsidx;
      lhsconstant = constant - rhs;
      *enforcebelow = TRUE; /* enforce expr <= rhs */
   }
   else
   {
      assert(lhsidx != -1);

      /* if lhs is infinity, it can't be soc */
      if( SCIPisInfinity(scip, -lhs) )
         goto CLEANUP;

      specialtermidx = lhsidx;
      lhsconstant = lhs - constant;

      /* negate all coefficients */
      for( i = 0; i < nchildren; ++i )
         childcoefs[i] = -childcoefs[i];
      *enforcebelow = FALSE; /* enforce lhs <= expr */
   }
   assert(childcoefs[specialtermidx] != 0.0);

   if( ishyperbolic )
   {
      SCIP_INTERVAL yactivity;
      SCIP_INTERVAL zactivity;

      assert(SCIPexprGetNChildren(children[specialtermidx]) == 2);

      SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(children[specialtermidx])[0]) );
      yactivity = SCIPexprGetActivity(SCIPexprGetChildren(children[specialtermidx])[0]);

      SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(children[specialtermidx])[1]) );
      zactivity = SCIPexprGetActivity(SCIPexprGetChildren(children[specialtermidx])[1]);

      if( SCIPisNegative(scip, yactivity.inf + zactivity.inf) )
      {
         /* the sum of the expressions in the bilinear term changes sign -> no SOC */
         if( SCIPisPositive(scip, yactivity.sup + zactivity.sup) )
            goto CLEANUP;

         rhssign = -1.0;
      }
      else
         rhssign = 1.0;

      lhsconstant *= 4.0 / -childcoefs[specialtermidx];
   }
   else if( SCIPisExprVar(scip, children[specialtermidx]) )
   {
      /* children[specialtermidx] can be a variable, in which case we treat it as if it is squared */
      rhssign = 1.0;
   }
   else
   {
      SCIP_INTERVAL rhsactivity;

      assert(SCIPexprGetNChildren(children[specialtermidx]) == 1);
      SCIP_CALL( SCIPevalExprActivity(scip, SCIPexprGetChildren(children[specialtermidx])[0]) );
      rhsactivity = SCIPexprGetActivity(SCIPexprGetChildren(children[specialtermidx])[0]);

      if( rhsactivity.inf < 0.0 )
      {
         /* rhs variable changes sign -> no SOC */
         if( rhsactivity.sup > 0.0 )
            goto CLEANUP;

         rhssign = -1.0;
      }
      else
         rhssign = 1.0;
   }

   if( SCIPisNegative(scip, lhsconstant) )
      goto CLEANUP;

   if( SCIPisZero(scip, lhsconstant) )
      lhsconstant = 0.0;

   /*
    * we have found an SOC-representable expression. Now build the nlhdlrexprdata
    *
    * in the non-hyperbolic case, c + sum_i coef_i expr_i^2 - coef_k expr_k^2 <= 0 is transformed to
    * sqrt( c + sum_i coef_i expr_i^2 ) <= coef_k expr_k
    * so there are nchildren many vars, nchildren (+ 1 if c != 0) many terms, nchildren many coefficients in the vs
    * in SOC representation
    *
    * in the hyperbolic case, c + sum_i coef_i expr_i^2 - coef_k expr_k expr_l <= 0 is transformed to
    * sqrt( 4(c + sum_i coef_i expr_i^2) + (expr_k - expr_l)^2 ) <= expr_k + expr_l
    * so there are nchildren + 1many vars, nchildren + 1(+ 1 if c != 0) many terms, nchildren + 3 many coefficients in
    * the vs in SOC representation
    */

   ntranscoefs = ishyperbolic ? nchildren + 3 : nchildren;
   nvars = ishyperbolic ? nchildren + 1 : nchildren;
   nterms = nvars;

   /* constant term */
   if( lhsconstant > 0.0 )
      nterms++;

   /* SOC was detected, allocate temporary memory for data to collect */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, nterms + 1) );

   *success = TRUE;
   nextentry = 0;

   /* collect all the v_i and beta_i */
   nnzinterms = 0;
   for( i = 0; i < nchildren; ++i )
   {
      /* variable and coef for rhs have to be set to the last entry */
      if( i == specialtermidx )
         continue;

      /* extract (unique) variable appearing in term */
      if( SCIPisExprVar(scip, children[i]) )
      {
         vars[nextentry] = children[i];

         assert(SCIPvarIsBinary(SCIPgetVarExprVar(vars[nextentry])));
      }
      else
      {
         assert(SCIPisExprPower(scip, children[i]));

         /* notify that we will require auxiliary variable */
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(children[i])[0], TRUE, FALSE, FALSE, FALSE) );
         vars[nextentry] = SCIPexprGetChildren(children[i])[0];
      }
      assert(vars[nextentry] != NULL);

      /* store v_i and beta_i */
      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = 0.0;

      transcoefsidx[nnzinterms] = nextentry;
      if( ishyperbolic )
      {
         /* we eliminate the coefficient of the bilinear term to arrive at standard form */
         assert(4.0 * childcoefs[i] / -childcoefs[specialtermidx] > 0.0);
         transcoefs[nnzinterms] = sqrt(4.0 * childcoefs[i] / -childcoefs[specialtermidx]);
      }
      else
      {
         assert(childcoefs[i] > 0.0);
         transcoefs[nnzinterms] = sqrt(childcoefs[i]);
      }

      /* finish adding nonzeros */
      ++nnzinterms;

      /* finish processing term */
      ++nextentry;
   }
   assert(nextentry == nchildren - 1);

   /* store term for constant (v_i = 0) */
   if( lhsconstant > 0.0 )
   {
      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = sqrt(lhsconstant);

      /* finish processing term; this term has 0 nonzero thus we do not increase nnzinterms */
      ++nextentry;
   }

   if( !ishyperbolic )
   {
      /* store rhs term */
      if( SCIPisExprVar(scip, children[specialtermidx]) )
      {
         /* this should be the "children[specialtermidx] can be a variable, in which case we treat it as if it is squared" case */
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, children[specialtermidx], TRUE, FALSE, FALSE, FALSE) );
         vars[nvars - 1] = children[specialtermidx];
      }
      else
      {
         assert(SCIPisExprPower(scip, children[specialtermidx]));
         assert(SCIPexprGetChildren(children[specialtermidx]) != NULL);
         SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(children[specialtermidx])[0], TRUE, FALSE, FALSE, FALSE) );
         vars[nvars - 1] = SCIPexprGetChildren(children[specialtermidx])[0];
      }

      assert(childcoefs[specialtermidx] < 0.0);

      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = 0.0;
      transcoefs[nnzinterms] = rhssign * sqrt(-childcoefs[specialtermidx]);
      transcoefsidx[nnzinterms] = nvars - 1;

      /* finish adding nonzeros */
      ++nnzinterms;

      /* finish processing term */
      ++nextentry;
   }
   else
   {
      /* store last lhs term and rhs term coming from the bilinear term */
      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(children[specialtermidx])[0], TRUE, FALSE, FALSE, FALSE) );
      vars[nvars - 2] = SCIPexprGetChildren(children[specialtermidx])[0];

      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, SCIPexprGetChildren(children[specialtermidx])[1], TRUE, FALSE, FALSE, FALSE) );
      vars[nvars - 1] = SCIPexprGetChildren(children[specialtermidx])[1];

      /* at this point, vars[nvars - 2] = expr_k and vars[nvars - 1] = expr_l;
       * on the lhs we have the term (expr_k - expr_l)^2
       */
      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = 0.0;

      /* expr_k */
      transcoefsidx[nnzinterms] = nvars - 2;
      transcoefs[nnzinterms] = 1.0;
      ++nnzinterms;

      /* - expr_l */
      transcoefsidx[nnzinterms] = nvars - 1;
      transcoefs[nnzinterms] = -1.0;
      ++nnzinterms;

      /* finish processing term */
      ++nextentry;

      /* on rhs we have +/-(expr_k + expr_l) */
      termbegins[nextentry] = nnzinterms;
      offsets[nextentry] = 0.0;

      /* rhssing * expr_k */
      transcoefsidx[nnzinterms] = nvars - 2;
      transcoefs[nnzinterms] = rhssign;
      ++nnzinterms;

      /* rhssing * expr_l */
      transcoefsidx[nnzinterms] = nvars - 1;
      transcoefs[nnzinterms] = rhssign;
      ++nnzinterms;

      /* finish processing term */
      ++nextentry;
   }
   assert(nextentry == nterms);
   assert(nnzinterms == ntranscoefs);

   /* sentinel value */
   termbegins[nextentry] = nnzinterms;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n  %g <= ", (void*)expr, lhs);
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, " <= %g\n", rhs);
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, offsets, transcoefs, transcoefsidx, termbegins, nvars, nterms,
            nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

CLEANUP:
   SCIPfreeBufferArrayNull(scip, &termbegins);
   SCIPfreeBufferArrayNull(scip, &transcoefsidx);
   SCIPfreeBufferArrayNull(scip, &transcoefs);
   SCIPfreeBufferArrayNull(scip, &offsets);
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &childcoefs);

   return SCIP_OKAY;
}

/** detects complex quadratic expressions that can be represented as SOC constraints
 *
 *  These are quadratic expressions with either exactly one positive or exactly one negative eigenvalue,
 *  in addition to some extra conditions. One needs to write the quadratic as
 *  sum eigval_i (eigvec_i . x)^2 + c &le; -eigval_k (eigvec_k . x)^2, where eigval_k is the negative eigenvalue,
 *  and c must be positive and (eigvec_k . x) must not change sign.
 *  This is described in more details in
 *  Mahajan, Ashutosh & Munson, Todd, Exploiting Second-Order Cone Structure for Global Optimization, 2010.
 *
 *  The eigen-decomposition is computed using Lapack.
 *  Binary linear variables are interpreted as quadratic terms.
 *
 * @todo: In the case -b <= a + x^2 - y^2 <= b, it is possible to represent both sides by SOC. Currently, the
 * datastructure can only handle one SOC. If this should appear more often, it could be worth to extend it,
 * such that both sides can be handled (see e.g. instance chp_partload).
 * FS: this shouldn't be possible. For a <= b + x^2 - y^2 <= c to be SOC representable on both sides, we would need
 * that a - b >= 0 and b -c >= 0, but this implies that a >= c and assuming the constraint is not trivially infeasible,
 * a <= b. Thus, a = b = c and the constraint is x^2 == y^2.
 *
 * @todo: Since cons_nonlinear multiplies as many terms out as possible during presolving, some SOC-representable
 * structures cannot be detected, (see e.g. instances bearing or wager). There is currently no obvious way
 * to handle this.
 */
static
SCIP_RETCODE detectSocQuadraticComplex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            enforcebelow,       /**< pointer to store whether we enforce <= (TRUE) or >= (FALSE); only
                                              *   valid when success is TRUE */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_EXPR** occurringexprs;
   SCIP_HASHMAP* expr2idx;
   SCIP_Real* offsets;
   SCIP_Real* transcoefs;
   SCIP_Real* eigvecmatrix;
   SCIP_Real* eigvals;
   SCIP_Real* lincoefs;
   SCIP_Real* bp;
   int* transcoefsidx;
   int* termbegins;
   SCIP_Real constant;
   SCIP_Real lhsconstant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_INTERVAL expractivity;
   int nvars;
   int nterms;
   int nchildren;
   int npos;
   int nneg;
   int ntranscoefs;
   int i;
   int j;
   SCIP_Bool rhsissoc;
   SCIP_Bool lhsissoc;
   SCIP_Bool isquadratic;

   assert(expr != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is a sum with at least 2 children */
   if( ! SCIPisExprSum(scip, expr) || SCIPexprGetNChildren(expr) < 2 )
   {
      return SCIP_OKAY;
   }

   /* we need Lapack to compute eigenvalues/vectors below */
   if( ! SCIPlapackIsAvailable() )
      return SCIP_OKAY;

   /* get children of the sum */
   nchildren = SCIPexprGetNChildren(expr);
   constant = SCIPgetConstantExprSum(expr);

   /* initialize data */
   offsets = NULL;
   transcoefs = NULL;
   transcoefsidx = NULL;
   termbegins = NULL;
   bp = NULL;

   SCIP_CALL( SCIPhashmapCreate(&expr2idx, SCIPblkmem(scip), 2 * nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &occurringexprs, 2 * nchildren) );

   /* check if the expression is quadratic and collect all occurring expressions */
   SCIP_CALL( checkAndCollectQuadratic(scip, expr, expr2idx, occurringexprs, &nvars, &isquadratic) );

   if( !isquadratic )
   {
      SCIPfreeBufferArray(scip, &occurringexprs);
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   /* check that nvars*nvars doesn't get too large, see also SCIPcomputeExprQuadraticCurvature() */
   if( nvars > 7000 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "nlhdlr_soc - number of quadratic variables is too large (%d) to check the curvature\n", nvars);
      SCIPfreeBufferArray(scip, &occurringexprs);
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   assert(SCIPhashmapGetNElements(expr2idx) == nvars);

   /* create datastructures for constaint defining matrix and vector */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &eigvecmatrix, nvars * nvars) ); /*lint !e647*/
   SCIP_CALL( SCIPallocClearBufferArray(scip, &lincoefs, nvars) );

   /* build constraint defining matrix (stored in eigvecmatrix) and vector (stored in lincoefs) */
   buildQuadExprMatrix(scip, expr, expr2idx, nvars, eigvecmatrix, lincoefs);

   SCIP_CALL( SCIPallocBufferArray(scip, &eigvals, nvars) );

   /* compute eigenvalues and vectors, A = PDP^t
    * note: eigvecmatrix stores P^t, i.e., P^t_{i,j} = eigvecmatrix[i*nvars+j]
    */
   if( SCIPlapackComputeEigenvalues(SCIPbuffer(scip), TRUE, nvars, eigvecmatrix, eigvals) != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Failed to compute eigenvalues and eigenvectors for expression:\n");

#ifdef SCIP_DEBUG
      SCIPdismantleExpr(scip, NULL, expr);
#endif

      goto CLEANUP;
   }

   SCIP_CALL( SCIPallocClearBufferArray(scip, &bp, nvars) );

   nneg = 0;
   npos = 0;
   ntranscoefs = 0;

   /* set small eigenvalues to 0 and compute b*P */
   for( i = 0; i < nvars; ++i )
   {
      for( j = 0; j < nvars; ++j )
      {
         bp[i] += lincoefs[j] * eigvecmatrix[i * nvars + j];

         /* count the number of transcoefs to be used later */
         if( !SCIPisZero(scip, eigvals[i]) && !SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
            ++ntranscoefs;
      }

      if( SCIPisZero(scip, eigvals[i]) )
      {
         /* if there is a purely linear variable, the constraint can't be written as a SOC */
         if( !SCIPisZero(scip, bp[i]) )
            goto CLEANUP;

         bp[i] = 0.0;
         eigvals[i] = 0.0;
      }
      else if( eigvals[i] > 0.0 )
         npos++;
      else
         nneg++;
   }

   /* a proper SOC constraint needs at least 2 variables */
   if( npos + nneg < 2 )
      goto CLEANUP;

   /* determine whether rhs or lhs of cons is potentially SOC, if any */
   rhsissoc = (nneg == 1 && SCIPgetExprNLocksPosNonlinear(expr) > 0);
   lhsissoc = (npos == 1 && SCIPgetExprNLocksNegNonlinear(expr) > 0);

   if( rhsissoc || lhsissoc )
   {
      if( conslhs == SCIP_INVALID || consrhs == SCIP_INVALID ) /*lint !e777*/
      {
         SCIP_CALL( SCIPevalExprActivity(scip, expr) );
         expractivity = SCIPexprGetActivity(expr);
         lhs = (conslhs == SCIP_INVALID ? expractivity.inf : conslhs); /*lint !e777*/
         rhs = (consrhs == SCIP_INVALID ? expractivity.sup : consrhs); /*lint !e777*/
      }
      else
      {
         lhs = conslhs;
         rhs = consrhs;
      }
   }
   else
   {
      /* if none of the sides is potentially SOC, stop */
      goto CLEANUP;
   }

   /* @TODO: what do we do if both sides are possible? */
   if( !rhsissoc )
   {
      assert(lhsissoc);

      /* lhs is potentially SOC, change signs */
      lhsconstant = lhs - constant;  /*lint !e644*/

      for( i = 0; i < nvars; ++i )
      {
         eigvals[i] = -eigvals[i];
         bp[i] = -bp[i];
      }
      *enforcebelow = FALSE; /* enforce lhs <= expr */
   }
   else
   {
      lhsconstant = constant - rhs;  /*lint !e644*/
      *enforcebelow = TRUE; /* enforce expr <= rhs */
   }

   /* initialize remaining datastructures for nonlinear handler */
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, npos + nneg + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, npos + nneg + 2) );

   /* try to fill the nlhdlrexprdata (at this point, it can still fail) */
   SCIP_CALL( tryFillNlhdlrExprDataQuad(scip, occurringexprs, eigvecmatrix, eigvals, bp, nvars, termbegins, transcoefs,
         transcoefsidx, offsets, &lhsconstant, &nterms, success) );

   if( !(*success) )
      goto CLEANUP;

   assert(0 < nterms && nterms <= npos + nneg + 1);
   assert(ntranscoefs == termbegins[nterms]);

   /*
    * at this point, the expression passed all checks and is SOC-representable
    */

   /* register all requests for auxiliary variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPregisterExprUsageNonlinear(scip, occurringexprs[i], TRUE, FALSE, FALSE, FALSE) );
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n%f <= ", (void*)expr, lhs);
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, "<= %f\n", rhs);
#endif

   /* finally, create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, occurringexprs, offsets, transcoefs, transcoefsidx, termbegins, nvars, nterms,
            nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

CLEANUP:
   SCIPfreeBufferArrayNull(scip, &termbegins);
   SCIPfreeBufferArrayNull(scip, &transcoefsidx);
   SCIPfreeBufferArrayNull(scip, &transcoefs);
   SCIPfreeBufferArrayNull(scip, &offsets);
   SCIPfreeBufferArrayNull(scip, &bp);
   SCIPfreeBufferArray(scip, &eigvals);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &eigvecmatrix);
   SCIPfreeBufferArray(scip, &occurringexprs);
   SCIPhashmapFree(&expr2idx);

   return SCIP_OKAY;
}

/** helper method to detect SOC structures
 *
 * The detection runs in 3 steps:
 *  1. check if expression is a norm of the form \f$\sqrt{\sum_i (\text{sqrcoef}_i\, \text{expr}_i^2 + \text{lincoef}_i\, \text{expr}_i) + c}\f$
 *  which can be transformed to the form \f$\sqrt{\sum_i (\text{coef}_i \text{expr}_i + \text{const}_i)^2 + c^*}\f$ with \f$c^* \geq 0\f$.\n
 *    -> this results in the SOC     expr &le; auxvar(expr)
 *
 *    TODO we should generalize and check for sqrt(positive-semidefinite-quadratic)
 *
 *  2. check if expression represents a quadratic function of one of the following forms (all coefs > 0)
 *     1. \f$(\sum_i   \text{coef}_i \text{expr}_i^2) - \text{coef}_k \text{expr}_k^2 \leq \text{RHS}\f$ or
 *     2. \f$(\sum_i - \text{coef}_i \text{expr}_i^2) + \text{coef}_k \text{expr}_k^2 \geq \text{LHS}\f$ or
 *     3. \f$(\sum_i   \text{coef}_i \text{expr}_i^2) - \text{coef}_k \text{expr}_k \text{expr}_l \leq \text{RHS}\f$ or
 *     4. \f$(\sum_i - \text{coef}_i \text{expr}_i^2) + \text{coef}_k \text{expr}_k \text{expr}_l \geq \text{LHS}\f$,
 *
 *     where RHS &ge; 0 or LHS &le; 0, respectively. For LHS and RHS we use the constraint sides if it is a root expr
 *     and the bounds of the auxiliary variable otherwise.
 *     The last two cases are called hyperbolic or rotated second order cone.\n
 *     -> this results in the SOC \f$\sqrt{(\sum_i \text{coef}_i \text{expr}_i^2) - \text{RHS}} \leq \sqrt{\text{coef}_k} \text{expr}_k\f$
 *                            or  \f$\sqrt{4(\sum_i \text{coef}_i \text{expr}_i^2) - 4\text{RHS} + (\text{expr}_k - \text{expr}_l)^2)} \leq \text{expr}_k + \text{expr}_l\f$.
 *                            (analogously for the LHS cases)
 *
 *  3. check if expression represents a quadratic inequality of the form \f$f(x) = x^TAx + b^Tx + c \leq 0\f$ such that \f$f(x)\f$
 *  has exactly one negative eigenvalue plus some extra conditions, see detectSocQuadraticComplex().
 *
 *  Note that step 3 is only performed if parameter `compeigenvalues` is set to TRUE.
 */
static
SCIP_RETCODE detectSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_NLHDLRDATA*      nlhdlrdata,         /**< nonlinear handler data */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_NLHDLREXPRDATA** nlhdlrexprdata,     /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            enforcebelow,       /**< pointer to store whether we enforce <= (TRUE) or >= (FALSE); only
                                              *   valid when success is TRUE */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   assert(expr != NULL);
   assert(nlhdlrdata != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is given as norm as described in case 1 above: if we have a constraint
    * sqrt(sum x_i^2) <= constant, then it might be better not to handle this here; thus, we only call detectSocNorm
    * when the expr is _not_ the root of a constraint
    */
   if( conslhs == SCIP_INVALID && consrhs == SCIP_INVALID ) /*lint !e777*/
   {
      SCIP_CALL( detectSocNorm(scip, expr, nlhdlrexprdata, success) );
      *enforcebelow = *success;
   }

   if( !(*success) )
   {
      /* check whether expression is a simple soc-respresentable quadratic expression as described in case 2 above */
      SCIP_CALL( detectSocQuadraticSimple(scip, expr, conslhs, consrhs, nlhdlrexprdata, enforcebelow, success) );
   }

   if( !(*success) && nlhdlrdata->compeigenvalues )
   {
      /* check whether expression is a more complex soc-respresentable quadratic expression as described in case 3 */
      SCIP_CALL( detectSocQuadraticComplex(scip, expr, conslhs, consrhs, nlhdlrexprdata, enforcebelow, success) );
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_NLHDLRCOPYHDLR(nlhdlrCopyhdlrSoc)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPnlhdlrGetName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeNlhdlrSoc(targetscip) );

   return SCIP_OKAY;
}

/** callback to free data of handler */
static
SCIP_DECL_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataSoc)
{ /*lint --e{715}*/
   assert(nlhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataSoc)
{  /*lint --e{715}*/
   assert(*nlhdlrexprdata != NULL);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );

   return SCIP_OKAY;
}


/** callback to be called in initialization */
#if 0
static
SCIP_DECL_NLHDLRINIT(nlhdlrInitSoc)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrInitSoc NULL
#endif


/** callback to be called in deinitialization */
#if 0
static
SCIP_DECL_NLHDLREXIT(nlhdlrExitSoc)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrExitSoc NULL
#endif


/** callback to detect structure in expression tree */
static
SCIP_DECL_NLHDLRDETECT(nlhdlrDetectSoc)
{ /*lint --e{715}*/
   SCIP_Real conslhs;
   SCIP_Real consrhs;
   SCIP_Bool enforcebelow;
   SCIP_Bool success;
   SCIP_NLHDLRDATA* nlhdlrdata;

   assert(expr != NULL);

   /* don't try if no sepa is required
    * TODO implement some bound strengthening
    */
   if( (*enforcing & SCIP_NLHDLR_METHOD_SEPABOTH) == SCIP_NLHDLR_METHOD_SEPABOTH )
      return SCIP_OKAY;

   assert(SCIPgetExprNAuxvarUsesNonlinear(expr) > 0);  /* since some sepa is required, there should have been demand for it */

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   conslhs = (cons == NULL ? SCIP_INVALID : SCIPgetLhsNonlinear(cons));
   consrhs = (cons == NULL ? SCIP_INVALID : SCIPgetRhsNonlinear(cons));

   SCIP_CALL( detectSOC(scip, nlhdlrdata, expr, conslhs, consrhs, nlhdlrexprdata, &enforcebelow, &success) );

   if( !success )
      return SCIP_OKAY;

   /* inform what we can do */
   *participating = enforcebelow ? SCIP_NLHDLR_METHOD_SEPABELOW : SCIP_NLHDLR_METHOD_SEPAABOVE;

   /* if we have been successful on sqrt(...) <= auxvar, then we enforce
    * otherwise, expr is quadratic and we separate for expr <= ub(auxvar) only
    * in that case, we enforce only if expr is the root of a constraint, since then replacing auxvar by up(auxvar) does not relax anything (auxvar <= ub(auxvar) is the only constraint on auxvar)
    */
   if( (SCIPisExprPower(scip, expr) && SCIPgetExponentExprPow(expr) == 0.5) || (cons != NULL) )
      *enforcing |= *participating;

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler
 * @todo: remember if we are in the original variables and avoid reevaluating
 */
static
SCIP_DECL_NLHDLREVALAUX(nlhdlrEvalauxSoc)
{ /*lint --e{715}*/
   int i;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->vars != NULL);
   assert(nlhdlrexprdata->transcoefs != NULL);
   assert(nlhdlrexprdata->transcoefsidx != NULL);
   assert(nlhdlrexprdata->nterms > 1);

   /* if the original expression is a norm, evaluate w.r.t. the auxiliary variables */
   if( SCIPisExprPower(scip, expr) )
   {
      assert(SCIPgetExponentExprPow(expr) == 0.5);

      updateVarVals(scip, nlhdlrexprdata, sol, FALSE);

      /* compute sum_i coef_i expr_i^2 */
      *auxvalue = 0.0;
      for( i = 0; i < nlhdlrexprdata->nterms - 1; ++i )
      {
         SCIP_Real termval;

         termval = evalSingleTerm(scip, nlhdlrexprdata, i);
         *auxvalue += SQR(termval);
      }

      assert(*auxvalue >= 0.0);

      /* compute sqrt(sum_i coef_i expr_i^2) */
      *auxvalue = sqrt(*auxvalue);
   }
   /* otherwise, evaluate the original quadratic expression w.r.t. the created auxvars of the children */
   else
   {
      SCIP_EXPR** children;
      SCIP_Real* childcoefs;
      int nchildren;

      assert(SCIPisExprSum(scip, expr));

      children = SCIPexprGetChildren(expr);
      childcoefs = SCIPgetCoefsExprSum(expr);
      nchildren = SCIPexprGetNChildren(expr);

      *auxvalue = SCIPgetConstantExprSum(expr);

      for( i = 0; i < nchildren; ++i )
      {
         if( SCIPisExprPower(scip, children[i]) )
         {
            SCIP_VAR* argauxvar;
            SCIP_Real solval;

            assert(SCIPgetExponentExprPow(children[i]) == 2.0);

            argauxvar = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(children[i])[0]);
            assert(argauxvar != NULL);

            solval = SCIPgetSolVal(scip, sol, argauxvar);
            *auxvalue += childcoefs[i] * SQR( solval );
         }
         else if( SCIPisExprProduct(scip, children[i]) )
         {
            SCIP_VAR* argauxvar1;
            SCIP_VAR* argauxvar2;

            assert(SCIPexprGetNChildren(children[i]) == 2);

            argauxvar1 = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(children[i])[0]);
            argauxvar2 = SCIPgetExprAuxVarNonlinear(SCIPexprGetChildren(children[i])[1]);
            assert(argauxvar1 != NULL);
            assert(argauxvar2 != NULL);

            *auxvalue += childcoefs[i] * SCIPgetSolVal(scip, sol, argauxvar1) * SCIPgetSolVal(scip, sol, argauxvar2);
         }
         else
         {
            SCIP_VAR* argauxvar;

            argauxvar = SCIPgetExprAuxVarNonlinear(children[i]);
            assert(argauxvar != NULL);

            *auxvalue += childcoefs[i] * SCIPgetSolVal(scip, sol, argauxvar);
         }
      }
   }

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSINITLP) */
static
SCIP_DECL_NLHDLRINITSEPA(nlhdlrInitSepaSoc)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;

   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   /* already needed for debug solution */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexprdata->varvals, nlhdlrexprdata->nvars) );

   /* if we have 3 or more terms in lhs create variable and row for disaggregation */
   if( nlhdlrexprdata->nterms > 3 )
   {
      /* create variables for cone disaggregation */
      SCIP_CALL( createDisaggrVars(scip, expr, nlhdlrexprdata) );

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real lhsval;
         SCIP_Real rhsval;
         SCIP_Real disvarval;
         int ndisvars;
         int nterms;
         int i;

         for( i = 0; i < nlhdlrexprdata->nvars; ++i )
         {
            SCIP_CALL( SCIPdebugGetSolVal(scip, SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[i]), &nlhdlrexprdata->varvals[i]) );
         }

         /*  the debug solution value of the disaggregation variables is set to
          *      (v_i^T x + beta_i)^2 / (v_{n+1}^T x + beta_{n+1})
          *  if (v_{n+1}^T x + beta_{n+1}) is different from 0.
          *  Otherwise, the debug solution value is set to 0.
          */

         nterms = nlhdlrexprdata->nterms;

         /* find value of rhs */
         rhsval = evalSingleTerm(scip, nlhdlrexprdata, nterms - 1);

         /* set value of disaggregation vars */
         ndisvars = nlhdlrexprdata->nterms - 1;

         if( SCIPisZero(scip, rhsval) )
         {
            for( i = 0; i < ndisvars; ++i )
            {
               SCIP_CALL( SCIPdebugAddSolVal(scip, nlhdlrexprdata->disvars[i], 0.0) );
            }
         }
         else
         {
            /* set value for each disaggregation variable corresponding to quadratic term */
            for( i = 0; i < ndisvars; ++i )
            {
               lhsval = evalSingleTerm(scip, nlhdlrexprdata, i);

               disvarval = SQR(lhsval) / rhsval;

               SCIP_CALL( SCIPdebugAddSolVal(scip, nlhdlrexprdata->disvars[i], disvarval) );
            }
         }
      }
#endif

      /* create the disaggregation row and store it in nlhdlrexprdata */
      SCIP_CALL( createDisaggrRow(scip, conshdlr, expr, nlhdlrexprdata) );
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "initlp for \n");
   printNlhdlrExprData(scip, nlhdlrexprdata);
#endif

   /* add some initial cuts on well-selected coordinates */
   if( nlhdlrexprdata->nterms == 2 )
   {
      /* we have |v_1^T x + \beta_1| \leq v_2^T x + \beta_2
       *
       * we should linearize at points where the first term is -1 or 1, so we can take
       *
       * x = v_1 / ||v_1||^2 (+/-1 - beta_1)
       */
      SCIP_Real plusminus1;
      SCIP_Real norm;
      int i;

      /* calculate ||v_1||^2 */
      norm = 0.0;
      for( i = nlhdlrexprdata->termbegins[0]; i < nlhdlrexprdata->termbegins[1]; ++i )
         norm += SQR(nlhdlrexprdata->transcoefs[i]);
      assert(norm > 0.0);

      BMSclearMemoryArray(nlhdlrexprdata->varvals, nlhdlrexprdata->nvars);

      for( plusminus1 = -1.0; plusminus1 <= 1.0; plusminus1 += 2.0 )
      {
         /* set x = v_1 / ||v_1||^2 (plusminus1 - beta_1) */
         for( i = nlhdlrexprdata->termbegins[0]; i < nlhdlrexprdata->termbegins[1]; ++i )
            nlhdlrexprdata->varvals[nlhdlrexprdata->transcoefsidx[i]] = nlhdlrexprdata->transcoefs[i] / norm * (plusminus1 - nlhdlrexprdata->offsets[0]);
         assert(SCIPisEQ(scip, evalSingleTerm(scip, nlhdlrexprdata, 0), plusminus1));

         /* compute gradient cut */
         SCIP_CALL( generateCutSolSOC(scip, &rowprep, expr, cons, nlhdlrexprdata, -SCIPinfinity(scip), evalSingleTerm(scip, nlhdlrexprdata, 1)) );

         if( rowprep != NULL )
         {
            SCIP_Bool success = FALSE;

            SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, NULL, SCIPgetHugeValue(scip), &success) );
            if( success )
            {
               SCIP_ROW* cut;
               SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, rowprep, cons) );
               SCIP_CALL( SCIPaddRow(scip, cut, FALSE, infeasible) );
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            }

            SCIPfreeRowprep(scip, &rowprep);
         }

         if( *infeasible )
            return SCIP_OKAY;
      }
   }
   else if( nlhdlrexprdata->nterms == 3 && nlhdlrexprdata->termbegins[0] != nlhdlrexprdata->termbegins[1] )
   {
      /* we have \sqrt{ (v_1^T x + \beta_1)^2 + (v_2^T x + \beta_2)^2 } \leq v_3^T x + \beta_3
       * with v_1 != 0
       *
       * we should linearize at points where the first and second term are (-1,0), (1,1), (1,-1), i.e.,
       *
       * v_1^T x + \beta_1 = -1   1   1
       * v_2^T x + \beta_2 =  0   1  -1
       *
       * Let i be the index of the first nonzero element in v_1.
       * Let j != i be an index of v_2, to be determined.
       * Assume all other entries of x will be set to 0.
       * Then we have
       *   (v_1)_i x_i + (v_1)_j x_j = c   with  c = -1 - beta_1
       *   (v_2)_i x_i + (v_2)_j x_j = d   with  d =  0 - beta_2
       *
       * Since (v_1)_i != 0, this gives
       *   x_i = 1/(v_1)_i (c - (v_1)_j x_j)
       * Substituting in the 2nd equation, we get
       *    (v_2)_i/(v_1)_i (c - (v_1)_j x_j) + (v_2)_j x_j = d
       * -> ((v_2)_j - (v_2)_i (v_1)_j / (v_1)_i) x_j = d - (v_2)_i/(v_1)_i c
       * Now find j such that (v_2)_j - (v_2)_i (v_1)_j / (v_1)_i != 0.
       *
       * If v_2 = 0, then linearize only for first term being -1 or 1 and don't care about value of second term.
       * We then set j arbitrary, x_i = 1/(v_1)_i c, other coordinates of x = zero.
       */
      static const SCIP_Real refpoints[3][2] = { {-1.0, 0.0}, {1.0, 1.0}, {1.0, -1.0} };
      SCIP_Real v1i, v1j = 0.0;
      SCIP_Real v2i, v2j = 0.0;
      SCIP_Bool v2zero;
      int i;
      int j = -1;
      int k;
      int pos;

      i = nlhdlrexprdata->transcoefsidx[nlhdlrexprdata->termbegins[0]];
      v1i = nlhdlrexprdata->transcoefs[nlhdlrexprdata->termbegins[0]];
      assert(v1i != 0.0);

      v2zero = nlhdlrexprdata->termbegins[1] == nlhdlrexprdata->termbegins[2];

      /* get (v_2)_i; as it is a sparse vector, we need to search for i in transcoefsidx */
      v2i = 0.0;
      if( !v2zero && SCIPsortedvecFindInt(nlhdlrexprdata->transcoefsidx + nlhdlrexprdata->termbegins[1], i, nlhdlrexprdata->termbegins[2] - nlhdlrexprdata->termbegins[1], &pos) )
      {
         assert(nlhdlrexprdata->transcoefsidx[nlhdlrexprdata->termbegins[1] + pos] == i);
         v2i = nlhdlrexprdata->transcoefs[nlhdlrexprdata->termbegins[1] + pos];
      }

      /* find a j, for now look only into indices with (v_2)_j != 0 */
      for( k = nlhdlrexprdata->termbegins[1]; k < nlhdlrexprdata->termbegins[2]; ++k )
      {
         /* check whether transcoefsidx[k] could be a good j */

         if( nlhdlrexprdata->transcoefsidx[k] == i ) /* i == j is not good */
            continue;

         /* get (v_1)_j; as it is a sparse vector, we need to search for j in transcoefsidx */
         v1j = 0.0;
         if( SCIPsortedvecFindInt(nlhdlrexprdata->transcoefsidx + nlhdlrexprdata->termbegins[0], nlhdlrexprdata->transcoefsidx[k], nlhdlrexprdata->termbegins[1] - nlhdlrexprdata->termbegins[0], &pos) )
         {
            assert(nlhdlrexprdata->transcoefsidx[nlhdlrexprdata->termbegins[0] + pos] == nlhdlrexprdata->transcoefsidx[k]);
            v1j = nlhdlrexprdata->transcoefs[nlhdlrexprdata->termbegins[0] + pos];
         }

         v2j = nlhdlrexprdata->transcoefs[k];

         if( SCIPisZero(scip, v2j - v2i * v1j / v1i) )  /* (v_2)_j - (v_2)_i (v_1)_j / (v_1)_i = 0 is also not good */
            continue;

         j = nlhdlrexprdata->transcoefsidx[k];
         break;
      }

      if( v2zero )
      {
         j = 0;
         v1j = 0.0;
         v2j = 0.0;
      }

      if( j != -1 )
      {
         SCIP_Real c, d;
         int point;

         BMSclearMemoryArray(nlhdlrexprdata->varvals, nlhdlrexprdata->nvars);

         for( point = 0; point < (v2zero ? 2 : 3); ++point )
         {
            c = refpoints[point][0] - nlhdlrexprdata->offsets[0];

            if( !v2zero )
            {
               /* set x_j and x_i */
               d = refpoints[point][1] - nlhdlrexprdata->offsets[1];
               nlhdlrexprdata->varvals[j] = (d - v2i/v1i*c) / (v2j - v2i * v1j / v1i);
               nlhdlrexprdata->varvals[i] = (c - v1j * nlhdlrexprdata->varvals[j]) / v1i;

               SCIPdebugMsg(scip, "<%s>(%d) = %g, <%s>(%d) = %g\n",
                  SCIPvarGetName(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[i])), i, nlhdlrexprdata->varvals[i],
                  SCIPvarGetName(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[j])), j, nlhdlrexprdata->varvals[j]);
            }
            else
            {
               /* set x_i */
               nlhdlrexprdata->varvals[i] = c / v1i;

               SCIPdebugMsg(scip, "<%s>(%d) = %g\n",
                  SCIPvarGetName(SCIPgetExprAuxVarNonlinear(nlhdlrexprdata->vars[i])), i, nlhdlrexprdata->varvals[i]);
            }

            assert(SCIPisEQ(scip, evalSingleTerm(scip, nlhdlrexprdata, 0), refpoints[point][0]));
            assert(v2zero || SCIPisEQ(scip, evalSingleTerm(scip, nlhdlrexprdata, 1), refpoints[point][1]));

            /* compute gradient cut */
            SCIP_CALL( generateCutSolSOC(scip, &rowprep, expr, cons, nlhdlrexprdata, -SCIPinfinity(scip), evalSingleTerm(scip, nlhdlrexprdata, 2)) );

            if( rowprep != NULL )
            {
               SCIP_Bool success = FALSE;

               SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, NULL, SCIPgetHugeValue(scip), &success) );
               if( success )
               {
                  SCIP_ROW* cut;
                  SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, rowprep, cons) );
                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, infeasible) );
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }

               SCIPfreeRowprep(scip, &rowprep);
            }

            if( *infeasible )
               return SCIP_OKAY;
         }
      }
   }
   else if( nlhdlrexprdata->nterms == 3 )
   {
      /* we have \sqrt{ \beta_1^2 + (v_2^T x + \beta_2)^2 } \leq v_3^T x + \beta_3
       * with v_2 != 0
       *
       * we should linearize at points where the second term is -1 or 1
       *
       * set x = v_2 / ||v_2||^2 (+/-1 - beta_2)
       */
      SCIP_Real plusminus1;
      SCIP_Real norm;
      int i;

      /* calculate ||v_2||^2 */
      norm = 0.0;
      for( i = nlhdlrexprdata->termbegins[1]; i < nlhdlrexprdata->termbegins[2]; ++i )
         norm += SQR(nlhdlrexprdata->transcoefs[i]);
      assert(norm > 0.0);

      BMSclearMemoryArray(nlhdlrexprdata->varvals, nlhdlrexprdata->nvars);

      for( plusminus1 = -1.0; plusminus1 <= 1.0; plusminus1 += 2.0 )
      {
         /* set x = v_2 / ||v_2||^2 (plusminus1 - beta_2) */
         for( i = nlhdlrexprdata->termbegins[1]; i < nlhdlrexprdata->termbegins[2]; ++i )
            nlhdlrexprdata->varvals[nlhdlrexprdata->transcoefsidx[i]] = nlhdlrexprdata->transcoefs[i] / norm * (plusminus1 - nlhdlrexprdata->offsets[1]);
         assert(SCIPisEQ(scip, evalSingleTerm(scip, nlhdlrexprdata, 1), plusminus1));

         /* compute gradient cut */
         SCIP_CALL( generateCutSolSOC(scip, &rowprep, expr, cons, nlhdlrexprdata, -SCIPinfinity(scip), evalSingleTerm(scip, nlhdlrexprdata, 2)) );

         if( rowprep != NULL )
         {
            SCIP_Bool success = FALSE;

            SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, NULL, SCIPgetHugeValue(scip), &success) );
            if( success )
            {
               SCIP_ROW* cut;
               SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, rowprep, cons) );
               SCIP_CALL( SCIPaddRow(scip, cut, FALSE, infeasible) );
               SCIP_CALL( SCIPreleaseRow(scip, &cut) );
            }

            SCIPfreeRowprep(scip, &rowprep);
         }

         if( *infeasible )
            return SCIP_OKAY;
      }
   }
   else
   {
      /* generate gradient cuts for the small rotated cones
       *  \f[
       *    \sqrt{4(v_k^T x + \beta_k)^2 + (v_n^T x + \beta_n - y_k)^2} - v_n^T x - \beta_n - y_k \leq 0.
       *  \f]
       *
       * we should linearize again at points where the first and second term (inside sqr) are (-1/2,0), (1/2,1), (1/2,-1).
       * Since we have y_k, we can achieve this more easily here via
       *   x = v_k/||v_k||^2 (+/-0.5 - beta_k)
       *   y_k = v_n^T x + beta_n + 0/1/-1
       *
       * If v_k = 0, then we use x = 0 and linearize for second term being 1 and -1 only
       */
      static const SCIP_Real refpoints[3][2] = { {-0.5, 0.0}, {0.5, 1.0}, {0.5, -1.0} };
      SCIP_Real rhsval;
      SCIP_Real norm;
      int point;
      int i;
      int k;
      SCIP_Bool vkzero;

      /* add disaggregation row to LP */
      SCIP_CALL( SCIPaddRow(scip, nlhdlrexprdata->disrow, FALSE, infeasible) );

      if( *infeasible )
         return SCIP_OKAY;

      for( k = 0; k < nlhdlrexprdata->nterms - 1; ++k )
      {
         vkzero = nlhdlrexprdata->termbegins[k+1] == nlhdlrexprdata->termbegins[k];
         assert(!vkzero || nlhdlrexprdata->offsets[k] != 0.0);

         /* calculate ||v_k||^2 */
         norm = 0.0;
         for( i = nlhdlrexprdata->termbegins[k]; i < nlhdlrexprdata->termbegins[k+1]; ++i )
            norm += SQR(nlhdlrexprdata->transcoefs[i]);
         assert(vkzero || norm > 0.0);

         BMSclearMemoryArray(nlhdlrexprdata->varvals, nlhdlrexprdata->nvars);

         for( point = vkzero ? 1 : 0; point < 3; ++point )
         {
            /* set x = v_k / ||v_k||^2 (refpoints[point][0] - beta_k) / 2 */
            for( i = nlhdlrexprdata->termbegins[k]; i < nlhdlrexprdata->termbegins[k+1]; ++i )
               nlhdlrexprdata->varvals[nlhdlrexprdata->transcoefsidx[i]] = nlhdlrexprdata->transcoefs[i] / norm * (refpoints[point][0] - nlhdlrexprdata->offsets[k]);  /*lint !e795*/
            assert(vkzero || SCIPisEQ(scip, evalSingleTerm(scip, nlhdlrexprdata, k), refpoints[point][0]));

            /* set y_k = v_n^T x + beta_n + 0/1/-1 */
            rhsval = evalSingleTerm(scip, nlhdlrexprdata, nlhdlrexprdata->nterms - 1);
            nlhdlrexprdata->disvarvals[k] = rhsval + refpoints[point][1];

            /* compute gradient cut */
            SCIP_CALL( generateCutSolDisagg(scip, &rowprep, expr, cons, nlhdlrexprdata, k, -SCIPinfinity(scip), rhsval) );

            if( rowprep != NULL )
            {
               SCIP_Bool success = FALSE;

               SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, NULL, SCIPgetHugeValue(scip), &success) );
               if( success )
               {
                  SCIP_ROW* cut;
                  SCIP_CALL( SCIPgetRowprepRowCons(scip, &cut, rowprep, cons) );
                  SCIP_CALL( SCIPaddRow(scip, cut, FALSE, infeasible) );
                  SCIP_CALL( SCIPreleaseRow(scip, &cut) );
               }

               SCIPfreeRowprep(scip, &rowprep);
            }

            if( *infeasible )
               return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
static
SCIP_DECL_NLHDLREXITSEPA(nlhdlrExitSepaSoc)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   /* free disaggreagation row */
   if( nlhdlrexprdata->disrow != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &nlhdlrexprdata->disrow) );
   }

   SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->varvals, nlhdlrexprdata->nvars);

   return SCIP_OKAY;
}


/** nonlinear handler separation callback */
static
SCIP_DECL_NLHDLRENFO(nlhdlrEnfoSoc)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_Real rhsval;
   int ndisaggrs;
   int k;
   SCIP_Bool infeasible;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nterms < 4 || nlhdlrexprdata->disrow != NULL);
   assert(nlhdlrexprdata->nterms > 1);

   *result = SCIP_DIDNOTFIND;

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* update varvals
    * set variables close to integer to integer, in particular when close to zero
    * for simple soc's (no large v_i, no offsets), variables close to zero would give coefficients close to zero in the cut,
    * which the cut cleanup may have problems to relax (and we end up with local or much relaxed cuts)
    * also when close to other integers, rounding now may prevent some relaxation in cut cleanup
    */
   updateVarVals(scip, nlhdlrexprdata, sol, TRUE);

   rhsval = evalSingleTerm(scip, nlhdlrexprdata, nlhdlrexprdata->nterms - 1);

   /* if there are three or two terms just compute gradient cut */
   if( nlhdlrexprdata->nterms < 4 )
   {
      SCIP_ROWPREP* rowprep;

      /* compute gradient cut */
      SCIP_CALL( generateCutSolSOC(scip, &rowprep, expr, cons, nlhdlrexprdata, SCIPgetLPFeastol(scip), rhsval) );

      if( rowprep != NULL )
      {
         SCIP_CALL( addCut(scip, nlhdlrdata, rowprep, sol, cons, allowweakcuts, result) );

         SCIPfreeRowprep(scip, &rowprep);
      }
      else
      {
         SCIPdebugMsg(scip, "failed to generate cut for SOC\n");
      }

      return SCIP_OKAY;
   }

   ndisaggrs = nlhdlrexprdata->nterms - 1;

   /* check whether the aggregation row is in the LP */
   if( !SCIProwIsInLP(nlhdlrexprdata->disrow) && -SCIPgetRowSolFeasibility(scip, nlhdlrexprdata->disrow, sol) > SCIPgetLPFeastol(scip) )
   {
      SCIP_CALL( SCIPaddRow(scip, nlhdlrexprdata->disrow, TRUE, &infeasible) );
      SCIPdebugMsg(scip, "added disaggregation row to LP, cutoff=%u\n", infeasible);

      if( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      *result = SCIP_SEPARATED;
   }

   for( k = 0; k < ndisaggrs && *result != SCIP_CUTOFF; ++k )
   {
      SCIP_ROWPREP* rowprep;

      /* compute gradient cut */
      SCIP_CALL( generateCutSolDisagg(scip, &rowprep, expr, cons, nlhdlrexprdata, k, SCIPgetLPFeastol(scip), rhsval) );

      if( rowprep != NULL )
      {
         SCIP_CALL( addCut(scip, nlhdlrdata, rowprep, sol, cons, allowweakcuts, result) );

         SCIPfreeRowprep(scip, &rowprep);
      }
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_NLHDLRSOLLINEARIZE(nlhdlrSollinearizeSoc)
{ /*lint --e{715}*/
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_Real rhsval;
   int k;

   assert(sol != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nterms < 4 || nlhdlrexprdata->disrow != NULL);
   assert(nlhdlrexprdata->nterms > 1);

   nlhdlrdata = SCIPnlhdlrGetData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* update varvals */
   updateVarVals(scip, nlhdlrexprdata, sol, TRUE);

   rhsval = evalSingleTerm(scip, nlhdlrexprdata, nlhdlrexprdata->nterms - 1);

   /* if there are three or two terms just compute gradient cut */
   if( nlhdlrexprdata->nterms < 4 )
   {
      SCIP_ROWPREP* rowprep;

      /* compute gradient cut */
      SCIP_CALL( generateCutSolSOC(scip, &rowprep, expr, cons, nlhdlrexprdata, -SCIPinfinity(scip), rhsval) );

      if( rowprep != NULL )
      {
         SCIP_CALL( addCutPool(scip, nlhdlrdata, rowprep, sol, cons) );

         SCIPfreeRowprep(scip, &rowprep);
      }

      return SCIP_OKAY;
   }

   for( k = 0; k < nlhdlrexprdata->nterms - 1; ++k )
   {
      SCIP_ROWPREP* rowprep;

      /* compute gradient cut */
      SCIP_CALL( generateCutSolDisagg(scip, &rowprep, expr, cons, nlhdlrexprdata, k, -SCIPinfinity(scip), rhsval) );

      if( rowprep != NULL )
      {
         SCIP_CALL( addCutPool(scip, nlhdlrdata, rowprep, sol, cons) );

         SCIPfreeRowprep(scip, &rowprep);
      }
   }

   return SCIP_OKAY;
}

/*
 * nonlinear handler specific interface methods
 */

/** includes SOC nonlinear handler in nonlinear constraint handler */
SCIP_RETCODE SCIPincludeNlhdlrSoc(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_NLHDLRDATA* nlhdlrdata;
   SCIP_NLHDLR* nlhdlr;

   assert(scip != NULL);

   /* create nonlinear handler data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &nlhdlrdata) );

   SCIP_CALL( SCIPincludeNlhdlrNonlinear(scip, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY, NLHDLR_ENFOPRIORITY, nlhdlrDetectSoc, nlhdlrEvalauxSoc, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPnlhdlrSetCopyHdlr(nlhdlr, nlhdlrCopyhdlrSoc);
   SCIPnlhdlrSetFreeHdlrData(nlhdlr, nlhdlrFreehdlrdataSoc);
   SCIPnlhdlrSetFreeExprData(nlhdlr, nlhdlrFreeExprDataSoc);
   SCIPnlhdlrSetInitExit(nlhdlr, nlhdlrInitSoc, nlhdlrExitSoc);
   SCIPnlhdlrSetSepa(nlhdlr, nlhdlrInitSepaSoc, nlhdlrEnfoSoc, NULL, nlhdlrExitSepaSoc);
   SCIPnlhdlrSetSollinearize(nlhdlr, nlhdlrSollinearizeSoc);

   /* add soc nlhdlr parameters */
   /* TODO should we get rid of this and use separating/mineffiacy(root) instead, which is 1e-4? */
   SCIP_CALL( SCIPaddRealParam(scip, "nlhdlr/" NLHDLR_NAME "/mincutefficacy",
         "Minimum efficacy which a cut needs in order to be added.",
         &nlhdlrdata->mincutefficacy, FALSE, DEFAULT_MINCUTEFFICACY, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "nlhdlr/" NLHDLR_NAME "/compeigenvalues",
         "Should Eigenvalue computations be done to detect complex cases in quadratic constraints?",
         &nlhdlrdata->compeigenvalues, FALSE, DEFAULT_COMPEIGENVALUES, NULL, NULL) );

   return SCIP_OKAY;
}

/** checks whether constraint is SOC representable in original variables and returns the SOC representation
 *
 * The SOC representation has the form:
 * \f$\sqrt{\sum_{i=1}^{n} (v_i^T x + \beta_i)^2} - v_{n+1}^T x - \beta_{n+1} \lessgtr 0\f$,
 * where \f$n+1 = \text{nterms}\f$ and the inequality type is given by sidetype (`SCIP_SIDETYPE_RIGHT` if inequality
 * is \f$\leq\f$, `SCIP_SIDETYPE_LEFT` if \f$\geq\f$).
 *
 * For each term (i.e. for each \f$i\f$ in the above notation as well as \f$n+1\f$), the constant \f$\beta_i\f$ is given by the
 * corresponding element `offsets[i-1]` and `termbegins[i-1]` is the starting position of the term in arrays
 * `transcoefs` and `transcoefsidx`. The overall number of nonzeros is `termbegins[nterms]`.
 *
 * Arrays `transcoefs` and `transcoefsidx` have size `termbegins[nterms]` and define the linear expressions \f$v_i^T x\f$
 * for each term. For a term \f$i\f$ in the above notation, the nonzeroes are given by elements
 * `termbegins[i-1]...termbegins[i]` of `transcoefs` and `transcoefsidx`. There may be no nonzeroes for some term (i.e.,
 * constant terms are possible). `transcoefs` contains the coefficients \f$v_i\f$ and `transcoefsidx` contains positions of
 * variables in the `vars` array.
 *
 * The `vars` array has size `nvars` and contains \f$x\f$ variables; each variable is included at most once.
 *
 * The arrays should be freed by calling SCIPfreeSOCArraysNonlinear().
 *
 * This function uses the methods that are used in the detection algorithm of the SOC nonlinear handler.
 */
SCIP_RETCODE SCIPisSOCNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< nonlinear constraint */
   SCIP_Bool             compeigenvalues,    /**< whether eigenvalues should be computed to detect complex cases */
   SCIP_Bool*            success,            /**< pointer to store whether SOC structure has been detected */
   SCIP_SIDETYPE*        sidetype,           /**< pointer to store which side of cons is SOC representable; only
                                              *   valid when success is TRUE */
   SCIP_VAR***           vars,               /**< variables (x) that appear on both sides; no duplicates are allowed */
   SCIP_Real**           offsets,            /**< offsets of both sides (beta_i) */
   SCIP_Real**           transcoefs,         /**< non-zeros of linear transformation vectors (v_i) */
   int**                 transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int**                 termbegins,         /**< starting indices of transcoefs for each term */
   int*                  nvars,              /**< total number of variables appearing (i.e. size of vars) */
   int*                  nterms              /**< number of summands in the SQRT +1 for RHS (n+1) */
   )
{
   SCIP_NLHDLRDATA nlhdlrdata;
   SCIP_NLHDLREXPRDATA *nlhdlrexprdata;
   SCIP_Real conslhs;
   SCIP_Real consrhs;
   SCIP_EXPR* expr;
   SCIP_Bool enforcebelow;
   int i;

   assert(cons != NULL);

   expr = SCIPgetExprNonlinear(cons);
   assert(expr != NULL);

   nlhdlrdata.mincutefficacy = 0.0;
   nlhdlrdata.compeigenvalues = compeigenvalues;

   conslhs = SCIPgetLhsNonlinear(cons);
   consrhs = SCIPgetRhsNonlinear(cons);

   SCIP_CALL( detectSOC(scip, &nlhdlrdata, expr, conslhs, consrhs, &nlhdlrexprdata, &enforcebelow, success) );

   /* the constraint must be SOC representable in original variables */
   if( *success )
   {
      assert(nlhdlrexprdata != NULL);

      for( i = 0; i < nlhdlrexprdata->nvars; ++i )
      {
         if( !SCIPisExprVar(scip, nlhdlrexprdata->vars[i]) )
         {
            *success = FALSE;
            break;
         }
      }
   }

   if( *success )
   {
      *sidetype = enforcebelow ? SCIP_SIDETYPE_RIGHT : SCIP_SIDETYPE_LEFT;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, vars, nlhdlrexprdata->nvars) );

      for( i = 0; i < nlhdlrexprdata->nvars; ++i )
      {
         (*vars)[i] = SCIPgetVarExprVar(nlhdlrexprdata->vars[i]);
         assert((*vars)[i] != NULL);
      }
      SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->vars, nlhdlrexprdata->nvars);
      *offsets = nlhdlrexprdata->offsets;
      *transcoefs = nlhdlrexprdata->transcoefs;
      *transcoefsidx = nlhdlrexprdata->transcoefsidx;
      *termbegins = nlhdlrexprdata->termbegins;
      *nvars = nlhdlrexprdata->nvars;
      *nterms = nlhdlrexprdata->nterms;
      SCIPfreeBlockMemory(scip, &nlhdlrexprdata);
   }
   else
   {
      if( nlhdlrexprdata != NULL )
      {
         SCIP_CALL( freeNlhdlrExprData(scip, &nlhdlrexprdata) );
      }
      *vars = NULL;
      *offsets = NULL;
      *transcoefs = NULL;
      *transcoefsidx = NULL;
      *termbegins = NULL;
      *nvars = 0;
      *nterms = 0;
   }

   return SCIP_OKAY;
}

/** frees arrays created by SCIPisSOCNonlinear() */
void SCIPfreeSOCArraysNonlinear(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR***           vars,               /**< variables that appear on both sides (x) */
   SCIP_Real**           offsets,            /**< offsets of both sides (beta_i) */
   SCIP_Real**           transcoefs,         /**< non-zeros of linear transformation vectors (v_i) */
   int**                 transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int**                 termbegins,         /**< starting indices of transcoefs for each term */
   int                   nvars,              /**< total number of variables appearing */
   int                   nterms              /**< number of summands in the SQRT +1 for RHS (n+1) */
   )
{
   int ntranscoefs;

   if( nvars == 0 )
      return;

   assert(vars != NULL);
   assert(offsets != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(termbegins != NULL);

   ntranscoefs = (*termbegins)[nterms];

   SCIPfreeBlockMemoryArray(scip, termbegins, nterms + 1);
   SCIPfreeBlockMemoryArray(scip, transcoefsidx, ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, transcoefs, ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, offsets, nterms);
   SCIPfreeBlockMemoryArray(scip, vars, nvars);
}
