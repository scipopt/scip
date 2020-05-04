/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_soc.c
 * @brief  nonlinear handler for second order cone constraints

 * @author Benjamin Mueller
 * @author Fabian Wegscheider
 *
 * This is a nonlinear handler for second order cone constraints of the form
 *
 * \f$\sqrt{\gamma + \sum_{i=1}^{n} (v_i^T x + \beta_i)^2} \leq v_{n+1}^T x + \beta_{n+1}\f$,
 *
 * where \f$\gamma \geq 0\f$ and \f$v_{n+1}^T x +\beta_{n+1} \geq 0\f$.
 *
 * @TODO: Don't disaggregate if n = 1. This requires a separate separation method for this case.
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_soc.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_var.h"
#include "scip/debug.h"
#include "nlpi/nlpi_ipopt.h"


/* fundamental nonlinear handler properties */
#define NLHDLR_NAME           "soc"
#define NLHDLR_DESC           "soc nonlinear handler"
#define NLHDLR_PRIORITY             100
#define DEFAULT_MINCUTEFFICACY     1e-5 /** default value for parameter mincutefficacy */
#define DEFAULT_ENFOFREQ              5 /** default value for parameter enfofreq */
#define DEFAULT_MAXENFOROUNDSROOT    -1 /** default value for parameter maxenforoundsroot */
#define DEFAULT_MAXENFOROUNDS         1 /** default value for parameter maxenforounds */
#define DEFAULT_COMPEIGENVALUES    TRUE /** default value for parameter compeigenvalues */

/*
 * Data structures
 */

/** nonlinear handler expression data. The data is structured in the following way:
 *
 *  A 'term' is one of the arguments of the  quadratic terms, i.e. v_i^T x + beta_i.
 *  The last term is always the one on the right-hand side. This means that nterms is
 *  equal to n+1 in the above description.
 *
 *  - vars contains a list of all variables that appear in the expression (no duplicates)
 *  - offsets contains the constants beta_i of each term
 *  - transcoefs contains the non-zero values of the transformation vectors v_i of each term
 *  - transcoefsidx contains for each entry of transcoefs the position of the respective variable in vars
 *  - termbegins contains the index at which the transcoefs of each term start
 *  - nnonzeroes contains the number of non-zeroes in v_i of each term
 *  - constant is gamma which is the constant inside the square root
 *  - nvars is the total number of unique variables appearing (length of vars)
 *  - nterms is the total number of terms appearing on both sides
 *  - ntranscoefs is the total number of entries in transcoefs and transcoefsidx
 *
 *  The disaggregation is implicitly stored in the variables disvars and disrow. An SOC as
 *  described above is replaced by n+1 (or n if gamma = 0) smaller SOCs
 *
 *              (v_i^T x + beta_i)^2 <= disvar_i     * (v_{n+1}^T x + beta_{n+1})
 *                             gamma <= disvar_{n+1} * (v_{n+1}^T x + beta_{n+1})
 *
 *  and the row       sum_i disvar_i <= v_{n+1}^T x + beta_{n+1}.
 *
 *
 *  Example: The constraint SQRT(5 + (3x - 4y + 2)^2 + y^2 + 7z^2) <= 5x - y - 1
 *           results in the following nlhdlrexprdata:
 *
 *           vars = {x, y, z}
 *           offsets = {2, 0, 0, -1}
 *           transcoefs = {3, -4, 1, SQRT(7), 5, -1}
 *           transcoefsidx = {0, 1, 1, 2, 0, 1}
 *           termbegins = {0, 2, 3, 4}
 *           nnonzeroes = {2, 1, 1, 2}
 *           constant = 5
 *           nvars = 3
 *           nterms = 4
 *           ntranscoefs = 6
 */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_VAR**            vars;               /**< variables appearing on both sides (x) */
   SCIP_Real*            offsets;            /**< offsets of both sides (beta_i) */
   SCIP_Real*            transcoefs;         /**< non-zeroes of linear transformation vectors (v_i) */
   int*                  transcoefsidx;      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  termbegins;         /**< starting indices of transcoefs for each term */
   int*                  nnonzeroes;         /**< number of non-zeroes in each v_i */
   SCIP_Real             constant;           /**< constant on left hand side (gamma) */
   int                   nvars;              /**< total number of variables appearing */
   int                   nterms;             /**< number of summands in the SQRT (excluding gamma) +1 for RHS (n+1) */
   int                   ntranscoefs;        /**< total number of entries in transcoefs */

   /* variables for cone disaggregation */
   SCIP_VAR**            disvars;            /**< disaggregation variables for each expression; the last entry
                                               *  corresponds to the constant term */
   SCIP_ROW*             disrow;             /**< disaggregation row */
};

struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_NODE*            prevnode;           /**< the node for which enforcement was last called */
   int                   nenfocalls;         /**< number of enforcement calls for the previous node */
   SCIP_Real             mincutefficacy;     /**< minimum efficacy a cut need to be added */
   int                   enfofreq;           /**< frequency of enforcement rounds (every x levels of depth) */
   int                   maxenforoundsroot;  /**< maximum number of enforcement rounds in the root round */
   int                   maxenforounds;      /**< maximum number of enforcement rounds in non-root rounds */
   SCIP_Bool             compeigenvalues;    /**< whether Eigenvalue computations should be done to detect complex cases */
};

/*
 * Local methods
 */

#ifdef SCIP_DEBUG
/** prints the  nlhdlr expression data */
static
void printNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata /**< pointer to store nonlinear handler expression data */
   )
{
   int nterms;
   int i;
   int j;

   nterms = nlhdlrexprdata->nterms;

   SCIPinfoMessage(scip, NULL, "SQRT( ", nlhdlrexprdata->constant);
   if( nlhdlrexprdata->constant != 0.0 )
      SCIPinfoMessage(scip, NULL, "%f + ", nlhdlrexprdata->constant);

   for( i = 0; i < nterms - 1; ++i )
   {
      int startidx;

      SCIPinfoMessage(scip, NULL, "(");

      startidx = nlhdlrexprdata->termbegins[i];

      for( j = startidx; j < startidx + nlhdlrexprdata->nnonzeroes[i]; ++j )
      {
         if( nlhdlrexprdata->transcoefs[j] != 1.0 )
            SCIPinfoMessage(scip, NULL, "%f*", nlhdlrexprdata->transcoefs[j]);
         SCIPinfoMessage(scip, NULL, "%s", SCIPvarGetName(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]));

         if( j < startidx + nlhdlrexprdata->nnonzeroes[i] - 1 )
            SCIPinfoMessage(scip, NULL, " + ");
         else if( nlhdlrexprdata->offsets[i] != 0.0 )
            SCIPinfoMessage(scip, NULL, " + %f", nlhdlrexprdata->offsets[i]);
      }

      SCIPinfoMessage(scip, NULL, ")^2");

      if( i < nterms - 2 )
         SCIPinfoMessage(scip, NULL, " + ");
   }

   SCIPinfoMessage(scip, NULL, " ) <= ");

   for( j = nlhdlrexprdata->termbegins[nterms-1]; j < nlhdlrexprdata->ntranscoefs; ++j )
   {
      if( nlhdlrexprdata->transcoefs[j] != 1.0 )
         SCIPinfoMessage(scip, NULL, "%f*", nlhdlrexprdata->transcoefs[j]);
      SCIPinfoMessage(scip, NULL, "%s", SCIPvarGetName(nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[j]]));

      if( j < nlhdlrexprdata->ntranscoefs - 1 )
         SCIPinfoMessage(scip, NULL, " + ");
      else if( nlhdlrexprdata->offsets[nterms-1] != 0.0 )
         SCIPinfoMessage(scip, NULL, " + %f", nlhdlrexprdata->offsets[nterms-1]);
   }

   SCIPinfoMessage(scip, NULL, "\n");
}
#endif

/** helper method to create variables for the cone disaggregation */
static
SCIP_RETCODE createDisaggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata  /**< nonlinear handler expression data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int ndisvars;
   int nterms;
   int i;

   assert(nlhdlrexprdata != NULL);

   nterms = nlhdlrexprdata->nterms;

   /* check whether constant has a separate entry */
   ndisvars = SCIPisZero(scip, nlhdlrexprdata->constant) ? nterms - 1 : nterms;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexprdata->disvars, ndisvars) );

   /* create disaggregation variables representing the epigraph of (v_i^T x + beta_i)^2 / (v_{n+1}^T x + beta_{n+1}) */
   for( i = 0; i < nterms - 1; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_%d", (void*) expr, i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &nlhdlrexprdata->disvars[i], name, 0.0, SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[i]) );

      SCIPvarMarkRelaxationOnly(nlhdlrexprdata->disvars[i]);
      SCIP_CALL( SCIPaddVarLocksType(scip, nlhdlrexprdata->disvars[i], SCIP_LOCKTYPE_MODEL, 1, 1) );
   }

   /* create disaggregation variable representing the epigraph of gamma / (v_{n+1}^T x + beta_{n+1}) */
   if( !SCIPisZero(scip, nlhdlrexprdata->constant) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_const", (void*) expr);
      SCIP_CALL( SCIPcreateVarBasic(scip, &nlhdlrexprdata->disvars[ndisvars - 1], name, 0.0, SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[ndisvars - 1]) );

      SCIPvarMarkRelaxationOnly(nlhdlrexprdata->disvars[ndisvars - 1]);
      SCIP_CALL( SCIPaddVarLocksType(scip, nlhdlrexprdata->disvars[ndisvars - 1], SCIP_LOCKTYPE_MODEL, 1, 1) );
   }

   return SCIP_OKAY;
}

/** helper method to free variables for the cone disaggregation */
static
SCIP_RETCODE freeDisaggrVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata  /**< nonlinear handler expression data */
   )
{
   int i;
   int size;

   assert(nlhdlrexprdata != NULL);

   /* check whether constant has a separate entry */
   size = SCIPisZero(scip, nlhdlrexprdata->constant) ? nlhdlrexprdata->nterms-1 : nlhdlrexprdata->nterms;

   /* release variables */
   for( i = 0; i < size; ++i )
   {
      SCIP_CALL( SCIPaddVarLocksType(scip, nlhdlrexprdata->disvars[i], SCIP_LOCKTYPE_MODEL, -1, -1) );
      SCIP_CALL( SCIPreleaseVar(scip, &nlhdlrexprdata->disvars[i]) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &nlhdlrexprdata->disvars, size);

   return SCIP_OKAY;
}

/** helper method to create the disaggregation row */
static
SCIP_RETCODE createDisaggrRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata  /**< nonlinear handler expression data */
   )
{
   SCIP_Real disrowrhs;
   char name[SCIP_MAXSTRLEN];
   int ndisvars;
   int nrhsvars;
   int nterms;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->disrow == NULL);

   nterms = nlhdlrexprdata->nterms;
   nrhsvars = nlhdlrexprdata->nnonzeroes[nterms-1];
   disrowrhs = nlhdlrexprdata->offsets[nterms - 1];

   /* check whether constant has a separate entry */
   ndisvars = SCIPisZero(scip, nlhdlrexprdata->constant) ? nterms - 1 : nterms;

   /* create row */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_row", (void*) expr);
   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &nlhdlrexprdata->disrow, conshdlr, name,
         -SCIPinfinity(scip), disrowrhs, FALSE, FALSE, TRUE) );

   /* add disvars to row */
   for( i = 0; i < ndisvars; ++i )
   {
      SCIP_CALL( SCIPaddVarToRow(scip, nlhdlrexprdata->disrow, nlhdlrexprdata->disvars[i], 1.0) );
   }

   /* add rhs vars to row */
   for( i = nlhdlrexprdata->ntranscoefs - nrhsvars; i < nlhdlrexprdata->ntranscoefs; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real coef;

      var = nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[i]];
      coef = -nlhdlrexprdata->transcoefs[i];

      SCIP_CALL( SCIPaddVarToRow(scip, nlhdlrexprdata->disrow, var, coef) );
   }

   return SCIP_OKAY;
}

/** helper method to create nonlinear handler expression data */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables appearing ob both sides (x) */
   SCIP_Real*            offsets,            /**< offsets of bot sides (beta_i) */
   SCIP_Real*            transcoefs,         /**< non-zeroes of linear transformation vectors (v_i) */
   int*                  transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  termbegins,         /**< starting indices of transcoefs for each term */
   int*                  nnonzeroes,         /**< number of non-zeroes in each v_i */
   SCIP_Real             constant,           /**< constant on left hand side (gamma) */
   int                   nvars,              /**< total number of variables appearing */
   int                   nterms,             /**< number of summands in the SQRT (excluding gamma) +1 for RHS */
   int                   ntranscoefs,        /**< total number of entries in transcoefs */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata /**< pointer to store nonlinear handler expression data */
   )
{
   int i;

   assert(vars != NULL);
   assert(offsets != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(termbegins != NULL);
   assert(nnonzeroes != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(constant >= 0.0);

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexprdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, vars, nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->offsets, offsets, nterms) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefs, transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefsidx, transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->termbegins, termbegins, nterms) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->nnonzeroes, nnonzeroes, nterms) );
   (*nlhdlrexprdata)->constant = constant;
   (*nlhdlrexprdata)->nvars = nvars;
   (*nlhdlrexprdata)->nterms = nterms;
   (*nlhdlrexprdata)->ntranscoefs = ntranscoefs;

   /* capture variables on LHS */
   for( i = 0; i < nvars; ++i )
   {
      assert(vars[i] != NULL);

      SCIP_CALL( SCIPcaptureVar(scip, vars[i]) );
   }

   (*nlhdlrexprdata)->disrow = NULL;
   (*nlhdlrexprdata)->disvars = NULL;

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "created nlhdlr data for the following soc expression:\n");
      printNlhdlrExprData(scip, *nlhdlrexprdata);
#endif

   return SCIP_OKAY;
}

/** helper method to free nonlinear handler expression data */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata /**< pointer to free nonlinear handler expression data */
   )
{
   int i;

   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   /* free variables and row for cone disaggregation */
   SCIP_CALL( freeDisaggrVars(scip, *nlhdlrexprdata) );

   /* release LHS variables */
   for( i = 0; i < (*nlhdlrexprdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*nlhdlrexprdata)->vars[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->nnonzeroes, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->termbegins, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefsidx, (*nlhdlrexprdata)->ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefs, (*nlhdlrexprdata)->ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->offsets, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, (*nlhdlrexprdata)->nvars);
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** evaluate a single term of the form \f$v_i^T x + \beta_i\f$ */
static
SCIP_Real evalSingleTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nonlinear handler expression data */
   SCIP_SOL*             sol,                /**< solution */
   int                   k                   /**< term to be evaluated */
   )
{
   SCIP_Real result;
   int termstart;
   int i;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(0 <= k && k < nlhdlrexprdata->nterms);

   termstart = nlhdlrexprdata->termbegins[k];
   result = nlhdlrexprdata->offsets[k];

   for( i = 0; i < nlhdlrexprdata->nnonzeroes[k]; ++i )
   {
      SCIP_Real varval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[termstart + i]]);
      result += nlhdlrexprdata->transcoefs[termstart + i] * varval;
   }

   return result;
}

/** helper method to compute and add a gradient cut for the k-th cone disaggregation
 *
 *  After the soc constraint \f$\sqrt{\gamma + \sum_i^n (v_i^T x + \beta_i)^2} \leq v_{n+1}^T x + \beta_{n+1}\f$
 *  is disggregated into the row \f$\sum_i^{n+1} y_i \leq v_{n+1}^T x + \beta_{n+1}\f$ and the smaller soc constraints
 *
 *  \f{align*}{
 *    (v_i^T x + \beta_i)^2 &\leq (v_{n+1}^T x + \beta_{n+1}) y_i & \text{for } i \in [n] \\
 *    \gamma                &\leq (v_{n+1}^T x + \beta_{n+1}) y_{n+1} &
 *  \f}
 *
 *  Now we want to separate one of the small rotated cones. We first transform it into standard form:
 *  \f[
 *    \sqrt{4(v_i^T x + \beta_i)^2 + (v_{n+1}^T x + \beta_{n+1} - y_i)^2} - v_{n+1}^T x - \beta_{n+1} - y_i \leq 0.
 *  \f]
 *
 *  Let the left-hand-side be called \f$f(x,y)\f$. We now compute the gradient by
 *  \f{align*}{
 *    \frac{\delta f}{\delta x_j} &= \frac{(v_i)_j(4v_i^T x + 4\beta_i) + (v_{n+1})_j(v_{n+1}^T x + \beta_{n+1} - y_i)}{\sqrt{4(v_i^T x + \beta_i)^2 + (v_{n+1}^T x + \beta_{n+1} - y_i)^2}} - (v_{n+1})_j \\
 *    \frac{\delta f}{\delta y_i} &= \frac{y_i - v_{n+1}^T x -\beta_{n+1}}{\sqrt{4(v_i^T x + \beta_i)^2 + (v_{n+1}^T x + \beta_{n+1} - y_i)^2}} - 1
 *  \f}
 *
 *  and the gradient cut is then \f$f(x^*, y^*) + \nabla f(x^*,y^*)((x,y) - (x^*, y^*)) \leq 0\f$.
 *  For the \f$\gamma\f$-cone the formula is similar.
 */
static
SCIP_RETCODE generateCutSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONS*            cons,               /**< the constraint that expr is part of */
   SCIP_SOL*             sol,                /**< solution to separate or NULL for the LP solution */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nonlinear handler expression data */
   int                   disaggidx,          /**< index of disaggregation to separate */
   SCIP_Real             mincutviolation,    /**< minimal required cut violation */
   SCIP_Real             rhsval,             /**< value of the rhs term */
   SCIP_ROW**            cut                 /**< pointer to store a cut */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** disvars;
   SCIP_Real* transcoefs;
   int* transcoefsidx;
   int* termbegins;
   int* nnonzeroes;
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* cutvar;
   SCIP_Real cutcoef;
   SCIP_Real fvalue;
   SCIP_Real disvarval;
   SCIP_Real lhsval;
   SCIP_Real constant;
   SCIP_Real denominator;
   int ncutvars;
   int nterms;
   int i;

   assert(expr != NULL);
   assert(cons != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(disaggidx < nlhdlrexprdata->nterms);
   assert(mincutviolation >= 0.0);
   assert(cut != NULL);

   vars = nlhdlrexprdata->vars;
   disvars = nlhdlrexprdata->disvars;
   transcoefs = nlhdlrexprdata->transcoefs;
   transcoefsidx = nlhdlrexprdata->transcoefsidx;
   termbegins = nlhdlrexprdata->termbegins;
   nnonzeroes = nlhdlrexprdata->nnonzeroes;
   nterms = nlhdlrexprdata->nterms;

   /* nterms is equal to n+1 in the description and disaggidx is in {0,...,n}
    * disaggidx = n corresponds to the disaggregation of gamma
    */

   *cut = NULL;

   disvarval = SCIPgetSolVal(scip, sol, disvars[disaggidx]);

   if( disaggidx < nterms - 1 )
   {
      lhsval = evalSingleTerm(scip, nlhdlrexprdata, sol, disaggidx);
      denominator = SQRT(4.0 * SQR(lhsval) + SQR(rhsval - disvarval));
   }
   else
   {
      lhsval = nlhdlrexprdata->constant;
      denominator = SQRT(4.0 * lhsval + SQR(rhsval - disvarval));
   }

   /* compute value of function to be separated (f(x*,y*)) */
   fvalue = denominator - rhsval - disvarval;

   /* if the soc is not violated don't compute cut */
   if( !SCIPisPositive(scip, fvalue) )
   {
      SCIPdebugMsg(scip, "skip cut on disaggregation index %d as violation=%g below feastol\n", disaggidx, fvalue);
      return SCIP_OKAY;
   }

   /* if the denominator is 0 -> the constraint can't be violated */
   assert(!SCIPisZero(scip, denominator));

   /* compute maximum number of variables in cut */
   ncutvars = nnonzeroes[nterms - 1] + 1 + (disaggidx < nterms - 1 ? nnonzeroes[disaggidx] : 0);

   /* create cut */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, ncutvars) );

   /* constant will be grad_f(x*,y*)^T  (x*, y*) */
   constant = 0.0;

   /* a variable could appear on the lhs and rhs, but we add the coefficients separately  */

   /* add terms for v_{disaggidx+1} (disaggidx corresponds to v_{disaggidx+1} in the description) */
   if( disaggidx < nterms - 1 && !SCIPisZero(scip, lhsval) )
   {
      for( i = 0; i < nnonzeroes[disaggidx]; ++i )
      {
         int idx;

         idx = termbegins[disaggidx] + i;
         cutvar = vars[transcoefsidx[idx]];

         /* cutcoef is (the first part of) the partial derivative w.r.t cutvar */
         cutcoef = 4.0 * lhsval * transcoefs[idx] / denominator;

         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

         constant += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
      }
   }

   /* add terms for v_{n+1} */
   for( i = 0; i < nnonzeroes[nterms - 1]; ++i )
   {
      int idx;

      idx = termbegins[nterms - 1] + i;
      cutvar = vars[transcoefsidx[idx]];

      /* cutcoef is the (second part of) the partial derivative w.r.t cutvar */
      cutcoef = (rhsval - disvarval) * transcoefs[idx] / denominator - transcoefs[idx];

      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

      constant += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
   }

   /* add term for disvar */

   /* cutcoef is the the partial derivative w.r.t. the disaggregation variable */
   cutcoef = (disvarval - rhsval) / denominator - 1.0;
   cutvar = disvars[disaggidx];

   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

   constant += cutcoef * SCIPgetSolVal(scip, sol, cutvar);

   /* add side */
   SCIPaddRowprepSide(rowprep, constant - fvalue);

   SCIP_CALL( SCIPcleanupRowprep2(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, SCIPinfinity(scip), NULL) );

   if( SCIPisGT(scip, SCIPgetRowprepViolation(scip, rowprep, sol, NULL), mincutviolation) )
   {
      (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "soc_%p_%d_%d", (void*) expr, disaggidx, SCIPgetNLPs(scip));
      SCIP_CALL( SCIPgetRowprepRowCons(scip, cut, rowprep, cons) );
   }
   else
   {
      SCIPdebugMsg(scip, "rowprep violation %g below mincutviolation %g\n", SCIPgetRowprepViolation(scip, rowprep, sol,
               NULL), mincutviolation);
      /* SCIPprintRowprep(scip, rowprep, NULL); */
   }

   /* free memory */
   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}

/** checks if an expression is quadratic and to collectall occurring expressions
 *
 * @pre @param expr2idx and @param occurringexprs need to be initialized with capacity 2 * nchildren
 *
 * @note: We assume that a linear term always appears before its corresponding
 * quadratic term in quadexpr; this should be ensured by canonicalize
 */
static
SCIP_RETCODE checkAndCollectQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   quadexpr,           /**< candidate for a quadratic expression */
   SCIP_HASHMAP*         expr2idx,           /**< hashmap to store expressions */
   SCIP_CONSEXPR_EXPR**  occurringexprs,     /**< array to store expressions */
   int*                  nexprs,             /**< buffer to store number of expressions */
   SCIP_Bool*            success             /**< buffer to store whether the check was successful */
   )
{
   SCIP_CONSEXPR_EXPR** children;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(quadexpr != NULL);
   assert(expr2idx != NULL);
   assert(occurringexprs != NULL);
   assert(nexprs != NULL);
   assert(success != NULL);

   *nexprs = 0;
   *success = FALSE;
   children = SCIPgetConsExprExprChildren(quadexpr);
   nchildren = SCIPgetConsExprExprNChildren(quadexpr);

   /* iterate in reverse order to ensure that quadratic terms are found before linear terms */
   for( i = nchildren - 1; i >= 0; --i )
   {
      SCIP_CONSEXPR_EXPR* child;

      child = children[i];
      if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         SCIP_CONSEXPR_EXPR* childarg;

         if( SCIPgetConsExprExprPowExponent(child) != 2.0 )
            return SCIP_OKAY;

         childarg = SCIPgetConsExprExprChildren(child)[0];

         if( !SCIPhashmapExists(expr2idx, (void*) childarg) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void*) childarg, *nexprs) );

            /* store the expression so we know it later */
            assert(*nexprs < 2 * nchildren);
            occurringexprs[*nexprs] = childarg;

            ++(*nexprs);
         }

      }
      else if( SCIPisConsExprExprVar(child) && SCIPvarIsBinary(SCIPgetConsExprExprVarVar(child)) )
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
      else if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrProduct(conshdlr) )
      {
         SCIP_CONSEXPR_EXPR* childarg1;
         SCIP_CONSEXPR_EXPR* childarg2;

         if( SCIPgetConsExprExprNChildren(child) != 2 )
            return SCIP_OKAY;

         childarg1 = SCIPgetConsExprExprChildren(child)[0];
         childarg2 = SCIPgetConsExprExprChildren(child)[1];

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
 * @pre @param quadmatrix and @param linvector need to be initialized with size @param nexprs^2 or @param nexprs, resp.
 */
static
void buildQuadExprMatrix(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   quadexpr,           /**< the quadratic expression */
   SCIP_HASHMAP*         expr2idx,           /**< hashmap mapping the occurring expressions to their index */
   int                   nexprs,             /**< number of occurring expressions */
   SCIP_Real*            quadmatrix,         /**< pointer to store (the lower-left triangle of) the quadratic matrix */
   SCIP_Real*            linvector           /**< pointer to store the linear vector */
   )
{
   SCIP_CONSEXPR_EXPR** children;
   SCIP_Real* childcoefs;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(quadexpr != NULL);
   assert(expr2idx != NULL);
   assert(quadmatrix != NULL);
   assert(linvector != NULL);
   assert(SCIPgetConsExprExprHdlrSum(conshdlr) == SCIPgetConsExprExprHdlr(quadexpr));

   children = SCIPgetConsExprExprChildren(quadexpr);
   nchildren = SCIPgetConsExprExprNChildren(quadexpr);
   childcoefs = SCIPgetConsExprExprSumCoefs(quadexpr);

   /* iterate over children to build the constraint defining matrix and vector */
   for( i = 0; i < nchildren; ++i )
   {
      int varpos;

      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         assert(SCIPgetConsExprExprPowExponent(children[i]) == 2.0);
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPgetConsExprExprChildren(children[i])[0]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPgetConsExprExprChildren(children[i])[0]);
         assert(0 <= varpos && varpos < nexprs);

         quadmatrix[varpos * nexprs + varpos] = childcoefs[i];
      }
      else if( SCIPisConsExprExprVar(children[i]) && SCIPvarIsBinary(SCIPgetConsExprExprVarVar(children[i])) )
      {
         assert(SCIPhashmapExists(expr2idx, (void*) children[i]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);
         assert(0 <= varpos && varpos < nexprs);

         quadmatrix[varpos * nexprs + varpos] = childcoefs[i];
      }
      else if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrProduct(conshdlr) )
      {
         int varpos2;

         assert(SCIPgetConsExprExprNChildren(children[i]) == 2);
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPgetConsExprExprChildren(children[i])[0]));
         assert(SCIPhashmapExists(expr2idx, (void*) SCIPgetConsExprExprChildren(children[i])[1]));

         varpos = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPgetConsExprExprChildren(children[i])[0]);
         assert(0 <= varpos && varpos < nexprs);

         varpos2 = SCIPhashmapGetImageInt(expr2idx, (void*) SCIPgetConsExprExprChildren(children[i])[1]);
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

/* tries to fill the nlhdlrexprdata for a potential quadratic SOC expression.
 * We say try because the expression might still turn out not to be an SOC at this point.
 */
static
void tryFillNlhdlrExprDataQuad(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  occurringexprs,     /**< array of all occurring expressions (nvars many) */
   SCIP_Real*            eigvecmatrix,       /**< array containing the Eigenvectors */
   SCIP_Real*            eigvals,            /**< array containing the Eigenvalues */
   SCIP_Real*            bp,                 /**< product of linear vector b * P (eigvecmatrix^t) */
   int                   nvars,              /**< number of variables */
   int                   nposeigvals,        /**< number of positive Eigenvalues */
   int                   nnegeigvals,        /**< number of negative Eigenvalues */
   int*                  nnonzeroes,         /**< pointer to store the nnonzeroes */
   int*                  termbegins,         /**< pointer to store the termbegins */
   SCIP_Real*            transcoefs,         /**< pointer to store the transcoefs */
   int*                  transcoefsidx,      /**< pointer to store the transcoefsidx */
   SCIP_Real*            offsets,            /**< pointer to store the offsets */
   int                   ntranscoefs,        /**< number of transcoefs */
   SCIP_Real*            lhsconstant,        /**< pointer to store the lhsconstant */
   SCIP_Bool*            success             /**< whether the expression is indeed a SOC */
   )
{
   int nextlhsterm = 0;
   int nexttranscoef = 0;
   int nrhstranscoefs = 0;
   int i;
   int j;

   assert(scip != NULL);
   assert(occurringexprs != NULL);
   assert(eigvecmatrix != NULL);
   assert(eigvals != NULL);
   assert(bp != NULL);
   assert(nnonzeroes != NULL);
   assert(termbegins != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(offsets != NULL);
   assert(lhsconstant != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* we have lhsconstant + x^t A x + b x <= 0 and A has a single negative eigenvalue; try to build soc */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real sqrteigval;

      if( SCIPisZero(scip, eigvals[i]) )
         continue;

      if( eigvals[i] > 0.0 )
      {
         sqrteigval = SQRT(eigvals[i]);

         offsets[nextlhsterm] = bp[i] / (2.0 * sqrteigval);
         *lhsconstant -= bp[i] * bp[i] / (4.0 * eigvals[i]);
         termbegins[nextlhsterm] = nexttranscoef;


         /* set transcoefs */
         for( j = 0; j < nvars; ++j )
         {
            if( !SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
            {
               transcoefs[nexttranscoef] = sqrteigval * eigvecmatrix[i * nvars + j];
               transcoefsidx[nexttranscoef] = j;
               ++nnonzeroes[nextlhsterm];

               ++nexttranscoef;
            }
         }

         ++nextlhsterm;
      }
      else
      {
         SCIP_Real rhsvarlb;
         SCIP_Real rhsvarub;
         SCIP_Real signfactor;

         sqrteigval = SQRT(-eigvals[i]);

         offsets[nposeigvals + nnegeigvals - 1] = -bp[i] / (2.0 * sqrteigval);

         rhsvarlb = 0.0;
         rhsvarub = 0.0;

         /* the expression can only be an soc if the resulting rhs term does not change sign;
          * the rhs term is a linear combination of variables, so estimate its bounds
          */
         for( j = 0; j < nvars; ++j )
         {
            SCIP_Real aux;

            if( SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
               continue;

            if( eigvecmatrix[i * nvars + j] > 0.0 )
            {
               aux = SCIPgetConsExprExprActivity(scip, occurringexprs[j]).inf;
               assert(!SCIPisInfinity(scip, aux));
            }
            else
            {
               aux = SCIPgetConsExprExprActivity(scip, occurringexprs[j]).sup;
               assert(!SCIPisInfinity(scip, -aux));
            }

            if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
            {
               rhsvarlb = -SCIPinfinity(scip);
               break;
            }
            else
               rhsvarlb += sqrteigval * eigvecmatrix[i * nvars + j] * aux;
         }
         rhsvarlb += offsets[nposeigvals + nnegeigvals - 1];

         for( j = 0; j < nvars; ++j )
         {
            SCIP_Real aux;

            if( SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
               continue;

            if( eigvecmatrix[i * nvars + j] > 0.0 )
            {
               aux =  SCIPgetConsExprExprActivity(scip, occurringexprs[j]).sup;
               assert(!SCIPisInfinity(scip, -aux));
            }
            else
            {
               aux =  SCIPgetConsExprExprActivity(scip, occurringexprs[j]).inf;
               assert(!SCIPisInfinity(scip, aux));
            }

            if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
            {
               rhsvarub = SCIPinfinity(scip);
               break;
            }
            else
               rhsvarub += sqrteigval * eigvecmatrix[i * nvars + j] * aux;
         }
         rhsvarub += offsets[nposeigvals + nnegeigvals - 1];

         /* since we are just interested in obtaining an interval that contains the real bounds
          * and is tight enough so that we can identify that the rhsvar does not change sign,
          * we swap the bounds in case of numerical troubles
          */
         if( rhsvarub < rhsvarlb )
         {
            assert(SCIPisEQ(scip, rhsvarub, rhsvarlb));
            SCIPswapReals(&rhsvarub, &rhsvarlb);
         }

         /* check whether rhsvar changes sign */
         if( SCIPisLT(scip, rhsvarlb, 0.0) && SCIPisGT(scip, rhsvarub, 0.0) )
            return;

         signfactor = SCIPisLE(scip, rhsvarub, 0.0) ? -1.0 : 1.0;

         offsets[nposeigvals + nnegeigvals - 1] *= signfactor;
         *lhsconstant -= bp[i] * bp[i] / (4.0 * eigvals[i]);

         nrhstranscoefs = 0;

         /* set transcoefs for rhs term */
         for( j = 0; j < nvars; ++j )
         {
            if( !SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
            {
               transcoefs[ntranscoefs - nrhstranscoefs - 1] = signfactor * sqrteigval * eigvecmatrix[i * nvars + j];
               transcoefsidx[ntranscoefs - nrhstranscoefs - 1] = j;
               ++nnonzeroes[nposeigvals + nnegeigvals - 1];

               ++nrhstranscoefs;
            }
         }
         assert(nrhstranscoefs > 0);

         termbegins[nposeigvals + nnegeigvals - 1] = ntranscoefs - nrhstranscoefs;
      }
   }
   assert(nextlhsterm == nposeigvals + nnegeigvals - 1);
   assert(nexttranscoef == ntranscoefs - nrhstranscoefs);

   /* if the lhs constant is nonnegative, it is a SOC... */
   if( !SCIPisNegative(scip, *lhsconstant) )
   {
      *success = TRUE;

      /* ... however, we need lhsconstant to be >= 0 */
      if( *lhsconstant < 0.0 )
         *lhsconstant = 0.0;
   }
}

/** detects if expr <= auxvar is of the form SQRT(sum_i coef_i (expr_i + shift_i)^2 + const) <= auxvar */
static
SCIP_RETCODE detectSocNorm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_CONSEXPR_EXPR** children;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR** vars;
   SCIP_HASHMAP* expr2idx;
   SCIP_HASHSET* linexprs;
   SCIP_Real* childcoefs;
   SCIP_Real* offsets;
   SCIP_Real* transcoefs;
   int* transcoefsidx;
   int* termbegins;
   int* nnonzeroes;
   SCIP_Real constant;
   int nchildren;
   int nvars;
   int nextentry;
   int i;

   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(auxvar != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* relation is not "<=" -> skip */
   if( SCIPgetConsExprExprNLocksPos(expr) == 0 )
      return SCIP_OKAY;

   assert(SCIPgetConsExprExprNChildren(expr) > 0);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   /* check whether expression is a SQRT and has a sum as child with at least 2 children and a non-negative constant */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrPower(conshdlr)
      || SCIPgetConsExprExprPowExponent(expr) != 0.5
      || SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrSum(conshdlr)
      || SCIPgetConsExprExprNChildren(child) < 2
      || SCIPgetConsExprExprSumConstant(child) < 0.0)
   {
      return SCIP_OKAY;
   }

   assert(SCIPvarGetLbLocal(auxvar) >= 0.0);

   /* get children of the sum */
   children = SCIPgetConsExprExprChildren(child);
   nchildren = SCIPgetConsExprExprNChildren(child);
   childcoefs = SCIPgetConsExprExprSumCoefs(child);

   /* TODO: should we initialize the hashmap with size SCIPgetNVars() so that it never has to be resized? */
   SCIP_CALL( SCIPhashmapCreate(&expr2idx, SCIPblkmem(scip), nchildren) );
   SCIP_CALL( SCIPhashsetCreate(&linexprs, SCIPblkmem(scip), nchildren) );

   /* we create coefs array here already, since we have to fill it in first loop in case of success */
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, nchildren) );

   nvars = 0;

   /* check if all children are squares or linear terms with matching square term:
    * if the i-th child is (pow, expr, 2) we store the association <|expr -> i|> in expr2idx and if expr was in
    * linexprs, we remove it from there.
    * if the i-th child is expr' (different from (pow, expr, 2)) and expr' is not a key of expr2idx, we add it
    * to linexprs.
    * if at the end there is any expr in linexpr -> we do not have a separable quadratic function.
    */
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

      exprhdlr = SCIPgetConsExprExprHdlr(children[i]);

      /* handle quadratic expressions children */
      if( exprhdlr == SCIPgetConsExprExprHdlrPower(conshdlr) && SCIPgetConsExprExprPowExponent(children[i]) == 2.0 )
      {
         SCIP_CONSEXPR_EXPR* squarearg = SCIPgetConsExprExprChildren(children[i])[0];

         if( !SCIPhashmapExists(expr2idx, (void*) squarearg) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void *) squarearg, nvars) );
         }

         transcoefs[nvars] = SQRT(childcoefs[i]);

         SCIP_CALL( SCIPhashsetRemove(linexprs, (void*) squarearg) );
         ++nvars;
      }
      /* handle binary variable children */
      else if( SCIPisConsExprExprVar(children[i]) && SCIPvarIsBinary(SCIPgetConsExprExprVarVar(children[i])) )
      {
         assert(!SCIPhashmapExists(expr2idx, (void*) children[i]));
         assert(!SCIPhashsetExists(linexprs, (void*) children[i]));

         SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void *) children[i], nvars) );

         transcoefs[nvars] = SQRT(childcoefs[i]);

         ++nvars;
      }
      else
      {
         if( !SCIPhashmapExists(expr2idx, (void*) children[i]) )
         {
            SCIP_CALL( SCIPhashsetInsert(linexprs, SCIPblkmem(scip), (void*) children[i]) );
         }
      }
   }

   /* there are linear terms without corresponding quadratic terms */
   if( SCIPhashsetGetNElements(linexprs) > 0 )
   {
      SCIPfreeBufferArray(scip, &transcoefs);
      SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   /* add one to variable counter for auxvar */
   ++nvars;

   constant = SCIPgetConsExprExprSumConstant(child);

   /* compute constant of possible soc expression to check its sign */
   for( i = 0; i < nchildren; ++i )
   {
      if( SCIPgetConsExprExprHdlr(children[i]) != SCIPgetConsExprExprHdlrPower(conshdlr)
         || SCIPgetConsExprExprPowExponent(children[i]) != 2.0 )
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

   /* allocate temporary memory to collect data */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nnonzeroes, nvars) );

   /* initialize arrays */
   for( i = 0; i < nvars; ++i )
   {
      transcoefsidx[i] = i;
      termbegins[i] = i;
      offsets[i] = 0.0;
      nnonzeroes[i] = 1;
   }

   /* add data for the auxiliary variable (RHS) */
   vars[nvars - 1] = auxvar;
   transcoefs[nvars - 1] = 1.0;


   /* create required auxiliary variables and fill offsets array */
   nextentry = 0;
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_VAR* argauxvar;

      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr)
         && SCIPgetConsExprExprPowExponent(children[i]) == 2.0 )
      {
         SCIP_CONSEXPR_EXPR* squarearg;

         squarearg = SCIPgetConsExprExprChildren(children[i])[0];
         assert(SCIPhashmapGetImageInt(expr2idx, (void*) squarearg) == nextentry);

         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, squarearg, &argauxvar) );
         assert(argauxvar != NULL);

         vars[nextentry] = argauxvar;
         ++nextentry;
      }
      else if( SCIPisConsExprExprVar(children[i]) && SCIPvarIsBinary(SCIPgetConsExprExprVarVar(children[i])) )
      {
         /* handle binary variable children: no need to create auxvar */
         assert(SCIPhashmapGetImageInt(expr2idx, (void*) children[i]) == nextentry);
         vars[nextentry] = SCIPgetConsExprExprVarVar(children[i]);
         ++nextentry;
      }
      else
      {
         int auxvarpos;

         assert(SCIPhashmapExists(expr2idx, (void*) children[i]));
         auxvarpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);

         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, children[i], &argauxvar) );
         assert(argauxvar != NULL);

         offsets[auxvarpos] = 0.5 * childcoefs[i] / transcoefs[auxvarpos];
      }
   }
   assert(nextentry == nvars - 1);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n", (void*)expr);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, " <= %s\n", SCIPvarGetName(auxvar));
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, offsets, transcoefs, transcoefsidx, termbegins, nnonzeroes, constant,
            nvars, nvars, nvars, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

   /* free memory */
   SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
   SCIPhashmapFree(&expr2idx);
   SCIPfreeBufferArray(scip, &nnonzeroes);
   SCIPfreeBufferArray(scip, &termbegins);
   SCIPfreeBufferArray(scip, &transcoefsidx);
   SCIPfreeBufferArray(scip, &offsets);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &transcoefs);

   return SCIP_OKAY;
}

/** helper method to detect c + sum_i coef_i expr_i^2 - coef_k expr_k^2 <= 0
 *  and c + sum_i coef_i expr_i^2 - coef_k expr_k expr_l <= 0
 *
 *  binary linear variables are interpreted as quadratic terms.
 */
static
SCIP_RETCODE detectSocQuadraticSimple(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_CONSEXPR_EXPR** children;
   SCIP_VAR** vars = NULL;
   SCIP_Real* childcoefs;
   SCIP_Real* offsets = NULL;
   SCIP_Real* transcoefs = NULL;
   int* transcoefsidx = NULL;
   int* termbegins = NULL;
   int* nnonzeroes = NULL;
   SCIP_Real constant;
   SCIP_Real lhsconstant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real rhssign;
   int nposquadterms;
   int nnegquadterms;
   int nposbilinterms;
   int nnegbilinterms;
   int rhsidx;
   int lhsidx;
   int specialtermidx;
   int nchildren;
   int ntranscoefs;
   int nterms;
   int nextentry;
   int i;
   SCIP_Bool ishyperbolic;

   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(auxvar != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is a sum with at least 2 children */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) || SCIPgetConsExprExprNChildren(expr) < 2 )
      return SCIP_OKAY;

   /* get children of the sum */
   children = SCIPgetConsExprExprChildren(expr);
   nchildren = SCIPgetConsExprExprNChildren(expr);
   constant = SCIPgetConsExprExprSumConstant(expr);

   /* we duplicate the child coefficients since we have to manipulate them */
   SCIP_CALL( SCIPduplicateBufferArray(scip, &childcoefs, SCIPgetConsExprExprSumCoefs(expr), nchildren) ); /*lint !e666*/

   /* initialize data */
   lhsidx = -1;
   rhsidx = -1;
   nposquadterms = 0;
   nnegquadterms = 0;
   nposbilinterms = 0;
   nnegbilinterms = 0;
   lhs = (conslhs == SCIP_INVALID ? SCIPvarGetLbGlobal(auxvar) : conslhs); /*lint !e777*/
   rhs = (consrhs == SCIP_INVALID ? SCIPvarGetUbGlobal(auxvar) : consrhs); /*lint !e777*/

   /* check if all children are quadratic or binary linear and count number of positive and negative terms */
   for( i = 0; i < nchildren; ++i )
   {
      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr) &&
            SCIPgetConsExprExprPowExponent(children[i]) == 2.0 )
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
      else if( SCIPisConsExprExprVar(children[i]) && SCIPvarIsBinary(SCIPgetConsExprExprVarVar(children[i])) )
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
      else if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrProduct(conshdlr) &&
            SCIPgetConsExprExprNChildren(children[i]) == 2 )
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

   /* detect case and store lhs/rhs information */
   if( (ishyperbolic && nnegbilinterms > 0) || (!ishyperbolic && nnegquadterms < 2) )
   {
      /* we have -x*y + z^2 ... -> we want to write  z^2 ... <= x*y;
       * or we have -x^2 + y^2  ... -> we want to write y^2 ... <= x^2;
       * in any case, we need a finite rhs
       */
      assert(nnegbilinterms == 1 || nnegquadterms == 1);
      assert(rhsidx != -1);

      /* if rhs is infinity, it can't be soc */
      if( SCIPisInfinity(scip, rhs) )
         goto CLEANUP;

      specialtermidx = rhsidx;
      lhsconstant = constant - rhs;
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
   }

   if( ishyperbolic )
   {
      SCIP_INTERVAL yactivity;
      SCIP_INTERVAL zactivity;

      yactivity = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(children[specialtermidx])[0]);
      zactivity = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(children[specialtermidx])[1]);

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
   else
   {
      SCIP_INTERVAL rhsactivity;

      rhsactivity = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(children[specialtermidx])[0]);

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
    * we have found an soc-representable expression. Now build the nlhdlrexprdata
    *
    * in the non-hyperbolic case, c + sum_i coef_i expr_i^2 - coef_k expr_k^2 <= 0 is transformed to
    * SQRT( c + sum_i coef_i expr_i^2 ) <= coef_k expr_k
    *
    * in the hyperbolic case, c + sum_i coef_i expr_i^2 - coef_k expr_k expr_l <= 0 is transformed to
    * SQRT( 4(c + sum_i coef_i expr_i^2) + (expr_k - expr_l)^2 ) <= expr_k + expr_l
    */

   nterms = ishyperbolic ? nchildren + 1 : nchildren;
   ntranscoefs = ishyperbolic ? nchildren + 3 : nchildren;

   /* SOC was detected, allocate temporary memory for data to collect */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nterms) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &offsets, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nnonzeroes, nterms) );

   *success = TRUE;
   nextentry = 0;

   for( i = 0; i < nchildren; ++i )
   {
      assert(childcoefs[specialtermidx] != 0.0);

      transcoefsidx[i] = i;
      termbegins[i] = i;
      nnonzeroes[i] = 1;

      /* variable and coef for rhs have to be set to the last entry */
      if( i == specialtermidx )
         continue;

      if( SCIPisConsExprExprVar(children[i]) )
      {
         vars[nextentry] = SCIPgetConsExprExprVarVar(children[i]);
         assert(SCIPvarIsBinary(vars[nextentry]));
      }
      else
      {
         assert(SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr));

         /* create the necessary auxiliary variable, if not existent yet */
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(children[i])[0],
                  &vars[nextentry]) );
      }

      if( ishyperbolic )
      {
         /* we eliminate the coefficient of the bilinear term to arrive at standard form */
         assert(4.0 * childcoefs[i] / -childcoefs[specialtermidx] > 0.0);
         transcoefs[nextentry] = SQRT(4.0 * childcoefs[i] / -childcoefs[specialtermidx]);
      }
      else
      {
         assert(childcoefs[i] > 0.0);
         transcoefs[nextentry] = SQRT(childcoefs[i]);
      }

      assert(vars[nextentry] != NULL);

      ++nextentry;
   }

   assert(nextentry == nchildren - 1);

   if( !ishyperbolic )
   {
      /* add data for the rhs variable */
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(children[specialtermidx])[0],
               &vars[nchildren - 1]) );
      assert(vars[nchildren - 1] != NULL);

      assert(childcoefs[specialtermidx] < 0.0);
      transcoefs[nchildren - 1] = rhssign * SQRT(-childcoefs[specialtermidx]);
   }
   else
   {
      /* add data for variables coming from bilinear term */
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(children[specialtermidx])[0],
               &vars[nchildren - 1]) );
      assert(vars[nchildren - 1] != NULL);

      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, SCIPgetConsExprExprChildren(children[specialtermidx])[1],
               &vars[nchildren]) );
      assert(vars[nchildren] != NULL);

      /* on the lhs we have the term (expr_k - expr_l)^2 and the rhs we have +/-(expr_k + expr_l)
       * at this point, vars[nchildren-1] = auxvar(expr_k) and vars[nchildren] = auxvar(expr_l)
       */

      termbegins[nterms - 1] = ntranscoefs - 2;

      nnonzeroes[nterms - 2] = 2;
      nnonzeroes[nterms - 1] = 2;

      transcoefsidx[ntranscoefs - 4] = nchildren - 1;
      transcoefsidx[ntranscoefs - 3] = nchildren;
      transcoefsidx[ntranscoefs - 2] = nchildren - 1;
      transcoefsidx[ntranscoefs - 1] = nchildren;

      transcoefs[ntranscoefs - 4] = 1.0;
      transcoefs[ntranscoefs - 3] = -1.0;
      transcoefs[ntranscoefs - 2] = rhssign;
      transcoefs[ntranscoefs - 1] = rhssign;
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n%f <= ", (void*)expr, lhs);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, "<= %f\n", rhs);
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, offsets, transcoefs, transcoefsidx, termbegins, nnonzeroes, lhsconstant,
            nterms, nterms, ntranscoefs, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

CLEANUP:
   SCIPfreeBufferArrayNull(scip, &nnonzeroes);
   SCIPfreeBufferArrayNull(scip, &termbegins);
   SCIPfreeBufferArrayNull(scip, &transcoefsidx);
   SCIPfreeBufferArrayNull(scip, &transcoefs);
   SCIPfreeBufferArrayNull(scip, &offsets);
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &childcoefs);

   return SCIP_OKAY;
}

/** detects complex quadratic expressions that can be represented by soc constraints.
 *  These are quadratic expressions with either exactly one positive or exactly one negative eigenvalue,
 *  in addition to some extra conditions. One needs to write the quadratic as
 *  sum eigval_i (eigvec_i . x)^2 + c <= -eigval_k (eigvec_k . x)^2, where eigval_k is the negative eigenvalue,
 *  and c must be positive and (eigvec_k . x) must not change sign.
 *  This is described in more details in:
 *
 *  Mahajan, Ashutosh & Munson, Todd. (2010). Exploiting Second-Order Cone Structure for Global Optimization.
 *
 *  The eigen-decomposition is computed using Lapack. Binary linear variables are interpreted as quadratic terms.
 *
 * @todo: In the case -b <= a + x^2 - y^2 <= b, it is possible to represent both sides by soc, Currently, the
 * datastructure can only handle one soc. If this should appear more often, it could be worth to extend it,
 * such that both sides can be handled (see e.g. instance chp_partload).
 *
 * @todo: Since consexpr multiplies as many terms out as possible during presolving, some soc-representable
 * structures cannot be detected, (see e.g. instances bearing or wager). There is currently no obvious way
 * to handle this.
 */
static
SCIP_RETCODE detectSocQuadraticComplex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_CONSEXPR_EXPR** occurringexprs;
   SCIP_VAR** vars;
   SCIP_HASHMAP* expr2idx;
   SCIP_Real* offsets;
   SCIP_Real* transcoefs;
   SCIP_Real* eigvecmatrix;
   SCIP_Real* eigvals;
   SCIP_Real* lincoefs;
   SCIP_Real* bp;
   int* transcoefsidx;
   int* termbegins;
   int* nnonzeroes;
   SCIP_Real constant;
   SCIP_Real lhsconstant;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int nvars;
   int nchildren;
   int npos;
   int nneg;
   int ntranscoefs;
   int i;
   int j;
   SCIP_Bool rhsissoc;
   SCIP_Bool lhsissoc;
   SCIP_Bool isquadratic;

   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(auxvar != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is a sum with at least 2 children */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) || SCIPgetConsExprExprNChildren(expr) < 2 )
   {
      return SCIP_OKAY;
   }

   /* get children of the sum */
   nchildren = SCIPgetConsExprExprNChildren(expr);
   constant = SCIPgetConsExprExprSumConstant(expr);

   /* initialize data */
   vars = NULL;
   offsets = NULL;
   transcoefs = NULL;
   transcoefsidx = NULL;
   nnonzeroes = NULL;
   termbegins = NULL;
   bp = NULL;
   lhs = (conslhs == SCIP_INVALID ? SCIPvarGetLbGlobal(auxvar) : conslhs); /*lint !e777*/
   rhs = (consrhs == SCIP_INVALID ? SCIPvarGetUbGlobal(auxvar) : consrhs); /*lint !e777*/

   SCIP_CALL( SCIPhashmapCreate(&expr2idx, SCIPblkmem(scip), 2 * nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &occurringexprs, 2 * nchildren) );

   /* check if the expression is quadratic and collect all occurring expressions */
   SCIP_CALL( checkAndCollectQuadratic(scip, conshdlr, expr, expr2idx, occurringexprs, &nvars, &isquadratic) );

   if( !isquadratic )
   {
      SCIPfreeBufferArray(scip, &occurringexprs);
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   assert(SCIPhashmapGetNElements(expr2idx) == nvars);

   /* create datastructures for constaint defining matrix and vector */
   SCIP_CALL( SCIPallocClearBufferArray(scip, &eigvecmatrix, nvars * nvars) ); /*lint !e647*/
   SCIP_CALL( SCIPallocClearBufferArray(scip, &lincoefs, nvars) );

   /* build constraint defining matrix (stored in eigvecmatrix) and vector (stored in lincoefs) */
   buildQuadExprMatrix(scip, conshdlr, expr, expr2idx, nvars, eigvecmatrix, lincoefs);

   SCIP_CALL( SCIPallocBufferArray(scip, &eigvals, nvars) );

   /* compute eigenvalues and vectors, A = PDP^t
    * note: eigvecmatrix stores P^t, i.e., P^t_{i,j} = eigvecmatrix[i*nvars+j]
    */
   if( LapackDsyev(TRUE, nvars, eigvecmatrix, eigvals) != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Failed to compute eigenvalues and eigenvectors for expression:\n");

#ifdef SCIP_DEBUG
      SCIPdismantleConsExprExpr(scip, expr);
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
   rhsissoc = (nneg == 1 && SCIPgetConsExprExprNLocksPos(expr) > 0);
   lhsissoc = (npos == 1 && SCIPgetConsExprExprNLocksNeg(expr) > 0);

   /* @TODO: what do we do if both sides are possible? */
   if( !rhsissoc )
   {
      /* if none of the sides is potentially SOC, stop */
      if( !lhsissoc )
         goto CLEANUP;

      /* lhs is potentially SOC, change signs */
      lhsconstant = lhs - constant;

      for( i = 0; i < nvars; ++i )
      {
         eigvals[i] = -eigvals[i];
         bp[i] = -bp[i];
      }
   }
   else
   {
      lhsconstant = constant - rhs;
   }

   /* initialize remaining datastructures for nonlinear handler */
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, npos + nneg) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, npos + nneg) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nnonzeroes, npos + nneg) );

   /* try to fill the nlhdlrexprdata (at this point, it can still fail) */
   tryFillNlhdlrExprDataQuad(scip, occurringexprs, eigvecmatrix, eigvals, bp, nvars, npos, nneg, nnonzeroes,
      termbegins, transcoefs, transcoefsidx, offsets, ntranscoefs, &lhsconstant, success);

   if( !(*success) )
      goto CLEANUP;

   /*
    * at this point, the expression passed all checks and is SOC-representable
    */

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   /* create or get all auxiliary variables */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, occurringexprs[i], &vars[i]) );
      assert(vars[i] != NULL);
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n%f <= ", (void*)expr, lhs);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, "<= %f\n", rhs);
#endif

   /* finally, create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, offsets, transcoefs, transcoefsidx,
         termbegins, nnonzeroes, lhsconstant, nvars, npos + nneg, ntranscoefs, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

CLEANUP:
   SCIPfreeBufferArrayNull(scip, &vars);
   SCIPfreeBufferArrayNull(scip, &nnonzeroes);
   SCIPfreeBufferArrayNull(scip, &termbegins);
   SCIPfreeBufferArrayNull(scip, &transcoefsidx);
   SCIPfreeBufferArrayNull(scip, &transcoefs);
   SCIPfreeBufferArrayNull(scip, &offsets);
   SCIPfreeBufferArrayNull(scip, &bp);
   SCIPfreeBufferArray(scip, &eigvals);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &occurringexprs);
   SCIPfreeBufferArray(scip, &eigvecmatrix);
   SCIPhashmapFree(&expr2idx);

   return SCIP_OKAY;
}

/** helper method to detect SOC structures. The dection runs in 3 steps:
 *
 *  1. check if expression is a norm of the form SQRT(sum_i (sqrcoef_i expr_i^2 + lincoef_i expr_i) + c)
 *  which can be transformed to the form SQRT(sum_i (coef_i expr_i + const_i)^2 + c*) with c* >= 0.
 *    -> this results in the SOC     expr <= auxvar(expr)
 *
 *  2. check if expression represents a quadratic function of one of the following forms (all coefs > 0)
 *     (sum_i coef_i expr_i^2) - coef_k expr_k^2      <= RHS or (sum_i - coef_i expr_i^2) + coef_k expr_k^2      >= LHS
 *  or (sum_i coef_i expr_i^2) - coef_k expr_k expr_l <= RHS or (sum_i - coef_i expr_i^2) + coef_k expr_k expr_l >= LHS
 *  where RHS >=0 or LHS <= 0, respectively. For LHS and RHS we use either the constraint sides if it is a root expr
 *  or the bounds of the auxiliary variable, otherwise.
 *  The cases in the second line above are called hyperbolic or rotated second order cone.
 *    -> this results in the SOC     SQRT((sum_i coef_i expr_i^2) - RHS) <= SQRT(coef_k) expr_k
 *                            or     SQRT(4*(sum_i coef_i expr_i^2) - 4*RHS + (expr_k - expr_l)^2) <= expr_k + expr_l
 *                            (analogously for the LHS cases)
 *
 *  3. check if expression represents a quadratic inequality of the form f(x) = x^TAx + b^Tx + c <= 0 such that f(x)
 *  has exactly one negative eigenvalue plus some extra conditions, see detectSocQuadraticComplex()
 *
 *  Note that step 3 is only performed if paramter compeigenvalues is set to TRUE.
 */
static
SCIP_RETCODE detectSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLR* nlhdlr,             /**< nonlinear handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(expr != NULL);
   assert(auxvar != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(success != NULL);
   assert(conshdlr != NULL);

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* check whether expression is given as norm as described in case 1 above */
   SCIP_CALL( detectSocNorm(scip, conshdlr, expr, auxvar, nlhdlrexprdata, success) );

   if( !(*success) )
   {
      /* check whether expression is a simple soc-respresentable quadratic expression as described in case 2 above */
      SCIP_CALL( detectSocQuadraticSimple(scip, conshdlr, expr, auxvar, conslhs, consrhs, nlhdlrexprdata, success) );
   }

   if( !(*success) && nlhdlrdata->compeigenvalues )
   {
      /* check whether expression is a more complex soc-respresentable quadratic expression as described in case 3 */
      SCIP_CALL( detectSocQuadraticComplex(scip, conshdlr, expr, auxvar, conslhs, consrhs, nlhdlrexprdata, success) );
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrSoc)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrCopyhdlrSoc NULL
#endif

/** callback to free data of handler */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataSoc)
{ /*lint --e{715}*/
   assert(nlhdlrdata != NULL);

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;

}

/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataSoc)
{  /*lint --e{715}*/
   assert(*nlhdlrexprdata != NULL);

   SCIP_CALL( freeNlhdlrExprData(scip, nlhdlrexprdata) );

   return SCIP_OKAY;
}


/** callback to be called in initialization */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINIT(nlhdlrInitSoc)
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
SCIP_DECL_CONSEXPR_NLHDLREXIT(nlhdlrExitSoc)
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
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectSoc)
{ /*lint --e{715}*/
   SCIP_VAR* auxvar;
   SCIP_Real conslhs;
   SCIP_Real consrhs;

   assert(expr != NULL);

   /* TODO is it worth to detect during presolving and then try to apply some bound strengthening? */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      return SCIP_OKAY;

   auxvar = SCIPgetConsExprExprAuxVar(expr);
   assert(auxvar != NULL);

   conslhs = (cons == NULL ? SCIP_INVALID : SCIPgetLhsConsExpr(scip, cons));
   consrhs = (cons == NULL ? SCIP_INVALID : SCIPgetRhsConsExpr(scip, cons));

   SCIP_CALL( detectSOC(scip, conshdlr, nlhdlr, expr, auxvar, conslhs, consrhs, nlhdlrexprdata, success) );

   if( *success )
   {
      /* create variables for cone disaggregation. @TODO: don't do this if nvars = 2 */
      SCIP_CALL( createDisaggrVars(scip, expr, (*nlhdlrexprdata)) );

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_Real lhsval;
         SCIP_Real rhsval;
         SCIP_Real disvarval;
         int termstart;
         int nterms;
         int i;
         int k;

         /*  the debug solution value of the disaggregation variables is set to
          *      (v_i^T x + beta_i)^2 / (v_{n+1}^T x + beta_{n+1})
          *  or                 gamma / (v_{n+1}^T x + beta_{n+1})
          *  if (v_{n+1}^T x + beta_{n+1}) is different from 0.
          *  Otherwise, the debug solution value is set to 0.
          */

         nterms = (*nlhdlrexprdata)->nterms;

         /* set value of rhs */
         termstart = (*nlhdlrexprdata)->termbegins[nterms - 1];
         rhsval = (*nlhdlrexprdata)->offsets[nterms - 1];

         for( i = termstart; i < termstart + (*nlhdlrexprdata)->nnonzeroes[nterms - 1]; ++i )
         {
            SCIP_VAR* var;
            SCIP_Real varval;

            var = (*nlhdlrexprdata)->vars[(*nlhdlrexprdata)->transcoefsidx[i]];

            SCIP_CALL( SCIPdebugGetSolVal(scip, var, &varval) );
            rhsval += (*nlhdlrexprdata)->transcoefs[i] * varval;
         }

         if( SCIPisZero(scip, rhsval) )
         {
            for( i = 0; i < nterms; ++i )
            {
               SCIP_CALL( SCIPdebugAddSolVal(scip, (*nlhdlrexprdata)->disvars[i], 0.0) );
            }
         }
         else
         {
            /* set value for each disaggregation variable corresponding to quadratic term */
            for( k = 0; k < nterms - 1; ++k )
            {
               termstart = (*nlhdlrexprdata)->termbegins[k];
               lhsval = (*nlhdlrexprdata)->offsets[k];

               for( i = termstart; i < termstart + (*nlhdlrexprdata)->nnonzeroes[k]; ++i )
               {
                  SCIP_VAR* var;
                  SCIP_Real varval;

                  var = (*nlhdlrexprdata)->vars[(*nlhdlrexprdata)->transcoefsidx[i]];

                  SCIP_CALL( SCIPdebugGetSolVal(scip, var, &varval) );
                  lhsval += (*nlhdlrexprdata)->transcoefs[i] * varval;
               }

               disvarval = SQR(lhsval) / rhsval;

               SCIP_CALL( SCIPdebugAddSolVal(scip, (*nlhdlrexprdata)->disvars[k], disvarval) );
            }

            /* set value for disaggregation variable corresponding to constant */
            if( !SCIPisZero(scip, (*nlhdlrexprdata)->constant) )
            {
               disvarval = (*nlhdlrexprdata)->constant / rhsval;

               SCIP_CALL( SCIPdebugAddSolVal(scip, (*nlhdlrexprdata)->disvars[nterms - 1], disvarval) );
            }
         }
      }
#endif
   }

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalauxSoc)
{ /*lint --e{715}*/
   SCIP_CONSHDLR* conshdlr;
   int i;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->vars != NULL);
   assert(nlhdlrexprdata->transcoefs != NULL);
   assert(nlhdlrexprdata->transcoefsidx != NULL);
   assert(nlhdlrexprdata->nnonzeroes != NULL);
   assert(nlhdlrexprdata->nterms > 1);

   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   /* if the original expression is a norm, evaluate w.r.t. the auxiliary variables */
   if( SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrPower(conshdlr) )
   {
      assert(SCIPgetConsExprExprPowExponent(expr) == 0.5);

      /* compute sum_i coef_i expr_i^2 + constant */
      *auxvalue = nlhdlrexprdata->constant;

      for( i = 0; i < nlhdlrexprdata->nterms - 1; ++i )
      {
         SCIP_Real termval;

         termval = evalSingleTerm(scip, nlhdlrexprdata, sol, i);
         *auxvalue += SQR(termval);
      }

      assert(*auxvalue >= 0.0);

      /* compute SQRT(sum_i coef_i expr_i^2 + constant) */
      *auxvalue = SQRT(*auxvalue);
   }
   /* otherwise, evaluate the original quadratic expression w.r.t. the created auxvars of the children */
   else
   {
      SCIP_CONSEXPR_EXPR** children;
      SCIP_Real* childcoefs;
      int nchildren;

      assert(SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr));

      children = SCIPgetConsExprExprChildren(expr);
      childcoefs = SCIPgetConsExprExprSumCoefs(expr);
      nchildren = SCIPgetConsExprExprNChildren(expr);

      *auxvalue = SCIPgetConsExprExprSumConstant(expr);

      for( i = 0; i < nchildren; ++i )
      {
         SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

         exprhdlr = SCIPgetConsExprExprHdlr(children[i]);
         if( exprhdlr == SCIPgetConsExprExprHdlrPower(conshdlr) )
         {
            SCIP_VAR* argauxvar;

            assert(SCIPgetConsExprExprPowExponent(children[i]) == 2.0);

            argauxvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[0]);
            assert(argauxvar != NULL);

            *auxvalue += childcoefs[i] * SCIPgetSolVal(scip, sol, argauxvar);
         }
         else if( exprhdlr == SCIPgetConsExprExprHdlrProduct(conshdlr) )
         {
            SCIP_VAR* argauxvar1;
            SCIP_VAR* argauxvar2;

            assert(SCIPgetConsExprExprNChildren(children[i]) == 2);

            argauxvar1 = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[0]);
            argauxvar2 = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[1]);
            assert(argauxvar1 != NULL);
            assert(argauxvar2 != NULL);

            *auxvalue += childcoefs[i] * SCIPgetSolVal(scip, sol, argauxvar1) * SCIPgetSolVal(scip, sol, argauxvar2);
         }
         else
         {
            SCIP_VAR* argauxvar;

            argauxvar = SCIPgetConsExprExprAuxVar(children[i]);
            assert(argauxvar != NULL);

            *auxvalue += childcoefs[i] * SCIPgetSolVal(scip, sol, argauxvar);
         }
      }
   }

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSINITLP) */
static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaSoc)
{ /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   /* create the disaggregation row and store it in nlhdlrexprdata */
   SCIP_CALL( createDisaggrRow(scip, conshdlr, expr, nlhdlrexprdata) );

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
static
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(nlhdlrExitSepaSoc)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   /* free disaggreagation row */
   if( nlhdlrexprdata->disrow != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &nlhdlrexprdata->disrow) );
   }

   return SCIP_OKAY;
}


/** nonlinear handler separation callback */
static
SCIP_DECL_CONSEXPR_NLHDLRENFO(nlhdlrEnfoSoc)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_Real rhsval;
   int depth;
   int naggrs;
   int k;
   SCIP_Bool infeasible;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->disrow != NULL);

   *result = SCIP_DIDNOTFIND;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( SCIPgetCurrentNode(scip) != nlhdlrdata->prevnode )
   {
      nlhdlrdata->nenfocalls = 0;
      nlhdlrdata->prevnode = SCIPgetCurrentNode(scip);
   }

   /* only call separator a given number of times at each node */
   depth = SCIPgetDepth(scip);
   if( (depth == 0 && nlhdlrdata->maxenforoundsroot >= 0 && nlhdlrdata->nenfocalls >= nlhdlrdata->maxenforoundsroot)
      || (depth > 0 && nlhdlrdata->maxenforounds >= 0 && nlhdlrdata->nenfocalls >= nlhdlrdata->maxenforounds)
      || (nlhdlrdata->enfofreq == 0 && depth != 0)
      || (nlhdlrdata->enfofreq > 0 && depth % nlhdlrdata->enfofreq != 0) )
   {
      SCIPdebugMsg(scip, "not running at depth=%d and nenfocalls=%d due to timing parameters (maxenforoundsroot=%d,\
         maxenforounds=%d, enfofreq=%d)\n", depth, nlhdlrdata->nenfocalls, nlhdlrdata->maxenforoundsroot,
            nlhdlrdata->maxenforounds, nlhdlrdata->enfofreq);
      return SCIP_OKAY;
   }

   ++nlhdlrdata->nenfocalls;

   naggrs = SCIPisZero(scip, nlhdlrexprdata->constant) ? nlhdlrexprdata->nterms-1 : nlhdlrexprdata->nterms;

   /* check whether the aggregation row is in the LP */
   if( !SCIProwIsInLP(nlhdlrexprdata->disrow) && SCIPisGE(scip, -SCIPgetRowSolFeasibility(scip, nlhdlrexprdata->disrow,
               sol), SCIPgetLPFeastol(scip) ) )
   {
      SCIP_CALL( SCIPaddRow(scip, nlhdlrexprdata->disrow, FALSE, &infeasible) );
      SCIPdebugMsg(scip, "added aggregation row to LP, cutoff=%d\n", infeasible);

      if( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      *result = SCIP_SEPARATED;
   }

   rhsval = evalSingleTerm(scip, nlhdlrexprdata, sol, nlhdlrexprdata->nterms - 1);

   for( k = 0; k < naggrs && *result != SCIP_CUTOFF; ++k )
   {
      SCIP_ROW* row;

      /* compute gradient cut */
      SCIP_CALL( generateCutSol(scip, expr, cons, sol, nlhdlrexprdata, k, SCIPgetLPFeastol(scip), rhsval, &row) );

      if( row != NULL )
      {
         SCIP_Real cutefficacy;

         cutefficacy = SCIPgetCutEfficacy(scip, sol, row);

         SCIPdebugMsg(scip, "generated row for aggregation %d, efficacy=%g, minefficacy=%g, allowweakcuts=%d\n",
            k, cutefficacy, nlhdlrdata->mincutefficacy, allowweakcuts);

         /* check whether cut is applicable */
         if( SCIPisCutApplicable(scip, row) && (allowweakcuts || cutefficacy >= nlhdlrdata->mincutefficacy) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
            SCIPdebugMsg(scip, "added cut with efficacy %g\n", cutefficacy);

            if( infeasible )
               *result = SCIP_CUTOFF;
            else
               *result = SCIP_SEPARATED;
         }

         /* release row */
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }

   return SCIP_OKAY;
}

/*
 * nonlinear handler specific interface methods
 */

/** includes SOC nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrSoc(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSEXPR_NLHDLR* nlhdlr;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   /* create nonlinear handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );

   nlhdlrdata->nenfocalls = 0;
   nlhdlrdata->prevnode = NULL;

   /* TODO: create and store nonlinear handler specific data here */

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectSoc, nlhdlrEvalauxSoc, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrSoc);
   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrFreehdlrdataSoc);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrFreeExprDataSoc);
   SCIPsetConsExprNlhdlrInitExit(scip, nlhdlr, nlhdlrInitSoc, nlhdlrExitSoc);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, nlhdlrInitSepaSoc, nlhdlrEnfoSoc, NULL, nlhdlrExitSepaSoc);

   /* add soc nlhdlr parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/enfofreq",
         "frequency for enforcement rounds (0: only in root node)",
         &nlhdlrdata->enfofreq, FALSE, DEFAULT_ENFOFREQ, 0, INT_MAX, NULL, NULL) );

   /* add soc nlhdlr parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/maxenforounds",
         "maximal number of enforcement rounds in non-root nodes (-1: unlimited)",
         &nlhdlrdata->maxenforounds, FALSE, DEFAULT_MAXENFOROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/maxenforoundsroot",
         "maximal number of enforcement rounds in the root node (-1: unlimited)",
         &nlhdlrdata->maxenforoundsroot, FALSE, DEFAULT_MAXENFOROUNDSROOT, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/mincutefficacy",
         "Minimum efficacy which a cut needs in order to be added.",
         &nlhdlrdata->mincutefficacy, FALSE, DEFAULT_MINCUTEFFICACY, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/compeigenvalues",
         "Should Eigenvalue computations be done to detect complex cases in quadratic constraints?",
         &nlhdlrdata->compeigenvalues, FALSE, DEFAULT_COMPEIGENVALUES, NULL, NULL) );

   return SCIP_OKAY;
}
