/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_soc.h
 * @brief  nonlinear handler for second order cone constraints \f$\sqrt{\gamma + \sum_{i=1}^{n} (\alpha_i\, (v_i^T x + \beta_i))^2} \leq \alpha_{n+1}\, (x_{n+1}+\beta_{n+1})\f$

 * @author Benjamin Mueller
 * @author Fabian Wegscheider
 *
 * @todo Add row that is stored in the nonlinear handler expression data to the LP if not happened so far.
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
#define NLHDLR_NAME         "soc"
#define NLHDLR_DESC         "soc nonlinear handler"
#define NLHDLR_PRIORITY     100
#define DEFAULT_ALLOWEIGENVALUECOMPS TRUE

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_VAR**            vars;               /**< variables appearing on both sides (x) */
   SCIP_Real*            coefs;              /**< coefficients of both sides (alpha_i) */
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
   SCIP_ROW*             row;                /**< disaggregation row */
};

struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_Bool alloweigenvaluecomps;           /**< whether Eigenvalue computations should be done to detect complex cases */
};

/*
 * Local methods
 */

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

      if( nlhdlrexprdata->coefs[i] != 1.0 )
         SCIPinfoMessage(scip, NULL, "%f*", nlhdlrexprdata->coefs[i]);
      SCIPinfoMessage(scip, NULL, "(", nlhdlrexprdata->coefs[i]);

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
   if( nlhdlrexprdata->coefs[nterms-1] != 1.0 )
      SCIPinfoMessage(scip, NULL, "%f*(", nlhdlrexprdata->coefs[nterms-1]);

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

   if( nlhdlrexprdata->coefs[nterms-1] != 1.0 )
      SCIPinfoMessage(scip, NULL, ")");

   SCIPinfoMessage(scip, NULL, "\n");
}

/** helper method to create nonlinear handler expression data */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            vars,               /**< variables appearing ob both sides (x) */
   SCIP_Real*            coefs,              /**< coefficients of both sides (alpha_i) */
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
   assert(coefs != NULL);
   assert(offsets != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(termbegins != NULL);
   assert(nnonzeroes != NULL);
   assert(coefs[nterms-1] != 0.0);
   assert(nlhdlrexprdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexprdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, vars, nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->coefs, coefs, nterms) );
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
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->coefs, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, (*nlhdlrexprdata)->nvars);
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** evaluate a single term of the form v_i^T x + \beta_i */
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
   assert(k >= 0);
   assert(k < nlhdlrexprdata->nterms);

   termstart = nlhdlrexprdata->termbegins[k];
   result = nlhdlrexprdata->offsets[k];

   for( i = 0; i < nlhdlrexprdata->nnonzeroes[k]; ++i )
   {
      SCIP_Real varval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[termstart + i]]);
      result += nlhdlrexprdata->transcoefs[termstart + i] * varval;
   }

   return result;
}

/** helper method to create variables for the cone disaggregation */
static
SCIP_RETCODE createDisaggr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata  /**< nonlinear handler expression data */
   )
{
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   char name[SCIP_MAXSTRLEN];
   int nvars;
   int size;
   int nrhsvars;
   int nterms;
   int i;

   assert(nlhdlrexprdata != NULL);

   nterms = nlhdlrexprdata->nterms;
   nrhsvars = nlhdlrexprdata->nnonzeroes[nterms-1];

   /* check whether constant has a separate entry */
   size = SCIPisZero(scip, nlhdlrexprdata->constant) ? nterms-1 : nterms;
   nvars = 0;

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexprdata->disvars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, size + nrhsvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, size + nrhsvars) );

   /* create disaggregation variables */
   for( i = 0; i < nterms-1; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_%d", (void*) expr, i);
      SCIP_CALL( SCIPcreateVar(scip, &nlhdlrexprdata->disvars[i], name, 0.0, SCIPinfinity(scip),
            0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[i]) );

      vars[nvars] = nlhdlrexprdata->disvars[i];
      coefs[nvars] = 1.0;

      SCIPvarMarkRelaxationOnly(vars[nvars]);
      SCIP_CALL( SCIPaddVarLocksType(scip, vars[nvars], SCIP_LOCKTYPE_MODEL, 1, 1) );

      ++nvars;
   }

   /* add constant <= rhscoef * rhvar * z_i */
   if( !SCIPisZero(scip, nlhdlrexprdata->constant) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_const", (void*) expr);
      SCIP_CALL( SCIPcreateVar(scip, &nlhdlrexprdata->disvars[size - 1], name, 0.0, SCIPinfinity(scip),
            0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[size - 1]) );

      vars[nvars] = nlhdlrexprdata->disvars[size - 1];
      coefs[nvars] = 1.0;

      SCIPvarMarkRelaxationOnly(vars[nvars]);
      SCIP_CALL( SCIPaddVarLocksType(scip, vars[nvars], SCIP_LOCKTYPE_MODEL, 1, 1) );

      ++nvars;
   }

   /* consider RHS variables */
   for( i = nlhdlrexprdata->ntranscoefs - nrhsvars; i < nlhdlrexprdata->ntranscoefs; ++i )
   {
      vars[nvars] = nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[i]];
      coefs[nvars] = -nlhdlrexprdata->transcoefs[i] * nlhdlrexprdata->coefs[nterms-1];
      ++nvars;
   }

   /* create row */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_row", (void*) expr);
   SCIP_CALL( SCIPcreateEmptyRowConshdlr(scip, &nlhdlrexprdata->row, conshdlr, name, -SCIPinfinity(scip),
         nlhdlrexprdata->offsets[nterms-1], FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarsToRow(scip, nlhdlrexprdata->row, nvars, vars, coefs) );

   /* free memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** helper method to free variables for the cone disaggregation */
static
SCIP_RETCODE freeDisaggr(
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

   if( nlhdlrexprdata->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &nlhdlrexprdata->row) );
   }

   /* free memory */
   SCIPfreeBlockMemoryArrayNull(scip, &nlhdlrexprdata->disvars, size);

   return SCIP_OKAY;
}

/** helper method to compute and add a gradient cut for the k-th cone disaggregation */
static
SCIP_RETCODE generateCutSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_SOL*             sol,                /**< solution to separate (might be NULL) */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nonlinear handler expression data */
   int                   k,                  /**< k-th disaggregation */
   SCIP_Real             mincutviolation,    /**< minimal required cut violation */
   SCIP_ROW**            cut                 /**< pointer to store a cut */
   )
{
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* cutvar;
   SCIP_Real cutcoef;
   SCIP_Real value;
   SCIP_Real disvarval;
   SCIP_Real rhsval;
   SCIP_Real lhsval;
   SCIP_Real sideval;
   SCIP_Real denominator;
   int termstartidx;
   int ncutvars;
   int nterms;
   int i;

   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(k < nlhdlrexprdata->nterms + 1);
   assert(mincutviolation >= 0.0);
   assert(cut != NULL);

   nterms = nlhdlrexprdata->nterms;

   *cut = NULL;

   disvarval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->disvars[k]);
   rhsval = evalSingleTerm(scip, nlhdlrexprdata, sol, nterms-1);

   if( k < nterms )
   {
      lhsval = evalSingleTerm(scip, nlhdlrexprdata, sol, k);
      value = nlhdlrexprdata->coefs[k] * SQR(lhsval);
   }
   else
   {
      lhsval = nlhdlrexprdata->constant;
      value = lhsval;
   }

   value -= rhsval * disvarval;
   SCIPdebugMsg(scip, "evaluate disaggregation: value=%g\n", value);

   /* if the cone is not violated or we would divide by 0, don't compute cut */
   if( value <= mincutviolation || SCIPisZero(scip, disvarval) || SCIPisZero(scip, rhsval) )
      return SCIP_OKAY;

   /* compute maximum number of variables in cut */
   ncutvars = (k < nterms ? nlhdlrexprdata->nnonzeroes[k] + nlhdlrexprdata->nnonzeroes[nterms-1] + 1 : 2);

   /* create cut */
   SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
   SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, ncutvars) );

   assert(nlhdlrexprdata->coefs[nterms-1] > 0);

   sideval = 0.0;
   termstartidx = nlhdlrexprdata->termbegins[nterms-1];
   denominator = 2.0 *  SQRT(nlhdlrexprdata->coefs[nterms-1]) * rhsval * disvarval;

   /* add terms for lhs */
   if( k < nterms  && !SCIPisZero(scip, lhsval) )
   {
      termstartidx = nlhdlrexprdata->termbegins[k];

      for( i = 0; i < nlhdlrexprdata->nnonzeroes[k]; ++i )
      {
         assert(nlhdlrexprdata->coefs[k] > 0);

         cutvar = nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[termstartidx + i]];
         cutcoef = SQRT(nlhdlrexprdata->coefs[k]) * nlhdlrexprdata->transcoefs[termstartidx + i];

         if( SCIPisNegative(scip, lhsval) )
            cutcoef = -cutcoef;

         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

         sideval += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
      }
   }

   /* add terms for rhs */
   for( i = 0; i < nlhdlrexprdata->nnonzeroes[nterms-1]; ++i )
   {
      cutvar = nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[termstartidx + i]];
      cutcoef = -nlhdlrexprdata->coefs[nterms-1] * nlhdlrexprdata->transcoefs[termstartidx + i] * disvarval;
      cutcoef /= denominator;

      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

      sideval += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
   }

   /* add term for disvar */
   cutvar = nlhdlrexprdata->disvars[k];
   cutcoef = -rhsval / denominator;

   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

   sideval += cutcoef * SCIPgetSolVal(scip, sol, cutvar);

   /* add side */
   SCIPaddRowprepSide(rowprep, sideval - value);

   if( SCIPisGT(scip, SCIPgetRowprepViolation(scip, rowprep, sol, NULL), mincutviolation) )
   {
      (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "soc_%p_%d", (void*) expr, k);
      SCIP_CALL( SCIPgetRowprepRowConshdlr(scip, cut, rowprep, conshdlr) );
   }

   /* free memory */
   SCIPfreeRowprep(scip, &rowprep);

   return SCIP_OKAY;
}


/** helper method to detect SQRT(sum_i coef_i (expr_i + shift_i)^2 + const) <= auxvar */
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
   SCIP_Real* coefs;
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
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nchildren) );

   nvars = 0;

   /* check if all children are squares or linear terms with matching square term */
   for( i = 0; i < nchildren; ++i )
   {
      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr)
         && SCIPgetConsExprExprPowExponent(children[i]) == 2.0 )
      {
         SCIP_CONSEXPR_EXPR* squarearg = SCIPgetConsExprExprChildren(children[i])[0];

         if( !SCIPhashmapExists(expr2idx, (void*) squarearg) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(expr2idx, (void *) squarearg, nvars) );
         }

         coefs[nvars] = childcoefs[i];

         SCIP_CALL( SCIPhashsetRemove(linexprs, (void*) squarearg) );
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

   if( SCIPhashsetGetNElements(linexprs) > 0 )
   {
      SCIPfreeBufferArray(scip, &coefs);
      SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
      SCIPhashmapFree(&expr2idx);
      return SCIP_OKAY;
   }

   ++nvars;

   /* allocate temporary memory for data to collect */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nnonzeroes, nvars) );

   for( i = 0; i < nvars; ++i )
   {
      transcoefs[i] = 1.0;
      transcoefsidx[i] = i;
      termbegins[i] = i;
      offsets[i] = 0.0;
      nnonzeroes[i] = 1;
   }

   /* add data for the auxiliary variable (RHS) */
   vars[nvars-1] = auxvar;
   coefs[nvars-1] = 1.0;

   nextentry = 0;
   constant = SCIPgetConsExprExprSumConstant(child);

   /* found SOC structure -> create required auxiliary variables */
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
      else
      {
         int auxvarpos;

         assert(SCIPhashmapExists(expr2idx, (void*) children[i]) );
         auxvarpos = SCIPhashmapGetImageInt(expr2idx, (void*) children[i]);

         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, children[i], &argauxvar) );
         assert(argauxvar != NULL);

         offsets[auxvarpos] = 0.5 * childcoefs[i] / coefs[auxvarpos];
         constant -= coefs[auxvarpos] * SQR(offsets[auxvarpos]);
      }
   }

   assert(nextentry == nvars-1);

   *success = TRUE;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n", (void*)expr);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, " <= %s\n", SCIPvarGetName(auxvar));
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, coefs, offsets, transcoefs, transcoefsidx,
         termbegins, nnonzeroes, constant, nvars, nvars, nvars, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

   /* free memory */
   SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
   SCIPhashmapFree(&expr2idx);
   SCIPfreeBufferArray(scip, &nnonzeroes);
   SCIPfreeBufferArray(scip, &termbegins);
   SCIPfreeBufferArray(scip, &transcoefsidx);
   SCIPfreeBufferArray(scip, &transcoefs);
   SCIPfreeBufferArray(scip, &offsets);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &coefs);

   return SCIP_OKAY;
}

/** Helper method to detect c + sum_i coef_i expr_i^2 <= coef_k expr_k^2.
 *  Binary linear variables are interpreted as quadratic terms.
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
   SCIP_VAR** vars;
   SCIP_Real* childcoefs;
   SCIP_Real* coefs;
   SCIP_Real* offsets;
   SCIP_Real* transcoefs;
   int* transcoefsidx;
   int* termbegins;
   int* nnonzeroes;
   SCIP_Real constant;
   SCIP_Real sideconstant;
   SCIP_Real lhs;
   SCIP_Real rhs;
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

   /* check whether expression is a sum with at least 2 quadratic children */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr)
      || SCIPgetConsExprExprNChildren(expr) < 2 )
      return SCIP_OKAY;

   /* get children of the sum */
   children = SCIPgetConsExprExprChildren(expr);
   nchildren = SCIPgetConsExprExprNChildren(expr);
   childcoefs = SCIPgetConsExprExprSumCoefs(expr);
   constant = SCIPgetConsExprExprSumConstant(expr);

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
      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr)
         && SCIPgetConsExprExprPowExponent(children[i]) == 2.0 )
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
      else if( SCIPisConsExprExprVar(children[i])
         && SCIPvarIsBinary(SCIPgetConsExprExprVarVar(children[i])) )
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
      else if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrProduct(conshdlr)
         && SCIPgetConsExprExprNChildren(children[i]) == 2 )
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
         return SCIP_OKAY;
      }

      if( nposquadterms > 1 && nnegquadterms > 1 )
         return SCIP_OKAY;

      if( nposbilinterms + nnegbilinterms > 1 )
         return SCIP_OKAY;

      if( nposbilinterms > 0 && nposquadterms > 0 )
         return SCIP_OKAY;

      if( nnegbilinterms > 0 && nnegquadterms > 0 )
         return SCIP_OKAY;
   }

   if( nposquadterms == nchildren || nnegquadterms == nchildren )
      return SCIP_OKAY;

   assert(nposquadterms <= 1 || nnegquadterms <= 1);
   assert(nposbilinterms + nnegbilinterms <= 1);
   assert(nposbilinterms == 0 || nposquadterms == 0);
   assert(nnegbilinterms == 0 || nnegquadterms == 0);

   /* if a bilinear term is involved, it is a hyperbolic expression */
   ishyperbolic = (nposbilinterms + nnegbilinterms > 0);

   /* detect case and store lhs/rhs information */
   if( nnegquadterms + nnegbilinterms <= 1)
   {
      assert(rhsidx != -1);

      /* if rhs is infinity, it can't be soc */
      if( SCIPgetConsExprExprNLocksPos(expr) == 0 )
         return SCIP_OKAY;

      specialtermidx = rhsidx;
      sideconstant = constant - rhs;
   }
   else
   {
      assert(lhsidx != -1);

      /* if lhs is infinity, it can't be soc */
      if( SCIPgetConsExprExprNLocksNeg(expr) == 0 )
         return SCIP_OKAY;

      specialtermidx = lhsidx;
      sideconstant = lhs - constant;

      /* negate all coefficients */
      for( i = 0; i < nchildren; ++i )
         childcoefs[i] = -childcoefs[i];
   }

   if( ishyperbolic )
   {
      /* one of the expressions in the bilinear term is not non-negative -> no SOC */
      if( SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(children[specialtermidx])[0]).inf < 0.0
         || SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(children[specialtermidx])[1]).inf < 0.0 )
      {
         return SCIP_OKAY;
      }

      sideconstant *= 4.0 / -childcoefs[specialtermidx];
   }
   else
   {
         /* rhs variable is not non-negative -> no SOC */
         if( SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(children[specialtermidx])[0]).inf < 0.0 )
            return SCIP_OKAY;
   }

   if( SCIPisNegative(scip, sideconstant) )
      return SCIP_OKAY;

   /**
    *  we have found an soc-representable expression
    */

   nterms = ishyperbolic ? nchildren + 1 : nchildren;
   ntranscoefs = ishyperbolic ? nchildren + 3 : nchildren;

   /* SOC was detected, allocate temporary memory for data to collect */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nterms) );
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

      transcoefs[i] = 1.0;
      transcoefsidx[i] = i;
      termbegins[i] = i;
      nnonzeroes[i] = 1;

      /* variable and coef for rhs has to be set to the last entry */
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
         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr,
               SCIPgetConsExprExprChildren(children[i])[0], &vars[nextentry]) );
      }

      coefs[nextentry] = ishyperbolic ? -4.0 * childcoefs[i] / childcoefs[specialtermidx]: childcoefs[i];

      assert(vars[nextentry] != NULL);
      assert(coefs[nextentry] > 0.0);

      ++nextentry;
   }

   assert(nextentry == nchildren-1);

   if( !ishyperbolic )
   {
      /* add data for the rhs variable */
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr,
            SCIPgetConsExprExprChildren(children[specialtermidx])[0], &vars[nchildren-1]) );
      assert(vars[nchildren-1] != NULL);

      assert(childcoefs[specialtermidx] < 0.0);
      coefs[nchildren-1] = SQRT(-childcoefs[specialtermidx]);
   }
   else
   {
      /* add data for variables coming from bilinear term */
      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr,
            SCIPgetConsExprExprChildren(children[specialtermidx])[0], &vars[nchildren-1]) );
      assert(vars[nchildren-1] != NULL);

      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr,
            SCIPgetConsExprExprChildren(children[specialtermidx])[1], &vars[nchildren]) );
      assert(vars[nchildren] != NULL);

      termbegins[nterms-1] = ntranscoefs-2;

      nnonzeroes[nterms-2] = 2;
      nnonzeroes[nterms-1] = 2;

      transcoefsidx[ntranscoefs-4] = nchildren-1;
      transcoefsidx[ntranscoefs-3] = nchildren;
      transcoefsidx[ntranscoefs-2] = nchildren-1;
      transcoefsidx[ntranscoefs-1] = nchildren;

      transcoefs[ntranscoefs-4] = 1.0;
      transcoefs[ntranscoefs-3] = -1.0;
      transcoefs[ntranscoefs-2] = 1.0;
      transcoefs[ntranscoefs-1] = 1.0;

      coefs[nterms-2] = 1.0;
      coefs[nterms-1] = 1.0;
   }

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n%f <= ", (void*)expr, lhs);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, "<= %f\n", rhs);
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, coefs, offsets, transcoefs, transcoefsidx,
         termbegins, nnonzeroes, sideconstant, nterms, nterms, ntranscoefs, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

   /* free memory */
   SCIPfreeBufferArray(scip, &nnonzeroes);
   SCIPfreeBufferArray(scip, &termbegins);
   SCIPfreeBufferArray(scip, &transcoefsidx);
   SCIPfreeBufferArray(scip, &transcoefs);
   SCIPfreeBufferArray(scip, &offsets);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** Helper method to detect quadratic expressions that can be represented by soc constraints.
 *  This is sdone by computing and analyzing the Eigenvalue decomposition.
 *  Binary linear variables are interpreted as quadratic terms.
 *
 * @TODO: In the case -b <= a + x^2 - y^2 <= b, it is possible to represent both sides by soc, Currently, the
 * datastructure can only handle one soc. If this should appear more often, it could be worth to extend it,
 * such that both sides can be handled (see e.g. instance chp_partload).
 *
 * @TODO: Since consexpr multiplies as many terms out as possible during presolving, some soc-representable
 * structured cannot be detected (see e.g. instances bearing or wager). There is currently no obvious way
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
   SCIP_CONSEXPR_EXPR** children;
   SCIP_VAR** vars;
   SCIP_HASHMAP* var2idx;
   SCIP_Real* childcoefs;
   SCIP_Real* coefs;
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
   int nextlhsterm;
   int nexttranscoef;
   int nrhstranscoefs;
   int ntranscoefs;
   int i;
   int j;
   SCIP_Bool rhsissoc;
   SCIP_Bool lhsissoc;

   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(auxvar != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* check whether expression is a sum with at least 2 quadratic children */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr)
      || SCIPgetConsExprExprNChildren(expr) < 2 )
   {
      return SCIP_OKAY;
   }

   /* get children of the sum */
   children = SCIPgetConsExprExprChildren(expr);
   nchildren = SCIPgetConsExprExprNChildren(expr);
   childcoefs = SCIPgetConsExprExprSumCoefs(expr);
   constant = SCIPgetConsExprExprSumConstant(expr);

   /* initialize data */
   coefs = NULL;
   offsets = NULL;
   transcoefs = NULL;
   transcoefsidx = NULL;
   nnonzeroes = NULL;
   termbegins = NULL;
   bp = NULL;
   lhs = (conslhs == SCIP_INVALID ? SCIPvarGetLbGlobal(auxvar) : conslhs); /*lint !e788*/
   rhs = (consrhs == SCIP_INVALID ? SCIPvarGetUbGlobal(auxvar) : consrhs); /*lint !e788*/

   /* TODO: should we initialize the hashmap with size SCIPgetNVars() so that it never has to be resized? */
   SCIP_CALL( SCIPhashmapCreate(&var2idx, SCIPblkmem(scip), nchildren) );
   nvars = 0;

   /* iterate over children once to collect variables that appear in quadratic/bilinear terms */
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_VAR* argvar;

      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         if( SCIPgetConsExprExprPowExponent(children[i]) != 2.0 )
         {
            SCIPhashmapFree(&var2idx);
            return SCIP_OKAY;
         }

         argvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[0]);

         if( !SCIPhashmapExists(var2idx, (void*) argvar) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(var2idx, (void*) argvar, nvars) );
            ++nvars;
         }
      }
      else if( SCIPisConsExprExprVar(children[i])
         && SCIPvarIsBinary(SCIPgetConsExprExprVarVar(children[i])) )
      {
         argvar = SCIPgetConsExprExprAuxVar(children[i]);

         if( !SCIPhashmapExists(var2idx, (void*) argvar) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(var2idx, (void*) argvar, nvars) );
            ++nvars;
         }
      }
      else if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrProduct(conshdlr) )
      {
         if( SCIPgetConsExprExprNChildren(children[i]) != 2 )
         {
            SCIPhashmapFree(&var2idx);
            return SCIP_OKAY;
         }

         argvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[0]);

         if( !SCIPhashmapExists(var2idx, (void*) argvar) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(var2idx, (void*) argvar, nvars) );
            ++nvars;
         }

         argvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[1]);

         if( !SCIPhashmapExists(var2idx, (void*) argvar) )
         {
            SCIP_CALL( SCIPhashmapInsertInt(var2idx, (void*) argvar, nvars) );
            ++nvars;
         }
      }
   }

   /* iterate over children a second time to check whether 'linear' terms also appear quadratically */
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_VAR* termauxvar;

      /* skip the already handled children */
      if( SCIPgetConsExprExprHdlr(children[i]) != SCIPgetConsExprExprHdlrPower(conshdlr)
         && SCIPgetConsExprExprHdlr(children[i]) != SCIPgetConsExprExprHdlrProduct(conshdlr) )
      {
         termauxvar = SCIPgetConsExprExprAuxVar(children[i]);

         /* if the auxiliary variable was not found in any quadratic term, it is not soc-representable  */
         if( !SCIPhashmapExists(var2idx, (void*) termauxvar) )
         {
            SCIPhashmapFree(&var2idx);
            return SCIP_OKAY;
         }
      }
   }

   SCIP_CALL( SCIPallocClearBufferArray(scip, &eigvecmatrix, nvars * nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &eigvals, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &lincoefs, nvars) );

   /* iterate over children a third time to build the constraint defining matrix and vector */
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_VAR* argvar;
      int varpos;

      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         assert(SCIPgetConsExprExprPowExponent(children[i]) == 2.0);

         argvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[0]);
         varpos = SCIPhashmapGetImageInt(var2idx, (void*) argvar);
         assert(varpos >= 0);
         assert(varpos < nvars);

         vars[varpos] = argvar;
         eigvecmatrix[varpos * nvars + varpos] = childcoefs[i];
      }
      else if( SCIPisConsExprExprVar(children[i]) && SCIPvarIsBinary(SCIPgetConsExprExprVarVar(children[i])) )
      {
         argvar = SCIPgetConsExprExprAuxVar(children[i]);
         varpos = SCIPhashmapGetImageInt(var2idx, (void*) argvar);
         assert(varpos >= 0);
         assert(varpos < nvars);

         vars[varpos] = argvar;
         eigvecmatrix[varpos * nvars + varpos] = childcoefs[i];
      }
      else if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrProduct(conshdlr) )
      {
         int varpos2;

         assert(SCIPgetConsExprExprNChildren(children[i]) == 2);

         argvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[0]);
         varpos = SCIPhashmapGetImageInt(var2idx, (void*) argvar);
         assert(varpos >= 0);
         assert(varpos < nvars);

         vars[varpos] = argvar;

         argvar = SCIPgetConsExprExprAuxVar(SCIPgetConsExprExprChildren(children[i])[1]);
         varpos2 = SCIPhashmapGetImageInt(var2idx, (void*) argvar);
         assert(varpos2 >= 0);
         assert(varpos2 < nvars);
         assert(varpos != varpos2);

         vars[varpos2] = argvar;
         eigvecmatrix[MIN(varpos, varpos2) * nvars + MAX(varpos, varpos2)] = childcoefs[i] / 2.0;
      }
      else
      {
         argvar = SCIPgetConsExprExprAuxVar(children[i]);
         varpos = SCIPhashmapGetImageInt(var2idx, (void*) argvar);
         assert(varpos >= 0);
         assert(varpos < nvars);

         lincoefs[varpos] = childcoefs[i];
      }
   }

   /* compute eigenvalues and vectors, A = PDP^t
    * note: eigvecmatrix stores P^t
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
      /* if non is potentially SOC, stop */
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
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, npos + nneg) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, npos + nneg) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &termbegins, npos + nneg) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &nnonzeroes, npos + nneg) );

   nextlhsterm = 0;
   nexttranscoef = 0;
   nrhstranscoefs = 0;

   /* we have lhsconstant + x^t A x + b x <= 0 and A has a single negative eigenvalue; try to build soc */
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPisZero(scip, eigvals[i]) )
         continue;

      if( eigvals[i] > 0.0 )
      {
         coefs[nextlhsterm] = sqrt(eigvals[i]);
         offsets[nextlhsterm] = bp[i] / (2 * eigvals[i]);
         lhsconstant -= bp[i] * bp[i] / (4 * eigvals[i]);

         termbegins[nextlhsterm] = nexttranscoef;

         /* set transcoefs */
         for( j = 0; j < nvars; ++j )
         {
            if( !SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
            {
               transcoefs[nexttranscoef] = eigvecmatrix[i * nvars + j];
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

         offsets[npos + nneg - 1] = bp[i] / (2 * eigvals[i]);
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
               aux = SCIPcomputeVarLbGlobal(scip, vars[j]);
               assert(!SCIPisInfinity(scip, aux));
            }
            else
            {
               aux = SCIPcomputeVarUbGlobal(scip, vars[j]);
               assert(!SCIPisInfinity(scip, -aux));
            }

            if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
            {
               rhsvarlb = -SCIPinfinity(scip);
               break;
            }
            else
               rhsvarlb += eigvecmatrix[i * nvars + j] * aux;
         }
         rhsvarlb += offsets[npos + nneg - 1];

         for( j = 0; j < nvars; ++j )
         {
            SCIP_Real aux;

            if( SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
               continue;

            if( eigvecmatrix[i * nvars + j] > 0.0 )
            {
               aux = SCIPcomputeVarUbGlobal(scip, vars[j]);
               assert(!SCIPisInfinity(scip, -aux));
            }
            else
            {
               aux = SCIPcomputeVarLbGlobal(scip, vars[j]);
               assert(!SCIPisInfinity(scip, aux));
            }

            if( SCIPisInfinity(scip, aux) || SCIPisInfinity(scip, -aux) )
            {
               rhsvarub = SCIPinfinity(scip);
               break;
            }
            else
               rhsvarub += eigvecmatrix[i * nvars + j] * aux;
         }
         rhsvarub += offsets[npos + nneg - 1];

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
         if( SCIPisGE(scip, rhsvarlb, 0.0) || SCIPisLE(scip, rhsvarub, 0.0) )
         {
            lhsconstant -= bp[i] * bp[i] / (4 * eigvals[i]);
            coefs[npos + nneg - 1] = SCIPisLE(scip, rhsvarub, 0.0) ? -sqrt(-eigvals[i]) : sqrt(-eigvals[i]);

            nrhstranscoefs = 0;

            /* set transcoefs for rhs term */
            for( j = 0; j < nvars; ++j )
            {
               if( !SCIPisZero(scip, eigvecmatrix[i * nvars + j]) )
               {
                  transcoefs[ntranscoefs - nrhstranscoefs - 1] = eigvecmatrix[i * nvars + j];
                  transcoefsidx[ntranscoefs - nrhstranscoefs - 1] = j;
                  ++nnonzeroes[npos + nneg - 1];

                  ++nrhstranscoefs;
               }
            }
            assert(nrhstranscoefs > 0);

            termbegins[npos + nneg - 1] = ntranscoefs - nrhstranscoefs;
         }
         else
         {
            SCIPdebugMsg(scip, "Failed because rhsterm [%g, %g] changes sign.\n", rhsvarlb, rhsvarub);
            goto CLEANUP;
         }
      }
   }
   assert(nextlhsterm == npos + nneg - 1);
   assert(nexttranscoef == ntranscoefs - nrhstranscoefs);

   /* if the lhs constant is negative, it is non an soc */
   if( SCIPisNegative(scip, lhsconstant) )
      goto CLEANUP;

   *success = TRUE;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n%f <= ", (void*)expr, lhs);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, "<= %f\n", rhs);
#endif

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, coefs, offsets, transcoefs, transcoefsidx,
         termbegins, nnonzeroes, lhsconstant, nvars, npos + nneg, ntranscoefs, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

CLEANUP:
   SCIPfreeBufferArrayNull(scip, &nnonzeroes);
   SCIPfreeBufferArrayNull(scip, &termbegins);
   SCIPfreeBufferArrayNull(scip, &transcoefsidx);
   SCIPfreeBufferArrayNull(scip, &transcoefs);
   SCIPfreeBufferArrayNull(scip, &offsets);
   SCIPfreeBufferArrayNull(scip, &coefs);
   SCIPfreeBufferArrayNull(scip, &bp);
   SCIPfreeBufferArray(scip, &lincoefs);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &eigvals);
   SCIPfreeBufferArray(scip, &eigvecmatrix);
   SCIPhashmapFree(&var2idx);

   return SCIP_OKAY;
}

/** helper method to detect SOC structures */
static
SCIP_RETCODE detectSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_Real             conslhs,            /**< lhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_Real             consrhs,            /**< rhs of the constraint that the expression defines (or SCIP_INVALID) */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(expr != NULL);
   assert(auxvar != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(success != NULL);

   /* no expression constraint handler -> skip */
   conshdlr = SCIPfindConshdlr(scip, "expr");
   if( conshdlr == NULL )
   {
      *success = FALSE;
      return SCIP_OKAY;
   }

   nlhdlrdata = SCIPgetConsExprNlhdlrData(SCIPfindConsExprNlhdlr(conshdlr, NLHDLR_NAME));
   assert(nlhdlrdata != NULL);

   /* check whether expression is given as norm */
   SCIP_CALL( detectSocNorm(scip, conshdlr, expr, auxvar, nlhdlrexprdata, success) );

   if( !(*success) )
   {
      /* check whether expression is a simple soc-respresentable quadratic expression */
      SCIP_CALL( detectSocQuadraticSimple(scip, conshdlr, expr, auxvar, conslhs, consrhs,
            nlhdlrexprdata, success) );
   }

   if( !(*success) && nlhdlrdata->alloweigenvaluecomps )
   {
      /* check whether expression is a more complex soc-respresentable quadratic expression */
      SCIP_CALL( detectSocQuadraticComplex(scip, conshdlr, expr, auxvar, conslhs, consrhs,
            nlhdlrexprdata, success) );
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

   SCIP_CALL( detectSOC(scip, expr, auxvar, conslhs, consrhs, nlhdlrexprdata, success) );

   if( *success )
   {
      /* create variables for cone disaggregation */
      SCIP_CALL( createDisaggr(scip, conshdlr, expr, (*nlhdlrexprdata)) );

#ifdef WITH_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         SCIP_SOL* debugsol;
         SCIP_Real rhsval;
         int i;
         int ndisaggvars;
         int nterms;

         SCIP_CALL( SCIPdebugGetSol(scip, &debugsol) );

         nterms = (*nlhdlrexprdata)->nterms;
         rhsval = (*nlhdlrexprdata)->coefs[nterms-1] * evalSingleTerm(scip, *nlhdlrexprdata, debugsol, nterms-1);

         if( debugsol != NULL )
         {
            ndisaggvars = SCIPisZero(scip, (*nlhdlrexprdata)->constant) ? nterms-1 : nterms;

            for( i = 0; i < ndisaggvars; ++i )
            {
               SCIP_Real disvarval;
               SCIP_Real lhsval;

               if( SCIPisZero(scip, rhsval) )
                  disvarval = 0.0;
               else
               {
                  lhsval = evalSingleTerm(scip, *nlhdlrexprdata, debugsol, i);

                  disvarval = (*nlhdlrexprdata)->coefs[i] * SQR(lhsval) / rhsval;
               }
               /* store debug solution value of disaggregation variable
               * assumes that expression has been evaluated in debug solution before
               */
               SCIP_CALL( SCIPdebugAddSolVal(scip, (*nlhdlrexprdata)->disvars[i], disvarval) );
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
   int i;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->vars != NULL);
   assert(nlhdlrexprdata->coefs != NULL);
   assert(nlhdlrexprdata->transcoefs != NULL);
   assert(nlhdlrexprdata->transcoefsidx != NULL);
   assert(nlhdlrexprdata->nnonzeroes != NULL);
   assert(nlhdlrexprdata->nterms > 1);

   /*
    * TODO the following code is valid if the detected expression is of the form || * || <= auxvar; however, it is not
    *      clear to me what needs to be evaluated if the original expression was quadratic
    */

   /* compute sum_i coef_i expr_i^2 + constant */
   *auxvalue = nlhdlrexprdata->constant;

   for( i = 0; i < nlhdlrexprdata->nterms-1; ++i )
   {
      SCIP_Real termval;

      termval = evalSingleTerm(scip, nlhdlrexprdata, sol, i);
      *auxvalue += nlhdlrexprdata->coefs[i] * SQR(termval);
   }
   assert(*auxvalue >= 0.0);

   /* compute SQRT(sum_i coef_i expr_i^2 + constant) */
   *auxvalue = SQRT(*auxvalue);

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSINITLP) */
static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaSoc)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
static
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(nlhdlrExitSepaSoc)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   /* free variable for cone disaggregation */
   SCIP_CALL( freeDisaggr(scip, nlhdlrexprdata) );

   return SCIP_OKAY;
}


/** nonlinear handler separation callback */
static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(nlhdlrSepaSoc)
{ /*lint --e{715}*/
   int naggrs;
   int k;
   SCIP_Bool infeasible;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->row != NULL);

   *result = SCIP_DIDNOTFIND;

   naggrs = SCIPisZero(scip, nlhdlrexprdata->constant) ? nlhdlrexprdata->nterms-1 : nlhdlrexprdata->nterms;

   /* check whether aggregation row is in the LP */
   if( !SCIProwIsInLP(nlhdlrexprdata->row) )
   {
      SCIP_CALL( SCIPaddRow(scip, nlhdlrexprdata->row, FALSE, &infeasible) );

      if( infeasible )
      {
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      *result = SCIP_SUCCESS;
   }

   for( k = 0; k < naggrs && *result != SCIP_CUTOFF; ++k )
   {
      SCIP_ROW* row;

      /* compute gradient cut */
      SCIP_CALL( generateCutSol(scip, expr, conshdlr, sol, nlhdlrexprdata, k, mincutviolation, &row) );

      if( row != NULL )
      {
         /* check whether cut is applicable */
         if( SCIPisCutApplicable(scip, row) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
            SCIPdebugMsg(scip, "added cut with efficacy %g\n", SCIPgetCutEfficacy(scip, sol, row));

            *ncuts += 1;

            if( infeasible )
               *result = SCIP_CUTOFF;
            else
               *result = SCIP_SUCCESS;
         }

         /* release row */
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }

   return SCIP_OKAY;
}


/** nonlinear handler under/overestimation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateSoc)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrEstimateSoc NULL
#endif


/** nonlinear handler interval evaluation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalSoc)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrIntevalSoc NULL
#endif


/** nonlinear handler callback for reverse propagation */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropSoc)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrReversepropSoc NULL
#endif


/** nonlinear handler callback for branching scores */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(nlhdlrBranchscoreSoc)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrBranchscoreSoc NULL
#endif


/** nonlinear handler callback for reformulation */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(nlhdlrReformulateSoc)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrReformulateSoc NULL
#endif

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

   /* TODO: create and store nonlinear handler specific data here */

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectSoc, nlhdlrEvalauxSoc, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrSoc);
   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrFreehdlrdataSoc);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrFreeExprDataSoc);
   SCIPsetConsExprNlhdlrInitExit(scip, nlhdlr, nlhdlrInitSoc, nlhdlrExitSoc);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, nlhdlrInitSepaSoc, nlhdlrSepaSoc, nlhdlrEstimateSoc, nlhdlrExitSepaSoc);
   SCIPsetConsExprNlhdlrProp(scip, nlhdlr, nlhdlrIntevalSoc, nlhdlrReversepropSoc);
   SCIPsetConsExprNlhdlrBranchscore(scip, nlhdlr, nlhdlrBranchscoreSoc);
   SCIPsetConsExprNlhdlrReformulate(scip, nlhdlr, nlhdlrReformulateSoc);

   /* add soc nlhdlr parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/alloweigenvaluecomps",
         "Should Eigenvalue computations be done to detect complex cases in quadratic constraints?",
         &nlhdlrdata->alloweigenvaluecomps, FALSE, DEFAULT_ALLOWEIGENVALUECOMPS, NULL, NULL) );

   return SCIP_OKAY;
}
