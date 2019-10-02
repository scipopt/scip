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
 *
 * @todo Add row that is stored in the nonlinear handler expression data to the LP if not happened so far.
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_soc.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_sum.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "soc"
#define NLHDLR_DESC         "soc nonlinear handler"
#define NLHDLR_PRIORITY     100

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_VAR**            vars;               /**< variables appearing ob both sides (x) */
   SCIP_Real*            coefs;              /**< coefficients of both sides (alpha_i) */
   SCIP_Real*            offsets;            /**< offsets of bot sides (beta_i) */
   SCIP_Real*            transcoefs;         /**< non-zeroes of linear transformation vectors (v_i) */
   int*                  transcoefsidx;      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  termbegins;         /**< starting indices of transcoefs for each term */
   int*                  nnonzeroes;         /**< number of non-zeroes in each v_i */
   SCIP_Real             constant;           /**< constant on left hand side (gamma) */
   int                   nvars;              /**< total number of variables appearing */
   int                   nterms;             /**< number of summands in the SQRT (excluding gamma) +1 for RHS (n+1)*/
   int                   ntranscoefs;        /**< total number of entries in transcoefs */

   /* variables for cone disaggregation */
   SCIP_VAR**            disvars;           /**< disaggregation variables for each expression; entry (nexprs + 1) corresponds to the constant term */
   SCIP_ROW*             row;               /**< disaggregation row */
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
};

/*
 * Local methods
 */

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
      SCIPcaptureVar(scip, vars[i]);
   }

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

   nrhsvars = nlhdlrexprdata->nnonzeroes[nlhdlrexprdata->nterms-1];
   nterms = nlhdlrexprdata->nterms;

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
      SCIP_CALL( SCIPcreateVar(scip, &nlhdlrexprdata->disvars[i], name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE,
            NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[i]) );

      vars[nvars] = nlhdlrexprdata->disvars[i];
      coefs[nvars] = 1.0;
      ++nvars;
   }

   /* add constant <= rhscoef * rhvar * z_i */
   if( !SCIPisZero(scip, nlhdlrexprdata->constant) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_const_%p", (void*) expr);
      SCIP_CALL( SCIPcreateVar(scip, &nlhdlrexprdata->disvars[size - 1], name, 0.0, SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[size - 1]) );

      vars[nvars] = nlhdlrexprdata->disvars[size - 1];
      coefs[nvars] = 1.0;
      ++nvars;
   }

   /* consider RHS variables */
   for( i = 0; i < nrhsvars; ++i)
   {
      int idx = nlhdlrexprdata->ntranscoefs - i - 1;
      vars[nvars] = nlhdlrexprdata->vars[idx];
      coefs[nvars] = -nlhdlrexprdata->transcoefs[idx] * nlhdlrexprdata->coefs[nterms-1];
      ++nvars;
   }

   /* create row */
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_row_%s", (void*) expr);
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
   SCIP_VAR* cutvar;
   SCIP_Real cutcoef;
   SCIP_Real value;
   SCIP_Real disvarval;
   SCIP_Real rhsval;
   SCIP_Real lhsval;
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

   if( value > mincutviolation )
   {
      SCIP_ROWPREP* rowprep;
      SCIP_Real sideval;
      int termstartidx;

      /* create cut */
      SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
      SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, k < nterms ? nterms+1 : 2) );

      sideval = 0;

      /* add terms for lhs */
      if( k < nterms )
      {
         termstartidx = nlhdlrexprdata->termbegins[k];

         for( i = 0; i < nlhdlrexprdata->nnonzeroes[k]; ++i )
         {
            cutvar = nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[termstartidx + i]];
            cutcoef = 2.0 * nlhdlrexprdata->coefs[k] * nlhdlrexprdata->transcoefs[termstartidx + i] * lhsval;

            SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

            sideval += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
         }
      }

      /* add terms for rhs */
      if( !SCIPisZero(scip, disvarval) )
      {
         termstartidx = nlhdlrexprdata->termbegins[nterms-1];

         for( i = 0; i < nlhdlrexprdata->nnonzeroes[nterms-1]; ++i )
         {
            cutvar = nlhdlrexprdata->vars[nlhdlrexprdata->transcoefsidx[termstartidx + i]];
            cutcoef = -nlhdlrexprdata->coefs[nterms-1] * nlhdlrexprdata->transcoefs[termstartidx + i] * disvarval;

            SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

            sideval += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
         }
      }

      /* add term for disvar */
      if( !SCIPisZero(scip, rhsval) )
      {
         cutvar = nlhdlrexprdata->disvars[k];
         cutcoef = -rhsval * nlhdlrexprdata->coefs[nterms-1];

         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, cutvar, cutcoef) );

         sideval += cutcoef * SCIPgetSolVal(scip, sol, cutvar);
      }


      /* add side */
      SCIPaddRowprepSide(rowprep, sideval - value);

      if( SCIPisGT(scip, SCIPgetRowprepViolation(scip, rowprep, sol, NULL), mincutviolation) )
      {
         (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "soc_%p_%d", (void*) expr, k);
         SCIP_CALL( SCIPgetRowprepRowConshdlr(scip, cut, rowprep, conshdlr) );
      }

      /* free memory */
      SCIPfreeRowprep(scip, &rowprep);
   }

   return SCIP_OKAY;
}


/** helper method to detect SQRT(sum_i coef_i (expr_i)^2 + const) <= auxvar */
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
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrPower(conshdlr) || SCIPgetConsExprExprPowExponent(expr) != 0.5
      || SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrSum(conshdlr) || SCIPgetConsExprExprNChildren(child) < 2
      || SCIPgetConsExprExprSumConstant(child) < 0.0)
      return SCIP_OKAY;

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
      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr) && SCIPgetConsExprExprPowExponent(children[i]) == 2.0 )
      {
         SCIP_CONSEXPR_EXPR* squarearg = SCIPgetConsExprExprChildren(children[i])[0];

         if( !SCIPhashmapExists(expr2idx, (void*) squarearg) )
         {
            SCIP_CALL(SCIPhashmapInsertInt(expr2idx, (void *) squarearg, nvars) );
         }

         coefs[nvars] = childcoefs[i];

         SCIPhashsetRemove(linexprs, (void*) squarearg);
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

      if( SCIPgetConsExprExprHdlr(children[i]) == SCIPgetConsExprExprHdlrPower(conshdlr) && SCIPgetConsExprExprPowExponent(children[i]) == 2.0 )
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

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, coefs, offsets, transcoefs, transcoefsidx, termbegins, nnonzeroes,
         constant, nvars, nvars, nvars, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n", (void*)expr);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, " <= %s\n", SCIPvarGetName(auxvar));
#endif

   /* free memory */
   SCIPhashsetFree(&linexprs, SCIPblkmem(scip) );
   SCIPhashmapFree(&expr2idx);
   SCIPfreeBufferArray(scip, &nnonzeroes);
   SCIPfreeBufferArray(scip, &termbegins);
   SCIPfreeBufferArray(scip, &transcoefsidx);
   SCIPfreeBufferArray(scip, &transcoefs);
   SCIPfreeBufferArray(scip, &offsets);
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** helper method to detect SOC structures */
static
SCIP_RETCODE detectSOC(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR*             auxvar,             /**< auxiliary variable */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nonlinear handler expression data */
   SCIP_Bool*            success             /**< pointer to store whether SOC structure has been detected */
   )
{
   SCIP_CONSHDLR* conshdlr;

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

   /* check whether expression is given as norm */
   SCIP_CALL( detectSocNorm(scip, conshdlr, expr, auxvar, nlhdlrexprdata, success) );

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
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataSoc)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of soc nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrFreehdlrdataSoc NULL
#endif


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

   assert(expr != NULL);

   /* TODO is it worth to detect during presolving and then try to apply some bound strengthening? */
   if( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING )
      return SCIP_OKAY;

   auxvar = SCIPgetConsExprExprAuxVar(expr);
   assert(auxvar != NULL);

   SCIP_CALL( detectSOC(scip, expr, auxvar, nlhdlrexprdata, success) );

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

   /* create variables for cone disaggregation */
   SCIP_CALL( createDisaggr(scip, conshdlr, expr, nlhdlrexprdata) );

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

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->row != NULL);

   *result = SCIP_DIDNOTRUN;

   naggrs = SCIPisZero(scip, nlhdlrexprdata->constant) ? nlhdlrexprdata->nterms-1 : nlhdlrexprdata->nterms;

   /* check whether aggregation row is in the LP */
   if( !SCIProwIsInLP(nlhdlrexprdata->row) )
   {
      SCIP_Bool infeasible;

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
      SCIP_Bool cutoff;

      /* compute gradient cut */
      SCIP_CALL( generateCutSol(scip, expr, conshdlr, sol, nlhdlrexprdata, k, mincutviolation, &row) );

      if( row != NULL )
      {
         /* check whether cut is applicable */
         if( SCIPisCutApplicable(scip, row) )
         {
            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &cutoff) );
            SCIPdebugMsg(scip, "added cut with efficacy %g\n", SCIPgetCutEfficacy(scip, sol, row));

            *ncuts += 1;

            if( cutoff )
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
   nlhdlrdata = NULL;

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

   return SCIP_OKAY;
}
