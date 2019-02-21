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
   SCIP_VAR**            vars;               /**< variables on left hand side (x) */
   SCIP_Real*            coefs;              /**< coefficients for left hand side (alpha_i) */
   SCIP_Real*            offsets;            /**< offsets for left hand side (beta_i) */
   SCIP_Real*            transcoefs;         /**< non-zeroes of linear transformation vectors on left hand side (v_i) */
   int*                  transcoefsidx;      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  nnonzeroes;         /**< number of non-zeroes in each v_i */
   SCIP_Real             constant;           /**< constant on left hand side (gamma) */
   int                   nvars;              /**< total number of variables appearing on left hand side */
   int                   nterms;             /**< number of summands in the SQRT (excluding gamma) */
   int                   ntranscoefs;        /**< total number of entries in transcoefs */

   SCIP_VAR*             rhsvar;             /**< right-hand side variable (x_{n+1}) */
   SCIP_Real             rhscoef;            /**< coefficient of right-hand side variable (alpha_{n+1}) */
   SCIP_Real             rhsoffset;          /**< offset for variable on right hand side (beta_{n+1}) */

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
   SCIP_VAR**            vars,               /**< variables on left hand side (x) */
   SCIP_Real*            coefs,              /**< coefficients for left hand side (alpha_i) */
   SCIP_Real*            offsets,            /**< offsets for left hand side (beta_i) */
   SCIP_Real*            transcoefs,         /**< non-zeros of linear transformation vectors on left hand side (v_i) */
   int*                  transcoefsidx,      /**< mapping of transformation coefficients to variable indices in vars */
   int*                  nnonzeroes,         /**< number of non-zeros in each v_i */
   SCIP_Real             constant,           /**< constant on left hand side (gamma) */
   int                   nvars,              /**< total number of variables appearing on left hand side */
   int                   nterms,             /**< number of summands in the SQRT (excluding gamma) */
   int                   ntranscoefs,        /**< total number of entries in transcoefs */
   SCIP_VAR*             rhsvar,             /**< right-hand side variable */
   SCIP_Real             rhscoef,            /**< coefficient of right-hand side variable */
   SCIP_Real             rhsoffset,          /**< offset for variable on right hand side (beta_{n+1}) */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata /**< pointer to store nonlinear handler expression data */
   )
{
   int i;

   assert(vars != NULL);
   assert(coefs != NULL);
   assert(offsets != NULL);
   assert(transcoefs != NULL);
   assert(transcoefsidx != NULL);
   assert(nnonzeroes != NULL);
   assert(rhsvar != NULL);
   assert(rhscoef != 0.0);
   assert(nlhdlrexprdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexprdata) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, vars, nvars) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->coefs, coefs, nterms) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->offsets, offsets, nterms) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefs, transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefsidx, transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*nlhdlrexprdata)->nnonzeroes, nnonzeroes, nterms) );
   (*nlhdlrexprdata)->constant = constant;
   (*nlhdlrexprdata)->nvars = nvars;
   (*nlhdlrexprdata)->nterms = nterms;
   (*nlhdlrexprdata)->ntranscoefs = ntranscoefs;
   (*nlhdlrexprdata)->rhsvar = rhsvar;
   (*nlhdlrexprdata)->rhscoef = rhscoef;
   (*nlhdlrexprdata)->rhsoffset = rhsoffset;

   /* capture variables on LHS */
   for( i = 0; i < nvars; ++i )
   {
      assert(vars[i] != NULL);
      SCIPcaptureVar(scip, vars[i]);
   }

   /* capture RHS variable */
   SCIP_CALL( SCIPcaptureVar(scip, rhsvar) );

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

   /* release RHS variable */
   SCIP_CALL( SCIPreleaseVar(scip, &(*nlhdlrexprdata)->rhsvar) );

   /* release LHS variables */
   for( i = 0; i < (*nlhdlrexprdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*nlhdlrexprdata)->vars[i]) );
   }

   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->nnonzeroes, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefsidx, (*nlhdlrexprdata)->ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->transcoefs, (*nlhdlrexprdata)->ntranscoefs);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->offsets, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->coefs, (*nlhdlrexprdata)->nterms);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->vars, (*nlhdlrexprdata)->nvars);
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** helper method to create variables for the cone disaggregation */
/*static
SCIP_RETCODE createDisaggr(
   SCIP*                 scip,               *//**< SCIP data structure *//*
   SCIP_CONSHDLR*        conshdlr,           *//**< expression constraint handler *//*
   SCIP_CONSEXPR_EXPR*   expr,               *//**< expression *//*
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata  *//**< nonlinear handler expression data *//*
   )
{
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   char name[SCIP_MAXSTRLEN];
   SCIP_Bool success;
   SCIP_Bool infeas;
   int nvars;
   int size;
   int i;

   assert(nlhdlrexprdata != NULL);

   *//* check whether constant has a separate entry *//*
   size = SCIPisZero(scip, nlhdlrexprdata->constant) ? nlhdlrexprdata->nexprs : nlhdlrexprdata->nexprs + 1;
   nvars = 0;

   *//* allocate memory *//*
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexprdata->disvars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, size + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, size + 1) );

   *//* create disaggregation variables *//*
   for( i = 0; i < nlhdlrexprdata->nexprs; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_%p_%d", (void*)expr, i);
      SCIP_CALL( SCIPcreateVar(scip, &nlhdlrexprdata->disvars[i], name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE,
            NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[i]) );

      vars[nvars] = nlhdlrexprdata->disvars[i];
      coefs[nvars] = 1.0;
      ++nvars;
   }

   *//* add constant <= rhscoef * rhvar * z_i *//*
   if( !SCIPisZero(scip, nlhdlrexprdata->constant) )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_const_%p", (void*)expr);
      SCIP_CALL( SCIPcreateVar(scip, &nlhdlrexprdata->disvars[size - 1], name, 0.0, SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, nlhdlrexprdata->disvars[size - 1]) );

      vars[nvars] = nlhdlrexprdata->disvars[size - 1];
      coefs[nvars] = 1.0;
      ++nvars;
   }

   *//* consider RHS variable *//*
   vars[nvars] = nlhdlrexprdata->rhsvar;
   coefs[nvars] = -nlhdlrexprdata->rhscoef;
   ++nvars;

   *//* create row *//*
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "conedis_row_%p", (void*)expr);
   SCIP_CALL( SCIPcreateEmptyRowCons(scip, &nlhdlrexprdata->row, conshdlr, name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, TRUE) );
   SCIP_CALL( SCIPaddVarsToRow(scip, nlhdlrexprdata->row, nvars, vars, coefs) );

   *//* free memory *//*
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}*/

/** helper method to free variables for the cone disaggregation */
/*static
SCIP_RETCODE freeDisaggr(
   SCIP*                 scip,               *//**< SCIP data structure *//*
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata  *//**< nonlinear handler expression data *//*
   )
{
   int i;
   int size;

   assert(nlhdlrexprdata != NULL);

   *//* check whether constant has a separate entry *//*
   size = SCIPisZero(scip, nlhdlrexprdata->constant) ? nlhdlrexprdata->nexprs : nlhdlrexprdata->nexprs + 1;

   if( nlhdlrexprdata->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &nlhdlrexprdata->row) );
   }

   *//* release variables *//*
   for( i = 0; i < size; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &nlhdlrexprdata->disvars[i]) );
   }

   *//* free memory *//*
   SCIPfreeBlockMemoryArrayNull(scip, &nlhdlrexprdata->disvars, size);

   return SCIP_OKAY;
}*/

/** helper method to evaluate a cone disaggregation term */
/*
static
SCIP_RETCODE evalDisaggr(
   SCIP*                 scip,               */
/**< SCIP data structure *//*

   SCIP_SOL*             sol,                */
/**< solution to evaluate (might be NULL) *//*

   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, */
/**< nonlinear handler expression data *//*

   int                   k,                  */
/**< k-th disaggregation term *//*

   SCIP_Real*            value,              */
/**< pointer to store the result *//*

   SCIP_Real*            gradient            */
/**< array to store the gradient *//*

   )
{
   SCIP_Real disvarval;
   SCIP_Real rhsval;
   SCIP_Real tmp;

   assert(nlhdlrexprdata != NULL);
   assert(k >= 0 && k < nlhdlrexprdata->nexprs + 1);
   assert(value != NULL);
   assert(gradient != NULL);

   disvarval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->disvars[k]);
   rhsval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->rhsvar);

   if( k < nlhdlrexprdata->nexprs )
   {
      SCIP_Real exprauxval = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(nlhdlrexprdata->exprs[k]));

      tmp = SQRT(4.0 * nlhdlrexprdata->coefs[k] * SQR(exprauxval) + SQR(nlhdlrexprdata->rhscoef * rhsval - disvarval));
      *value = tmp - disvarval - nlhdlrexprdata->rhscoef * rhsval;

      */
/* gradient w.r.t. auxiliary variable for the expression *//*

      gradient[0] = (4.0 * nlhdlrexprdata->coefs[k] * exprauxval) / tmp;
   }
   else
   {
      assert(!SCIPisZero(scip, nlhdlrexprdata->constant));

      tmp = SQRT(4.0 * nlhdlrexprdata->constant + SQR(nlhdlrexprdata->rhscoef * rhsval - disvarval));
      *value = tmp - disvarval - nlhdlrexprdata->rhscoef * rhsval;

      */
/* gradient w.r.t. auxiliary variable for the expression *//*

      gradient[0] = 0.0;
   }

   */
/* gradient w.r.t. the disaggregation variable *//*

   gradient[1] = -nlhdlrexprdata->rhscoef - (nlhdlrexprdata->rhscoef * rhsval - disvarval) / tmp;

   */
/* gradient w.r.t. the RHS variable *//*

   gradient[2] = -1.0 - nlhdlrexprdata->rhscoef * (nlhdlrexprdata->rhscoef * rhsval - disvarval) / tmp;

   return SCIP_OKAY;
}
*/

/** helper method to compute and add a gradient cut for the k-th cone disaggregation */
/*
static
SCIP_RETCODE generateCutSol(
   SCIP*                 scip,               */
/**< SCIP data structure *//*

   SCIP_CONSEXPR_EXPR*   expr,               */
/**< expression *//*

   SCIP_CONSHDLR*        conshdlr,           */
/**< expression constraint handler *//*

   SCIP_SOL*             sol,                */
/**< solution to separate (might be NULL) *//*

   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, */
/**< nonlinear handler expression data *//*

   int                   k,                  */
/**< k-th disaggregation *//*

   SCIP_Real             mincutviolation,    */
/**< minimal required cut violation *//*

   SCIP_ROW**            row                 */
/**< pointer to store a cut *//*

   )
{
   SCIP_Real gradient[3];
   SCIP_Real value;

   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(k < nlhdlrexprdata->nexprs + 1);
   assert(mincutviolation >= 0.0);
   assert(row != NULL);

   *row = NULL;

   SCIP_CALL( evalDisaggr(scip, sol, nlhdlrexprdata, k, &value, gradient) );
   SCIPdebugMsg(scip, "evaluate disaggregation: value=%g gradient=(%g,%g,%g)\n", value, gradient[0], gradient[1], gradient[2]);

   if( value > mincutviolation )
   {
      SCIP_ROWPREP* rowprep;
      SCIP_Real disvarval;
      SCIP_Real exprauxval;
      SCIP_Real rhsval;

      disvarval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->disvars[k]);
      rhsval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->rhsvar);
      exprauxval = k < nlhdlrexprdata->nexprs ? SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(nlhdlrexprdata->exprs[k])) : 0.0;

      */
/* create cut *//*

      SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, SCIP_SIDETYPE_RIGHT, FALSE) );
      SCIP_CALL( SCIPensureRowprepSize(scip, rowprep, 3) );

      */
/* add terms *//*

      if( exprauxval != 0.0 )
      {
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, SCIPgetConsExprExprAuxVar(nlhdlrexprdata->exprs[k]), gradient[0]) );
      }
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, nlhdlrexprdata->disvars[k], gradient[1]) );
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, nlhdlrexprdata->rhsvar, gradient[2]) );

      */
/* add side *//*

      SCIPaddRowprepSide(rowprep, exprauxval * gradient[0] + disvarval * gradient[1] + rhsval * gradient[2] - value);

      if( SCIPisGT(scip, SCIPgetRowprepViolation(scip, rowprep, sol, NULL), mincutviolation) )
      {
         (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "soc_%p_%d", (void*)expr, k);
         SCIP_CALL( SCIPgetRowprepRowCons(scip, row, rowprep, conshdlr) );
      }

      */
/* free memory *//*

      SCIPfreeRowprep(scip, &rowprep);
   }

   return SCIP_OKAY;
}
*/

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
   SCIP_HASHMAP* vars2idx;
   SCIP_Real* offsets;
   SCIP_Real* transcoefs;
   int* transcoefsidx;
   int* nnonzeroes;
   int ntranscoefs;
   int nchildren;
   int nvars;
   int nexttranscoef;
   int i;

   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(auxvar != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* relation is not "<=" -> skip */
   if( SCIPisInfinity(scip, SCIPvarGetUbLocal(auxvar)) )
      return SCIP_OKAY;

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   /* check whether expression is a SQRT and has a sum as children with at least 2 children */
   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrPow(conshdlr) || SCIPgetConsExprExprPowExponent(expr) != 0.5
      || SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrSum(conshdlr) || SCIPgetConsExprExprNChildren(child) < 2 )
      return SCIP_OKAY;

   /* get children of the sum */
   children = SCIPgetConsExprExprChildren(child);
   nchildren = SCIPgetConsExprExprNChildren(child);

   /* check if all children are squares */
   for( i = 0; i < nchildren; ++i )
   {
      if( SCIPgetConsExprExprHdlr(children[i]) != SCIPgetConsExprExprHdlrPow(conshdlr) || SCIPgetConsExprExprPowExponent(children[i]) != 2.0 )
         return SCIP_OKAY;
   }

   ntranscoefs = 0;

   /* iterate over children and count number of summands (1 for non-sum-expressions) */
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CONSEXPR_EXPR* squarearg;
      assert(SCIPgetConsExprExprNChildren(children[i]) == 1);

      squarearg = SCIPgetConsExprExprChildren(children[i])[0];

      if( SCIPgetConsExprExprHdlr(squarearg) == SCIPgetConsExprExprHdlrSum(conshdlr) )
         ntranscoefs += SCIPgetConsExprExprNChildren(squarearg);
      else
         ++ntranscoefs;
   }

   /* allocate temporary memory for data to collect */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &offsets, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefs, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &transcoefsidx, ntranscoefs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nnonzeroes, nchildren) );
   SCIP_CALL( SCIPhashmapCreate(&vars2idx, SCIPblkmem(scip), nvars) );

   nvars = 0;
   nexttranscoef = 0;

   /* found SOC structure -> create required auxiliary variables */
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CONSEXPR_EXPR* squarearg;
      assert(SCIPgetConsExprExprNChildren(children[i]) == 1);

      squarearg = SCIPgetConsExprExprChildren(children[i])[0];

      /* for all but some expression, just make sure the auxiliary variable exists */
      if( SCIPgetConsExprExprHdlr(squarearg) != SCIPgetConsExprExprHdlrSum(conshdlr) )
      {
         SCIP_VAR* childauxvar;
         int auxvarpos;

         SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, squarearg, &childauxvar) );
         assert(childauxvar != NULL);

         /* check if childauxvar is already in vars array (possible for common subexpressions */
         if( SCIPhashmapExists(vars2idx, (void*) childauxvar) )
         {
            auxvarpos = SCIPhashmapGetImageInt(vars2idx, (void*) childauxvar);
         }
         else
         {
            SCIP_CALL( SCIPhashmapInsertInt(vars2idx, (void*) childauxvar, nvars) );
            vars[nvars] = childauxvar;
            auxvarpos = nvars;
            ++nvars;
         }

         transcoefs[nexttranscoef] = 1.0;
         transcoefsidx[nexttranscoef] = auxvarpos;
         offsets[i] = 0.0;
         nnonzeroes[i] = 1;
         ++nexttranscoef;
      }
      /* for sum expressions, do the same for all children of the sum */
      else
      {
         SCIP_CONSEXPR_EXPR** argchildren;
         SCIP_Real* argcoefs;
         int nargchildren;
         int j;

         /* get data of children of square argument (sum expr) */
         argchildren = SCIPgetConsExprExprChildren(squarearg);
         nargchildren = SCIPgetConsExprExprNChildren(squarearg);
         argcoefs = SCIPgetConsExprExprSumCoefs(squarearg);

         for( j = 0; j < nargchildren; ++j )
         {
            SCIP_VAR* childauxvar;
            int auxvarpos;

            /* if the summand is not a variable, create/get the childauxvar, otherwise this just gives the resp. variable */
            SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, argchildren[j], &childauxvar) );
            assert(childauxvar != NULL);

            /* check if childauxvar is already in vars array */
            if( SCIPhashmapExists(vars2idx, (void*) childauxvar) )
            {
               auxvarpos = SCIPhashmapGetImageInt(vars2idx, (void*) childauxvar);
            }
            else
            {
               SCIP_CALL( SCIPhashmapInsertInt(vars2idx, (void*) childauxvar, nvars) );
               vars[nvars] = childauxvar;
               auxvarpos = nvars;
               ++nvars;
            }

            transcoefs[nexttranscoef] = argcoefs[j];
            transcoefsidx[nexttranscoef] = auxvarpos;
            ++nexttranscoef;
         }

         nnonzeroes[i] = nargchildren;
         offsets[i] = SCIPgetConsExprExprSumConstant(squarearg);
      }
   }

   assert(nexttranscoef == ntranscoefs);

   *success = TRUE;

   /* create and store nonlinear handler expression data */
   SCIP_CALL( createNlhdlrExprData(scip, vars, SCIPgetConsExprExprSumCoefs(child), offsets, transcoefs, transcoefsidx,
         nnonzeroes, SCIPgetConsExprExprSumConstant(child), nvars, nchildren, ntranscoefs, auxvar, 1.0, 0.0, nlhdlrexprdata) );
   assert(*nlhdlrexprdata != NULL);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "found SOC structure for expression %p\n", (void*)expr);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, " <= %s\n", SCIPvarGetName(auxvar));
#endif

   /* free memory */
   SCIPhashmapFree(&vars2idx);
   SCIPfreeBufferArray(scip, &nnonzeroes);
   SCIPfreeBufferArray(scip, &transcoefsidx);
   SCIPfreeBufferArray(scip, &transcoefs);
   SCIPfreeBufferArray(scip, &offsets);
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

   /* check whether expression is given as */
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

   if( isroot )
   {
      SCIP_CALL( detectSOC(scip, expr, auxvar, nlhdlrexprdata, success) );
   }
   else
   {
      *success = FALSE;
   }

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalauxSoc)
{ /*lint --e{715}*/
   /*int i;

   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->exprs != NULL);
   assert(nlhdlrexprdata->coefs != NULL);
   assert(nlhdlrexprdata->nexprs > 1);

   *//*
    * TODO the following code is valid if the detected expression is of the form || * || <= auxvar; however, it is not
    *      clear to me what needs to be evaluated if the original expression was quadratic
    *//*

   *//* compute sum_i coef_i expr_i^2 + constant *//*
   *auxvalue = nlhdlrexprdata->constant;*/

   /*for( i = 0; i < nlhdlrexprdata->nexprs; ++i )
   {
      SCIP_CONSEXPR_EXPR* child;
      SCIP_VAR* var;

      child = nlhdlrexprdata->exprs[i];
      assert(child != NULL);

      var = SCIPgetConsExprExprAuxVar(child);
      assert(var != NULL);

      *auxvalue += nlhdlrexprdata->coefs[i] * SQR(SCIPgetSolVal(scip, sol, var));
   }
   assert(*auxvalue >= 0.0);

   *//* compute SQRT(sum_i coef_i expr_i^2 + constant) *//*
   *auxvalue = SQRT(*auxvalue);

   return SCIP_OKAY;*/
}


/** separation deinitialization method of a nonlinear handler (called during CONSINITLP) */
static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaSoc)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   /* create variables for cone disaggregation */
/*
   SCIP_CALL( createDisaggr(scip, conshdlr, expr, nlhdlrexprdata) );
*/

   return SCIP_OKAY;
}


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
static
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(nlhdlrExitSepaSoc)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);

   /* free variable for cone disaggregation */
/*
   SCIP_CALL( freeDisaggr(scip, nlhdlrexprdata) );
*/

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

   /*naggrs = SCIPisZero(scip, nlhdlrexprdata->constant) ? nlhdlrexprdata->nexprs : nlhdlrexprdata->nexprs + 1;

   *//* check whether aggregation row is in the LP *//*
   if( SCIProwIsInLP(nlhdlrexprdata->row) )
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

      *//* compute gradient cut *//*
      SCIP_CALL( generateCutSol(scip, expr, conshdlr, sol, nlhdlrexprdata, k, mincutviolation, &row) );

      if( row != NULL )
      {
         *//* check whether cut is applicable *//*
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

         *//* release row *//*
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }
   }*/

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