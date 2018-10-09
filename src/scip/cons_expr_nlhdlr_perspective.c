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

/**@file   cons_expr_nlhdlr_perspective.h
 * @brief  perspective nonlinear handler
 * @author Benjamin Mueller
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_perspective.h"
#include "scip/cons_expr.h"
#include "scip/cons_linear.h"
#include "scip/struct_cons.h"
#include "scip/cons_expr_pow.h"
#include "scip/cons_expr_var.h"
#include "cons_expr_sum.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "perspective"
#define NLHDLR_DESC         "perspective handler for expressions"
#define NLHDLR_PRIORITY     0

/*
 * Data structures
 */

/* TODO: fill in the necessary nonlinear handler data */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_Bool   detected;
   SCIP_CONS** vbdconss;       /**< varbound constraints */
   int         nvbdconss;      /**< number of varbound constraints */
};

/*
 * Local methods
 */


static
SCIP_RETCODE findUnivariateVar(
   SCIP* scip,                             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR* expr,               /**< expression to check */
   SCIP_VAR**          expr_var            /**< the variable of a proper univariate function */
   )
{
   int c;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* var;

   SCIPinfoMessage(scip, NULL, "\nExpr handler: %s", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)));

   if( strcmp("sum", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr))) == 0 )
   {
      for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
      {
         child = SCIPgetConsExprExprChildren(expr)[c];

         /* check that the continuous var is the same in all children */
         if( strcmp("pow", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) == 0 )
         {
            if( SCIPgetConsExprExprPowExponent(child) != 2 )
            {
               SCIPinfoMessage(scip, NULL, "\nThe sum expression is not quadratic");
               *expr_var = NULL;
               return SCIP_OKAY;
            }
            var = SCIPgetConsExprExprVarVar(SCIPgetConsExprExprChildren(child)[0]);
            if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
            {
               SCIPinfoMessage(scip, NULL, "\nSquared binary variable, not reformulating");
               *expr_var = NULL;
               return SCIP_OKAY;
            }
            if( *expr_var == NULL )
               *expr_var = var;
            else if( *expr_var != var )
            {
               *expr_var = NULL;
               return SCIP_OKAY;
            }
         }
         else if( strcmp("var", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) == 0 )
         {
            var = SCIPgetConsExprExprVarVar(child);
            if( *expr_var != NULL )
            {
               if( SCIPvarGetType(var) != SCIP_VARTYPE_BINARY && SCIPvarGetType(*expr_var) != SCIP_VARTYPE_BINARY && *expr_var != var )
               {
                  SCIPinfoMessage(scip, NULL, "\nFound two different non-binary variables");
                  *expr_var = NULL;
                  return SCIP_OKAY;
               }
            }
            else
               *expr_var = var;
         }
         else
         {
            if( strcmp("value", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) != 0 )
            {
               SCIPinfoMessage(scip, NULL, "\nThe sum expression is not quadratic");
               *expr_var = NULL;
               return SCIP_OKAY;
            }
         }
      }
   }

   /* TODO: check other constraint handlers, detect other types of nonlinear constraints: pow, exp, etc? */

   if( *expr_var != NULL )
      SCIPinfoMessage(scip, NULL, "\nThe function is quadratic univariate");
   else
      SCIPinfoMessage(scip, NULL, "\nThe function is not of the proper form");

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrCopyhdlrPerspective NULL
#endif

/** callback to free data of handler */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataPerspective)
{ /*lint --e{715}*/
   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}


/** callback to free expression specific data */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataPerspective)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrFreeExprDataPerspective NULL
#endif


/** callback to be called in initialization */
static
SCIP_DECL_CONSEXPR_NLHDLRINIT(nlhdlrInitPerspective)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   return SCIP_OKAY;
}


/** callback to be called in deinitialization */
static
SCIP_DECL_CONSEXPR_NLHDLREXIT(nlhdlrExitPerspective)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   /* TODO free the memory that has been allocated in HDLRINIT */
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrdata)->vbdconss, nlhdlrdata->nvbdconss);

   return SCIP_OKAY;
}


/** callback to detect structure in expression tree */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectPerspective)
{ /*lint --e{715}*/
   *success = FALSE;

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalauxPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** callback to detect structure in expression tree */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINITSEPA(nlhdlrInitSepaPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrInitSepaPerspective NULL
#endif


/** separation deinitialization method of a nonlinear handler (called during CONSEXITSOL) */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLREXITSEPA(nlhdlrExitSepaPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrExitSepaPerspective NULL
#endif


/** nonlinear handler separation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(nlhdlrSepaPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrSepaPerspective NULL
#endif


/** nonlinear handler under/overestimation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimatePerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrEstimatePerspective NULL
#endif


/** nonlinear handler interval evaluation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrIntevalPerspective NULL
#endif


/** nonlinear handler callback for reverse propagation */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREVERSEPROP(nlhdlrReversepropPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrReversepropPerspective NULL
#endif


/** nonlinear handler callback for branching scores */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(nlhdlrBranchscorePerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define nlhdlrBranchscorePerspective NULL
#endif


/** nonlinear handler callback for reformulation */
static
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(nlhdlrReformulatePerspective)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_VAR* bvar;
   SCIP_Real pmin = -SCIPinfinity(scip), pmax = SCIPinfinity(scip), bvarcoef = 0;
   int c;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if( !nlhdlrdata->detected )
   {
      SCIP_CONSHDLR* varboundconshdlr;
      SCIP_CONSHDLR* linconshdlr;

      nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
      assert(nlhdlrdata != NULL);

      /* look for important constraint handlers */
      linconshdlr = SCIPfindConshdlr(scip, "linear");
      varboundconshdlr = SCIPfindConshdlr(scip, "varbound");

      for( c = 0; c < SCIPgetNConss(scip); ++c )
      {
         SCIP_CONS* cons = SCIPgetConss(scip)[c];
         assert(cons != NULL);

         if( SCIPconsGetHdlr(cons) == linconshdlr )
         {
            if( SCIPgetNVarsLinear(scip, cons) == 2 )
            {
               SCIP_VARTYPE vartype1 = SCIPvarGetType(SCIPgetVarsLinear(scip, cons)[0]);
               SCIP_VARTYPE vartype2 = SCIPvarGetType(SCIPgetVarsLinear(scip, cons)[1]);

               if( (vartype1 == SCIP_VARTYPE_BINARY && vartype2 == SCIP_VARTYPE_CONTINUOUS) || (vartype2 == SCIP_VARTYPE_BINARY && vartype1 == SCIP_VARTYPE_CONTINUOUS)  )
               {
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrdata->vbdconss, nlhdlrdata->nvbdconss, nlhdlrdata->nvbdconss+1) );
                  nlhdlrdata->vbdconss[nlhdlrdata->nvbdconss] = cons;
                  nlhdlrdata->nvbdconss++;
               }
            }
         }
         else if( SCIPconsGetHdlr(cons) == varboundconshdlr )
         {
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrdata->vbdconss, nlhdlrdata->nvbdconss, nlhdlrdata->nvbdconss+1) );
            nlhdlrdata->vbdconss[nlhdlrdata->nvbdconss] = cons;
            nlhdlrdata->nvbdconss++;
         }
         else
         {
            /* TODO check whether other constraint handlers are important */
         }
      }

      SCIPinfoMessage(scip, NULL, "\nVar bound constraints found:\n");
      for( c = 0; c < nlhdlrdata->nvbdconss; c++ )
      {
         SCIP_CALL( SCIPprintCons(scip, nlhdlrdata->vbdconss[c], NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
      }

      nlhdlrdata->detected = TRUE;
   }

   /* TODO detect structure */

   /* is the expression univariate (disregarding adding binary variables)? if yes, expr_var stores the continuous var */
   SCIP_VAR* expr_var = NULL;
   SCIP_CALL( findUnivariateVar(scip, expr, &expr_var) );

   /* no reformulation is possible */
   if( expr_var == NULL )
   {
      *refexpr = expr;
      SCIPcaptureConsExprExpr(*refexpr);
      return SCIP_OKAY;
   }

   /* check if expr_var is semicontinuous */
   for( c = 0; c < nlhdlrdata->nvbdconss; ++c )
   {
      SCIP_CONS* vbcons = nlhdlrdata->vbdconss[c];
      int bpos, cpos;
      if( SCIPgetVarsLinear(scip, vbcons)[0] == expr_var )
      {
         bpos = 1;
         cpos = 0;
      }
      else if( SCIPgetVarsLinear(scip, vbcons)[1] == expr_var )
      {
         bpos = 0;
         cpos = 1;
      }
      else
         continue;


      bvar = SCIPgetVarsLinear(scip, vbcons)[bpos]; /* TODO what if there are varbounds with different bvars? */

      SCIP_Bool leq0 = SCIPgetRhsLinear(scip, vbcons) == 0 && SCIPgetLhsLinear(scip, vbcons) == -SCIPinfinity(scip);
      SCIP_Bool geq0 = SCIPgetLhsLinear(scip, vbcons) == 0 && SCIPgetRhsLinear(scip, vbcons) == SCIPinfinity(scip);
      SCIP_Real ccoef = SCIPgetValsLinear(scip, vbcons)[cpos];
      if( (leq0 && ccoef > 0) || (geq0 && ccoef < 0) )
         pmax = -SCIPgetValsLinear(scip, vbcons)[bpos] / ccoef;
      else if( (geq0 && ccoef > 0) || (leq0 && ccoef < 0) )
         pmin = -SCIPgetValsLinear(scip, vbcons)[bpos] / ccoef;
      else /* currently handling only one-sided constraints with 0 lhs or rhs, could be generalised later */
         SCIPinfoMessage(scip, NULL, "\nThe varbound constraint of this form can't be handled now");
      if( pmin != -SCIPinfinity(scip) && pmax != SCIPinfinity(scip) )
         break;
   }
   if( pmin == -SCIPinfinity(scip) || pmax == SCIPinfinity(scip) )
   {
      SCIPinfoMessage(scip, NULL, "\nThe variable does not have proper on/off bounds, %f, %f", pmin, pmax);
      *refexpr = expr;
      SCIPcaptureConsExprExpr(*refexpr);
      return SCIP_OKAY;
   }
   else
      SCIPinfoMessage(scip, NULL, "\nFound a semicont variable, pmin = %f, pmax = %f", pmin, pmax);

   /* find the coefficient of the binary variable, if not found leave it at 0 */
   if( strcmp("sum", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr))) == 0 )
   {
      for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
      {
         SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[c];

         if( strcmp("var", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) == 0 )
         {
            if( SCIPgetConsExprExprVarVar(child) == bvar )
               bvarcoef = SCIPgetConsExprExprSumCoefs(expr)[c];
         }
      }
   }

   SCIPinfoMessage(scip, NULL, "\nc = %f", bvarcoef);

   /* TODO create expression and store it in refexpr */


   /* set refexpr to expr and capture it if no reformulation is possible */
   *refexpr = expr;
   SCIPcaptureConsExprExpr(*refexpr);

   return SCIP_OKAY;
}

/*
 * nonlinear handler specific interface methods
 */

/** includes Perspective nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrPerspective(
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
   BMSclearMemory(nlhdlrdata);

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectPerspective, nlhdlrEvalauxPerspective, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrPerspective);
   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrFreehdlrdataPerspective);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrFreeExprDataPerspective);
   SCIPsetConsExprNlhdlrInitExit(scip, nlhdlr, nlhdlrInitPerspective, nlhdlrExitPerspective);
   SCIPsetConsExprNlhdlrReformulate(scip, nlhdlr, nlhdlrReformulatePerspective);

   return SCIP_OKAY;
}
