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

#include "scip/cons_varbound.h"
#include "scip/cons_expr_nlhdlr_perspective.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sum.h"
#include "scip/type_cons_expr.h"


/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "perspective"
#define NLHDLR_DESC         "perspective handler for expressions"
#define NLHDLR_PRIORITY     0

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_HASHMAP* scvars;
   SCIP_CONSEXPR_EXPR** perspterms;
   SCIP_CONSEXPR_EXPR** linterms;
   SCIP_VAR* perspbvars;
   SCIP_Real* perspcoefs;
   SCIP_Real* lincoefs;
   int nperspterms;
   int nlinterms;
   int persptermssize;
   int lintermssize;
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
};

/*
 * Local methods
 */

/* a variable is semicontinuous if its bounds depend on the binary variable bvar
 * and bvar == 0 => var = v_off for some real constant v_off.
 * If the bvar is not specified, find the first binary variable that var depends on.
 */
static
SCIP_Bool varIsSemicontinuous(
   SCIP*         scip,       /**< SCIP data structure */
   SCIP_VAR*     var,        /**< the variable to check */
   SCIP_VAR**    bvar,       /**< the binary variable */
   SCIP_HASHMAP* scvars      /**< semicontinuous variable information */
   )
{
   SCIP_SCVARDATA* scv;
   SCIP_CONSHDLR* varboundconshdlr;
   SCIP_CONS** vbconss;
   SCIP_SCVARDATA* olddata;
   int nvbconss, c;
   SCIP_Real lb0, ub0, lb1, ub1;

   assert(scip != NULL);
   assert(var != NULL);

   olddata = SCIPhashmapGetImage(scvars, (void*)var);
   if( olddata != NULL )
   {
      if( olddata->bvar == *bvar )
      {
         return TRUE;
      }
      else
      {
         return FALSE;
      }
   }

   lb0 = -SCIPinfinity(scip);
   ub0 = SCIPinfinity(scip);
   lb1 = -SCIPinfinity(scip);
   ub1 = SCIPinfinity(scip);

   varboundconshdlr = SCIPfindConshdlr(scip, "varbound");
   nvbconss = SCIPconshdlrGetNConss(varboundconshdlr);
   vbconss = SCIPconshdlrGetConss(varboundconshdlr);

   for( c = 0; c < nvbconss; ++c )
   {
      /* look for varbounds containing var and, if specified, bvar */
      if( SCIPgetVarVarbound(scip, vbconss[c]) != var || (*bvar != NULL && SCIPgetVbdvarVarbound(scip, vbconss[c]) != *bvar) )
         continue;

      *bvar = SCIPgetVbdvarVarbound(scip, vbconss[c]);
      SCIP_CALL( SCIPprintCons(scip, vbconss[c], NULL) );
      SCIPinfoMessage(scip, NULL, "\n");

      if( SCIPgetLhsVarbound(scip, vbconss[c]) > lb0 )
         lb0 = SCIPgetLhsVarbound(scip, vbconss[c]);
      if( SCIPgetRhsVarbound(scip, vbconss[c]) < ub0 )
         ub0 = SCIPgetRhsVarbound(scip, vbconss[c]);
      if( SCIPgetLhsVarbound(scip, vbconss[c]) - SCIPgetVbdcoefVarbound(scip, vbconss[c]) > lb1 )
         lb1 = SCIPgetLhsVarbound(scip, vbconss[c]) - SCIPgetVbdcoefVarbound(scip, vbconss[c]);
      if( SCIPgetRhsVarbound(scip, vbconss[c]) - SCIPgetVbdcoefVarbound(scip, vbconss[c]) > ub1 )
         ub1 = SCIPgetRhsVarbound(scip, vbconss[c]) - SCIPgetVbdcoefVarbound(scip, vbconss[c]);
   }

   if( SCIPvarGetLbGlobal(var) >= lb0 )
      lb0 = SCIPvarGetLbGlobal(var);
   if( SCIPvarGetUbGlobal(var) <= ub0 )
      ub0 = SCIPvarGetUbGlobal(var);
   if( SCIPvarGetLbGlobal(var) >= lb1 )
      lb1 = SCIPvarGetLbGlobal(var);
   if( SCIPvarGetUbGlobal(var) <= ub1 )
      ub1 = SCIPvarGetUbGlobal(var);

   if( lb0 == ub0 && (lb0 != lb1 || ub0 != ub1) )
   {
      SCIPallocBlockMemory(scip, &scv);
      scv->bvar = *bvar;
      scv->lb0 = lb0;
      scv->ub0 = ub0;
      scv->lb1 = lb1;
      scv->ub1 = ub1;
      SCIPhashmapInsert(scvars, var, scv);
      return TRUE;
   }

   return FALSE;
}

static
SCIP_RETCODE nlhdlrexprdataAddPerspTerm(SCIP* scip, SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, SCIP_Real coef, SCIP_CONSEXPR_EXPR* expr)
{
   int newsize;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);

   if( nlhdlrexprdata->nperspterms + 1 > nlhdlrexprdata->persptermssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, nlhdlrexprdata->nperspterms + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->perspterms,  nlhdlrexprdata->persptermssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->perspcoefs, nlhdlrexprdata->persptermssize, newsize) );
      nlhdlrexprdata->persptermssize = newsize;
   }
   assert(nlhdlrexprdata->nperspterms + 1 <= nlhdlrexprdata->persptermssize);

   nlhdlrexprdata->perspcoefs[nlhdlrexprdata->nperspterms] = coef;
   nlhdlrexprdata->perspterms[nlhdlrexprdata->nperspterms] = expr;
   nlhdlrexprdata->nperspterms++;

   return SCIP_OKAY;
}

static
SCIP_RETCODE nlhdlrexprdataAddLinTerm(SCIP* scip, SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, SCIP_Real coef, SCIP_CONSEXPR_EXPR* expr)
{
   int newsize;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);

   if( nlhdlrexprdata->nlinterms + 1 > nlhdlrexprdata->lintermssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, nlhdlrexprdata->nlinterms + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->linterms,  nlhdlrexprdata->lintermssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->lincoefs, nlhdlrexprdata->lintermssize, newsize) );
      nlhdlrexprdata->lintermssize = newsize;
   }
   assert(nlhdlrexprdata->nlinterms + 1 <= nlhdlrexprdata->lintermssize);

   nlhdlrexprdata->lincoefs[nlhdlrexprdata->nlinterms] = coef;
   nlhdlrexprdata->linterms[nlhdlrexprdata->nlinterms] = expr;
   nlhdlrexprdata->nlinterms++;

   return SCIP_OKAY;
}

/** frees nlhdlrexprdata structure */
static
void freeNlhdlrExprData(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata   /**< nlhdlr expression data */
)
{
   SCIP_HASHMAPENTRY* entry;
   SCIP_SCVARDATA* data;

   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->linterms), nlhdlrexprdata->lintermssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->lincoefs), nlhdlrexprdata->lintermssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->perspterms), nlhdlrexprdata->persptermssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->perspcoefs), nlhdlrexprdata->persptermssize);
   if(nlhdlrexprdata->scvars != NULL)
   {
      for( int c = 0; c < SCIPhashmapGetNEntries(nlhdlrexprdata->scvars); ++c )
      {
         entry = SCIPhashmapGetEntry(nlhdlrexprdata->scvars, c);
         if( entry != NULL )
         {
            data = (SCIP_SCVARDATA*) SCIPhashmapEntryGetImage(entry);
            SCIPfreeBlockMemory(scip, &data);
         }
      }
      SCIPhashmapFree(&(nlhdlrexprdata->scvars));
   }
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
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINIT(nlhdlrInitPerspective)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSHDLR* linconshdlr;
   SCIP_CONSHDLR* varboundconshdlr;
   int c;

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
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");

         /* TODO extract information with interface functions from cons_linear.h */
      }
      else if( SCIPconsGetHdlr(cons) == varboundconshdlr )
      {
         SCIP_CALL( SCIPprintCons(scip, cons, NULL) );
         SCIPinfoMessage(scip, NULL, "\n");

         /* TODO extract information with interface functions from cons_varbound.h */
      }
      else
      {
         /* TODO check whether other constraint handlers are important */
      }
   }

   /* TODO store data in nlhdlrdata */

   return SCIP_OKAY;
}
#else
#define nlhdlrInitPerspective NULL
#endif


/** callback to be called in deinitialization */
static
SCIP_DECL_CONSEXPR_NLHDLREXIT(nlhdlrExitPerspective)
{  /*lint --e{715}*/

   /* TODO free the memory that has been allocated in HDLRINIT */

   return SCIP_OKAY;
}


/** callback to detect structure in expression tree */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectPerspective)
{ /*lint --e{715}*/
   SCIP_EXPRCURV curvature;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_VAR* var;
   SCIP_VAR* bvar;
   int nvars, v, c, nchildren;
   SCIP_CONSEXPR_EXPR** children;
   SCIP_Real* coefs;
   SCIP_Bool expr_is_onoff;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   SCIPinfoMessage(scip, NULL, "\nCalled perspective detect, expr = ");
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);

   *success = FALSE;
   /* curvature = SCIPgetConsExprExprCurvature(expr); */

   SCIP_CALL(SCIPallocClearBlockMemory(scip, nlhdlrexprdata));
   SCIPgetConsExprExprNVars(scip, conshdlr, expr, &nvars);
   SCIPhashmapCreate(&((*nlhdlrexprdata)->scvars), SCIPblkmem(scip), nvars);

   if( strcmp("sum", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr))) == 0 )
   {
      children = SCIPgetConsExprExprChildren(expr);
      nchildren = SCIPgetConsExprExprNChildren(expr);
      coefs = SCIPgetConsExprExprSumCoefs(expr);
   }
   else
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &children, 1) );
      *children = expr;
      nchildren = 1;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &coefs, 1) );
      *coefs = 1.0;
   }

   for( c = 0; c < nchildren; ++c )
   {
      /* check if variables are semicontinuous */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip)) );
      SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, children[c], varexprs, &nvars) );
      bvar = NULL;

      expr_is_onoff = TRUE;
      if( nvars == 0 ) expr_is_onoff = FALSE;

      /* all variables should be semicontinuous */
      for( v = 0; v < nvars; ++v )
      {
         var = SCIPgetConsExprExprVarVar(varexprs[v]);
         if( !varIsSemicontinuous(scip, var, &bvar, (*nlhdlrexprdata)->scvars))
         {
            expr_is_onoff = FALSE;
            break;
         }
      }
      if( !expr_is_onoff )
      {
         SCIP_CALL( nlhdlrexprdataAddLinTerm(scip, *nlhdlrexprdata, coefs[c], children[c]) );
      }
      else
      {
      SCIPinfoMessage(scip, NULL, "\nFound on/off expression: ");
      SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
         /* TODO: check curvature */
      }

      for( v = 0; v < nvars; ++v )
      {
         SCIPreleaseConsExprExpr(scip, &varexprs[v]);
      }
      SCIPfreeBlockMemoryArray(scip, &varexprs, SCIPgetNTotalVars(scip));
   }

   /* check whether expression is nonlinear, convex or concave, and is not handled by another nonlinear handler */
   if( curvature == SCIP_EXPRCURV_CONVEX && !*enforcedbelow )
   {
      *enforcedbelow = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
      SCIPinfoMessage(scip, NULL, "detected expr to be convex -> can enforce expr <= auxvar\n");
   }
   else if( curvature == SCIP_EXPRCURV_CONCAVE && !*enforcedabove )
   {
      *enforcedabove = TRUE;
      *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
      SCIPinfoMessage(scip, NULL, "detected expr to be concave -> can enforce expr >= auxvar\n");
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "\nExpr is neither convex nor concave");
      *success = FALSE;
   }

   if( strcmp("sum", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr))) != 0 )
   {
      SCIPfreeBlockMemoryArrayNull(scip, &children, 1);
      SCIPfreeBlockMemoryArrayNull(scip, &coefs, 1);
   }

   if( !*success )
   {
      freeNlhdlrExprData(scip, *nlhdlrexprdata);
      SCIPfreeBlockMemory(scip, nlhdlrexprdata);
      *nlhdlrexprdata = NULL;
   }

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
static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(nlhdlrSepaPerspective)
{ /*lint --e{715}*/
   SCIPinfoMessage(scip, NULL, "method of perspective nonlinear handler not implemented yet\n");
   return SCIP_OKAY;
}


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
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(nlhdlrReformulatePerspective)
{ /*lint --e{715}*/

   /* TODO detect structure */

   /* TODO create expression and store it in refexpr */

   /* set refexpr to expr and capture it if no reformulation is possible */
   *refexpr = expr;
   SCIPcaptureConsExprExpr(*refexpr);

   return SCIP_OKAY;
}
#else
#define nlhdlrReformulatePerspective NULL
#endif

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
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, nlhdlrSepaPerspective, NULL, NULL);

   return SCIP_OKAY;
}
