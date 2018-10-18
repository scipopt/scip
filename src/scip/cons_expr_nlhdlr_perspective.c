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
#include "scip/cons_varbound.h"
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
   SCIP_HASHMAP* oldSemiContVars;  /**< previously detected semicontinuous variables */
   SCIP_HASHMAP* semiContVars;     /**< newly detected semicontinuous variables */
   SCIP_Bool     detected;
   SCIP_CONS**   vbdconss;         /**< varbound constraints */
   int           nvbdconss;        /**< number of varbound constraints */
};

/** data of a semicontinuous variable */
struct SCIP_SCVarData
{
   SCIP_Real             xmin;             /**< lower bound in the 'on' state */
   SCIP_Real             xmax;             /**< upper bound in the 'on' state */
   SCIP_VAR*             bvar;             /**< associated binary variable */
};
typedef struct SCIP_SCVarData SCIP_SCVARDATA;

/*
 * Local methods
 */

static
SCIP_Bool varIsSemicontinuous(
   SCIP*                       scip,                  /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLRDATA*   nlhdlrdata,            /**< nonlinear handler data */
   SCIP_VAR*                   var,                   /**< the variable to be checked */
   SCIP_VAR**                  bvar                   /**< the binary variable assosiated with var, NULL if var is not semicontinuous */
)
{
   if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      return FALSE;

   SCIP_SCVARDATA* data = (SCIP_SCVARDATA*) SCIPhashmapGetImage(nlhdlrdata->semiContVars, var);
   if( data != NULL && data->xmin != -SCIPinfinity(scip) && data->xmax != SCIPinfinity(scip) )
   {
      *bvar = data->bvar;
      return TRUE;
   }

   return FALSE;
}

#if 0
static
SCIP_Bool varIsSemicontinuous(
   SCIP*                       scip,                  /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLRDATA*   nlhdlrdata,            /**< nonlinear handler data */
   SCIP_VAR*                   var,                   /**< the variable to be checked */
   SCIP_VAR**                  bvar                   /**< the binary variable assosiated with var, NULL if var is not semicontinuous */
   )
{
   if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
      return FALSE;

   int c;
   SCIP_Real pmin = -SCIPinfinity(scip), pmax = SCIPinfinity(scip);
   *bvar = NULL;

   /* check if var is semicontinuous */
   for( c = 0; c < nlhdlrdata->nvbdconss; ++c )
   {
      SCIP_CONS* vbcons = nlhdlrdata->vbdconss[c];
      if(vbcons->consdata == NULL) continue;
      int bpos, cpos;
      if( SCIPgetVarsLinear(scip, vbcons)[0] == var )
      {
         bpos = 1;
         cpos = 0;
      }
      else if( SCIPgetVarsLinear(scip, vbcons)[1] == var )
      {
         bpos = 0;
         cpos = 1;
      }
      else
         continue;


      *bvar = SCIPgetVarsLinear(scip, vbcons)[bpos]; /* TODO what if there are varbounds with different bvars? */

      SCIP_Bool leq0 = SCIPgetRhsLinear(scip, vbcons) == 0 && SCIPgetLhsLinear(scip, vbcons) == -SCIPinfinity(scip);
      SCIP_Bool geq0 = SCIPgetLhsLinear(scip, vbcons) == 0 && SCIPgetRhsLinear(scip, vbcons) == SCIPinfinity(scip);
      SCIP_Real ccoef = SCIPgetValsLinear(scip, vbcons)[cpos];
      /* TODO currently handling only one-sided constraints with 0 lhs or rhs, could be generalised later */
      if( (leq0 && ccoef > 0) || (geq0 && ccoef < 0) )
         pmax = -SCIPgetValsLinear(scip, vbcons)[bpos] / ccoef;
      else if( (geq0 && ccoef > 0) || (leq0 && ccoef < 0) )
         pmin = -SCIPgetValsLinear(scip, vbcons)[bpos] / ccoef;
      if( pmin != -SCIPinfinity(scip) && pmax != SCIPinfinity(scip) )
         break;
   }
   if( pmin == -SCIPinfinity(scip) || pmax == SCIPinfinity(scip) )
   {
      *bvar = NULL;
      return FALSE;
   }

   SCIPinfoMessage(scip, NULL, "\nFound a semicontinuous variable, pmin = %f, pmax = %f", pmin, pmax);
   return TRUE;
}
#endif


static
SCIP_RETCODE addRefterm(
   SCIP*                  scip,
   SCIP_CONSHDLR*         conshdlr,
   SCIP_CONSEXPR_EXPR*    expr,
   SCIP_Real              coef,
   SCIP_CONSEXPR_EXPR**   refexpr,
   SCIP_Bool              reformulate
   )
{
   SCIP_CONSEXPR_EXPR* refterm;
   if( reformulate )
   {
      /* TODO: reformulate */
#if 0
     refterm = reformulateTerm(expr);
#endif
      SCIPinfoMessage(scip, NULL, "\nExpr can be reformulated: ");
      SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
      refterm = expr;
   }
   else
      refterm = expr;
   if( *refexpr == NULL )
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, refexpr, 1, &refterm, &coef, 0.0) );
   else
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *refexpr, refterm, coef) );
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
   if(nlhdlrdata->semiContVars != NULL)
      SCIPhashmapFree(&(nlhdlrdata)->semiContVars);

   if(nlhdlrdata->oldSemiContVars != NULL)
   {
      for( int c = 0; c < SCIPhashmapGetNEntries(nlhdlrdata->oldSemiContVars); ++c )
      {
         SCIP_HASHMAPENTRY* entry;
         SCIP_SCVARDATA* data;
         entry = SCIPhashmapGetEntry(nlhdlrdata->oldSemiContVars, c);

         if( entry != NULL )
         {
            data = (SCIP_SCVARDATA*) SCIPhashmapEntryGetImage(entry);
            SCIPfreeBlockMemory(scip, &data);
         }
      }
      SCIPhashmapFree(&(nlhdlrdata)->oldSemiContVars);
   }

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


static
SCIP_DECL_CONSEXPR_NLHDLRUPDATE(nlhdlrUpdatePerspective)
{
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONS** conss;
   SCIP_HASHMAP* newSCvars;
   SCIP_VAR* x;
   SCIP_Bool found = FALSE;
   int c, nconss;

   SCIP_CONSHDLR* varboundconshdlr;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* look for important constraint handlers */
   varboundconshdlr = SCIPfindConshdlr(scip, "varbound");
   if( varboundconshdlr == NULL )
      return SCIP_OKAY; /* nothing to do if there are no varbounds */

   nconss = SCIPconshdlrGetNConss(varboundconshdlr);
   conss = SCIPconshdlrGetConss(varboundconshdlr);

   if( nlhdlrdata->semiContVars!= NULL )
      SCIPhashmapFree(&(nlhdlrdata)->semiContVars);

   SCIPhashmapCreate(&(nlhdlrdata)->semiContVars, SCIPblkmem(scip), nconss);
   newSCvars = nlhdlrdata->semiContVars;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
      x = SCIPgetVarVarbound(scip, conss[c]);
      if( SCIPgetLhsVarbound(scip, conss[c]) != 0 && SCIPgetRhsVarbound(scip, conss[c]) != 0 )
         continue;
      SCIP_Real xmin = SCIPgetLhsVarbound(scip, conss[c]) == 0 ? -SCIPgetVbdcoefVarbound(scip, conss[c]) : -SCIPinfinity(scip);
      SCIP_Real xmax = SCIPgetRhsVarbound(scip, conss[c]) == 0 ? -SCIPgetVbdcoefVarbound(scip, conss[c]) : SCIPinfinity(scip);
      SCIP_VAR* bvar = SCIPgetVbdvarVarbound(scip, conss[c]);

      SCIP_SCVARDATA* oldvardata;
      if( nlhdlrdata->oldSemiContVars == NULL )
         oldvardata = NULL;
      else
         oldvardata = (SCIP_SCVARDATA *) SCIPhashmapGetImage(nlhdlrdata->oldSemiContVars, x);

      /* the variable was already detected as semicontinuous at a previous round */
      if( oldvardata != NULL && oldvardata->xmin != -SCIPinfinity(scip) && oldvardata->xmax != SCIPinfinity(scip) )
         continue;

      /* the constraint does not add anything to existing semicontinuous var data */
      if(oldvardata != NULL && oldvardata->bvar == bvar && oldvardata->xmin == xmin && oldvardata->xmax == xmax)
         continue;

      if( !SCIPhashmapExists(nlhdlrdata->semiContVars, x) )
      {
         SCIP_SCVARDATA* newvardata;
         SCIP_CALL( SCIPallocBlockMemory(scip, &newvardata) );
         newvardata->bvar = bvar;
         newvardata->xmin = xmin;
         newvardata->xmax = xmax;
         SCIP_CALL( SCIPhashmapInsert(nlhdlrdata->semiContVars, x, newvardata) );
         found = TRUE;
      }
      else
      {
         SCIP_SCVARDATA* newvardata;
         newvardata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(nlhdlrdata->semiContVars, x);
         if( bvar != newvardata->bvar )
            continue;
         if( xmin > newvardata->xmin )
            newvardata->xmin = xmin;
         if( xmax < newvardata->xmax )
            newvardata->xmax = xmax;
      }
   }


   if( newSCvars != NULL )
   {
      int nentries = SCIPhashmapGetNEntries(newSCvars);
      for( c = 0; c < nentries; ++c )
      {
         SCIP_HASHMAPENTRY* entry;
         entry = SCIPhashmapGetEntry(newSCvars, c);

         if( entry != NULL )
         {
            SCIP_SCVARDATA* data = (SCIP_SCVARDATA*)SCIPhashmapEntryGetImage(entry);
            SCIP_VAR* origin = (SCIP_VAR*)SCIPhashmapEntryGetOrigin(entry);
            if( nlhdlrdata->oldSemiContVars == NULL )
               SCIPhashmapCreate(&(nlhdlrdata)->oldSemiContVars, SCIPblkmem(scip), nentries);
            SCIPhashmapInsert(nlhdlrdata->oldSemiContVars, origin, data);
         }
      }
   }

#if 0
   SCIPinfoMessage(scip, NULL, "\nNew semicontinuous variables: ");
   if( !found )
      return SCIP_OKAY;
   for( c = 0; c < SCIPhashmapGetNEntries(nlhdlrdata->semiContVars); ++c )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_VAR* newx;
      SCIP_SCVARDATA* newdata;
      entry = SCIPhashmapGetEntry(nlhdlrdata->semiContVars, c);

      if( entry != NULL )
      {
         newx = (SCIP_VAR*) SCIPhashmapEntryGetOrigin(entry);
         assert(newx != NULL);
         newdata = (SCIP_SCVARDATA*) SCIPhashmapEntryGetImage(entry);

         SCIPinfoMessage(scip, NULL, " \n%s: bvar = %s, xmin = %f, xmax = %f", SCIPvarGetName(newx), SCIPvarGetName(newdata->bvar), newdata->xmin, newdata->xmax);
         assert(SCIPhashmapExists(nlhdlrdata->semiContVars, newx));
      }
   }
   SCIPinfoMessage(scip, NULL, "\n");

   SCIPinfoMessage(scip, NULL, "\nOld semicontinuous variables: ");
   for( c = 0; c < SCIPhashmapGetNEntries(nlhdlrdata->oldSemiContVars); ++c )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_VAR* newx;
      SCIP_SCVARDATA* newdata;
      entry = SCIPhashmapGetEntry(nlhdlrdata->oldSemiContVars, c);

      if( entry != NULL )
      {
         newx = (SCIP_VAR*) SCIPhashmapEntryGetOrigin(entry);
         assert(newx != NULL);
         newdata = (SCIP_SCVARDATA*) SCIPhashmapEntryGetImage(entry);

         SCIPinfoMessage(scip, NULL, " \n%s: bvar = %s, xmin = %f, xmax = %f", SCIPvarGetName(newx), SCIPvarGetName(newdata->bvar), newdata->xmin, newdata->xmax);
      }
   }
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   return SCIP_OKAY;
}

/** nonlinear handler callback for reformulation */
static
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(nlhdlrReformulatePerspective)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_VAR* bvar;
   SCIP_CONSEXPR_EXPR* child;
   int c;

#if 0
   SCIPinfoMessage(scip, NULL, "\n------------------\nCalled for expr: ");
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
#endif

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   if(SCIPhashmapIsEmpty(nlhdlrdata->semiContVars) )
   {
      /* No new semicontinuous variables have been found by update */
      *refexpr = expr;
      SCIPcaptureConsExprExpr(*refexpr);
      return SCIP_OKAY;
   }

#if 0
   if( !nlhdlrdata->detected )
   {
      SCIP_CONSHDLR* varboundconshdlr;
      SCIP_CONSHDLR* linconshdlr;
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

      SCIPinfoMessage(scip, NULL, "\nVar bound constraints:\n");
      for( c = 0; c < nlhdlrdata->nvbdconss; c++ )
      {
         SCIP_CALL( SCIPprintCons(scip, nlhdlrdata->vbdconss[c], NULL) );
         SCIPinfoMessage(scip, NULL, "\n");
      }
            nlhdlrdata->detected = TRUE;
   }
#endif

   /* TODO detect structure (assumes that the terms are sorted) */

   SCIP_CONSEXPR_EXPR** children;
   int nchildren;

   if( strcmp("sum", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr))) == 0 )
   {
      children = SCIPgetConsExprExprChildren(expr);
      nchildren = SCIPgetConsExprExprNChildren(expr);
   }
   else
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &children, 0, 1) );
      *children = expr;
      nchildren = 1;
   }


   SCIP_VAR* qvar = NULL;
   int qvarpos = -1;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      child = SCIPgetConsExprExprChildren(expr)[c];

      if( strcmp("var", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) == 0)
      {
         if( varIsSemicontinuous(scip, nlhdlrdata, SCIPgetConsExprExprVarVar(child), &bvar) )
         {
            /* save the variable, do nothing */
            qvar = SCIPgetConsExprExprVarVar(child);
            qvarpos = c;
            continue;
         }
         else
            addRefterm(scip, conshdlr, child, SCIPgetConsExprExprSumCoefs(expr)[c], refexpr, FALSE);
      }
      else if( strcmp("pow", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) == 0 )
      {
         SCIP_CONSEXPR_EXPR* powbase;
         powbase = SCIPgetConsExprExprChildren(child)[0];
         if( strcmp("var", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(powbase))) != 0  )
            continue;
         SCIP_VAR* powvar;
         powvar = SCIPgetConsExprExprVarVar(powbase);
         if(SCIPvarGetType(powvar) != SCIP_VARTYPE_CONTINUOUS)
            continue;
         if( SCIPgetConsExprExprPowExponent(child) == 2) /* quadratic univariate */
         {
            if( qvar == powvar )
            {
               SCIP_CONSEXPR_EXPR* quadsum;
               SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &quadsum, 1, &SCIPgetConsExprExprChildren(expr)[qvarpos], &SCIPgetConsExprExprSumCoefs(expr)[qvarpos], 0.0) );
               SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, quadsum, child, SCIPgetConsExprExprSumCoefs(expr)[c]) );
               addRefterm(scip, conshdlr, quadsum, 1.0, refexpr, TRUE);
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &quadsum) );
            }
            else
            {
               if( varIsSemicontinuous(scip, nlhdlrdata, powvar, &bvar) )
                  addRefterm(scip, conshdlr, child, SCIPgetConsExprExprSumCoefs(expr)[c], refexpr, TRUE);
               else
                  addRefterm(scip, conshdlr, child, SCIPgetConsExprExprSumCoefs(expr)[c], refexpr, FALSE);
               if( qvar != NULL )
               {
                  addRefterm(scip, conshdlr, SCIPgetConsExprExprChildren(expr)[qvarpos],
                             SCIPgetConsExprExprSumCoefs(expr)[qvarpos], refexpr, FALSE);
                  qvar = NULL;
                  qvarpos = -1;
               }
            }

         }
         else /* univariate rational power (non-quadratic) */
         {
            if( varIsSemicontinuous(scip, nlhdlrdata, powvar, &bvar) )
               addRefterm(scip, conshdlr, child, SCIPgetConsExprExprSumCoefs(expr)[c], refexpr, TRUE);
            else
               addRefterm(scip, conshdlr, child, SCIPgetConsExprExprSumCoefs(expr)[c], refexpr, FALSE);
         }
      }
      else if( strcmp("exp", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child))) == 0 )
      {
         SCIP_CONSEXPR_EXPR* exppow;
         exppow = SCIPgetConsExprExprChildren(child)[0];
         if( strcmp("var", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(exppow))) != 0  )
            continue;
         SCIP_VAR* expvar;
         expvar = SCIPgetConsExprExprVarVar(exppow);
         if(SCIPvarGetType(expvar) != SCIP_VARTYPE_CONTINUOUS)
            continue;
         if( varIsSemicontinuous(scip, nlhdlrdata, expvar, &bvar) )
               addRefterm(scip, conshdlr, child, SCIPgetConsExprExprSumCoefs(expr)[c], refexpr, TRUE);
         else
               addRefterm(scip, conshdlr, child, SCIPgetConsExprExprSumCoefs(expr)[c], refexpr, FALSE);
      }
      else
         addRefterm(scip, conshdlr, child, SCIPgetConsExprExprSumCoefs(expr)[c], refexpr, FALSE);
   }


   if( strcmp("sum", SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr))) == 0 )
   {
      SCIPsetConsExprExprSumConstant(*refexpr, SCIPgetConsExprExprSumConstant(expr));
   }
   else
   {
      SCIPfreeBlockMemoryArrayNull(scip, &children, nchildren);
   }


   SCIPinfoMessage(scip, NULL, "\n\nrefexpr:\n");
   SCIPprintConsExprExpr(scip, conshdlr, *refexpr, NULL);
   SCIPinfoMessage(scip, NULL, "\n\n");


   /* TODO create expression and store it in refexpr */

   if( *refexpr != NULL )
   {
      SCIPreleaseConsExprExpr(scip,refexpr);
      *refexpr = NULL;
   }

   /* set refexpr to expr and capture it if no reformulation is possible */
   if(*refexpr == NULL)
   {
      *refexpr = expr;
      SCIPcaptureConsExprExpr(*refexpr);
   }

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
   SCIPsetConsExprNlhdlrUpdate(scip, nlhdlr, nlhdlrUpdatePerspective);

   return SCIP_OKAY;
}
