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

/**@file   cons_expr_nlhdlr_perspective.c
 * @brief  perspective nonlinear handler
 * @author Ksenia Bestuzheva
 */

#include <string.h>

#include "scip/cons_varbound.h"
#include "scip/cons_expr_nlhdlr_perspective.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/scip_sol.h"
#include "scip/cons_expr_iterator.h"
#include "struct_cons_expr.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "perspective"
#define NLHDLR_DESC               "perspective handler for expressions"
#define NLHDLR_DETECTPRIORITY     -20 /* detect last so that to make use of what other handlers detected */
#define NLHDLR_ENFOPRIORITY       125 /* enforce first because perspective cuts are always stronger */

#define DEFAULT_DETECTSUM    FALSE
#define DEFAULT_MULTCUTS     TRUE

/*
 * Data structures
 */

/** data structure to store information of a semicontinuous variable */
struct SCVarData
{
   SCIP_Real*            vals0;              /**< values of the variable when the corresponding bvars[i] = 0 */
   SCIP_Real*            lbs1;               /**< lower bounds of the variable when the corresponding bvars[i] = 0 */
   SCIP_Real*            ubs1;               /**< upper bounds of the variable when the corresponding bvars[i] = 0 */
   SCIP_VAR**            bvars;              /**< the binary variables on which the variable domain depends */
   int                   nbnds;              /**< number of suitable on/off bounds the var has */
   int                   bndssize;           /**< size of the arrays */
};
typedef struct SCVarData SCVARDATA;

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
#ifdef SCIP_DISABLED_CODE
   SCIP_EXPRCURV         curvature;          /**< curvature of the expression */

   SCIP_CONSEXPR_EXPR**  onoffterms;         /**< on/off terms for which we apply perspective cuts */
   SCIP_Real*            onoffcoefs;         /**< coefficients of onoffterms */
   int                   nonoffterms;        /**< number of on/off expressions */
   int                   onofftermssize;     /**< size of arrays describing on/off terms */
   SCIP_VAR***           termbvars;          /**< binary vars associated with onoffterms */
   int*                  ntermbvars;         /**< number of binary variables for each term */

   SCIP_CONSEXPR_EXPR**  convterms;          /**< convex terms for which we apply gradient cuts */
   SCIP_Real*            convcoefs;          /**< coefficients of convterms */
   int                   nconvterms;         /**< number of convterms */
   int                   convtermssize;      /**< size of the convterms array */
#endif

   SCIP_Real*            exprvals0;          /**< 'off' values of the expression for each indicator variable */

   SCIP_CONSEXPR_EXPR**  varexprs;           /**< variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */

   SCIP_VAR**            indicators;         /**< all indicator variables for the expression */
   int                   nindicators;        /**< number of indicator variables */
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_Bool             detectsum;          /**< whether to run detection when the root of an expression is a sum */
   SCIP_Bool             multcuts;           /**< whether to add cuts for all suitable indicator variables */ /* TODO do we need this? */
   SCIP_HASHMAP*         scvars;             /**< maps semicontinuous variables to their on/off bounds */
};

/*
 * Local methods
 */

/** adds an indicator to the data of a semicontinuous variable */
static
SCIP_RETCODE addSCVarIndicator(
   SCIP*                 scip,               /**< SCIP data structure */
   SCVARDATA*            scvdata,            /**< semicontinuous variable data */
   SCIP_VAR*             indicator,          /**< indicator to be added */
   SCIP_Real             val0,               /**< value of the variable when indicator = 0 */
   SCIP_Real             lb1,                /**< lower bound of the variable when indicator = 1 */
   SCIP_Real             ub1                 /**< upper bound of the variable when indicator = 1 */
   )
{
   int newsize;
   int i;
   SCIP_Bool found;
   int pos;

   assert(scvdata != NULL);
   assert(indicator != NULL);

   /* ensure sizes */
   if( scvdata->nbnds + 1 > scvdata->bndssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, scvdata->nbnds + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->bvars, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->vals0, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->lbs1, scvdata->bndssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->ubs1, scvdata->bndssize, newsize) );
      scvdata->bndssize = newsize;
   }
   assert(scvdata->nbnds + 1 <= scvdata->bndssize);

   /* find the position where to insert */
   found = SCIPsortedvecFindPtr((void**)scvdata->bvars, SCIPvarComp, (void*)indicator, scvdata->nbnds, &pos);

   if( found )
      return SCIP_OKAY;

   /* move entries if needed */
   for( i = scvdata->nbnds; i > pos; --i )
   {
      scvdata->bvars[i] = scvdata->bvars[i-1];
      scvdata->vals0[i] = scvdata->vals0[i-1];
      scvdata->lbs1[i] = scvdata->lbs1[i-1];
      scvdata->ubs1[i] = scvdata->ubs1[i-1];
   }

   scvdata->bvars[pos] = indicator;
   scvdata->vals0[pos] = val0;
   scvdata->lbs1[pos] = lb1;
   scvdata->ubs1[pos] = ub1;
   ++scvdata->nbnds;

   return SCIP_OKAY;
}

/** checks if a variable is semicontinuous and, if needed, updates the hashmap
 *
 * A variable is semicontinuous if its bounds depend on the binary variable bvar and bvar == 0 => var = v_off for some
 * real constant v_off. If the bvar is not specified, find the first binary variable that var depends on.
 */
static
SCIP_RETCODE varIsSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< the variable to check */
   SCIP_HASHMAP*         scvars,             /**< semicontinuous variable information */
   SCIP_VAR*             indicator,          /**< indicator variable which var should depend on (NULL if doesn't matter) */
   SCIP_Real*            val0,               /**< buffer to store value of var when indicator == 0 (NULL if not interested) */
   SCIP_Bool*            result              /**< buffer to store whether var is semicontinuous */
   )
{
   SCIP_Real lb0;
   SCIP_Real ub0;
   SCIP_Real lb1;
   SCIP_Real ub1;
   SCIP_Real glb;
   SCIP_Real gub;
   SCIP_Bool exists;
   int c;
   int pos;
   SCIP_VAR** vlbvars;
   SCIP_VAR** vubvars;
   SCIP_Real* vlbcoefs;
   SCIP_Real* vubcoefs;
   SCIP_Real* vlbconstants;
   SCIP_Real* vubconstants;
   int nvlbs;
   int nvubs;
   SCVARDATA* scvdata;
   SCIP_VAR* bvar;

   assert(scip != NULL);
   assert(var != NULL);
   assert(scvars != NULL);
   assert(result != NULL);

   scvdata = (SCVARDATA*) SCIPhashmapGetImage(scvars, (void*)var);
   if( scvdata != NULL )
   {
      *result = TRUE;

      if( indicator != NULL )
      { /* if the indicator variable matters, look for it */
         exists = SCIPsortedvecFindPtr((void**)scvdata->bvars, SCIPvarComp, (void*)indicator, scvdata->nbnds, &pos);
         if( !exists )
            *result = FALSE;
         else if( val0 != NULL )
            *val0 = scvdata->vals0[pos];
      }

      return SCIP_OKAY;
   }

   vlbvars = SCIPvarGetVlbVars(var);
   vubvars = SCIPvarGetVubVars(var);
   vlbcoefs = SCIPvarGetVlbCoefs(var);
   vubcoefs = SCIPvarGetVubCoefs(var);
   vlbconstants = SCIPvarGetVlbConstants(var);
   vubconstants = SCIPvarGetVubConstants(var);
   nvlbs = SCIPvarGetNVlbs(var);
   nvubs = SCIPvarGetNVubs(var);
   glb = SCIPvarGetLbGlobal(var);
   gub = SCIPvarGetUbGlobal(var);

   *result = FALSE;

   /* Scan through lower bounds; for each binary vlbvar save the corresponding lb0 and lb1.
    * Then check if there is an upper bound with this vlbvar and save ub0 and ub1.
    * If the found bounds imply that the var value is fixed to some val0 when vlbvar = 0,
    * save vlbvar and val0 to scvdata.
    */
   for( c = 0; c < nvlbs; ++c )
   {
      if( SCIPvarGetType(vlbvars[c]) != SCIP_VARTYPE_BINARY )
         continue;

      SCIPdebugMsg(scip, "var <%s>[%f, %f] lower bound: %f <%s> %+f", SCIPvarGetName(var), glb, gub, vlbcoefs[c], SCIPvarGetName(vlbvars[c]), vlbconstants[c]);

      bvar = vlbvars[c];

      lb0 = MAX(vlbconstants[c], glb);
      lb1 = MAX(vlbconstants[c] + vlbcoefs[c], glb);

      /* look for bvar in vubvars */
      if( vubvars != NULL )
         exists = SCIPsortedvecFindPtr((void**)vubvars, SCIPvarComp, bvar, nvubs, &pos);
      else
         exists = FALSE;
      if( exists )
      { /*lint --e{644}*/
         SCIPdebugMsgPrint(scip, ", upper bound: %f <%s> %+f", vubcoefs[pos], SCIPvarGetName(vubvars[pos]), vubconstants[pos]); /*lint !e613*/

         /* save the upper bounds */
         ub0 = MIN(vubconstants[pos], gub);
         ub1 = MIN(vubconstants[pos] + vubcoefs[pos], gub);
      }
      else
      {
         /* if there is no upper bound with vubvar = bvar, use global var bounds */
         ub0 = gub;
         ub1 = gub;
      }

      /* the 'off' domain of a semicontinuous var should reduce to a single point and be different from the 'on' domain */
      SCIPdebugMsgPrint(scip, " -> <%s> in [%f, %f] (off), [%f, %f] (on)\n", SCIPvarGetName(var), lb0, ub0, lb1, ub1);
      if( lb0 == ub0 && (lb0 != lb1 || ub0 != ub1) ) /*lint !e777*/
      {
         if( scvdata == NULL )
         {
            SCIP_CALL( SCIPallocClearBlockMemory(scip, &scvdata) );
         }

         addSCVarIndicator(scip, scvdata, bvar, lb0, lb1, ub1);
      }
   }

   /* look for vubvars that have not been processed yet */
   assert(vubvars != NULL || nvubs == 0);
   for( c = 0; c < nvubs; ++c )
   {
      if( SCIPvarGetType(vubvars[c]) != SCIP_VARTYPE_BINARY)  /*lint !e613*/
         continue;

      bvar = vubvars[c];  /*lint !e613*/

      /* skip vars that are in vlbvars */
      if( vlbvars != NULL && SCIPsortedvecFindPtr((void**)vlbvars, SCIPvarComp, bvar, nvlbs, &pos) )
         continue;

      SCIPdebugMsg(scip, "var <%s>[%f, %f] upper bound: %f <%s> %+f",
         SCIPvarGetName(var), glb, gub, vubcoefs[c], SCIPvarGetName(vubvars[c]), vubconstants[c]);  /*lint !e613*/

      lb0 = glb;
      lb1 = glb;
      ub0 = MIN(vubconstants[c], gub);
      ub1 = MIN(vubconstants[c] + vubcoefs[c], gub);

      /* the 'off' domain of a semicontinuous var should reduce to a single point and be different from the 'on' domain */
      SCIPdebugMsgPrint(scip, " -> <%s> in [%f, %f] (off), [%f, %f] (on)\n", SCIPvarGetName(var), lb0, ub0, lb1, ub1);
      if( lb0 == ub0 && (lb0 != lb1 || ub0 != ub1) ) /*lint !e777*/
      {
         if( scvdata == NULL )
         {
            SCIP_CALL( SCIPallocClearBlockMemory(scip, &scvdata) );
         }

         addSCVarIndicator(scip, scvdata, bvar, lb0, lb1, ub1);
      }
   }

   if( scvdata != NULL )
   {
#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "var <%s> has global bounds [%f, %f] and the following on/off bounds:\n", SCIPvarGetName(var), glb, gub);
      for( c = 0; c < scvdata->nbnds; ++c )
      {
         SCIPdebugMsg(scip, " c = %d, bvar <%s>: val0 = %f\n", c, SCIPvarGetName(scvdata->bvars[c]), scvdata->vals0[c]);
      }
#endif
      SCIP_CALL( SCIPhashmapInsert(scvars, var, scvdata) );
      *result = TRUE;
   }

   return SCIP_OKAY;
}

/* TODO keeping the unused functions for now, will decide if they are needed later */
#ifdef SCIP_DISABLED_CODE
/** adds an expression to the array of on/off expressions */
static
SCIP_RETCODE addOnoffTerm(
   SCIP*                         scip,                     /**< SCIP data structure */
   SCIP_CONSHDLR*                conshdlr,                 /**< constraint handler */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata,           /**< expression data */
   SCIP_Real                     coef,                     /**< coef of the added term */
   SCIP_CONSEXPR_EXPR*           expr,                     /**< expr to add */
   SCIP_VAR**                    bvars,                    /**< binary variables */
   int                           nbvars                    /**< number of binary variables */
)
{
   int newsize;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(bvars != NULL);

   if( nlhdlrexprdata->nonoffterms + 1 > nlhdlrexprdata->onofftermssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, nlhdlrexprdata->nonoffterms + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->onoffterms,  nlhdlrexprdata->onofftermssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->onoffcoefs, nlhdlrexprdata->onofftermssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->termbvars, nlhdlrexprdata->onofftermssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->ntermbvars, nlhdlrexprdata->onofftermssize, newsize) );
      nlhdlrexprdata->onofftermssize = newsize;
   }
   assert(nlhdlrexprdata->nonoffterms + 1 <= nlhdlrexprdata->onofftermssize);

   nlhdlrexprdata->onoffcoefs[nlhdlrexprdata->nonoffterms] = coef;
   nlhdlrexprdata->onoffterms[nlhdlrexprdata->nonoffterms] = expr;
   nlhdlrexprdata->termbvars[nlhdlrexprdata->nonoffterms] = bvars;
   nlhdlrexprdata->ntermbvars[nlhdlrexprdata->nonoffterms] = nbvars;
   nlhdlrexprdata->nonoffterms++;

   return SCIP_OKAY;
}

/** adds an expression to the array of convex expressions */
static
SCIP_RETCODE addConvTerm(
   SCIP*                           scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA*   nlhdlrexprdata,   /**< nonlinear handler expression data */
   SCIP_Real                       coef,             /**< coefficient of expr in the original expression */
   SCIP_CONSEXPR_EXPR*             expr              /**< expression to be added */
)
{
   int newsize;

   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(expr != NULL);

   if( nlhdlrexprdata->nconvterms + 1 > nlhdlrexprdata->convtermssize )
   {
      newsize = SCIPcalcMemGrowSize(scip, nlhdlrexprdata->nconvterms + 1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->convterms,  nlhdlrexprdata->convtermssize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &nlhdlrexprdata->convcoefs, nlhdlrexprdata->convtermssize, newsize) );
      nlhdlrexprdata->convtermssize = newsize;
   }
   assert(nlhdlrexprdata->nconvterms + 1 <= nlhdlrexprdata->convtermssize);

   nlhdlrexprdata->convcoefs[nlhdlrexprdata->nconvterms] = coef;
   nlhdlrexprdata->convterms[nlhdlrexprdata->nconvterms] = expr;
   nlhdlrexprdata->nconvterms++;

   return SCIP_OKAY;
}

/** constructs gradient linearization of a given expression and adds it to rowprep */
static
SCIP_RETCODE addGradientLinearisation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_ROWPREP*         rowprep,            /**< a rowprep where the linearization is stored */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be linearized */
   SCIP_Real             coef,               /**< coefficient of expr in the original expression */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_Bool*            success             /**< indicates whether the linearization could be computed */
)
{
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_Real constant;
   int i, v, nvars;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(rowprep != NULL);
   assert(expr != NULL);
   assert(success != NULL);

   /* compute gradient */
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );

   /* gradient evaluation error -> skip */
   if( SCIPgetConsExprExprDerivative(expr) == SCIP_INVALID ) /*lint !e777*/
   {
      *success = FALSE;
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* get g(x*) */
   constant = SCIPgetConsExprExprValue(expr);

   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(constant)) )
   {
      *success = FALSE;
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", constant, (void*)expr);
      return SCIP_OKAY;
   }

   /* compute gradient cut */
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, SCIPgetNTotalVars(scip)) );
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, varexprs, &nvars) );
   for( i = 0; i < nvars; ++i )
   {
      SCIP_VAR* var;
      SCIP_Real derivative;
      SCIP_Real val;

      assert(varexprs[i] != NULL);
      assert(SCIPisConsExprExprVar(varexprs[i]));

      /* get the variable of the variable expression */
      var = SCIPgetConsExprExprVarVar(varexprs[i]);
      assert(var != NULL);

      /* get solution value */
      val = SCIPgetSolVal(scip, sol, var);

      /* avoid overhead of SCIPgetConsExprExprPartialDiff by accessing the derivative directly */
      derivative = SCIPgetConsExprExprDerivative(varexprs[i]);
      assert(SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var) == derivative); /*lint !e777*/

      /* evaluation error or too large values -> skip */
      if( SCIPisInfinity(scip, REALABS(derivative * val)) )
      {
         *success = FALSE;
         SCIPdebugMsg(scip, "evaluation error / too large values (%g %g) for %s in %p\n", derivative, val,
                      SCIPvarGetName(var), (void*)expr);
         goto TERMINATE;
      }

      /* - grad(g(x*))_i x*_i */
      constant -= derivative * val;

      /* grad(g(x*))_i x_i */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef*derivative) );
   }

   /* add constant */
   SCIPaddRowprepConstant(rowprep, coef*constant);

 TERMINATE:
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
   }
   SCIPfreeBufferArray(scip, &varexprs);

   return SCIP_OKAY;
}

/** constructs perspective linearization of a given expression and adds it to rowprep */
static
SCIP_RETCODE addPerspectiveLinearisation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_HASHMAP*         scvars,             /**< hashmap linking semicontinuous vars to their data */
   SCIP_ROWPREP*         rowprep,            /**< a rowprep where the linearization is stored */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be linearized */
   SCIP_Real             coef,               /**< coefficient of expr */
   SCIP_VAR*             bvar,               /**< binary variable */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_Bool*            success             /**< indicates whether the linearization could be computed */
   )
{
   SCIP_SOL* sol0;
   SCIP_Real* vals0;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_VAR** vars;
   SCIP_Real scalar_prod, fval, fval0;
   int nvars, v, pos;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(scvars != NULL);
   assert(rowprep != NULL);
   assert(expr != NULL);
   assert(bvar != NULL);
   assert(success != NULL);

   /* add the cut: auxvar >= (x - x0) \nabla f(sol) + (f(sol) - f(x0) - (sol - x0) \nabla f(sol)) z + f(x0),
    * where x is semicontinuous, z is binary and x0 is the value of x when z = 0 */

   SCIP_CALL( SCIPcreateSol(scip, &sol0, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, SCIPgetNTotalVars(scip)) );
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, varexprs, &nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals0, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "bvar = \n");
   SCIPprintVar(scip, bvar, NULL);
   SCIPdebugMsg(scip, NULL, "pexpr = \n");
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
#endif

   /* get x0 */
   for( v = 0; v < nvars; ++v )
   {
      SCVARDATA* vardata;
      vars[v] = SCIPgetConsExprExprVarVar(varexprs[v]);
      vardata = (SCVARDATA*)SCIPhashmapGetImage(scvars, (void*)vars[v]);

      /* find bvar in vardata->bvars */
      (void) SCIPsortedvecFindPtr((void**)vardata->bvars, SCIPvarComp, (void*)bvar, vardata->nbnds, &pos);
      assert(pos < vardata->nbnds);
      assert(vardata->bvars[pos] == bvar);

      vals0[v] = vardata->vals0[pos];
   }

   /* set x to x0 in sol0 */
   SCIP_CALL( SCIPsetSolVals(scip, sol0, nvars, vars, vals0) );

   /* get f(x0) */
   SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr, sol0, 0) );
   fval0 = SCIPgetConsExprExprValue(expr);
   SCIP_CALL( SCIPfreeSol(scip, &sol0) );

   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(fval0)) )
   {
      *success = FALSE;
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", fval0, (void*)expr);
      goto TERMINATE;
   }

   /* TODO it should not be necessary to reevaluate in sol, cons_expr should have done that already */
   /* get f(sol) */
   SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr, sol, 0) );
   fval = SCIPgetConsExprExprValue(expr);

   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(fval)) )
   {
      *success = FALSE;
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", fval, (void*)expr);
      goto TERMINATE;
   }

   /* add (f(sol) - f(x0))z + f(x0) */
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, bvar, coef*(fval - fval0)) );
   SCIPaddRowprepConstant(rowprep, coef*fval0);

   /* compute gradient */
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, expr, sol, 0) );

   /* gradient evaluation error -> skip */
   if( SCIPgetConsExprExprDerivative(expr) == SCIP_INVALID ) /*lint !e777*/
   {
      *success = FALSE;
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)expr);
      goto TERMINATE;
   }

   scalar_prod = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      var = SCIPgetConsExprExprVarVar(varexprs[v]);

      /* add xi f'xi(sol) */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, coef*SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var)) );
      /* add -x0i f'xi(sol) */
      SCIPaddRowprepConstant(rowprep, -coef*vals0[v]*SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var));

      /* compute -(soli - x0i) f'xi(sol) */
      scalar_prod -= (SCIPgetSolVal(scip, sol, var) - vals0[v])*SCIPgetConsExprExprPartialDiff(scip, conshdlr, expr, var);
   }

   /* add -(sol - x0) \nabla f(sol)) z */
   SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, bvar, coef*scalar_prod) );

 TERMINATE:
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &vals0);
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
   }
   SCIPfreeBufferArray(scip, &varexprs);

   return SCIP_OKAY;
}
#endif

/* checks if an expression is semicontinuous
 *
 * An expression is semicontinuous if all of its variables are semicontinuous
 * and share at least one common indicator variable
 */
static
SCIP_RETCODE exprIsSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata,     /**< nonlinear handler data */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nlhdlr expression data */
   SCIP_Bool*            res                 /**< buffer to store whether the expression is semicontinuous */
   )
{
   int v;
   int nvarexprs;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_Bool var_is_sc;
   SCVARDATA* scvdata;
   SCIP_VAR* var;
   int nindicators;
   int nbnds0;
   SCIP_VAR** indicators;

   *res = FALSE;

   nvarexprs = nlhdlrexprdata->nvarexprs;
   varexprs = nlhdlrexprdata->varexprs;

   /* constant expression is not semicontinuous */
   if( nvarexprs == 0 )
   {
      return SCIP_OKAY;
   }

   /* all variables of an on/off term should be semicontinuous */
   for( v = 0; v < nvarexprs; ++v )
   {
      var = SCIPgetConsExprExprVarVar(varexprs[v]);
      SCIP_CALL( varIsSemicontinuous(scip, var, nlhdlrdata->scvars, NULL, NULL, &var_is_sc) );
      if( !var_is_sc )
         return SCIP_OKAY;
   }

   /* look for common binary variables for all variables of the expression */
   scvdata = (SCVARDATA*)SCIPhashmapGetImage(nlhdlrdata->scvars, (void*)SCIPgetConsExprExprVarVar(varexprs[0]));
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &indicators, scvdata->bvars, scvdata->nbnds) );

   nbnds0 = scvdata->nbnds;
   nindicators = nbnds0;

   SCIPdebugMsg(scip, "\nArray intersection for vars %s, *nbvars = %d", SCIPvarGetName(SCIPgetConsExprExprVarVar(varexprs[0])), nindicators);
   for( v = 1; v < nvarexprs; ++v )
   {
#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "\n%s; ", SCIPvarGetName(SCIPgetConsExprExprVarVar(varexprs[v])));
#endif
      scvdata = (SCVARDATA*)SCIPhashmapGetImage(nlhdlrdata->scvars, (void*)SCIPgetConsExprExprVarVar(varexprs[v]));

      SCIPcomputeArraysIntersectionPtr((void**)indicators, nindicators, (void**)scvdata->bvars, scvdata->nbnds, SCIPvarComp, (void**)indicators, &nindicators);

      /* if we have found out that the intersection is empty, expr is not semicontinuous */
      if( nindicators == 0 )
      {
         SCIPfreeBlockMemoryArray(scip, &indicators, nbnds0);
         return SCIP_OKAY;
      }
   }

   assert(nindicators > 0 && nindicators <= nbnds0);

   if( nindicators < nbnds0 )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &indicators, nbnds0, nindicators) );
   }

   nlhdlrexprdata->indicators = indicators;
   nlhdlrexprdata->nindicators = nindicators;
   *res = TRUE;

   return SCIP_OKAY;
}

/** add the cut given by rowprep to sepastore */
static
SCIP_RETCODE addCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< expression constraint */
   SCIP_ROWPREP*         rowprep,            /**< cut to be added */
   SCIP_SOL*             sol,                /**< solution to be separated */
   SCIP_RESULT*          result              /**< pointer to store result */
   )
{
   SCIP_Bool success;
   SCIP_ROW* row;

   /* merge coefficients that belong to same variable */
   SCIPmergeRowprepTerms(scip, rowprep);

   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, SCIPgetLPFeastol(scip), NULL, &success) );

   /* if cut looks good (numerics ok and cutting off solution), then turn into row and add to sepastore */
   if( success )
   {
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

      if( infeasible )
      {
         *result = SCIP_CUTOFF;
      }
      else
      {
         *result = SCIP_SEPARATED;
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** frees nlhdlrexprdata structure */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata /**< nlhdlr expression data */
   )
{
   int c;

   if( nlhdlrexprdata->nindicators != 0 )
   {
      assert(nlhdlrexprdata->indicators != NULL);
      SCIPfreeBlockMemoryArray(scip, &(nlhdlrexprdata->indicators), nlhdlrexprdata->nindicators);
      SCIPfreeBlockMemoryArray(scip, &(nlhdlrexprdata->exprvals0), nlhdlrexprdata->nindicators);
   }

   if( nlhdlrexprdata->varexprs != NULL )
   {
      for( c = 0; c < nlhdlrexprdata->nvarexprs; ++c )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(nlhdlrexprdata->varexprs[c])) );
      }
      SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->varexprs, nlhdlrexprdata->nvarexprs);
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

/** nonlinear handler copy callback */
static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrPerspective)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrPerspective(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}


/** callback to free data of handler */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrFreehdlrdataPerspective)
{ /*lint --e{715}*/
   SCIP_HASHMAPENTRY* entry;
   SCVARDATA* data;
   int c;

   if( (*nlhdlrdata)->scvars != NULL )
   {
      for( c = 0; c < SCIPhashmapGetNEntries((*nlhdlrdata)->scvars); ++c )
      {
         entry = SCIPhashmapGetEntry((*nlhdlrdata)->scvars, c);
         if( entry != NULL )
         {
            data = (SCVARDATA*) SCIPhashmapEntryGetImage(entry);
            SCIPfreeBlockMemoryArray(scip, &data->ubs1, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->lbs1, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->vals0, data->bndssize);
            SCIPfreeBlockMemoryArray(scip, &data->bvars, data->bndssize);
            SCIPfreeBlockMemory(scip, &data);
         }
      }
      SCIPhashmapFree(&((*nlhdlrdata)->scvars));
   }

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}


/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrFreeExprDataPerspective)
{  /*lint --e{715}*/
   SCIP_CALL( freeNlhdlrExprData(scip, *nlhdlrexprdata) );
   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}


/** callback to be called in initialization */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINIT(nlhdlrInitPerspective)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}
#endif


/** callback to be called in deinitialization */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLREXIT(nlhdlrExitPerspective)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}
#endif

/* remove an indicator from nonlinear expression data */
static
void removeExprIndicator(
  SCIP_CONSEXPR_NLHDLREXPRDATA* nlexprdata,  /**< nonlinear expression data */
  int                    pos                 /**< position of the indicator */
  )
{
   int i;

   assert(pos >= 0 && pos < nlexprdata->nindicators);

   for( i = pos; i < nlexprdata->nindicators - 1; ++i )
   {
      nlexprdata->indicators[i] = nlexprdata->indicators[i+1];
   }

   --nlexprdata->nindicators;
}

/** computes the 'off' value of the expression and the 'off' values of
  * semicontinuous auxiliary variables for each indicator variable
  */
static
SCIP_RETCODE computeOffValues(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_NLHDLRDATA* hdlrdata,       /**< nonlinear handler data */
   SCIP_CONSEXPR_NLHDLREXPRDATA* exprdata,   /**< nonlinear expression data */
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   SCIP_CONSEXPR_ITERATOR* it;
   SCIP_SOL* sol;
   int i;
   int v;
   SCIP_Real* vals0;
   SCIP_VAR** vars;
   SCIP_Bool var_is_sc;
   SCVARDATA* scvdata;
   SCIP_VAR* auxvar;
   SCIP_CONSEXPR_EXPR* curexpr;

   assert(expr != NULL);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(exprdata->exprvals0), exprdata->nindicators) );

   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, exprdata->nvarexprs) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals0, exprdata->nvarexprs) );
   SCIP_CALL( SCIPexpriteratorCreate(&it, conshdlr, SCIPblkmem(scip)) );

   for( v = 0; v < exprdata->nvarexprs; ++v )
   {
      vars[v] = SCIPgetConsExprExprVarVar(exprdata->varexprs[v]);
   }

   for( i = 0; i < exprdata->nindicators; ++i )
   {
      /* set sol to the off value of all expr vars for this indicator */
      for( v = 0; v < exprdata->nvarexprs; ++v )
      {
         SCIP_CALL( varIsSemicontinuous(scip, vars[v], hdlrdata->scvars, exprdata->indicators[i], &vals0[v], &var_is_sc) );
         assert(var_is_sc);
      }
      SCIPsetSolVals(scip, sol, exprdata->nvarexprs, vars, vals0);
      SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, expr, sol, 0) );

      if( SCIPgetConsExprExprValue(expr) == SCIP_INVALID )
      {
         SCIPdebugMsg(scip, "expression evaluation failed for %p, removing the indicator\n", (void*)expr);
         removeExprIndicator(exprdata, i);
         continue;
      }

      exprdata->exprvals0[i] = SCIPgetConsExprExprValue(expr);

      /* iterate through the expression and create scvdata for all aux vars */
      SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_CONSEXPRITERATOR_DFS, FALSE) );
      curexpr = SCIPexpriteratorGetCurrent(it);

      while( !SCIPexpriteratorIsEnd(it) )
      {
         if( curexpr->auxvar != NULL )
         {
            /* we know that all vars are sc with respect to exprdata->indicators; it remains to:
             * - get or create the scvdata structure
             * - add it to scvars hashmap
             * - find the expr's off value
             * - add the indicator and off value to scvdata
             */
            auxvar = SCIPgetConsExprExprAuxVar(curexpr);

            scvdata = (SCVARDATA*) SCIPhashmapGetImage(hdlrdata->scvars, (void*)auxvar);
            if( scvdata == NULL )
            {
               SCIP_CALL( SCIPallocClearBlockMemory(scip, &scvdata) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scvdata->bvars,  exprdata->nindicators) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scvdata->vals0, exprdata->nindicators) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scvdata->lbs1, exprdata->nindicators) );
               SCIP_CALL( SCIPallocBlockMemoryArray(scip, &scvdata->ubs1, exprdata->nindicators) );
               scvdata->bndssize = exprdata->nindicators;
               SCIP_CALL( SCIPhashmapInsert(hdlrdata->scvars, auxvar, scvdata) );
            }

            SCIP_CALL( addSCVarIndicator(scip, scvdata, exprdata->indicators[i], SCIPgetConsExprExprValue(curexpr),
                     SCIPvarGetLbGlobal(auxvar), SCIPvarGetUbGlobal(auxvar)) );
         }

         curexpr = SCIPexpriteratorGetNext(it);
      }
   }

   SCIPexpriteratorFree(&it);
   SCIPfreeBufferArray(scip, &vals0);
   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeSol(scip, &sol);

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree
 *
 *  We are looking for expressions g(x), where x is a vector of semicontinuous variables that all share at least one
 *  indicator variable.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectPerspective)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrdata != NULL);

   *success = FALSE;

   /* do not run in presolve, as we only do separation */
   if( SCIPgetStage(scip) <= SCIP_STAGE_INITSOLVE )
   {
      return SCIP_OKAY;
   }

   if( SCIPgetNBinVars(scip) == 0 )
   {
      SCIPdebugMsg(scip, "\nproblem has no binary variables, not running perspective detection");
      return SCIP_OKAY;
   }

   /* some other nonlinear handler should be able to separate */
   if( !(*enforcemethods & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) && !(*enforcemethods & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) )
   {
      SCIPdebugMsg(scip, "\nno enforcement method, exiting detect");
      return SCIP_OKAY;
   }

#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "Called perspective detect, expr = %p: ", expr);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPdebugMsgPrint(scip, "\n");
#endif

   /* allocate memory */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   if( nlhdlrdata->scvars == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&(nlhdlrdata->scvars), SCIPblkmem(scip), SCIPgetNVars(scip)) );
   }

   /* save varexprs to nlhdlrexprdata */
   SCIP_CALL( SCIPgetConsExprExprNVars(scip, conshdlr, expr, &(*nlhdlrexprdata)->nvarexprs) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->varexprs, (*nlhdlrexprdata)->nvarexprs) );
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, (*nlhdlrexprdata)->varexprs, &(*nlhdlrexprdata)->nvarexprs) );

   /* check if expr is semicontinuous and save indicator variables */
   SCIP_CALL( exprIsSemicontinuous(scip, nlhdlrdata, *nlhdlrexprdata, success) );

   if( *success )
   {
      int sindicators;

      sindicators = (*nlhdlrexprdata)->nindicators;

      /* compute 'off' values and propagate semicontinuity */
      SCIP_CALL( computeOffValues(scip, conshdlr, nlhdlrdata, *nlhdlrexprdata, expr) );

      /* some indicator variables might have been removed if evaluation failed, check how many remain */
      if( (*nlhdlrexprdata)->nindicators == 0 )
      {
         *success = FALSE;
      }
      else if( (*nlhdlrexprdata)->nindicators < sindicators )
      {
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->indicators, sindicators,
               (*nlhdlrexprdata)->nindicators) );
         SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->exprvals0, sindicators,
                                                (*nlhdlrexprdata)->nindicators) );
      }
   }

   if( *success )
   {
      assert((*nlhdlrexprdata)->nindicators > 0);

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "\ndetected an on/off expr: ");
      SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
#endif

      /* if we get here, enforcemethods should have already been set by other handler(s) */
      assert(((*enforcemethods & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) && *enforcedbelow)
         || ((*enforcemethods & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) && *enforcedabove));

      assert(*nlhdlrexprdata != NULL);
   }
   else
   {
      SCIP_CALL( freeNlhdlrExprData(scip, *nlhdlrexprdata) );
      SCIPfreeBlockMemory(scip, nlhdlrexprdata);
      *nlhdlrexprdata = NULL;
   }

   return SCIP_OKAY;
}


/** auxiliary evaluation callback of nonlinear handler */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalauxPerspective)
{ /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(auxvalue != NULL);

   *auxvalue = SCIPgetConsExprExprValue(expr);

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
#endif


/** nonlinear handler estimation callback
 *
 * "Perspectivies" cuts produced by other handlers. Suppose that we want to separate x from the set g(x) <= 0.
 * If g(x) = g0 if indicator z = 0, and a cut is given by sum aixi + c <= aux, where xi = xi0 if z = 0 for all i,
 * then the "perspectivied" cut is sum aixi + c + (1 - z)*(g0 - c - sum aix0i) <= aux. This ensures that at z = 1,
 * the new cut is equivalent to the given cut, and at z = 0 it reduces to g0 <= aux.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimatePerspective)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* auxvar;
   int i;
   int j;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_Real cst0;
   SCIP_VAR* indicator;
   SCIP_PTRARRAY* rowpreps2;
   int nrowpreps;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "enforcement method of perspective nonlinear handler called for expr %p: ", expr);
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   assert(scip != NULL);
   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(rowpreps != NULL);
   assert(nlhdlrdata != NULL);

   auxvar = SCIPgetConsExprExprAuxVar(expr);
   assert(auxvar != NULL);

   nrowpreps = 0;
   *success = FALSE;

   SCIP_CALL( SCIPcreatePtrarray(scip, &rowpreps2) );

   /* build cuts for every indicator variable */
   for( i = 0; i < nlhdlrexprdata->nindicators; ++i )
   {
      indicator = nlhdlrexprdata->indicators[i];

      /* use cuts from every suitable nlhdlr */
      for( j = 0; j < expr->nenfos; ++j )
      {
         SCIP_Bool addedbranchscores2;
         SCIP_CONSEXPR_NLHDLR* nlhdlr2;
         int minidx;
         int maxidx;
         int r;
         int success2;

         nlhdlr2 = expr->enfos[j]->nlhdlr;

         if( !SCIPhasConsExprNlhdlrEstimate(nlhdlr2) || nlhdlr2 == nlhdlr )
            continue;

         SCIPdebugMsg(scip, "\nasking handler %s to %sestimate", SCIPgetConsExprNlhdlrName(nlhdlr2), overestimate ? "over" : "under");

         /* evaluate auxiliary before calling estimate */
         SCIP_CALL( SCIPevalauxConsExprNlhdlr(scip, nlhdlr2, expr, expr->enfos[j]->nlhdlrexprdata, &expr->enfos[j]->auxvalue, sol) );

         /* ask the handler for an estimator */
         SCIP_CALL( SCIPestimateConsExprNlhdlr(scip, conshdlr, nlhdlr2, expr, expr->enfos[j]->nlhdlrexprdata, sol, expr->enfos[j]->auxvalue,
               overestimate, SCIPgetSolVal(scip, sol, auxvar), rowpreps2, &success2, FALSE, &addedbranchscores2) );

         minidx = SCIPgetPtrarrayMinIdx(scip, rowpreps2);
         maxidx = SCIPgetPtrarrayMaxIdx(scip, rowpreps2);

         assert((success2 && minidx <= maxidx) || (!success2 && minidx > maxidx));

         for( r = minidx; r <= maxidx; ++r )
         {
            int v;
            SCIP_VAR *var;
            SCIP_Bool var_is_sc;

            rowprep = (SCIP_ROWPREP*) SCIPgetPtrarrayVal(scip, rowpreps2, r);

            assert(rowprep != NULL);

            SCIPdebugMsg(scip, "\nrowprep before perspectivy is: ");
#ifdef SCIP_DEBUG
            SCIPprintRowprep(scip, rowprep, NULL);
#endif

            /* perspectivy the estimator by adding (1-z)(g0 - c - sum aix0i),
             * where sum aixi + c = rowprep */

            /* we want cst0 = g0 - c - sum aix0i; first add g0 - c */
            cst0 = nlhdlrexprdata->exprvals0[i] + rowprep->side;

            for( v = 0; v < rowprep->nvars; ++v )
            {
               SCIP_Real val0;

               var = rowprep->vars[v];

               /* is var sc with respect to this indicator? */
               SCIP_CALL(varIsSemicontinuous(scip, var, nlhdlrdata->scvars, indicator, &val0, &var_is_sc));

               assert(var_is_sc);

               cst0 -= rowprep->coefs[v] * val0;
            }

            /* update the rowprep by adding cst0 - cst0*z */
            SCIPaddRowprepConstant(rowprep, cst0);
            SCIP_CALL(SCIPaddRowprepTerm(scip, rowprep, indicator, -cst0));

            SCIP_CALL(SCIPaddRowprepTerm(scip, rowprep, auxvar, -1.0));

            SCIPdebugMsg(scip, "\nrowprep after perspectivy is: ");
#ifdef SCIP_DEBUG
            SCIPprintRowprep(scip, rowprep, NULL);
#endif

            *success = TRUE;
            SCIP_CALL( SCIPsetPtrarrayVal(scip, rowpreps, nrowpreps, rowprep) );
            ++nrowpreps;
         }

         SCIP_CALL( SCIPclearPtrarray(scip, rowpreps2) );
      }
   }

   SCIP_CALL( SCIPfreePtrarray(scip, &rowpreps2) );

   return SCIP_OKAY;
}


/** nonlinear handler interval evaluation callback */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRINTEVAL(nlhdlrIntevalPerspective)
{ /*lint --e{715}*/
   SCIPerrorMessage("method of perspective nonlinear handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
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
#endif


/** nonlinear handler callback for reformulation */
#if 0
static
SCIP_DECL_CONSEXPR_NLHDLRREFORMULATE(nlhdlrReformulatePerspective)
{ /*lint --e{715}*/

   /* set refexpr to expr and capture it if no reformulation is possible */
   *refexpr = expr;
   SCIPcaptureConsExprExpr(*refexpr);

   return SCIP_OKAY;
}
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

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_DETECTPRIORITY,
      NLHDLR_ENFOPRIORITY, nlhdlrDetectPerspective, nlhdlrEvalauxPerspective, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/detectsum",
      "whether to run convexity detection when the root of an expression is a sum",
      &nlhdlrdata->detectsum, FALSE, DEFAULT_DETECTSUM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/multcuts",
      "whether to add cuts for all suitable indicator variables",
      &nlhdlrdata->multcuts, FALSE, DEFAULT_MULTCUTS, NULL, NULL) );

   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrPerspective);
   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrFreehdlrdataPerspective);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrFreeExprDataPerspective);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimatePerspective, NULL);

   return SCIP_OKAY;
}
