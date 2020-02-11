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
#include "scip/cons_expr_sum.h"
#include "scip/scip_sol.h"
#include "struct_cons_expr.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME               "perspective"
#define NLHDLR_DESC               "perspective handler for expressions"
#define NLHDLR_DETECTPRIORITY     -20 /* detect last so that to make use of what other handlers detected */
#define NLHDLR_ENFOPRIORITY       125 /* enforce first because perspective cuts are always stronger */

#define DEFAULT_DETECTSUM    TRUE
#define DEFAULT_MULTCUTS     TRUE

/*
 * Data structures
 */

/** data structure to store information of a semicontinuous variable */
struct SCIP_SCVarData
{
   SCIP_Real*            vals0;              /**< values of the variable when the corresponding bvars[i] = 0 */
   SCIP_VAR**            bvars;              /**< the binary variables on which the variable domain depends */
   int                   nbnds;              /**< number of suitable on/off bounds the var has */
   int                   bndssize;           /**< size of the arrays */
};
typedef struct SCIP_SCVarData SCIP_SCVARDATA;

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_EXPRCURV         curvature;          /**< curvature of the expression */

   SCIP_CONSEXPR_EXPR**  onoffterms;         /**< on/off terms for which we apply perspective cuts */
   SCIP_Real*            onoffcoefs;         /**< coefficients of onoffterms */
   SCIP_Real*            termvals0;          /**< 'off' values of the term for each indicator variable */
   SCIP_VAR***           termbvars;          /**< binary vars associated with onoffterms */
   int*                  ntermbvars;         /**< number of binary variables for each term */
   int                   nonoffterms;        /**< number of on/off expressions */
   int                   onofftermssize;     /**< size of arrays describing on/off terms */

   SCIP_CONSEXPR_EXPR**  convterms;          /**< convex terms for which we apply gradient cuts */
   SCIP_Real*            convcoefs;          /**< coefficients of convterms */
   int                   nconvterms;         /**< number of convterms */
   int                   convtermssize;      /**< size of the convterms array */

   SCIP_CONSEXPR_EXPR**  varexprs;           /**< variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */

   SCIP_VAR**            indicators;         /**< all indicator variables for the expression */
   int                   nindicators;        /**< number of indicator variables */
   SCIP_VAR***           indscvars;          /**< sorted arrays of suitable semicontinuous variables for each indicator */
   int*                  nindscvars;         /**< numbers of suitable semicontinuous variables for each indicator */
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_Bool             detectsum;          /**< whether to run detection when the root of an expression is a sum */
   SCIP_Bool             multcuts;           /**< whether to add cuts for all suitable indicator variables */
   SCIP_HASHMAP*         scvars;             /**< maps semicontinuous variables to their on/off bounds */
};

/*
 * Local methods
 */

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
   int newsize;
   SCIP_VAR** vlbvars;
   SCIP_VAR** vubvars;
   SCIP_Real* vlbcoefs;
   SCIP_Real* vubcoefs;
   SCIP_Real* vlbconstants;
   SCIP_Real* vubconstants;
   int nvlbs;
   int nvubs;
   SCIP_SCVARDATA* scvdata;
   SCIP_VAR* bvar;

   assert(scip != NULL);
   assert(var != NULL);
   assert(scvars != NULL);
   assert(result != NULL);

   scvdata = (SCIP_SCVARDATA*) SCIPhashmapGetImage(scvars, (void*)var);
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
      {/*lint --e{644}*/
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

         if( scvdata->nbnds + 1 > scvdata->bndssize )
         {
            newsize = SCIPcalcMemGrowSize(scip, scvdata->nbnds + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->bvars,  scvdata->bndssize, newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->vals0, scvdata->bndssize, newsize) );
            scvdata->bndssize = newsize;
         }
         assert(scvdata->nbnds + 1 <= scvdata->bndssize);

         scvdata->bvars[scvdata->nbnds] = bvar;
         scvdata->vals0[scvdata->nbnds] = lb0;
         ++scvdata->nbnds;
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

         if( scvdata->nbnds + 1 > scvdata->bndssize )
         {
            newsize = SCIPcalcMemGrowSize(scip, scvdata->nbnds + 1);
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->bvars, scvdata->bndssize, newsize) );
            SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &scvdata->vals0, scvdata->bndssize, newsize) );
            scvdata->bndssize = newsize;
         }
         assert(scvdata->nbnds + 1 <= scvdata->bndssize);

         scvdata->bvars[scvdata->nbnds] = bvar;
         scvdata->vals0[scvdata->nbnds] = lb0;
         ++scvdata->nbnds;
      }
   }

   if( scvdata != NULL )
   {
      /* sort bvars and vals0 */
      SCIPsortPtrReal((void**)scvdata->bvars, scvdata->vals0, SCIPvarComp, scvdata->nbnds);
      SCIPdebugMsg(scip, "var <%s> has global bounds [%f, %f] and the following on/off bounds:\n", SCIPvarGetName(var), glb, gub);
      for( c = 0; c < scvdata->nbnds; ++c )
      {
         SCIPdebugMsg(scip, " c = %d, bvar <%s>: val0 = %f\n", c, SCIPvarGetName(scvdata->bvars[c]), scvdata->vals0[c]);
      }
      SCIP_CALL( SCIPhashmapInsert(scvars, var, scvdata) );
      *result = TRUE;
   }

   return SCIP_OKAY;
}

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
      SCIP_SCVARDATA* vardata;
      vars[v] = SCIPgetConsExprExprVarVar(varexprs[v]);
      vardata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(scvars, (void*)vars[v]);

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

/* checks if an expression is semicontinuous
 *
 * An expression is semicontinuous if all of its nonlinear variables are semicontinuous
 * and share at least one common indicator variable */
static
SCIP_RETCODE exprIsSemicontinuous(
   SCIP*                         scip,            /**< SCIP data structure */
   SCIP_CONSHDLR*                conshdlr,        /**< constraint handler */
   SCIP_CONSEXPR_NLHDLRDATA*     nlhdlrdata,      /**< nonlinear handler data */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata,  /**< nlhdlr expression data */
   SCIP_CONSEXPR_EXPR*           expr,            /**< expression term to be added */
   SCIP_Bool                     isroot,          /**< is expr the root expression? */
   int*                          nbvars,          /**< buffer to store the number of indicator variables */
   SCIP_VAR**                    expr_bvars,      /**< buffer to store indicator variables of the expression */
   SCIP_Bool*                    res              /**< buffer to store whether the expression is semicontinuous */
   )
{
   int v;
   int v0;
   int nvars;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_Bool var_is_sc;
   SCIP_SCVARDATA* scvdata;
   SCIP_VAR* var;

   /* TODO this needs sc prop */
   /* we are happy with linear terms anyway */
   var = SCIPgetConsExprExprAuxVar(expr);
   if( var != NULL && !isroot )
   {
      SCIP_CALL( varIsSemicontinuous(scip, var, nlhdlrdata->scvars, NULL, NULL, &var_is_sc) );
      *res = TRUE;
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, nlhdlrexprdata->nvarexprs) );
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, varexprs, &nvars) );

   /* expression with a constant is fine */
   if( nvars == 0 )
   {
      assert(!isroot);
      *res = TRUE;
      goto TERMINATE;
   }

   *res = FALSE;

   /* all nonlinear variables of a nonlinear on/off term should be semicontinuous */
   /* TODO propagate auxiliary variables */
   for( v = 0; v < nvars; ++v )
   {
      var = SCIPgetConsExprExprVarVar(varexprs[v]);
      SCIP_CALL( varIsSemicontinuous(scip, var, nlhdlrdata->scvars, NULL, NULL, &var_is_sc) );
      if( !var_is_sc )
         goto TERMINATE;
   }

   /* find common binary variables for all variables of children[c] */
   scvdata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(nlhdlrdata->scvars, (void*)SCIPgetConsExprExprVarVar(varexprs[0]));

   if( *nbvars == 0 )
      BMScopyMemoryArray(expr_bvars, scvdata->bvars, scvdata->nbnds);

   v0 = *nbvars == 0 ? 1 : 0;
   *nbvars = scvdata->nbnds;

   SCIPdebugMsg(scip, "\nArray intersection for vars %s, *nbvars = %d", SCIPvarGetName(SCIPgetConsExprExprVarVar(varexprs[0])), *nbvars);
   for( v = v0; v < nvars; ++v )
   {
#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "\n%s; ", SCIPvarGetName(SCIPgetConsExprExprVarVar(varexprs[v])));
#endif
      scvdata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(nlhdlrdata->scvars, (void*)SCIPgetConsExprExprVarVar(varexprs[v]));
      SCIPcomputeArraysIntersectionPtr((void**)expr_bvars, *nbvars, (void**)scvdata->bvars, scvdata->nbnds, SCIPvarComp, (void**)expr_bvars, nbvars);

      /* if we have found out that the intersection is empty, expr is not semicontinuous */
      if( *nbvars == 0 )
         goto TERMINATE;
   }

   *res = TRUE;

   TERMINATE:
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
   }
   SCIPfreeBufferArray(scip, &varexprs);

   return SCIP_OKAY;
}

/** adds an expression term either to convterms or to onoffterms */
static
SCIP_RETCODE addTerm(
   SCIP*                         scip,            /**< SCIP data structure */
   SCIP_CONSHDLR*                conshdlr,        /**< constraint handler */
   SCIP_CONSEXPR_NLHDLRDATA*     nlhdlrdata,      /**< nonlinear handler data */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata,  /**< nlhdlr expression data */
   SCIP_CONSEXPR_EXPR*           term,            /**< expression term to be added */
   SCIP_Real                     coef             /**< coefficient of the term in the parent expression */
   )
{
   SCIP_CONSEXPR_EXPR** varexprs;
   int nvars, v, nbvars;
   SCIP_Bool expr_is_sc;
   SCIP_VAR** expr_bvars;

   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, nlhdlrexprdata->nvarexprs) );
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, term, varexprs, &nvars) );

   SCIP_CALL( exprIsSemicontinuous(scip, conshdlr, nlhdlrdata, nlhdlrexprdata, term, FALSE, &nbvars, expr_bvars, &expr_is_sc) );

   if( expr_is_sc )
   {
#ifdef SCIP_DEBUG
      SCIPdebugMsg(scip, "Adding on/off term: ");
      SCIPprintConsExprExpr(scip, conshdlr, term, NULL);
#endif
      /* if the term satisfies the requirements for g_i(x_i), add it to onoffterms */
      SCIP_CALL( addOnoffTerm(scip, conshdlr, nlhdlrexprdata, coef, term, expr_bvars, nbvars) );
   }
   else
      SCIP_CALL( addConvTerm(scip, nlhdlrexprdata, coef, term) );

   /* free varexprs */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
   }
   SCIPfreeBufferArray(scip, &varexprs);

   return SCIP_OKAY;
}

/** add the cut given by rowprep to sepastore */
static
SCIP_RETCODE addCut(
   SCIP*         scip,
   SCIP_CONS*    cons,
   SCIP_ROWPREP* rowprep,
   SCIP_SOL*     sol,
   double        mincutviolation,
   int*          ncuts,
   SCIP_RESULT*  result
   )
{
   SCIP_Bool success;
   SCIP_ROW* row;

   /* merge coefficients that belong to same variable */
   SCIPmergeRowprepTerms(scip, rowprep);

   SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, sol, SCIP_CONSEXPR_CUTMAXRANGE, mincutviolation, NULL, &success) );

   /* if cut looks good (numerics ok and cutting off solution), then turn into row and add to sepastore */
   if( success )
   {
      SCIP_Bool infeasible;

      SCIPdebugMsg(scip, "Separating sol point by perspective cut\n");
      SCIPdebug( SCIPprintRowprepSol(scip, rowprep, sol, NULL) );

      SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

      SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );

      if( infeasible )
      {
         *result = SCIP_CUTOFF;
      }
      else
      {
         *result = SCIP_SEPARATED;
         ++*ncuts;
      }

      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

   return SCIP_OKAY;
}

/** frees nlhdlrexprdata structure */
static
SCIP_RETCODE freeNlhdlrExprData(
   SCIP*                         scip,            /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata   /**< nlhdlr expression data */
   )
{
   int c;

   if( nlhdlrexprdata->indicators != NULL )
   {
      for( c = 0; c < nlhdlrexprdata->nindicators; ++c )
      {
         SCIPfreeBlockMemoryArray(scip, &(nlhdlrexprdata->indscvars[c]), nlhdlrexprdata->nvarexprs);
      }
      SCIPfreeBlockMemoryArray(scip, &(nlhdlrexprdata->nindscvars), nlhdlrexprdata->nindicators);
      SCIPfreeBlockMemoryArray(scip, &(nlhdlrexprdata->indscvars), nlhdlrexprdata->nindicators);
      SCIPfreeBlockMemoryArray(scip, &(nlhdlrexprdata->indicators), nlhdlrexprdata->nindicators);
   }

   if( nlhdlrexprdata->nonoffterms != 0 )
      SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->termvals0), nlhdlrexprdata->ntermbvars[0]);
   for( c = 0; c < nlhdlrexprdata->nonoffterms; ++c )
   {
      SCIPfreeBlockMemoryArray(scip, &(nlhdlrexprdata->termbvars[c]), nlhdlrexprdata->ntermbvars[c]);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->termbvars), nlhdlrexprdata->onofftermssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->ntermbvars), nlhdlrexprdata->onofftermssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->onoffcoefs), nlhdlrexprdata->onofftermssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->onoffterms), nlhdlrexprdata->onofftermssize);

   if( nlhdlrexprdata->varexprs != NULL )
   {
      for( c = 0; c < nlhdlrexprdata->nvarexprs; ++c )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(nlhdlrexprdata->varexprs[c])) );
      }
      SCIPfreeBlockMemoryArray(scip, &nlhdlrexprdata->varexprs, nlhdlrexprdata->nvarexprs);
   }

   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->convcoefs), nlhdlrexprdata->convtermssize);
   SCIPfreeBlockMemoryArrayNull(scip, &(nlhdlrexprdata->convterms), nlhdlrexprdata->convtermssize);

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
   SCIP_SCVARDATA* data;
   int c;

   if( (*nlhdlrdata)->scvars != NULL )
   {
      for( c = 0; c < SCIPhashmapGetNEntries((*nlhdlrdata)->scvars); ++c )
      {
         entry = SCIPhashmapGetEntry((*nlhdlrdata)->scvars, c);
         if( entry != NULL )
         {
            data = (SCIP_SCVARDATA*) SCIPhashmapEntryGetImage(entry);
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

/** considering a sum expression as a graph where variables are nodes and an edge exists if the two variables are
 *  contained in the same term, find the connected components of this graph
 */
SCIP_RETCODE getConnectedComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_VAR***           components,         /**< array to store the connected components */
   int*                  ncomponents,        /**< buffer to store number of connected components */
   int*                  ncompvars,          /**< array to store number of vars in each component */
   SCIP_CONSEXPR_EXPR**  varexprs,           /**< variable expressions of expr */
   int                   nvars               /**< number of variable expressions of expr */
   )
{
   int* varmarks;
   SCIP_HASHMAP* vartopos;
   int i;
   int j;
   int k;
   int nchildren;
   SCIP_CONSEXPR_EXPR** children;
   int insertpos;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(components != NULL);
   assert(ncomponents != NULL);
   assert(ncompvars != NULL);
   assert(varexprs != NULL);

   SCIP_CALL( SCIPhashmapCreate(&vartopos, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &varmarks, nvars) );

   for( i = 0; i < nvars; ++i )
   {
      int idx;

      idx = SCIPvarGetIndex(SCIPgetConsExprExprVarVar(varexprs[i]));
      SCIP_CALL( SCIPhashmapInsertInt(vartopos, (void*)(size_t) idx, i) );
      varmarks[i] = -1;
   }

   /* prepare the list of terms */
   if( SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) )
   {
      children = SCIPgetConsExprExprChildren(expr);
      nchildren = SCIPgetConsExprExprNChildren(expr);
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &children, 1) );
      *children = expr;
      nchildren = 1;
   }

   *ncomponents = 0;

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CONSEXPR_EXPR** varexprs_child;
      int nvars_child;
      int pos;
      int pos2;
      SCIP_VAR* childvar;
      int compid;

      SCIP_CALL( SCIPallocBufferArray(scip, &varexprs_child, nvars) );
      SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, children[i], varexprs_child, &nvars_child) );

      /* does any var already belong to a component? */
      compid = -1;
      for( j = 0; j < nvars_child; ++j )
      {
         childvar = SCIPgetConsExprExprVarVar(varexprs_child[j]);

         pos = SCIPhashmapGetImageInt(vartopos, (void*)(size_t) SCIPvarGetIndex(childvar));
         SCIPinfoMessage(scip, NULL, "\nvar encountered %s", SCIPvarGetName(childvar));
         if( varmarks[pos] != -1 && compid == -1 )
         {
            /* this variable already belongs to some component */
            compid = varmarks[pos];
            SCIPinfoMessage(scip, NULL, "\nset compid to %d", compid);
         }

         if( compid != -1 && varmarks[pos] != -1 && varmarks[pos] != compid )
         {
            SCIPinfoMessage(scip, NULL, "\nvar %s belongs to a different comp already! compid %d with %d vars, var compid %d with %d vars",
                  SCIPvarGetName(childvar), compid, ncompvars[compid], varmarks[pos], ncompvars[varmarks[pos]]);
            SCIP_VAR* vartomove;
            int deletedcomp;

            deletedcomp = varmarks[pos];

            /* two variables belong to different connected components -> merge the components */
            for( k = 0; k < ncompvars[deletedcomp]; ++k )
            {
               vartomove = components[deletedcomp][k];
               SCIPinfoMessage(scip, NULL, "\nset comps[%d][%d] to %s", compid, ncompvars[compid]+k, SCIPvarGetName(vartomove));
               components[compid][ncompvars[compid]+k] = vartomove;
               pos2 = SCIPhashmapGetImageInt(vartopos, (void*)(size_t) SCIPvarGetIndex(vartomove));
               varmarks[pos2] = compid;
            }

            SCIPinfoMessage(scip, NULL, "\nncompvars[%d] = %d", compid, ncompvars[compid]);
            ncompvars[compid]+= ncompvars[deletedcomp];
            ncompvars[deletedcomp] = 0;
            SCIPinfoMessage(scip, NULL, "\nncompvars[%d] = %d", compid, ncompvars[compid]);
         }
      }

      if( compid == -1 )
      {
         /* we have a new connected component */
         compid = *ncomponents;
         SCIPinfoMessage(scip, NULL, "\nnew comp, set compid to %d", compid);
         (*ncomponents)++;
         SCIP_CALL( SCIPallocBufferArray(scip, &components[compid], nvars) );
         ncompvars[compid] = 0;
      }

      /* add vars to the component */
      for( j = 0; j < nvars_child; ++j )
      {
         childvar = SCIPgetConsExprExprVarVar(varexprs_child[j]);

         pos = SCIPhashmapGetImageInt(vartopos, (void*)(size_t) SCIPvarGetIndex(childvar));
         if( varmarks[pos] == -1 )
         {
            /* this variable hasn't yet been added to the component */
            varmarks[pos] = compid;
            components[compid][ncompvars[compid]] = childvar;
            ncompvars[compid]++;
         }
         assert(varmarks[pos] == compid);
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs_child[j]) );
      }
      SCIPfreeBufferArray(scip, &varexprs_child);
   }

   /* remove empty components */
   insertpos = -1;
   for( i = 0; i < *ncomponents; ++i )
   {
      if( ncompvars[i] == 0 )
      {
         SCIPfreeBufferArray(scip, &components[i]);
         if( insertpos == -1 )
            insertpos = i;
      }
      else
      {
         SCIPsortPtr((void**)components[i], SCIPvarComp, ncompvars[i]);
         if( insertpos != -1 )
         {
            components[insertpos] = components[i];
            ncompvars[insertpos] = ncompvars[i];
            components[i] = NULL;
            ncompvars[i] = 0;
            insertpos++;
         }
      }
   }
   *ncomponents = insertpos;


   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(conshdlr) )
   {
      SCIPfreeBufferArray(scip, &children);
   }
   SCIPfreeBufferArray(scip, &varmarks);
   SCIPhashmapFree(&vartopos);

   return SCIP_OKAY;
}

/** adds an array of semicontinuous variables to the indscvar array for indicator */
SCIP_RETCODE addIndicatorSCvars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlexprdata, /**< nlhdlr expression data */
   SCIP_VAR*             indicator,          /**< indicator var */
   SCIP_VAR**            scvars,             /**< SC vars to be added to an indicator */
   int                   nscvars             /**< number of SC vars to be added to an indicator */
   )
{
   SCIP_Bool foundind;
   int pos;
   int i;
   int j;
   int k;

   SCIPinfoMessage(scip, NULL, "\nAdding scvars for indicator %s: ", SCIPvarGetName(indicator));
   for( i = 0; i < nscvars; ++i )
   {
      SCIPinfoMessage(scip, NULL, "\n%s; ", SCIPvarGetName(scvars[i]));
   }

   if( nlexprdata->indicators == NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(nlexprdata->indicators), SCIPgetNBinVars(scip)) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(nlexprdata->nindscvars), SCIPgetNBinVars(scip)) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(nlexprdata->indscvars), SCIPgetNBinVars(scip)) );
      nlexprdata->nindicators = 0;
      pos = 0;
      foundind = FALSE;
   }
   else
   {
      /* look for indicator in nlexprdata */
      foundind = SCIPsortedvecFindPtr((void**)nlexprdata->indicators, SCIPvarComp, (void*)indicator, nlexprdata->nindicators, &pos);
   }

   if( !foundind )
   {
      /* move all indicators */
      for( i = nlexprdata->nindicators; i > pos; --i )
      {
         nlexprdata->indicators[i] = nlexprdata->indicators[i-1];
         nlexprdata->indscvars[i] = nlexprdata->indscvars[i-1];
         nlexprdata->nindscvars[i] = nlexprdata->nindscvars[i-1];
      }

      /* store indicator */
      nlexprdata->indicators[pos] = indicator;
      nlexprdata->nindicators++;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(nlexprdata->indscvars[pos]), nlexprdata->nvarexprs) );

      /* just add all vars */
      BMScopyMemoryArray(nlexprdata->indscvars[pos], scvars, nscvars);
      nlexprdata->nindscvars[pos] = nscvars;
   }
   else
   {
      /* add scvars that haven't been added yet */
      j = 0;
      for( i = 0; i < nscvars; ++i )
      {
         while( SCIPvarComp(scvars[i], nlexprdata->indscvars[pos][j]) < 0 )
            ++j;

         if( SCIPvarComp(scvars[i], nlexprdata->indscvars[pos][j]) < 0 == 0 )
            continue; /* scvar i has already been added, go on to next scvar */

         /* move all indscvars */
         for( k = nlexprdata->nindscvars[pos]; k > j; --k )
         {
            nlexprdata->indscvars[pos][k] = nlexprdata->indscvars[pos][k-1];
         }

         nlexprdata->indscvars[pos][j] = scvars[i];
         nlexprdata->nindscvars[pos]++;
      }
   }

   return SCIP_OKAY;
}

/** collects the indicators and creates an array of suitable semicontinuous variables for each indicator.
 *
 *  This is done by going over every connected component and for each component, getting the indicator variables that
 *  are linked to all variables in the component. For each indicator thus chosen, all the variables of the component
 *  are added to the SC vars array.
 */
SCIP_RETCODE fillSCvarArrays(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata,     /**< nonlinear handler data */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlexprdata, /**< nonlinear expression data */
   SCIP_VAR***           comps,              /**< connected components */
   int                   ncomps,             /**< number of connected components */
   int*                  ncompvars           /**< number of variabels in each connected component */
   )
{
   int c;
   SCIP_Bool var_is_sc;
   int v;
   int i;
   SCIP_VAR* var;
   SCIP_VAR** indicators;
   SCIP_SCVARDATA* scvdata;

   SCIP_CALL( SCIPallocBufferArray(scip, &indicators, SCIPgetNBinVars(scip)) );

   /* get a list of indicators for each component */
   for( c = 0; c < ncomps; ++c )
   {
      SCIP_Bool compissc;
      int nindicators;

      compissc = TRUE;
      nindicators = 0;

      /* get an intersection of all indicators */
      for( v = 0; v < ncompvars[c]; ++v )
      {
         var = comps[c][v];
         SCIP_CALL( varIsSemicontinuous(scip, var, nlhdlrdata->scvars, NULL, NULL, &var_is_sc) );
         if( !var_is_sc )
         {
            compissc = FALSE;
            break; /* this is not an sc component for any indicator */
         }

         scvdata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(nlhdlrdata->scvars, (void*)var);

         if( nindicators == 0 )
         {
            BMScopyMemoryArray(indicators, scvdata->bvars, scvdata->nbnds);
            nindicators = scvdata->nbnds;
         }
         else
         {
            SCIPcomputeArraysIntersectionPtr((void**)indicators, nindicators, (void**)scvdata->bvars, scvdata->nbnds, SCIPvarComp, (void**)indicators, &nindicators);

            /* if the intersection is empty, component[c] is not semicontinuous */
            if( nindicators == 0 )
            {
               compissc = FALSE;
               break;
            }
         }
      }

      if( !compissc )
         break;

      /* save indicators and their sc vars */
      for( i = 0; i < nindicators; ++i )
      {
         SCIP_CALL( addIndicatorSCvars(scip, nlexprdata, indicators[i], comps[c], ncompvars[c]) );
      }
   }

   SCIPinfoMessage(scip, NULL, "\nIndicator and sc vars arrays: ");
   for( i = 0; i < nlexprdata->nindicators; ++i )
   {
      SCIPinfoMessage(scip, NULL, "\n%s (nscvars = %d): ", SCIPvarGetName(nlexprdata->indicators[i]), nlexprdata->nindscvars[i]);
      for( v = 0; v < nlexprdata->nindscvars[i]; ++v )
      {
         SCIPinfoMessage(scip, NULL, "%s; ", SCIPvarGetName(nlexprdata->indscvars[i][v]));
      }
   }

   if( nlexprdata->nindicators > 0 && nlexprdata->nindicators < SCIPgetNBinVars(scip) )
   {
      assert(nlexprdata->indicators != NULL);
      assert(nlexprdata->indscvars != NULL);
      assert(nlexprdata->nindscvars != NULL);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(nlexprdata->indicators), SCIPgetNBinVars(scip), nlexprdata->nindicators) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(nlexprdata->indscvars), SCIPgetNBinVars(scip), nlexprdata->nindicators) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(nlexprdata->nindscvars), SCIPgetNBinVars(scip), nlexprdata->nindicators) );
   }

   SCIPfreeBufferArray(scip, &indicators);

   return SCIP_OKAY;
}

/** callback to detect structure in expression tree
 *
 * We are looking for expressions of the form: \sum\limits_{i=1}^p g_i(x_i) + g_0(x_0), where:
 *  each vector x_i has a single fixed value x^{off}_i when a binary var b_i is 0;
 *  g_i, i=1,..,p are nonlinear and either all convex or all concave;
 *  g_0 is either linear or has the same curvature as g_i, i=1,..,p;
 *  p != 0.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectPerspective)
{ /*lint --e{715}*/
   int nbvars;
   int c;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_VAR** expr_bvars;
   SCIP_Real exprval0;
   SCIP_SOL* sol0;
//   SCIP_Bool isroot;
   SCIP_VAR*** components;
   SCIP_CONSEXPR_EXPR** varexprs;
   int nvars;
   int v;
   int ncomps;
   int* ncompvars;

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

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "Called perspective detect, expr = %p: ", expr);
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPdebugMsgPrint(scip, "\n");
#endif

   /* ignore sums */
   if( !nlhdlrdata->detectsum && SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) ) /*lint !e506 !e774*/
   {
      return SCIP_OKAY;
   }

   /* some other nonlinear handler should be able to separate */
   if( !(*enforcemethods & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) && !(*enforcemethods & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) )
   {
      SCIPdebugMsg(scip, "\nno enforcement method, exiting detect");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );

   (*nlhdlrexprdata)->curvature = SCIPgetConsExprExprCurvature(expr);
   SCIPdebugMsg(scip, "expr %p is %s\n", expr, (*nlhdlrexprdata)->curvature == SCIP_EXPRCURV_CONVEX ? "convex" : "concave");

   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, SCIPgetNTotalVars(scip)) );
   SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, varexprs, &(*nlhdlrexprdata)->nvarexprs) );
   nvars = (*nlhdlrexprdata)->nvarexprs;

   if( nlhdlrdata->scvars == NULL )
   {
      SCIP_CALL( SCIPhashmapCreate(&(nlhdlrdata->scvars), SCIPblkmem(scip), nvars) );
   }

   /* allocate memory for indicator variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &expr_bvars, SCIPgetNBinVars(scip)) );

   /* check that all nonlinear terms that satisfy the conditions for g_i(x_i) and the corresponding binary variables;
    * also collect information on semicontinuous variables */

   nbvars = 0;
   *success = TRUE;

   SCIP_CALL( SCIPallocBufferArray(scip, &components, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ncompvars, nvars) );

   SCIP_CALL( getConnectedComponents(scip, conshdlr, expr, components, &ncomps, ncompvars, varexprs, nvars) );
   SCIP_CALL( fillSCvarArrays(scip, nlhdlrdata, *nlhdlrexprdata, components, ncomps, ncompvars) );

   /* check if expr is semicontinuous */
//   SCIP_CALL( exprIsSemicontinuous(scip, conshdlr, nlhdlrdata, *nlhdlrexprdata, children[c], isroot, &nbvars, expr_bvars, success) );
   SCIP_CALL( exprIsSemicontinuous(scip, conshdlr, nlhdlrdata, *nlhdlrexprdata, expr, TRUE, &nbvars, expr_bvars, success) );

//   /* this needs to be checked since non-semicontinuous linear terms are allowed,
//    * so success could remain TRUE, but there could be no semicontinuous terms */
//   if( nbvars == 0 )
//      *success = FALSE;

   if( *success )
   {
      SCIP_Real* vals0;
      SCIP_VAR** vars;

      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &expr_bvars, SCIPgetNBinVars(scip), nbvars) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*nlhdlrexprdata)->termvals0), nbvars) );

      /* find the 'off' value of the expression for each indicator var */

      /* allocate memory for the solution */
      SCIP_CALL( SCIPcreateSol(scip, &sol0, NULL) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals0, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

      /* loop through indicator variables */
      for( v = 0; v < nbvars; ++v )
      {
         int v1;

         exprval0 = 0.0;

         /* get x0 */
         for( v1 = 0; v1 < nvars; ++v1 )
         {
            SCIP_Real val0;
            SCIP_Bool var_is_sc;

            vars[v1] = SCIPgetConsExprExprVarVar(varexprs[v1]);

            /* is the variable semicontinuous with respect to this indicator? */
            SCIP_CALL( varIsSemicontinuous(scip, vars[v1], nlhdlrdata->scvars, expr_bvars[v], &val0, &var_is_sc) );

            vals0[v1] = val0;
         }

         /* set x to x0 in sol0 */
         SCIP_CALL( SCIPsetSolVals(scip, sol0, nvars, vars, vals0) );

//         /* evaluate each semicontinuous term */
//         for( c = 0; c < nchildren; ++c )
//         {
//            SCIP_VAR* auxvar;
//            SCIP_Bool var_is_sc;
//
//            /* constant goes into exprval0 too */
//            if( SCIPgetConsExprExprHdlr(children[c]) == SCIPgetConsExprExprHdlrValue(conshdlr) )
//               exprval0 += SCIPgetConsExprExprValue(children[c]);
//
//            /* TODO this works best with sc prop */
//            auxvar = SCIPgetConsExprExprAuxVar(children[c]);
//
//            /* skip non-semicontinuous (with respect to this indicator) linear terms */
//            if( auxvar != NULL )
//            {
//               SCIP_CALL( varIsSemicontinuous(scip, auxvar, nlhdlrdata->scvars, expr_bvars[v], NULL, &var_is_sc) );
//               if( !var_is_sc )
//                  continue;
//            }
//
//            /* get exprval0 */
//            SCIP_CALL( SCIPevalConsExprExpr(scip, conshdlr, children[c], sol0, 0) );
//            exprval0 += SCIPgetConsExprExprValue(children[c]);
//         }

         /* save exprval0 */
         (*nlhdlrexprdata)->termvals0[v] = exprval0;
      }

      SCIPfreeBufferArray(scip, &vars);
      SCIPfreeBufferArray(scip, &vals0);
      SCIP_CALL( SCIPfreeSol(scip, &sol0) );


      SCIPdebugMsg(scip, "\ndetected an on/off expr");

      addOnoffTerm(scip, conshdlr, *nlhdlrexprdata, 1.0, expr, expr_bvars, nbvars);

      /* if we get here, enforcemethods should have already been set by other handler(s) */
      assert(((*enforcemethods & SCIP_CONSEXPR_EXPRENFO_SEPABELOW) && *enforcedbelow)
         || ((*enforcemethods & SCIP_CONSEXPR_EXPRENFO_SEPAABOVE) && *enforcedabove));

      /* save varexprs to nlhdlrexprdata */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->varexprs, (*nlhdlrexprdata)->nvarexprs) ); /* TODO this is repeating */
      SCIP_CALL( SCIPgetConsExprExprVarExprs(scip, conshdlr, expr, (*nlhdlrexprdata)->varexprs, &(*nlhdlrexprdata)->nvarexprs) );
      assert(*nlhdlrexprdata != NULL);
   }
   else
   {
      SCIPfreeBlockMemoryArray(scip, &expr_bvars, SCIPgetNBinVars(scip));
      SCIP_CALL( freeNlhdlrExprData(scip, *nlhdlrexprdata) );
      SCIPfreeBlockMemory(scip, nlhdlrexprdata);
      *nlhdlrexprdata = NULL;
   }

   for( c = 0; c < ncomps; ++c )
   {
      SCIPfreeBufferArray(scip, &components[c]);
   }
   SCIPfreeBufferArray(scip, &components);

   for( v = 0; v < nvars; ++v )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[v]) );
   }
   SCIPfreeBufferArray(scip, &varexprs);

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


/** nonlinear handler separation callback
 *
 * Applies perspective linearization to on/off terms and gradient linearization to everything else.
 */
static
SCIP_DECL_CONSEXPR_NLHDLRSEPA(nlhdlrSepaPerspective)
{ /*lint --e{715}*/
   SCIP_ROWPREP* rowprep;
   SCIP_VAR* auxvar;
   int i, j;
   SCIP_CONSEXPR_EXPR* pexpr;
   SCIP_Bool success;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_Real pcoef;
   SCIP_VAR** bvars;
   SCIP_CONSEXPR_NLHDLR* nlhdlr2;
   SCIP_Real cst0;

   *result = SCIP_DIDNOTFIND;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "sepa method of perspective nonlinear handler called for expr %p: ", expr);
   SCIP_CALL( SCIPprintConsExprExpr(scip, conshdlr, expr, NULL) );
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   assert(scip != NULL);
   assert(expr != NULL);
   assert(conshdlr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(result != NULL);
   assert(ncuts != NULL);
   assert(nlhdlrdata != NULL);

   *ncuts = 0;

   auxvar = SCIPgetConsExprExprAuxVar(expr);
   assert(auxvar != NULL);

   /* build cuts for every indicator variable */
   for( i = 0; i < nlhdlrexprdata->ntermbvars[0]; ++i )
   {
      SCIP_VAR* indicator;

      indicator = nlhdlrexprdata->termbvars[0][i];

      /* use cuts from every suitable nlhdlr */
      for( j = 0; j < expr->nenfos; ++j )
      {
         nlhdlr2 = expr->enfos[j]->nlhdlr;

         if( nlhdlr2->estimate == NULL )
            continue;

         SCIP_CALL( SCIPcreateRowprep(scip, &rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, FALSE) );

         /* ask the handler for an estimator */
         SCIP_CALL( nlhdlr2->estimate(scip, conshdlr, nlhdlr2, expr, expr->enfos[j]->nlhdlrexprdata, sol, auxvalue,
            overestimate, SCIPgetSolVal(scip, sol, auxvar), rowprep, &success) );

         if( success )
         {
            int v;
            SCIP_VAR* var;
            SCIP_Bool var_is_sc;

            SCIPinfoMessage(scip, NULL, "\nrowprep before perspectivy is: ");
            SCIPprintRowprep(scip, rowprep, NULL);

            /* perspectivy the estimator by adding (1-z)(f0 - sum aix0i)
             * (i are the indices of semicontinuous variables) */

            cst0 = nlhdlrexprdata->termvals0[i] + rowprep->side;

            SCIPinfoMessage(scip, NULL, "\nexpr = %g at %s = 0 ", cst0, SCIPvarGetName(indicator));
            SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);

            for( v = 0; v < rowprep->nvars; ++v )
            {
               SCIP_Real val0;

               var = rowprep->vars[v];

               /* is var sc with respect to this indicator? */
               SCIP_CALL( varIsSemicontinuous(scip, var, nlhdlrdata->scvars, indicator, &val0, &var_is_sc) );

               if( !var_is_sc )
                  continue;

               cst0 -= rowprep->coefs[v]*val0;
               SCIPinfoMessage(scip, NULL, "\nvar %s = %g at 0 ", SCIPvarGetName(var), val0);
            }

            /* update the rowprep by adding cst0 - cst0*z */
            SCIPaddRowprepConstant(rowprep, cst0);
            SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, indicator, -cst0) );

            SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, auxvar, -1.0) );

            SCIPinfoMessage(scip, NULL, "\nrowprep after perspectivy is: ");
            SCIPprintRowprep(scip, rowprep, NULL);

            SCIP_CALL( addCut(scip, cons, rowprep, sol, mincutviolation, ncuts, result) );
         }

         SCIPfreeRowprep(scip, &rowprep);
      }
   }

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


/** nonlinear handler callback for branching scores */
static
SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(nlhdlrBranchscorePerspective)
{ /*lint --e{715}*/
   SCIP_Real violation;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX || SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONCAVE);
   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(auxvalue == SCIPgetConsExprExprValue(expr)); /* given auxvalue should have been computed by nlhdlrEvalAuxConvex */  /*lint !e777*/
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->varexprs != NULL);
   assert(nlhdlrexprdata->nvarexprs > 0);
   assert(success != NULL);

   *success = FALSE;

   /* we separate only convex functions here, so there should be little use for branching
    * if violations are small or there are numerical issues, then we will not have generated a cut, though
    * in that case, we will still branch, that is, register branchscores for all depending var exprs
    */

   /* compute violation */
   if( auxvalue == SCIP_INVALID ) /*lint !e777*/
      violation = SCIPinfinity(scip); /* evaluation error -> we should branch */
   else if( SCIPgetConsExprExprCurvature(expr) == SCIP_EXPRCURV_CONVEX  )
      violation = auxvalue - SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));
   else
      violation = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr)) - auxvalue;

   /* if violation is not on the side that we need to enforce, then no need for branching */
   if( violation <= 0.0 )
      return SCIP_OKAY;

   /* TODO try to figure out which variables appear linear and skip them here */
   for( i = 0; i < nlhdlrexprdata->nvarexprs; ++i )
   {
      assert(nlhdlrexprdata->varexprs[i] != NULL);
      assert(SCIPisConsExprExprVar(nlhdlrexprdata->varexprs[i]));

      SCIPaddConsExprExprBranchScore(scip, nlhdlrexprdata->varexprs[i], brscoretag, violation);
   }

   *success = TRUE;

   return SCIP_OKAY;
}

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
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, nlhdlrSepaPerspective, NULL, NULL);
   SCIPsetConsExprNlhdlrBranchscore(scip, nlhdlr, nlhdlrBranchscorePerspective);

   return SCIP_OKAY;
}
