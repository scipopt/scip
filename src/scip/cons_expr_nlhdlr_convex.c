/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_nlhdlr_convex.c
 * @brief  nonlinear handler for convex expressions
 * @author Benjamin Mueller
 * @author Stefan Vigerske
 * @author Ksenia Bestuzheva
 *
 * TODO curvature information that have been computed during the detection
 *      of other nonlinear handler can not be used right now
 *
 * TODO perturb reference point if separation fails due to too large numbers
 * TODO if univariate integer, then do secant on 2 nearest integers instead of tangent
 */

#include <string.h>

#include "scip/cons_expr_nlhdlr_convex.h"
#include "scip/cons_expr.h"
#include "scip/cons_expr_iterator.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_pow.h"

/* fundamental nonlinear handler properties */
#define NLHDLR_NAME         "convex"
#define NLHDLR_DESC         "convex handler for expressions"
#define NLHDLR_PRIORITY     50

#define DEFAULT_DETECTSUM   FALSE
#define DEFAULT_PERSPECTIVE TRUE
#define DEFAULT_CVXSIGNOMIAL TRUE
#define DEFAULT_CREATEAUX   TRUE

/*
 * Data structures
 */

/** expression tree */
typedef struct NlHdlr_Expr NLHDLR_EXPR;
struct NlHdlr_Expr
{
   SCIP_CONSEXPR_EXPR*   origexpr;           /**< the original expression that is represented */
   NLHDLR_EXPR**         children;           /**< nlhdlr-expressions for children, or NULL if no children or to use auxvar */
   int                   varidx;             /**< if expression without children, then index of auxvar in nlhdlr exprdata (auxvarexprs) */

   NLHDLR_EXPR*          parent;             /**< parent in expression tree TODO: remove as unused? */

   SCIP_EXPRCURV         curv;               /**< required curvature */
   int                   depth;              /**< distance from root, root as distance 1 */

   SCIP_Real*            childrenval;        /**< buffer to store value of children */
   SCIP_Real             val;                /**< value from last eval */
   SCIP_Real             deriv;              /**< partial derivative w.r.t. root from last gradient cut computation */
};

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
   NLHDLR_EXPR*          nlexpr;             /**< subexpression for which this nlhdlr estimates */

   int                   nauxvars;           /**< number of distinct (auxiliary) variables handled */
   SCIP_CONSEXPR_EXPR**  auxvarexprs;        /**< original expressions which auxiliary variable we use */

   /* perspective cuts */
   int                   nbinvars;           /**< number of binary variables for which nlexpr is an on/off term */
   SCIP_VAR**            binvars;            /**< binary variables that make every auxvar an on/off term */
   SCIP_Real*            nlexprvalx0;        /**< value of nlexpr when binvar is 0, for each binvar */
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
   SCIP_HASHMAP*         scvars;             /***< maps semicontinuous variables to their on/off bounds */

   /* parameters */
   SCIP_Bool             detectsum;          /**< whether to run detection when the root of an expression is a sum */
   SCIP_Bool             perspective;        /**< whether to generate perspective cuts */

   /* parameters that problaby should be removed again */
   SCIP_Bool             cvxsignomial;       /**< whether to use convexity check on signomials */
   SCIP_Bool             createaux;          /**< whether to allow creating auxvars, i.e., do not require convexity in original variables */
};

/*
 * static methods
 */

/** create nlhdlr-expression
 *
 * does not create children, i.e., assumes that this will be a leaf
 */
static
SCIP_RETCODE nlhdlrExprCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR**         nlhdlrexpr,         /**< buffer to store created nlhdlr-expr */
   SCIP_CONSEXPR_EXPR*   origexpr,           /**< original expression to be represented */
   NLHDLR_EXPR*          parent,             /**< parent nlhdlr-expr, or NULL if root */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
)
{
   assert(scip != NULL);
   assert(nlhdlrexpr != NULL);
   assert(origexpr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, nlhdlrexpr) );
   (*nlhdlrexpr)->origexpr = origexpr;
   (*nlhdlrexpr)->children = NULL;
   (*nlhdlrexpr)->childrenval = NULL;
   (*nlhdlrexpr)->parent = parent;
   (*nlhdlrexpr)->curv = curv;
   (*nlhdlrexpr)->depth = parent != NULL ? parent->depth + 1 : 0;

   SCIPcaptureConsExprExpr(origexpr);

   return SCIP_OKAY;
}

/** expand nlhdlr-expression by adding children according to original expression */
static
SCIP_RETCODE nlhdlrExprGrowChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR*          nlhdlrexpr,         /**< nlhdlr-expr for which to create children */
   SCIP_EXPRCURV*        childrencurv        /**< curvature required for children, or NULL if to set to UNKNOWN */
   )
{
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(nlhdlrexpr != NULL);

   nchildren = SCIPgetConsExprExprNChildren(nlhdlrexpr->origexpr);
   if( nchildren == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexpr->children, nchildren) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexpr->childrenval, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( nlhdlrExprCreate(scip, &nlhdlrexpr->children[i], SCIPgetConsExprExprChildren(nlhdlrexpr->origexpr)[i], nlhdlrexpr,
         childrencurv != NULL ? childrencurv[i] : SCIP_EXPRCURV_UNKNOWN) );
   }

   return SCIP_OKAY;
}

/** free nlhdlr-expression, incl children */
static
void nlhdlrExprFree(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR**         nlhdlrexpr          /**< pointer to nlhdlr-expr to be freed */
)
{
   assert(nlhdlrexpr != NULL);

   if( *nlhdlrexpr == NULL )
      return;

   assert((*nlhdlrexpr)->origexpr != NULL);

   if( (*nlhdlrexpr)->children != NULL )
   {
      int nchildren = SCIPgetConsExprExprNChildren((*nlhdlrexpr)->origexpr);
      int i;

      for( i = 0; i < nchildren; ++i )
         nlhdlrExprFree(scip, &(*nlhdlrexpr)->children[i]);

      SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexpr)->childrenval, nchildren);
      SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexpr)->children, nchildren);
   }

   SCIPreleaseConsExprExpr(scip, &(*nlhdlrexpr)->origexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexpr);

   assert(*nlhdlrexpr == NULL);
}

/** evaluate nlhdlr-expression for given solution */
static
SCIP_RETCODE nlhdlrExprEval(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL */
   )
{
   int nchildren;
   int i;

   assert(nlexpr != NULL);

   if( nlexpr->children == NULL )
   {
      SCIP_VAR* var = SCIPgetConsExprExprAuxVar(nlexpr->origexpr);

      /* if do not have an auxvar, then this should be a value expression, which value is already stored in nlexpr */
      assert(var != NULL || SCIPgetConsExprExprHdlr(nlexpr->origexpr) == SCIPgetConsExprExprHdlrValue(SCIPfindConshdlr(scip, "expr")));
      assert(var != NULL || nlexpr->val == SCIPgetConsExprExprValueValue(nlexpr->origexpr));

      if( var != NULL )
         nlexpr->val = SCIPgetSolVal(scip, sol, var);
      return SCIP_OKAY;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( nlhdlrExprEval(scip, nlexpr->children[i], sol) );
      if( nlexpr->children[i]->val == SCIP_INVALID )
      {
         nlexpr->val = SCIP_INVALID;
         return SCIP_OKAY;
      }
      nlexpr->childrenval[i] = nlexpr->children[i]->val;
   }

   SCIP_CALL( SCIPevalConsExprExprHdlr(scip, nlexpr->origexpr, &nlexpr->val, nlexpr->childrenval, sol) );

   return SCIP_OKAY;
}

/** evaluate nlhdlr-expression for X0 in semi-continuous setting */
static
SCIP_RETCODE nlhdlrExprEvalX0(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   SCIP_HASHMAP*         scvars,             /**< info on semi-continuous variables */
   SCIP_VAR*             bvar                /**< binary variable to consider 0 */
   )
{
   int nchildren;
   int i;

   assert(nlexpr != NULL);

   if( nlexpr->children == NULL )
   {
      SCIP_VAR* var;
      SCIP_SCVARDATA* vardata;
      int pos;

      var = SCIPgetConsExprExprAuxVar(nlexpr->origexpr);
      /* if do not have an auxvar, then this should be a value expression, which value is already stored in nlexpr */
      assert(var != NULL || SCIPgetConsExprExprHdlr(nlexpr->origexpr) == SCIPgetConsExprExprHdlrValue(SCIPfindConshdlr(scip, "expr")));
      assert(var != NULL || nlexpr->val == SCIPgetConsExprExprValueValue(nlexpr->origexpr));
      if( var == NULL )
         return SCIP_OKAY;

      vardata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(scvars, (void*)var);
      assert(vardata != NULL);

      /* find bvar in vardata->bvars */
      (void) SCIPsortedvecFindPtr((void**)vardata->bvars, SCIPvarComp, (void*)bvar, vardata->nbnds, &pos);
      assert(pos < vardata->nbnds);
      assert(vardata->bvars[pos] == bvar);

      nlexpr->val = vardata->vals0[pos];

      return SCIP_OKAY;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( nlhdlrExprEvalX0(scip, nlexpr->children[i], scvars, bvar) );
      if( nlexpr->children[i]->val == SCIP_INVALID )
      {
         nlexpr->val = SCIP_INVALID;
         return SCIP_OKAY;
      }
      nlexpr->childrenval[i] = nlexpr->children[i]->val;
   }

   SCIP_CALL( SCIPevalConsExprExprHdlr(scip, nlexpr->origexpr, &nlexpr->val, nlexpr->childrenval, NULL) );

   return SCIP_OKAY;
}

/* differentiate nlhdlr-expression and store contribution to gradient or perspective cut in rowprep */
static
SCIP_RETCODE nlhdlrExprGradientCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   SCIP_VAR*             perspvar,           /**< binary variable to use for perspective cut, or NULL if ordinary gradient cut */
   SCIP_HASHMAP*         scvars,             /**< info on semi-continuous variables */
   SCIP_ROWPREP*         rowprep,            /**< row where to add coefficients and constants */
   SCIP_Bool*            success             /**< set to FALSE if gradient evaluation error */
   )
{
   int nchildren;
   int i;

   assert(nlexpr != NULL);
   assert(success != NULL);

   *success = TRUE;

   if( nlexpr->children == NULL )
   {
      SCIP_VAR* var;

      var = SCIPgetConsExprExprAuxVar(nlexpr->origexpr);
      /* if do not have an auxvar, then this should be a value expression, which value is already stored in nlexpr */
      assert(var != NULL || SCIPgetConsExprExprHdlr(nlexpr->origexpr) == SCIPgetConsExprExprHdlrValue(SCIPfindConshdlr(scip, "expr")));
      assert(var != NULL || nlexpr->val == SCIPgetConsExprExprValueValue(nlexpr->origexpr));
      if( var == NULL )
         return SCIP_OKAY;

      if( perspvar == NULL )
      {
         SCIPdebugMsg(scip, "add %g * (<%s> - %g) to rowprep\n", nlexpr->deriv, SCIPvarGetName(var), nlexpr->val);

         /* add deriv * (var - varval) to rowprep */
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, nlexpr->deriv) );
         SCIPaddRowprepConstant(rowprep, -nlexpr->deriv * nlexpr->val);
      }
      else
      {
         /* add deriv (x - x0) - deriv (sol - x0) perspvar to rowprep */
         SCIP_SCVARDATA* vardata;
         SCIP_Real varval0;
         int pos;

         vardata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(scvars, (void*)var);
         assert(vardata != NULL);

         /* find perspvar in vardata->bvars */
         (void) SCIPsortedvecFindPtr((void**)vardata->bvars, SCIPvarComp, (void*)perspvar, vardata->nbnds, &pos);
         assert(pos < vardata->nbnds);
         assert(vardata->bvars[pos] == perspvar);
         varval0 = vardata->vals0[pos];

         SCIPdebugMsg(scip, "add %g (<%s> %+g) %+g (%g %+g) <%s> to rowprep\n",
            nlexpr->deriv, SCIPvarGetName(var), -varval0, -nlexpr->deriv, nlexpr->val, -varval0, SCIPvarGetName(perspvar));

         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, nlexpr->deriv) );
         SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, perspvar, -nlexpr->deriv * (nlexpr->val - varval0)) );
         SCIPaddRowprepConstant(rowprep, -nlexpr->deriv * varval0);
      }

      return SCIP_OKAY;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);

   /* assemble children values
    * should still be uptodate from last eval
   for( i = 0; i < nchildren; ++i )
      nlexpr->childrenval[i] = nlexpr->children[i]->val;
    */

   /* differentiate w.r.t. each child */
   for( i = 0; i < nchildren; ++i )
   {
      /* compute partial derivative w.r.t. child i */
      SCIP_CALL( SCIPbwdiffConsExprExprHdlr(scip, nlexpr->origexpr, i, &nlexpr->children[i]->deriv, nlexpr->childrenval, nlexpr->val) );
      if( nlexpr->children[i]->deriv == SCIP_INVALID )
      {
         *success = FALSE;
         break;
      }
      nlexpr->children[i]->deriv *= nlexpr->deriv;

      /* push derivatives further down */
      SCIP_CALL( nlhdlrExprGradientCut(scip, conshdlr, nlexpr->children[i], perspvar, scvars, rowprep, success) );
      if( !*success )
         break;
   }

   return SCIP_OKAY;
}

/** register violation as branchscore in leafs of nlhdlr-expression */
static
void nlhdlrExprBranchscore(
   SCIP*                 scip,               /**< SCIP data structure */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expression */
   unsigned int          brscoretag,         /**< branchscore tag */
   SCIP_Real             violation           /**< violation */
   )
{
   int nchildren;
   int c;

   if( nlexpr->children == NULL )
   {
      SCIPaddConsExprExprBranchScore(scip, nlexpr->origexpr, brscoretag, violation);
      return;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);
   for( c = 0; c < nchildren; ++c )
      nlhdlrExprBranchscore(scip, nlexpr->children[c], brscoretag, violation);
}

/** collect the original expressions which are leafs in the nlexpr, thus which auxvar we use, and assign an index to them
 *
 * Also ensure that for every leaf there is an auxiliary variable.
 */
static
SCIP_RETCODE collectAuxVarExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expr */
   SCIP_HASHMAP*         origexpr2index,     /**< mapping from original expression to index */
   int*                  nindices            /**< number of indices */
   )
{
   int nchildren;
   int i;

   assert(nlexpr != NULL);
   assert(origexpr2index != NULL);
   assert(nindices != NULL);

   if( nlexpr->children == NULL )
   {
      /* do not add auxvars for value expressions, we can handle them as such */
      if( SCIPgetConsExprExprHdlr(nlexpr->origexpr) == SCIPgetConsExprExprHdlrValue(conshdlr) )
      {
         assert(SCIPgetConsExprExprAuxVar(nlexpr->origexpr) == NULL);  /* we use no-auxvar as an indication for value-expressions in eval, etc */
         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, nlexpr->origexpr, NULL) );

      if( SCIPhashmapExists(origexpr2index, (void*)nlexpr->origexpr) )
      {
         nlexpr->varidx = SCIPhashmapGetImageInt(origexpr2index, (void*)nlexpr->origexpr);
      }
      else
      {
         nlexpr->varidx = (*nindices)++;
         SCIP_CALL( SCIPhashmapInsertInt(origexpr2index, (void*)nlexpr->origexpr, nlexpr->varidx) );
      }

      return SCIP_OKAY;
   }

   nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);
   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( collectAuxVarExprs(scip, conshdlr, nlexpr->children[i], origexpr2index, nindices) );
   }

   return SCIP_OKAY;
}

/** looks whether top of given expression looks like a monomial that can have a given curvature
 * e.g., sqrt(x)*sqrt(y) is convex if x,y >= 0 and x and y are convex
 * unfortunately, doesn't work for tls, because i) it's originally sqrt(x*y), and ii) it expanded into some sqrt(z*y+y)
 * but works for cvxnonsep_nsig
 */
static
SCIP_RETCODE constructExprCheckMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   NLHDLR_EXPR*          nlexpr,             /**< nlhdlr-expr to check */
   NLHDLR_EXPR***        stack,              /**< pointer to stack where to add generated leafs */
   int*                  stackpos,           /**< current top position of stack */
   int*                  stacksize,          /**< length of *stack */
   SCIP_Bool*            success             /**< whether we found something */
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_Real* exponents;
   SCIP_INTERVAL* bounds;
   SCIP_EXPRCURV* curv;
   int nfactors;
   int i;

   assert(nlexpr != NULL);
   assert(stack != NULL);
   assert(*stack != NULL);
   assert(stackpos != NULL);
   assert(*stackpos >= -1);
   assert(stacksize != NULL);
   assert(success != NULL);

   *success = FALSE;

   expr = nlexpr->origexpr;
   assert(expr != NULL);

   if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrProduct(conshdlr) )
      return SCIP_OKAY;

   nfactors = SCIPgetConsExprExprNChildren(expr);
   if( nfactors <= 1 )  /* boooring */
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &exponents, nfactors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nfactors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &curv, nfactors) );

   for( i = 0; i < nfactors; ++i )
   {
      child = SCIPgetConsExprExprChildren(expr)[i];
      assert(child != NULL);

      if( SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         exponents[i] = 1.0;
         bounds[i] = SCIPgetConsExprExprActivity(scip, child);
      }
      else
      {
         exponents[i] = SCIPgetConsExprExprPowExponent(child);
         bounds[i] = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(child)[0]);
      }
   }

   if( !SCIPexprcurvMonomialInv(SCIPexprcurvMultiply(SCIPgetConsExprExprProductCoef(expr), nlexpr->curv), nfactors, exponents, bounds, curv) )
      goto TERMINATE;

   /* add immediate children to nlexpr
    * some entries in curv actually apply to arguments of pow's, will correct this next
    */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr, curv) );
   assert(nlexpr->children != NULL);

   /* make sure there is enough space on the stack for all children (or their children)  */
   if( *stackpos+1 + nfactors > *stacksize )
   {
      *stacksize = SCIPcalcMemGrowSize(scip, *stackpos+1 + nfactors);
      SCIP_CALL( SCIPreallocBufferArray(scip, stack, *stacksize) );
   }

   /* put children that are not power on stack
    * grow child for children that are power and put this child on stack
    */
   for( i = 0; i < nfactors; ++i )
   {
      child = SCIPgetConsExprExprChildren(expr)[i];
      assert(child != NULL);
      assert(nlexpr->children[i] != NULL);

      if( SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         assert(*stackpos+1 < *stacksize);
         (*stack)[++*stackpos] = nlexpr->children[i];
      }
      else
      {
         SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr->children[i], &curv[i]) );
         assert(nlexpr->children[i]->children != NULL);
         (*stack)[++*stackpos] = nlexpr->children[i]->children[0];
      }
   }

   *success = TRUE;

TERMINATE:
   SCIPfreeBufferArray(scip, &curv);
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &exponents);

   return SCIP_OKAY;
}

/** construct a subexpression (as nlhdlr-expression) of maximal size that has a given curvature
 *
 * If the curvature cannot be achieved for an expression in the original expression graph,
 * then this expression becomes a leaf in the nlhdlr-expression.
 *
 * Sets *rootnlexpr to NULL if failed.
 */
static
SCIP_RETCODE constructExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata,     /**< nonlinear handler data */
   NLHDLR_EXPR**         rootnlexpr,         /**< buffer to store created expression */
   int*                  maxdepth,           /**< buffer to store created depth of expression */
   SCIP_CONSEXPR_EXPR*   rootexpr,           /**< expression */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
   )
{
   NLHDLR_EXPR* nlexpr;
   SCIP_EXPRCURV* childcurv;
   int childcurvsize;
   int nchildren;
   SCIP_Bool success;

   /* to do list: expressions where to check whether they can have the desired curvature when taking their children into account */
   NLHDLR_EXPR** stack;
   int stacksize;
   int stackpos;

   assert(scip != NULL);
   assert(rootnlexpr != NULL);
   assert(maxdepth != NULL);
   assert(rootexpr != NULL);
   assert(curv == SCIP_EXPRCURV_CONVEX || curv == SCIP_EXPRCURV_CONCAVE);

   childcurvsize = SCIPgetConsExprExprNChildren(rootexpr);
   SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, childcurvsize) );

   /* create root nlhdlr-expression */
   SCIP_CALL( nlhdlrExprCreate(scip, rootnlexpr, rootexpr, NULL, curv) );
   *maxdepth = 0;

   stacksize = 20;
   SCIP_CALL( SCIPallocBufferArray(scip, &stack, stacksize) );

   stack[0] = *rootnlexpr;
   stackpos = 0;
   while( stackpos >= 0 )
   {
      /* take expression from stack */
      nlexpr = stack[stackpos--];
      assert(nlexpr != NULL);
      assert(nlexpr->origexpr != NULL);

      *maxdepth = MAX(*maxdepth, nlexpr->depth);
      nchildren = SCIPgetConsExprExprNChildren(nlexpr->origexpr);

      /* TODO if sum with > N children, then always turn into variable?
       * especially if we have a big linear term, this would help to consider less vars for vertex-polyhedral funcs
       */

      /* TODO if bwdiff not implemented, then proceed further? */

      if( childcurvsize < nchildren )
      {
         childcurvsize = SCIPcalcMemGrowSize(scip, nchildren);
         SCIP_CALL( SCIPreallocBufferArray(scip, &childcurv, childcurvsize) );
      }

      success = TRUE;
      if( nlexpr->curv != SCIP_EXPRCURV_UNKNOWN )
      {
         /* TODO here should be a nice system where a number of convexity-detection rules can be tried */
         if( nlhdlrdata->cvxsignomial )
         {
            SCIP_CALL( constructExprCheckMonomial(scip, conshdlr, nlexpr, &stack, &stackpos, &stacksize, &success) );
            if( success )
               continue;  /* constructExprCheckMonomial will have updated stack */
         }

         /* check whether and under which conditions nlexpr->origexpr can have desired curvature */
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, nlexpr->origexpr, nlexpr->curv, &success, childcurv) );
         /* SCIPprintConsExprExpr(scip, conshdlr, nlexpr->origexpr, NULL);
         SCIPinfoMessage(scip, NULL, " is %s? %d\n", SCIPexprcurvGetName(nlexpr->curv), success); */
         if( success )
         {
            /* if origexpr can have curvature curv, then don't treat it as leaf, but include its children */
            SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr, childcurv) );
         }
         else if( nlexpr == *rootnlexpr )
         {
            /* if there is no way to ensure curv for root expression, then we failed */
            nlhdlrExprFree(scip, rootnlexpr);
            break;
         }
         else if( !nlhdlrdata->createaux )
         {
            /* require that desired curvature can be achieved w.r.t. original variables, i.e., without introducing auxvars */
            nlhdlrExprFree(scip, rootnlexpr);
            break;
         }
      }
      else if( nchildren > 0 )
      {
         /* if we don't care about curvature in this subtree anymore (very unlikely),
          * then only continue iterating this subtree to assemble leaf expressions
          */
         SCIP_CALL( nlhdlrExprGrowChildren(scip, nlexpr, NULL) );
      }

      if( success && nchildren > 0 )
      {
         /* add children expressions to to do list (stack) */
         assert(nlexpr->children != NULL);

         if( stackpos+1 + nchildren > stacksize )
         {
            stacksize = SCIPcalcMemGrowSize(scip, stackpos+1 + nchildren);
            SCIP_CALL( SCIPreallocBufferArray(scip, &stack, stacksize) );
         }
         memcpy(stack + (stackpos+1), nlexpr->children, nchildren * sizeof(NLHDLR_EXPR*));
         stackpos += nchildren;
      }
   }

   SCIPfreeBufferArray(scip, &stack);
   SCIPfreeBufferArray(scip, &childcurv);

   return SCIP_OKAY;
}

/** creates nonlinear handler expression data structure */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nlhdlr expression data */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   NLHDLR_EXPR*          nlexpr              /**< nlhdlr expression */
   )
{
   SCIP_HASHMAP* origexpr2index;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata == NULL);
   assert(nlexpr != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   (*nlhdlrexprdata)->nlexpr = nlexpr;

#ifdef SCIP_DEBUG
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, " is handled as %s\n", SCIPexprcurvGetName(nlexpr->curv));
#endif

   /* make sure there are auxvars and collect them */
   SCIP_CALL( SCIPhashmapCreate(&origexpr2index, SCIPblkmem(scip), 10) );   /* TODO replace 10 by number of leafs, which we could count in constructExpr */
   (*nlhdlrexprdata)->nauxvars = 0;
   SCIP_CALL( collectAuxVarExprs(scip, conshdlr, nlexpr, origexpr2index, &(*nlhdlrexprdata)->nauxvars) );

   /* assemble auxvarexprs array */
   assert((*nlhdlrexprdata)->nauxvars > 0);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->auxvarexprs, (*nlhdlrexprdata)->nauxvars) );
   for( i = 0; i < SCIPhashmapGetNEntries(origexpr2index); ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_CONSEXPR_EXPR* origexpr;
      int idx;

      entry = SCIPhashmapGetEntry(origexpr2index, i);
      if( entry == NULL )
         continue;

      origexpr = (SCIP_CONSEXPR_EXPR*) SCIPhashmapEntryGetOrigin(entry);
      assert(expr != NULL);

      idx = SCIPhashmapEntryGetImageInt(entry);
      assert(idx >= 0);
      assert(idx < (*nlhdlrexprdata)->nauxvars);

      (*nlhdlrexprdata)->auxvarexprs[idx] = origexpr;

      SCIPdebugMsg(scip, "auxvar %d: <%s>\n", idx, SCIPvarGetName(SCIPgetConsExprExprAuxVar(origexpr)));
   }

   SCIPhashmapFree(&origexpr2index);

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
   SCIP_Bool*            result              /**< buffer to store whether var is semicontinuous */
   )
{
   SCIP_Real lb0, ub0, lb1, ub1, glb, gub;
   SCIP_Bool exists;
   int c, pos, newsize;
   SCIP_VAR** vlbvars;
   SCIP_VAR** vubvars;
   SCIP_Real* vlbcoefs;
   SCIP_Real* vubcoefs;
   SCIP_Real* vlbconstants;
   SCIP_Real* vubconstants;
   int nvlbs, nvubs;
   SCIP_SCVARDATA* scvdata;
   SCIP_VAR* bvar;

   assert(scip != NULL);
   assert(var != NULL);
   assert(scvars != NULL);
   assert(result != NULL);

   *result = FALSE;

   scvdata = (SCIP_SCVARDATA*) SCIPhashmapGetImage(scvars, (void*)var);
   if( scvdata != NULL )
   {
      *result = TRUE;
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
         SCIPdebugMsgPrint(scip, ", upper bound: %f <%s> %+f", vubcoefs[pos], SCIPvarGetName(vubvars[pos]), vubconstants[pos]);

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

/** check whether (aux)variables are semicontinuous w.r.t. one or several common binary variables */
static
SCIP_RETCODE checkSemicontinuous(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_NLHDLREXPRDATA* nlhdlrexprdata, /**< nlhdlr expression data */
   SCIP_HASHMAP*         scvars              /**< maps semicontinuous variables to their on/off bounds */
   )
{
   SCIP_VAR* var;
   SCIP_SCVARDATA* scvdata;
   SCIP_VAR** expr_bvars = NULL;
   SCIP_Bool var_is_sc;
   int nbvars;
   int nbvars0;
   int nvars;
   int v;

   assert(scip != NULL);
   assert(scvars != NULL);
   assert(nlhdlrexprdata != NULL);

   nvars = nlhdlrexprdata->nauxvars;
   assert(nvars > 0);

   /* check whether all variables are semicontinuous, so nlexpr can be handled as an on/off term */
   for( v = 0; v < nvars; ++v )
   {
      var = SCIPgetConsExprExprAuxVar(nlhdlrexprdata->auxvarexprs[v]);
      SCIP_CALL( varIsSemicontinuous(scip, var, scvars, &var_is_sc) );
      if( !var_is_sc )
         return SCIP_OKAY;
   }

   /* find common binary variables for all variables */
   scvdata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(scvars, (void*)SCIPgetConsExprExprAuxVar(nlhdlrexprdata->auxvarexprs[0]));
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &expr_bvars, scvdata->bvars, scvdata->nbnds) );
   nbvars0 = nbvars = scvdata->nbnds;

   SCIPdebugMsg(scip, "Array intersection for vars %s", SCIPvarGetName(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->auxvarexprs[0])));
   for( v = 1; v < nvars; ++v )
   {
      SCIPdebugMsgPrint(scip, ", %s", SCIPvarGetName(SCIPgetConsExprExprAuxVar(nlhdlrexprdata->auxvarexprs[v])));
      scvdata = (SCIP_SCVARDATA*)SCIPhashmapGetImage(scvars, (void*)SCIPgetConsExprExprAuxVar(nlhdlrexprdata->auxvarexprs[v]));
      SCIPcomputeArraysIntersectionPtr((void**)expr_bvars, nbvars, (void**)scvdata->bvars, scvdata->nbnds, SCIPvarComp, (void**)expr_bvars, &nbvars);

      /* if we have found out that the intersection is empty, we have no on/off term */
      if( nbvars == 0 )
      {
         SCIPfreeBlockMemoryArrayNull(scip, &expr_bvars, nbvars0);
         return SCIP_OKAY;
      }
   }
   SCIPdebugMsgPrint(scip, " is: ");
   for( v = 0; v < nbvars; ++v )
   {
      SCIPdebugMsgPrint(scip, "%s; ", SCIPvarGetName(expr_bvars[v]));
   }
   SCIPdebugMsgPrint(scip, "\n");

   /* store binary variables in nlhdlr's expr data */
   SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &expr_bvars, nbvars0, nbvars) );
   nlhdlrexprdata->nbinvars = nbvars;
   nlhdlrexprdata->binvars = expr_bvars;

   /* evaluate expression in the point defined by the binary variable being 0, for each binary variable */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nlhdlrexprdata->nlexprvalx0, nbvars) );
   for( v = 0; v < nbvars; ++v )
   {
      SCIP_CALL( nlhdlrExprEvalX0(scip, nlhdlrexprdata->nlexpr, scvars, expr_bvars[v]) );
      nlhdlrexprdata->nlexprvalx0[v] = nlhdlrexprdata->nlexpr->val;
   }

   return SCIP_OKAY;
}

/*
 * Callback methods of nonlinear handler
 */

static
SCIP_DECL_CONSEXPR_NLHDLRFREEHDLRDATA(nlhdlrfreeHdlrDataConvex)
{
   assert(scip != NULL);
   assert(nlhdlrdata != NULL);
   assert(*nlhdlrdata != NULL);

   if( (*nlhdlrdata)->scvars != NULL )
   {
      SCIP_SCVARDATA* scvdata;
      int c;

      for( c = 0; c < SCIPhashmapGetNEntries((*nlhdlrdata)->scvars); ++c )
      {
         SCIP_HASHMAPENTRY* entry = SCIPhashmapGetEntry((*nlhdlrdata)->scvars, c);
         if( entry != NULL )
         {
            scvdata = (SCIP_SCVARDATA*) SCIPhashmapEntryGetImage(entry);
            SCIPfreeBlockMemoryArray(scip, &scvdata->vals0, scvdata->bndssize);
            SCIPfreeBlockMemoryArray(scip, &scvdata->bvars, scvdata->bndssize);
            SCIPfreeBlockMemory(scip, &scvdata);
         }
      }
      SCIPhashmapFree(&(*nlhdlrdata)->scvars);
   }

   SCIPfreeBlockMemory(scip, nlhdlrdata);

   return SCIP_OKAY;
}

/** callback to free expression specific data */
static
SCIP_DECL_CONSEXPR_NLHDLRFREEEXPRDATA(nlhdlrfreeExprDataConvex)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexprdata)->nlexprvalx0, (*nlhdlrexprdata)->nbinvars);
   SCIPfreeBlockMemoryArrayNull(scip, &(*nlhdlrexprdata)->binvars, (*nlhdlrexprdata)->nbinvars);
   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->auxvarexprs, (*nlhdlrexprdata)->nauxvars);
   nlhdlrExprFree(scip, &(*nlhdlrexprdata)->nlexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** the detection assumes that the curvature information of the expression has been computed already */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectConvex)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   NLHDLR_EXPR* nlexpr = NULL;
   int depth;

   assert(scip != NULL);
   assert(nlhdlr != NULL);
   assert(expr != NULL);
   assert(enforcemethods != NULL);
   assert(enforcedbelow != NULL);
   assert(enforcedabove != NULL);
   assert(success != NULL);
   assert(nlhdlrexprdata != NULL);

   *success = FALSE;

   /* we currently cannot contribute in presolve */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* ignore sums */
   if( !nlhdlrdata->detectsum && SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) ) /*lint !e506 !e774*/
      return SCIP_OKAY;

   /* ignore pure constants and variables */
   if( SCIPgetConsExprExprNChildren(expr) == 0 )
      return SCIP_OKAY;

   /* TODO we are also interested in handling the concave side? (-> nlhdlr_vertexpolyhedral?) */

   if( !*enforcedbelow )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, nlhdlrdata, &nlexpr, &depth, expr, SCIP_EXPRCURV_CONVEX) );
      assert(nlexpr == NULL || depth >= 1);
      if( depth < 1 )  /* TODO change back to <= 1, i.e. free if only immediate children, but no grand-daughters; but what if we can do perspective? */
      {
         nlhdlrExprFree(scip, &nlexpr);
      }
      else if( nlexpr != NULL )
      {
         *enforcedbelow = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be convex -> can enforce expr <= auxvar\n", (void*)expr);
      }
   }

   if( !*enforcedabove && nlexpr == NULL )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, nlhdlrdata, &nlexpr, &depth, expr, SCIP_EXPRCURV_CONCAVE) );
      assert(nlexpr == NULL || depth >= 1);
      if( depth < 1 )  /* TODO change back to <= 1, see above */
      {
         nlhdlrExprFree(scip, &nlexpr);
      }
      else if( nlexpr != NULL )
      {
         *enforcedabove = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be concave -> can enforce expr >= auxvar\n", (void*)expr);
      }
   }

   assert(*success || nlexpr == NULL);
   if( !*success )
      return SCIP_OKAY;


   /* store variable expressions into the expression data of the nonlinear handler */
   SCIP_CALL( createNlhdlrExprData(scip, conshdlr, nlhdlrexprdata, expr, nlexpr) );

   if( nlhdlrdata->perspective )
   {
      if( nlhdlrdata->scvars == NULL )
      {
         SCIP_CALL( SCIPhashmapCreate(&nlhdlrdata->scvars, SCIPblkmem(scip), (*nlhdlrexprdata)->nauxvars) );
      }

      SCIP_CALL( checkSemicontinuous(scip, *nlhdlrexprdata, nlhdlrdata->scvars) );
   }

   return SCIP_OKAY;
}

/** auxiliary evaluation callback */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalAuxConvex)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlexpr != NULL);
   assert(nlhdlrexprdata->nlexpr->origexpr == expr);
   assert(auxvalue != NULL);

   SCIP_CALL( nlhdlrExprEval(scip, nlhdlrexprdata->nlexpr, sol) );
   *auxvalue = nlhdlrexprdata->nlexpr->val;

   return SCIP_OKAY;
}

/** estimator callback */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateConvex)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   NLHDLR_EXPR* nlexpr;
   SCIP_VAR* perspvar;
   SCIP_Real perspvarval;
   SCIP_Real nlexprvalx0;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);
   assert(nlexpr->origexpr == expr);
   assert(nlexpr->curv == SCIP_EXPRCURV_CONVEX || nlexpr->curv == SCIP_EXPRCURV_CONCAVE);
   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* if estimating on non-convex side, then do nothing
    * TODO we are vertex-polyhedral and so should compute something
    */
   if( ( overestimate && nlexpr->curv == SCIP_EXPRCURV_CONVEX) ||
       (!overestimate && nlexpr->curv == SCIP_EXPRCURV_CONCAVE) )
      return SCIP_OKAY;

   /* we can skip eval as nlhdlrEvalAux should have been called for same solution before */
   /* SCIP_CALL( nlhdlrExprEval(scip, nlexpr, sol) ); */
   assert(auxvalue == nlexpr->val); /* given value (originally from nlhdlrEvalAuxConvex) should coincide with the one stored in nlexpr */  /*lint !e777*/
   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(nlexpr->val)) )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", nlexpr->val, (void*)expr);
      return SCIP_OKAY;
   }

   nlhdlrdata = SCIPgetConsExprNlhdlrData(nlhdlr);
   assert(nlhdlrdata != NULL);

   /* heuristically choose the most promising binary variable (one closest to 0) for perspective formulation */  /*TODO? nlhdlr_perspective is generating one cut for each binary var by default */
   perspvar = NULL;
   perspvarval = 2.0;
   for( i = 0; i < nlhdlrexprdata->nbinvars; ++i )
   {
      SCIP_Real bval = SCIPgetSolVal(scip, sol, nlhdlrexprdata->binvars[i]);
      if( bval < perspvarval && nlhdlrexprdata->nlexprvalx0[i] != SCIP_INVALID )
      {
         perspvarval = bval;
         perspvar = nlhdlrexprdata->binvars[i];
         nlexprvalx0 = nlhdlrexprdata->nlexprvalx0[i];
      }
   }

   /* compute gradient or perspective cut, add to rowprep
    *
    * gradient cut adds (x - sol) \nabla f(sol) + f(sol)
    * perspective cut adds (x - x0) \nabla f(sol) + (f(sol) - f(x0) - (sol - x0) \nabla f(sol)) perspvar + f(x0)
    */
   nlexpr->deriv = 1.0;
   *success = TRUE;
   SCIP_CALL( nlhdlrExprGradientCut(scip, conshdlr, nlexpr, perspvar, nlhdlrdata->scvars, rowprep, success) );

   if( !*success )
      return SCIP_OKAY;

   if( perspvar == NULL )
   {
      /* add f(sol) */
      SCIPaddRowprepConstant(rowprep, nlexpr->val);
   }
   else
   {
      /* add (f(sol)-f(x0)) perspvar + f(x0) */
      SCIPdebugMsg(scip, "add (%g %+g) <%s> %+g to rowprep\n", nlexpr->val, -nlexprvalx0, SCIPvarGetName(perspvar), nlexprvalx0);
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, perspvar, nlexpr->val - nlexprvalx0) );
      SCIPaddRowprepConstant(rowprep, nlexprvalx0);
   }
   rowprep->local = FALSE;

   (void) SCIPsnprintf(rowprep->name, SCIP_MAXSTRLEN, "%sestimate_convex%p_%s%d",
      overestimate ? "over" : "under",
      (void*)expr,
      sol != NULL ? "sol" : "lp",
      sol != NULL ? SCIPsolGetIndex(sol) : SCIPgetNLPs(scip));

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRBRANCHSCORE(nlhdlrBranchscoreConvex)
{ /*lint --e{715}*/
   NLHDLR_EXPR* nlexpr;
   SCIP_Real violation;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(success != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);

   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(auxvalue == nlexpr->val); /* given auxvalue should have been computed by nlhdlrEvalAuxConvex */  /*lint !e777*/

   *success = FALSE;

   /* we separate only convex functions here, so there should be little use for branching
    * if violations are small or there are numerical issues, then we will not have generated a cut, though
    * in that case, we will still branch, that is, register branchscores for all depending var exprs
    */

   /* compute violation */
   if( auxvalue == SCIP_INVALID ) /*lint !e777*/
      violation = SCIPinfinity(scip); /* evaluation error -> we should branch */
   else if( nlexpr->curv == SCIP_EXPRCURV_CONVEX  )
      violation = auxvalue - SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));
   else
      violation = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr)) - auxvalue;

   /* if violation is not on the side that we need to enforce, then no need for branching */
   if( violation <= 0.0 )
      return SCIP_OKAY;

   /* TODO try to figure out which variables appear linear and skip them here */
   nlhdlrExprBranchscore(scip, nlhdlrexprdata->nlexpr, brscoretag, violation);

   *success = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_NLHDLRCOPYHDLR(nlhdlrCopyhdlrConvex)
{ /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetconsexprhdlr != NULL);
   assert(sourcenlhdlr != NULL);
   assert(strcmp(SCIPgetConsExprNlhdlrName(sourcenlhdlr), NLHDLR_NAME) == 0);

   SCIP_CALL( SCIPincludeConsExprNlhdlrConvex(targetscip, targetconsexprhdlr) );

   return SCIP_OKAY;
}

/** includes convex nonlinear handler to consexpr */
SCIP_RETCODE SCIPincludeConsExprNlhdlrConvex(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_NLHDLR* nlhdlr;
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;

   assert(scip != NULL);
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &nlhdlrdata) );
   nlhdlrdata->scvars = NULL;

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectConvex, nlhdlrEvalAuxConvex, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/detectsum",
      "whether to run convexity detection when the root of an expression is a sum",
      &nlhdlrdata->detectsum, FALSE, DEFAULT_DETECTSUM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/perspective",
      "whether to check for semicontinuous variables and use perspective cuts",
      &nlhdlrdata->perspective, FALSE, DEFAULT_PERSPECTIVE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/cvxsignomial",
      "whether to use convexity check on signomials",
      &nlhdlrdata->cvxsignomial, FALSE, DEFAULT_CVXSIGNOMIAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/createaux",
      "whether to allow creating auxvars, i.e., do not require convexity in original variables",
      &nlhdlrdata->createaux, FALSE, DEFAULT_CREATEAUX, NULL, NULL) );

   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrfreeHdlrDataConvex);
   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrConvex);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataConvex);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateConvex, NULL);
   SCIPsetConsExprNlhdlrBranchscore(scip, nlhdlr, nlhdlrBranchscoreConvex);

   return SCIP_OKAY;
}
