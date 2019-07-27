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
 *
 * TODO curvature information that have been computed during the detection
 *      of other nonlinear handler can not be used right now
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
#define DEFAULT_PREFEREXTENDED TRUE
#define DEFAULT_CVXSIGNOMIAL TRUE

/*
 * Data structures
 */

/** nonlinear handler expression data */
struct SCIP_ConsExpr_NlhdlrExprData
{
   SCIP_CONSEXPR_EXPR*   nlexpr;             /**< expression (copy) for which this nlhdlr estimates */
   SCIP_HASHMAP*         nlexpr2origexpr;    /**< mapping of our copied expression to original expression */

   int                   nauxvars;           /**< number of distinct (auxiliary) variables handled */
   SCIP_VAR**            auxvars;            /**< variables we use in expression */
};

/** nonlinear handler data */
struct SCIP_ConsExpr_NlhdlrData
{
   /* parameters */
   SCIP_Bool             detectsum;          /**< whether to run detection when the root of an expression is a sum */
   SCIP_Bool             preferextended;     /**< whether to prefer extended formulations */

   /* parameters that probably should be removed again */
   SCIP_Bool             cvxsignomial;       /**< whether to use convexity check on signomials */
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
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from copied to original expression */
   SCIP_CONSEXPR_EXPR**  nlhdlrexpr,         /**< buffer to store created expr */
   SCIP_CONSEXPR_EXPR*   origexpr,           /**< original expression to be copied */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
)
{
   assert(scip != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(nlhdlrexpr != NULL);
   assert(origexpr != NULL);

   if( SCIPgetConsExprExprNChildren(origexpr) == 0 )
   {
      /* for leaves, do not copy */
      *nlhdlrexpr = origexpr;
      SCIPcaptureConsExprExpr(*nlhdlrexpr);
      SCIP_CALL( SCIPhashmapInsert(nlexpr2origexpr, (void*)*nlhdlrexpr, (void*)origexpr) );
      return SCIP_OKAY;
   }

   /* create copy of expression, but without children */
   SCIP_CALL( SCIPduplicateConsExprExpr(scip, conshdlr, origexpr, nlhdlrexpr, FALSE) );
   assert(*nlhdlrexpr != NULL);  /* copies within the same SCIP must always work */

   /* store the curvature we want to get in the curvature flag of the copied expression
    * it's a bit of a misuse, but once we are done with everything, this is actually correct
    */
   SCIPsetConsExprExprCurvature(*nlhdlrexpr, curv);

   /* remember which the original expression was */
   SCIP_CALL( SCIPhashmapInsert(nlexpr2origexpr, (void*)*nlhdlrexpr, (void*)origexpr) );

   return SCIP_OKAY;
}

/** expand nlhdlr-expression by adding children according to original expression */
static
SCIP_RETCODE nlhdlrExprGrowChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from copied to original expression */
   SCIP_CONSEXPR_EXPR*   nlhdlrexpr,         /**< expression for which to create children */
   SCIP_EXPRCURV*        childrencurv        /**< curvature required for children, or NULL if to set to UNKNOWN */
   )
{
   SCIP_CONSEXPR_EXPR* origexpr;
   SCIP_CONSEXPR_EXPR* child;
   int nchildren;
   int i;

   assert(scip != NULL);
   assert(nlhdlrexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(nlhdlrexpr) == 0);

   origexpr = SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlhdlrexpr);

   nchildren = SCIPgetConsExprExprNChildren(origexpr);
   if( nchildren == 0 )
      return SCIP_OKAY;

   for( i = 0; i < nchildren; ++i )
   {
      SCIP_CALL( nlhdlrExprCreate(scip, conshdlr, nlexpr2origexpr, &child, SCIPgetConsExprExprChildren(origexpr)[i],
         childrencurv != NULL ? childrencurv[i] : SCIP_EXPRCURV_UNKNOWN) );
      SCIP_CALL( SCIPappendConsExprExpr(scip, nlhdlrexpr, child) );
      /* append captures child, so we can release the capture from nlhdlrExprCreate */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &child) );
   }

   assert(SCIPgetConsExprExprNChildren(nlhdlrexpr) == SCIPgetConsExprExprNChildren(origexpr));

   return SCIP_OKAY;
}

/** register violation as branchscore in leafs of nlhdlr-expression */
static
void nlhdlrExprBranchscore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   nlexpr,             /**< nlhdlr-expression */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< map expression copy to original */
   unsigned int          brscoretag,         /**< branchscore tag */
   SCIP_Real             violation           /**< violation */
   )
{
   SCIP_CONSEXPR_EXPR* origexpr;
   int nchildren;
   int c;

   /* TODO expression iterator */

   nchildren = SCIPgetConsExprExprNChildren(nlexpr);
   if( nchildren == 0 )
   {
      origexpr = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr);
      assert(origexpr != NULL);

      SCIPaddConsExprExprBranchScore(scip, origexpr, brscoretag, violation);
      return;
   }

   for( c = 0; c < nchildren; ++c )
      nlhdlrExprBranchscore(scip, SCIPgetConsExprExprChildren(nlexpr)[c], nlexpr2origexpr, brscoretag, violation);
}

/** collect and possibly create variables used by our expression
 *
 * For children where we could not achieve the desired curvature, introduce an auxvar and replace the child by a var-expression that points to this auxvar.
 * Collect all used variables and index them.
 */
static
SCIP_RETCODE collectVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   nlexpr,             /**< nlhdlr-expr */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from our expression copy to original */
   SCIP_HASHMAP*         var2index,          /**< mapping from variables to index */
   int*                  nindices            /**< number of indices */
   )
{
   SCIP_VAR* var;
   SCIP_CONSEXPR_EXPR* origexpr;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_CONSEXPR_EXPR* newchild;
   int i;

   assert(nlexpr != NULL);
   assert(nlexpr2origexpr != NULL);
   assert(var2index != NULL);
   assert(nindices != NULL);

   assert(SCIPgetConsExprExprNChildren(nlexpr) > 0);
   assert(SCIPgetConsExprExprChildren(nlexpr) != NULL);

   /* TODO use iterator */

   for( i = 0; i < SCIPgetConsExprExprNChildren(nlexpr); ++i )
   {
      child = SCIPgetConsExprExprChildren(nlexpr)[i];
      origexpr = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(nlexpr2origexpr, (void*)child);
      assert(origexpr != NULL);

      if( SCIPgetConsExprExprNChildren(child) == 0 )
      {
         if( SCIPgetConsExprExprNChildren(origexpr) > 0 )
         {
            /* having a child that had children in original but not in copy means that we could not achieve the desired curvature
             * thus, replace by a new child that points to the auxvar of the original expression
             */
            SCIP_CALL( SCIPcreateConsExprExprAuxVar(scip, conshdlr, origexpr, &var) );
            assert(var != NULL);
            SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &newchild, var) );  /* this captures newchild once */

            SCIP_CALL( SCIPreplaceConsExprExprChild(scip, nlexpr, i, newchild) );  /* this captures newchild again */

            SCIP_CALL( SCIPhashmapRemove(nlexpr2origexpr, (void*)child) );
            SCIP_CALL( SCIPhashmapInsert(nlexpr2origexpr, (void*)newchild, (void*)origexpr) );

            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &newchild) );  /* because it was captured by both create and replace */

            if( !SCIPhashmapExists(var2index, (void*)var) )
            {
               /* new variable -> new index and remember in hashmap */
               SCIP_CALL( SCIPhashmapInsertInt(var2index, (void*)var, (*nindices)++) );
            }
         }
         else if( SCIPisConsExprExprVar(child) )
         {
            /* if variable, then add to hashmap, if not already there */
            var = SCIPgetConsExprExprVarVar(child);
            if( !SCIPhashmapExists(var2index, (void*)var) )
            {
               SCIP_CALL( SCIPhashmapInsertInt(var2index, (void*)var, (*nindices)++) );
            }
         }
         /* else: it's probably a value-expression, nothing to do */
      }
      else
      {
         SCIP_CALL( collectVars(scip, conshdlr, child, nlexpr2origexpr, var2index, nindices) );
      }
   }

   return SCIP_OKAY;
}

/** looks whether top of given expression looks like a monomial that can have a given curvature
 * e.g., sqrt(x)*sqrt(y) is convex if x,y >= 0 and x and y are convex
 * unfortunately, doesn't work for tls, because i) it's originally sqrt(x*y), and ii) it is expanded into some sqrt(z*y+y)
 * but works for cvxnonsep_nsig
 */
static
SCIP_RETCODE constructExprCheckMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONSEXPR_EXPR*   nlexpr,             /**< nlhdlr-expr to check */
   SCIP_CONSEXPR_EXPR*** stack,              /**< pointer to stack where to add generated leafs */
   int*                  stackpos,           /**< current top position of stack */
   int*                  stacksize,          /**< length of *stack */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from our expression copy to original expression */
   SCIP_Bool             preferextended,     /**< whether we prefer extended formulations */
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
   assert(nlexpr2origexpr != NULL);
   assert(success != NULL);

   *success = FALSE;

   if( SCIPgetConsExprExprHdlr(nlexpr) != SCIPgetConsExprExprHdlrProduct(conshdlr) )
      return SCIP_OKAY;

   expr = SCIPhashmapGetImage(nlexpr2origexpr, (void*)nlexpr);
   assert(expr != NULL);

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

   if( !SCIPexprcurvMonomialInv(SCIPexprcurvMultiply(SCIPgetConsExprExprProductCoef(expr), SCIPgetConsExprExprCurvature(nlexpr)), nfactors, exponents, bounds, curv) )
      goto TERMINATE;

   /* add immediate children to nlexpr
    * some entries in curv actually apply to arguments of pow's, will correct this next
    */
   SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, nlexpr, curv) );
   assert(SCIPgetConsExprExprNChildren(nlexpr) == nfactors);

   /* make sure there is enough space on the stack for all children (or their children) */
   if( *stackpos+1 + nfactors > *stacksize )
   {
      *stacksize = SCIPcalcMemGrowSize(scip, *stackpos+1 + nfactors);
      SCIP_CALL( SCIPreallocBufferArray(scip, stack, *stacksize) );
   }

   /* put children that are not power on stack
    * grow child for children that are power and put this child on stack
    * if preferextended, then require children to be linear
    * unless they are linear, an auxvar will be introduced for them and thus they will be handled as var here
    */
   for( i = 0; i < nfactors; ++i )
   {
      child = SCIPgetConsExprExprChildren(nlexpr)[i];
      assert(child != NULL);

      if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrPower(conshdlr) )
      {
         SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, child, &curv[i]) );
         assert(SCIPgetConsExprExprNChildren(child) == 1);
         child = SCIPgetConsExprExprChildren(child)[0];
      }
      assert(SCIPgetConsExprExprNChildren(child) == 0);

      if( preferextended )
      {
         SCIPsetConsExprExprCurvature(child, SCIP_EXPRCURV_LINEAR);
#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "Extendedform: Require linearity for ");
         SCIPprintConsExprExpr(scip, conshdlr, child, NULL);
         SCIPinfoMessage(scip, NULL, "\n");
#endif
      }

      assert(*stackpos+1 < *stacksize);
      (*stack)[++*stackpos] = child;
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
   SCIP_CONSEXPR_EXPR**  rootnlexpr,         /**< buffer to store created expression */
   SCIP_HASHMAP*         nlexpr2origexpr,    /**< mapping from our expression copy to original expression */
   SCIP_CONSEXPR_EXPR*   rootexpr,           /**< expression */
   SCIP_EXPRCURV         curv                /**< curvature to achieve */
   )
{
   SCIP_CONSEXPR_EXPR* nlexpr;
   SCIP_CONSEXPR_EXPR* origexpr;
   SCIP_EXPRCURV* childcurv;
   int childcurvsize;
   int nchildren;
   SCIP_Bool success;

   /* to do list: expressions where to check whether they can have the desired curvature when taking their children into account */
   SCIP_CONSEXPR_EXPR** stack;
   int stacksize;
   int stackpos;

   assert(scip != NULL);
   assert(rootnlexpr != NULL);
   assert(rootexpr != NULL);
   assert(curv == SCIP_EXPRCURV_CONVEX || curv == SCIP_EXPRCURV_CONCAVE);

   childcurvsize = SCIPgetConsExprExprNChildren(rootexpr);
   SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, childcurvsize) );

   /* create root expression */
   SCIP_CALL( nlhdlrExprCreate(scip, conshdlr, nlexpr2origexpr, rootnlexpr, rootexpr, curv) );

   stacksize = 20;
   SCIP_CALL( SCIPallocBufferArray(scip, &stack, stacksize) );

   stack[0] = *rootnlexpr;
   stackpos = 0;
   while( stackpos >= 0 )
   {
      /* take expression from stack */
      nlexpr = stack[stackpos--];
      assert(nlexpr != NULL);
      assert(SCIPgetConsExprExprNChildren(nlexpr) == 0);

      /* SCIPprintConsExprExpr(scip, conshdlr, nlexpr, NULL);
      SCIPinfoMessage(scip, NULL, "\n"); */

      origexpr = SCIPhashmapGetImage(nlexpr2origexpr, nlexpr);
      assert(origexpr != NULL);
      nchildren = SCIPgetConsExprExprNChildren(origexpr);

      success = TRUE;
      if( !SCIPhasConsExprExprHdlrBwdiff(SCIPgetConsExprExprHdlr(nlexpr)) )
      {
         /* if bwdiff is not implemented, then we could not generate cuts, so stop */
         success = FALSE;
      }
      else if( SCIPgetConsExprExprCurvature(nlexpr) != SCIP_EXPRCURV_UNKNOWN )
      {
         /* TODO here should be a nice system where a number of convexity-detection rules can be tried */
         if( nlhdlrdata->cvxsignomial )
         {
            SCIP_CALL( constructExprCheckMonomial(scip, conshdlr, nlexpr, &stack, &stackpos, &stacksize, nlexpr2origexpr, nlhdlrdata->preferextended, &success) );
            if( success )
               continue;  /* constructExprCheckMonomial will have updated stack */
         }

         if( childcurvsize < nchildren )
         {
            childcurvsize = SCIPcalcMemGrowSize(scip, nchildren);
            SCIP_CALL( SCIPreallocBufferArray(scip, &childcurv, childcurvsize) );
         }

         /* check whether and under which conditions origexpr can have desired curvature */
         SCIP_CALL( SCIPcurvatureConsExprExprHdlr(scip, conshdlr, origexpr, SCIPgetConsExprExprCurvature(nlexpr), &success, childcurv) );
         /* SCIPprintConsExprExpr(scip, conshdlr, nlexpr->origexpr, NULL);
         SCIPinfoMessage(scip, NULL, " is %s? %d\n", SCIPexprcurvGetName(nlexpr->curv), success); */
         if( success )
         {
            /* if origexpr can have curvature curv, then don't treat it as leaf, but include its children */
            SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, nlexpr, childcurv) );
         }
         else if( nlexpr == *rootnlexpr )
         {
            /* if there is no way to ensure curvature for root expression, then we failed */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, rootnlexpr) );
            break;
         }
      }
      else if( nchildren > 0 )
      {
         /* if we don't care about curvature in this subtree anymore (very unlikely),
          * then only continue iterating this subtree to assemble leaf expressions
          */
         SCIP_CALL( nlhdlrExprGrowChildren(scip, conshdlr, nlexpr2origexpr, nlexpr, NULL) );
      }

      /* If more than one child and we prefer extended formulations, then require all children to be linear.
       * Unless they are, auxvars will be introduced and they will be handles as variables, which can be an advantage in the context of extended formulations.
       */
      if( success && nchildren > 1 && nlhdlrdata->preferextended )
      {
         int i;
         for( i = 0; i < nchildren; ++i )
            SCIPsetConsExprExprCurvature(SCIPgetConsExprExprChildren(nlexpr)[i], SCIP_EXPRCURV_LINEAR);
#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "require linearity for children of ");
         SCIPprintConsExprExpr(scip, conshdlr, origexpr, NULL);
         SCIPinfoMessage(scip, NULL, " within ");
         SCIPprintConsExprExpr(scip, conshdlr, rootexpr, NULL);
         SCIPinfoMessage(scip, NULL, "\n");
#endif
      }

      if( success && nchildren > 0 )
      {
         /* add children expressions to to-do list (stack) */
         assert(SCIPgetConsExprExprChildren(nlexpr) != NULL);

         if( stackpos+1 + nchildren > stacksize )
         {
            stacksize = SCIPcalcMemGrowSize(scip, stackpos+1 + nchildren);
            SCIP_CALL( SCIPreallocBufferArray(scip, &stack, stacksize) );
         }
         memcpy(stack + (stackpos+1), SCIPgetConsExprExprChildren(nlexpr), nchildren * sizeof(SCIP_CONSEXPR_EXPR*));
         stackpos += nchildren;
      }
   }
   assert(stackpos == -1);

   SCIPfreeBufferArray(scip, &stack);
   SCIPfreeBufferArray(scip, &childcurv);

   if( *rootnlexpr != NULL )
   {
      SCIP_Bool istrivial = TRUE;
      int i;

      /* if all children do not have children, i.e., are variables, or will be replaced by auxvars, then free
       * also if rootnlexpr has no children, then free
       * TODO could also figure this out by looking at number of expressions we looked at above
       */
      for( i = 0; i < SCIPgetConsExprExprNChildren(*rootnlexpr); ++i )
      {
         if( SCIPgetConsExprExprNChildren(SCIPgetConsExprExprChildren(*rootnlexpr)[i]) > 0 )
         {
            istrivial = FALSE;
            break;
         }
      }

      if( istrivial )
      {
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, rootnlexpr) );
      }
   }

   return SCIP_OKAY;
}

/** creates nonlinear handler expression data structure */
static
SCIP_RETCODE createNlhdlrExprData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_NLHDLREXPRDATA** nlhdlrexprdata, /**< pointer to store nlhdlr expression data */
   SCIP_CONSEXPR_EXPR*   expr,               /**< original expression */
   SCIP_CONSEXPR_EXPR*   nlexpr,             /**< our copy of expression */
   SCIP_HASHMAP*         nlexpr2origexpr     /**< mapping of expression copy to original */
   )
{
   SCIP_HASHMAP* var2index;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(*nlhdlrexprdata == NULL);
   assert(nlexpr != NULL);
   assert(nlexpr2origexpr != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, nlhdlrexprdata) );
   (*nlhdlrexprdata)->nlexpr = nlexpr;
   (*nlhdlrexprdata)->nlexpr2origexpr = nlexpr2origexpr;

   /* make sure there are auxvars and collect all variables */
   SCIP_CALL( SCIPhashmapCreate(&var2index, SCIPblkmem(scip), 10) );   /* TODO replace 10 by number of leafs, which we could count in constructExpr */
   (*nlhdlrexprdata)->nauxvars = 0;
   SCIP_CALL( collectVars(scip, conshdlr, nlexpr, nlexpr2origexpr, var2index, &(*nlhdlrexprdata)->nauxvars) );

   /* assemble auxvars array */
   assert((*nlhdlrexprdata)->nauxvars > 0);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*nlhdlrexprdata)->auxvars, (*nlhdlrexprdata)->nauxvars) );
   for( i = 0; i < SCIPhashmapGetNEntries(var2index); ++i )
   {
      SCIP_HASHMAPENTRY* entry;
      SCIP_VAR* auxvar;
      int idx;

      entry = SCIPhashmapGetEntry(var2index, i);
      if( entry == NULL )
         continue;

      auxvar = (SCIP_VAR*) SCIPhashmapEntryGetOrigin(entry);
      assert(auxvar != NULL);

      idx = SCIPhashmapEntryGetImageInt(entry);
      assert(idx >= 0);
      assert(idx < (*nlhdlrexprdata)->nauxvars);

      (*nlhdlrexprdata)->auxvars[idx] = auxvar;

      SCIPdebugMsg(scip, "auxvar %d: <%s>\n", idx, SCIPvarGetName(auxvar));
   }

   SCIPhashmapFree(&var2index);

#ifdef SCIP_DEBUG
   SCIPprintConsExprExpr(scip, conshdlr, nlexpr, NULL);
   SCIPinfoMessage(scip, NULL, " is handled as %s\n", SCIPexprcurvGetName(SCIPgetConsExprExprCurvature(nlexpr)));
#endif

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

   SCIPfreeBlockMemoryArray(scip, &(*nlhdlrexprdata)->auxvars, (*nlhdlrexprdata)->nauxvars);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*nlhdlrexprdata)->nlexpr) );
   SCIPhashmapFree(&(*nlhdlrexprdata)->nlexpr2origexpr);

   SCIPfreeBlockMemory(scip, nlhdlrexprdata);

   return SCIP_OKAY;
}

/** the detection assumes that the curvature information of the expression has been computed already */
static
SCIP_DECL_CONSEXPR_NLHDLRDETECT(nlhdlrDetectConvex)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_NLHDLRDATA* nlhdlrdata;
   SCIP_CONSEXPR_EXPR* nlexpr = NULL;
   SCIP_HASHMAP* nlexpr2origexpr;

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

   /* ignore sums if > 1 children
    * NOTE: this means we may treat 1+f(x) with f begin a trivial expression here; probably that's ok, just thought to mention it anyway
    */
   if( !nlhdlrdata->detectsum && SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr) && SCIPgetConsExprExprNChildren(expr) > 1 )
      return SCIP_OKAY;

   /* ignore pure constants and variables */
   if( SCIPgetConsExprExprNChildren(expr) == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPhashmapCreate(&nlexpr2origexpr, SCIPblkmem(scip), 10) );  /* TODO 10? */

   if( !*enforcedbelow )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, nlhdlrdata, &nlexpr, nlexpr2origexpr, expr, SCIP_EXPRCURV_CONVEX) );
      if( nlexpr != NULL )
      {
         assert(SCIPgetConsExprExprNChildren(nlexpr) > 0);  /* should not be trivial */

         *enforcedbelow = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPABELOW;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be convex -> can enforce expr <= auxvar\n", (void*)expr);
      }
      else
         SCIPhashmapRemoveAll(nlexpr2origexpr);
   }

   if( !*enforcedabove && nlexpr == NULL )
   {
      SCIP_CALL( constructExpr(scip, conshdlr, nlhdlrdata, &nlexpr, nlexpr2origexpr, expr, SCIP_EXPRCURV_CONCAVE) );
      if( nlexpr != NULL )
      {
         assert(SCIPgetConsExprExprNChildren(nlexpr) > 0);  /* should not be trivial */

         *enforcedabove = TRUE;
         *enforcemethods |= SCIP_CONSEXPR_EXPRENFO_SEPAABOVE;
         *success = TRUE;

         SCIPdebugMsg(scip, "detected expr %p to be concave -> can enforce expr >= auxvar\n", (void*)expr);
      }
   }

   assert(*success || nlexpr == NULL);
   if( !*success )
   {
      SCIPhashmapFree(&nlexpr2origexpr);
      return SCIP_OKAY;
   }

   /* store variable expressions into the expression data of the nonlinear handler */
   SCIP_CALL( createNlhdlrExprData(scip, conshdlr, nlhdlrexprdata, expr, nlexpr, nlexpr2origexpr) );

   return SCIP_OKAY;
}

/** auxiliary evaluation callback */
static
SCIP_DECL_CONSEXPR_NLHDLREVALAUX(nlhdlrEvalAuxConvex)
{ /*lint --e{715}*/
   assert(nlhdlrexprdata != NULL);
   assert(nlhdlrexprdata->nlexpr != NULL);
   assert(auxvalue != NULL);

   SCIP_CALL( SCIPevalConsExprExpr(scip, SCIPfindConshdlr(scip, "expr"), nlhdlrexprdata->nlexpr, sol, 0) );
   *auxvalue = SCIPgetConsExprExprValue(nlhdlrexprdata->nlexpr);

   return SCIP_OKAY;
}

/** estimator callback */
static
SCIP_DECL_CONSEXPR_NLHDLRESTIMATE(nlhdlrEstimateConvex)
{ /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* nlexpr;
   SCIP_EXPRCURV curvature;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);
   assert(SCIPhashmapGetImage(nlhdlrexprdata->nlexpr2origexpr, (void*)nlexpr) == expr);
   assert(rowprep != NULL);
   assert(success != NULL);

   *success = FALSE;

   /* if estimating on non-convex side, then do nothing
    * TODO we are vertex-polyhedral and so should compute something
    */
   curvature = SCIPgetConsExprExprCurvature(nlexpr);
   assert(curvature == SCIP_EXPRCURV_CONVEX || curvature == SCIP_EXPRCURV_CONCAVE);
   if( ( overestimate && curvature == SCIP_EXPRCURV_CONVEX) ||
       (!overestimate && curvature == SCIP_EXPRCURV_CONCAVE) )
      return SCIP_OKAY;

   /* we can skip eval as nlhdlrEvalAux should have been called for same solution before */
   /* SCIP_CALL( nlhdlrExprEval(scip, nlexpr, sol) ); */
   assert(auxvalue == SCIPgetConsExprExprValue(nlexpr)); /* given value (originally from nlhdlrEvalAuxConvex) should coincide with the one stored in nlexpr */  /*lint !e777*/
   /* evaluation error or a too large constant -> skip */
   if( SCIPisInfinity(scip, REALABS(auxvalue)) )
   {
      SCIPdebugMsg(scip, "evaluation error / too large value (%g) for %p\n", auxvalue, (void*)expr);
      return SCIP_OKAY;
   }

   /* compute gradient (TODO: this also reevaluates (soltag=0), which shouldn't be necessary) */
   SCIP_CALL( SCIPcomputeConsExprExprGradient(scip, conshdlr, nlexpr, sol, 0) );

   /* gradient evaluation error -> skip */
   if( SCIPgetConsExprExprDerivative(nlexpr) == SCIP_INVALID ) /*lint !e777*/
   {
      SCIPdebugMsg(scip, "gradient evaluation error for %p\n", (void*)expr);
      return SCIP_OKAY;
   }

   /* add gradient underestimator to rowprep: first contribution of each variable, (x - sol) \nabla f(sol) */
   *success = TRUE;
   for( i = 0; i < nlhdlrexprdata->nauxvars; ++i )
   {
      SCIP_VAR* var = nlhdlrexprdata->auxvars[i];
      SCIP_Real deriv;
      SCIP_Real varval;

      deriv = SCIPgetConsExprExprPartialDiff(scip, conshdlr, nlexpr, var);
      if( deriv == SCIP_INVALID )
      {
         *success = FALSE;
         break;
      }

      varval = SCIPgetSolVal(scip, sol, var);

      SCIPdebugMsg(scip, "add %g * (<%s> - %g) to rowprep\n", deriv, SCIPvarGetName(var), varval);

      /* add deriv * (var - varval) to rowprep */
      SCIP_CALL( SCIPaddRowprepTerm(scip, rowprep, var, deriv) );
      SCIPaddRowprepConstant(rowprep, -deriv * varval);
   }

   if( !*success )
      return SCIP_OKAY;

   /* next add f(sol) */
   SCIPaddRowprepConstant(rowprep, auxvalue);
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
   SCIP_CONSEXPR_EXPR* nlexpr;
   SCIP_Real violation;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(nlhdlrexprdata != NULL);
   assert(success != NULL);

   nlexpr = nlhdlrexprdata->nlexpr;
   assert(nlexpr != NULL);

   assert(SCIPgetConsExprExprAuxVar(expr) != NULL);
   assert(auxvalue == SCIPgetConsExprExprValue(nlexpr)); /* given auxvalue should have been computed by nlhdlrEvalAuxConvex */  /*lint !e777*/

   *success = FALSE;

   /* we separate only convex functions here, so there should be little use for branching
    * if violations are small or there are numerical issues, then we will not have generated a cut, though
    * in that case, we will still branch, that is, register branchscores for all depending var exprs
    */

   /* compute violation */
   if( auxvalue == SCIP_INVALID ) /*lint !e777*/
      violation = SCIPinfinity(scip); /* evaluation error -> we should branch */
   else if( SCIPgetConsExprExprCurvature(nlexpr) == SCIP_EXPRCURV_CONVEX  )
      violation = auxvalue - SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr));
   else
      violation = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprAuxVar(expr)) - auxvalue;

   /* if violation is not on the side that we need to enforce, then no need for branching */
   if( violation <= 0.0 )
      return SCIP_OKAY;

   /* TODO try to figure out which variables appear linear and skip them here */
   nlhdlrExprBranchscore(scip, nlhdlrexprdata->nlexpr, nlhdlrexprdata->nlexpr2origexpr, brscoretag, violation);

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

   SCIP_CALL( SCIPincludeConsExprNlhdlrBasic(scip, consexprhdlr, &nlhdlr, NLHDLR_NAME, NLHDLR_DESC, NLHDLR_PRIORITY, nlhdlrDetectConvex, nlhdlrEvalAuxConvex, nlhdlrdata) );
   assert(nlhdlr != NULL);

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/detectsum",
      "whether to run convexity detection when the root of an expression is a sum",
      &nlhdlrdata->detectsum, FALSE, DEFAULT_DETECTSUM, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/preferextended",
      "whether to prefer extended formulations",
      &nlhdlrdata->preferextended, FALSE, DEFAULT_PREFEREXTENDED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/expr/nlhdlr/" NLHDLR_NAME "/cvxsignomial",
      "whether to use convexity check on signomials",
      &nlhdlrdata->cvxsignomial, FALSE, DEFAULT_CVXSIGNOMIAL, NULL, NULL) );

   SCIPsetConsExprNlhdlrFreeHdlrData(scip, nlhdlr, nlhdlrfreeHdlrDataConvex);
   SCIPsetConsExprNlhdlrCopyHdlr(scip, nlhdlr, nlhdlrCopyhdlrConvex);
   SCIPsetConsExprNlhdlrFreeExprData(scip, nlhdlr, nlhdlrfreeExprDataConvex);
   SCIPsetConsExprNlhdlrSepa(scip, nlhdlr, NULL, NULL, nlhdlrEstimateConvex, NULL);
   SCIPsetConsExprNlhdlrBranchscore(scip, nlhdlr, nlhdlrBranchscoreConvex);

   return SCIP_OKAY;
}
