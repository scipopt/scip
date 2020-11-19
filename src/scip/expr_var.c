/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_var.c
 * @brief  variable expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include "scip/expr_var.h"
#include "scip/expr_sum.h"

#define EXPRHDLR_NAME         "var"
#define EXPRHDLR_DESC         "variable expression"
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(22153.0)

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))


/** simplifies a variable expression
 *
 * We replace the variable when fixed by its value
 * If a variable is fixed, (multi)aggregated or more generally, inactive, we replace it with its active counterpart
 * IMPLEMENTATION NOTE:
 * - we follow the general approach of the simplify, where we replace the var expression for its
 *   simplified expression only in the current parent. So if we see that there is any performance issue in the simplify
 *   we might have to revisit this decision.
 * - we build the sum expression by appending variable expressions one at a time. This may be
 *   speed-up if we allocate memory for all the variable expressions and build the sum directly.
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifyVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real constant;
   int nvars;
   int varssize;
   int requsize;
   int i;
   SCIP_EXPR* sumexpr;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);

   var = SCIPgetVarExprVar(expr);
   assert(var != NULL);

   /* if var is active then there is nothing to simplify */
   if( SCIPvarIsActive(var) )
   {
      *simplifiedexpr = expr;
      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureExpr(*simplifiedexpr);
      return SCIP_OKAY;
   }

   /* var is not active; obtain active representation var = constant + sum coefs_i vars_i */
   varssize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &vars,  varssize) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, varssize) );

   vars[0]  = var;
   coefs[0] = 1.0;
   constant = 0.0;
   nvars = 1;
   SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, coefs, &nvars, varssize, &constant, &requsize, TRUE) );

   if( requsize > varssize )
   {
      SCIP_CALL( SCIPreallocBufferArray(scip, &vars,  requsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &coefs, requsize) );
      varssize = requsize;
      SCIP_CALL( SCIPgetProbvarLinearSum(scip, vars, coefs, &nvars, varssize, &constant, &requsize, TRUE) );
      assert(requsize <= nvars);
   }

   /* create expression for constant + sum coefs_i vars_i */
   SCIP_CALL( SCIPcreateExprSum(scip, &sumexpr, 0, NULL, NULL, constant, ownerdatacreate, ownerdatacreatedata) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_EXPR* child;

      SCIP_CALL( SCIPcreateExprVar(scip, &child, vars[i], ownerdatacreate, ownerdatacreatedata) );
      SCIP_CALL( SCIPappendExprSumExpr(scip, sumexpr, child, coefs[i]) );
      SCIP_CALL( SCIPreleaseExpr(scip, &child) );
   }

   /* simplify since it might not really be a sum */
   SCIP_CALL( SCIPsimplifyExprShallow(scip, sumexpr, simplifiedexpr, ownerdatacreate, ownerdatacreatedata) );

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "expr_var simplify: <%s> := ", SCIPvarGetName(var));
   SCIPprintExpr(scip, *simplifiedexpr, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   /* we cannot handle fixings to infinity at the moment (TODO we should) */
   assert(!SCIPisInfinity(scip, REALABS(constant)));

   /* release no longer used sumexpr */
   SCIP_CALL( SCIPreleaseExpr(scip, &sumexpr) );

   /* free memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** the order of two variable is given by their indices
 * @note: this is affected by permutations in the problem! */
static
SCIP_DECL_EXPRCOMPARE(compareVar)
{  /*lint --e{715}*/
   int index1;
   int index2;

   index1 = SCIPvarGetIndex(SCIPgetVarExprVar(expr1));
   index2 = SCIPvarGetIndex(SCIPgetVarExprVar(expr2));

   return index1 < index2 ? -1 : index1 == index2 ? 0 : 1;
}

static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrVar)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprHdlrVar(scip) );

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPRCOPYDATA(copydataVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   /* copying into a different SCIP should be handled on the SCIPexprCopy() level (via mapexpr) */
   assert(targetscip == sourcescip);

   var = SCIPgetVarExprVar(sourceexpr);
   assert(var != NULL);

   *targetexprdata = (SCIP_EXPRDATA*)var;

   SCIP_CALL( SCIPcaptureVar(targetscip, var) );

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPRFREEDATA(freedataVar)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   SCIP_CALL( SCIPreleaseVar(scip, (SCIP_VAR**)&exprdata) );

   SCIPexprSetData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_EXPRPRINT(printVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetVarExprVar(expr) != NULL);

   if( stage == SCIP_EXPRITER_ENTEREXPR )
   {
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName(SCIPgetVarExprVar(expr)));
   }

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetVarExprVar(expr) != NULL);

   *val = SCIPgetSolVal(scip, sol, SCIPgetVarExprVar(expr));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetVarExprVar(expr) != NULL);

   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRFWDIFF(fwdiffVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);

   *dot = SCIPexprGetDot(expr);

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);

   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert(expr != NULL);

   var = SCIPgetVarExprVar(expr);
   assert(var != NULL);

   if( intevalvar != NULL )
      *interval = intevalvar(scip, var, intevalvardata);
   else
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);

      SCIPintervalSetBounds(interval,  /*lint !e666*/
         -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb),    /*lint !e666*/
          infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  ub));   /*lint !e666*/
   }

   return SCIP_OKAY;
}

/** variable hash callback */
static
SCIP_DECL_EXPRHASH(hashVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);
   assert(hashkey != NULL);

   var = SCIPgetVarExprVar(expr);
   assert(var != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash((SCIP_Real)SCIPvarGetIndex(var));

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);

   /* x -> x is linear, convex, and concave */
   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicityVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(SCIPexprGetNChildren(expr) == 0);

   *result = SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralityVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   *isintegral = SCIPvarIsIntegral(SCIPgetVarExprVar(expr));

   return SCIP_OKAY;
}

/** creates the handler for variable expression and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprHdlrVar(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeExprHdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, 0, evalVar, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrVar, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataVar, freedataVar);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifyVar);
   SCIPexprhdlrSetCompare(exprhdlr, compareVar);
   SCIPexprhdlrSetPrint(exprhdlr, printVar);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalVar);
   SCIPexprhdlrSetHash(exprhdlr, hashVar);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffVar, fwdiffVar, bwfwdiffVar);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureVar);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicityVar);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralityVar);

   return SCIP_OKAY;
}

/** creates a variable expression */
SCIP_RETCODE SCIPcreateExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_VAR*             var,                /**< variable to be stored */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata  /**< data to pass to ownerdatacreate */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(var != NULL);

   /* capture the variable so that it doesn't disappear while the expr still points to it */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   exprdata = (SCIP_EXPRDATA*)var;

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPgetExprHdlrVar(scip), exprdata, 0, NULL, ownerdatacreate, ownerdatacreatedata) );

   return SCIP_OKAY;
}

/* from pub_expr.h */

/** gets the variable of a variable expression */
SCIP_VAR* SCIPgetVarExprVar(
   SCIP_EXPR*            expr                /**< variable expression */
   )
{
   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(SCIPexprGetData(expr) != NULL);

   return (SCIP_VAR*)SCIPexprGetData(expr);
}

#if !1  // FIXME move into cons_nonlinear
/** registers event handler to catch variable events on variable
 *
 * Additionally, the given constraint is stored in the data of the variable-expression.
 * When an event occurs, all stored constraints are notified.
 */
SCIP_RETCODE SCIPcatchConsExprExprVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< variable expression */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< expr constraint */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

#ifndef NDEBUG
   /* assert that constraint does not double-catch variable */
   {
      int i;
      for( i = 0; i < exprdata->nconss; ++i )
         assert(exprdata->conss[i] != cons);
   }
#endif

   /* append cons to exprdata->conss */
   if( exprdata->consssize < exprdata->nconss + 1 )
   {
      int newsize = SCIPcalcMemGrowSize(scip, exprdata->nconss+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &exprdata->conss, exprdata->consssize, newsize) );
      exprdata->consssize = newsize;
   }
   exprdata->conss[exprdata->nconss++] = cons;
   /* we're not capturing the constraint here to avoid circular references */

   /* updated sorted flag */
   if( exprdata->nconss <= 1 )
      exprdata->consssorted = TRUE;
   else if( exprdata->consssorted )
      exprdata->consssorted = SCIPcompareIndexConsNonlinear(exprdata->conss[exprdata->nconss-2], exprdata->conss[exprdata->nconss-1]) > 0;

   /* catch variable events, if not done so yet (first constraint) */
   if( exprdata->filterpos < 0 )
   {
      SCIP_EVENTTYPE eventtype;

      assert(exprdata->nconss == 1);

      eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED;

      SCIP_CALL( SCIPcatchVarEvent(scip, exprdata->var, eventtype, eventhdlr, (SCIP_EVENTDATA*)expr, &exprdata->filterpos) );
      assert(exprdata->filterpos >= 0);
   }

   return SCIP_OKAY;
}

/** unregisters event handler to catch variable events on variable
 *
 * The given constraint is removed from the constraints array in the data of the variable-expression.
 * If this was the last constraint, then the event handler is unregistered for this variable.
 */
SCIP_RETCODE SCIPdropConsExprExprVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*   expr,               /**< variable expression */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< expr constraint */
   )
{
   SCIP_EXPRDATA* exprdata;
   int pos;

   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);
   assert(exprdata->nconss > 0);

   if( exprdata->conss[exprdata->nconss-1] == cons )
   {
      pos = exprdata->nconss-1;
   }
   else
   {
      if( !exprdata->consssorted )
      {
         SCIPsortPtr((void**)exprdata->conss, SCIPcompareIndexConsNonlinear, exprdata->nconss);
         exprdata->consssorted = TRUE;
      }

      if( !SCIPsortedvecFindPtr((void**)exprdata->conss, SCIPcompareIndexConsNonlinear, cons, exprdata->nconss, &pos) )
      {
         SCIPerrorMessage("Constraint <%s> not in constraint array of expression for variable <%s>\n", SCIPconsGetName(cons), SCIPvarGetName(exprdata->var));
         return SCIP_ERROR;
      }
      assert(pos >= 0 && pos < exprdata->nconss);
   }
   assert(exprdata->conss[pos] == cons);

   /* move last constraint into position of removed constraint */
   if( pos < exprdata->nconss-1 )
   {
      exprdata->conss[pos] = exprdata->conss[exprdata->nconss-1];
      exprdata->consssorted = FALSE;
   }
   --exprdata->nconss;

   /* drop variable events if that was the last constraint */
   if( exprdata->nconss == 0 )
   {
      SCIP_EVENTTYPE eventtype;

      assert(exprdata->filterpos >= 0);

      eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED;

      SCIP_CALL( SCIPdropVarEvent(scip, exprdata->var, eventtype, eventhdlr, (SCIP_EVENTDATA*)expr, exprdata->filterpos) );
      exprdata->filterpos = -1;
   }

   return SCIP_OKAY;
}

/** returns whether the variable events on variable are catched */
SCIP_Bool SCIPisConsExprExprVarEventCatched(
   SCIP_EXPR*   expr                /**< variable expression */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   return exprdata->filterpos >= 0;
}

/** gives number of constraints for which the expression catches bound change events on the variable */
int SCIPgetConsExprExprVarNConss(
   SCIP_EXPR*   expr                /**< variable expression */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   return exprdata->nconss;
}

/** gives constraints for which the expression catches bound change events on the variable */
SCIP_CONS** SCIPgetConsExprExprVarConss(
   SCIP_EXPR*   expr                /**< variable expression */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   return exprdata->conss;
}
#endif
