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
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_sum.h"

#define EXPRHDLR_NAME         "var"
#define EXPRHDLR_DESC         "variable expression"
#define EXPRHDLR_HASHKEY     SCIPcalcFibHash(22153.0)

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) ((val) >= (infty1) ? (infty2) : (val))

/** expression data */
struct SCIP_ConsExpr_ExprData
{
   SCIP_VAR*             var;                /**< SCIP variable that this constraint represents */

   SCIP_CONS**           conss;              /**< constraints in which this variable appears */
   int                   nconss;             /**< current number of constraints in conss */
   int                   consssize;          /**< length of conss array */
   SCIP_Bool             consssorted;        /**< is the array of constraints sorted */

   int                   filterpos;          /**< position of eventdata in SCIP's event filter, -1 if not catching events */
};

/** simplifies a variable expression.
 * We replace the variable when fixed by its value
 * If a variable is fixed, (multi)aggregated or more generally, inactive, we replace it with its active counterpart
 * IMPLEMENTATION NOTE: - we follow the general approach of the simplify, where we replace the var expression for its
 * simplified expression only in the current parent. So if we see that there is any performance issue in the simplify
 * we might have to revisit this decision.
 *                      - we build the sum expression by appending variable expressions one at a time. This may be
 * speed-up if we allocate memory for all the variable expressions and build the sum directly.
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_Real constant;
   int nvars;
   int varssize;
   int requsize;
   int i;
   SCIP_CONSEXPR_EXPR* sumexpr;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrVar(conshdlr));

   var = SCIPgetConsExprExprVarVar(expr);
   assert(var != NULL);

   /* if var is active then there is nothing to simplify */
   if( SCIPvarIsActive(var) )
   {
      *simplifiedexpr = expr;
      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureConsExprExpr(*simplifiedexpr);
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
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sumexpr, 0, NULL, NULL, constant) );

   for( i = 0; i < nvars; ++i )
   {
      SCIP_CONSEXPR_EXPR* child;

      SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, &child, vars[i]) );
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, sumexpr, child, coefs[i]) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &child) );
   }

   /* simplify since it might not really be a sum */
   SCIP_CALL( SCIPsimplifyConsExprExprHdlr(scip, conshdlr, sumexpr, simplifiedexpr) );

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "expr_var simplify: <%s> := ", SCIPvarGetName(var));
   SCIPprintConsExprExpr(scip, conshdlr, *simplifiedexpr, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   /* we cannot handle fixings to infinity at the moment (TODO we should) */
   assert(!SCIPisInfinity(scip, REALABS(constant)));

   /* release no longer used sumexpr */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sumexpr) );

   /* free memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** the order of two variable is given by their indices
 * @note: this is affected by permutations in the problem! */
static
SCIP_DECL_CONSEXPR_EXPRCOMPARE(compareVar)
{  /*lint --e{715}*/
   int index1;
   int index2;

   index1 = SCIPvarGetIndex(SCIPgetConsExprExprVarVar(expr1));
   index2 = SCIPvarGetIndex(SCIPgetConsExprExprVarVar(expr2));

   return index1 < index2 ? -1 : index1 == index2 ? 0 : 1;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrVar)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrVar(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** expression handler free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEHDLR(freehdlrVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(exprhdlr != NULL);
   assert(exprhdlrdata != NULL);

   /* free variable to variable expression map */
   assert(SCIPhashmapGetNElements((SCIP_HASHMAP*) (*exprhdlrdata)) == 0);
   SCIPhashmapFree((SCIP_HASHMAP**) exprhdlrdata);
   *exprhdlrdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataVar)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   SCIP_CALL( SCIPallocClearBlockMemory(targetscip, targetexprdata) );
   (*targetexprdata)->filterpos = -1;

   if( mapvar == NULL )
   {
      /* identical mapping */
      assert(targetscip == sourcescip);

      (*targetexprdata)->var = SCIPgetConsExprExprVarVar(sourceexpr);

      SCIP_CALL( SCIPcaptureVar(targetscip, (*targetexprdata)->var) );
   }
   else
   {
      /* call mapvar callback (captures targetvar) */
      SCIP_CALL( (*mapvar)(targetscip, &(*targetexprdata)->var, sourcescip, SCIPgetConsExprExprVarVar(sourceexpr), mapvardata) );
      assert((*targetexprdata)->var != NULL);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataVar)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_HASHMAP* var2expr;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);
   assert(exprdata->var != NULL);
   /* there should be no constraints left that still use this variable */
   assert(exprdata->nconss == 0);
   /* thus, there should also be no variable event catched (via this exprhdlr) */
   assert(exprdata->filterpos == -1);

   var2expr = (SCIP_HASHMAP*) SCIPgetConsExprExprHdlrData(SCIPgetConsExprExprHdlr(expr));
   assert(var2expr != NULL);

   assert(SCIPhashmapExists(var2expr, (void*) exprdata->var));

   /* remove variable expression from the hashmap */
   SCIP_CALL( SCIPhashmapRemove(var2expr, (void*) exprdata->var) );

   SCIP_CALL( SCIPreleaseVar(scip, &exprdata->var) );

   /* free expression data */
   SCIPfreeBlockMemoryArrayNull(scip, &exprdata->conss, exprdata->consssize);
   SCIPfreeBlockMemory(scip, &exprdata);
   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprVarVar(expr) != NULL);

   if( stage == SCIP_CONSEXPRITERATOR_ENTEREXPR )
   {
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName(SCIPgetConsExprExprVarVar(expr)));
   }

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprVarVar(expr) != NULL);

   *val = SCIPgetSolVal(scip, sol, (SCIP_VAR*)SCIPgetConsExprExprVarVar(expr));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffVar)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprVarVar(expr) != NULL);

   /* this should never happen because variable expressions do not have children */
   return SCIP_INVALIDCALL;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert(expr != NULL);

   var = SCIPgetConsExprExprVarVar(expr);
   assert(var != NULL);

   if( intevalvar != NULL )
      *interval = intevalvar(scip, var, intevalvardata);
   else
   {
      SCIP_Real lb;
      SCIP_Real ub;

      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);

      SCIPintervalSetBounds(interval,  /*lint !e666*/
         -infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY, -lb),    /*lint !e666*/
          infty2infty(SCIPinfinity(scip), SCIP_INTERVAL_INFINITY,  ub));   /*lint !e666*/
   }
   return SCIP_OKAY;
}

/** variable hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashVar)
{  /*lint --e{715}*/
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 0);
   assert(hashkey != NULL);

   var = SCIPgetConsExprExprVarVar(expr);
   assert(var != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash((SCIP_Real)SCIPvarGetIndex(var));

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_CONSEXPR_EXPRCURVATURE(curvatureVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(success != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 0);

   /* x -> x is linear, convex, and concave */
   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(monotonicityVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 0);

   *result = SCIP_MONOTONE_INC;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEGRALITY(integralityVar)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   *isintegral = SCIPvarIsIntegral(SCIPgetConsExprExprVarVar(expr));

   return SCIP_OKAY;
}

/** creates the handler for variable expression and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_HASHMAP* var2expr;

   /* initialize hash map to reuse variable expressions for the same variables */
   SCIP_CALL( SCIPhashmapCreate(&var2expr, SCIPblkmem(scip), 100) );

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, 0, evalVar, (SCIP_CONSEXPR_EXPRHDLRDATA*) var2expr) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrVar, freehdlrVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataVar, freedataVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCurvature(scip, consexprhdlr, exprhdlr, curvatureVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrMonotonicity(scip, consexprhdlr, exprhdlr, monotonicityVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntegrality(scip, consexprhdlr, exprhdlr, integralityVar) );

   return SCIP_OKAY;
}

/** creates a variable expression */
SCIP_RETCODE SCIPcreateConsExprExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_VAR*             var                 /**< variable to be stored */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_HASHMAP* var2expr;

   assert(consexprhdlr != NULL);
   assert(expr != NULL);
   assert(var != NULL);

   var2expr = (SCIP_HASHMAP*) SCIPgetConsExprExprHdlrData(SCIPgetConsExprExprHdlrVar(consexprhdlr));
   assert(var2expr != NULL);

   /* check if we have already created a variable expression representing the given variable */
   if( SCIPhashmapExists(var2expr, (void*) var) )
   {
      *expr = (SCIP_CONSEXPR_EXPR*) SCIPhashmapGetImage(var2expr, (void*) var);
      assert(*expr != NULL);

      /* we need to capture the variable expression */
      SCIPcaptureConsExprExpr(*expr);
   }
   else
   {
      /* important to capture variable once since there will be only one variable expression representing this variable */
      SCIP_CALL( SCIPcaptureVar(scip, var) );

      SCIP_CALL( SCIPallocClearBlockMemory(scip, &exprdata) );
      exprdata->var = var;
      exprdata->filterpos = -1;

      SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPgetConsExprExprHdlrVar(consexprhdlr), exprdata, 0, NULL) );

      /* store the variable expression */
      SCIP_CALL( SCIPhashmapInsert(var2expr, (void*) var, (void*) *expr) );
   }

   return SCIP_OKAY;
}

/** gets the variable of a variable expression */
SCIP_VAR* SCIPgetConsExprExprVarVar(
   SCIP_CONSEXPR_EXPR*   expr                /**< variable expression */
   )
{
   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(SCIPgetConsExprExprData(expr) != NULL);

   return SCIPgetConsExprExprData(expr)->var;
}

/** registers event handler to catch variable events on variable
 *
 * Additionally, the given constraint is stored in the data of the variable-expression.
 * When an event occurs, all stored constraints are notified.
 */
SCIP_RETCODE SCIPcatchConsExprExprVarEvent(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< variable expression */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< expr constraint */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
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
      exprdata->consssorted = SCIPcompareConsExprIndex(exprdata->conss[exprdata->nconss-2], exprdata->conss[exprdata->nconss-1]) > 0;

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
   SCIP_CONSEXPR_EXPR*   expr,               /**< variable expression */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< expr constraint */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int pos;

   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(eventhdlr != NULL);
   assert(cons != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
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
         SCIPsortPtr((void**)exprdata->conss, SCIPcompareConsExprIndex, exprdata->nconss);
         exprdata->consssorted = TRUE;
      }

      if( !SCIPsortedvecFindPtr((void**)exprdata->conss, SCIPcompareConsExprIndex, cons, exprdata->nconss, &pos) )
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


/** gives number of constraints for which the expression catches bound change events on the variable */
int SCIPgetConsExprExprVarNConss(
   SCIP_CONSEXPR_EXPR*   expr                /**< variable expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   return exprdata->nconss;
}

/** gives constraints for which the expression catches bound change events on the variable */
SCIP_CONS** SCIPgetConsExprExprVarConss(
   SCIP_CONSEXPR_EXPR*   expr                /**< variable expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   return exprdata->conss;
}
