/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_var.c
 * @brief  variable expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_var.h"

#define VAR_HASHKEY     SCIPcalcFibHash(22153)

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrVar)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrVar(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** expression handler free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEHDLR(freehdlrVar)
{
   assert(scip != NULL);
   assert(consexprhdlr != NULL);
   assert(exprhdlr != NULL);
   assert(exprhdlrdata != NULL);

   /* free variable to variable expression map */
   assert(SCIPhashmapGetNEntries((SCIP_HASHMAP*) (*exprhdlrdata)) == 0);
   SCIPhashmapFree((SCIP_HASHMAP**) exprhdlrdata);
   *exprhdlrdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataVar)
{
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   if( mapvar == NULL )
   {
      /* identical mapping: just copy data pointer */
      assert(targetscip == sourcescip);

      *targetexprdata = SCIPgetConsExprExprData(sourceexpr);
      assert(*targetexprdata != NULL);

      SCIP_CALL( SCIPcaptureVar(targetscip, (SCIP_VAR*)*targetexprdata) );
   }
   else
   {
      /* call mapvar callback (captures targetvar) */
      SCIP_CALL( (*mapvar)(targetscip, (SCIP_VAR**)targetexprdata, sourcescip, (SCIP_VAR*)SCIPgetConsExprExprData(sourceexpr), mapvardata) );
      assert(*targetexprdata != NULL);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataVar)
{
   SCIP_HASHMAP* var2expr;
   SCIP_VAR* var;

   assert(expr != NULL);

   var2expr = (SCIP_HASHMAP*) SCIPgetConsExprExprHdlrData(SCIPgetConsExprExprHdlr(expr));
   assert(var2expr != NULL);

   var = (SCIP_VAR*)SCIPgetConsExprExprData(expr);
   assert(var != NULL);
   assert(SCIPhashmapExists(var2expr, (void*) var));

   /* remove variable expression from the hashmap */
   SCIP_CALL( SCIPhashmapRemove(var2expr, (void*) var) );

   SCIP_CALL( SCIPreleaseVar(scip, &var) );

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printVar)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);

   if( stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR )
   {
      SCIPinfoMessage(scip, file, "<%s>", SCIPvarGetName((SCIP_VAR*)SCIPgetConsExprExprData(expr)));
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPREVAL(evalVar)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);

   *val = SCIPgetSolVal(scip, sol, (SCIP_VAR*)SCIPgetConsExprExprData(expr));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalVar)
{
   SCIP_VAR* var;

   assert(expr != NULL);

   var = (SCIP_VAR*) SCIPgetConsExprExprData(expr);
   assert(var != NULL);

   SCIPintervalSetBounds(interval, MIN(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)),
      MAX(SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   return SCIP_OKAY;
}

/** expression reverse propagaton callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropVar)
{
   SCIPerrorMessage("Unexpected call of reverse propagation callback of variable expression handler.");
   SCIPABORT();

   return SCIP_OKAY;
}

/** variable hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashVar)
{
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 0);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   var = (SCIP_VAR*) SCIPgetConsExprExprData(expr);
   assert(var != NULL);

   *hashkey = VAR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(SCIPvarGetIndex(var));

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

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "var", "variable", 0, evalVar, (SCIP_CONSEXPR_EXPRHDLRDATA*) var2expr) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrVar, freehdlrVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataVar, freedataVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashVar) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropVar) );

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
      /* it is very important to capture the variable only once since there will be only one variable expression
       * representing this variable
       */
      SCIP_CALL( SCIPcaptureVar(scip, var) );

      SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPgetConsExprExprHdlrVar(consexprhdlr), (SCIP_CONSEXPR_EXPRDATA*)var, 0, NULL) );

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

   return (SCIP_VAR*)SCIPgetConsExprExprData(expr);
}
