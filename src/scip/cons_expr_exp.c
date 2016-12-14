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

/**@file   cons_expr_exp.c
 * @brief  exponential expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 *
 * @todo initsepaExp
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_value.h"
#include "scip/cons_expr_exp.h"

#define EXP_PRECEDENCE  85000
#define EXP_HASHKEY     SCIPcalcFibHash(10181.0)

/*
 * Data structures
 */


/*
 * Local methods
 */

/*
 * Callback methods of expression handler
 */

/** simplifies an exp expression.
 * Evaluates the exponential function when its child is a value expression
 * TODO: exp(log(*)) = *
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyExp)
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_CONSHDLR* conshdlr;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   conshdlr = SCIPfindConshdlr(scip, "expr");
   assert(conshdlr != NULL);

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   /* check for value expression */
   if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrValue(conshdlr) )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, exp(SCIPgetConsExprExprValueValue(child))) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrExp)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrExp(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataExp)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataExp)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printExp)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print function with opening parenthesis */
         SCIPinfoMessage(scip, file, "exp(");
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         assert(SCIPgetConsExprExprWalkCurrentChild(expr) == 0);
         break;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      {
         /* print closing parenthesis */
         SCIPinfoMessage(scip, file, ")");
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      default: ;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseExp)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseConsExprExpr(scip, consexprhdlr, string, endstring, &childexpr) );
   assert(childexpr != NULL);

   /* create exponential expression */
   SCIP_CALL( SCIPcreateConsExprExprExp(scip, consexprhdlr, expr, childexpr) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the exponential expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalExp)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = exp(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]));

   return SCIP_OKAY;
}


/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffExp)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(idx >= 0 && idx < SCIPgetConsExprExprNChildren(expr));
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(expr)[idx])), "val") != 0);
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID);

   *val = SCIPgetConsExprExprValue(expr);

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalExp)
{
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

   SCIPintervalExp(SCIPinfinity(scip), interval, childinterval);

   return SCIP_OKAY;
}

/** helper function to separate a given point; needed for proper unittest */
static
SCIP_RETCODE separatePointExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_ROW**            cut                 /**< pointer to store the row */
   )
{
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;
   SCIP_VAR* childvar;
   SCIP_Bool overestimate;
   SCIP_Real violation;
   SCIP_Real refpoint;
   SCIP_Real lincoef;
   SCIP_Real linconstant;
   SCIP_Bool islocal;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "exp") == 0);
   assert(cut != NULL);

   *cut = NULL;

   /* get expression data */
   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);
   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);
   childvar = SCIPgetConsExprExprLinearizationVar(child);
   assert(childvar != NULL);

   /* compute the violation; this determines whether we need to over- or underestimate */
   violation = exp(SCIPgetSolVal(scip, sol, childvar)) - SCIPgetSolVal(scip, sol, auxvar);

   /* check if there is a violation */
   if( SCIPisEQ(scip, violation, 0.0) )
      return SCIP_OKAY;

   /* determine if we need to under- or overestimate */
   overestimate = SCIPisLT(scip, violation, 0.0);
   refpoint = SCIPgetSolVal(scip, sol, childvar);
   lincoef = 0.0;
   linconstant = 0.0;
   success = TRUE;

   /* adjust the reference points */
   refpoint = SCIPisLT(scip, refpoint, SCIPvarGetLbLocal(childvar)) ? SCIPvarGetLbLocal(childvar) : refpoint;
   refpoint = SCIPisGT(scip, refpoint, SCIPvarGetUbLocal(childvar)) ? SCIPvarGetUbLocal(childvar) : refpoint;
   assert(SCIPisLE(scip, refpoint, SCIPvarGetUbLocal(childvar)) && SCIPisGE(scip, refpoint, SCIPvarGetLbLocal(childvar)));

   if( overestimate )
   {
      SCIPaddExpSecant(scip, SCIPvarGetLbLocal(childvar), SCIPvarGetUbLocal(childvar), &lincoef, &linconstant,
         &success);
      islocal = TRUE; /* secants are only valid locally */
   }
   else
   {
      SCIPaddExpLinearization(scip, SCIPvarGetLbLocal(childvar), SCIPvarGetUbLocal(childvar), refpoint, &lincoef, &linconstant,
         &success);
      islocal = FALSE; /* linearization are globally valid */
   }

   /* create cut if it was successful */
   if( success )
   {
      SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, "exp_cut", 0, NULL, NULL,
            overestimate ? -linconstant : -SCIPinfinity(scip),
            overestimate ? SCIPinfinity(scip) : -linconstant,
            islocal, FALSE, FALSE) );

      SCIP_CALL( SCIPaddVarToRow(scip, *cut, auxvar, -1.0) );
      SCIP_CALL( SCIPaddVarToRow(scip, *cut, childvar, lincoef) );
   }

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaExp)
{
   SCIP_ROW* cut;
   SCIP_Bool infeasible;

   cut = NULL;
   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( separatePointExp(scip, conshdlr, expr, sol, &cut) );

   /* failed to compute a cut */
   if( cut == NULL )
      return SCIP_OKAY;

   SCIP_CALL( SCIPmassageConsExprExprCut(scip, &cut, sol, minviolation) );

   /* cut violation or numerics were too bad */
   if( cut == NULL )
      return SCIP_OKAY;

   /* add cut */
   SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
   *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
   *ncuts += 1;

#ifdef SCIP_DEBUG
   SCIPdebugMsg(scip, "add cut with violation %e\n", violation);
   SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** expression reverse propagaton callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropExp)
{
   SCIP_INTERVAL childbound;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(nreductions != NULL);
   assert(SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)) >= 0.0);

   *nreductions = 0;

   if( SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)) <= 0.0 )
   {
      *infeasible = TRUE;
      return SCIP_OKAY;
   }

   /* f = exp(c0) -> c0 = log(f) */
   SCIPintervalLog(SCIPinfinity(scip), &childbound, SCIPgetConsExprExprInterval(expr));

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, SCIPgetConsExprExprChildren(expr)[0], childbound, force, infeasible,
         nreductions) );

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashExp)
{
   unsigned int childhash;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   *hashkey = EXP_HASHKEY;

   assert(SCIPhashmapExists(expr2key, (void*)SCIPgetConsExprExprChildren(expr)[0]));
   childhash = (unsigned int)(size_t)SCIPhashmapGetImage(expr2key, SCIPgetConsExprExprChildren(expr)[0]);

   *hashkey ^= childhash;

   return SCIP_OKAY;
}

/** creates the handler for exponential expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "exp", "exponential function",
         EXP_PRECEDENCE, evalExp, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrExp, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataExp, freedataExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashExp) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffExp) );

   return SCIP_OKAY;
}

/** creates an exponential expression */
SCIP_RETCODE SCIPcreateConsExprExprExp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< single child */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindConsExprExprHdlr(consexprhdlr, "exp") != NULL);

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPfindConsExprExprHdlr(consexprhdlr, "exp"), NULL, 1, &child) );

   return SCIP_OKAY;
}
