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

/**@file   cons_expr_abs.c
 * @brief  absolute expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_abs.h"

#define ABS_PRECEDENCE  70000
#define ABS_HASHKEY     SCIPcalcFibHash(7187)
/*
 * Data structures
 */


/*
 * Local methods
 */


/*
 * Callback methods of expression handler
 */


static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrAbs)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrAbs(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataAbs)
{
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataAbs)
{
   assert(expr != NULL);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printAbs)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print function with opening parenthesis */
         SCIPinfoMessage(scip, file, "abs(");
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

      default: ;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseAbs)
{
   SCIP_CONSEXPR_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseConsExprExpr(scip, consexprhdlr, string, endstring, &childexpr) );
   assert(childexpr != NULL);

   /* create absolute expression */
   SCIP_CALL( SCIPcreateConsExprExprAbs(scip, consexprhdlr, expr, childexpr) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the absolute expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPREVAL(evalAbs)
{
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID);

   *val = REALABS(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalAbs)
{
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) == NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

   SCIPintervalAbs(SCIPinfinity(scip), interval, childinterval);

   return SCIP_OKAY;
}

/** helper function to separate a given point; needed for proper unittest */
static
SCIP_RETCODE separatePointAbs(
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
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "abs") == 0);
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
   violation = REALABS(SCIPgetSolVal(scip, sol, childvar)) - SCIPgetSolVal(scip, sol, auxvar);

   /* check if there is a violation */
   if( SCIPisEQ(scip, violation, 0.0) )
      return SCIP_OKAY;

   /* determine whether we need to under- or overestimate */
   overestimate = SCIPisLT(scip, violation, 0.0);

   refpoint = SCIPgetSolVal(scip, sol, childvar);
   lincoef = 0.0;
   linconstant = 0.0;
   success = FALSE;
   lb = SCIPvarGetLbLocal(childvar);
   ub = SCIPvarGetUbLocal(childvar);

   if( overestimate )
   {
      printf("OVERESTIMATE\n");

      if( !SCIPisInfinity(scip, REALABS(lb)) && !SCIPisInfinity(scip, REALABS(ub)) )
      {
         printf("lb = %e ub = %e\n", lb, ub);
         success = TRUE;
         lincoef = (REALABS(ub) - REALABS(lb)) / (ub - lb);
         linconstant = REALABS(ub) - lincoef * ub;

         printf("lincoef = %e linconst = %e\n", lincoef, linconstant);
      }
   }
   else
   {
      printf("UNDERESTIMATE\n");

      /* this case can be removed after implementing the INITLP callback of the expression handler */
      success = TRUE;
      lincoef = SCIPisLE(scip, refpoint, 0.0) ? -1.0 : 1.0;
      linconstant = 0.0;
   }

   /* create cut */
   if( success )
   {
      /* create cut */
      SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, "abs_cut", 0, NULL, NULL,
            overestimate ? -linconstant : -SCIPinfinity(scip),
            overestimate ? SCIPinfinity(scip) : -linconstant,
            FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarToRow(scip, *cut, childvar, lincoef) );
      SCIP_CALL( SCIPaddVarToRow(scip, *cut, auxvar, -1) );
   }

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaAbs)
{
   SCIP_ROW* cut;
   SCIP_Real efficacy;

   cut = NULL;
   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( separatePointAbs(scip, conshdlr, expr, sol, &cut) );

   /* failed to compute a cut */
   if( cut == NULL )
      return SCIP_OKAY;

   efficacy = -SCIPgetRowSolFeasibility(scip, cut, sol);

   /* add cut if it is violated */
   if( SCIPisGE(scip, efficacy, minefficacy) )
   {
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPaddCut(scip, NULL, cut, FALSE, &infeasible) );
      *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;

#ifdef SCIP_DEBUG
      SCIPdebugMessage("add cut with efficacy %e\n", efficacy);
      SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif
   }

   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropAbs)
{
   SCIP_INTERVAL childbound;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(nreductions != NULL);
   assert(SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)) >= 0.0);

   *nreductions = 0;

   /* handle absolute expression as identity if child expression is already non-negative */
   if( SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]).inf >= 0.0 )
   {
      SCIPintervalSetBounds(&childbound, SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)),
         SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)));
   }
   /* handle absolute expression as -identity if child expression is already non-positive */
   else if( SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]).sup <= 0.0 )
   {
      assert(-SCIPgetConsExprExprInterval(expr).sup <= -SCIPgetConsExprExprInterval(expr).inf);
      SCIPintervalSetBounds(&childbound, -SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)),
         -SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)));
   }
   /* f = abs(c0) => c0 = -f union f */
   else
   {
      SCIPintervalSetBounds(&childbound, -SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)),
         SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)));
   }

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, SCIPgetConsExprExprChildren(expr)[0], childbound, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** expression hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashAbs)
{
   unsigned int childhash;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   *hashkey = ABS_HASHKEY;

   assert(SCIPhashmapExists(expr2key, (void*) SCIPgetConsExprExprChildren(expr)[0]));
   childhash = (unsigned int)(size_t) SCIPhashmapGetImage(expr2key, SCIPgetConsExprExprChildren(expr)[0]);

   *hashkey ^= childhash;

   return SCIP_OKAY;
}

/** creates the handler for absolute expression and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrAbs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "abs", "absolute expression",
         ABS_PRECEDENCE, evalAbs, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrAbs, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataAbs, freedataAbs) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printAbs) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseAbs) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalAbs) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaAbs) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashAbs) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropAbs) );

   return SCIP_OKAY;
}

/** creates an absolute expression */
SCIP_RETCODE SCIPcreateConsExprExprAbs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< single child */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindConsExprExprHdlr(consexprhdlr, "abs") != NULL);

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPfindConsExprExprHdlr(consexprhdlr, "abs"), NULL, 1, &child) );

   return SCIP_OKAY;
}
