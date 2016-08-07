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
#define ABS_HASHKEY     SCIPcalcFibHash(7187.0)

/*
 * Data structures
 */

struct SCIP_ConsExpr_ExprData
{
   SCIP_ROW*  rowneg;  /**< left tangent z >= -x */
   SCIP_ROW*  rowpos;  /**< right tangent z <= x */
};

/*
 * Local methods
 */


/*
 * Callback methods of expression handler
 */


static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrAbs)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrAbs(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataAbs)
{  /*lint --e{715}*/
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataAbs)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPfreeBlockMemory(scip, &exprdata);
   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printAbs)
{  /*lint --e{715}*/
   assert(expr != NULL);

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

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
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
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = REALABS(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalAbs)
{
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

   SCIPintervalAbs(SCIPinfinity(scip), interval, childinterval);

   return SCIP_OKAY;
}

/** computes both tangent underestimates and secant */
static
SCIP_RETCODE computeCutsAbs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< absolute expression */
   SCIP_ROW**            rowneg,             /**< buffer to store first tangent (might be NULL) */
   SCIP_ROW**            rowpos,             /**< buffer to store second tangent (might be NULL) */
   SCIP_ROW**            secant              /**< buffer to store secant (might be NULL) */
   )
{
   SCIP_VAR* x;
   SCIP_VAR* z;
   SCIP_VAR* vars[2];
   SCIP_Real coefs[2];
   char name[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "abs") == 0);

   x = SCIPgetConsExprExprLinearizationVar(SCIPgetConsExprExprChildren(expr)[0]);
   z = SCIPgetConsExprExprLinearizationVar(expr);
   assert(x != NULL);
   assert(z != NULL);

   vars[0] = z;
   vars[1] = x;
   coefs[0] = -1.0;

   /* compute left tangent -z -x <= 0 */
   if( rowneg != NULL )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "abs_neg_%s", SCIPvarGetName(x));
      coefs[1] = -1.0;
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowneg, conshdlr, name, -SCIPinfinity(scip), 0.0, FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarsToRow(scip, *rowneg, 2, vars, coefs) );
   }

   /* compute right tangent 0 <= -z + x */
   if( rowpos != NULL )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "abs_pos_%s", SCIPvarGetName(x));
      coefs[1] = 1.0;
      SCIP_CALL( SCIPcreateEmptyRowCons(scip, rowpos, conshdlr, name, 0.0, SCIPinfinity(scip), FALSE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddVarsToRow(scip, *rowpos, 2, vars, coefs) );
   }

   /* compute secant */
   if( secant != NULL )
   {
      SCIP_Real lb;
      SCIP_Real ub;

      *secant = NULL;
      lb = SCIPvarGetLbLocal(x);
      ub = SCIPvarGetUbLocal(x);

      /* it does not make sense to add a cut if child variable is unbounded, fixed, non-positive, or non-negative */
      if( !SCIPisInfinity(scip, -lb) && !SCIPisInfinity(scip, ub)
         && SCIPisLT(scip, lb, 0.0) && SCIPisGT(scip, ub, 0.0)
         && !SCIPisEQ(scip, lb, ub) )
      {
         SCIP_Real alpha;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "abs_secant_%s", SCIPvarGetName(x));

         /* let alpha = (|ub|-|lb|) / (ub-lb) then the resulting secant looks like
          *
          * z - |ub| <= alpha * (x - ub)  <=> alpha * ub - |ub| <= -z + alpha * x
          */
         alpha = (REALABS(ub) - REALABS(lb)) / (ub - lb);
         coefs[1] = alpha;

         /* secants are only valid locally */
         SCIP_CALL( SCIPcreateEmptyRowCons(scip, secant, conshdlr, name, (alpha * ub - REALABS(ub)), SCIPinfinity(scip),
               TRUE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddVarsToRow(scip, *secant, 2, vars, coefs) );
      }
   }

   return SCIP_OKAY;
}

/** expression separation initialization callback */
static
SCIP_DECL_CONSEXPR_EXPRINITSEPA(initSepaAbs)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_ROW* secant;

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);
   assert(exprdata->rowneg == NULL);
   assert(exprdata->rowpos == NULL);

   *infeasible = FALSE;
   secant = NULL;

   /* compute initial cuts; do no store the secant in the expression data */
   SCIP_CALL( computeCutsAbs(scip, conshdlr, expr, &exprdata->rowneg, &exprdata->rowpos, &secant) );
   assert(exprdata->rowneg != NULL && exprdata->rowpos != NULL);

   SCIP_CALL( SCIPaddCut(scip, NULL, exprdata->rowneg, FALSE, infeasible) );
   if( !*infeasible )
   {
      SCIP_CALL( SCIPaddCut(scip, NULL, exprdata->rowpos, FALSE, infeasible) );
   }

   /* it might happen that we could not compute a secand (because of fixed or unbounded variables) */
   if( !*infeasible && secant != NULL )
   {
      SCIP_CALL( SCIPaddCut(scip, NULL, secant, FALSE, infeasible) );
   }

   /* release secant */
   if( secant != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &secant) );
   }
   assert(secant == NULL);

   return SCIP_OKAY;
}

/** expression separation deinitialization callback */
static
SCIP_DECL_CONSEXPR_EXPREXITSEPA(exitSepaAbs)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   if( exprdata->rowneg != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &exprdata->rowneg) );
   }

   if( exprdata->rowpos != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &exprdata->rowpos) );
   }

   assert(exprdata->rowneg == NULL);
   assert(exprdata->rowpos == NULL);

   return SCIP_OKAY;
}


/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaAbs)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_ROW* rows[3];
   SCIP_Real efficacy;
   SCIP_Bool infeasible;
   int i;

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   infeasible = FALSE;
   *ncuts = 0;
   *result = SCIP_DIDNOTFIND;

   /* create tangents if it not happened so far (might be possible if the constraint is not 'initial') */
   if( exprdata->rowneg == NULL || exprdata->rowpos == NULL )
   {
      assert(exprdata->rowneg == NULL && exprdata->rowpos == NULL);

      SCIP_CALL( computeCutsAbs(scip, conshdlr, expr, &exprdata->rowneg, &exprdata->rowpos, NULL) );
   }
   assert(exprdata->rowneg != NULL && exprdata->rowpos != NULL);

   rows[0] = exprdata->rowneg;
   rows[1] = exprdata->rowpos;
   SCIP_CALL( computeCutsAbs(scip, conshdlr, expr, NULL, NULL, &rows[2]) );

   for( i = 0; i < 3; ++i )
   {
      if( rows[i] == NULL || SCIProwIsInLP(rows[i]) )
         continue;

      efficacy = -SCIPgetRowSolFeasibility(scip, rows[i], sol);
      if( SCIPisGE(scip, efficacy, minefficacy) )
      {
         SCIP_CALL( SCIPaddCut(scip, sol, rows[i], FALSE, &infeasible) );
         *ncuts += 1;

         if( infeasible )
         {
            *result = SCIP_CUTOFF;
            break;
         }
         else
            *result = SCIP_SEPARATED;
      }
   }

   /* release the secant */
   if( rows[2] != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &rows[2]) );
   }

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
   SCIP_CALL( SCIPsetConsExprExprHdlrInitSepa(scip, consexprhdlr, exprhdlr, initSepaAbs) );
   SCIP_CALL( SCIPsetConsExprExprHdlrExitSepa(scip, consexprhdlr, exprhdlr, exitSepaAbs) );
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
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindConsExprExprHdlr(consexprhdlr, "abs") != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &exprdata) );
   assert(exprdata != NULL);

   BMSclearMemory(exprdata);

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPfindConsExprExprHdlr(consexprhdlr, "abs"), exprdata, 1, &child) );

   return SCIP_OKAY;
}
