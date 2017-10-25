/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  CIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cons_expr_cos.c
 * @brief  handler for cosine expressions
 * @author Fabian Wegscheider
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define _USE_MATH_DEFINES   /* to get M_PI on Windows */  /*lint !750 */

#include <string.h>
#include <math.h>
#include "scip/cons_expr_cos.h"
#include "cons_expr_value.h"

/* fundamental expression handler properties */
#define EXPRHDLR_NAME         "cos"
#define EXPRHDLR_DESC         "cosine expression"
#define EXPRHDLR_PRECEDENCE   9200
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(82463.0)

/*
 * Local methods
 */

/* TODO: put your local methods here, and declare them static */

/*
 * Callback methods of expression handler
 */

/** expression handler copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrCos)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrCos(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

/** simplifies a cos expression
 * Evaluates the sine value function when its child is a value expression
 * TODO: add further simplifications
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyCos)
{  /*lint --e{715}*/
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
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr,
            COS(SCIPgetConsExprExprValueValue(child))) );
   }
   else
   {
      *simplifiedexpr = expr;

      /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataCos)
{  /*lint --e{715}*/
   assert(targetscip != NULL);
   assert(targetexprdata != NULL);
   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);
   assert(SCIPgetConsExprExprData(sourceexpr) == NULL);

   *targetexprdata = NULL;

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataCos)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_CONSEXPR_EXPRPRINT(printCos)
{  /*lint --e{715}*/
   assert(expr != NULL);

   switch( stage )
   {
   case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
   {
      /* print function with opening parenthesis */
      SCIPinfoMessage(scip, file, "%s(", EXPRHDLR_NAME);
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

/** expression parse callback */
static
SCIP_DECL_CONSEXPR_EXPRPARSE(parseCos)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* childexpr;

   assert(expr != NULL);

   /* parse child expression from remaining string */
   SCIP_CALL( SCIPparseConsExprExpr(scip, consexprhdlr, string, endstring, &childexpr) );
   assert(childexpr != NULL);

   /* create cosine expression */
   SCIP_CALL( SCIPcreateConsExprExprCos(scip, consexprhdlr, expr, childexpr) );
   assert(*expr != NULL);

   /* release child expression since it has been captured by the cosine expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &childexpr) );

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression (point-) evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalCos)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]) != SCIP_INVALID); /*lint !e777*/

   *val = COS(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[0]));

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffCos)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(idx >= 0 && idx < SCIPgetConsExprExprNChildren(expr));
   assert(SCIPgetConsExprExprValue(expr) != SCIP_INVALID); /*lint !e777*/

   child = SCIPgetConsExprExprChildren(expr)[idx];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);

   *val = -SIN(SCIPgetConsExprExprValue(child));

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalCos)
{  /*lint --e{715}*/
   SCIP_INTERVAL childinterval;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);

   childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[0]);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

   SCIPintervalCos(SCIPinfinity(scip), interval, childinterval);

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_CONSEXPR_EXPRINITSEPA(initSepaCos)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaCos)
{  /*lint --e{715}*/
   assert(expr != NULL);

   SCIPerrorMessage("method of xyz constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropCos)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;
   SCIP_INTERVAL interval;
   SCIP_INTERVAL childbound;
   SCIP_Real newinf;
   SCIP_Real newsup;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(nreductions != NULL);
   assert(SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)) >= -1.0);
   assert(SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)) <= 1.0);

   *nreductions = 0;

   child = SCIPgetConsExprExprChildren(expr)[0];
   assert(child != NULL);

   childbound = SCIPgetConsExprExprInterval(child);
   newinf = childbound.inf;
   newsup = childbound.sup;

   interval = SCIPgetConsExprExprInterval(expr);

   if( !SCIPisInfinity(scip, -newinf) )
   {
      /* l(x) and u(x) are lower/upper bound of child, l(s) and u(s) are lower/upper bound of cos expr
       *
       * if cos(l(x)) < l(s), we are looking for k minimal s.t. a + 2k*pi > l(x) where a = acos(l(s))
       * then the new lower bound is a + 2k*pi
       */
      if( SCIPisLT(scip, COS(newinf), interval.inf) )
      {
         SCIP_Real a = ACOS(interval.inf);
         int k = (int) ceil((newinf - a) / (2.0*M_PI));
         newinf = a + 2.0*M_PI * k;
      }

      /* if cos(l(x)) > u(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) where a = acos(u(s))
       * then the new lower bound is pi - a + 2k*pi
       */
      else if( SCIPisGT(scip, COS(newinf), interval.sup) )
      {
         SCIP_Real a = ACOS(interval.sup);
         int k = (int) ceil((newinf + a) / (2.0*M_PI) - 0.5);
         newinf = M_PI * (2.0*k + 1.0) - a;
      }
   }

   if( !SCIPisInfinity(scip, newsup) )
   {
      /* if cos(u(x)) > u(s), we are looking for k minimal s.t. a + 2k*pi > u(x) - 2*pi where a = acos(u(s))
       * then the new upper bound is a + 2k*pi
       */
      if ( SCIPisGT(scip, COS(newsup), interval.sup) )
      {
         SCIP_Real a = ACOS(interval.sup);
         int k = (int) ceil((newsup - a ) / (2.0*M_PI)) - 1;
         newsup = a + 2.0*M_PI * k;
      }

      /* if cos(u(x)) < l(s), we are looking for k minimal s.t. pi - a + 2k*pi > l(x) - 2*pi where a = acos(l(s))
       * then the new upper bound is pi - a + 2k*pi
       */
      if( SCIPisLT(scip, COS(newsup), interval.inf) )
      {
         SCIP_Real a = ACOS(interval.inf);
         int k = (int) ceil((newsup + a) / (2*M_PI) - 0.5) - 1;
         newsup = M_PI * (2.0*k + 1.0) - a;
      }
   }

   assert(newinf >= childbound.inf);
   assert(newsup <= childbound.sup);
   assert(newinf <= newsup);
   assert(SCIPisGE(scip, COS(newinf), interval.inf));
   assert(SCIPisLE(scip, COS(newinf), interval.sup));
   assert(SCIPisGE(scip, COS(newsup), interval.inf));
   assert(SCIPisLE(scip, COS(newsup), interval.sup));

   SCIPintervalSetBounds(&childbound, newinf, newsup);

   /* try to tighten the bounds of the child node */
   SCIP_CALL( SCIPtightenConsExprExprInterval(scip, child, childbound, force, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** cos hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashCos)
{  /*lint --e{715}*/
   unsigned int childhash;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) == 1);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   *hashkey = EXPRHDLR_HASHKEY;

   assert(SCIPhashmapExists(expr2key, (void*) SCIPgetConsExprExprChildren(expr)[0]));
   childhash = (unsigned int)(size_t) SCIPhashmapGetImage(expr2key, SCIPgetConsExprExprChildren(expr)[0]);

   *hashkey ^= childhash;

   return SCIP_OKAY;
}

/** creates the handler for cos expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrCos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   /* include expression handler */
   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalCos, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrCos, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataCos, freedataCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrParse(scip, consexprhdlr, exprhdlr, parseCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrInitSepa(scip, consexprhdlr, exprhdlr, initSepaCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashCos) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffCos) );

   return SCIP_OKAY;
}

/** creates a cos expression */
SCIP_RETCODE SCIPcreateConsExprExprCos(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPR*   child               /**< single child */
   )
{
   assert(expr != NULL);
   assert(child != NULL);
   assert(SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME) != NULL);

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPfindConsExprExprHdlr(consexprhdlr, EXPRHDLR_NAME), NULL, 1,
         &child) );

   return SCIP_OKAY;
}
