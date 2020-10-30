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

/**@file   scip_expr.c
 * @brief  public functions to work with algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#include "scip/scip_expr.h"
#include "scip/expr.h"
#include "scip/set.h"

/**@name Expression Handler Methods */
/**@{ */

/** creates the handler for an expression handler and includes it into SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeExprHdlr(
   SCIP*                 scip,         /**< SCIP data structure */
   SCIP_EXPRHDLR**       exprhdlr,     /**< buffer where to store created expression handler */
   const char*           name,         /**< name of expression handler (must not be NULL) */
   const char*           desc,         /**< description of expression handler (can be NULL) */
   unsigned int          precedence,   /**< precedence of expression operation (used for printing) */
   SCIP_DECL_EXPREVAL((*eval)),        /**< point evaluation callback (must not be NULL) */
   SCIP_EXPRHDLRDATA*    data          /**< data of expression handler (can be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(exprhdlr != NULL);
   assert(name != NULL);

   SCIP_CALL( SCIPallocClearMemory(scip, exprhdlr) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*exprhdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*exprhdlr)->desc, desc, strlen(desc)+1) );
   }

   (*exprhdlr)->precedence = precedence;
   (*exprhdlr)->eval = eval;
   (*exprhdlr)->data = data;

   /* create clocks */
   SCIP_CALL( SCIPcreateClock(scip, &(*exprhdlr)->estimatetime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*exprhdlr)->proptime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*exprhdlr)->intevaltime) );
   SCIP_CALL( SCIPcreateClock(scip, &(*exprhdlr)->simplifytime) );

   SCIP_CALL( SCIPsetIncludeExprhdlr(scip->set, exprhdlr) );

   return SCIP_OKAY;
}

/** gives expression handlers */
SCIP_EXPRHDLR** SCIPgetExprHdlrs(
   SCIP*                      scip           /**< SCIP data structure */
)
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->exprhdlrs;
}

/** gives number of expression handlers */
int SCIPgetNExprHdlrs(
   SCIP*                      scip           /**< SCIP data structure */
)
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->nexprhdlrs;
}

/** returns an expression handler of a given name (or NULL if not found) */
SCIP_EXPRHDLR* SCIPfindExprHdlr(
   SCIP*                      scip,          /**< SCIP data structure */
   const char*                name           /**< name of expression handler */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return SCIPsetFindExprhdlr(scip->set, name);
}

/** returns expression handler for variable expressions (or NULL if not included) */
SCIP_EXPRHDLR* SCIPgetExprHdlrVar(
   SCIP*                      scip           /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->exprhdlrvar;
}

/** returns expression handler for constant value expressions (or NULL if not included) */
SCIP_EXPRHDLR* SCIPgetExprHdlrValue(
   SCIP*                      scip           /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->exprhdlrval;
}

/** returns expression handler for sum expressions (or NULL if not included) */
SCIP_EXPRHDLR* SCIPgetExprHdlrSum(
   SCIP*                      scip           /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->exprhdlrsum;
}

/** returns expression handler for product expressions (or NULL if not included) */
SCIP_EXPRHDLR* SCIPgetExprHdlrProduct(
   SCIP*                      scip           /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->exprhdlrproduct;
}

/** returns expression handler for power expressions (or NULL if not included) */
SCIP_EXPRHDLR* SCIPgetExprHdlrPower(
   SCIP*                      scip           /**< SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   return scip->set->exprhdlrpow;
}

/**@} */


/**@name Expression Methods */
/**@{ */

/** creates and captures an expression with given expression data and children */
SCIP_RETCODE SCIPcreateExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data (expression assumes ownership) */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children            /**< children (can be NULL if nchildren is 0) */
   )
{
   assert(scip != NULL);
   assert(scip->set != NULL);

   SCIP_CALL( SCIPexprCreate(scip->set, scip->mem->probmem, expr, exprhdlr, exprdata, nchildren, children) );

   return SCIP_OKAY;
}

/** creates and captures an expression with given expression data and up to two children */
SCIP_RETCODE SCIPcreateExpr2(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data */
   SCIP_EXPR*            child1,             /**< first child (can be NULL) */
   SCIP_EXPR*            child2              /**< second child (can be NULL) */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(exprhdlr != NULL);

   if( child1 != NULL && child2 != NULL )
   {
      SCIP_EXPR* pair[2];
      pair[0] = child1;
      pair[1] = child2;

      SCIP_CALL( SCIPcreateExpr(scip, expr, exprhdlr, exprdata, 2, pair) );
   }
   else if( child2 == NULL )
   {
      SCIP_CALL( SCIPcreateExpr(scip, expr, exprhdlr, exprdata, child1 == NULL ? 0 : 1, &child1) );
   }
   else
   {
      /* child2 != NULL, child1 == NULL */
      SCIP_CALL( SCIPcreateExpr(scip, expr, exprhdlr, exprdata, 1, &child2) );
   }

   return SCIP_OKAY;
}

/** creates and captures an expression representing a quadratic function */
SCIP_RETCODE SCIPcreateExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nlinvars,           /**< number of linear terms */
   SCIP_VAR**            linvars,            /**< array with variables in linear part */
   SCIP_Real*            lincoefs,           /**< array with coefficients of variables in linear part */
   int                   nquadterms,         /**< number of quadratic terms */
   SCIP_VAR**            quadvars1,          /**< array with first variables in quadratic terms */
   SCIP_VAR**            quadvars2,          /**< array with second variables in quadratic terms */
   SCIP_Real*            quadcoefs           /**< array with coefficients of quadratic terms */
   )
{
   SCIP_EXPR** children;
   SCIP_Real* coefs;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(nlinvars == 0 || (linvars != NULL && lincoefs != NULL));
   assert(nquadterms == 0 || (quadvars1 != NULL && quadvars2 != NULL && quadcoefs != NULL));

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &children, nquadterms + nlinvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nquadterms + nlinvars) );

   /* create children for quadratic terms */
   for( i = 0; i < nquadterms; ++i )
   {
      assert(quadvars1 != NULL && quadvars1[i] != NULL);
      assert(quadvars2 != NULL && quadvars2[i] != NULL);

      /* quadratic term */
      if( quadvars1[i] == quadvars2[i] )
      {
         SCIP_EXPR* xexpr;

         /* create variable expression; intentionally not using createExprVar here,
          * since expression created here is not part of a constraint (they will be copied when a constraint is created)
          */
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, &xexpr, quadvars1[i]) );

         /* create pow expression */
         SCIP_CALL( SCIPcreateConsExprExprPow(scip, &children[i], xexpr, 2.0) );

         /* release variable expression; note that the variable expression is still captured by children[i] */
         SCIP_CALL( SCIPreleaseExpr(scip, &xexpr) );
      }
      else /* bilinear term */
      {
         SCIP_EXPR* exprs[2];

         /* create variable expressions; intentionally not using createExprVar here,
          * since expression created here is not part of a constraint (they will be copied when a constraint is created)
          */
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, &exprs[0], quadvars1[i]) );
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, &exprs[1], quadvars2[i]) );

         /* create product expression */
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, &children[i], 2, exprs, 1.0) );

         /* release variable expressions; note that the variable expressions are still captured by children[i] */
         SCIP_CALL( SCIPreleaseExpr(scip, &exprs[1]) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exprs[0]) );
      }

      /* store coefficient */
      coefs[i] = quadcoefs[i];
   }

   /* create children for linear terms */
   for( i = 0; i < nlinvars; ++i )
   {
      assert(linvars != NULL && linvars[i] != NULL);

      /* create variable expression; intentionally not using createExprVar here,
       * since expression created here is not part of a constraint (they will be copied when a constraint is created);
       * release variable expression after the sum expression has been created
       */
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, &children[nquadterms + i], linvars[i]) );

      /* store coefficient */
      coefs[nquadterms + i] = lincoefs[i];
   }

   /* create sum expression */
   SCIP_CALL( SCIPcreateConsExprExprSum(scip, expr, nquadterms + nlinvars, children, coefs, 0.0) );

   /* release children */
   for( i = 0; i < nquadterms + nlinvars; ++i )
   {
      assert(children[i] != NULL);
      SCIP_CALL( SCIPreleaseExpr(scip, &children[i]) );
   }

   /* free memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &children);

   return SCIP_OKAY;
}

/** creates and captures an expression representing a monomial */
SCIP_RETCODE SCIPcreateExprMonomial(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nfactors,           /**< number of factors in monomial */
   SCIP_VAR**            vars,               /**< variables in the monomial */
   SCIP_Real*            exponents           /**< exponent in each factor, or NULL if all 1.0 */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(nfactors >= 0);

   /* return 1 as constant expression if there are no factors */
   if( nfactors == 0 )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, expr, 1.0) );
   }
   else if( nfactors == 1 )
   {
      /* only one factor and exponent is 1 => return factors[0] */
      if( exponents == NULL || exponents[0] == 1.0 )
      {
         /* intentionally not using createExprVar here, since expression created here is not part of
          * a constraint (they will be copied when a constraint is created)
          */
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, expr, vars[0]) );
      }
      else
      {
         SCIP_EXPR* varexpr;

         /* create variable and power expression; intentionally not using createExprVar here,
          * since expression created here is not part of a constraint (they will be copied when a constraint is created)
          */
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, &varexpr, vars[0]) );
         SCIP_CALL( SCIPcreateConsExprExprPow(scip, expr, varexpr, exponents[0]) );
         SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );
      }
   }
   else
   {
      SCIP_EXPR** children;
      int i;

      /* allocate memory to store the children */
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nfactors) );

      /* create children */
      for( i = 0; i < nfactors; ++i )
      {
         /* check whether to create a power expression or not, i.e., exponent == 1 */
         if( exponents == NULL || exponents[i] == 1.0 )
         {
            SCIP_CALL( SCIPcreateConsExprExprVar(scip, &children[i], vars[i]) );
         }
         else
         {
            SCIP_EXPR* varexpr;

            /* create variable and pow expression */
            SCIP_CALL( SCIPcreateConsExprExprVar(scip, &varexpr, vars[i]) );
            SCIP_CALL( SCIPcreateConsExprExprPow(scip, &children[i], varexpr, exponents[i]) );
            SCIP_CALL( SCIPreleaseExpr(scip, &varexpr) );
         }
      }

      /* create product expression */
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, expr, nfactors, children, 1.0) );

      /* release children */
      for( i = 0; i < nfactors; ++i )
      {
         assert(children[i] != NULL);
         SCIP_CALL( SCIPreleaseExpr(scip, &children[i]) );
      }

      /* free memory */
      SCIPfreeBufferArray(scip, &children);
   }

   return SCIP_OKAY;
}

/** appends child to the children list of expr */
SCIP_RETCODE SCIPappendExprChild(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPR*            child               /**< expression to be appended */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprAppendChild(scip->set, scip->mem->probmem, expr, child) );

   return SCIP_OKAY;
}

/** overwrites/replaces a child of an expressions
 *
 * @note the old child is released and the newchild is captured, unless they are the same (=same pointer)
 */
SCIP_RETCODE SCIPreplaceExprChild(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPR*              expr,             /**< expression which is going to replace a child */
   int                     childidx,         /**< index of child being replaced */
   SCIP_EXPR*              newchild          /**< the new child */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprReplaceChild(scip->set, scip->stat, scip->mem->probmem, expr, childidx, newchild) );

   return SCIP_OKAY;
}

/** remove all children of expr
 *
 * @attention only use if you really know what you are doing
 */
SCIP_RETCODE SCIPremoveExprChildren(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprRemoveChildren(scip->set, scip->stat, scip->mem->probmem, expr) );

   return SCIP_OKAY;
}

/** duplicates the given expression (including children) */
SCIP_RETCODE SCIPcopyExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_DECL_EXPR_MAPVAR((*mapvar)),         /**< variable mapping function, or NULL for identity mapping */
   void*                 mapvardata,         /**< data of variable mapping function */
   SCIP_DECL_EXPR_MAPEXPR((*mapexpr)),       /**< expression mapping function, or NULL for creating new expressions */
   void*                 mapexprdata,        /**< data of expression mapping function */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call on expression copy to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata, /**< data to pass to ownerdatacreate */
   SCIP_DECL_EXPR_OWNERDATAFREE((*ownerdatafree)),     /**< function to call when freeing expression, e.g., to free ownerdata */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprCopy(scip->set, scip->stat, scip->mem->probmem, scip->set, scip->mem->probmem, expr, copyexpr, mapvar, mapvardata, mapexpr, mapexprdata, ownerdatacreate, ownerdatacreatedata, ownerdatafree) );

   return SCIP_OKAY;
}

/** duplicates the given expression without its children */
SCIP_RETCODE SCIPcopyExprShallow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr            /**< buffer to store (shallow) duplicate of expr */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   /* copy expression data */
   SCIP_EXPRDATA* exprdatacopy = NULL;
   if( SCIPexprGetData(expr) != NULL )
   {
      assert(expr->exprhdlr->copydata != NULL);
      SCIP_CALL( expr->exprhdlr->copydata(scip, expr->exprhdlr, &exprdatacopy, scip, expr, NULL, NULL) );
   }

   /* create expression with same handler and copied data, but without children */
   SCIP_CALL( SCIPexprCreate(scip->set, scip->mem->probmem, copyexpr, expr->exprhdlr, exprdatacopy, 0, NULL) );

   return SCIP_OKAY;
}

/** creates an expression from a string
 *
 * We specify the grammar that defines the syntax of an expression. Loosely speaking, a Base will be any "block",
 * a Factor is a Base to a power, a Term is a product of Factors and an Expression is a sum of terms
 * The actual definition:
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * Term       -> Factor { ("*" | "/" ) Factor }
 * Factor     -> Base [ "^" "number" | "^(" "number" ")" ]
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 * where [a|b] means a or b or none, (a|b) means a or b, {a} means 0 or more a.
 *
 * Note that Op and OpExpression are undefined. Op corresponds to the name of an expression handler and
 * OpExpression to whatever string the expression handler accepts (through its parse method).
 *
 * See also @ref parseExpr in expr.c.
 */
SCIP_RETCODE SCIPparseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer to store the expr parsed */
   const char*           exprstr,            /**< string with the expr to parse */
   const char**          finalpos            /**< buffer to store the position of exprstr where we finished reading, or NULL if not of interest */
)
{
   const char* finalpos_;
   SCIP_RETCODE retcode;
   SCIP_HASHMAP* vartoexprvarmap;

   SCIP_CALL( SCIPhashmapCreate(&vartoexprvarmap, SCIPblkmem(scip), 5 * SCIPgetNVars(scip)) );

   /* if parseExpr fails, we still want to free hashmap */
   retcode = SCIPexprParse(scip->set, scip->mem->probmem, vartoexprvarmap, exprstr, &finalpos_, expr);

   SCIPhashmapFree(&vartoexprvarmap);

   if( finalpos != NULL )
      *finalpos = finalpos_;

   return retcode;
}

/** captures an expression (increments usage count) */
void SCIPcaptureExpr(
   SCIP_EXPR*            expr                /**< expression to be captured */
   )
{
   SCIPexprCapture(expr);
}

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_RETCODE SCIPreleaseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr                /**< pointer to expression to be released */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprRelease(scip->set, scip->stat, scip->mem->probmem, expr) );

   return SCIP_OKAY;
}

/** returns whether an expression is a variable expression */
SCIP_Bool SCIPisExprVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(scip != NULL);

   return SCIPexprIsVar(scip->set, expr);
}

/** returns whether an expression is a value expression */
SCIP_Bool SCIPisExprValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(scip != NULL);

   return SCIPexprIsValue(scip->set, expr);
}

/** print an expression as info-message */
SCIP_RETCODE SCIPprintExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be printed */
   FILE*                 file                /**< file to print to, or NULL for stdout */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprPrint(scip->set, scip->stat, scip->mem->blkmem, scip->messagehdlr, expr, file) );

   return SCIP_OKAY;
}

/** initializes printing of expressions in dot format to a give FILE* pointer */
SCIP_RETCODE SCIPprintExprDotInit(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata,        /**< buffer to store dot printing data */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_EXPRPRINT_WHAT     whattoprint       /**< info on what to print for each expression */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprPrintDotInit(scip->set, scip->stat, scip->mem->probmem, printdata, file, whattoprint) );

   return SCIP_OKAY;
}

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_RETCODE SCIPprintExprDotInit2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata,        /**< buffer to store dot printing data */
   const char*             filename,         /**< name of file to print to */
   SCIP_EXPRPRINT_WHAT     whattoprint       /**< info on what to print for each expression */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprPrintDotInit2(scip->set, scip->stat, scip->mem->probmem, printdata, filename, whattoprint) );

   return SCIP_OKAY;
}

/** main part of printing an expression in dot format */
SCIP_RETCODE SCIPprintExprDot(
   SCIP*                  scip,              /**< SCIP data structure */
   SCIP_EXPRPRINTDATA*    printdata,         /**< data as initialized by \ref SCIPprintExprDotInit() */
   SCIP_EXPR*             expr               /**< expression to be printed */
   )
{
   assert(scip != NULL);

   SCIP_CALL( SCIPexprPrintDot(scip->set, scip->messagehdlr, printdata, expr) );

   return SCIP_OKAY;
}

/** finishes printing of expressions in dot format */
SCIP_RETCODE SCIPprintExprDotFinal(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_EXPRPRINTDATA**    printdata         /**< buffer where dot printing data has been stored */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprPrintDotFinal(scip->set, scip->stat, scip->mem->probmem, printdata) );

   return SCIP_OKAY;
}

/** shows a single expression by use of dot and gv
 *
 * This function is meant for debugging purposes.
 * It's signature is kept as simple as possible to make it
 * easily callable from gdb, for example.
 *
 * It prints the expression into a temporary file in dot format, then calls dot to create a postscript file, then calls ghostview (gv) to show the file.
 * SCIP will hold until ghostscript is closed.
 */
SCIP_RETCODE SCIPshowExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression to be printed */
   )
{
   /* this function is for developers, so don't bother with C variants that don't have popen() */
#if _POSIX_C_SOURCE < 2
   SCIPerrorMessage("No POSIX version 2. Try http://distrowatch.com/.");
   return SCIP_ERROR;
#else
   SCIP_EXPRPRINTDATA* dotdata;
   FILE* f;

   assert(scip != NULL);
   assert(expr != NULL);

   /* call dot to generate postscript output and show it via ghostview */
   f = popen("dot -Tps | gv --media=a3 -", "w");
   if( f == NULL )
   {
      SCIPerrorMessage("Calling popen() failed");
      return SCIP_FILECREATEERROR;
   }

   /* print all of the expression into the pipe */
   SCIP_CALL( SCIPprintExprDotInit(scip, &dotdata, f, SCIP_EXPRPRINT_ALL) );
   SCIP_CALL( SCIPprintExprDot(scip, dotdata, expr) );
   SCIP_CALL( SCIPprintExprDotFinal(scip, &dotdata) );

   /* close the pipe */
   (void) pclose(f);

   return SCIP_OKAY;
#endif
}

/** prints structure of an expression a la Maple's dismantle */
SCIP_RETCODE SCIPdismantleExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< file to print to, or NULL for stdout */
   SCIP_EXPR*            expr                /**< expression to dismantle */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprDismantle(scip->set, scip->stat, scip->mem->blkmem, scip->messagehdlr, file, expr) );

   return SCIP_OKAY;
}

/**@} */

/**@name Expression Iterator Methods */
/**@{ */

/** creates an expression iterator */
SCIP_RETCODE SCIPcreateExpriter(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRITER**       iterator            /**< buffer to store expression iterator */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexpriterCreate(scip->set, scip->stat, scip->mem->probmem, iterator) );

   return SCIP_OKAY;
}

/** frees an expression iterator */
void SCIPfreeExpriter(
   SCIP_EXPRITER**       iterator            /**< pointer to the expression iterator */
   )
{
   SCIPexpriteratorFree(iterator);
}

/**@} */
