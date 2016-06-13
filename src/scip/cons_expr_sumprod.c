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

/**@file   cons_expr_sumprod.c
 * @brief  sum and product expression handlers
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 *
 * Implementation of the sum expression, representing a summation of a constant
 * and the arguments, each multiplied by a coefficients, i.e., sum_i a_i*x_i + constant.
 * Implementation of the product expression, representing a signomial term,
 * i.e., coef * prod_i x_i^e_i.
 * As both expressions store similar data, we implement them in the same C file.
 * The data (a_i and constant, or e_i and coef) is currently stored as a SCIP_Real
 * array of length nchildren + 1, storing the constant/coef in the first position,
 * and the a_i/e_i afterwards.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_sumprod.h"

#define SUM_PRECEDENCE      40000
#define PRODUCT_PRECEDENCE  50000

/** ensures that a block memory array has at least a given size
 *
 *  if cursize is 0, then *array1 can be NULL
 */
#define ENSUREBLOCKMEMORYARRAYSIZE(scip, array1, cursize, minsize)      \
   do {                                                                 \
      int __newsize;                                                    \
      assert((scip)  != NULL);                                          \
      if( (cursize) >= (minsize) )                                      \
         break;                                                         \
      __newsize = SCIPcalcMemGrowSize(scip, minsize);                   \
      assert(__newsize >= (minsize));                                   \
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(array1), cursize, __newsize) ); \
      (cursize) = __newsize;                                            \
   } while( FALSE )

/*
 * Data structures
 */

struct SCIP_ConsExpr_ExprData
{
   SCIP_Real  constant;     /**< constant coefficient */
   SCIP_Real* coefficients; /**< coefficients / exponents of children */
   int        ncoefs;       /**< number of coefficients (usually also number of children, but don't count on it) */
   int        coefssize;    /**< size of the coefficients array */
};

/*
 * Local methods
 */

static
SCIP_RETCODE createData(
   SCIP*                    scip,            /**< SCIP data structure */
   SCIP_CONSEXPR_EXPRDATA** exprdata,        /**< pointer where to store expression data */
   int                      ncoefficients,   /**< number of coefficients (i.e., number of children) */
   SCIP_Real*               coefficients,    /**< array with coefficients for all children (or NULL if all 1.0) */
   SCIP_Real                constant         /**< constant term of sum */
   )
{
   SCIP_Real* edata;

   assert(exprdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, exprdata) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &edata, ncoefficients) );

   if( coefficients != NULL )
   {
      memcpy(edata, coefficients, ncoefficients * sizeof(SCIP_Real));
   }
   else
   {
      int i;
      for( i = 0; i < ncoefficients; ++i )
         edata[i] = 1.0;
   }

   (*exprdata)->coefficients = edata;
   (*exprdata)->coefssize    = ncoefficients;
   (*exprdata)->constant     = constant;

   return SCIP_OKAY;
}

/* helper function to compare two sums or two products; see compareSum or compareProduct */
static
int compareSumProduct(
   SCIP_Real                const1,          /**< expr1's constant/factor */
   SCIP_Real*               coefs1,          /**< expr1's coefficients/exponents */
   SCIP_CONSEXPR_EXPR**     children1,       /**< expr1's children */
   int                      nchildren1,      /**< number of expr1's children */
   SCIP_Real                const2,          /**< expr2's constant/factor */
   SCIP_Real*               coefs2,          /**< expr2's coefficients/exponents */
   SCIP_CONSEXPR_EXPR**     children2,       /**< expr2's children */
   int                      nchildren2       /**< number of expr2's children */
   )
{
   int compareresult;
   int i;
   int j;

   for( i = nchildren1 - 1, j = nchildren2 - 1; i >= 0 && j >= 0; --i, --j )
   {
      compareresult = SCIPcompareConsExprExprs(children1[i], children2[j]);
      if( compareresult != 0 )
         return compareresult;
      else
      {
         /* expressions are equal, compare coefficient/exponent */
         if( coefs1[i] < coefs2[j] )
            return -1;
         if( coefs1[i] > coefs2[j] )
            return 1;

         /* coefficients are equal, continue */
      }
   }

   /* all children of one expression are children of the other expression, use number of children as a tie-breaker */
   if( i < j )
   {
      assert(i == -1);
      /* expr1 has less elements, hence expr1 < expr2 */
      return -1;
   }
   if( i > j )
   {
      assert(j == -1);
      /* expr1 has more elements, hence expr1 > expr2 */
      return 1;
   }

   /* everything is equal, use constant/coefficient as tie-breaker */
   assert(i == -1 && j == -1);
   if( const1 < const2 )
      return -1;
   if( const1 > const2 )
      return 1;

   /* they are equal */
   return 0;
}

/*
 * Callback methods of expression handler
 */

/** the order of two sum expressions is a lexicographical order on the terms.
 *  So, starting from the *last*, we find the first children where they differ.
 *  Suppose that is children is the i-th children, then u < v <=> u_i < v_i.
 *  If case there is no such children and they have different number of children, then u < v <=> nchildren(u) < nchildren(v)
 *  Otherwise, they are the same
 *  Note: we are assuming expression are simplified, so within u, we have u_1 < u_2, etc
 *  Example: y + z < x + y + z, 2*x + 3*y < 3*x + 3*y
 */
static
SCIP_DECL_CONSEXPR_EXPRCMP(compareSum)
{
   return compareSumProduct(
         SCIPgetConsExprExprSumConstant(expr1), SCIPgetConsExprExprSumCoefs(expr1),
         SCIPgetConsExprExprChildren(expr1), SCIPgetConsExprExprNChildren(expr1),
         SCIPgetConsExprExprSumConstant(expr2), SCIPgetConsExprExprSumCoefs(expr2),
         SCIPgetConsExprExprChildren(expr2), SCIPgetConsExprExprNChildren(expr2));
}

/** the order of two product expressions is a lexicographical order on the terms. See compareSum
 *  Example: y * z < x * y * z, x^2 * y^3 < x^3 * y^3
 */
static
SCIP_DECL_CONSEXPR_EXPRCMP(compareProduct)
{
   return compareSumProduct(
         SCIPgetConsExprExprProductCoef(expr1), SCIPgetConsExprExprProductExponents(expr1),
         SCIPgetConsExprExprChildren(expr1), SCIPgetConsExprExprNChildren(expr1),
         SCIPgetConsExprExprProductCoef(expr2), SCIPgetConsExprExprProductExponents(expr2),
         SCIPgetConsExprExprChildren(expr2), SCIPgetConsExprExprNChildren(expr2));
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrSum)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrSum(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrProduct)
{
   SCIP_CALL( SCIPincludeConsExprExprHdlrProduct(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataSumProduct)
{
   SCIP_CONSEXPR_EXPRDATA* sourceexprdata;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   sourceexprdata = SCIPgetConsExprExprData(sourceexpr);
   assert(sourceexprdata != NULL);

   SCIP_CALL( createData(targetscip, targetexprdata, SCIPgetConsExprExprNChildren(sourceexpr),
            sourceexprdata->coefficients, sourceexprdata->constant) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataSumProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPfreeBlockMemoryArray(scip, &(exprdata->coefficients), exprdata->coefssize);
   SCIPfreeBlockMemory(scip, &exprdata);

   SCIPsetConsExprExprData(expr, NULL);

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRPRINT(printSum)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print opening parenthesis, if necessary */
         if( SUM_PRECEDENCE <= SCIPgetConsExprExprWalkParentPrecedence(expr) )
         {
            SCIPinfoMessage(scip, file, "(");
         }

         /* print constant, if nonzero */
         if( exprdata->constant != 0.0 )
         {
            SCIPinfoMessage(scip, file, "%g", exprdata->constant);
         }
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         int childidx;
         SCIP_Real coef;

         childidx = SCIPgetConsExprExprWalkCurrentChild(expr);
         coef = exprdata->coefficients[childidx];

         /* print coefficient, if necessary */
         if( coef == 1.0 )
         {
            /* if coefficient is 1.0, then print only "+" if not the first term */
            if( exprdata->constant != 0.0 || childidx > 0 )
            {
               SCIPinfoMessage(scip, file, "+");
            }
         }
         else if( coef == -1.0 )
         {
            /* if coefficient is -1.0, then print only "-" */
            SCIPinfoMessage(scip, file, "-");
         }
         else
         {
            /* force "+" sign on positive coefficient if not the first term */
            SCIPinfoMessage(scip, file, (exprdata->constant != 0.0 || childidx > 0) ? "%+g*" : "%g*", coef);
         }

         break;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      {
         /* print closing parenthesis, if necessary */
         if( SUM_PRECEDENCE <= SCIPgetConsExprExprWalkParentPrecedence(expr) )
         {
            SCIPinfoMessage(scip, file, ")");
         }
         break;
      }

      default: ;
   }

   return SCIP_OKAY;
}


static
SCIP_DECL_CONSEXPR_EXPRPRINT(printProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* print opening parenthesis, if necessary */
         if( PRODUCT_PRECEDENCE <= SCIPgetConsExprExprWalkParentPrecedence(expr) )
         {
            SCIPinfoMessage(scip, file, "(");
         }

         /* print constant factor, if not one */
         if( exprdata->constant != 1.0 )
         {
            if( exprdata->constant < 0.0 && PRODUCT_PRECEDENCE > SCIPgetConsExprExprWalkParentPrecedence(expr) )
            {
               SCIPinfoMessage(scip, file, "(%g)", exprdata->constant);
            }
            else
            {
               SCIPinfoMessage(scip, file, "%g", exprdata->constant);
            }
         }
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         int childidx = SCIPgetConsExprExprWalkCurrentChild(expr);

         if( exprdata->coefficients[childidx] >= 0.0 )
         {
            /* print multiplication sign, if not first factor */
            if( exprdata->constant != 1.0 || childidx > 0 )
            {
               SCIPinfoMessage(scip, file, "*");
            }
         }
         else
         {
            if( exprdata->constant != 1.0 || childidx > 0 )
            {
               /* print division sign, if not first factor */
               SCIPinfoMessage(scip, file, "/");
            }
            else
            {
               /* print "1/", if first factor */
               SCIPinfoMessage(scip, file, "1/");
            }
         }
         break;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      {
         SCIP_Real exponent;
         exponent = exprdata->coefficients[SCIPgetConsExprExprWalkCurrentChild(expr)];

         /* print absolute value of exponent, if not 1.0 (sign is taken care of in VISITINGCHILD) */
         exponent = REALABS(exponent);
         if( exponent != 1.0 )
         {
            SCIPinfoMessage(scip, file, "^%g", exponent);
         }

         break;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      {
         /* print closing parenthesis, if necessary */
         if( PRODUCT_PRECEDENCE <= SCIPgetConsExprExprWalkParentPrecedence(expr) )
         {
            SCIPinfoMessage(scip, file, ")");
         }
         break;
      }
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPREVAL(evalSum)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *val = exprdata->constant;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[c]) != SCIP_INVALID);

      *val += exprdata->coefficients[c] * SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[c]);
   }

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalSum)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_INTERVAL suminterval;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPintervalSet(interval, exprdata->constant);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_INTERVAL childinterval;

      childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[c]);
      assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

      /* compute coefficients[c] * childinterval and add the result to the so far computed interval */
      if( exprdata->coefficients[c] == 1.0 )
      {
         SCIPintervalAdd(SCIPinfinity(scip), interval, *interval, childinterval);
      }
      else
      {
         SCIPintervalMulScalar(SCIPinfinity(scip), &suminterval, childinterval, exprdata->coefficients[c]);
         SCIPintervalAdd(SCIPinfinity(scip), interval, *interval, suminterval);
      }
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPREVAL(evalProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_Real childval;
   SCIP_Real powval;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *val = exprdata->constant;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      childval = SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[c]);
      assert(childval != SCIP_INVALID);

      /* according to the man page of pow(), this should also work fine for cases like pow(<negative>, <integer>) */
      powval = pow(childval, exprdata->coefficients[c]);

      /* if there is a domain, pole, or range error, pow() should return some kind of NaN, infinity, or HUGE_VAL
       * we could also work with floating point exceptions or errno, but I am not sure this would be thread-safe
       */
      if( !SCIPisFinite(powval) || powval == HUGE_VAL || powval == -HUGE_VAL )
      {
         *val = SCIP_INVALID;
         return SCIP_OKAY;
      }

      *val *= powval;
   }

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_INTERVAL powinterval;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPintervalSet(interval, exprdata->constant);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_INTERVAL childinterval;

      childinterval = SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[c]);
      assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), childinterval));

      /* compute interval resulting from childinterval^exprdata->coefficients[c] */
      SCIPintervalPowerScalar(SCIPinfinity(scip), &powinterval, childinterval, exprdata->coefficients[c]);

      if( SCIPintervalIsEmpty(SCIPinfinity(scip), powinterval) )
      {
         SCIPintervalSetEmpty(interval);
         return SCIP_OKAY;
      }

      /* multiply powinterval with the so far computed interval */
      SCIPintervalMul(SCIPinfinity(scip), interval, *interval, powinterval);
   }

   return SCIP_OKAY;
}

/** creates the handler for sum expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "sum", "summation with coefficients and a constant",
         SUM_PRECEDENCE, evalSum, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrSum, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataSumProduct, freedataSumProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalSum) );

   return SCIP_OKAY;
}

/** creates a sum expression */
SCIP_RETCODE SCIPcreateConsExprExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children,           /**< children */
   SCIP_Real*            coefficients,       /**< array with coefficients for all children (or NULL if all 1.0) */
   SCIP_Real             constant            /**< constant term of sum */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   SCIP_CALL( createData(scip, &exprdata, nchildren, coefficients, constant) );

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPgetConsExprExprHdlrSum(consexprhdlr), exprdata, nchildren, children) );

   return SCIP_OKAY;
}

/** gets the coefficients of a summation expression */
SCIP_Real* SCIPgetConsExprExprSumCoefs(
   SCIP_CONSEXPR_EXPR*   expr                /**< sum expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   return exprdata->coefficients;
}

/** gets the constant of a summation expression */
SCIP_Real SCIPgetConsExprExprSumConstant(
   SCIP_CONSEXPR_EXPR*   expr                /**< sum expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   return exprdata->constant;
}

/** sets the constant of a summation expression */
void SCIPsetConsExprExprSumConstant(
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   exprdata->constant = constant;
}



/** creates the handler for product expressions and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr        /**< expression constraint handler */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, "prod",
         "product of children with exponents (actually a signomial)", PRODUCT_PRECEDENCE, evalProduct, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrProduct, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataSumProduct, freedataSumProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalProduct) );

   return SCIP_OKAY;
}

/** creates a product expression */
SCIP_RETCODE SCIPcreateConsExprExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children,           /**< children */
   SCIP_Real*            exponents,          /**< array with exponents for all children (or NULL if all 1.0) */
   SCIP_Real             constant            /**< constant coefficient of product */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   SCIP_CALL( createData(scip, &exprdata, nchildren, exponents, constant) );

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPgetConsExprExprHdlrProduct(consexprhdlr), exprdata, nchildren, children) );

   return SCIP_OKAY;
}

/** gets the exponents of a product expression */
SCIP_Real* SCIPgetConsExprExprProductExponents(
   SCIP_CONSEXPR_EXPR*   expr                /**< product expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   return exprdata->coefficients;
}

/** gets the constant coefficient of a product expression */
SCIP_Real SCIPgetConsExprExprProductCoef(
   SCIP_CONSEXPR_EXPR*   expr                /**< product expression */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   return exprdata->constant;
}

/** appends an expression to a sum expression */
SCIP_RETCODE SCIPappendConsExprExprSumExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_CONSEXPR_EXPR*   child,              /**< expression to be appended */
   SCIP_Real             childcoef           /**< child's coefficient */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int nchildren;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   nchildren = SCIPgetConsExprExprNChildren(expr);

   ENSUREBLOCKMEMORYARRAYSIZE(scip, exprdata->coefficients, exprdata->coefssize, nchildren + 1);

   assert(exprdata->coefssize > nchildren);
   exprdata->coefficients[nchildren] = childcoef;

   SCIP_CALL( SCIPappendConsExprExpr(scip, expr, child) );

   return SCIP_OKAY;
}

/** appends an expression to a product expression */
SCIP_RETCODE SCIPappendConsExprExprProductExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< product expression */
   SCIP_CONSEXPR_EXPR*   child,              /**< expression to be appended */
   SCIP_Real             childexponent       /**< child's exponent */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int nchildren;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   nchildren = SCIPgetConsExprExprNChildren(expr);

   ENSUREBLOCKMEMORYARRAYSIZE(scip, exprdata->coefficients, exprdata->coefssize, nchildren + 1);

   assert(exprdata->coefssize > nchildren);
   exprdata->coefficients[nchildren] = childexponent;

   SCIP_CALL( SCIPappendConsExprExpr(scip, expr, child) );

   return SCIP_OKAY;
}

/** multiplies given sum expr by a constant */
void SCIPmultiplyConsExprExprSumByConstant(
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant that multiplies sum expression */
   )
{
   int i;
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
   {
      exprdata->coefficients[i] *= constant;
   }
   exprdata->constant *= constant;
}

/** exponentiate given product expr by a constant
 * TODO: should this function create abs children when exponent is fractional and resulting exponent is odd?
 * FIXME: No! Somebody (who?) should have enforced (add a constraint?, change bound, etc) that the child
 * with fractional exponent is positive */
void SCIPexponentiateConsExprExprProductByConstant(
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_Real             exponent            /**< exponent */
   )
{
   int i;
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
   {
      exprdata->coefficients[i] *= exponent;
   }
   exprdata->constant = pow(exprdata->constant, exponent);
}
