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

/**@file   cons_expr_sum.c
 * @brief  sum expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 *
 * Implementation of the sum expression, representing a summation of a constant
 * and the arguments, each multiplied by a coefficients, i.e., sum_i a_i*x_i + constant.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include <stddef.h>

#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_product.h"
#include "scip/cons_expr_exp.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_rowprep.h"

#define EXPRHDLR_NAME         "sum"
#define EXPRHDLR_DESC         "summation with coefficients and a constant"
#define EXPRHDLR_PRECEDENCE      40000
#define EXPRHDLR_HASHKEY         SCIPcalcFibHash(47161.0)

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

/** macro to activate/deactivate debugging information of simplify method */
#ifdef SIMPLIFY_DEBUG
#define debugSimplify                   printf
#else
#define debugSimplify                   while( FALSE ) printf
#endif

/*
 * Data structures
 */

struct SCIP_ConsExpr_ExprData
{
   SCIP_Real  constant;     /**< constant coefficient */
   SCIP_Real* coefficients; /**< coefficients of children */
   int        coefssize;    /**< size of the coefficients array */
};

/*
 * Local methods
 */

/** simplifies the `idx`-th child of the sum expression `duplicate` in order for it to be able to be a child of a
 * simplified sum; for example, this means that the `idx`-th child cannot be itself a sum; if it is, we have to flatten
 * it, i.e., take all its children and make them children of `duplicate`
 */
static
SCIP_RETCODE simplifyTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< consexpr handler */
   SCIP_CONSEXPR_EXPR*   duplicate,          /**< expression to be simplified */
   int                   idx,                /**< idx of children to be simplified */
   SCIP_Bool*            changed             /**< pointer to store if some term actually got simplified */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_CONSEXPR_EXPR** children;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Real* coefs;
   SCIP_Real constant ;
   SCIP_Real coef;

   children  = SCIPgetConsExprExprChildren(duplicate);
   coefs     = SCIPgetConsExprExprSumCoefs(duplicate);
   constant  = SCIPgetConsExprExprSumConstant(duplicate);

   assert(0 <= idx && idx < SCIPgetConsExprExprNChildren(duplicate));

   coef = coefs[idx];
   expr = children[idx];

   assert(expr != NULL);

   exprhdlr = SCIPgetConsExprExprHdlr(expr);

   /* enforces SS3 */
   if( exprhdlr == SCIPgetConsExprExprHdlrValue(conshdlr) )
   {
      *changed = TRUE;
      constant += coef * SCIPgetConsExprExprValueValue(expr);
      SCIPsetConsExprExprSumConstant(duplicate, constant);

      /* TODO: remove child? */
      coefs[idx] = 0.0;

      return SCIP_OKAY;
   }

   /* enforces SS2 */
   if( exprhdlr == SCIPgetConsExprExprHdlrSum(conshdlr) )
   {
      *changed = TRUE;

      /* pass constant to parent */
      constant += coef * SCIPgetConsExprExprSumConstant(expr);
      SCIPsetConsExprExprSumConstant(duplicate, constant);

      /* append all children of expr on parent except the first one */
      if( SCIPgetConsExprExprNChildren(expr) > 1 )
      {
         int i;

         for( i = 1; i < SCIPgetConsExprExprNChildren(expr); ++i )
         {
            assert(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(expr)[i]) != SCIPgetConsExprExprHdlrSum(conshdlr));
            SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, duplicate, SCIPgetConsExprExprChildren(expr)[i], coef *
                     SCIPgetConsExprExprSumCoefs(expr)[i]) );
         }
      }

      /* replace expr with first child; need to get data again since it might be re-allocated */
      assert(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(expr)[0]) != SCIPgetConsExprExprHdlrSum(conshdlr));

      coefs = SCIPgetConsExprExprSumCoefs(duplicate);

      coefs[idx] = coef * SCIPgetConsExprExprSumCoefs(expr)[0];
      SCIP_CALL( SCIPreplaceConsExprExprChild(scip, duplicate, idx, SCIPgetConsExprExprChildren(expr)[0]) );

      return SCIP_OKAY;
   }

   /* enforce SS9 */
   if( REALABS(coef) != 1.0 && exprhdlr == SCIPgetConsExprExprHdlrProduct(conshdlr) )
   {
      SCIP_CONSEXPR_EXPR* expchild = NULL;
      int i;

      for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
      {
         SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[i];
         assert(child != NULL);

         if( SCIPgetConsExprExprHdlr(child) == SCIPgetConsExprExprHdlrExponential(conshdlr) )
         {
            expchild = child;
            break;
         }
      }

      /* coef != +- 1, term is product and one factor is an exponential -> enforce SS9 */
      if( expchild != NULL )
      {
         SCIP_CONSEXPR_EXPR* sum;
         SCIP_CONSEXPR_EXPR* prod;
         SCIP_CONSEXPR_EXPR* simplifiedprod;
         SCIP_CONSEXPR_EXPR* simplifiedsum;
         SCIP_CONSEXPR_EXPR* exponential;
         SCIP_CONSEXPR_EXPR* simplifiedexp;
         SCIP_Real expconstant;

         /* inform that expression will change */
         *changed = TRUE;

         /* compute expchild's coefficient as +- 1.0 * exp(log(abs(coef))) */
         if( coef > 0.0 )
         {
            expconstant = log(coef);
            coefs[idx] = 1.0;
         }
         else
         {
            expconstant = log(-coef);
            coefs[idx] = -1.0;
         }

         /* add constant to exponential's child */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, &sum, 1, SCIPgetConsExprExprChildren(expchild), NULL,
                  expconstant) );

         /* simplify sum */
         SCIP_CALL( SCIPsimplifyConsExprExprHdlr(scip, conshdlr, sum, &simplifiedsum) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &sum) );

         /* create exponential with new child */
         SCIP_CALL( SCIPcreateConsExprExprExp(scip, conshdlr, &exponential, simplifiedsum) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplifiedsum) );

         /* simplify exponential */
         SCIP_CALL( SCIPsimplifyConsExprExprHdlr(scip, conshdlr, exponential, &simplifiedexp) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &exponential) );

         /* create product with new child */
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &prod, 0, NULL, 1.0) );

         for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
         {
            if( SCIPgetConsExprExprChildren(expr)[i] == expchild )
            {
               SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, prod, simplifiedexp) );
            }
            else
            {
               SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, prod, SCIPgetConsExprExprChildren(expr)[i]) );
            }
         }
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplifiedexp) );

         /* simplify product */
         SCIP_CALL( SCIPsimplifyConsExprExprHdlr(scip, conshdlr, prod, &simplifiedprod) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prod) );

         /* replace current child with simplified product */
         SCIP_CALL( SCIPreplaceConsExprExprChild(scip, duplicate, idx, simplifiedprod) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplifiedprod) );

         /* since the simplified product can be a sum ( exp(-1)*exp(log(x+y)+1) -> x+y ), we call the function we are in
          * again; this is not recursive, since the coef is now +- 1
          */
         SCIP_CALL( simplifyTerm(scip, conshdlr, duplicate, idx, changed) );

         return SCIP_OKAY;
      }
   }

   /* other types of (simplified) expressions can be a child of a simplified sum */
   assert(exprhdlr != SCIPgetConsExprExprHdlrSum(conshdlr));
   assert(exprhdlr != SCIPgetConsExprExprHdlrValue(conshdlr));

   return SCIP_OKAY;
}


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

/*
 * Callback methods of expression handler
 */

static
SCIP_DECL_SORTPTRCOMP(sortExprComp)
{
   SCIP_CONSEXPR_EXPR* expr1 = (SCIP_CONSEXPR_EXPR*) elem1;
   SCIP_CONSEXPR_EXPR* expr2 = (SCIP_CONSEXPR_EXPR*) elem2;

   return SCIPcompareConsExprExprs(expr1, expr2);
}

/** simplifies a sum expression
 *
 * goes through each child and simplifies it; then sorts the simplified children; then sum the children that are equal;
 * finally creates a sum expression with all the children that do not have a 0 coefficient and post-process so that SS6
 * and SS7 are satisfied
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifySum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR** children;
   SCIP_CONSEXPR_EXPR* duplicate = NULL;
   SCIP_CONSEXPR_EXPR** newchildren = NULL;
   SCIP_Real* newcoefs = NULL;
   int nnewchildren;
   SCIP_Real newconstant;

   SCIP_Real* coefs;
   int i;
   int nchildren;
   SCIP_Bool changed;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPgetConsExprExprHdlr(expr) == SCIPgetConsExprExprHdlrSum(conshdlr));

   changed = FALSE;

   /* TODO: maybe have a flag to know if it is simplified ? */
   SCIP_CALL( SCIPduplicateConsExprExpr(scip, conshdlr, expr, &duplicate, TRUE) );
   assert(duplicate != NULL);

   nchildren = SCIPgetConsExprExprNChildren(duplicate);

   for( i = 0; i < nchildren; i++ )
   {
      /* enforces SS8 TODO: remove child? */
      /* we have to ask for the coefs everytime, since it might get realloced in simpifyTerm */
      if( SCIPgetConsExprExprSumCoefs(duplicate)[i] == 0.0 )
      {
         changed = TRUE;
         continue;
      }

      /* enforces SS2, SS3 and SS9 */
      SCIP_CALL( simplifyTerm(scip, conshdlr, duplicate, i, &changed) );
   }

   /* simplifyTerm can add new children to duplicate and realloc them; so get them again */
   nchildren = SCIPgetConsExprExprNChildren(duplicate);

   /* enforces SS5: sort children; if nothing has changed so far, we need to find to find out if sorting changes
    * anything
    */
   if( changed )
   {
      SCIPsortPtrPtr((void**)SCIPgetConsExprExprChildren(duplicate), (void**)SCIPgetConsExprExprSumCoefs(duplicate),
            sortExprComp, nchildren);
   }
   else
   {
      int* order;

      SCIP_CALL( SCIPallocBufferArray(scip, &order, nchildren) );
      for( i = 0; i < nchildren; i++ )
         order[i] = i;

      SCIPsortPtrPtrInt((void**)SCIPgetConsExprExprChildren(duplicate), (void**)SCIPgetConsExprExprSumCoefs(duplicate),
            order, sortExprComp, nchildren);

      for( i = 0; i < nchildren; i++ )
      {
         if( order[i] != i )
         {
            changed = TRUE;
            break;
         }
      }
      SCIPfreeBufferArray(scip, &order);
   }

   /* post-process */
   children = SCIPgetConsExprExprChildren(duplicate);
   coefs    = SCIPgetConsExprExprSumCoefs(duplicate);

   /* treat zero term case */
   if( nchildren == 0 )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, SCIPgetConsExprExprSumConstant(duplicate)) );
      goto CLEANUP;
   }

   /* treat one term case */
   if( nchildren == 1 )
   {
      if( coefs[0] == 0.0 )
      {
         SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, SCIPgetConsExprExprSumConstant(duplicate)) );
         goto CLEANUP;
      }

      if( coefs[0] == 1.0 && SCIPgetConsExprExprSumConstant(duplicate) == 0.0 )
         *simplifiedexpr = children[0]; /* SS7 */
      else
         *simplifiedexpr = changed ? duplicate : expr;

      SCIPcaptureConsExprExpr(*simplifiedexpr);

      goto CLEANUP;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &newchildren, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newcoefs, nchildren) );

   /* enforces SS4 */
   nnewchildren = 0;
   for( i = 0; i < nchildren - 1; i++ )
   {
      if( coefs[i] == 0.0 )
      {
         changed = TRUE;
         continue;
      }

      /* sum the expressions */
      if( SCIPcompareConsExprExprs(children[i], children[i + 1]) == 0 )
      {
         changed = TRUE;
         /* if we substract two almost equal not-so-small numbers, then set new coefficient to 0.0
          * instead of some tiny value that is likely the result of some random round-off error
          * E.g., on instance ex1221, we have x1^2 + b3 = 1.25.
          *   Probing finds an aggregation x1 = 1.11803 - 0.618034 b3.
          *   Simplify would then produce 1.25 + 1e-16 x1 = 1.25.
          */
         if( SCIPisEQ(scip, coefs[i], -coefs[i+1]) && REALABS(coefs[i]) >= 1.0 )
            coefs[i+1] = 0.0;
         else
            coefs[i+1] += coefs[i];
         continue;
      }

      /* insert expression to newchildren */
      newchildren[nnewchildren] = children[i];
      newcoefs[nnewchildren] = coefs[i];
      nnewchildren++;
   }
   /* treat last one */
   assert(i == nchildren - 1);
   if( coefs[i] == 0.0 )
   {
      changed = TRUE;
   }
   else
   {
      /* insert expression to newchildren */
      newchildren[nnewchildren] = children[i];
      newcoefs[nnewchildren] = coefs[i];
      nnewchildren++;
   }

   /* build sum expression from finalchildren and post-simplify */
   newconstant = SCIPgetConsExprExprSumConstant(duplicate);

   debugSimplify("what to do? finalchildren has length %d\n", nnewchildren); /*lint !e506 !e681*/

   /* enforces SS6: if they are no children, return value */
   if( nnewchildren == 0 )
   {
      debugSimplify("[sum] got empty list, return value %g\n", newconstant); /*lint !e506 !e681*/
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, simplifiedexpr, newconstant) );

      goto CLEANUP;
   }

   /* enforces SS7
    * if list consists of one expr with coef 1.0 and constant is 0, return that expr */
   if( nnewchildren == 1 && newcoefs[0] == 1.0 && newconstant == 0.0 )
   {
      *simplifiedexpr = newchildren[0];
      SCIPcaptureConsExprExpr(*simplifiedexpr);

      goto CLEANUP;
   }

   /* build sum expression from children */
   if( changed )
   {
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, simplifiedexpr, nnewchildren, newchildren, newcoefs,
               newconstant) );

      goto CLEANUP;
   }

   *simplifiedexpr = expr;

   /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
   SCIPcaptureConsExprExpr(*simplifiedexpr);

   /* free memory */
CLEANUP:
   SCIPfreeBufferArrayNull(scip, &newcoefs);
   SCIPfreeBufferArrayNull(scip, &newchildren);
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &duplicate) );

   return SCIP_OKAY;
}

/** the order of two sum expressions is a lexicographical order on the terms.
 *  Starting from the *last*, we find the first child where they differ, say, the i-th.
 *  Then u < v <=> u_i < v_i.
 *  If there is no such children and they have different number of children, then u < v <=> nchildren(u) < nchildren(v)
 *  If there is no such children and they have the same number of children, then u < v <=> const(u) < const(v)
 *  Otherwise, they are the same
 *  Note: we are assuming expression are simplified, so within u, we have u_1 < u_2, etc
 *  Example: y + z < x + y + z, 2*x + 3*y < 3*x + 3*y
 */
static
SCIP_DECL_CONSEXPR_EXPRCOMPARE(compareSum)
{  /*lint --e{715}*/
   SCIP_Real const1;
   SCIP_Real* coefs1;
   SCIP_CONSEXPR_EXPR** children1;
   int nchildren1;
   SCIP_Real const2;
   SCIP_Real* coefs2;
   SCIP_CONSEXPR_EXPR** children2;
   int nchildren2;
   int compareresult;
   int i;
   int j;

   nchildren1 = SCIPgetConsExprExprNChildren(expr1);
   nchildren2 = SCIPgetConsExprExprNChildren(expr2);
   children1 = SCIPgetConsExprExprChildren(expr1);
   children2 = SCIPgetConsExprExprChildren(expr2);
   coefs1 = SCIPgetConsExprExprSumCoefs(expr1);
   coefs2 = SCIPgetConsExprExprSumCoefs(expr2);
   const1 = SCIPgetConsExprExprSumConstant(expr1);
   const2 = SCIPgetConsExprExprSumConstant(expr2);

   for( i = nchildren1 - 1, j = nchildren2 - 1; i >= 0 && j >= 0; --i, --j )
   {
      compareresult = SCIPcompareConsExprExprs(children1[i], children2[j]);
      if( compareresult != 0 )
         return compareresult;
      else
      {
         /* expressions are equal, compare coefficient */
         if( (coefs1 ? coefs1[i] : 1.0) < (coefs2 ? coefs2[j] : 1.0) )
            return -1;
         if( (coefs1 ? coefs1[i] : 1.0) > (coefs2 ? coefs2[j] : 1.0) )
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

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrSum)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrSum(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataSum)
{  /*lint --e{715}*/
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
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataSum)
{  /*lint --e{715}*/
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
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   switch( stage )
   {
      case SCIP_CONSEXPRITERATOR_ENTEREXPR :
      {
         /* print opening parenthesis, if necessary */
         if( EXPRHDLR_PRECEDENCE <= parentprecedence )
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

      case SCIP_CONSEXPRITERATOR_VISITINGCHILD :
      {
         SCIP_Real coef;

         coef = exprdata->coefficients[currentchild];

         /* print coefficient, if necessary */
         if( coef == 1.0 )
         {
            /* if coefficient is 1.0, then print only "+" if not the first term */
            if( exprdata->constant != 0.0 || currentchild > 0 )
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
            SCIPinfoMessage(scip, file, (exprdata->constant != 0.0 || currentchild > 0) ? "%+g*" : "%g*", coef);
         }

         break;
      }

      case SCIP_CONSEXPRITERATOR_LEAVEEXPR :
      {
         /* print closing parenthesis, if necessary */
         if( EXPRHDLR_PRECEDENCE <= parentprecedence )
         {
            SCIPinfoMessage(scip, file, ")");
         }
         break;
      }

      case SCIP_CONSEXPRITERATOR_VISITEDCHILD :
      default: ;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPREVAL(evalSum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *val = exprdata->constant;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      assert(SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[c]) != SCIP_INVALID); /*lint !e777*/

      *val += exprdata->coefficients[c] * SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[c]);
   }

   return SCIP_OKAY;
}

/** expression forward derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRFWDIFF(fwdiffSum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int c;

   assert(expr != NULL);
   assert(dot != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *dot = 0.0;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      assert(SCIPgetConsExprExprDot(SCIPgetConsExprExprChildren(expr)[c]) != SCIP_INVALID); /*lint !e777*/

      *dot += exprdata->coefficients[c] * SCIPgetConsExprExprDot(SCIPgetConsExprExprChildren(expr)[c]);
   }

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffSum)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);
   assert(childidx >= 0 && childidx < SCIPgetConsExprExprNChildren(expr));
   assert(SCIPgetConsExprExprChildren(expr)[childidx] != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(expr)[childidx])), "val") != 0);

   *val = SCIPgetConsExprExprSumCoefs(expr)[childidx];

   return SCIP_OKAY;
}

/** expression backward forward derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWFWDIFF(bwfwdiffSum)
{  /*lint --e{715}*/
   assert(bardot != NULL);

   *bardot = 0.0;

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEVAL(intevalSum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_INTERVAL suminterval;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIPintervalSet(interval, exprdata->constant);

   SCIPdebugMsg(scip, "inteval %p with %d children: %.20g", (void*)expr, SCIPgetConsExprExprNChildren(expr), exprdata->constant);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_INTERVAL childinterval;

      childinterval = SCIPgetConsExprExprActivity(scip, SCIPgetConsExprExprChildren(expr)[c]);
      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childinterval) )
      {
         SCIPintervalSetEmpty(interval);
         break;
      }

      /* compute coefficients[c] * childinterval and add the result to the so far computed interval */
      if( exprdata->coefficients[c] == 1.0 )
      {
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, interval, *interval, childinterval);
      }
      else
      {
         SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &suminterval, childinterval, exprdata->coefficients[c]);
         SCIPintervalAdd(SCIP_INTERVAL_INFINITY, interval, *interval, suminterval);
      }

      SCIPdebugMsgPrint(scip, " %+.20g*[%.20g,%.20g]", exprdata->coefficients[c], childinterval.inf, childinterval.sup);

   }
   SCIPdebugMsgPrint(scip, " = [%.20g,%.20g]\n", interval->inf, interval->sup);

   return SCIP_OKAY;
}

/** helper function to separate a given point; needed for proper unittest */
static
SCIP_RETCODE separatePointSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_Bool             overestimate,       /**< should the expression be overestimated? */
   SCIP_ROWPREP**        rowprep             /**< pointer to store the row preparation */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_VAR* auxvar;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(rowprep != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   auxvar = SCIPgetConsExprExprAuxVar(expr);
   assert(auxvar != NULL);

   *rowprep = NULL;

   /* create rowprep */
   SCIP_CALL( SCIPcreateRowprep(scip, rowprep, overestimate ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT, FALSE));
   SCIP_CALL( SCIPensureRowprepSize(scip, *rowprep, SCIPgetConsExprExprNChildren(expr) + 1) );

   /* compute  w = sum_i alpha_i z_i + const */
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIP_CONSEXPR_EXPR* child;
      SCIP_VAR* var;

      child = SCIPgetConsExprExprChildren(expr)[c];
      assert(child != NULL);

      /* value expressions should have been removed during simplification */
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "value") != 0);

      var = SCIPgetConsExprExprAuxVar(child);
      assert(var != NULL);

      SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, var, exprdata->coefficients[c]) );
   }

   /* add -1 * auxvar and set side */
   SCIP_CALL( SCIPaddRowprepTerm(scip, *rowprep, auxvar, -1.0) );
   SCIPaddRowprepSide(*rowprep, -exprdata->constant);

   (void) SCIPsnprintf((*rowprep)->name, SCIP_MAXSTRLEN, "sum");  /* @todo make cutname unique, e.g., add LP number */

   return SCIP_OKAY;
}

/* TODO we should get rid of this somehow
 * maybe initsepa of nlhdlr_default should call estimate,
 * but it currently doesn't know for which expr's the estimator can be computed without a solution
 */
/** separation initialization callback */
static
SCIP_DECL_CONSEXPR_EXPRINITSEPA(initSepaSum)
{  /*lint --e{715}*/
   int i;

   assert(overestimate || underestimate);

   *infeasible = FALSE;

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "initSepaSum %d children: ", SCIPgetConsExprExprNChildren(expr));
   SCIPprintConsExprExpr(scip, conshdlr, expr, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   /* i = 0 for overestimation; i = 1 for underestimation */
   for( i = 0; i < 2 && !*infeasible; ++i )
   {
      SCIP_ROWPREP* rowprep;
      SCIP_Bool success;

      if( (i == 0 && !overestimate) || (i == 1 && !underestimate) )
         continue;

      /* create rowprep */
      SCIP_CALL( separatePointSum(scip, expr, i == 0, &rowprep) );
      assert(rowprep != NULL);

      /* first try scale-up rowprep to try to get rid of within-epsilon of integer coefficients */
      (void) SCIPscaleupRowprep(scip, rowprep, 1.0, &success);

      if( success && underestimate && overestimate )
      {
         SCIP_ROW* row;

         assert(i == 0);

         SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

         /* since we did not relax the overestimator (i=0), we can turn the row into an equality if we need an underestimator, too */
         if( rowprep->sidetype == SCIP_SIDETYPE_LEFT )
         {
            SCIP_CALL( SCIPchgRowRhs(scip, row, rowprep->side) );
         }
         else
         {
            SCIP_CALL( SCIPchgRowLhs(scip, row, rowprep->side) );
         }

#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "adding row ");
         SCIPprintRow(scip, row, NULL);
#endif

         SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );

         /* free rowprep */
         SCIPfreeRowprep(scip, &rowprep);

         break;
      }

      if( !success )
      {
         /* if scale-up is not sufficient, then do clean-up
          * this might relax the row, so we only get a bounding cut
          */
         SCIP_CALL( SCIPcleanupRowprep(scip, rowprep, NULL, SCIP_CONSEXPR_CUTMAXRANGE, 0.0, NULL, &success) );
      }

      /* create a SCIP_ROW and add it to the initial LP */
      if( success )
      {
         SCIP_ROW* row;

         SCIP_CALL( SCIPgetRowprepRowCons(scip, &row, rowprep, cons) );

#ifdef SCIP_DEBUG
         SCIPinfoMessage(scip, NULL, "adding row ");
         SCIPprintRow(scip, row, NULL);
         SCIPinfoMessage(scip, NULL, "\n");
#endif

         SCIP_CALL( SCIPaddRow(scip, row, FALSE, infeasible) );
         SCIP_CALL( SCIPreleaseRow(scip, &row) );
      }

      /* free rowprep */
      SCIPfreeRowprep(scip, &rowprep);
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRESTIMATE(estimateSum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(success != NULL);
   assert(branchcand != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   memcpy(coefs, exprdata->coefficients, SCIPgetConsExprExprNChildren(expr) * sizeof(SCIP_Real));
   *constant = exprdata->constant;
   *success = TRUE;

   /* for none of our children, branching would improve the underestimator, so set branchcand[i]=FALSE everywhere
    * if we branch for numerical reasons, then cons-expr-core should figure out what the candidates are
    */
   memset(branchcand, 0, SCIPgetConsExprExprNChildren(expr) * sizeof(SCIP_Bool));

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_EXPRREVERSEPROP(reversepropSum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) > 0);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   SCIP_CALL( SCIPreverseConsExprExprPropagateWeightedSum(scip, conshdlr, SCIPgetConsExprExprNChildren(expr),
            SCIPgetConsExprExprChildren(expr), exprdata->coefficients, exprdata->constant,
            bounds, infeasible, nreductions) );

   return SCIP_OKAY;
}

/** sum hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashSum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(exprdata->constant);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
      *hashkey ^= SCIPcalcFibHash(exprdata->coefficients[c]) ^ childrenhashes[c];

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_CONSEXPR_EXPRCURVATURE(curvatureSum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(childcurv != NULL);
   assert(success != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
      childcurv[i] = SCIPexprcurvMultiply(exprdata->coefficients[i], exprcurvature);

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_CONSEXPR_EXPRMONOTONICITY(monotonicitySum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx >= 0 && childidx < SCIPgetConsExprExprNChildren(expr));

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *result = exprdata->coefficients[childidx] >= 0.0 ? SCIP_MONOTONE_INC : SCIP_MONOTONE_DEC;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_CONSEXPR_EXPRINTEGRALITY(integralitySum)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *isintegral = EPSISINT(exprdata->constant, 0.0); /*lint !e835*/

   for( i = 0; i < SCIPgetConsExprExprNChildren(expr) && *isintegral; ++i )
   {
      SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[i];
      assert(child != NULL);

      *isintegral = EPSISINT(exprdata->coefficients[i], 0.0) && SCIPisConsExprExprIntegral(child); /*lint !e835*/
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

   SCIP_CALL( SCIPincludeConsExprExprHdlrBasic(scip, consexprhdlr, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC,
         EXPRHDLR_PRECEDENCE, evalSum, NULL) );
   assert(exprhdlr != NULL);

   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeHdlr(scip, consexprhdlr, exprhdlr, copyhdlrSum, NULL) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataSum, freedataSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifySum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, initSepaSum, NULL, estimateSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrFwdiff(scip, consexprhdlr, exprhdlr, fwdiffSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwfwdiff(scip, consexprhdlr, exprhdlr, bwfwdiffSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCurvature(scip, consexprhdlr, exprhdlr, curvatureSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrMonotonicity(scip, consexprhdlr, exprhdlr, monotonicitySum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntegrality(scip, consexprhdlr, exprhdlr, integralitySum) );

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

/** multiplies given sum expression by a constant */
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

/** reverse propagate a weighted sum of expressions in the given interval */
SCIP_RETCODE SCIPreverseConsExprExprPropagateWeightedSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   int                   nexprs,             /**< number of expressions to propagate */
   SCIP_CONSEXPR_EXPR**  exprs,              /**< expressions to propagate */
   SCIP_Real*            weights,            /**< weights of expressions in sum */
   SCIP_Real             constant,           /**< constant in sum */
   SCIP_INTERVAL         interval,           /**< constant + sum weight_i expr_i in interval */
   SCIP_Bool*            infeasible,         /**< buffer to store if propagation produced infeasibility */
   int*                  nreductions         /**< buffer to store the number of interval reductions */
   )
{
   SCIP_INTERVAL* bounds;
   SCIP_ROUNDMODE prevroundmode;
   SCIP_INTERVAL childbounds;
   SCIP_Real minlinactivity;
   SCIP_Real maxlinactivity;
   int minlinactivityinf;
   int maxlinactivityinf;
   int c;

   assert(scip != NULL);
   assert(exprs != NULL || nexprs == 0);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   *infeasible = FALSE;
   *nreductions = 0;

   /* not possible to conclude finite bounds if the interval is [-inf,inf] */
   if( SCIPintervalIsEntire(SCIP_INTERVAL_INFINITY, interval) )
      return SCIP_OKAY;

   prevroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();

   minlinactivity = constant;
   maxlinactivity = -constant; /* use -constant because of the rounding mode */
   minlinactivityinf = 0;
   maxlinactivityinf = 0;

   SCIPdebugMsg(scip, "reverse prop with %d children: %.20g", constant);

   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, nexprs) );

   /* shift coefficients into the intervals of the children; compute the min and max activities */
   for( c = 0; c < nexprs; ++c )
   {
      childbounds = SCIPgetConsExprExprBounds(scip, conshdlr, exprs[c]);
      SCIPdebugMsgPrint(scip, " %+.20g*[%.20g,%.20g]", weights[c], childbounds.inf, childbounds.sup); /*lint !e613 */

      if( SCIPintervalIsEmpty(SCIP_INTERVAL_INFINITY, childbounds) )
      {
         *infeasible = TRUE;
         goto TERMINATE;
      }

      SCIPintervalMulScalar(SCIP_INTERVAL_INFINITY, &bounds[c], childbounds, weights[c]);  /*lint !e613 */

      if( SCIPisInfinity(scip, SCIPintervalGetSup(bounds[c])) )
         ++maxlinactivityinf;
      else
      {
         assert(SCIPintervalGetSup(bounds[c]) > -SCIP_INTERVAL_INFINITY);
         maxlinactivity -= SCIPintervalGetSup(bounds[c]);
      }

      if( SCIPisInfinity(scip, -SCIPintervalGetInf(bounds[c])) )
         ++minlinactivityinf;
      else
      {
         assert(SCIPintervalGetInf(bounds[c]) < SCIP_INTERVAL_INFINITY);
         minlinactivity += SCIPintervalGetInf(bounds[c]);
      }
   }
   maxlinactivity = -maxlinactivity; /* correct sign */

   SCIPdebugMsgPrint(scip, " = [%.20g,%.20g] in rhs = [%.20g,%.20g]\n",
      minlinactivityinf ? -SCIP_INTERVAL_INFINITY : minlinactivity,
      maxlinactivityinf ?  SCIP_INTERVAL_INFINITY : maxlinactivity,
      interval.inf, interval.sup);

   /* if there are too many unbounded bounds, then could only compute infinite bounds for children, so give up */
   if( (minlinactivityinf >= 2 || SCIPisInfinity(scip, SCIPintervalGetSup(interval))) &&
      ( maxlinactivityinf >= 2 || SCIPisInfinity(scip, -SCIPintervalGetInf(interval)))
      )
   {
      goto TERMINATE;
   }

   for( c = 0; c < nexprs && !(*infeasible); ++c )
   {
      /* upper bounds of c_i is
       *   node->bounds.sup - (minlinactivity - c_i.inf), if c_i.inf > -infinity and minlinactivityinf == 0
       *   node->bounds.sup - minlinactivity, if c_i.inf == -infinity and minlinactivityinf == 1
       */
      SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &childbounds);
      if( !SCIPisInfinity(scip, SCIPintervalGetSup(interval)) )
      {
         /* we are still in downward rounding mode, so negate and negate to get upward rounding */
         if( bounds[c].inf <= -SCIP_INTERVAL_INFINITY && minlinactivityinf <= 1 )
         {
            assert(minlinactivityinf == 1);
            childbounds.sup = SCIPintervalNegateReal(minlinactivity - interval.sup);
         }
         else if( minlinactivityinf == 0 )
         {
            childbounds.sup = SCIPintervalNegateReal(minlinactivity - interval.sup - bounds[c].inf);
         }
      }

      /* lower bounds of c_i is
       *   node->bounds.inf - (maxlinactivity - c_i.sup), if c_i.sup < infinity and maxlinactivityinf == 0
       *   node->bounds.inf - maxlinactivity, if c_i.sup == infinity and maxlinactivityinf == 1
       */
      if( interval.inf > -SCIP_INTERVAL_INFINITY )
      {
         if( bounds[c].sup >= SCIP_INTERVAL_INFINITY && maxlinactivityinf <= 1 )
         {
            assert(maxlinactivityinf == 1);
            childbounds.inf = interval.inf - maxlinactivity;
         }
         else if( maxlinactivityinf == 0 )
         {
            childbounds.inf = interval.inf - maxlinactivity + bounds[c].sup;
         }
      }

      SCIPdebugMsg(scip, "child %d: %.20g*x in [%.20g,%.20g]", c, weights[c], childbounds.inf, childbounds.sup);

      /* divide by the child coefficient */
      SCIPintervalDivScalar(SCIP_INTERVAL_INFINITY, &childbounds, childbounds, weights[c]);

      SCIPdebugMsgPrint(scip, " -> x = [%.20g,%.20g]\n", childbounds.inf, childbounds.sup);

      /* try to tighten the bounds of the expression */
      SCIP_CALL( SCIPtightenConsExprExprInterval(scip, conshdlr, exprs[c], childbounds, infeasible, nreductions) );  /*lint !e613 */
   }

TERMINATE:
   SCIPintervalSetRoundingMode(prevroundmode);
   SCIPfreeBufferArray(scip, &bounds);

   return SCIP_OKAY;
}
