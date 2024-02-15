/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   expr_sum.c
 * @ingroup DEFPLUGINS_EXPR
 * @brief  sum expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>
#include <stddef.h>

#include "scip/expr_sum.h"
#include "scip/expr_value.h"
#include "scip/expr_product.h"
#include "scip/expr_exp.h"
#include "scip/expr_pow.h"
#include "symmetry/struct_symmetry.h"

#define EXPRHDLR_NAME         "sum"
#define EXPRHDLR_DESC         "summation with coefficients and a constant"
#define EXPRHDLR_PRECEDENCE   40000
#define EXPRHDLR_HASHKEY      SCIPcalcFibHash(47161.0)

/** macro to activate/deactivate debugging information of simplify method */
/*lint -emacro(774,debugSimplify) */
#ifdef SIMPLIFY_DEBUG
#define debugSimplify                   printf
#else
#define debugSimplify                   while( FALSE ) printf
#endif

/*
 * Data structures
 */

/** expression data */
struct SCIP_ExprData
{
   SCIP_Real             constant;           /**< constant coefficient */
   SCIP_Real*            coefficients;       /**< coefficients of children */
   int                   coefssize;          /**< size of the coefficients array */
};

/*
 * Local methods
 */

/** creates expression data */
static
SCIP_RETCODE createData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRDATA**       exprdata,           /**< pointer where to store expression data */
   int                   ncoefficients,      /**< number of coefficients (i.e., number of children) */
   SCIP_Real*            coefficients,       /**< array with coefficients for all children (or NULL if all 1.0) */
   SCIP_Real             constant            /**< constant term of sum */
   )
{
   assert(exprdata != NULL);
   assert(ncoefficients >= 0);

   SCIP_CALL( SCIPallocBlockMemory(scip, exprdata) );

   if( coefficients != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*exprdata)->coefficients, coefficients, ncoefficients) );
   }
   else
   {
      int i;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*exprdata)->coefficients, ncoefficients) );
      for( i = 0; i < ncoefficients; ++i )
         (*exprdata)->coefficients[i] = 1.0;
   }

   (*exprdata)->coefssize = ncoefficients;
   (*exprdata)->constant  = constant;

   return SCIP_OKAY;
}

/** simplifies the `idx`-th child of the sum expression `duplicate` in order for it to be able to be a child of a simplified sum
 *
 * for example, this means that the `idx`-th child cannot be itself a sum
 * if it is, we have to flatten it, i.e., take all its children and make them children of `duplicate`
 */
static
SCIP_RETCODE simplifyTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            duplicate,          /**< expression to be simplified */
   int                   idx,                /**< idx of children to be simplified */
   SCIP_Bool*            changed,            /**< pointer to store if some term actually got simplified */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPR** children;
   SCIP_EXPR* expr;
   SCIP_Real* coefs;
   SCIP_Real constant ;
   SCIP_Real coef;

   assert(duplicate != NULL);
   assert(idx >= 0);
   assert(idx < SCIPexprGetNChildren(duplicate));
   assert(changed != NULL);

   children  = SCIPexprGetChildren(duplicate);
   coefs     = SCIPgetCoefsExprSum(duplicate);
   constant  = SCIPgetConstantExprSum(duplicate);

   coef = coefs[idx];
   expr = children[idx];
   assert(expr != NULL);

   /* enforces SS3 */
   if( SCIPisExprValue(scip, expr) )
   {
      *changed = TRUE;
      constant += coef * SCIPgetValueExprValue(expr);
      SCIPsetConstantExprSum(duplicate, constant);

      /* TODO: remove child? */
      coefs[idx] = 0.0;

      return SCIP_OKAY;
   }

   /* enforces SS2 */
   if( SCIPisExprSum(scip, expr) )
   {
      *changed = TRUE;

      /* pass constant to parent */
      constant += coef * SCIPgetConstantExprSum(expr);
      SCIPsetConstantExprSum(duplicate, constant);

      /* append all children of expr on parent except the first one */
      if( SCIPexprGetNChildren(expr) > 1 )
      {
         int i;

         for( i = 1; i < SCIPexprGetNChildren(expr); ++i )
         {
            assert(!SCIPisExprSum(scip, SCIPexprGetChildren(expr)[i]));
            SCIP_CALL( SCIPappendExprSumExpr(scip, duplicate, SCIPexprGetChildren(expr)[i],
                  coef * SCIPgetCoefsExprSum(expr)[i]) );
         }
      }

      /* replace expr with first child; need to get data again since it might be re-allocated */
      assert(!SCIPisExprSum(scip, SCIPexprGetChildren(expr)[0]));

      coefs = SCIPgetCoefsExprSum(duplicate);

      coefs[idx] = coef * SCIPgetCoefsExprSum(expr)[0];
      SCIP_CALL( SCIPreplaceExprChild(scip, duplicate, idx, SCIPexprGetChildren(expr)[0]) );

      return SCIP_OKAY;
   }

   /* enforce SS9 */
   if( REALABS(coef) != 1.0 && SCIPisExprProduct(scip, expr) )
   {
      SCIP_EXPR* expchild = NULL;
      int i;

      for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
      {
         SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
         assert(child != NULL);

         if( SCIPisExprExp(scip, child) )
         {
            expchild = child;
            break;
         }
      }

      /* coef != +- 1, term is product and one factor is an exponential -> enforce SS9 */
      if( expchild != NULL )
      {
         SCIP_EXPR* sum;
         SCIP_EXPR* prod;
         SCIP_EXPR* simplifiedprod;
         SCIP_EXPR* simplifiedsum;
         SCIP_EXPR* exponential;
         SCIP_EXPR* simplifiedexp;
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
         SCIP_CALL( SCIPcreateExprSum(scip, &sum, 1, SCIPexprGetChildren(expchild), NULL, expconstant, ownercreate,
               ownercreatedata) );

         /* simplify sum */
         SCIP_CALL( SCIPcallExprSimplify(scip, sum, &simplifiedsum, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &sum) );

         /* create exponential with new child */
         SCIP_CALL( SCIPcreateExprExp(scip, &exponential, simplifiedsum, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedsum) );

         /* simplify exponential */
         SCIP_CALL( SCIPcallExprSimplify(scip, exponential, &simplifiedexp, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &exponential) );

         /* create product with new child */
         SCIP_CALL( SCIPcreateExprProduct(scip, &prod, 0, NULL, 1.0, ownercreate, ownercreatedata) );

         for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
         {
            if( SCIPexprGetChildren(expr)[i] == expchild )
            {
               SCIP_CALL( SCIPappendExprChild(scip, prod, simplifiedexp) );
            }
            else
            {
               SCIP_CALL( SCIPappendExprChild(scip, prod, SCIPexprGetChildren(expr)[i]) );
            }
         }
         SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedexp) );

         /* simplify product */
         SCIP_CALL( SCIPcallExprSimplify(scip, prod, &simplifiedprod, ownercreate, ownercreatedata) );
         SCIP_CALL( SCIPreleaseExpr(scip, &prod) );

         /* replace current child with simplified product */
         SCIP_CALL( SCIPreplaceExprChild(scip, duplicate, idx, simplifiedprod) );
         SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedprod) );

         /* since the simplified product can be a sum ( exp(-1)*exp(log(x+y)+1) -> x+y ),
          * we call the function we are in again
          * this is no endless recursion, since the coef is now +- 1
          */
         SCIP_CALL( simplifyTerm(scip, duplicate, idx, changed, ownercreate, ownercreatedata) );

         return SCIP_OKAY;
      }
   }

   /* enforce SS10 */
   if( REALABS(coef) != 1.0 && SCIPisExprExp(scip, expr) )
   {
      /* coef != +- 1, term is exponential -> enforce SS10 by moving |coef| into argument of exponential */

      SCIP_EXPR* sum;
      SCIP_EXPR* simplifiedsum;
      SCIP_EXPR* exponential;
      SCIP_EXPR* simplifiedexp;
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
      SCIP_CALL( SCIPcreateExprSum(scip, &sum, 1, SCIPexprGetChildren(expr), NULL, expconstant, ownercreate,
            ownercreatedata) );  /* expconstant+expchild */

      /* simplify sum */
      SCIP_CALL( SCIPcallExprSimplify(scip, sum, &simplifiedsum, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &sum) );

      /* create exponential with new child */
      SCIP_CALL( SCIPcreateExprExp(scip, &exponential, simplifiedsum, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedsum) );

      /* simplify exponential */
      SCIP_CALL( SCIPcallExprSimplify(scip, exponential, &simplifiedexp, ownercreate, ownercreatedata) );
      SCIP_CALL( SCIPreleaseExpr(scip, &exponential) );

      /* replace current child with simplified exponential */
      SCIP_CALL( SCIPreplaceExprChild(scip, duplicate, idx, simplifiedexp) );
      SCIP_CALL( SCIPreleaseExpr(scip, &simplifiedexp) );

      return SCIP_OKAY;
   }

   /* other types of (simplified) expressions can be a child of a simplified sum */
   assert(!SCIPisExprSum(scip, expr));
   assert(!SCIPisExprValue(scip, expr));

   return SCIP_OKAY;
}

/** helper struct for expressions sort */
typedef struct
{
   SCIP*                 scip;               /**< SCIP data structure */
   SCIP_EXPR**           exprs;              /**< expressions */
} SORTEXPRDATA;

static
SCIP_DECL_SORTINDCOMP(sortExprComp)
{
   SORTEXPRDATA* data = (SORTEXPRDATA*)dataptr;

   return SCIPcompareExpr(data->scip, data->exprs[ind1], data->exprs[ind2]);
}

/*
 * Callback methods of expression handler
 */

/** simplifies a sum expression
 *
 * goes through each child and simplifies it; then sorts the simplified children; then sum the children that are equal;
 * finally creates a sum expression with all the children that do not have a 0 coefficient and post-process so that SS6
 * and SS7 are satisfied
 */
static
SCIP_DECL_EXPRSIMPLIFY(simplifySum)
{  /*lint --e{715}*/
   SCIP_EXPR** children;
   SCIP_EXPR* duplicate = NULL;
   SCIP_EXPR** newchildren = NULL;
   SCIP_Real* newcoefs = NULL;
   int nnewchildren;
   SCIP_Real newconstant;
   SCIP_Real* coefs;
   int i;
   int nchildren;
   SCIP_Bool changed;
   SORTEXPRDATA sortdata;
   int* order = NULL;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(SCIPexprGetHdlr(expr) == SCIPgetExprhdlrSum(scip));

   changed = FALSE;

   /* TODO: maybe have a flag to know if it is simplified ? */
   /* TODO: can we do this with a shallow duplicate + copy of children pointer? currently simplifyTerm may modify children,
    * so one would need to be careful
    */
   SCIP_CALL( SCIPduplicateExpr(scip, expr, &duplicate, NULL, NULL, ownercreate, ownercreatedata) );
   assert(duplicate != NULL);

   nchildren = SCIPexprGetNChildren(duplicate);
   for( i = 0; i < nchildren; i++ )
   {
      /* enforces SS8 TODO: remove child? */
      /* we have to ask for the coefs everytime, since it might get realloced in simpifyTerm */
      if( SCIPgetCoefsExprSum(duplicate)[i] == 0.0 )
      {
         changed = TRUE;
         continue;
      }

      /* enforces SS2, SS3, SS9, and SS10 */
      SCIP_CALL( simplifyTerm(scip, duplicate, i, &changed, ownercreate, ownercreatedata) );
   }

   /* simplifyTerm can add new children to duplicate and realloc them; so get them again */
   nchildren = SCIPexprGetNChildren(duplicate);
   children  = SCIPexprGetChildren(duplicate);
   coefs     = SCIPgetCoefsExprSum(duplicate);

   /* treat zero term case */
   if( nchildren == 0 )
   {
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, SCIPgetConstantExprSum(duplicate), ownercreate, ownercreatedata) );
      goto CLEANUP;
   }

   /* treat one term case */
   if( nchildren == 1 )
   {
      if( coefs[0] == 0.0 )
      {
         SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, SCIPgetConstantExprSum(duplicate), ownercreate, ownercreatedata) );
         goto CLEANUP;
      }

      if( coefs[0] == 1.0 && SCIPgetConstantExprSum(duplicate) == 0.0 )
         *simplifiedexpr = children[0]; /* SS7 */
      else
         *simplifiedexpr = changed ? duplicate : expr;

      SCIPcaptureExpr(*simplifiedexpr);

      goto CLEANUP;
   }

   /* enforces SS5: sort children */
   SCIP_CALL( SCIPallocBufferArray(scip, &order, nchildren) );
   for( i = 0; i < nchildren; i++ )
      order[i] = i;
   sortdata.scip = scip;
   sortdata.exprs = children;
   SCIPsortInd(order, sortExprComp, (void*)&sortdata, nchildren);

   /* create sorted variant of children and coefs */
   SCIP_CALL( SCIPallocBufferArray(scip, &newchildren, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &newcoefs, nchildren) );
   for( i = 0; i < nchildren; ++i )
   {
      newchildren[i] = children[order[i]];
      newcoefs[i] = coefs[order[i]];
      if( order[i] != i )
         changed = TRUE;
   }

   /* post-process */

   /* enforces SS4 */
   nnewchildren = 0;
   for( i = 0; i < nchildren; i++ )
   {
      /* eliminate zero-coefficients */
      if( newcoefs[i] == 0.0 )
      {
         changed = TRUE;
         continue;
      }

      /* sum equal expressions */
      if( i < nchildren-1 && SCIPcompareExpr(scip, newchildren[i], newchildren[i+1]) == 0 )
      {
         changed = TRUE;
         /* if we substract two almost equal not-so-small numbers, then set new coefficient to 0.0
          * instead of some tiny value that is likely the result of some random round-off error
          * E.g., on instance ex1221, we have x1^2 + b3 = 1.25.
          *   Probing finds an aggregation x1 = 1.11803 - 0.618034 b3.
          *   Simplify would then produce 1.25 + 1e-16 x1 = 1.25.
          */
         if( SCIPisEQ(scip, newcoefs[i], -newcoefs[i+1]) && REALABS(newcoefs[i]) >= 1.0 )
            newcoefs[i+1] = 0.0;
         else
            newcoefs[i+1] += newcoefs[i];
         continue;
      }

      /* move i-th child to new position */
      newchildren[nnewchildren] = newchildren[i];
      newcoefs[nnewchildren] = newcoefs[i];
      nnewchildren++;
   }

   /* build sum expression from finalchildren and post-simplify */
   newconstant = SCIPgetConstantExprSum(duplicate);

   debugSimplify("what to do? finalchildren has length %d\n", nnewchildren); /*lint !e506 !e681*/

   /* enforces SS6: if they are no children, return value */
   if( nnewchildren == 0 )
   {
      debugSimplify("[sum] got empty list, return value %g\n", newconstant); /*lint !e506 !e681*/
      SCIP_CALL( SCIPcreateExprValue(scip, simplifiedexpr, newconstant, ownercreate, ownercreatedata) );

      goto CLEANUP;
   }

   /* enforces SS7: if list consists of one expr with coef 1.0 and constant is 0, return that expr */
   if( nnewchildren == 1 && newcoefs[0] == 1.0 && newconstant == 0.0 )
   {
      *simplifiedexpr = newchildren[0];
      SCIPcaptureExpr(*simplifiedexpr);

      goto CLEANUP;
   }

   /* build sum expression from children */
   if( changed )
   {
      SCIP_CALL( SCIPcreateExprSum(scip, simplifiedexpr, nnewchildren, newchildren, newcoefs, newconstant,
            ownercreate, ownercreatedata) );

      goto CLEANUP;
   }

   *simplifiedexpr = expr;

   /* we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created */
   SCIPcaptureExpr(*simplifiedexpr);

   /* free memory */
 CLEANUP:
   SCIPfreeBufferArrayNull(scip, &newcoefs);
   SCIPfreeBufferArrayNull(scip, &newchildren);
   SCIPfreeBufferArrayNull(scip, &order);
   SCIP_CALL( SCIPreleaseExpr(scip, &duplicate) );

   return SCIP_OKAY;
}

/** expression callback to get information for symmetry detection */
static
SCIP_DECL_EXPRGETSYMDATA(getSymDataSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, symdata) );

   (*symdata)->nconstants = 1;
   (*symdata)->ncoefficients = exprdata->coefssize;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*symdata)->constants, 1) );
   (*symdata)->constants[0] = exprdata->constant;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*symdata)->coefficients, exprdata->coefssize) );
   for( i = 0; i < exprdata->coefssize; ++i )
      (*symdata)->coefficients[i] = exprdata->coefficients[i];

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*symdata)->children, exprdata->coefssize) );
   for( i = 0; i < exprdata->coefssize; ++i )
      (*symdata)->children[i] = SCIPexprGetChildren(expr)[i];

   return SCIP_OKAY;
}

/** compares two sum expressions
 *
 *  The order of two sum expressions is a lexicographical order on the terms.
 *
 *  Starting from the *last*, we find the first child where they differ, say, the i-th.
 *  Then u < v <=> u_i < v_i.
 *  If there are no such children and they have different number of children, then u < v <=> nchildren(u) < nchildren(v).
 *  If there are no such children and they have the same number of children, then u < v <=> const(u) < const(v).
 *  Otherwise, they are the same.
 *
 *  Note: we are assuming expression are simplified, so within u, we have u_1 < u_2, etc
 *
 *  Example: y + z < x + y + z, 2*x + 3*y < 3*x + 3*y
 */
static
SCIP_DECL_EXPRCOMPARE(compareSum)
{  /*lint --e{715}*/
   SCIP_Real const1;
   SCIP_Real* coefs1;
   SCIP_EXPR** children1;
   int nchildren1;
   SCIP_Real const2;
   SCIP_Real* coefs2;
   SCIP_EXPR** children2;
   int nchildren2;
   int compareresult;
   int i;
   int j;

   nchildren1 = SCIPexprGetNChildren(expr1);
   nchildren2 = SCIPexprGetNChildren(expr2);
   children1 = SCIPexprGetChildren(expr1);
   children2 = SCIPexprGetChildren(expr2);
   coefs1 = SCIPgetCoefsExprSum(expr1);
   coefs2 = SCIPgetCoefsExprSum(expr2);
   const1 = SCIPgetConstantExprSum(expr1);
   const2 = SCIPgetConstantExprSum(expr2);

   for( i = nchildren1 - 1, j = nchildren2 - 1; i >= 0 && j >= 0; --i, --j )
   {
      compareresult = SCIPcompareExpr(scip, children1[i], children2[j]);
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

/** expression handler copy callback */
static
SCIP_DECL_EXPRCOPYHDLR(copyhdlrSum)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeExprhdlrSum(scip) );

   return SCIP_OKAY;
}

/** expression data copy callback */
static
SCIP_DECL_EXPRCOPYDATA(copydataSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* sourceexprdata;

   assert(targetexprdata != NULL);
   assert(sourceexpr != NULL);

   sourceexprdata = SCIPexprGetData(sourceexpr);
   assert(sourceexprdata != NULL);

   SCIP_CALL( createData(targetscip, targetexprdata, SCIPexprGetNChildren(sourceexpr),
            sourceexprdata->coefficients, sourceexprdata->constant) );

   return SCIP_OKAY;
}

/** expression data free callback */
static
SCIP_DECL_EXPRFREEDATA(freedataSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   SCIPfreeBlockMemoryArray(scip, &(exprdata->coefficients), exprdata->coefssize);
   SCIPfreeBlockMemory(scip, &exprdata);

   SCIPexprSetData(expr, NULL);

   return SCIP_OKAY;
}

/** expression print callback */
static
SCIP_DECL_EXPRPRINT(printSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   /**! [SnippetExprPrintSum] */
   switch( stage )
   {
      case SCIP_EXPRITER_ENTEREXPR :
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

      case SCIP_EXPRITER_VISITINGCHILD :
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

      case SCIP_EXPRITER_LEAVEEXPR :
      {
         /* print closing parenthesis, if necessary */
         if( EXPRHDLR_PRECEDENCE <= parentprecedence )
         {
            SCIPinfoMessage(scip, file, ")");
         }
         break;
      }

      case SCIP_EXPRITER_VISITEDCHILD :
      default: ;
   }
   /**! [SnippetExprPrintSum] */

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_EXPREVAL(evalSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int c;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   /**! [SnippetExprEvalSum] */
   *val = exprdata->constant;
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      assert(SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[c]) != SCIP_INVALID); /*lint !e777*/

      *val += exprdata->coefficients[c] * SCIPexprGetEvalValue(SCIPexprGetChildren(expr)[c]);
   }
   /**! [SnippetExprEvalSum] */

   return SCIP_OKAY;
}

/** expression forward derivative evaluation callback */
static
SCIP_DECL_EXPRFWDIFF(fwdiffSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int c;

   assert(expr != NULL);
   assert(dot != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   *dot = 0.0;
   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      assert(SCIPexprGetDot(SCIPexprGetChildren(expr)[c]) != SCIP_INVALID); /*lint !e777*/

      *dot += exprdata->coefficients[c] * SCIPexprGetDot(SCIPexprGetChildren(expr)[c]);
   }

   return SCIP_OKAY;
}

/** expression derivative evaluation callback */
static
SCIP_DECL_EXPRBWDIFF(bwdiffSum)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPexprGetData(expr) != NULL);
   assert(childidx >= 0 && childidx < SCIPexprGetNChildren(expr));
   assert(SCIPexprGetChildren(expr)[childidx] != NULL);
   assert(!SCIPisExprValue(scip, SCIPexprGetChildren(expr)[childidx]));

   *val = SCIPgetCoefsExprSum(expr)[childidx];

   return SCIP_OKAY;
}

/** expression backward forward derivative evaluation callback */
static
SCIP_DECL_EXPRBWFWDIFF(bwfwdiffSum)
{  /*lint --e{715}*/
   assert(bardot != NULL);

   *bardot = 0.0;

   return SCIP_OKAY;
}

/** expression interval evaluation callback */
static
SCIP_DECL_EXPRINTEVAL(intevalSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   SCIP_INTERVAL suminterval;
   int c;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   SCIPintervalSet(interval, exprdata->constant);

   SCIPdebugMsg(scip, "inteval %p with %d children: %.20g", (void*)expr, SCIPexprGetNChildren(expr), exprdata->constant);

   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
   {
      SCIP_INTERVAL childinterval;

      childinterval = SCIPexprGetActivity(SCIPexprGetChildren(expr)[c]);
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

/** initial estimators callback */
static
SCIP_DECL_EXPRINITESTIMATES(initEstimatesSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "initEstimatesSum %d children: ", SCIPexprGetNChildren(expr));
   SCIPprintExpr(scip, expr, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif
   assert(scip != NULL);
   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);

   assert(coefs[0] != NULL);
   assert(constant != NULL);
   assert(nreturned != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   BMScopyMemoryArray(coefs[0], exprdata->coefficients, SCIPexprGetNChildren(expr));
   *constant = exprdata->constant;
   *nreturned = 1;

   return SCIP_OKAY;
}

/** expression estimate callback */
static
SCIP_DECL_EXPRESTIMATE(estimateSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(strcmp(SCIPexprhdlrGetName(SCIPexprGetHdlr(expr)), EXPRHDLR_NAME) == 0);
   assert(islocal != NULL);
   assert(success != NULL);
   assert(branchcand != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   /* NOTE: nlhdlr_default assumes in nlhdlrInitSepaDefault that this estimator can be used for both under- and overestimation */

   BMScopyMemoryArray(coefs, exprdata->coefficients, SCIPexprGetNChildren(expr));
   *constant = exprdata->constant;
   *islocal = FALSE;
   *success = TRUE;

   /* for none of our children, branching would improve the underestimator, so set branchcand[i]=FALSE everywhere
    * if we branch for numerical reasons, then cons-expr-core should figure out what the candidates are
    */
   BMSclearMemoryArray(branchcand, SCIPexprGetNChildren(expr));

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_EXPRREVERSEPROP(reversepropSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   SCIP_INTERVAL* newbounds;
   int nchildren;
   int nreductions;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(infeasible != NULL);

   nchildren = SCIPexprGetNChildren(expr);
   assert(nchildren > 0);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &newbounds, nchildren) );

   nreductions = SCIPintervalPropagateWeightedSum(SCIP_INTERVAL_INFINITY, nchildren, childrenbounds,
         exprdata->coefficients, exprdata->constant, bounds, newbounds, infeasible);

   if( !*infeasible && nreductions > 0 )
      BMScopyMemoryArray(childrenbounds, newbounds, nchildren);

   SCIPfreeBufferArray(scip, &newbounds);

   return SCIP_OKAY;
}

/** sum hash callback */
static
SCIP_DECL_EXPRHASH(hashSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(hashkey != NULL);
   assert(childrenhashes != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   /**! [SnippetExprHashSum] */
   *hashkey = EXPRHDLR_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(exprdata->constant);

   for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
      *hashkey ^= SCIPcalcFibHash(exprdata->coefficients[c]) ^ childrenhashes[c];
   /**! [SnippetExprHashSum] */

   return SCIP_OKAY;
}

/** expression curvature detection callback */
static
SCIP_DECL_EXPRCURVATURE(curvatureSum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(childcurv != NULL);
   assert(success != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
      childcurv[i] = SCIPexprcurvMultiply(exprdata->coefficients[i], exprcurvature);

   *success = TRUE;

   return SCIP_OKAY;
}

/** expression monotonicity detection callback */
static
SCIP_DECL_EXPRMONOTONICITY(monotonicitySum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(result != NULL);
   assert(childidx >= 0 && childidx < SCIPexprGetNChildren(expr));

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   *result = exprdata->coefficients[childidx] >= 0.0 ? SCIP_MONOTONE_INC : SCIP_MONOTONE_DEC;

   return SCIP_OKAY;
}

/** expression integrality detection callback */
static
SCIP_DECL_EXPRINTEGRALITY(integralitySum)
{  /*lint --e{715}*/
   SCIP_EXPRDATA* exprdata;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(isintegral != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   /**! [SnippetExprIntegralitySum] */
   *isintegral = EPSISINT(exprdata->constant, 0.0); /*lint !e835*/

   for( i = 0; i < SCIPexprGetNChildren(expr) && *isintegral; ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
      assert(child != NULL);

      *isintegral = EPSISINT(exprdata->coefficients[i], 0.0) && SCIPexprIsIntegral(child); /*lint !e835*/
   }
   /**! [SnippetExprIntegralitySum] */

   return SCIP_OKAY;
}

/** creates the handler for sum expressions and includes it into SCIP */
SCIP_RETCODE SCIPincludeExprhdlrSum(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EXPRHDLR* exprhdlr;

   SCIP_CALL( SCIPincludeExprhdlr(scip, &exprhdlr, EXPRHDLR_NAME, EXPRHDLR_DESC, EXPRHDLR_PRECEDENCE, evalSum, NULL) );
   assert(exprhdlr != NULL);

   SCIPexprhdlrSetCopyFreeHdlr(exprhdlr, copyhdlrSum, NULL);
   SCIPexprhdlrSetCopyFreeData(exprhdlr, copydataSum, freedataSum);
   SCIPexprhdlrSetSimplify(exprhdlr, simplifySum);
   SCIPexprhdlrSetCompare(exprhdlr, compareSum);
   SCIPexprhdlrSetPrint(exprhdlr, printSum);
   SCIPexprhdlrSetIntEval(exprhdlr, intevalSum);
   SCIPexprhdlrSetEstimate(exprhdlr, initEstimatesSum, estimateSum);
   SCIPexprhdlrSetReverseProp(exprhdlr, reversepropSum);
   SCIPexprhdlrSetHash(exprhdlr, hashSum);
   SCIPexprhdlrSetDiff(exprhdlr, bwdiffSum, fwdiffSum, bwfwdiffSum);
   SCIPexprhdlrSetCurvature(exprhdlr, curvatureSum);
   SCIPexprhdlrSetMonotonicity(exprhdlr, monotonicitySum);
   SCIPexprhdlrSetIntegrality(exprhdlr, integralitySum);
   SCIPexprhdlrSetGetSymdata(exprhdlr, getSymDataSum);

   return SCIP_OKAY;
}

/** creates a sum expression */
SCIP_RETCODE SCIPcreateExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children */
   SCIP_Real*            coefficients,       /**< array with coefficients for all children (or NULL if all 1.0) */
   SCIP_Real             constant,           /**< constant term of sum */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_EXPRDATA* exprdata;

   SCIP_CALL( createData(scip, &exprdata, nchildren, coefficients, constant) );

   SCIP_CALL( SCIPcreateExpr(scip, expr, SCIPgetExprhdlrSum(scip), exprdata, nchildren, children, ownercreate, ownercreatedata) );

   return SCIP_OKAY;
}

/** sets the constant of a summation expression */
void SCIPsetConstantExprSum(
   SCIP_EXPR*            expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   exprdata->constant = constant;
}

/** appends an expression to a sum expression */
SCIP_RETCODE SCIPappendExprSumExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< sum expression */
   SCIP_EXPR*            child,              /**< expression to be appended */
   SCIP_Real             childcoef           /**< child's coefficient */
   )
{
   SCIP_EXPRDATA* exprdata;
   int nchildren;

   assert(expr != NULL);
   assert(SCIPisExprSum(scip, expr));

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   nchildren = SCIPexprGetNChildren(expr);

   SCIP_CALL( SCIPensureBlockMemoryArray(scip, &exprdata->coefficients, &exprdata->coefssize, nchildren + 1) );

   assert(exprdata->coefssize > nchildren);
   exprdata->coefficients[nchildren] = childcoef;

   SCIP_CALL( SCIPappendExprChild(scip, expr, child) );

   return SCIP_OKAY;
}

/** multiplies given sum expression by a constant */
void SCIPmultiplyByConstantExprSum(
   SCIP_EXPR*            expr,               /**< sum expression */
   SCIP_Real             constant            /**< constant that multiplies sum expression */
   )
{
   int i;
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
      exprdata->coefficients[i] *= constant;
   exprdata->constant *= constant;
}

/** constructs the expanded product of two sum expressions */
SCIP_RETCODE SCIPmultiplyBySumExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           product,            /**< buffer where to store multiplied sums (expanded as sum) */
   SCIP_EXPR*            factor1,            /**< first sum */
   SCIP_EXPR*            factor2,            /**< second sum */
   SCIP_Bool             simplify,           /**< whether to simplify created terms and sum */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   SCIP_Real constant1;
   SCIP_Real constant2;
   int nchildren1;
   int nchildren2;
   int i1;
   int i2;

   assert(scip != NULL);
   assert(product != NULL);
   assert(factor1 != NULL);
   assert(SCIPisExprSum(scip, factor1));
   assert(factor2 != NULL);
   assert(SCIPisExprSum(scip, factor2));

   constant1 = SCIPgetConstantExprSum(factor1);
   constant2 = SCIPgetConstantExprSum(factor2);
   nchildren1 = SCIPexprGetNChildren(factor1);
   nchildren2 = SCIPexprGetNChildren(factor2);

   /* TODO might be nice to integrate more with simplify and construct a simplified sum right away */

   SCIP_CALL( SCIPcreateExprSum(scip, product, 0, NULL, NULL, constant1 * constant2, ownercreate, ownercreatedata) );

   /* first add constant1 * factor2
    * constant * constant2 already added above
    */
   if( constant1 != 0.0 )
   {
      for( i2 = 0; i2 < nchildren2; ++i2 )
      {
         SCIP_CALL( SCIPappendExprSumExpr(scip, *product, SCIPexprGetChildren(factor2)[i2], constant1 * SCIPgetCoefsExprSum(factor2)[i2]) );
      }
   }

   for( i1 = 0; i1 < nchildren1; ++i1 )
   {
      SCIP_EXPR* child1;
      SCIP_EXPR* child2;
      SCIP_Real coef1;
      SCIP_Real coef2;

      coef1 = SCIPgetCoefsExprSum(factor1)[i1];
      child1 = SCIPexprGetChildren(factor1)[i1];

      if( constant2 != 0.0 )
      {
         /* add coef1 * child1 * constant2 */
         SCIP_CALL( SCIPappendExprSumExpr(scip, *product, child1, coef1 * constant2) );
      }

      for( i2 = 0; i2 < nchildren2; ++i2 )
      {
         /* add coef1 * child1 * coef2 * child2 */
         SCIP_EXPR* termprod;
         SCIP_EXPR* termprodsimplified;
         SCIP_EXPR* termfactors[2];

         coef2 = SCIPgetCoefsExprSum(factor2)[i2];
         child2 = SCIPexprGetChildren(factor2)[i2];

         /* create child1 * child2 expr */
         termfactors[0] = child1;
         termfactors[1] = child2;
         SCIP_CALL( SCIPcreateExprProduct(scip, &termprod, 2, termfactors, 1.0, NULL, NULL) );

         if( simplify )
         {
            SCIP_Bool changed;
            SCIP_Bool infeasible;

            /* simplify child1*child2 */
            SCIP_CALL( SCIPsimplifyExpr(scip, termprod, &termprodsimplified, &changed, &infeasible, ownercreate, ownercreatedata) );
            assert(!infeasible);  /* simplify of products should never be infeasible */
            SCIP_CALL( SCIPreleaseExpr(scip, &termprod) );
         }
         else
            termprodsimplified = termprod;

         /* append to product */
         SCIP_CALL( SCIPappendExprSumExpr(scip, *product, termprodsimplified, coef1 * coef2) );
         SCIP_CALL( SCIPreleaseExpr(scip, &termprodsimplified) );
      }
   }

   if( simplify )
   {
      SCIP_EXPR* prodsimplified;
      SCIP_Bool changed;
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPsimplifyExpr(scip, *product, &prodsimplified, &changed, &infeasible, ownercreate, ownercreatedata) );
      assert(!infeasible);
      SCIP_CALL( SCIPreleaseExpr(scip, product) );
      *product = prodsimplified;
   }

   return SCIP_OKAY;
}

/** constructs the expanded power of a sum expression
 *
 * @attention The number of terms in the expansion grows exponential with the exponent. Be aware of what you wish for.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPpowerExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR**           result,             /**< buffer where to store expanded power of sum */
   SCIP_EXPR*            base,               /**< sum */
   int                   exponent,           /**< exponent > 1 */
   SCIP_Bool             simplify,           /**< whether to simplify created terms and sum */
   SCIP_DECL_EXPR_OWNERCREATE((*ownercreate)), /**< function to call to create ownerdata */
   void*                 ownercreatedata     /**< data to pass to ownercreate */
   )
{
   /* for sum_i alpha_i expr_i  and exponent p, this constructs
   *   sum_{beta} (p over beta) prod_i (alpha_i expr_i)^beta_i   over all multiindex beta such that sum_i beta_i = p
   * See also https://en.wikipedia.org/wiki/Multinomial_theorem
   * Calculation of multinomials adapted from https://github.com/m-j-w/MultinomialSeries.jl
   *
   * TODO might be nice to integrate more with simplify and construct a simplified sum right away
   */

   SCIP_EXPR** children;
   SCIP_EXPR*** childrenpower;
   SCIP_Real* coefs;
   SCIP_Real constant;
   SCIP_Bool haveconst;
   int nchildren;
   int nterms;
   int* factorials;
   int* beta;
   int betapos;
   int restsum;
   int multinomialcoef;
   int i;

   SCIP_EXPR* newterm;
   SCIP_Real newtermcoef;

   assert(scip != NULL);
   assert(result != NULL);
   assert(base != NULL);
   assert(exponent > 1);
   assert(SCIPisExprSum(scip, base));
   assert(exponent > 1.0);

   children = SCIPexprGetChildren(base);
   nchildren = SCIPexprGetNChildren(base);
   coefs = SCIPgetCoefsExprSum(base);
   constant = SCIPgetConstantExprSum(base);
   haveconst = constant != 0.0;
   nterms = nchildren + (haveconst ? 1 : 0);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "expanding (");
   SCIPprintExpr(scip, base, NULL);
   SCIPinfoMessage(scip, NULL, ")^%d\n", exponent);
#endif

   SCIP_CALL( SCIPcreateExprSum(scip, result, 0, NULL, NULL, 0.0, ownercreate, ownercreatedata) );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &beta, nterms) );
   SCIP_CALL( SCIPallocBufferArray(scip, &factorials, exponent+1) );

   /* precompute factorials k!, k = 0...exponent */
   factorials[0] = 1;
   for( i = 1; i <= exponent; ++i )
     factorials[i] = factorials[i-1] * i;

   /* precompute children^k, k=1..exponent */
   SCIP_CALL( SCIPallocBufferArray(scip, &childrenpower, nchildren) );
   for( i = 0; i < nchildren; ++i )
   {
      int k;
      SCIP_CALL( SCIPallocBufferArray(scip, &childrenpower[i], exponent+1) );
      childrenpower[i][1] = children[i];
      for( k = 2; k <= exponent; ++k )
      {
         SCIP_CALL( SCIPcreateExprPow(scip, &childrenpower[i][k], children[i], (SCIP_Real)k, NULL, NULL) );
         if( simplify )
         {
            SCIP_Bool changed;
            SCIP_Bool infeasible;
            SCIP_EXPR* simplified;

            SCIP_CALL( SCIPsimplifyExpr(scip, childrenpower[i][k], &simplified, &changed, &infeasible, ownercreate, ownercreatedata) );
            assert(!infeasible);
            SCIP_CALL( SCIPreleaseExpr(scip, &childrenpower[i][k]) );
            childrenpower[i][k] = simplified;
         }
      }
   }

   /* first multinomial is (exponent,0,0,...) */
   beta[0] = exponent;
   betapos = 0;
   do
   {
      /* compute multinomial coef exponent! /  (beta[0]! * ... * beta[nterms-1]!) */
      multinomialcoef = factorials[exponent];
      for( i = 0; i < nterms; ++i )
      {
         assert(beta[i] >= 0);
         assert(beta[i] <= exponent);
         multinomialcoef /= factorials[beta[i]];
      }

      SCIPdebugMsg(scip, "multinomial (");
      for( i = 0; i < nterms; ++i )
         SCIPdebugPrintf("%d ", beta[i]);
      SCIPdebugPrintf(")  with coef %d\n", multinomialcoef);

      /* construct new term for expanded sum */
      SCIP_CALL( SCIPcreateExprProduct(scip, &newterm, 0, NULL, 1.0, ownercreate, ownercreatedata) );
      newtermcoef = multinomialcoef;
      for( i = 0; i < nterms; ++i )
      {
         if( beta[i] == 0 )
            continue;

         if( i == nterms-1 && haveconst )
         {
            /* if constant term, then update newtermcoef */
            newtermcoef *= pow(constant, (double)beta[i]);
            continue;
         }

         /* alpha_i^beta_i */
         newtermcoef *= pow(coefs[i], (double)beta[i]);

         /* expr_i^beta_i*/
         SCIP_CALL( SCIPappendExprChild(scip, newterm, childrenpower[i][beta[i]]) );
      }

      /* append newterm to sum */
      switch( SCIPexprGetNChildren(newterm) )
      {
         case 0:
         {
            /* no factor in product, so it is a constant -> update constant in sum */
            SCIPsetConstantExprSum(*result, SCIPgetConstantExprSum(*result) + newtermcoef);
            break;
         }

         case 1:
         {
            /* only one factor in product -> add this factor itself to sum */
            SCIP_CALL( SCIPappendExprSumExpr(scip, *result, SCIPexprGetChildren(newterm)[0], newtermcoef) );
            break;
         }

         default:
         {
            if( simplify )
            {
               /* simplify product */
               SCIP_Bool changed;
               SCIP_Bool infeasible;
               SCIP_EXPR* simplified;

               SCIP_CALL( SCIPsimplifyExpr(scip, newterm, &simplified, &changed, &infeasible, ownercreate, ownercreatedata) );
               assert(!infeasible);
               SCIP_CALL( SCIPreleaseExpr(scip, &newterm) );
               newterm = simplified;
            }

            /* append new term to sum */
            SCIP_CALL( SCIPappendExprSumExpr(scip, *result, newterm, newtermcoef) );
            break;
         }
      }
      SCIP_CALL( SCIPreleaseExpr(scip, &newterm) );

      /* determine next beta */
      while( beta[betapos] == 0 )
         --betapos;

      if( betapos == nterms-1 )
      {
         do
         {
            if( betapos == 0 )
               goto TERMINATE;
            --betapos;
         }
         while( beta[betapos] == 0 );

         restsum = 0;
         for( i = betapos+1; i < nterms; ++i )
         {
            restsum += beta[i];
            beta[i] = 0;
         }
         beta[betapos+1] = restsum;
      }

      if( beta[betapos] > 0 )
      {
         --beta[betapos];
         ++beta[++betapos];
      }
   }
   while( TRUE );  /*lint !e506*/

 TERMINATE:

   if( simplify )
   {
      SCIP_Bool changed;
      SCIP_Bool infeasible;
      SCIP_EXPR* simplified;

      SCIP_CALL( SCIPsimplifyExpr(scip, *result, &simplified, &changed, &infeasible, ownercreate, ownercreatedata) );
      assert(!infeasible);
      SCIP_CALL( SCIPreleaseExpr(scip, result) );
      *result = simplified;
   }

   for( i = nchildren-1; i >= 0; --i )
   {
      int k;
      for( k = exponent; k >= 2; --k )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &childrenpower[i][k]) );
      }
      SCIPfreeBufferArray(scip, &childrenpower[i]);
   }
   SCIPfreeBufferArray(scip, &childrenpower);
   SCIPfreeBufferArray(scip, &factorials);
   SCIPfreeBufferArray(scip, &beta);

#ifdef SCIP_DEBUG
   SCIPinfoMessage(scip, NULL, "-> ");
   SCIPprintExpr(scip, *result, NULL);
   SCIPinfoMessage(scip, NULL, "\n");
#endif

   return SCIP_OKAY;
}

/* from pub_expr.h */

/** gets the coefficients of a summation expression */
SCIP_Real* SCIPgetCoefsExprSum(
   SCIP_EXPR*            expr                /**< sum expression */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   return exprdata->coefficients;
}

/** gets the constant of a summation expression */
SCIP_Real SCIPgetConstantExprSum(
   SCIP_EXPR*            expr                /**< sum expression */
   )
{
   SCIP_EXPRDATA* exprdata;

   assert(expr != NULL);

   exprdata = SCIPexprGetData(expr);
   assert(exprdata != NULL);

   return exprdata->constant;
}
