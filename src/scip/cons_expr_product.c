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

/**@file   cons_expr_product.c
 * @brief  product expression handler
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 *
 * Implementation of the product expression, representing a product of expressions
 * and a constant, i.e., coef * prod_i x_i.
 *
 * @todo initsepaProduct
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "scip/cons_expr_product.h"
#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_pow.h"

#define PRODUCT_PRECEDENCE  50000
#define PRODUCT_HASHKEY     SCIPcalcFibHash(54949.0)

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
   SCIP_Real* coefficients; /**< coefficients / exponents of children */
   int        coefssize;    /**< size of the coefficients array */

   SCIP_ROW*  row;          /**< row created during initLP() */
};

/** node for linked list of expressions */
struct exprnode
{
   SCIP_CONSEXPR_EXPR*   expr;               /**< expression in node */
   SCIP_Real             coef;               /**< coefficient or exponent of expr*/
   struct exprnode*      next;               /**< next node */
};

typedef struct exprnode EXPRNODE;


/*
 * Local methods
 */

/*  methods for handling linked list of expressions */
/** inserts newnode at beginning of list */
static
void insertFirstList(
   EXPRNODE*             newnode,            /**< node to insert */
   EXPRNODE**            list                /**< list */
   )
{
   assert(list != NULL);
   assert(newnode != NULL);

   newnode->next = *list;
   *list = newnode;
}

/** removes first element of list and returns it */
static
EXPRNODE* listPopFirst(
   EXPRNODE**            list                /**< list */
   )
{
   EXPRNODE* first;

   assert(list != NULL);

   if( *list == NULL )
      return NULL;

   first = *list;
   *list = (*list)->next;
   first->next = NULL;

   return first;
}

/** returns length of list */
static
int listLength(
   EXPRNODE*             list                /**< list */
   )
{
   int length;

   if( list == NULL )
      return 0;

   length = 1;
   while( (list=list->next) != NULL )
      ++length;

   return length;
}

/** creates expression node and capture expression */
static
SCIP_RETCODE createExprNode(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression stored at node */
   SCIP_Real             coef,               /**< coefficient/exponent of expression */
   EXPRNODE**            newnode             /**< pointer to store node */
   )
{
   SCIP_CALL( SCIPallocBuffer(scip, newnode) );

   (*newnode)->expr = expr;
   (*newnode)->coef = coef;
   (*newnode)->next = NULL;
   SCIPcaptureConsExprExpr(expr);

   return SCIP_OKAY;
}

/** creates expression list from expressions */
static
SCIP_RETCODE createExprlistFromExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  exprs,              /**< expressions stored in list */
   SCIP_Real*            coefs,              /**< coefficients/exponents of expression */
   SCIP_Real             coef,               /**< coefficient/exponent to multiply coefs (distributing) */
   int                   nexprs,             /**< number of expressions */
   EXPRNODE**            list                /**< pointer to store list */
   )
{
   int i;

   assert(*list == NULL);
   assert(nexprs > 0);
   assert(coef != 0);

   debugSimplify("building expr list from %d expressions\n", nexprs); /*lint !e506 !e681*/
   for( i = nexprs - 1; i >= 0; --i )
   {
      EXPRNODE* newnode;

      SCIP_CALL( createExprNode(scip, exprs[i], coef*(coefs ? coefs[i] : 1.0), &newnode) );
      insertFirstList(newnode, list);
   }
   assert(nexprs > 1 || (*list)->next == NULL);

   return SCIP_OKAY;
}

/** frees expression node and release expressions */
static
SCIP_RETCODE freeExprNode(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE**            node                /**< node to be freed */
   )
{
   assert(node != NULL && *node != NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*node)->expr) );
   SCIPfreeBuffer(scip, node);

   return SCIP_OKAY;
}

/** frees an expression list */
static
SCIP_RETCODE freeExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE**            exprlist            /**< list */
   )
{
   EXPRNODE* current;

   if( *exprlist == NULL )
      return SCIP_OKAY;

   current = *exprlist;
   while( current != NULL )
   {
      EXPRNODE* tofree;

      tofree = current;
      current = current->next;
      SCIP_CALL( freeExprNode(scip, &tofree) );
   }
   assert(current == NULL);
   *exprlist = NULL;

   return SCIP_OKAY;
}

/* helper functions for simplifying expressions */

/** merges tomerge into finalchildren
 *
 * Both, tomerge and finalchildren contain expressions that could be the children of a simplified sum
 * (except for SS6 and SS7 which are enforced later).
 * However, the concatenation of both lists, will not in general yield a simplified sum expression,
 * because both SS4 and SS5 could be violated. So the purpose of this method is to enforce SS4 and SS5.
 * In the process of enforcing SS4, it could happen that SS8 is violated, but this is easy to fix.
 * note: if tomerge has more than one element, then they are the children of a simplified sum expression
 * so no values, nor sum expressions, but products, variable or function expressions
 */
static
SCIP_RETCODE mergeSumExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE*             tomerge,            /**< list to merge */
   EXPRNODE**            finalchildren       /**< pointer to store the result of merge between tomerge and *finalchildren */
   )
{
   EXPRNODE* tomergenode;
   EXPRNODE* current;
   EXPRNODE* previous;

   if( tomerge == NULL )
      return SCIP_OKAY;

   if( *finalchildren == NULL )
   {
      *finalchildren = tomerge;
      return SCIP_OKAY;
   }

   tomergenode = tomerge;
   current = *finalchildren;
   previous = NULL;

   while( tomergenode != NULL && current != NULL )
   {
      int compareres;
      EXPRNODE* aux;

      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(tomergenode->expr)), "sum") != 0);
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(tomergenode->expr)), "val") != 0);
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(current->expr)), "sum") != 0);
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(current->expr)), "val") != 0);
      assert(previous == NULL || previous->next == current);

      compareres = SCIPcompareConsExprExprs(current->expr, tomergenode->expr);

      debugSimplify("comparing exprs:\n"); /*lint !e506 !e681*/
      #ifdef SIMPLIFY_DEBUG
      SCIP_CALL( SCIPprintConsExprExpr(scip, current->expr, NULL) );
      SCIPinfoMessage(scip, NULL, " vs ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, tomergenode->expr, NULL) );
      SCIPinfoMessage(scip, NULL, ": won %d\n", compareres);
      #endif

      if( compareres == 0 )
      {
         /* enforces SS4 and SS8 */
         current->coef += tomergenode->coef;

         /* destroy tomergenode (since will not use the node again) */
         aux = tomergenode;
         tomergenode = tomergenode->next;
         SCIP_CALL( freeExprNode(scip, &aux) );

         /* if coef is 0, remove term: if current is the first node, pop; if not, use previous and current to remove */
         if( current->coef == 0.0 )
         {
            debugSimplify("GOT 0 WHILE ADDING UP\n"); /*lint !e506 !e681*/
            if( current == *finalchildren )
            {
               assert(previous == NULL);
               aux = listPopFirst(finalchildren);
               assert(aux == current);
               current = *finalchildren;
            }
            else
            {
               assert(previous != NULL);
               aux = current;
               current = current->next;
               previous->next = current;
            }

            SCIP_CALL( freeExprNode(scip, &aux) );
         }
      }
      /* enforces SS5 */
      else if( compareres == -1 )
      {
         /* current < tomergenode => move current */
         previous = current;
         current = current->next;
      }
      else
      {
         assert(compareres == 1);

         /* insert: if current is the first node, then insert at beginning; otherwise, insert between previous and current */
         if( current == *finalchildren )
         {
            assert(previous == NULL);
            aux = tomergenode;
            tomergenode = tomergenode->next;
            insertFirstList(aux, finalchildren);
            previous = *finalchildren;
         }
         else
         {
            assert(previous != NULL);
            /* extract */
            aux = tomergenode;
            tomergenode = tomergenode->next;
            /* insert */
            previous->next = aux;
            aux->next = current;
            previous = aux;
         }
      }
   }

   /* if all nodes of tomerge were merged, we are done */
   if( tomergenode == NULL )
      return SCIP_OKAY;

   assert(current == NULL);

   /* if all nodes of finalchildren were cancelled by nodes of tomerge, then the rest of tomerge is finalchildren */
   if( *finalchildren == NULL )
   {
      assert(previous == NULL);
      *finalchildren = tomergenode;
      return SCIP_OKAY;
   }

   /* there are still nodes of tomerge unmerged; these nodes are larger than finalchildren, so append at end */
   assert(previous != NULL && previous->next == NULL);
   previous->next = tomergenode;

   return SCIP_OKAY;
}

/** creates a sum expression with the elements of exprlist as its children */
static
SCIP_RETCODE createExprSumFromExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE*             exprlist,           /**< list containing the children of expr */
   SCIP_Real             constant,           /**< constant of expr */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to store the sum expression */
   )
{
   int i;
   int nchildren;
   SCIP_CONSEXPR_EXPR** children;
   SCIP_Real* coefs;

   nchildren = listLength(exprlist);

   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      children[i] = exprlist->expr;
      coefs[i] = exprlist->coef;
      exprlist = exprlist->next;
   }
   assert(exprlist == NULL);

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, SCIPfindConshdlr(scip, "expr"), expr, nchildren, children, coefs, constant) );

   SCIPfreeBufferArray(scip, &children);
   SCIPfreeBufferArray(scip, &coefs);

   return SCIP_OKAY;
}

/** creates a product expression with the elements of exprlist as its children */
static
SCIP_RETCODE createExprProductFromExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE*             exprlist,           /**< list containing the children of expr */
   SCIP_Real             coef,               /**< coef of expr */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to store the product expression */
   )
{
   int i;
   int nchildren;
   SCIP_CONSEXPR_EXPR** children;

   /* asserts SP8 */
   assert(coef == 1.0);
   nchildren = listLength(exprlist);

   SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      children[i] = exprlist->expr;
      assert(exprlist->coef == 1.0);
      exprlist = exprlist->next;
   }

   assert(exprlist == NULL);

   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, SCIPfindConshdlr(scip, "expr"), expr, nchildren, children, coef) );

   SCIPfreeBufferArray(scip, &children);

   return SCIP_OKAY;
}

/** simplifies a term of a sum expression: constant * expr, so that it is a valid child of a simplified sum expr.
 * @note: in contrast to other simplify methods, this does *not* return a simplified expression.
 * Instead, the method is intended to be called only when simplifying a sum expression,
 * Since in general, constant*expr is not a simplified child of a sum expression, this method returns
 * a list of expressions L, such that (sum L) = constant * expr *and* each expression in L
 * is a valid child of a simplified sum expression.
 */
static
SCIP_RETCODE simplifyTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be simplified */
   SCIP_Real             coef,               /**< coefficient of expressions to be simplified */
   SCIP_Real*            simplifiedconstant, /**< constant of parent sum expression */
   EXPRNODE**            simplifiedterm      /**< pointer to store the resulting expression node/list of nodes */
   )
{
   const char* exprtype;

   assert(simplifiedterm != NULL);
   assert(*simplifiedterm == NULL);
   assert(expr != NULL);

   exprtype = SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr));

   /* enforces SS3 */
   if( strcmp(exprtype, "val") == 0 )
   {
      *simplifiedconstant += coef * SCIPgetConsExprExprValueValue(expr);
      return SCIP_OKAY;
   }

   /* enforces SS2
    * We do not need to modify `expr`: we still need to distribute `coef` over `expr`
    * and this operation can still render `expr` unsimplified, e.g., (sum 0 2 (sum 0 1/2 x)) -> (sum 0 1 (sum 0 1 x)),
    * which is unsimplified. However, this is the only case. To see this, notice that we can regard `expr` as a sum
    * with constant 0 (because the constant will be passed to the parent), so `expr` = (sum 0 coef1 expr1 coef2 expr2...)
    * and after distributing `coef`, expr' = (sum coef1' expr1 coef2' expr2...) which will clearly satisfy SS1-SS4, SS6 and
    * SS8. SS5 is satisfied, because if coef1 expr1 < coef2 expr2 are children in a simplified sum, then expr1 != expr2.
    * Therefore expr1 < expr2, which implies that C1 * expr1 < C2 * expr2 for any C1, C2 different from 0
    * So the only condition that can fail is SS7. In that case, expr = (sum coef1 expr1) and expr' = (sum 1 expr1)
    * and so simplifying expr' gives expr1.
    * All this can be done and checked without modifying expr
    */
   if( strcmp(exprtype, "sum") == 0 )
   {
      /* pass constant to parent */
      *simplifiedconstant += coef * SCIPgetConsExprExprSumConstant(expr);

      /* check if SS7 could fail after distributing */
      if( SCIPgetConsExprExprNChildren(expr) == 1 && coef * SCIPgetConsExprExprSumCoefs(expr)[0] == 1.0 )
      {
         /* this is a one level recursion: expr is sum -> no child of expr is sum */
         assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(expr)[0])), "sum") != 0);
         SCIP_CALL( simplifyTerm(scip, SCIPgetConsExprExprChildren(expr)[0], 1.0, simplifiedconstant, simplifiedterm) );
         return SCIP_OKAY;
      }

      /* distribute */
      SCIP_CALL( createExprlistFromExprs(scip, SCIPgetConsExprExprChildren(expr),
               SCIPgetConsExprExprSumCoefs(expr), coef, SCIPgetConsExprExprNChildren(expr), simplifiedterm) );
   }
   else
   {
      /* other types of (simplified) expressions can be a child of a simplified sum */
      assert(strcmp(exprtype, "sum") != 0);
      assert(strcmp(exprtype, "val") != 0);

      SCIP_CALL( createExprNode(scip, expr, coef, simplifiedterm) );
   }

   return SCIP_OKAY;
}

/* signatures are needed by simplifyFactor */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyProduct);
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifySum);
/** simplifies a factor of a product expression: (base)^exponent, so that it is a valid children of a simplified product expr
 * @note: in contrast to other simplify methods, this does *not* return a simplified expression.
 * Instead, the method is intended to be called only when simplifying a product expression,
 * Since in general, (base)^exponent is not a simplified child of a product expression, this method returns
 * a list of expressions (with exponents) L, such that (prod L) = (base)^exponent *and* each expression in L
 * is a valid child of a simplified product expression.
 * TODO: handle more cases when exponent is non-integer and base >= 0 (of course one has to adapt the definition
 * of a simplified product)
 */
static
SCIP_RETCODE simplifyFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   base,               /**< expression to be simplified */
   SCIP_Real             exponent,           /**< exponent of expressions to be simplified */
   SCIP_Real*            simplifiedcoef,     /**< coefficient of parent product expression */
   EXPRNODE**            simplifiedfactor    /**< pointer to store the resulting expression node/list of nodes */
   )
{
   const char* basetype;

   assert(simplifiedfactor != NULL);
   assert(*simplifiedfactor == NULL);
   assert(base != NULL);

   basetype = SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(base));

   /* enforces SP7 */
   if( strcmp(basetype, "val") == 0 )
   {
      /* TODO: if val < 0 and exponent non integer -> domain error/undefined etc */
      debugSimplify("[simplifyFactor] seeing value %g, exponent %g -> include in coef\n", SCIPgetConsExprExprValueValue(base), exponent); /*lint !e506 !e681*/
      *simplifiedcoef *= pow(SCIPgetConsExprExprValueValue(base), exponent);
      return SCIP_OKAY;
   }

   /* enforces SP2 */
   if( strcmp(basetype, "prod") == 0 )
   {
      assert(SCIPgetConsExprExprProductCoef(base) == 1.0);
      debugSimplify("[simplifyPow] seeing a product: include its children\n"); /*lint !e506 !e681*/

      SCIP_CALL( createExprlistFromExprs(scip, SCIPgetConsExprExprChildren(base),
               NULL, 1.0, SCIPgetConsExprExprNChildren(base), simplifiedfactor) );

      return SCIP_OKAY;
   }

   /* the given (simplified) expressions can be children of a simplified product */
   assert(strcmp(basetype, "prod") != 0);
   assert(strcmp(basetype, "val") != 0);
   SCIP_CALL( createExprNode(scip, base, exponent, simplifiedfactor) );

   return SCIP_OKAY;
}

/** merges tomerge into finalchildren
 * Both, tomerge and finalchildren contain expressions that could be the children of a simplified product
 * (except for SP8-SP10 which are enforced later).
 * However, the concatenation of both lists will not in general yield a simplified product expression,
 * because both SP4 and SP5 could be violated. So the purpose of this method is to enforce SP4 and SP5.
 * In the process of enforcing SP4, it could happen that SP2 and SP3 get violated. Since enforcing those
 * could in principle generate further violations, we remove the affected children from finalchildren
 * and include them in unsimplifiedchildren for further processing.
 * note: if tomerge has more than one element, then they are the children of a simplified product expression;
 * it can contain products, but only because they are acting as non-integer powers!
 * TODO: this function and mergeSumExprlist are very similar... merge them?
 */
static
SCIP_RETCODE mergeProductExprlist(
   SCIP*                 scip,               /**< SCIP data structure */
   EXPRNODE*             tomerge,            /**< list to merge */
   EXPRNODE**            finalchildren,      /**< pointer to store the result of merge between tomerge and *finalchildren */
   EXPRNODE**            unsimplifiedchildren/**< the list of children that should go to the product expression; they are
                                                  unsimplified when seen as children of a simplified product */
   )
{
   EXPRNODE* tomergenode;
   EXPRNODE* current;
   EXPRNODE* previous;

   if( tomerge == NULL )
      return SCIP_OKAY;

   if( *finalchildren == NULL )
   {
      *finalchildren = tomerge;
      return SCIP_OKAY;
   }

   tomergenode = tomerge;
   current = *finalchildren;
   previous = NULL;

   while( tomergenode != NULL && current != NULL )
   {
      int compareres;
      EXPRNODE* aux;
      SCIP_CONSEXPR_EXPR* base1;
      SCIP_CONSEXPR_EXPR* base2;
      SCIP_Real expo1;
      SCIP_Real expo2;

      /* assert invariants */
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(tomergenode->expr)), "val") != 0);
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(current->expr)), "val") != 0);
      assert(previous == NULL || previous->next == current);

      /* in general the base of an expression is itself if type(expr) != pow, otherwise it is child of pow */
      /* TODO: better documentation
       *       clean code */
      if( strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(current->expr)), "pow") == 0 )
      {
         base1 = SCIPgetConsExprExprChildren(current->expr)[0];
         expo1 = SCIPgetConsExprExprPowExponent(current->expr);
      }
      else
      {
         base1 = current->expr;
         expo1 = 1.0;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(tomergenode->expr)), "pow") == 0 )
      {
         base2 = SCIPgetConsExprExprChildren(tomergenode->expr)[0];
         expo2 = SCIPgetConsExprExprPowExponent(tomergenode->expr);
      }
      else
      {
         base2 = tomergenode->expr;
         expo2 = 1.0;
      }

      /* if both bases are the same: have to build simplifiy(base^(expo1 + expo2)) */
      if( SCIPcompareConsExprExprs(base1, base2) == 0 )
      {
         SCIP_CONSEXPR_EXPR* power;
         SCIP_CONSEXPR_EXPR* simplifiedpower;

         SCIP_CALL( SCIPcreateConsExprExprPow(scip, SCIPfindConshdlr(scip, "expr"), &power, base1, expo1 + expo2) );
         SCIP_CALL( SCIPsimplifyConsExprExpr(scip, power, &simplifiedpower) ); /* FIXME: call simplifyPow */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &power) );

         /* replace tomergenode's expression with simplifiedpower */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &tomergenode->expr) );
         tomergenode->expr = simplifiedpower;
         /* move tomergenode to unsimplifiedchildren */
         aux = tomergenode;
         tomergenode = tomergenode->next;
         insertFirstList(aux, unsimplifiedchildren);

         /* destroy current */
         if( current == *finalchildren )
         {
            assert(previous == NULL);
            aux = listPopFirst(finalchildren);
            assert(aux == current);
            current = *finalchildren;
         }
         else
         {
            assert(previous != NULL);
            aux = current;
            current = current->next;
            previous->next = current;
         }
         SCIP_CALL( freeExprNode(scip, &aux) );

         /* continue */
         continue;
      }

      /* bases are not the same, then expressions cannot be the same */
      compareres = SCIPcompareConsExprExprs(current->expr, tomergenode->expr);
      if( compareres == -1 )
      {
         /* current < tomergenode => move current */
         previous = current;
         current = current->next;
      }
      else
      {
         assert(compareres == 1);

         /* insert: if current is the first node, then insert at beginning; otherwise, insert between previous and current */
         if( current == *finalchildren )
         {
            assert(previous == NULL);
            aux = tomergenode;
            tomergenode = tomergenode->next;
            insertFirstList(aux, finalchildren);
            previous = *finalchildren;
         }
         else
         {
            assert(previous != NULL);
            /* extract */
            aux = tomergenode;
            tomergenode = tomergenode->next;
            /* insert */
            previous->next = aux;
            aux->next = current;
            previous = aux;
         }
      }
   }

   /* if all nodes of tomerge were merged, we are done */
   if( tomergenode == NULL )
      return SCIP_OKAY;

   assert(current == NULL);

   /* if all nodes of finalchildren were cancelled by nodes of tomerge (ie, transfered to unsimplifiedchildren),
    * then the rest of tomerge is finalchildren */
   if( *finalchildren == NULL )
   {
      assert(previous == NULL);
      *finalchildren = tomergenode;
      return SCIP_OKAY;
   }

   /* there are still nodes of tomerge unmerged; these nodes are larger than finalchildren, so append at end */
   assert(previous != NULL && previous->next == NULL);
   previous->next = tomergenode;

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
   (*exprdata)->row          = NULL;

   return SCIP_OKAY;
}

/*
 * Callback methods of expression handler
 */

/** simplifies a sum expression
 *
 * Summary: we first build a list of expressions (called finalchildren) which will be the children of the simplified sum
 * and then we process this list in order to enforce SS6 and SS7.
 * Description: To build finalchildren, each child of sum is manipulated (see simplifyTerm) in order to satisfy
 * SS2, SS3 and SS8 as follows
 * SS8: if the child's coefficient is 0, ignore it
 * SS3: if the child is a value, add the value to the sum's constant
 * SS2: if the child is a sum, we distribution that child's coefficient to its children and then build a list with the
 *      child's children. Note that distributing will not render the child unsimplified.
 * Otherwise (if it satisfies SS2, SS3 and SS8) we build a list with that child.
 * Then, we merge the built list into finalchildren (see mergeSumExprlist).
 * After finalchildren is done, we build the simplified sum expression out of it, taking care that SS6 and SS7 are satisfied
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifySum)
{
   SCIP_CONSEXPR_EXPR** children;
   EXPRNODE* finalchildren;
   SCIP_Real simplifiedconstant;
   SCIP_Real* coefs;
   int i;
   int nchildren;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "sum") == 0);

   children  = SCIPgetConsExprExprChildren(expr);
   nchildren = SCIPgetConsExprExprNChildren(expr);
   coefs     = SCIPgetConsExprExprSumCoefs(expr);

   /* while there are still children to process */
   finalchildren  = NULL;
   simplifiedconstant = SCIPgetConsExprExprSumConstant(expr);
   for( i = 0; i < nchildren; i++ )
   {
      EXPRNODE* tomerge;

      /* enforces SS8 */
      if( coefs[i] == 0.0 )
         continue;

      /* enforces SS2 and SS3 */
      tomerge = NULL;
      SCIP_CALL( simplifyTerm(scip, children[i], coefs[i], &simplifiedconstant, &tomerge) );

      /* enforces SS4 and SS5
       * note: merge frees (or uses) the nodes of the list tomerge */
      SCIP_CALL( mergeSumExprlist(scip, tomerge, &finalchildren) );
   }

   /* build sum expression from finalchildren and post-simplify */
   debugSimplify("what to do? finalchildren = %p has length %d\n", (void *)finalchildren, listLength(finalchildren)); /*lint !e506 !e681*/

   /* enforces SS6: if list is empty, return value */
   if( finalchildren == NULL )
   {
      debugSimplify("[sum] got empty list, return value %g\n", simplifiedconstant); /*lint !e506 !e681*/
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, SCIPfindConshdlr(scip, "expr"), simplifiedexpr, simplifiedconstant) );
   }
   /* enforces SS7
    * if list contains one expr with coef 1.0 and constant is 0, return that expr */
   else if( finalchildren->next == NULL && finalchildren->coef == 1.0 && simplifiedconstant == 0.0 )
   {
      *simplifiedexpr = finalchildren->expr;
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }
   /* build sum expression from list */
   else
   {
      SCIP_CALL( createExprSumFromExprlist(scip, finalchildren, simplifiedconstant, simplifiedexpr) );
   }

   /* free memory */
   SCIP_CALL( freeExprlist(scip, &finalchildren) );
   assert(finalchildren == NULL);

   return SCIP_OKAY;
}

/** simplifies a product expression
 *
 * Summary: we first build a list of expressions (called finalchildren) which will be the children of the simplified product
 * and then we process this list in order to enforce SP8-10
 * Description: In order to build finalchildren, we first build list of unsimplified children (called unsimplifiedchildren)
 * with the children of the product. Each node of the list is manipulated (see simplifyFactor) in order to satisfy
 * SP2, SP3, SP6 and SP7 as follows
 * SP7: if the node's expression is a value, multiply the value^exponent to the products's coef
 * SP6: if the node's exponent is 0, ignore it
 * SP2: if the node's expression is a product and its exponent is integer != 1, distribute exponent among its children
 *      and simplify again. If its exponent equals 1, then build a list with the child's children
 * SP3: if the node's expression is a sum with constant 0 with a unique child, multiply child's coef^exponent to the
 *      products coef and consider the sum's child as the child of product (simplify again)
 * Then, we merge the built list (or the simplified node) into finalchildren. While merging, nodes from finalchildren can
 * go back to unsimplifiedchildren for further processing (see mergeProductExprlist for more details)
 * After finalchildren is done, we build the simplified product expression out of it, taking care that SP8-10 are satisfied
 */
static
SCIP_DECL_CONSEXPR_EXPRSIMPLIFY(simplifyProduct)
{
   EXPRNODE* unsimplifiedchildren;
   EXPRNODE* finalchildren;
   SCIP_Real simplifiedcoef;

   assert(expr != NULL);
   assert(simplifiedexpr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "prod") == 0);

   /* set up list of current children (when looking at each of them individually, they are simplified, but as
    * children of a product expression they might be unsimplified) */
   unsimplifiedchildren = NULL;
   SCIP_CALL( createExprlistFromExprs(scip, SCIPgetConsExprExprChildren(expr),
            NULL, 1.0, SCIPgetConsExprExprNChildren(expr), &unsimplifiedchildren) );

   /* while there are still children to process */
   finalchildren  = NULL;
   simplifiedcoef = SCIPgetConsExprExprProductCoef(expr);
   while( unsimplifiedchildren != NULL )
   {
      EXPRNODE* tomerge;
      EXPRNODE* first;

      /* if the simplified coefficient is 0, we can return value 0
       * TODO: there might be some domain error in the unprocessed expressions. should we take care of this?
       */
      if( simplifiedcoef == 0.0 )
      {
         SCIP_CALL( freeExprlist(scip, &finalchildren) );
         SCIP_CALL( freeExprlist(scip, &unsimplifiedchildren) );
         assert(finalchildren == NULL);
         break;
      }

      first = listPopFirst(&unsimplifiedchildren);
      assert(first != NULL);

      /* enforces SP2 and SP7 */
      tomerge = NULL;
      assert(first->coef == 1.0);
      SCIP_CALL( simplifyFactor(scip, first->expr, first->coef, &simplifiedcoef, &tomerge) );

      /* enforces SP4 and SP5
       * note: merge frees (or uses) the nodes of the tomerge list */
      SCIP_CALL( mergeProductExprlist(scip, tomerge, &finalchildren, &unsimplifiedchildren) );

      /* free first */
      SCIP_CALL( freeExprlist(scip, &first) );
   }

   /* build product expression from finalchildren and post-simplify */
   debugSimplify("what to do? finalchildren = %p has length %d\n", (void *)finalchildren, listLength(finalchildren)); /*lint !e506 !e681*/

   /* enforces SP10: if list is empty, return value */
   if( finalchildren == NULL )
   {
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, SCIPfindConshdlr(scip, "expr"), simplifiedexpr, simplifiedcoef) );
   }
   /* enforces SP9: if finalchildren has only one expr with exponent 1.0 and coef 1.0, return that expr */
   else if( finalchildren->next == NULL && finalchildren->coef == 1.0 && simplifiedcoef == 1.0 )
   {
      *simplifiedexpr = finalchildren->expr;
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }
   /* enforces SP8 and SP9: if finalchildren has only one expr with exponent 1.0 and coef != 1.0, return (sum 0 coef expr) */
   else if( finalchildren->next == NULL && finalchildren->coef == 1.0 )
   {
      SCIP_CONSEXPR_EXPR* aux;

      SCIP_CALL( SCIPcreateConsExprExprSum(scip, SCIPfindConshdlr(scip, "expr"), &aux,
               1, &(finalchildren->expr), &simplifiedcoef, 0.0) );

      /* simplifying here is necessary, the product could have had sums as children
       * e.g., (prod 2 1 (sum 1 <x>)) -> (sum 0 2 (sum 1 <x>)) and that needs to be simplified
       */
      SCIP_CALL( simplifySum(scip, aux, simplifiedexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );
   }
   /* enforces SP8: if simplifiedcoef != 1.0, transform it into a sum with the (simplified) product as child */
   else if( simplifiedcoef != 1.0 )
   {
      SCIP_CONSEXPR_EXPR* aux;

      SCIP_CALL( createExprProductFromExprlist(scip, finalchildren, 1.0, &aux) );
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, SCIPfindConshdlr(scip, "expr"), simplifiedexpr,
               1, &aux, &simplifiedcoef, 0.0) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );
   }
   /* build product expression from list */
   else
   {
      SCIP_CALL( createExprProductFromExprlist(scip, finalchildren, simplifiedcoef, simplifiedexpr) );
   }

   /* free memory */
   SCIP_CALL( freeExprlist(scip, &finalchildren) );
   assert(finalchildren == NULL);

   assert(*simplifiedexpr != NULL);
   return SCIP_OKAY;
}

/** the order of two product expressions is a lexicographical order on the terms. See compareSum
 *  Example: y * z < x * y * z, x^2 * y^3 < x^3 * y^3
 */
static
SCIP_DECL_CONSEXPR_EXPRCMP(compareProduct)
{  /*lint --e{715}*/
   int compareresult;
   int i;
   int j;
   int nchildren1;
   int nchildren2;
   SCIP_CONSEXPR_EXPR** children1;
   SCIP_CONSEXPR_EXPR** children2;

   nchildren1 = SCIPgetConsExprExprNChildren(expr1);
   nchildren2 = SCIPgetConsExprExprNChildren(expr2);
   children1 = SCIPgetConsExprExprChildren(expr1);
   children2 = SCIPgetConsExprExprChildren(expr2);

   for( i = nchildren1 - 1, j = nchildren2 - 1; i >= 0 && j >= 0; --i, --j )
   {
      compareresult = SCIPcompareConsExprExprs(children1[i], children2[j]);
      if( compareresult != 0 )
         return compareresult;
      /* expressions are equal, continue */
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

   /* everything is equal, use coefficient as tie-breaker */
   assert(i == -1 && j == -1);
   if( SCIPgetConsExprExprProductCoef(expr1) < SCIPgetConsExprExprProductCoef(expr2) )
      return -1;
   if( SCIPgetConsExprExprProductCoef(expr1) > SCIPgetConsExprExprProductCoef(expr2) )
      return 1;

   /* they are equal */
   return 0;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYHDLR(copyhdlrProduct)
{  /*lint --e{715}*/
   SCIP_CALL( SCIPincludeConsExprExprHdlrProduct(scip, consexprhdlr) );
   *valid = TRUE;

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA(copydataProduct)
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
SCIP_DECL_CONSEXPR_EXPRFREEDATA(freedataProduct)
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

/** sum and products hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *hashkey = PRODUCT_HASHKEY;
   *hashkey ^= SCIPcalcFibHash(exprdata->constant);

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      unsigned int childhash;

      assert(SCIPhashmapExists(expr2key, (void*)SCIPgetConsExprExprChildren(expr)[c]));
      childhash = SCIPcalcFibHash(exprdata->coefficients[c]) ^ (unsigned int)(size_t)SCIPhashmapGetImage(expr2key, SCIPgetConsExprExprChildren(expr)[c]);

      *hashkey ^= childhash;
   }

   return SCIP_OKAY;
}

/** expression point evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPREVAL(evalProduct)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_Real childval;
   SCIP_Real powval;
   int c;

   assert(expr != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *val = exprdata->constant;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr) && (*val != 0.0); ++c )
   {
      childval = SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[c]);
      assert(childval != SCIP_INVALID); /*lint !e777*/

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

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffProduct)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);
   assert(idx >= 0 && idx < SCIPgetConsExprExprNChildren(expr));

   child = SCIPgetConsExprExprChildren(expr)[idx];
   assert(child != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(child)), "val") != 0);
   assert(SCIPgetConsExprExprValue(child) != SCIP_INVALID);

   if( !SCIPisZero(scip, SCIPgetConsExprExprValue(child)) )
      *val = SCIPgetConsExprExprValue(expr) / SCIPgetConsExprExprValue(child);
   else
   {
      int i;

      *val = SCIPgetConsExprExprData(expr)->constant;
      for( i = 0; i < SCIPgetConsExprExprNChildren(expr) && (*val != 0.0); ++i )
      {
         if( i == idx )
            continue;

         *val *= SCIPgetConsExprExprValue(SCIPgetConsExprExprChildren(expr)[i]);
      }
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

/** helper function to separate a given point; needed for proper unittest */
static
SCIP_RETCODE separatePointProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_SOL*             sol,                /**< solution to be separated (NULL for the LP solution) */
   SCIP_ROW**            cut                 /**< pointer to store the row */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR* auxvar;
   SCIP_VAR* var;
   SCIP_Real violation;
   SCIP_Bool overestimate;
   SCIP_Bool success;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "prod") == 0);
   assert(cut != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);
   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);

   *cut = NULL;

   /* compute violation of the expression by evaluating auxiliary variables */
   violation = exprdata->constant;
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      child = SCIPgetConsExprExprChildren(expr)[c];

      /* value expressions should have been removed during simplification */
      assert(SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrValue(conshdlr));

      var = SCIPgetConsExprExprLinearizationVar(child);
      assert(var != NULL);

      violation *= pow(SCIPgetSolVal(scip, sol, var), exprdata->coefficients[c]);
   }
   violation -= SCIPgetSolVal(scip, sol, auxvar);

   /* no violation in this sub-expression */
   if( SCIPisEQ(scip, violation, 0.0) )
      return SCIP_OKAY;

   overestimate = SCIPisLT(scip, violation, 0.0);
   success = FALSE;

   /* square term */
   if( SCIPgetConsExprExprNChildren(expr) == 1 && exprdata->coefficients[0] == 2.0 )
   {
      SCIP_VAR* x;
      SCIP_Real lincoef;
      SCIP_Real linconstant;
      SCIP_Real refpoint;
      SCIP_Bool islocal;

      /* collect variable */
      child = SCIPgetConsExprExprChildren(expr)[0];
      x = SCIPgetConsExprExprLinearizationVar(child);
      assert(x != NULL);

      lincoef = 0.0;
      linconstant = 0.0;
      refpoint = SCIPgetSolVal(scip, sol, x);
      success = TRUE;

      /* adjust the reference points */
      refpoint = SCIPisLT(scip, refpoint, SCIPvarGetLbLocal(x)) ? SCIPvarGetLbLocal(x) : refpoint;
      refpoint = SCIPisGT(scip, refpoint, SCIPvarGetUbLocal(x)) ? SCIPvarGetUbLocal(x) : refpoint;
      assert(SCIPisLE(scip, refpoint, SCIPvarGetUbLocal(x)) && SCIPisGE(scip, refpoint, SCIPvarGetLbLocal(x)));

      /* decide whether to use linearization or secant */
      if( (exprdata->constant < 0 && overestimate) || (exprdata->constant > 0 && !overestimate) )
      {
         SCIPaddSquareLinearization(scip, exprdata->constant, refpoint, SCIPvarIsIntegral(x), &lincoef, &linconstant, &success);
         islocal = FALSE; /* linearizations are globally valid */
      }
      else
      {
         SCIPaddSquareSecant(scip, exprdata->constant, SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), refpoint, &lincoef, &linconstant, &success);
         islocal = TRUE; /* secants are only valid locally */
      }

      /* @todo  allow lhs/rhs of +/- infinity? */
      if( success && !SCIPisInfinity(scip, REALABS(linconstant)) )
      {
         SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, islocal ? "square_secant" : "square_linearization", 0, NULL, NULL, -SCIPinfinity(scip),
               SCIPinfinity(scip), islocal, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *cut, x, lincoef) );
         SCIP_CALL( SCIPaddVarToRow(scip, *cut, auxvar, -1.0) );

         if( overestimate )
         {
            SCIP_CALL( SCIPchgRowLhs(scip, *cut, -linconstant) );
         }
         else
         {
            SCIP_CALL( SCIPchgRowRhs(scip, *cut, -linconstant) );
         }
      }
   }
   /* bilinear term */
   else if( SCIPgetConsExprExprNChildren(expr) == 2 && exprdata->coefficients[0] == 1.0 && exprdata->coefficients[1] == 1.0 )
   {
      SCIP_VAR* x;
      SCIP_VAR* y;
      SCIP_Real lincoefx;
      SCIP_Real lincoefy;
      SCIP_Real linconstant;
      SCIP_Real refpointx;
      SCIP_Real refpointy;

      /* collect first variable */
      child = SCIPgetConsExprExprChildren(expr)[0];
      x = SCIPgetConsExprExprLinearizationVar(child);
      assert(x != NULL);

      /* collect second variable */
      child = SCIPgetConsExprExprChildren(expr)[1];
      y = SCIPgetConsExprExprLinearizationVar(child);
      assert(y != NULL);

      lincoefx = 0.0;
      lincoefy = 0.0;
      linconstant = 0.0;
      refpointx = SCIPgetSolVal(scip, sol, x);
      refpointy = SCIPgetSolVal(scip, sol, y);
      success = TRUE;

      /* adjust the reference points */
      refpointx = SCIPisLT(scip, refpointx, SCIPvarGetLbLocal(x)) ? SCIPvarGetLbLocal(x) : refpointx;
      refpointx = SCIPisGT(scip, refpointx, SCIPvarGetUbLocal(x)) ? SCIPvarGetUbLocal(x) : refpointx;
      refpointy = SCIPisLT(scip, refpointy, SCIPvarGetLbLocal(y)) ? SCIPvarGetLbLocal(y) : refpointy;
      refpointy = SCIPisGT(scip, refpointy, SCIPvarGetUbLocal(y)) ? SCIPvarGetUbLocal(y) : refpointy;
      assert(SCIPisLE(scip, refpointx, SCIPvarGetUbLocal(x)) && SCIPisGE(scip, refpointx, SCIPvarGetLbLocal(x)));
      assert(SCIPisLE(scip, refpointy, SCIPvarGetUbLocal(y)) && SCIPisGE(scip, refpointy, SCIPvarGetLbLocal(y)));

      SCIPaddBilinMcCormick(scip, exprdata->constant, SCIPvarGetLbLocal(x), SCIPvarGetUbLocal(x), refpointx,
         SCIPvarGetLbLocal(y), SCIPvarGetUbLocal(y), refpointy, overestimate, &lincoefx, &lincoefy, &linconstant,
         &success);

      /* @todo allow lhs/rhs of +/- infinity? */
      if( success && !SCIPisInfinity(scip, REALABS(linconstant)) )
      {
         /* McCormicks are only valid locally */
         SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, "mccormick", 0, NULL, NULL, -SCIPinfinity(scip),
               SCIPinfinity(scip), TRUE, FALSE, FALSE) );

         SCIP_CALL( SCIPaddVarToRow(scip, *cut, x, lincoefx) );
         SCIP_CALL( SCIPaddVarToRow(scip, *cut, y, lincoefy) );
         SCIP_CALL( SCIPaddVarToRow(scip, *cut, auxvar, -1.0) );

         if( overestimate )
         {
            SCIP_CALL( SCIPchgRowLhs(scip, *cut, -linconstant) );
         }
         else
         {
            SCIP_CALL( SCIPchgRowRhs(scip, *cut, -linconstant) );
         }
      }
   }
   else
   {
      /* @todo can not be handled so far */
   }

   return SCIP_OKAY;
}


/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaProduct)
{
   SCIP_ROW* cut;
   SCIP_Bool infeasible;

   *result = SCIP_DIDNOTFIND;
   *ncuts = 0;

   /* try to find a cut */
   SCIP_CALL( separatePointProduct(scip, conshdlr, expr, sol, &cut) );

   if( cut == NULL )
      return SCIP_OKAY;

   /* check whether its violation and numerical properties are ok (and maybe improve) */
   SCIP_CALL( SCIPmassageConsExprExprCut(scip, &cut, sol, minviolation) );

   if( cut == NULL )
      return SCIP_OKAY;

   assert(-SCIPgetRowSolFeasibility(scip, cut, sol) >= minviolation);

   /* add cut */
   SCIP_CALL( SCIPaddCut(scip, sol, cut, FALSE, &infeasible) );
   *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
   *ncuts += 1;

#ifdef SCIP_DEBUG
   SCIPdebugMessage("add cut with violation %e\n", violation);
   SCIP_CALL( SCIPprintRow(scip, cut, NULL) );
#endif

   /* release cut */
   SCIP_CALL( SCIPreleaseRow(scip, &cut) );

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropProduct)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_INTERVAL childbounds;
   SCIP_INTERVAL tmp;
   int i;
   int j;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) > 0);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   *nreductions = 0;
   *infeasible = FALSE;

   /* too expensive (runtime here is quadratic in number of children) */
   if( SCIPgetConsExprExprNChildren(expr) > 10 )
      return SCIP_OKAY;

   /* not possible to learn bounds if expression interval is unbounded in both directions */
   if( SCIPintervalIsEntire(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)) )
      return SCIP_OKAY;

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   /* f = const * prod_i c_i ^ n_i => c_i = (f / (const * prod_{j:j!=i} c_j ^ n_j) )^ (1/n_i) */
   for( i = 0; i < SCIPgetConsExprExprNChildren(expr) && !(*infeasible); ++i )
   {
      SCIPintervalSet(&childbounds, exprdata->constant);

      /* compute prod_{j:j!=i} c_j */
      for( j = 0; j < SCIPgetConsExprExprNChildren(expr); ++j )
      {
         if( i == j )
            continue;

         SCIPintervalPowerScalar(SCIPinfinity(scip), &tmp, SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[j]), exprdata->coefficients[j]);
         SCIPintervalMul(SCIPinfinity(scip), &childbounds, childbounds, tmp);

         /* if there is 0.0 in the product, then later division will hardly give useful bounds, so give up for this i */
         if( childbounds.inf <= 0.0 && childbounds.sup >= 0.0 )
            break;
      }

      if( j == SCIPgetConsExprExprNChildren(expr) )
      {
         /* f / (const * prod_{j:j!=i} c_j ^ n_j) */
         SCIPintervalDiv(SCIPinfinity(scip), &childbounds, SCIPgetConsExprExprInterval(expr), childbounds);

         /* ^(1/n_i) */
         SCIPintervalPowerScalarInverse(SCIPinfinity(scip), &childbounds,
            SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[i]), exprdata->coefficients[i], childbounds);

         /* try to tighten the bounds of the expression */
         SCIP_CALL( SCIPtightenConsExprExprInterval(scip, SCIPgetConsExprExprChildren(expr)[i], childbounds, infeasible, nreductions) );
      }
   }

   return SCIP_OKAY;
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
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataProduct, freedataProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifyProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashProduct) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffProduct) );

   return SCIP_OKAY;
}

/** creates a product expression */
SCIP_RETCODE SCIPcreateConsExprExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**  expr,               /**< pointer where to store expression */
   int                   nchildren,          /**< number of children */
   SCIP_CONSEXPR_EXPR**  children,           /**< children */
   SCIP_Real             constant            /**< constant coefficient of product */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   SCIP_CALL( createData(scip, &exprdata, nchildren, NULL, constant) );

   SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, SCIPgetConsExprExprHdlrProduct(consexprhdlr), exprdata, nchildren, children) );

   return SCIP_OKAY;
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

/** appends an expression to a product expression */
SCIP_RETCODE SCIPappendConsExprExprProductExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< product expression */
   SCIP_CONSEXPR_EXPR*   child               /**< expression to be appended */
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
   exprdata->coefficients[nchildren] = 1.0;

   SCIP_CALL( SCIPappendConsExprExpr(scip, expr, child) );

   return SCIP_OKAY;
}

