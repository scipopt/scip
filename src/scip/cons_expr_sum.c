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

#include "scip/cons_expr_sum.h"
#include "scip/cons_expr_value.h"

#define SUM_PRECEDENCE      40000
#define SUM_HASHKEY         SCIPcalcFibHash(47161.0)

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
 * so no values nor sum expressions, but products, variable or function expressions
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

/* helper function to compare two sums; see compareSum */
static
int compareSumHelper(
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
{  /*lint --e{715}*/
   return compareSumHelper(
         SCIPgetConsExprExprSumConstant(expr1), SCIPgetConsExprExprSumCoefs(expr1),
         SCIPgetConsExprExprChildren(expr1), SCIPgetConsExprExprNChildren(expr1),
         SCIPgetConsExprExprSumConstant(expr2), SCIPgetConsExprExprSumCoefs(expr2),
         SCIPgetConsExprExprChildren(expr2), SCIPgetConsExprExprNChildren(expr2));
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

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
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

/** expression derivative evaluation callback */
static
SCIP_DECL_CONSEXPR_EXPRBWDIFF(bwdiffSum)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(SCIPgetConsExprExprData(expr) != NULL);
   assert(idx >= 0 && idx < SCIPgetConsExprExprNChildren(expr));
   assert(SCIPgetConsExprExprChildren(expr)[idx] != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(expr)[idx])), "val") != 0);

   *val = SCIPgetConsExprExprSumCoefs(expr)[idx];

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

/** helper function to separate a given point; needed for proper unittest */
static
SCIP_RETCODE separatePointSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR*   expr,               /**< sum expression */
   SCIP_ROW**            cut                 /**< pointer to store the row */
   )
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_CONSEXPR_EXPR* child;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_VAR* auxvar;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), "expr") == 0);
   assert(expr != NULL);
   assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "sum") == 0);
   assert(cut != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   auxvar = SCIPgetConsExprExprLinearizationVar(expr);
   assert(auxvar != NULL);

   *cut = NULL;

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, SCIPgetConsExprExprNChildren(expr) + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coefs, SCIPgetConsExprExprNChildren(expr) + 1) );

   /* compute  w = sum_i alpha_i z_i + const */
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      child = SCIPgetConsExprExprChildren(expr)[c];

      /* value expressions should have been removed during simplification */
      assert(SCIPgetConsExprExprHdlr(child) != SCIPgetConsExprExprHdlrValue(conshdlr));

      vars[c] = SCIPgetConsExprExprLinearizationVar(child);
      assert(vars[c] != NULL);
      coefs[c] = exprdata->coefficients[c];
   }

   vars[SCIPgetConsExprExprNChildren(expr)] = auxvar;
   coefs[SCIPgetConsExprExprNChildren(expr)] = -1.0;

   /* create cut; it is globally valid because it is linear */
   SCIP_CALL( SCIPcreateRowCons(scip, cut, conshdlr, "sum_cut", 0, NULL, NULL, -exprdata->constant, -exprdata->constant,
         FALSE, FALSE, FALSE) );
   SCIP_CALL( SCIPaddVarsToRow(scip, *cut, SCIPgetConsExprExprNChildren(expr) + 1, vars, coefs) );

   /* release memory */
   SCIPfreeBufferArray(scip, &coefs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** separation initialization callback */
static
SCIP_DECL_CONSEXPR_EXPRINITSEPA(initSepaSum)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);
   assert(exprdata->row == NULL);

   *infeasible = FALSE;

   /* compute initial cut */
   SCIP_CALL( separatePointSum(scip, conshdlr, expr, &exprdata->row) );
   assert(exprdata->row != NULL);

   SCIP_CALL( SCIPaddCut(scip, NULL, exprdata->row, FALSE, infeasible) );

   return SCIP_OKAY;
}

/** separation deinitialization callback */
static
SCIP_DECL_CONSEXPR_EXPREXITSEPA(exitSepaSum)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   if( exprdata->row != NULL )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &exprdata->row) );
   }
   assert(exprdata->row == NULL);

   return SCIP_OKAY;
}

/** expression separation callback */
static
SCIP_DECL_CONSEXPR_EXPRSEPA(sepaSum)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_Real violation;

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);
   assert(exprdata->row != NULL);

   *result = SCIP_DIDNOTFIND;

   /* row should not be violated if it is still in the LP */
   if( exprdata->row != NULL && SCIProwIsInLP(exprdata->row) )
   {
      /* We expect that the LP solution is feasible w.r.t. absolute tolerances,
       * that is, that the row feasibility is >= -feastol.
       * This is also checked by lp.c:SCIPlpGetSol (see also #1041).
       */
      assert(!SCIPisFeasNegative(scip, SCIPgetRowSolFeasibility(scip, exprdata->row, sol)));

      return SCIP_OKAY;
   }

   /* compute and store cut (might happen if the constraint is not 'initial') */
   if( exprdata->row == NULL )
   {
      SCIP_CALL( separatePointSum(scip, conshdlr, expr, &exprdata->row) );
   }
   assert(exprdata->row != NULL);

   /* compute violation */
   violation = -SCIPgetRowSolFeasibility(scip, exprdata->row, sol);

   if( SCIPisGE(scip, violation, minviolation) )
   {
      SCIP_Bool infeasible;

      SCIP_CALL( SCIPaddCut(scip, sol, exprdata->row, FALSE, &infeasible) );
      *result = infeasible ? SCIP_CUTOFF : SCIP_SEPARATED;
      *ncuts += 1;

#ifdef SCIP_DEBUG
      SCIPdebugMessage("add cut with violation %e\n", violation);
      SCIP_CALL( SCIPprintRow(scip, exprdata->row, NULL) );
#endif
   }

   return SCIP_OKAY;
}

/** expression reverse propagation callback */
static
SCIP_DECL_CONSEXPR_REVERSEPROP(reversepropSum)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   SCIP_INTERVAL* bounds;
   SCIP_ROUNDMODE prevroundmode;
   SCIP_INTERVAL childbounds;
   SCIP_Real minlinactivity;
   SCIP_Real maxlinactivity;
   int minlinactivityinf;
   int maxlinactivityinf;
   int c;
   assert(scip != NULL);
   assert(expr != NULL);
   assert(SCIPgetConsExprExprNChildren(expr) > 0);
   assert(infeasible != NULL);
   assert(nreductions != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *infeasible = FALSE;
   *nreductions = 0;

   /* not possible to conclude finite bounds if the interval of the expression is [-inf,inf] */
   if( SCIPintervalIsEntire(SCIPinfinity(scip), SCIPgetConsExprExprInterval(expr)) )
      return SCIP_OKAY;

   /* TODO handle cases nchildren = 1, 2 separately to be more efficient */

   prevroundmode = SCIPintervalGetRoundingMode();
   SCIPintervalSetRoundingModeDownwards();

   minlinactivity = exprdata->constant;
   maxlinactivity = -exprdata->constant; /* use -constant because of the rounding mode */
   minlinactivityinf = 0;
   maxlinactivityinf = 0;

   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, SCIPgetConsExprExprNChildren(expr)) );

   /* shift coefficients into the intervals of the children; compute the min and max activities */
   for( c = 0; c < SCIPgetConsExprExprNChildren(expr); ++c )
   {
      SCIPintervalMulScalar(SCIPinfinity(scip), &bounds[c], SCIPgetConsExprExprInterval(SCIPgetConsExprExprChildren(expr)[c]),
         exprdata->coefficients[c]);

      if( SCIPisInfinity(scip, SCIPintervalGetSup(bounds[c])) )
         ++maxlinactivityinf;
      else
      {
         assert(SCIPintervalGetSup(bounds[c]) > -SCIPinfinity(scip));
         maxlinactivity -= SCIPintervalGetSup(bounds[c]);
      }

      if( SCIPisInfinity(scip, -SCIPintervalGetInf(bounds[c])) )
         ++minlinactivityinf;
      else
      {
         assert(SCIPintervalGetInf(bounds[c]) < SCIPinfinity(scip));
         minlinactivity += SCIPintervalGetInf(bounds[c]);
      }
   }
   maxlinactivity = -maxlinactivity; /* correct sign */

   /* if there are too many unbounded bounds, then could only compute infinite bounds for children, so give up */
   if( (minlinactivityinf >= 2 || SCIPisInfinity(scip, SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)))) &&
      ( maxlinactivityinf >= 2 || SCIPisInfinity(scip, -SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr))))
      )
   {
      goto TERMINATE;
   }

   for( c = 0; c < SCIPgetConsExprExprNChildren(expr) && !(*infeasible); ++c )
   {
      /* upper bounds of c_i is
       *   node->bounds.sup - (minlinactivity - c_i.inf), if c_i.inf > -infinity and minlinactivityinf == 0
       *   node->bounds.sup - minlinactivity, if c_i.inf == -infinity and minlinactivityinf == 1
       */
      SCIPintervalSetEntire(SCIPinfinity(scip), &childbounds);
      if( !SCIPisInfinity(scip, SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr))) )
      {
         /* we are still in downward rounding mode, so negate and negate to get upward rounding */
         if( bounds[c].inf <= -SCIPinfinity(scip) && minlinactivityinf <= 1 )
         {
            assert(minlinactivityinf == 1);
            childbounds.sup = SCIPintervalNegateReal(minlinactivity - SCIPgetConsExprExprInterval(expr).sup);
         }
         else if( minlinactivityinf == 0 )
         {
            childbounds.sup = SCIPintervalNegateReal(minlinactivity - SCIPgetConsExprExprInterval(expr).sup - bounds[c].inf);
         }
      }

      /* lower bounds of c_i is
       *   node->bounds.inf - (maxlinactivity - c_i.sup), if c_i.sup < infinity and maxlinactivityinf == 0
       *   node->bounds.inf - maxlinactivity, if c_i.sup == infinity and maxlinactivityinf == 1
       */
      if( SCIPgetConsExprExprInterval(expr).inf > -SCIPinfinity(scip) )
      {
         if( bounds[c].sup >= SCIPinfinity(scip) && maxlinactivityinf <= 1 )
         {
            assert(maxlinactivityinf == 1);
            childbounds.inf = SCIPgetConsExprExprInterval(expr).inf - maxlinactivity;
         }
         else if( maxlinactivityinf == 0 )
         {
            childbounds.inf = SCIPgetConsExprExprInterval(expr).inf - maxlinactivity + bounds[c].sup;
         }
      }

      /* divide by the child coefficient */
      SCIPintervalDivScalar(SCIPinfinity(scip), &childbounds, childbounds, exprdata->coefficients[c]);

      /* try to tighten the bounds of the expression */
      SCIP_CALL( SCIPtightenConsExprExprInterval(scip, SCIPgetConsExprExprChildren(expr)[c], childbounds, infeasible, nreductions) );
   }

TERMINATE:
   SCIPintervalSetRoundingMode(prevroundmode);
   SCIPfreeBufferArray(scip, &bounds);

   return SCIP_OKAY;
}

/** sum hash callback */
static
SCIP_DECL_CONSEXPR_EXPRHASH(hashSum)
{
   SCIP_CONSEXPR_EXPRDATA* exprdata;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(expr2key != NULL);
   assert(hashkey != NULL);

   exprdata = SCIPgetConsExprExprData(expr);
   assert(exprdata != NULL);

   *hashkey = SUM_HASHKEY;
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
   SCIP_CALL( SCIPsetConsExprExprHdlrCopyFreeData(scip, consexprhdlr, exprhdlr, copydataSum, freedataSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSimplify(scip, consexprhdlr, exprhdlr, simplifySum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrCompare(scip, consexprhdlr, exprhdlr, compareSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrPrint(scip, consexprhdlr, exprhdlr, printSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrIntEval(scip, consexprhdlr, exprhdlr, intevalSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrInitSepa(scip, consexprhdlr, exprhdlr, initSepaSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrExitSepa(scip, consexprhdlr, exprhdlr, exitSepaSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrSepa(scip, consexprhdlr, exprhdlr, sepaSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrReverseProp(scip, consexprhdlr, exprhdlr, reversepropSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrHash(scip, consexprhdlr, exprhdlr, hashSum) );
   SCIP_CALL( SCIPsetConsExprExprHdlrBwdiff(scip, consexprhdlr, exprhdlr, bwdiffSum) );

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

