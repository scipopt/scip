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

/*
 * local functions
 */

/** variable mapping data passed on during copying expressions when copying SCIP instances */
typedef struct
{
   SCIP_HASHMAP*         varmap;             /**< SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP */
   SCIP_HASHMAP*         consmap;            /**< SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP */
   SCIP_Bool             global;             /**< should a global or a local copy be created */
   SCIP_Bool             valid;              /**< indicates whether every variable copy was valid */
} COPY_MAPEXPR_DATA;

/** variable expression mapping callback to call when copying expressions (within same or different SCIPs) */
static
SCIP_DECL_EXPR_MAPEXPR(copyVarExpr)
{
   COPY_MAPEXPR_DATA* data;
   SCIP_Bool valid;
   SCIP_VAR* targetvar;

   assert(sourcescip != NULL);
   assert(sourceexpr != NULL);
   assert(targetscip != NULL);
   assert(targetexpr != NULL);
   assert(mapexprdata != NULL);

   *targetexpr = NULL;

   if( !SCIPisExprVar(sourcescip, sourceexpr) )
      return SCIP_OKAY;

   data = (COPY_MAPEXPR_DATA*)mapexprdata;

   SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, SCIPgetConsExprExprVarVar(sourceexpr), &targetvar, data->varmap, data->consmap, data->global, &valid) );
   assert(*targetvar != NULL);

   /* if copy was not valid, store so in mapvar data */
   if( !valid )
      data->valid = FALSE;

   SCIP_CALL( SCIPcreateConsExprExprVar(targetscip, NULL, targetexpr, targetvar) );

   return SCIP_OKAY;
}

/** @name Simplifying expressions (hashing, common subexpressions, simplify)
 * @{

/** returns an equivalent expression for a given expression if possible
 *
 * it adds the expression to key2expr if the map does not contain the key
 */
static
SCIP_RETCODE findEqualExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_EXPR*            expr,               /**< expression to replace */
   SCIP_MULTIHASH*       key2expr,           /**< mapping of hashes to expressions */
   SCIP_EXPR**           newexpr             /**< pointer to store an equivalent expression (NULL if there is none) */
   )
{  /*lint --e{438}*/
   SCIP_MULTIHASHLIST* multihashlist;

   assert(set != NULL);
   assert(expr != NULL);
   assert(key2expr != NULL);
   assert(newexpr != NULL);

   *newexpr = NULL;
   multihashlist = NULL;
   do
   {
      /* search for an equivalent expression */
      *newexpr = (SCIP_EXPR*)(SCIPmultihashRetrieveNext(key2expr, &multihashlist, (void*)expr));

      if( *newexpr == NULL )
      {
         /* processed all expressions like expr from hash table, so insert expr */
         SCIP_CALL( SCIPmultihashInsert(key2expr, (void*) expr) );
         break;
      }
      else if( expr != *newexpr )
      {
         assert(SCIPexprCompare(set, expr, *newexpr) == 0);
         break;
      }
      else
      {
         /* can not replace expr since it is already contained in the hashtablelist */
         assert(expr == *newexpr);
         *newexpr = NULL;
         break;
      }
   }
   while( TRUE ); /*lint !e506*/

   return SCIP_OKAY;
}

/** get key of hash element */
static
SCIP_DECL_HASHGETKEY(hashCommonSubexprGetKey)
{
   return elem;
}  /*lint !e715*/

/** checks if two expressions are structurally the same */
static
SCIP_DECL_HASHKEYEQ(hashCommonSubexprEq)
{
   SCIP_EXPR* expr1;
   SCIP_EXPR* expr2;

   expr1 = (SCIP_EXPR*)key1;
   expr2 = (SCIP_EXPR*)key2;
   assert(expr1 != NULL);
   assert(expr2 != NULL);

   return expr1 == expr2 || SCIPexprCompare(expr1, expr2) == 0;
}  /*lint !e715*/

/** get value of hash element when comparing with another expression */
static
SCIP_DECL_HASHKEYVAL(hashCommonSubexprKeyval)
{
   SCIP_EXPR* expr;
   SCIP_EXPRITER* hashiterator;

   expr = (SCIP_EXPR*) key;
   assert(expr != NULL);

   hashiterator = (SCIP_EXPRITER*) userptr;
   assert(hashiterator != NULL);

   return SCIPexpriterGetExprUserData(hashiterator, expr).uintval;
}  /*lint !e715*/

/** hashes an expression using an already existing iterator
 *
 * The iterator must by of type DFS with allowrevisit=FALSE and only the leaveexpr stage enabled.
 * The hashes of all visited expressions will be stored in the iterators expression data.
 */
static
SCIP_RETCODE hashExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_EXPR*            expr,               /**< expression to hash */
   SCIP_EXPRITER*        hashiterator,       /**< iterator to use for hashing */
   int*                  nvisitedexprs       /**< counter to increment by the number of expressions visited, or NULL */
   )
{
   SCIP_EXPRITER_USERDATA iterdata;
   unsigned int* childrenhashes;
   int childrenhashessize;
   int i;

   assert(set != NULL);
   assert(expr != NULL);
   assert(hashiterator != NULL);

   childrenhashessize = 5;
   SCIP_ALLOC( BMSallocBufferMemoryArray(bufmem, &childrenhashes, childrenhashessize) );

   for( expr = SCIPexpriterRestartDFS(hashiterator, expr); !SCIPexpriterIsEnd(hashiterator); expr = SCIPexpriterGetNext(hashiterator) ) /*lint !e441*/
   {
      assert(SCIPexpriterGetStageDFS(hashiterator) == SCIP_EXPRITER_LEAVEEXPR);

      if( nvisitedexprs != NULL )
         ++*nvisitedexprs;

      /* collect hashes of children */
      if( childrenhashessize < expr->nchildren )
      {
         childrenhashessize = SCIPsetCalcMemGrowSize(set, expr->nchildren);
         SCIP_ALLOC( BMSreallocBufferMemoryArray(bufmem, &childrenhashes, childrenhashessize) );
      }
      for( i = 0; i < expr->nchildren; ++i )
         childrenhashes[i] = SCIPexpriterGetExprUserData(hashiterator, expr->children[i]).uintval;

      SCIP_CALL( SCIPcallExprhdlrHash(set->scip, expr, &iterdata.uintval, childrenhashes) );

      SCIPexpriterSetCurrentUserData(hashiterator, iterdata);
   }

   BMSfreeBufferMemoryArray(bufmem, &childrenhashes);

   return SCIP_OKAY;
}

/** replaces common sub-expressions in a given expression graph by using a hash key for each expression
 *
 *  The algorithm consists of two steps:
 *
 *  1. traverse through all given expressions and compute for each of them a (not necessarily unique) hash
 *
 *  2. initialize an empty hash table and traverse through all expression; check for each of them if we can find a
 *     structural equivalent expression in the hash table; if yes we replace the expression by the expression inside the
 *     hash table, otherwise we add it to the hash table
 *
 *  @note the hash keys of the expressions are used for the hashing inside the hash table; to compute if two expressions
 *  (with the same hash) are structurally the same we use the function SCIPexprCompare()
 */
static
SCIP_RETCODE replaceCommonSubexpressions(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_EXPR**           exprs,              /**< expressions (possibly replaced by equivalent on output) */
   int                   nexprs,             /**< total number of expressions */
   SCIP_Bool*            replacedroot        /**< buffer to store whether any root expression (expression in exprs) was replaced */
   )
{
   SCIP_EXPRITER* hashiterator;
   SCIP_EXPRITER* repliterator;
   SCIP_MULTIHASH* key2expr;
   SCIP_CONSDATA* consdata;
   int i;
   int nvisitedexprs = 0;

   assert(set != NULL);
   assert(stat != NULL);
   assert(exprs != NULL);
   assert(nexprs >= 0);
   assert(replacedroot != NULL);

   if( nexprs == 0 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &hashiterator) );
   SCIP_CALL( SCIPexpriterInit(hashiterator, NULL, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(hashiterator, SCIP_EXPRITER_LEAVEEXPR);

   /* compute all hashes for each sub-expression */
   for( i = 0; i < nexprs; ++i )
   {
      assert(exprs[i] != NULL);
      SCIP_CALL( hashExpr(set, bufmem, exprs[i], hashiterator, &nvisitedexprs) );
   }

   /* replace equivalent sub-expressions */
   SCIP_CALL( SCIPmultihashCreate(&key2expr, blkmem, nvisitedexprs,
         hashCommonSubexprGetKey, hashCommonSubexprEq, hashCommonSubexprKeyval, (void*)hashiterator) );

   SCIP_CALL( SCIPexpriterCreate(stat, blkmem, &repliterator) );

   for( i = 0; i < nexprs; ++i )
   {
      SCIP_EXPR* newroot;
      SCIP_EXPR* newchild;
      SCIP_EXPR* child;

      /* check the root for equivalence separately first */
      SCIP_CALL( findEqualExpr(set, exprs[i], key2expr, &newroot) );

      if( newroot != NULL )
      {
         assert(newroot != exprs[i]);
         assert(SCIPexprCompare(exprs[i], newroot) == 0);

         SCIPdebugMsg(scip, "replacing common root expression of %dth expr: %p -> %p\n", i, (void*)exprs[i], (void*)newroot);

         SCIP_CALL( SCIPreleaseExpr(set->scip, &exprs[i]) );

         exprs[i] = newroot;
         SCIPexprCapture(newroot);

         *replacedroot = TRUE;

         continue;
      }

      /* replace equivalent sub-expressions in the tree */
      SCIP_CALL( SCIPexpriterInit(repliterator, exprs[i], SCIP_EXPRITER_DFS, FALSE) );
      SCIPexpriterSetStagesDFS(repliterator, SCIP_EXPRITER_VISITINGCHILD);

      while( !SCIPexpriterIsEnd(repliterator) )
      {
         child = SCIPexpriterGetChildExprDFS(repliterator);
         assert(child != NULL);

         /* try to find an equivalent expression */
         SCIP_CALL( findEqualExpr(set, child, key2expr, &newchild) );

         /* replace child with newchild */
         if( newchild != NULL )
         {
            assert(child != newchild);
            assert(SCIPexprCompare(child, newchild) == 0);

            SCIPdebugMsg(scip, "replacing common child expression %p -> %p\n", (void*)child, (void*)newchild);

            SCIP_CALL( SCIPreplaceExprChild(scip, SCIPexpriterGetCurrent(repliterator), SCIPexpriterGetChildIdxDFS(repliterator), newchild) );

            (void) SCIPexpriterSkipDFS(repliterator);
         }
         else
         {
            (void) SCIPexpriterGetNext(repliterator);
         }
      }
   }

   /* free memory */
   SCIPexpriteratorFree(&repliterator);
   SCIPmultihashFree(&key2expr);
   SCIPexpriteratorFree(&hashiterator);

   return SCIP_OKAY;
}

/** @} */  /* end of simplify methods */

/*
 * public functions
 */

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
SCIP_RETCODE SCIPduplicateExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_DECL_EXPR_MAPEXPR((*mapexpr)),       /**< expression mapping function, or NULL for creating new expressions */
   void*                 mapexprdata,        /**< data of expression mapping function */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call on expression copy to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata, /**< data to pass to ownerdatacreate */
   SCIP_DECL_EXPR_OWNERDATAFREE((*ownerdatafree)),     /**< function to call when freeing expression, e.g., to free ownerdata */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprCopy(scip->set, scip->stat, scip->mem->probmem, scip->set, scip->stat, scip->mem->probmem, expr, copyexpr, mapvar, mapvardata, mapexpr, mapexprdata, ownerdatacreate, ownerdatacreatedata, ownerdatafree) );

   return SCIP_OKAY;
}

/** duplicates the given expression without its children */
SCIP_RETCODE SCIPduplicateExprShallow(
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

/** copies an expression to use in a (possibly different) SCIP instance (including children) */
SCIP_RETCODE SCIPcopyExpr(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 targetscip,         /**< target SCIP data structure */
   SCIP_EXPR*            expr,               /**< original expression */
   SCIP_EXPR**           copyexpr,           /**< buffer to store duplicate of expr */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call on expression copy to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata, /**< data to pass to ownerdatacreate */
   SCIP_DECL_EXPR_OWNERDATAFREE((*ownerdatafree)),     /**< function to call when freeing expression, e.g., to free ownerdata */
   SCIP_HASHMAP*         varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to the corresponding
                                              *   variables of the target SCIP, or NULL */
   SCIP_HASHMAP*         consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
   SCIP_Bool             global,             /**< create a global or a local copy? */
   SCIP_Bool*            valid,              /**< pointer to store whether all checked or enforced constraints were validly copied */
   )
{
   assert(sourcescip != NULL);
   assert(sourcescip->mem != NULL);
   assert(targetscip != NULL);
   assert(targetscip->mem != NULL);

   COPY_MAPEXPR_DATA copydata = {
      .varmap = varmap,
      .consmap = consmap,
      .global = global,
      .valid = TRUE
   };

   SCIP_CALL( SCIPexprCopy(sourcescip->set, sourcescip->stat, sourcescip->mem->probmem,
      targetscip->set, targetscip->stat, targetscip->mem->probmem,
      expr, copyexpr, copyVarExpr, &copydata, ownerdatacreate, ownerdatacreatedata, ownerdatafree) );

   *valid = copydata.valid;

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

   SCIP_CALL( SCIPexprPrint(scip->set, scip->stat, scip->mem->probmem, scip->messagehdlr, expr, file) );

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

   SCIP_CALL( SCIPexprDismantle(scip->set, scip->stat, scip->mem->probmem, scip->messagehdlr, file, expr) );

   return SCIP_OKAY;
}

/** evaluate an expression in a point
 *
 * Iterates over expressions to also evaluate children, if necessary.
 * Value can be received via SCIPexprGetEvalValue().
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 *
 * If a nonzero \p soltag is passed, then only (sub)expressions are
 * reevaluated that have a different solution tag. If a soltag of 0
 * is passed, then subexpressions are always reevaluated.
 * The tag is stored together with the value and can be received via
 * SCIPexprGetEvalTag().
 */
SCIP_RETCODE SCIPevalExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprEval(scip->set, scip->stat, scip->mem->probmem, expr, sol, soltag) );

   return SCIP_OKAY;
}

/** evaluates gradient of an expression for a given point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_RETCODE SCIPevalExprGradient(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   SCIP_SOL*             sol,                /**< solution to be evaluated (NULL for the current LP solution) */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprEvalGradient(scip->set, scip->stat, scip->mem->probmem, expr, sol, soltag) );

   return SCIP_OKAY;
}

/** evaluates Hessian-vector product of an expression for a given point and direction
 *
 * Evaluates children, if necessary.
 * Value can be received via SCIPgetExprPartialDiffGradientDirNonlinear().
 * If an error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 */
SCIP_RETCODE SCIPevalExprHessianDir(
   SCIP*                 scip,             /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression to be differentiated */
   SCIP_SOL*             sol,              /**< solution to be evaluated (NULL for the current LP solution) */
   unsigned int          soltag,           /**< tag that uniquely identifies the solution (with its values), or 0. */
   SCIP_SOL*             direction         /**< direction */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprEvalHessianDir(scip->set, scip->stat, scip->mem->probmem, expr, sol, soltag, direction) );

   return SCIP_OKAY;
}

/** evaluates activity of expression w.r.t. current local variable bounds
 *
 * Value van be received via SCIPexprGetActivity().
 *
 * Reevaluate activity if any variable has changed bounds since last eval call.
 */
SCIP_RETCODE SCIPevalExprActivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprEvalActivity(scip->set, scip->stat, scip->mem->probmem, expr) );

   return SCIP_OKAY;
}

/** compute the hash value of an expression */
SCIP_RETCODE SCIPhashExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   unsigned int*         hashval             /**< pointer to store the hash value */
   )
{
   SCIP_EXPRITER* it;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(expr != NULL);
   assert(hashval != NULL);

   SCIP_CALL( SCIPexpriterCreate(scip->stat, scip->mem->probmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   SCIP_CALL( hashExpr(scip->set, scip->mem->buffer, expr, it, NULL) );

   *hashval = SCIPexpriterGetExprUserData(it, expr).uintval;

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/* simplifies an expression (duplication of long doxygen comment omitted here) */
SCIP_RETCODE SCIPsimplifyExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            rootexpr,           /**< expression to be simplified */
   SCIP_EXPR**           simplified,         /**< buffer to store simplified expression */
   SCIP_Bool*            changed,            /**< buffer to store if rootexpr actually changed */
   SCIP_Bool*            infeasible          /**< buffer to store whether infeasibility has been detected */
   )
{
   SCIP_EXPR* expr;
   SCIP_EXPRITER* it;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(rootexpr != NULL);
   assert(simplified != NULL);
   assert(changed != NULL);
   assert(infeasible != NULL);

   /* simplify bottom up
    * when leaving an expression it simplifies it and stores the simplified expr in its iterators expression data
    * after the child was visited, it is replaced with the simplified expr
    */
   SCIP_CALL( SCIPexpriterCreate(scip->stat, scip->mem->probmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, rootexpr, SCIP_EXPRITER_DFS, TRUE) );  /* TODO can we set allowrevisited to FALSE?*/
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_VISITEDCHILD | SCIP_EXPRITER_LEAVEEXPR);

   *changed = FALSE;
   *infeasible = FALSE;
   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) ) /*lint !e441*/
   {
      switch( SCIPexpriterGetStageDFS(it) )
      {
         case SCIP_EXPRITER_VISITEDCHILD:
         {
            SCIP_EXPR* newchild;
            SCIP_EXPR* child;

            newchild = (SCIP_EXPR*)SCIPexpriterGetChildUserDataDFS(it).ptrval;
            child = SCIPexpriterGetChildExprDFS(it);
            assert(newchild != NULL);

            /* if child got simplified, replace it with the new child */
            if( newchild != child )
            {
               SCIP_CALL( SCIPreplaceExprChild(scip, expr, SCIPexpriterGetChildIdxDFS(it), newchild) );
            }

            /* we do not need to hold newchild anymore */
            SCIP_CALL( SCIPreleaseExpr(scip, &newchild) );

            break;
         }

         case SCIP_EXPRITER_LEAVEEXPR:
         {
            SCIP_EXPR* refexpr = NULL;
            SCIP_EXPRITER_USERDATA iterdata;

            /* TODO we should do constant folding (handle that all children are value-expressions) here in a generic way
             * instead of reimplementing it in every handler
             */

            /* use simplification of expression handlers */
            SCIP_CALL( SCIPexprhdlrSimplifyExpr(expr->exprhdlr, scip->set, expr, &refexpr) );
            assert(refexpr != NULL);
            if( expr != refexpr )
               *changed = TRUE;

            iterdata.ptrval = (void*) refexpr;
            SCIPexpriterSetCurrentUserData(it, iterdata);

            break;
         }

         default:
            SCIPABORT(); /* we should never be called in this stage */
            break;
      }
   }

   *simplified = (SCIP_EXPR*)SCIPexpriterGetExprUserData(it, rootexpr).ptrval;
   assert(*simplified != NULL);

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** computes the curvature of a given expression and all its subexpressions
 *
 *  @note this function also evaluates all subexpressions w.r.t. current variable bounds
 */
SCIP_RETCODE SCIPcomputeExprCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPRITER* it;
   SCIP_EXPRCURV curv;
   SCIP_EXPRCURV* childcurv;
   int childcurvsize;
   SCIP_Bool success;
   SCIP_EXPRCURV trialcurv[3] = { SCIP_EXPRCURV_LINEAR, SCIP_EXPRCURV_CONVEX, SCIP_EXPRCURV_CONCAVE };
   int i, c;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(expr != NULL);

   childcurvsize = 5;
   SCIP_CALL( SCIPallocBufferArray(scip, &childcurv, childcurvsize) );

   SCIP_CALL( SCIPexpriterCreate(scip->stat, scip->mem->probmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      curv = SCIP_EXPRCURV_UNKNOWN;

      if( expr->exprhdlr->curvature == NULL )
      {
         /* set curvature in expression */
         SCIPexprSetCurvature(expr, curv);
         continue;
      }

      if( SCIPexprGetNChildren(expr) > childcurvsize )
      {
         childcurvsize = SCIPcalcMemGrowSize(scip, SCIPexprGetNChildren(expr));
         SCIP_CALL( SCIPreallocBufferArray(scip, &childcurv, childcurvsize) );
      }

      for( i = 0; i < 3; ++i )
      {
         /* check if expression can have a curvature trialcurv[i] */
         SCIP_CALL( SCIPexprhdlrCurvatureExpr(expr->exprhdlr, set, expr, trialcurv[i], &success, childcurv) );
         if( !success )
            continue;

         /* check if conditions on children are satisfied */
         for( c = 0; c < SCIPexprGetNChildren(expr); ++c )
         {
            if( (childcurv[c] & SCIPexprGetCurvature(SCIPexprGetChildren(expr)[c])) != childcurv[c] )
            {
               success = FALSE;
               break;
            }
         }

         if( success )
         {
            curv = trialcurv[i];
            break;
         }
      }

      /* set curvature in expression */
      SCIPexprSetCurvature(expr, curv);
   }

   SCIPexpriteratorFree(&it);

   SCIPfreeBufferArray(scip, &childcurv);

   return SCIP_OKAY;
}

/** get the monotonicity of an expression w.r.t. to a given child */
SCIP_RETCODE SCIPgetExprMonotonicity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   int                   childidx,           /**< index of child */
   SCIP_MONOTONE*        monotonicity        /**< buffer to store monotonicity */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);

   SCIP_CALL( SCIPexprhdlrMonotonicityExpr(expr->exprhdlr, scip->set, expr, childidx, monotonicity) );

   return SCIP_OKAY;
}

/** computes integrality information of a given expression and all its subexpressions
 *
 * the integrality information can be accessed via SCIPexprIsIntegral()
 */
SCIP_RETCODE SCIPcomputeExprIntegrality(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   SCIP_EXPRITER* it;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(expr != NULL);

   /* shortcut for expr without children */
   if( expr->nchildren == 0 )
   {
      /* compute integrality information */
      SCIP_CALL( SCIPexprhdlrIntegralityExpr(expr->exprhdlr, scip->set, expr, &expr->isintegral) );

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPexpriterCreate(scip->stat, scip->mem->probmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );
   SCIPexpriterSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   for( expr = SCIPexpriterGetCurrent(it); !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      /* compute integrality information */
      SCIP_CALL( SCIPexprhdlrIntegralityExpr(expr->exprhdlr, scip->set, expr, &expr->isintegral) );
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** returns the total number of variables in an expression
 *
 * The function counts variables in common sub-expressions only once.
 */
SCIP_RETCODE SCIPgetExprNVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   int*                  nvars               /**< buffer to store the total number of variables */
   )
{
   SCIP_EXPRITER* it;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(expr != NULL);
   assert(nvars != NULL);

   SCIP_CALL( SCIPexpriterCreate(scip->stat, scip->mem->probmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );

   *nvars = 0;
   for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
      if( SCIPexprIsVar(scip->set, expr) )
         ++(*nvars);

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}

/** returns all variable expressions contained in a given expression
 *
 * the array to store all variable expressions needs
 * to be at least of size the number of unique variables in the expression which is given by SCIPgetExprNVars()
 * and can be bounded by SCIPgetNVars().
 *
 * @note function captures variable expressions
 */
SCIP_RETCODE SCIPgetExprVarExprs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_EXPR**           varexprs,           /**< array to store all variable expressions */
   int*                  nvarexprs           /**< buffer to store the total number of variable expressions */
   )
{
   SCIP_EXPRITER* it;

   assert(scip != NULL);
   assert(scip->mem != NULL);
   assert(expr != NULL);
   assert(varexprs != NULL);
   assert(nvarexprs != NULL);

   SCIP_CALL( SCIPexpriterCreate(scip->stat, scip->mem->probmem, &it) );
   SCIP_CALL( SCIPexpriterInit(it, expr, SCIP_EXPRITER_DFS, FALSE) );

   *nvarexprs = 0;
   for( ; !SCIPexpriterIsEnd(it); expr = SCIPexpriterGetNext(it) )
   {
      assert(expr != NULL);

      if( SCIPexprIsVar(scip->set, expr) )
      {
         /* add variable expression to array and capture expr */
         assert(SCIPgetNTotalVars(scip) >= *nvarexprs + 1);

         varexprs[(*nvarexprs)++] = expr;

         /* capture expression */
         SCIPcaptureExpr(expr);
      }
   }

   /* @todo sort variable expressions here? */

   SCIPexpriteratorFree(&it);

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

   SCIP_CALL( SCIPexpriterCreate(scip->stat, scip->mem->probmem, iterator) );

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


/**@name Quadratic expression functions */
/**@{ */

/** checks whether an expression is quadratic
 *
 * An expression is quadratic if it is either a square (of some expression), a product (of two expressions),
 * or a sum of terms where at least one is a square or a product.
 *
 * Use \ref SCIPexprGetQuadraticData to get data about the representation as quadratic.
 */
SCIP_RETCODE SCIPcheckExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< expression */
   SCIP_Bool*            isquadratic         /**< buffer to store result */
   )
{
   assert(scip != NULL);
   assert(scip->mem != NULL);

   SCIP_CALL( SCIPexprCheckQuadratic(scip->set, scip->mem->probmem, expr, isquadratic) );

   return SCIP_OKAY;
}

/** evaluates quadratic term in a solution
 *
 * \note This requires that every expr used in the quadratic data is a variable expression.
 */
SCIP_Real SCIPevalExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_SOL*             sol                 /**< solution to evaluate, or NULL for LP solution */
   )
{
   SCIP_QUADEXPR* quaddata;
   SCIP_Real auxvalue;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);

   quaddata = expr->quaddata;
   assert(expr->quaddata != NULL);

   auxvalue = quaddata->constant;

   for( i = 0; i < quaddata->nlinexprs; ++i ) /* linear exprs */
   {
      assert(SCIPexprIsVar(scip->set, quaddata->linexprs[i]));
      auxvalue += quaddata->lincoefs[i] * SCIPgetSolVal(scip, sol, SCIPgetConsExprExprVarVar(quaddata->linexprs[i]));
   }

   for( i = 0; i < quaddata->nquadexprs; ++i ) /* quadratic terms */
   {
      SCIP_QUADEXPR_QUADTERM* quadexprterm;
      SCIP_Real solval;

      quadexprterm = &quaddata->quadexprterms[i];
      assert(SCIPexprIsVar(scip->set, quadexprterm->expr));

      solval = SCIPgetSolVal(scip, sol, SCIPgetConsExprExprVarVar(quadexprterm->expr));
      auxvalue += (quadexprterm->lincoef + quadexprterm->sqrcoef * solval) * solval;
   }

   for( i = 0; i < quaddata->nbilinexprterms; ++i ) /* bilinear terms */
   {
      SCIP_QUADEXPR_BILINTERM* bilinexprterm;

      bilinexprterm = &quaddata->bilinexprterms[i];
      assert(SCIPexprIsVar(scip->set, bilinexprterm->expr1));
      assert(SCIPexprIsVar(scip->set, bilinexprterm->expr2));
      auxvalue += bilinexprterm->coef *
         SCIPgetSolVal(scip, sol, SCIPgetConsExprExprVarVar(bilinexprterm->expr1)) *
         SCIPgetSolVal(scip, sol, SCIPgetConsExprExprVarVar(bilinexprterm->expr2));
   }

   return auxvalue;
}

/** prints quadratic expression */
SCIP_RETCODE SCIPprintExprQuadratic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr                /**< quadratic expression */
   )
{
   SCIP_QUADEXPR* quaddata;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);

   quaddata = expr->quaddata;
   assert(quaddata != NULL);

   SCIPinfoMessage(scip, NULL, "Constant: %g\n", quaddata->constant);
   SCIPinfoMessage(scip, NULL, "Linear: ");
   for( c = 0; c < quaddata->nlinexprs; ++c )
   {
      SCIPinfoMessage(scip, NULL, "%g * ", quaddata->lincoefs[c]);
      SCIP_CALL( SCIPprintExpr(scip, quaddata->linexprs[c], NULL) );
      if( c < quaddata->nlinexprs-1 )
         SCIPinfoMessage(scip, NULL, " + ");
   }
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "Quadratic: ");
   for( c = 0; c < quaddata->nquadexprs; ++c )
   {
      SCIPinfoMessage(scip, NULL, "(%g * sqr(", quaddata->quadexprterms[c].sqrcoef);
      SCIP_CALL( SCIPprintExpr(scip, quaddata->quadexprterms[c].expr, NULL) );
      SCIPinfoMessage(scip, NULL, ") + %g) * ", quaddata->quadexprterms[c].lincoef);
      SCIP_CALL( SCIPprintExpr(scip, quaddata->quadexprterms[c].expr, NULL) );
      if( c < quaddata->nquadexprs-1 )
         SCIPinfoMessage(scip, NULL, " + ");
   }
   SCIPinfoMessage(scip, NULL, "\n");
   if( quaddata->nbilinexprterms > 0 )
   {
      SCIPinfoMessage(scip, NULL, "Bilinear: ");
      for( c = 0; c < quaddata->nbilinexprterms; ++c )
      {
         SCIPinfoMessage(scip, NULL, "%g * ", quaddata->bilinexprterms[c].coef);
         SCIP_CALL( SCIPprintExpr(scip, quaddata->bilinexprterms[c].expr1, NULL) );
         SCIPinfoMessage(scip, NULL, " * ");
         SCIP_CALL( SCIPprintExpr(scip, quaddata->bilinexprterms[c].expr2, NULL) );
         if( c < quaddata->nbilinexprterms-1 )
            SCIPinfoMessage(scip, NULL, " + ");
      }
      SCIPinfoMessage(scip, NULL, "\n");
      SCIPinfoMessage(scip, NULL, "Bilinear of quadratics: \n");
      for( c = 0; c < quaddata->nquadexprs; ++c )
      {
         int i;
         SCIPinfoMessage(scip, NULL, "  For ");
         SCIP_CALL( SCIPprintExpr(scip, quaddata->quadexprterms[c].expr, NULL) );
         SCIPinfoMessage(scip, NULL, " we see: ");
         for( i = 0; i < quaddata->quadexprterms[c].nadjbilin; ++i )
         {
            SCIPinfoMessage(scip, NULL, "%g * ", quaddata->bilinexprterms[quaddata->quadexprterms[c].adjbilin[i]].coef);
            SCIP_CALL( SCIPprintExpr(scip, quaddata->bilinexprterms[quaddata->quadexprterms[c].adjbilin[i]].expr1, NULL) );
            SCIPinfoMessage(scip, NULL, " * ");
            SCIP_CALL( SCIPprintExpr(scip, quaddata->bilinexprterms[quaddata->quadexprterms[c].adjbilin[i]].expr2, NULL) );
            if( i < quaddata->quadexprterms[c].nadjbilin - 1 )
               SCIPinfoMessage(scip, NULL, " + ");
         }
         SCIPinfoMessage(scip, NULL, "\n");
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "Bilinear: none\n");
   }

   return SCIP_OKAY;
}

/** Checks the curvature of the quadratic function, x^T Q x + b^T x stored in quaddata
 *
 * For this, it builds the matrix Q and computes its eigenvalues using LAPACK; if Q is
 * - semidefinite positive -> provided is set to sepaunder
 * - semidefinite negative -> provided is set to sepaover
 * - otherwise -> provided is set to none
 *
 * If assumevarfixed is given and some entries of x correspond to variables present in
 * this hashmap, then the corresponding rows and columns are ignored in the matrix Q.
 */
SCIP_RETCODE SCIPcomputeExprQuadraticCurvature(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPR*            expr,               /**< quadratic expression */
   SCIP_EXPRCURV*        curv,               /**< pointer to store the curvature of quadratics */
   SCIP_HASHMAP*         assumevarfixed,     /**< hashmap containing variables that should be assumed to be fixed, or NULL */
   SCIP_Bool             storeeigeninfo      /**< whether the eigenvalues and eigenvectors should be stored */
   )
{
   SCIP_QUADEXPR* quaddata;
   SCIP_HASHMAP* expr2matrix;
   double* matrix;
   double* alleigval;
   int nvars;
   int nn;
   int n;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(curv != NULL);

   quaddata = expr->quaddata;
   assert(quaddata != NULL);

   /* do not store eigen information if we are not considering full matrix */
   if( assumevarfixed != NULL )
      storeeigeninfo = FALSE;

   if( quaddata->eigeninfostored || (quaddata->curvaturechecked && ! storeeigeninfo) )
   {
      *curv = quaddata->curvature;
      /* if we are convex or concave on the full set of variables, then we will also be so on a subset */
      if( assumevarfixed == NULL || quaddata->curvature != SCIP_EXPRCURV_UNKNOWN )
         return SCIP_OKAY;
   }
   assert(quaddata->curvature == SCIP_EXPRCURV_UNKNOWN || assumevarfixed != NULL || (storeeigeninfo && !quaddata->eigeninfostored));

   *curv = SCIP_EXPRCURV_UNKNOWN;

   n  = quaddata->nquadexprs;

   /* do not check curvature if nn will be too large
    * we want nn * sizeof(real) to fit into an unsigned int, so n must be <= sqrt(unit_max/sizeof(real))
    * sqrt(2*214748364/8) = 7327.1475350234
    */
   if( n > 7000 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "number of quadratic variables is too large (%d) to check the curvature\n", n);
      return SCIP_OKAY;
   }

   /* TODO do some simple tests first; like diagonal entries don't change sign, etc */

   if( !SCIPisIpoptAvailableIpopt() )
      return SCIP_OKAY;

   nn = n * n;
   assert(nn > 0);
   assert((unsigned)nn < UINT_MAX / sizeof(SCIP_Real));

   if( storeeigeninfo )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &quaddata->eigenvalues, n));
      SCIP_CALL( SCIPallocClearBlockMemoryArray(scip, &quaddata->eigenvectors, nn));

      alleigval = quaddata->eigenvalues;
      matrix = quaddata->eigenvectors;
   }
   else
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &alleigval, n) );
      SCIP_CALL( SCIPallocClearBufferArray(scip, &matrix, nn) );
   }


   SCIP_CALL( SCIPhashmapCreate(&expr2matrix, SCIPblkmem(scip), n) );

   /* fill matrix's diagonal */
   nvars = 0;
   for( i = 0; i < n; ++i )
   {
      SCIP_QUADEXPR_QUADTERM quadexprterm;

      quadexprterm = quaddata->quadexprterms[i];

      assert(!SCIPhashmapExists(expr2matrix, (void*)quadexprterm.expr));

      /* skip expr if it is a variable mentioned in assumevarfixed */
      if( assumevarfixed != NULL && SCIPexprIsVar(scip->set, quadexprterm.expr) && SCIPhashmapExists(assumevarfixed, (void*)SCIPgetConsExprExprVarVar(quadexprterm.expr)) )
         continue;

      if( quadexprterm.sqrcoef == 0.0 && ! storeeigeninfo )
      {
         assert(quadexprterm.nadjbilin > 0);
         /* SCIPdebugMsg(scip, "var <%s> appears in bilinear term but is not squared --> indefinite quadratic\n", SCIPvarGetName(quadexprterm.var)); */
         goto CLEANUP;
      }

      matrix[nvars * n + nvars] = quadexprterm.sqrcoef;

      /* remember row of variable in matrix */
      SCIP_CALL( SCIPhashmapInsert(expr2matrix, (void *)quadexprterm.expr, (void *)(size_t)nvars) );
      nvars++;
   }

   /* fill matrix's upper-diagonal */
   for( i = 0; i < quaddata->nbilinexprterms; ++i )
   {
      SCIP_QUADEXPR_BILINTERM bilinexprterm;
      int col;
      int row;

      bilinexprterm = quaddata->bilinexprterms[i];

      /* each factor should have been added to expr2matrix unless it corresponds to a variable mentioned in assumevarfixed */
      assert(SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr1) || (assumevarfixed != NULL && SCIPisExprVar(bilinexprterm.expr1) && SCIPhashmapExists(assumevarfixed, (void*)SCIPgetConsExprExprVarVar(bilinexprterm.expr1))));
      assert(SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr2) || (assumevarfixed != NULL && SCIPisExprVar(bilinexprterm.expr2) && SCIPhashmapExists(assumevarfixed, (void*)SCIPgetConsExprExprVarVar(bilinexprterm.expr2))));

      /* skip bilinear terms where at least one of the factors should be assumed to be fixed (i.e., not present in expr2matrix map) */
      if( !SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr1) || !SCIPhashmapExists(expr2matrix, (void*)bilinexprterm.expr2) )
         continue;

      row = (int)(size_t)SCIPhashmapGetImage(expr2matrix, bilinexprterm.expr1);
      col = (int)(size_t)SCIPhashmapGetImage(expr2matrix, bilinexprterm.expr2);

      assert(row != col);

      if( row < col )
         matrix[row * n + col] = bilinexprterm.coef / 2.0;
      else
         matrix[col * n + row] = bilinexprterm.coef / 2.0;
   }

   /* compute eigenvalues */
   if( LapackDsyev(storeeigeninfo, n, matrix, alleigval) != SCIP_OKAY )
   {
      SCIPwarningMessage(scip, "Failed to compute eigenvalues of quadratic coefficient matrix --> don't know curvature\n");
      goto CLEANUP;
   }

   /* check convexity */
   if( !SCIPisNegative(scip, alleigval[0]) )
      *curv = SCIP_EXPRCURV_CONVEX;
   else if( !SCIPisPositive(scip, alleigval[n-1]) )
      *curv = SCIP_EXPRCURV_CONCAVE;

CLEANUP:
   SCIPhashmapFree(&expr2matrix);

   if( ! storeeigeninfo )
   {
      SCIPfreeBufferArray(scip, &matrix);
      SCIPfreeBufferArray(scip, &alleigval);
   }
   else
   {
      assert(!quaddata->eigeninfostored);
      quaddata->eigeninfostored = TRUE;
   }

   /* if checked convexity on full Q matrix, then remember it
    * if indefinite on submatrix, then it will also be indefinite on full matrix, so can remember that, too */
   if( assumevarfixed == NULL || (*curv == SCIP_EXPRCURV_UNKNOWN) )
   {
      quaddata->curvature = *curv;
      quaddata->curvaturechecked = TRUE;
   }

   return SCIP_OKAY;
}

/**@} */
