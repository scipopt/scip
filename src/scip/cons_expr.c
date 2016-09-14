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

/**@file   cons_expr.c
 * @brief  constraint handler for expression constraints (in particular, nonlinear constraints)
 * @author Stefan Vigerske
 * @author Benjamin Mueller
 * @author Felipe Serrano
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_expr.h"
#include "scip/struct_cons_expr.h"
#include "scip/cons_expr_var.h"
#include "scip/cons_expr_value.h"
#include "scip/cons_expr_sumprod.h"
#include "scip/cons_expr_exp.h"
#include "scip/cons_expr_log.h"
#include "scip/cons_expr_abs.h"
#include "scip/cons_expr_pow.h"
#include "scip/debug.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "expr"
#define CONSHDLR_DESC          "constraint handler for expressions"
#define CONSHDLR_ENFOPRIORITY       -60 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY -4000010 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

/* optional constraint handler properties */
#define CONSHDLR_SEPAPRIORITY        10 /**< priority of the constraint handler for separation */
#define CONSHDLR_SEPAFREQ             1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */

#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING     SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler*/

#define CONSHDLR_PRESOLTIMING    SCIP_PRESOLTIMING_ALWAYS /**< presolving timing of the constraint handler (fast, medium, or exhaustive) */
#define CONSHDLR_MAXPREROUNDS        -1 /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */



/* enable nonlinear constraint upgrading */
#include "scip/cons_nonlinear.h"
#define NONLINCONSUPGD_PRIORITY   100000 /**< priority of the constraint handler for upgrading of nonlinear constraints */

/* enable quadratic constraint upgrading */
#include "scip/cons_quadratic.h"
#define QUADCONSUPGD_PRIORITY     100000 /**< priority of the constraint handler for upgrading of quadratic constraints */



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

/** eventdata for variable bound change events in constraints */
typedef struct
{
   SCIP_CONS*            cons;               /**< constraint */
   SCIP_CONSEXPR_EXPR*   varexpr;            /**< variable expression */
   int                   filterpos;          /**< position of eventdata in SCIP's event filter */
} SCIP_VAREVENTDATA;

/** constraint data for expr constraints */
struct SCIP_ConsData
{
   SCIP_CONSEXPR_EXPR**  varexprs;           /**< array containing all variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */
   SCIP_VAREVENTDATA**   vareventdata;       /**< array containing eventdata for bound change of variables */

   SCIP_CONSEXPR_EXPR*   expr;               /**< expression that represents this constraint */
   SCIP_Real             lhs;                /**< left-hand side */
   SCIP_Real             rhs;                /**< right-hand side */

   SCIP_Real             lhsviol;            /**< violation of left-hand side by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of right-hand side by current solution (used temporarily inside constraint handler) */

   unsigned int          ispropagated:1;     /**< did we propagate the current bounds already? */

   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */

   int                   nlockspos;          /**< number of positive locks */
   int                   nlocksneg;          /**< number of negative locks */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_CONSEXPR_EXPRHDLR** exprhdlrs;       /**< expression handlers */
   int                      nexprhdlrs;      /**< number of expression handlers */
   int                      exprhdlrssize;   /**< size of exprhdlrs array */

   SCIP_CONSEXPR_EXPRHDLR*  exprvarhdlr;     /**< variable expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprvalhdlr;     /**< value expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprsumhdlr;     /**< summation expression handler */
   SCIP_CONSEXPR_EXPRHDLR*  exprprodhdlr;    /**< product expression handler */

   int                      auxvarid;        /**< unique id for the next auxiliary variable */

   unsigned int             lastsoltag;      /**< last solution tag used to evaluate current solution */
   unsigned int             lastsepatag;     /**< last separation tag used to compute cuts */
   unsigned int             lastinitsepatag; /**< last separation initialization flag used */


   SCIP_EVENTHDLR*          eventhdlr;       /**< handler for variable bound change events */

   int                      maxproprounds;   /**< limit on number of propagation rounds for a set of constraints within one round of SCIP propagation */

   /* separation parameters */
   SCIP_Real                mincutviolationsepa;    /**< minimal violation of a cut in order to add it to relaxation during separation */
   SCIP_Real                mincutviolationenfofac; /**< minimal target violation of a cut in order to add it to relaxation during enforcement as factor of feasibility tolerance (may be ignored) */

   unsigned int             restart;         /**< whether we are running after a restart */
};

/** data passed on during expression evaluation in a point */
typedef struct
{
   SCIP_SOL*             sol;                /**< solution that is evaluated */
   unsigned int          soltag;             /**< solution tag */
   SCIP_Bool             aborted;            /**< whether the evaluation has been aborted due to an evaluation error */
} EXPREVAL_DATA;

/** data passed on during expression evaluation over a box */
typedef struct
{
   unsigned int          boxtag;             /**< box tag */
   SCIP_Bool             aborted;            /**< whether the evaluation has been aborted due to an empty interval */
   SCIP_Bool             intersect;          /**< should the computed expression interval be intersected with the existing one? */
   SCIP_Real             varboundrelax;      /**< by how much to relax variable bounds (at most) */
} EXPRINTEVAL_DATA;

/** data passed on during variable locking  */
typedef struct
{
   SCIP_CONSEXPR_EXPRHDLR* exprvarhdlr;      /**< handler for variable expressions (to recognize variable expressions) */
   int                     nlockspos;        /**< number of positive locks */
   int                     nlocksneg;        /**< number of negative locks */
} EXPRLOCK_DATA;

/** data passed on during collecting all expression variables */
typedef struct
{
   SCIP_CONSEXPR_EXPR**  varexprs;           /**< array to store variable expressions */
   int                   nvarexprs;          /**< total number of variable expressions */
   SCIP_HASHMAP*         varexprsmap;        /**< map containing all visited variable expressions */
} GETVARS_DATA;

/** data passed on during copying expressions */
typedef struct
{
   SCIP*                   targetscip;                 /**< target SCIP pointer */
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA_MAPVAR((*mapvar));  /**< variable mapping function, or NULL for identity mapping (used in handler for var-expressions) */
   void*                   mapvardata;                 /**< data of variable mapping function */
   SCIP_CONSEXPR_EXPR*     targetexpr;                 /**< pointer that will hold the copied expression after the walk */
} COPY_DATA;

/** variable mapping data passed on during copying expressions when copying SCIP instances */
typedef struct
{
   SCIP_HASHMAP*           varmap;           /**< SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP */
   SCIP_HASHMAP*           consmap;          /**< SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP */
   SCIP_Bool               global;           /**< should a global or a local copy be created */
   SCIP_Bool               valid;            /**< indicates whether every variable copy was valid */
} COPY_MAPVAR_DATA;

/** data passed on during separation initialization */
typedef struct
{
   SCIP_CONSHDLR*          conshdlr;         /**< expression constraint handler */
   SCIP_Bool               infeasible;       /**< pointer to store whether the problem is infeasible */
   unsigned int            initsepatag;      /**< tag used for the separation initialization round */
} INITSEPA_DATA;

/** data passed on during separation */
typedef struct
{
   SCIP_CONSHDLR*          conshdlr;         /**< expression constraint handler */
   SCIP_SOL*               sol;              /**< solution to separate (NULL for separating the LP solution) */
   SCIP_Real               minviolation;     /**< minimal violation of a cut if it should be added to the LP */
   SCIP_RESULT             result;           /**< buffer to store a result */
   int                     ncuts;            /**< buffer to store the total number of added cuts */
   unsigned int            sepatag;          /**< separation tag */
} SEPA_DATA;

struct SCIP_ConsExpr_PrintDotData
{
   FILE*                   file;             /**< file to print to */
   SCIP_Bool               closefile;        /**< whether file need to be closed when finished printing */
   SCIP_HASHMAP*           visitedexprs;     /**< hashmap storing expressions that have been printed already */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint;  /**< flags that indicate what to print for each expression */
};

/** data passed on during creating of auxiliary variables */
typedef struct
{
   SCIP_CONSHDLR*          conshdlr;         /**< expression constraint handler */
#ifdef SCIP_DEBUG_SOLUTION
   SCIP_SOL*               debugsol;         /**< debug solution (or NULL if not debugging) */
#endif
} CREATE_AUXVARS_DATA;


/*
 * Local methods
 */

/** create and include conshdlr to SCIP and set everything except for expression handlers */
static
SCIP_RETCODE includeConshdlrExprBasic(SCIP* scip);

/** copy expression handlers from sourceconshdlr to (target's) scip consexprhdlr */
static
SCIP_RETCODE copyConshdlrExprExprHdlr(
   SCIP*                 scip,               /**< (target) SCIP data structure */
   SCIP_CONSHDLR*        sourceconshdlr,     /**< source constraint expression handler */
   SCIP_Bool*            valid               /**< was the copying process valid? */
   )
{
   int                i;
   SCIP_CONSHDLR*     conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLRDATA* sourceconshdlrdata;

   assert(strcmp(SCIPconshdlrGetName(sourceconshdlr), CONSHDLR_NAME) == 0);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   assert(conshdlr != sourceconshdlr);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   sourceconshdlrdata = SCIPconshdlrGetData(sourceconshdlr);
   assert(sourceconshdlrdata != NULL);

   /* copy expression handlers */
   *valid = TRUE;
   for( i = 0; i < sourceconshdlrdata->nexprhdlrs; i++ )
   {
      SCIP_Bool localvalid;
      SCIP_CONSEXPR_EXPRHDLR* sourceexprhdlr;

      sourceexprhdlr = sourceconshdlrdata->exprhdlrs[i];

      if( sourceexprhdlr->copyhdlr != NULL )
      {
         SCIP_CALL( sourceexprhdlr->copyhdlr(scip, conshdlr, sourceconshdlr, sourceexprhdlr, &localvalid) );
         *valid &= localvalid;
      }
      else
      {
         *valid = FALSE;
      }
   }

   /* set pointer to important expression handlers in conshdlr of target SCIP */
   conshdlrdata->exprvarhdlr = SCIPfindConsExprExprHdlr(conshdlr, "var");
   conshdlrdata->exprvalhdlr = SCIPfindConsExprExprHdlr(conshdlr, "val");
   conshdlrdata->exprsumhdlr = SCIPfindConsExprExprHdlr(conshdlr, "sum");
   conshdlrdata->exprprodhdlr = SCIPfindConsExprExprHdlr(conshdlr, "prod");

   return SCIP_OKAY;
}

/** returns an equivalent expression for a given expression if possible; it adds the expression to key2expr if the map
 *  does not contain the key
 */
static
SCIP_RETCODE findEqualExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR *  expr,               /**< expression to replace */
   SCIP_HASHTABLE*       key2expr,           /**< mapping of hashes to expressions */
   SCIP_CONSEXPR_EXPR**  newexpr             /**< pointer to store an equivalent expression (NULL if there is none) */
   )
{  /*lint --e{438}*/
   SCIP_HASHTABLELIST* hashtablelist;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(key2expr != NULL);
   assert(newexpr != NULL);

   *newexpr = NULL;
   hashtablelist = NULL;

   do
   {
      /* search for an equivalent expression */
      *newexpr = (SCIP_CONSEXPR_EXPR*)(SCIPhashtableRetrieveNext(key2expr, &hashtablelist, (void*)expr));

      if( *newexpr == NULL )
      {
         /* processed all expressions like expr from hash table, so insert expr */
         SCIP_CALL( SCIPhashtableInsert(key2expr, (void*) expr) );
         break;
      }
      else if( expr != *newexpr )
      {
         assert(SCIPcompareConsExprExprs(expr, *newexpr) == 0);
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

/** @name Walking methods
 *
 * Several operations need to traverse the whole expression tree: print, evaluate, free, etc.
 * These operations have a very natural recursive implementation. However, deep recursion can raise stack overflows.
 * To avoid this issue, the method SCIPwalkConsExprExprDF is introduced to traverse the tree and execute callbacks
 * at different places. Here are the callbacks needed for performing the mentioned operations.
 *
 * @{
 */

/** expression walk callback to copy an expression
 *
 * In expr->walkio is given the targetexpr which is expected to hold the copy of expr.
 */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(copyExpr)
{
   assert(expr != NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* create expr that will hold the copy */
         SCIP_CONSEXPR_EXPR*     targetexpr;
         SCIP_CONSEXPR_EXPRHDLR* targetexprhdlr;
         SCIP_CONSEXPR_EXPRDATA* targetexprdata;
         COPY_DATA* copydata;

         copydata = (COPY_DATA*)data;

         /* get the exprhdlr of the target scip */
         if( copydata->targetscip != scip )
         {
            SCIP_CONSHDLR* targetconsexprhdlr;

            targetconsexprhdlr = SCIPfindConshdlr(copydata->targetscip, "expr");
            assert(targetconsexprhdlr != NULL);

            targetexprhdlr = SCIPfindConsExprExprHdlr(targetconsexprhdlr,
                  SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)));

            if( targetexprhdlr == NULL )
            {
               /* expression handler not in target scip (probably did not have a copy callback) -> abort */
               expr->walkio.ptrval = NULL;
               *result = SCIP_CONSEXPREXPRWALK_SKIP;
               return SCIP_OKAY;
            }
         }
         else
         {
            targetexprhdlr = SCIPgetConsExprExprHdlr(expr);
         }
         assert(targetexprhdlr != NULL);

         /* if the source is a variable expression create a variable expression directly; otherwise copy the expression data */
         if( strcmp(expr->exprhdlr->name, "var") == 0 )
         {
            SCIP_VAR* sourcevar;
            SCIP_VAR* targetvar;

            sourcevar = SCIPgetConsExprExprVarVar(expr);
            assert(sourcevar != NULL);
            targetvar = NULL;

            /* get the corresponding variable in the target SCIP */
            if( copydata->mapvar != NULL )
            {
               SCIP_CALL( copydata->mapvar(copydata->targetscip, &targetvar, scip, sourcevar, copydata->mapvardata) );
               SCIP_CALL( SCIPcreateConsExprExprVar(copydata->targetscip, SCIPfindConshdlr(copydata->targetscip, "expr"), &targetexpr, targetvar) );

               /* we need to release once since it has been captured by the mapvar() and SCIPcreateConsExprExprVar() call */
               SCIP_CALL( SCIPreleaseVar(copydata->targetscip, &targetvar) );
            }
            else
            {
               targetvar = sourcevar;
               SCIP_CALL( SCIPcreateConsExprExprVar(copydata->targetscip, SCIPfindConshdlr(copydata->targetscip, "expr"), &targetexpr, targetvar) );
            }
         }
         else
         {
            /* copy expression data */
            if( expr->exprhdlr->copydata != NULL )
            {
               SCIP_CALL( expr->exprhdlr->copydata(
                     copydata->targetscip,
                     targetexprhdlr,
                     &targetexprdata,
                     scip,
                     expr,
                     copydata->mapvar,
                     copydata->mapvardata) );
            }
            else if( expr->exprdata != NULL )
            {
               /* no copy callback for expression data implemented -> abort
                * (we could also just copy the exprdata pointer, but for now let's say that
                *  an expression handler should explicitly implement this behavior, if desired)
                */
               expr->walkio.ptrval = NULL;
               *result = SCIP_CONSEXPREXPRWALK_SKIP;
               return SCIP_OKAY;
            }
            else
            {
               targetexprdata = NULL;
            }

            /* create in targetexpr an expression of the same type as expr, but without children for now */
            SCIP_CALL( SCIPcreateConsExprExpr(copydata->targetscip, &targetexpr, targetexprhdlr, targetexprdata, 0, NULL) );
         }

         /* store targetexpr */
         expr->walkio.ptrval = targetexpr;

         /* continue */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }


      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      {
         /* just visited child so a copy of himself should be available; append it */
         SCIP_CONSEXPR_EXPR* child;
         SCIP_CONSEXPR_EXPR* targetchild;
         SCIP_CONSEXPR_EXPR* targetexpr;
         COPY_DATA* copydata;

         assert(expr->walkcurrentchild < expr->nchildren);

         child = expr->children[expr->walkcurrentchild];
         copydata = (COPY_DATA*)data;

         /* get copy of child */
         targetchild = (SCIP_CONSEXPR_EXPR*)child->walkio.ptrval;

         if( targetchild == NULL )
         {
            /* release targetexpr (should free also the already copied children) */
            SCIP_CALL( SCIPreleaseConsExprExpr(copydata->targetscip, (SCIP_CONSEXPR_EXPR**)&expr->walkio.ptrval) );

            /* abort */
            *result = SCIP_CONSEXPREXPRWALK_SKIP;
            return SCIP_OKAY;
         }

         /* append child to copyexpr */
         targetexpr = (SCIP_CONSEXPR_EXPR*)expr->walkio.ptrval;
         SCIP_CALL( SCIPappendConsExprExpr(copydata->targetscip, targetexpr, targetchild) );

         /* release targetchild (captured by targetexpr) */
         SCIP_CALL( SCIPreleaseConsExprExpr(copydata->targetscip, &targetchild) );

         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      {
         COPY_DATA* copydata;

         /* store the copied expression in the copydata, in case this expression was the root, for which the walkio will be cleared */
         copydata = (COPY_DATA*)data;
         copydata->targetexpr = (SCIP_CONSEXPR_EXPR*)expr->walkio.ptrval;

         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      default:
      {
         SCIPABORT(); /* we should never be called in this stage */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE; /*lint !e527*/
         return SCIP_OKAY;
      }
   }
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA_MAPVAR(transformVar)
{   /*lint --e{715}*/
   assert(sourcevar != NULL);
   assert(targetvar != NULL);
   assert(sourcescip == targetscip);

   /* transform variable (does not capture target variable) */
   SCIP_CALL( SCIPgetTransformedVar(sourcescip, sourcevar, targetvar) );
   assert(*targetvar != NULL);

   /* caller assumes that target variable has been captured */
   SCIP_CALL( SCIPcaptureVar(sourcescip, *targetvar) );

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA_MAPVAR(copyVar)
{
   COPY_MAPVAR_DATA* data;
   SCIP_Bool valid;

   assert(sourcevar != NULL);
   assert(targetvar != NULL);
   assert(mapvardata != NULL);

   data = (COPY_MAPVAR_DATA*)mapvardata;

   SCIP_CALL( SCIPgetVarCopy(sourcescip, targetscip, sourcevar, targetvar, data->varmap, data->consmap, data->global, &valid) );
   assert(*targetvar != NULL);

   /* if copy was not valid, store so in mapvar data */
   if( !valid )
      data->valid = FALSE;

   /* caller assumes that target variable has been captured */
   SCIP_CALL( SCIPcaptureVar(targetscip, *targetvar) );

   return SCIP_OKAY;
}

/** expression walk callback to free an expression including its children (if not used anywhere else) */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(freeExpr)
{  /*lint --e{715}*/
   assert(expr != NULL);

   /* expression should be used by its parent and maybe by the walker (only the root!)
    * in VISITEDCHILD we assert that expression is only used by its parent
    */
   assert(0 <= expr->nuses && expr->nuses <= 2);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      {
         /* check whether a child needs to be visited (nuses == 1)
          * if not, then we still have to release it
          */
         SCIP_CONSEXPR_EXPR* child;

         assert(expr->walkcurrentchild < expr->nchildren);
         child = expr->children[expr->walkcurrentchild];
         if( child->nuses > 1 )
         {
            /* child is not going to be freed: just release it */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &child) );
            *result = SCIP_CONSEXPREXPRWALK_SKIP;
         }
         else
         {
            assert(child->nuses == 1);

            if( child->exprdata != NULL )
            {
               /* free child's expression data when entering child */
               if( child->exprhdlr->freedata != NULL )
               {
                  SCIP_CALL( child->exprhdlr->freedata(scip, child) );
                  assert(child->exprdata == NULL);
               }
               else
               {
                  child->exprdata = NULL;
               }
            }

            *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         }

         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      {
         /* free child after visiting it */
         SCIP_CONSEXPR_EXPR* child;

         assert(expr->walkcurrentchild < expr->nchildren);

         child = expr->children[expr->walkcurrentchild];
         /* child should only be used by its parent */
         assert(child->nuses == 1);

         /* child should have no data associated */
         assert(child->exprdata == NULL);

         /* free child's children array, if any */
         SCIPfreeBlockMemoryArrayNull(scip, &child->children, child->childrensize);

         /* expression should not be locked anymore */
         assert(child->nlockspos == 0);
         assert(child->nlocksneg == 0);

         SCIPfreeBlockMemory(scip, &child);
         assert(child == NULL);
         expr->children[expr->walkcurrentchild] = NULL;

         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      default:
      {
         SCIPABORT(); /* we should never be called in this stage */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE; /*lint !e527*/
         return SCIP_OKAY;
      }
   }
}

/** expression walk callback to print an expression */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(printExpr)
{
   FILE* file;

   assert(expr != NULL);
   assert(expr->exprhdlr != NULL);

   file = (FILE*)data;

   if( expr->exprhdlr->print == NULL )
   {
      /* default: <hdlrname>(<child1>, <child2>, ...) */
      switch( stage )
      {
         case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
         {
            SCIPinfoMessage(scip, file, SCIPgetConsExprExprHdlrName(expr->exprhdlr));
            if( expr->nchildren > 0 )
            {
               SCIPinfoMessage(scip, file, "(");
            }
            break;
         }

         case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
         {
            if( SCIPgetConsExprExprWalkCurrentChild(expr) < expr->nchildren-1 )
            {
               SCIPinfoMessage(scip, file, ", ");
            }
            else
            {
               SCIPinfoMessage(scip, file, ")");
            }

            break;
         }

         case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
         case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
         default: ;
      }
   }
   else
   {
      /* redirect to expression callback */
      SCIP_CALL( (*expr->exprhdlr->print)(scip, expr, stage, file) );
   }

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

/** expression walk callback to print an expression in dot format */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(printExprDot)
{  /*lint --e{715}*/
   SCIP_CONSEXPR_PRINTDOTDATA* dotdata;
   SCIP_CONSEXPR_EXPR* parentbackup;
   SCIP_Real color;
   int c;

   assert(expr != NULL);
   assert(expr->exprhdlr != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);
   assert(data != NULL);

   dotdata = (SCIP_CONSEXPR_PRINTDOTDATA*)data;

   /* skip expressions that have been printed already */
   if( SCIPhashmapExists(dotdata->visitedexprs, (void*)expr) )
   {
      *result = SCIP_CONSEXPREXPRWALK_SKIP;
      return SCIP_OKAY;
   }

   /* print expression as dot node */

   /* make up some color from the expression type (it's name) */
   color = 0.0;
   for( c = 0; expr->exprhdlr->name[c] != '\0'; ++c )
      color += (tolower(expr->exprhdlr->name[c]) - 'a') / 26.0;
   color = SCIPfrac(scip, color);
   SCIPinfoMessage(scip, dotdata->file, "n%p [fillcolor=\"%g,%g,%g\", label=\"", expr, color, color, color);

   if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_EXPRHDLR )
   {
      SCIPinfoMessage(scip, dotdata->file, "%s\\n", SCIPgetConsExprExprHdlrName(expr->exprhdlr));
   }

   if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_EXPRSTRING )
   {
      /* print expression string as label */
      parentbackup = expr->walkparent;
      expr->walkparent = NULL;
      assert(expr->walkcurrentchild == 0); /* as we are in enterexpr */

      SCIP_CALL( printExpr(scip, expr, SCIP_CONSEXPREXPRWALK_ENTEREXPR, (void*)dotdata->file, result) );
      for( c = 0; c < expr->nchildren; ++c )
      {
         expr->walkcurrentchild = c;
         SCIP_CALL( printExpr(scip, expr, SCIP_CONSEXPREXPRWALK_VISITINGCHILD, (void*)dotdata->file, result) );
         SCIPinfoMessage(scip, dotdata->file, "c%d", c);
         SCIP_CALL( printExpr(scip, expr, SCIP_CONSEXPREXPRWALK_VISITEDCHILD, (void*)dotdata->file, result) );
      }
      SCIP_CALL( printExpr(scip, expr, SCIP_CONSEXPREXPRWALK_LEAVEEXPR, (void*)dotdata->file, result) );
      SCIPinfoMessage(scip, dotdata->file, "\\n");

      expr->walkcurrentchild = 0;
      expr->walkparent = parentbackup;
   }

   if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_NUSES )
   {
      /* print number of uses */
      SCIPinfoMessage(scip, dotdata->file, "%d uses\\n", expr->nuses);
   }

   if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_NUSES )
   {
      /* print number of locks */
      SCIPinfoMessage(scip, dotdata->file, "%d,%d +,-locks\\n", expr->nlockspos, expr->nlocksneg);
   }

   if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_EVALVALUE )
   {
      /* print eval value */
      SCIPinfoMessage(scip, dotdata->file, "val=%g", expr->evalvalue);

      if( (dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_EVALTAG) == SCIP_CONSEXPR_PRINTDOT_EVALTAG )
      {
         /* print also eval tag */
         SCIPinfoMessage(scip, dotdata->file, " (%u)", expr->evaltag);
      }
      SCIPinfoMessage(scip, dotdata->file, "\\n");
   }

   if( dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_INTERVAL )
   {
      /* print interval value */
      SCIPinfoMessage(scip, dotdata->file, "[%g,%g]", expr->interval.inf, expr->interval.sup);

      if( (dotdata->whattoprint & SCIP_CONSEXPR_PRINTDOT_INTERVALTAG) == SCIP_CONSEXPR_PRINTDOT_INTERVALTAG )
      {
         /* print also interval eval tag */
         SCIPinfoMessage(scip, dotdata->file, " (%u)", expr->intevaltag);
      }
      SCIPinfoMessage(scip, dotdata->file, "\\n");
   }

   SCIPinfoMessage(scip, dotdata->file, "\"]\n");  /* end of label and end of node */

   /* add edges from expr to its children */
   for( c = 0; c < expr->nchildren; ++c )
      SCIPinfoMessage(scip, dotdata->file, "n%p -> n%p [label=\"c%d\"]\n", (void*)expr, (void*)expr->children[c], c);

   /* remember that we have printed this expression */
   SCIP_CALL( SCIPhashmapInsert(dotdata->visitedexprs, (void*)expr, NULL) );

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

/** expression walk callback when evaluating expression, called before child is visited */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(evalExprVisitChild)
{  /*lint --e{715}*/
   EXPREVAL_DATA* evaldata;

   assert(expr != NULL);
   assert(data != NULL);

   evaldata = (EXPREVAL_DATA*)data;

   /* skip child if it has been evaluated for that solution already */
   if( evaldata->soltag != 0 && evaldata->soltag == expr->children[expr->walkcurrentchild]->evaltag )
   {
      if( expr->children[expr->walkcurrentchild]->evalvalue == SCIP_INVALID ) /*lint !e777*/
      {
         evaldata->aborted = TRUE;
         *result = SCIP_CONSEXPREXPRWALK_ABORT;
      }
      else
      {
         *result = SCIP_CONSEXPREXPRWALK_SKIP;
      }
   }
   else
   {
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   }

   return SCIP_OKAY;
}

/** expression walk callback when evaluating expression, called when expression is left */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(evalExprLeaveExpr)
{  /*lint --e{715}*/
   EXPREVAL_DATA* evaldata;

   assert(expr != NULL);
   assert(data != NULL);
   assert(expr->exprhdlr->eval != NULL);

   evaldata = (EXPREVAL_DATA*)data;

   SCIP_CALL( (*expr->exprhdlr->eval)(scip, expr, &expr->evalvalue, evaldata->sol) );
   expr->evaltag = evaldata->soltag;

   if( expr->evalvalue == SCIP_INVALID ) /*lint !e777*/
   {
      evaldata->aborted = TRUE;
      *result = SCIP_CONSEXPREXPRWALK_ABORT;
   }
   else
   {
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   }

   return SCIP_OKAY;
}

/** expression walk callback when evaluating expression on intervals, called before child is visited */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(intevalExprVisitChild)
{  /*lint --e{715}*/
   EXPRINTEVAL_DATA* propdata;

   assert(expr != NULL);
   assert(data != NULL);

   propdata = (EXPRINTEVAL_DATA*)data;

   /* skip child if it has been evaluated already */
   if( propdata->boxtag != 0 && propdata->boxtag == expr->children[expr->walkcurrentchild]->intevaltag )
   {
      if( SCIPintervalIsEmpty(SCIPinfinity(scip), expr->children[expr->walkcurrentchild]->interval) )
      {
         propdata->aborted = TRUE;
         *result = SCIP_CONSEXPREXPRWALK_ABORT;
      }
      else
      {
         *result = SCIP_CONSEXPREXPRWALK_SKIP;
      }
   }
   else
   {
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   }

   return SCIP_OKAY;
}

/** expression walk callback when evaluating expression on intervals, called when expression is left */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(intevalExprLeaveExpr)
{  /*lint --e{715}*/
   EXPRINTEVAL_DATA* propdata;
   SCIP_INTERVAL interval;

   assert(expr != NULL);
   assert(data != NULL);

   propdata = (EXPRINTEVAL_DATA*)data;

   /* set tag in any case */
   expr->intevaltag = propdata->boxtag;

   /* mark expression as not tightened if we do not intersect expression intervals; this happens onces before calling
    * the reverse propagation
    */
   if( !propdata->intersect )
      expr->hastightened = FALSE;

   /* set interval to [-inf,+inf] if interval evaluation callback is not implemented */
   if( expr->exprhdlr->inteval == NULL )
   {
      SCIPintervalSetEntire(SCIPinfinity(scip), &expr->interval);
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

      return SCIP_OKAY;
   }

   /* evaluate current expression and move on */
   SCIP_CALL( (*expr->exprhdlr->inteval)(scip, expr, &interval, propdata->varboundrelax) );

   /* update expression interval */
   if( !SCIPintervalIsEmpty(SCIPinfinity(scip), interval) )
   {
      /* intersect new interval with the previous one */
      if( propdata->intersect )
         SCIPintervalIntersect(&expr->interval, expr->interval, interval);
      else
         SCIPintervalSetBounds(&expr->interval, interval.inf, interval.sup);

      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   }

   /* stop if the computed or resulting interval is empty */
   if( SCIPintervalIsEmpty(SCIPinfinity(scip), interval) || SCIPintervalIsEmpty(SCIPinfinity(scip), expr->interval) )
   {
      SCIPintervalSetEmpty(&expr->interval);
      propdata->aborted = TRUE;
      *result = SCIP_CONSEXPREXPRWALK_ABORT;
   }

   return SCIP_OKAY;
}

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(lockVar)
{
   EXPRLOCK_DATA* lockdata;

   assert(expr != NULL);
   assert(data != NULL);
   assert(result != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   lockdata = (EXPRLOCK_DATA*)data;

   expr->nlockspos += lockdata->nlockspos;
   expr->nlocksneg += lockdata->nlocksneg;

   if( SCIPgetConsExprExprHdlr(expr) == lockdata->exprvarhdlr )
   {
      /* if a variable, then also add nlockspos/nlocksneg from lockdata via SCIPaddVarLocks() */
      SCIP_CALL( SCIPaddVarLocks(scip, SCIPgetConsExprExprVarVar(expr), lockdata->nlockspos, lockdata->nlocksneg) );
   }

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

/** prints structure a la Maple's dismantle */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(dismantleExpr)
{
   assert(expr != NULL);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR:
      {
         int* depth;
         int nspaces;
         const char* type;

         depth = (int*)data;
         ++*depth;
         nspaces = 3 * *depth;
         type = SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr));

         /* use depth of expression to align output */
         SCIPinfoMessage(scip, NULL, "%*s[%s]: ", nspaces, "", type);

         if(strcmp(type, "var") == 0)
            SCIPinfoMessage(scip, NULL, "%s\n", SCIPvarGetName(SCIPgetConsExprExprVarVar(expr)));
         else if(strcmp(type, "sum") == 0)
            SCIPinfoMessage(scip, NULL, "%g\n", SCIPgetConsExprExprSumConstant(expr));
         else if(strcmp(type, "prod") == 0)
            SCIPinfoMessage(scip, NULL, "%g\n", SCIPgetConsExprExprProductCoef(expr));
         else if(strcmp(type, "val") == 0)
            SCIPinfoMessage(scip, NULL, "%g\n", SCIPgetConsExprExprValueValue(expr));
         else if(strcmp(type, "exp") == 0)
            SCIPinfoMessage(scip, NULL, "\n");
         else if(strcmp(type, "log") == 0)
            SCIPinfoMessage(scip, NULL, "\n");
         else if(strcmp(type, "abs") == 0)
            SCIPinfoMessage(scip, NULL, "\n");
         else
            SCIPinfoMessage(scip, NULL, "NOT IMPLEMENTED YET\n");
         break;
      }
      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD:
      {
         int* depth;
         int nspaces;
         const char* type;

         depth = (int*)data;
         nspaces = 3 * *depth;
         type = SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr));

         if( strcmp(type, "sum") == 0 )
         {
            SCIPinfoMessage(scip, NULL, "%*s   ", nspaces, "");
            SCIPinfoMessage(scip, NULL, "[coef]: %g\n", SCIPgetConsExprExprSumCoefs(expr)[SCIPgetConsExprExprWalkCurrentChild(expr)]);
         }
         else if( strcmp(type, "prod") == 0 )
         {
            SCIPinfoMessage(scip, NULL, "%*s   ", nspaces, "");
            SCIPinfoMessage(scip, NULL, "[expo]: %g\n", SCIPgetConsExprExprProductExponents(expr)[SCIPgetConsExprExprWalkCurrentChild(expr)]);
         }
         break;
      }
      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR:
      {
         int* depth;

         depth = (int*)data;
         --*depth;
         break;
      }
      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD:
      default:
      {
         /* shouldn't be here */
         SCIPABORT();
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE; /*lint !e527*/
         return SCIP_OKAY;
      }
   }
   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

/** expression walk callback to skip expression which have already been hashed */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(hashExprVisitingExpr)
{
   SCIP_HASHMAP* expr2key;
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_VISITINGCHILD);

   expr2key = (SCIP_HASHMAP*) data;
   assert(expr2key != NULL);

   assert(expr->walkcurrentchild < expr->nchildren);
   child = expr->children[expr->walkcurrentchild];
   assert(child != NULL);

   /* skip child if the expression is already in the map */
   *result = SCIPhashmapExists(expr2key, (void*) child) ? SCIP_CONSEXPREXPRWALK_SKIP : SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

/** expression walk callback to compute an hash value for an expression */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(hashExprLeaveExpr)
{
   SCIP_HASHMAP* expr2key;
   unsigned int hashkey;
   int i;

   assert(expr != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_LEAVEEXPR);

   expr2key = (SCIP_HASHMAP*) data;
   assert(expr2key != NULL);

   hashkey = 0;
   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   if( expr->exprhdlr->hash != NULL )
   {
      SCIP_CALL( (*expr->exprhdlr->hash)(scip, expr, expr2key, &hashkey) );
   }
   else
   {
      /* compute hash from expression handler name if callback is not implemented
       * this can lead to more collisions and thus a larger number of expensive expression compare calls
       */
      for( i = 0; expr->exprhdlr->name[i] != '\0'; i++ )
         hashkey += (unsigned int) expr->exprhdlr->name[i]; /*lint !e571*/

      hashkey = SCIPcalcFibHash((SCIP_Real)hashkey);
   }

   /* put the hash key into expr2key map */
   SCIP_CALL( SCIPhashmapInsert(expr2key, (void*)expr, (void*)(size_t)hashkey) );

   return SCIP_OKAY;
}

/** expression walk callback to replace common sub-expression */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(commonExprVisitingExpr)
{
   SCIP_HASHTABLE* key2expr;
   SCIP_CONSEXPR_EXPR* newchild;
   SCIP_CONSEXPR_EXPR* child;

   assert(expr != NULL);
   assert(data != NULL);
   assert(result != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_VISITINGCHILD);

   key2expr = (SCIP_HASHTABLE*)data;
   assert(key2expr != NULL);

   assert(expr->walkcurrentchild < expr->nchildren);
   child = expr->children[expr->walkcurrentchild];
   assert(child != NULL);

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   /* try to find an equivalent expression */
   SCIP_CALL( findEqualExpr(scip, child, key2expr, &newchild) );

   /* replace child with newchild */
   if( newchild != NULL )
   {
      assert(child != newchild);
      assert(SCIPcompareConsExprExprs(child, newchild) == 0);

      SCIPdebugMessage("replacing common child expression %p -> %p\n", (void*)child, (void*)newchild);

      SCIP_CALL( SCIPreplaceConsExprExprChild(scip, expr, expr->walkcurrentchild, newchild) );

      *result = SCIP_CONSEXPREXPRWALK_SKIP;
   }

   return SCIP_OKAY;
}

/** expression walk callback to collect all variable expressions */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(getVarExprsLeaveExpr)
{
   GETVARS_DATA* getvarsdata;

   assert(expr != NULL);
   assert(result != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_LEAVEEXPR);

   getvarsdata = (GETVARS_DATA*) data;
   assert(getvarsdata != NULL);

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   /* add variable expression if not seen so far; there is only one variable expression representing a variable */
   if( strcmp(expr->exprhdlr->name, "var") == 0 && !SCIPhashmapExists(getvarsdata->varexprsmap, (void*) expr) )
   {
      assert(SCIPgetNTotalVars(scip) >= getvarsdata->nvarexprs + 1);

      getvarsdata->varexprs[ getvarsdata->nvarexprs ] = expr;
      assert(getvarsdata->varexprs[getvarsdata->nvarexprs] != NULL);
      ++(getvarsdata->nvarexprs);
      SCIP_CALL( SCIPhashmapInsert(getvarsdata->varexprsmap, (void*) expr, NULL) );

      /* capture expression */
      SCIPcaptureConsExprExpr(expr);
   }

   return SCIP_OKAY;
}

/**@} */  /* end of walking methods */

/** @name Simplifying methods
 *
 * This is largely inspired in Joel Cohen's
 * Computer algebra and symbolic computation: Mathematical methods
 * In particular Chapter 3
 * The other fountain of inspiration is the current simplifying methods in expr.c.
 *
 * Definition of simplified expressions
 * ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 * An expression is simplified if it
 * - is a value expression
 * - is a var expression
 * - is a product expression such that
 *    SP1:  every child is simplified
 *    SP2:  no child with integer exponent is a product
 *    SP3:  no child with integer exponent is a sum with a single term ((2*x)^2 -> 4*x^2)
 *    SP4:  no two children are the same expression (those should be multiplied)
 *    SP5:  the children are sorted [commutative rule]
 *    SP6:  no exponent is 0
 *    SP7:  no child is a value
 *    SP8:  its coefficient is 1.0 (otherwise should be written as sum)
 *    SP9:  if it consists of a single child, then its exponent != 1.0
 *    SP10: it has at least one child
 *    ? at most one child is an exp
 *    ? at most one child is an abs
 * - is a sum expression such that
 *    SS1: every child is simplified
 *    SS2: no child is a sum
 *    SS3: no child is a value (values should go in the constant of the sum)
 *    SS4: no two children are the same expression (those should be summed up)
 *    SS5: the children are sorted [commutative rule]
 *    SS6: it has at least one child
 *    SS7: if it consists of a single child, then either constant is != 0.0 or coef != 1
 *    SS8: no child has coefficient 0
 *    x if it consists of a single child, then its constant != 0.0 (otherwise, should be written as a product)
 * - it is a function with simplified arguments, but not all of them can be values
 * ? a logarithm doesn't have a product as a child
 * ? the exponent of an exponential is always 1
 *
 * ORDERING RULES
 * ^^^^^^^^^^^^^^
 * These rules define a total order on *simplified* expressions.
 * There are two groups of rules, when comparing equal type expressions and different type expressions
 * Equal type expressions:
 * OR1: u,v value expressions: u < v <=> val(u) < val(v)
 * OR2: u,v var expressions: u < v <=> SCIPvarGetIndex(var(u)) < SCIPvarGetIndex(var(v))
 * OR3: u,v are both sum or product expression: < is a lexicographical order on the terms
 * OR4: u,v are u = FUN(u_1, ..., u_n), v = FUN(v_1, ..., v_m): u < v <=> For the first k such that u_k != v_k, u_k < v_k,
 *      or if such a k doesn't exist, then n < m.
 *
 * Different type expressions:
 * OR5: u value, v other: u < v always
 * OR6: u sum, v product, var or func: u < v <=> u < 0+v
 *      In other words, u = \sum_{i = 1}^n \alpha_i u_i, then u < v <=> u_n < v or if u_n = v and \alpha_n < 1
 * OR7: u product, v var or func: u < v <=> u < 1*v
 *      In other words, u = \Pi_{i = 1}^n u_i ^ \alpha_i, then u < v <=> u_n < v or if u_n = v and \alpha_n < 1
 *      @note: since this applies only to simplified expressions, the form of the product is correct. Simplified products
 *             do *not* have constant coefficients
 * OR8: u var, v func: u < v always
 * OR9: u func, v other type of func: u < v <=> name(type(u)) < name(type(v))
 * OR10: none of the rules apply: u < v <=> ! v < u
 * Examples:
 * OR10: x < x^2 ?:  x is var and x^2 product, so none applies.
 *       Hence, we try to answer x^2 < x ?: x^2 < x <=> x < x or if x = x and 2 < 1 <=> 2 < 1 <=> False, so x < x^2 is True
 *
 * Algorithm
 * ^^^^^^^^^
 * The recursive version of the algorithm is
 *
 * EXPR simplify(expr)
 *    for c in 1..expr->nchildren
 *       expr->children[c] = simplify(expr->children[c])
 *    end
 *    return expr->exprhdlr->simplify(expr)
 * end
 *
 * Important: Whatever is returned by a simplify callback **has** to be simplified.
 * Also, all children of the given expression **are** already simplified
 *
 * Here is an outline of the algorithm for simplifying sum expressions:
 * The idea is to create a list of all the children that the simplified expr must containt.
 * We use a linked list to construct it
 *
 * INPUT:  expr  :: sum expression to be simplified
 * OUTPUT: sexpr :: simplified expression
 * NOTE: it *can* modify expr
 *
 * simplified_coefficient <- expr->coefficient
 * expr_list <- empty list (list containing the simplified children of the final simplified expr)
 * For each child in expr->children:
 *    1. if child's coef is 0: continue
 *    2. if child is value: add it to simplified_coefficient and continue
 *    3. if child is not a sum: build list L = [(coef,child)]
 *    4. if child is sum:
 *       4.1. if coef is not 1.0: multiply child by coef (*)
 *       4.2. build list with the children of child, L = [(val, expr) for val in child->coeffs, expr in child->children)]
 *    5. mergeSum(L, expr_list)
 * if expr_list is empty, return value expression with value simplified_coefficient
 * if expr_list has only one child and simplified_coefficient is 1, return child
 * otherwise, build sum expression using the exprs in expr_list as children
 *
 * The method mergeSum simply inserts the elements of L into expr_list. Note that both lists are sorted.
 * While inserting, collisions can happen. A collision means that we have to add the two expressions.
 * However, after adding them, we need to simplify the resulting expression (e.g., the coefficient may become 0.0).
 * Fortunately, the coefficient being 0 is the only case we have to handle.
 * PROOF: all expressions in expr_list are simplified wrt to the sum, meaning that if we would build a sum
 * expression from them, it would yield a simplified sum expression. If there is a collision, then the expression
 * in L has to be in expr_list. The sum yields coef*expr and from the previous one can easily verify that it is
 * a valid child of a simplified sum (it is not a sum, etc), except for the case coef = 0.
 * Note: the context where the proof works is while merging (adding) children. Before this step, the children
 * go through a "local" simplification (i.e., 1-4 above). There, we *do* have to take care of other cases.
 * But, in contrast to products, after this steps, no child in finalchildren is a sum and the proof works.
 *
 * The algorithm for simplifying a product is basically the same. One extra difficulty occurs at (*):
 * The distribution of the exponent over a product children can only happen if the exponent is integral.
 * Also, in that case, the resulting new child could be unsimplified, so it must be re-simplified.
 * While merging, multiplying similar product expressions can render them unsimplified. So to handle it
 * one basically needs to simulate (the new) (*) while merging. Hence, a further merge might be necessary
 * (and then all the book-keeping information to perform the original merge faster is lost)
 *
 * @{
 */

/** expression walk callback to simplify an expression
 * simplifies bottom up; when leaving an expression it simplifies it and stores the simplified expr in its walkio ptr
 * and the walk data;
 * after the child was visited, it is replaced with the simplified expr
 */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(simplifyExpr)
{
   assert(expr != NULL);
   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD:
      {
         SCIP_CONSEXPR_EXPR* newchild;
         int currentchild;

         currentchild = SCIPgetConsExprExprWalkCurrentChild(expr);
         newchild = (SCIP_CONSEXPR_EXPR*)expr->children[currentchild]->walkio.ptrval;

         SCIP_CALL( SCIPreplaceConsExprExprChild(scip, expr, currentchild, newchild) );

         /* SCIPreplaceConsExprExprChild has captured the new child and we don't need it anymore */
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &newchild) );
         expr->children[currentchild]->walkio.ptrval = NULL;

         /* continue */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }
      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR:
      {
         SCIP_CONSEXPR_EXPR* simplifiedexpr;

         if( *expr->exprhdlr->simplify != NULL )
         {
            SCIP_CALL( (*expr->exprhdlr->simplify)(scip, expr, &simplifiedexpr) );
         }
         else
         {
            assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "sum")  != 0);
            assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "prod") != 0);
            assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "var") != 0);
            assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "abs") != 0);
            assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "log") != 0);
            assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "exp") != 0);

            /* if an expression handler doesn't implement simplify, we assume all those type of expressions are simplified
             * we have to capture it, since it must simulate a "normal" simplified call in which a new expression is created
             */
            simplifiedexpr = expr;
            SCIPcaptureConsExprExpr(simplifiedexpr);
         }
         assert(simplifiedexpr != NULL);
         expr->walkio.ptrval = (void *)simplifiedexpr;

         *(SCIP_CONSEXPR_EXPR**)data = simplifiedexpr;

         /* continue */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      default:
      {
         SCIPABORT(); /* we should never be called in this stage */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE; /*lint !e527*/
         return SCIP_OKAY;
      }
   }
}

/** Implements OR4: default comparison method of expressions of the same type:
 * expr1 < expr2 if and only if expr1_i = expr2_i for all i < k and expr1_k < expr2_k.
 * if there is no such k, use number of children to decide
 * if number of children is equal, both expressions are equal
 * @note: Warning, this method doesn't know about expression data. So if your expressions have special data,
 * you must implement the compare callback: SCIP_DECL_CONSEXPR_EXPRCMP
 */
static
int compareConsExprExprsDefault(
   SCIP_CONSEXPR_EXPR*   expr1,              /**< first expression */
   SCIP_CONSEXPR_EXPR*   expr2               /**< second expression */
   )
{
   int i;
   int nchildren1;
   int nchildren2;
   int compareresult;

   nchildren1 = SCIPgetConsExprExprNChildren(expr1);
   nchildren2 = SCIPgetConsExprExprNChildren(expr2);

   for( i = 0; i < nchildren1 && i < nchildren2; ++i )
   {
      compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[i], SCIPgetConsExprExprChildren(expr2)[i]);
      if( compareresult != 0 )
         return compareresult;
   }

   return nchildren1 == nchildren2 ? 0 : nchildren1 < nchildren2 ? -1 : 1;
}

/** compare expressions
 * @return -1, 0 or 1 if expr1 <, =, > expr2, respectively
 * @note: The given expressions are assumed to be simplified.
 */
int SCIPcompareConsExprExprs(
   SCIP_CONSEXPR_EXPR*   expr1,              /**< first expression */
   SCIP_CONSEXPR_EXPR*   expr2               /**< second expression */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr1;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr2;
   int retval;

   exprhdlr1 = SCIPgetConsExprExprHdlr(expr1);
   exprhdlr2 = SCIPgetConsExprExprHdlr(expr2);

   /* expressions are of the same kind/type; use compare callback or default method */
   if( exprhdlr1 == exprhdlr2 )
   {
      if( exprhdlr1->compare != NULL )
         /* enforces OR1-OR3 */
         return exprhdlr1->compare(expr1, expr2);
      else
         /* enforces OR4 */
         return compareConsExprExprsDefault(expr1, expr2);
   }

   /* expressions are of different kind/type */
   /* enforces OR5 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "val") == 0 )
   {
      return -1;
   }
   /* enforces OR10 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "val") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR6 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "sum") == 0 )
   {
      int compareresult;
      int nchildren;

      nchildren = SCIPgetConsExprExprNChildren(expr1);
      compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[nchildren-1], expr2);

      if( compareresult != 0 )
         return compareresult;

      /* "base" of the largest expression of the sum is equal to expr2, coefficient might tell us that expr2 is larger */
      if( SCIPgetConsExprExprSumCoefs(expr1)[nchildren-1] < 1.0 )
         return -1;

      /* largest expression of sum is larger or equal than expr2 => expr1 > expr2 */
      return 1;
   }
   /* enforces OR10 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "sum") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR7 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "prod") == 0 )
   {
      int compareresult;
      int nchildren;

      nchildren = SCIPgetConsExprExprNChildren(expr1);
      compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[nchildren-1], expr2);

      if( compareresult != 0 )
         return compareresult;

      /* base of the largest expression of the product is equal to expr2, exponent might tell us that expr2 is larger */
      if( SCIPgetConsExprExprProductExponents(expr1)[nchildren-1] < 1.0 )
         return -1;

      /* largest expression of product is larger or equal than expr2 => expr1 > expr2 */
      return 1;
   }
   /* enforces OR10 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "prod") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR8 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "var") == 0 )
      return -1;
   /* enforces OR10 */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "var") == 0 )
      return -SCIPcompareConsExprExprs(expr2, expr1);

   /* enforces OR9 */
   retval = strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), SCIPgetConsExprExprHdlrName(exprhdlr2));
   return retval == 0 ? 0 : retval < 0 ? -1 : 1;
}

/**@} */  /* end of simplifying methods */

/* export this function here, so it can be used by unittests but is not really part of the API */
/** propagates bounds for each sub-expression in the constraint by using variable bounds; the resulting bounds for the
 *  root expression will be intersected with the [lhs,rhs] which might lead to an empty interval
 */
static
SCIP_RETCODE forwardPropCons(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONS*              cons,             /**< constraint to propagate */
   SCIP_Bool               intersect,        /**< should the new expr. bounds be intersected with the previous ones? */
   SCIP_Bool*              infeasible,       /**< buffer to store whether an expression's bounds were propagated to an empty interval */
   SCIP_Bool*              redundant,        /**< buffer to store whether the constraint is redundant */
   int*                    ntightenings      /**< buffer to store the number of auxiliary variable tightenings */
   )
{
   SCIP_INTERVAL interval;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(infeasible != NULL);
   assert(redundant != NULL);
   assert(ntightenings != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *infeasible = FALSE;
   *redundant = FALSE;
   *ntightenings = 0;

   /* propagate active and non-deleted constraints only */
   if( SCIPconsIsDeleted(cons) || !SCIPconsIsActive(cons) )
      return SCIP_OKAY;

   /* use 0 tag to recompute intervals
    * we cannot trust variable bounds from SCIP, so relax them a little bit (a.k.a. epsilon)
    */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, consdata->expr, intersect, 0, SCIPepsilon(scip)) );

   /* @todo delete constraint locally if they are redundant w.r.t. bounds used by the LP solver; the LP solution might
    * violate variable bounds by more than SCIPfeastol() because of relative comparisons
    */
#ifdef SCIP_DISABLED_CODE
   /* if the root expression interval could not be tightened by constraint sides, then the constraint is redundant and
    * should be deleted (locally)
    *
    * @todo how to check this even if we have used the constraint sides during propagation, i.e. intersect is TRUE
    */
   if( !intersect && (SCIPisInfinity(scip, -consdata->lhs) || SCIPisLE(scip, consdata->lhs, consdata->expr->interval.inf))
      && (SCIPisInfinity(scip, consdata->rhs) || SCIPisGE(scip, consdata->rhs, consdata->expr->interval.sup)) )
   {
      SCIPdebugMessage("removing redundant constraint %s activity=[%e,%e] sides=[%e,%e]\n", SCIPconsGetName(cons),
         consdata->expr->interval.inf, consdata->expr->interval.sup, consdata->lhs, consdata->rhs);
      SCIP_CALL( SCIPdelConsLocal(scip, cons) );
      *redundant = TRUE;
      return SCIP_OKAY;
   }
#endif

   /* it may happen that we detect infeasibility during forward propagation if we use previously computed intervals */
   if( SCIPintervalIsEmpty(SCIPinfinity(scip), SCIPgetConsExprExprInterval(consdata->expr)) )
   {
      *infeasible = TRUE;
   }
   else
   {
      /* compare root expression interval with constraint sides; store the result in the root expression */
      SCIPintervalSetBounds(&interval, consdata->lhs, consdata->rhs);

      /* consider auxiliary variable stored in the root expression
       * it might happen that some other plug-ins tighten the bounds of these variables
       * we don't trust these bounds, so relax by epsilon
       */
      if( consdata->expr->auxvar != NULL )
      {
         SCIP_INTERVAL auxvarinterval;
         assert(SCIPvarGetLbLocal(consdata->expr->auxvar) <= SCIPvarGetUbLocal(consdata->expr->auxvar));  /* can SCIP ensure this by today? */

         SCIPintervalSetBounds(&auxvarinterval, SCIPvarGetLbLocal(consdata->expr->auxvar) - SCIPepsilon(scip), SCIPvarGetUbLocal(consdata->expr->auxvar) + SCIPepsilon(scip));
         SCIPintervalIntersect(&interval, interval, auxvarinterval);
      }

      SCIP_CALL( SCIPtightenConsExprExprInterval(scip, consdata->expr, interval, infeasible, ntightenings) );
   }

#ifdef SCIP_DEBUG
   if( *infeasible )
   {
      SCIPdebugMessage(" -> found empty bound for an expression during forward propagation of constraint %s\n",
         SCIPconsGetName(cons));
   }
#endif

   return SCIP_OKAY;
}

/* export this function here, so it can be used by unittests but is not really part of the API */
/** propagates bounds for each sub-expression of a given set of constraints by starting from the root expressions; the
 *  expression will be traversed in breadth first search by using a queue
 *
 *  @note calling this function requires feasible intervals for each sub-expression; this is guaranteed by calling
 *  forwardPropCons() before calling this function
 */
static
SCIP_RETCODE reversePropConss(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONS**             conss,            /**< constraints to propagate */
   int                     nconss,           /**< total number of constraints to propagate */
   SCIP_Bool*              infeasible,       /**< buffer to store whether an expression's bounds were propagated to an empty interval */
   int*                    ntightenings      /**< buffer to store the number of (variable) tightenings */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_QUEUE* queue;
   int i;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(infeasible != NULL);
   assert(ntightenings != NULL);

   *infeasible = FALSE;
   *ntightenings = 0;

   if( nconss == 0 )
      return SCIP_OKAY;

   /* create queue */
   SCIP_CALL( SCIPqueueCreate(&queue, SCIPgetNVars(scip), 2.0) );

   /* add root expressions to the queue */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* propagate active and non-deleted constraints only */
      if( SCIPconsIsDeleted(conss[i]) || !SCIPconsIsActive(conss[i]) )
         continue;

      /* skip expressions that could not have been tightened or do not implement the reverseprop callback; */
      if( !consdata->expr->hastightened || consdata->expr->exprhdlr->reverseprop == NULL )
         continue;

      /* add expressions which are not in the queue so far */
      if( !consdata->expr->inqueue )
      {
         SCIP_CALL( SCIPqueueInsert(queue, (void*) consdata->expr) );
         consdata->expr->inqueue = TRUE;
      }
   }

   /* main loop */
   while( !SCIPqueueIsEmpty(queue) && !(*infeasible) )
   {
      SCIP_CONSEXPR_EXPR* expr;
      int nreds;

      expr = (SCIP_CONSEXPR_EXPR*) SCIPqueueRemove(queue);
      assert(expr != NULL);
      assert(expr->exprhdlr->reverseprop != NULL);

      nreds = 0;

      /* mark that the expression is not in the queue anymore */
      expr->inqueue = FALSE;

#ifdef SCIP_DEBUG
      SCIPdebugMessage("call reverse propagation for ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "\n");
#endif

      /* call reverse propagation callback */
      SCIP_CALL( (*expr->exprhdlr->reverseprop)(scip, expr, infeasible, &nreds) );
      assert(nreds >= 0);
      *ntightenings += nreds;

      /* stop propagation if the problem is infeasible */
      if( *infeasible )
         break;

      /* add tightened children with at least one child to the queue */
      for( i = 0; i < expr->nchildren; ++i )
      {
         SCIP_CONSEXPR_EXPR* child;

         child = expr->children[i];
         assert(child != NULL);

         /* add child to the queue */
         /* @todo put children which are in the queue to the end of it! */
         if( !child->inqueue && child->hastightened && child->nchildren > 0 && child->exprhdlr->reverseprop != NULL )
         {
            SCIP_CALL( SCIPqueueInsert(queue, (void*) child) );
            child->inqueue = TRUE;
         }
      }
   }

   /* free the queue */
   SCIPqueueFree(&queue);

   return SCIP_OKAY;
}

/* export this function here, so it can be used by unittests but is not really part of the API */
/** calls domain propagation for a given set of constraints; the algorithm alternates calls of forward and reverse
 *  propagation; the latter only for nodes which have been tightened during the propagation loop;
 *
 *  the propagation algorithm works as follows:
 *
 *   0.) mark all expressions as non-tightened
 *
 *   1.) apply forward propagation and intersect the root expressions with the constraint sides; mark root nodes which
 *       have been changed after intersecting with the constraint sides
 *
 *   2.) apply reverse propagation to each root expression which has been marked as tightened; don't explore
 *       sub-expressions which have not changed since the beginning of the propagation loop
 *
 *   3.) if we have found enough tightenings go to 1.) otherwise leave propagation loop
 *
 *  @note after calling forward propagation for a constraint we mark this constraint as propagated; this flag might be
 *  reset during the reverse propagation when we find a bound tightening of a variable expression contained in the
 *  constraint; resetting this flag is done in the EVENTEXEC callback of the event handler
 *
 *  @note when using forward and reverse propagation alternatingly we reuse expression intervals computed in previous
 *  iterations
 */
static
SCIP_RETCODE propConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  nchgbds,            /**< buffer to add the number of changed bounds */
   int*                  ndelconss           /**< buffer to add the number of deleted constraints */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
   SCIP_Bool redundant;
   SCIP_Bool success;
   int ntightenings;
   int roundnr;
   int i;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);
   assert(result != NULL);
   assert(nchgbds != NULL);
   assert(*nchgbds >= 0);
   assert(ndelconss != NULL);

   /* no constraints to propagate */
   if( nconss == 0 )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTFIND;
   roundnr = 0;
   cutoff = FALSE;

   /* main propagation loop */
   do
   {
      SCIPdebugMessage("start propagation round %d\n", roundnr);

      /* apply forward propagation; recompute expression intervals if it is called for the first time (this also marks
       * all expressions as non-tightened)
       */
      for( i = 0; i < nconss; ++i )
      {
         consdata = SCIPconsGetData(conss[i]);
         assert(consdata != NULL);

         if( SCIPconsIsActive(conss[i]) && !consdata->ispropagated )
         {
            SCIPdebugMessage("call forwardPropCons() for constraint <%s>\n", SCIPconsGetName(conss[i]));
            SCIPdebugPrintCons(scip, conss[i], NULL);

            cutoff = FALSE;
            redundant = FALSE;
            ntightenings = 0;

            SCIP_CALL( forwardPropCons(scip, conss[i], (roundnr != 0), &cutoff, &redundant, &ntightenings) );
            assert(ntightenings >= 0);
            *nchgbds += ntightenings;

            if( cutoff )
            {
               SCIPdebugMessage(" -> cutoff\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            if( ntightenings > 0 )
               *result = SCIP_REDUCEDDOM;
            if( redundant )
               *ndelconss += 1;

            /* mark constraint as propagated; this will be reset via the event system when we find a variable tightening */
            consdata->ispropagated = TRUE;
         }
      }

      /* apply backward propagation; mark constraint as propagated */
      SCIP_CALL( reversePropConss(scip, conss, nconss, &cutoff, &ntightenings) );

      /* @todo add parameter for the minimum number of tightenings to trigger a new propagation round */
      success = ntightenings > 0;

      if( nchgbds != NULL )
         *nchgbds += ntightenings;

      if( cutoff )
      {
         SCIPdebugMessage(" -> cutoff\n");
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }

      if( success )
         *result = SCIP_REDUCEDDOM;
   }
   while( success && ++roundnr < conshdlrdata->maxproprounds );

   return SCIP_OKAY;
}

/* export this function here, so it can be used by unittests but is not really part of the API */
/** returns all variable expressions contained in a given expression; the array to store all variable expressions needs
 * to be at least of size the number of variables in the expression which is bounded by SCIPgetNVars() since there are
 * no two different variable expression sharing the same variable
 *
 * @note function captures variable expressions
 */
SCIP_RETCODE getVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPR**    varexprs,         /**< array to store all variable expressions */
   int*                    nvarexprs         /**< buffer to store the total number of variable expressions */
   );
SCIP_RETCODE getVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPR**    varexprs,         /**< array to store all variable expressions */
   int*                    nvarexprs         /**< buffer to store the total number of variable expressions */
   )
{
   GETVARS_DATA getvarsdata;

   assert(expr != NULL);
   assert(varexprs != NULL);
   assert(nvarexprs != NULL);

   getvarsdata.nvarexprs = 0;
   getvarsdata.varexprs = varexprs;

   /* use a hash map to dicide whether we have stored a variable expression already */
   SCIP_CALL( SCIPhashmapCreate(&getvarsdata.varexprsmap, SCIPblkmem(scip),
         SCIPcalcHashtableSize(SCIPgetNTotalVars(scip))) );

   /* collect all variable expressions */
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, NULL, NULL, getVarExprsLeaveExpr, (void*)&getvarsdata) );
   *nvarexprs = getvarsdata.nvarexprs;

   /* @todo sort variable expressions here? */

   SCIPhashmapFree(&getvarsdata.varexprsmap);

   return SCIP_OKAY;
}

/** stores all variable expressions into a given constraint */
static
SCIP_RETCODE storeVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSDATA*          consdata          /**< constraint data */
   )
{
   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already*/
   if( consdata->varexprs != NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs == NULL);
   assert(consdata->nvarexprs == 0);

   /* create array to store all variable expressions; the number of variable expressions is bounded by SCIPgetNVars() */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->varexprs, SCIPgetNTotalVars(scip)) );

   SCIP_CALL( getVarExprs(scip, consdata->expr, consdata->varexprs, &(consdata->nvarexprs)) );
   assert(SCIPgetNTotalVars(scip) >= consdata->nvarexprs);

   /* realloc array if there are less variable expression than variables */
   if( SCIPgetNTotalVars(scip) > consdata->nvarexprs )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->varexprs, SCIPgetNTotalVars(scip), consdata->nvarexprs) );
   }

   return SCIP_OKAY;
}

/** frees all variable expression stored in storeVarExprs() */
static
SCIP_RETCODE freeVarExprs(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSDATA*          consdata          /**< constraint data */
   )
{
   int i;

   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already*/
   if( consdata->varexprs == NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   /* release variable expressions */
   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->varexprs[i]) );
      assert(consdata->varexprs[i] == NULL);
   }

   /* free variable expressions */
   SCIPfreeBlockMemoryArrayNull(scip, &consdata->varexprs, consdata->nvarexprs);
   consdata->varexprs = NULL;
   consdata->nvarexprs = 0;

   return SCIP_OKAY;
}

/** computes violation of a constraint */
static
SCIP_RETCODE computeViolation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_SOL*             sol,                /**< solution or NULL if LP solution should be used */
   unsigned int          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real activity;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   SCIP_CALL( SCIPevalConsExprExpr(scip, consdata->expr, sol, soltag) );
   activity = SCIPgetConsExprExprValue(consdata->expr);

   /* consider constraint as violated if it is undefined in the current point */
   if( activity == SCIP_INVALID ) /*lint !e777*/
   {
      consdata->lhsviol = SCIPinfinity(scip);
      consdata->rhsviol = SCIPinfinity(scip);
      return SCIP_OKAY;
   }

   /* compute violations */
   consdata->lhsviol = SCIPisInfinity(scip, -consdata->lhs) ? -SCIPinfinity(scip) : consdata->lhs  - activity;
   consdata->rhsviol = SCIPisInfinity(scip,  consdata->rhs) ? -SCIPinfinity(scip) : activity - consdata->rhs;

   return SCIP_OKAY;
}

/** catch variable events */
static
SCIP_RETCODE catchVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to catch bound change events */
   )
{
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int i;

   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

   /* check if we have catched variable events already */
   if( consdata->vareventdata != NULL )
      return SCIP_OKAY;

   assert(consdata->vareventdata == NULL);

   SCIPdebugMessage("catchVarEvents for %s\n", SCIPconsGetName(cons));

   eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED;

   /* allocate enough memory to store all event data structs */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->vareventdata, consdata->nvarexprs) );

   for( i = 0; i < consdata->nvarexprs; ++i )
   {
      assert(consdata->varexprs[i] != NULL);
      assert(strcmp(consdata->varexprs[i]->exprhdlr->name, "var") == 0);

      var = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
      assert(var != NULL);

      SCIP_CALL( SCIPallocBlockMemory(scip, &(consdata->vareventdata[i])) ); /*lint !e866*/
      consdata->vareventdata[i]->cons = cons;
      consdata->vareventdata[i]->varexpr = consdata->varexprs[i];

      SCIP_CALL( SCIPcatchVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*) consdata->vareventdata[i],
            &(consdata->vareventdata[i]->filterpos)) );
   }

   return SCIP_OKAY;
}

/** drop variable events */
static
SCIP_RETCODE dropVarEvents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EVENTHDLR*       eventhdlr,          /**< event handler */
   SCIP_CONS*            cons                /**< constraint for which to drop bound change events */
   )
{
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int i;

   assert(eventhdlr != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* check if we have catched variable events already */
   if( consdata->vareventdata == NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);
   assert(consdata->vareventdata != NULL);

   eventtype = SCIP_EVENTTYPE_BOUNDCHANGED | SCIP_EVENTTYPE_VARFIXED;

   SCIPdebugMessage("dropVarEvents for %s\n", SCIPconsGetName(cons));

   for( i = consdata->nvarexprs - 1; i >= 0; --i )
   {
      var = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
      assert(var != NULL);

      assert(SCIPgetConsExprExprVarVar(consdata->vareventdata[i]->varexpr) == var);
      assert(consdata->vareventdata[i]->cons == cons);
      assert(consdata->vareventdata[i]->varexpr == consdata->varexprs[i]);
      assert(consdata->vareventdata[i]->filterpos >= 0);

      SCIP_CALL( SCIPdropVarEvent(scip, var, eventtype, eventhdlr, (SCIP_EVENTDATA*) consdata->vareventdata[i], consdata->vareventdata[i]->filterpos) );

      SCIPfreeBlockMemory(scip, &consdata->vareventdata[i]); /*lint !e866*/
      consdata->vareventdata[i] = NULL;
   }

   SCIPfreeBlockMemoryArray(scip, &consdata->vareventdata, consdata->nvarexprs);
   consdata->vareventdata = NULL;

   return SCIP_OKAY;
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{  /*lint --e{715}*/
   SCIP_EVENTTYPE eventtype;
   SCIP_CONSEXPR_EXPR* varexpr;
   SCIP_CONSDATA* consdata;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   assert(eventdata != NULL);

   cons = ((SCIP_VAREVENTDATA*) eventdata)->cons;
   assert(cons != NULL);
   consdata = SCIPconsGetData(cons);
   assert(cons != NULL);

   varexpr = ((SCIP_VAREVENTDATA*) eventdata)->varexpr;
   assert(varexpr != NULL);
   assert(strcmp(varexpr->exprhdlr->name, "var") == 0);

   var = SCIPgetConsExprExprVarVar(varexpr);
   assert(var != NULL);

   eventtype = SCIPeventGetType(event);
   assert((eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != 0 || (eventtype & SCIP_EVENTTYPE_VARFIXED) != 0);

   SCIPdebugMessage("  exec event %u for %s in %s\n", eventtype, SCIPvarGetName(var), SCIPconsGetName(cons));

   /* mark constraint to be propagated again */
   /* TODO: we only need to re-propagate if SCIP_EVENTTYPE_BOUNDTIGHTENED, but we need to reevaluate
    * the intervals (forward-propagation) when SCIP_EVENTTYPE_BOUNDRELAXED
    * at some point we should start using the intevaltag for this
    */
   if( (eventtype & SCIP_EVENTTYPE_BOUNDCHANGED) != (unsigned int) 0 )
   {
      SCIPdebugMessage("  propagate %s again\n", SCIPconsGetName(cons));
      consdata->ispropagated = FALSE;
   }

   return SCIP_OKAY;
}

/** propagates variable locks through expression and adds lock to variables */
static
SCIP_RETCODE propagateLocks(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   int                   nlockspos,          /**< number of positive locks */
   int                   nlocksneg           /**< number of negative locks */
   )
{
   SCIP_CONSHDLR* conshdlr;
   EXPRLOCK_DATA lockdata;

   assert(expr != NULL);

   /* if no locks, then nothing to do, then do nothing */
   if( nlockspos == 0 && nlocksneg == 0 )
      return SCIP_OKAY;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   lockdata.exprvarhdlr = SCIPgetConsExprExprHdlrVar(conshdlr);
   lockdata.nlockspos = nlockspos;
   lockdata.nlocksneg = nlocksneg;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, lockVar, NULL, NULL, NULL, &lockdata) );

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
   SCIP_CONSEXPR_EXPR* expr1;
   SCIP_CONSEXPR_EXPR* expr2;

   expr1 = (SCIP_CONSEXPR_EXPR*)key1;
   expr2 = (SCIP_CONSEXPR_EXPR*)key2;
   assert(expr1 != NULL);
   assert(expr2 != NULL);

   return expr1 == expr2 || SCIPcompareConsExprExprs(expr1, expr2) == 0;
}  /*lint !e715*/

/** get value of hash element when comparing with another expression */
static
SCIP_DECL_HASHKEYVAL(hashCommonSubexprKeyval)
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_HASHMAP* expr2key;

   expr = (SCIP_CONSEXPR_EXPR*) key;
   assert(expr != NULL);

   expr2key = (SCIP_HASHMAP*) userptr;
   assert(expr2key != NULL);
   assert(SCIPhashmapExists(expr2key, (void*)expr));

   return (unsigned int)(size_t)SCIPhashmapGetImage(expr2key, (void*)expr);
}  /*lint !e715*/

/* export this function here, so it can be used by unittests but is not really part of the API */
/** replaces common sub-expressions in the current expression graph by using a hash key for each expression; the
 *  algorithm consists of two steps:
 *
 *  1. traverse through all expressions trees of given constraints and compute for each of them a (not necessarily
 *     unique) hash
 *
 *  2. initialize an empty hash table and traverse through all expression; check for each of them if we can find a
 *     structural equivalent expression in the hash table; if yes we replace the expression by the expression inside the
 *     hash table, otherwise we add it to the hash table
 *
 *  @note the hash keys of the expressions are used for the hashing inside the hash table; to compute if two expressions
 *  (with the same hash) are structurally the same we use the function SCIPcompareConsExprExprs()
 */
SCIP_RETCODE replaceCommonSubexpressions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< total number of constraints */
   );
SCIP_RETCODE replaceCommonSubexpressions(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< total number of constraints */
   )
{
   SCIP_HASHMAP* expr2key;
   SCIP_HASHTABLE* key2expr;
   SCIP_CONSDATA* consdata;
   int i;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);

   /* create empty map to store all sub-expression hashes */
   SCIP_CALL( SCIPhashmapCreate(&expr2key, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip))) );

   /* compute all hashes for each sub-expression */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( consdata->expr != NULL )
      {
         SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, NULL, hashExprVisitingExpr, NULL, hashExprLeaveExpr, (void*)expr2key) );
      }
   }

   /* replace equivalent sub-expressions */
   SCIP_CALL( SCIPhashtableCreate(&key2expr, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPhashmapGetNEntries(expr2key)),
         hashCommonSubexprGetKey, hashCommonSubexprEq, hashCommonSubexprKeyval, (void*)expr2key) );

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CONSEXPR_EXPR* newroot;

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( consdata->expr == NULL )
         continue;

      /* since the root has not been checked for equivalence, it has to be checked separately */
      SCIP_CALL( findEqualExpr(scip, consdata->expr, key2expr, &newroot) );

      if( newroot != NULL )
      {
         assert(newroot != consdata->expr);
         assert(SCIPcompareConsExprExprs(consdata->expr, newroot) == 0);

         SCIPdebugMessage("replacing common root expression of constraint <%s>: %p -> %p\n", SCIPconsGetName(conss[i]), (void*)consdata->expr, (void*)newroot);

         /* remove locks on old expression */
         SCIP_CALL( propagateLocks(scip, consdata->expr, -consdata->nlockspos, -consdata->nlocksneg) );

         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );

         consdata->expr = newroot;
         SCIPcaptureConsExprExpr(newroot);

         /* add locks on new expression */
         SCIP_CALL( propagateLocks(scip, consdata->expr, consdata->nlockspos, consdata->nlocksneg) );
      }
      else
      {
         /* replace equivalent sub-expressions in the tree */
         SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, NULL, commonExprVisitingExpr, NULL, NULL, (void*)key2expr) );
      }
   }

   /* free memory */
   SCIPhashtableFree(&key2expr);
   SCIPhashmapFree(&expr2key);

   return SCIP_OKAY;
}

/** simplifies expressions in constraints */
/* @todo put the constant to the constraint sides
 * @todo call removeFixedAndBoundConstraints() from here and remove it from CONSPRESOL
 */
static
SCIP_RETCODE simplifyConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< total number of constraints */
   )
{
   int i;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(nconss >= 0);

   /* simplify each constraint's expression */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( consdata->expr != NULL )
      {
         SCIP_CONSEXPR_EXPR* simplified;

         SCIP_CALL( SCIPsimplifyConsExprExpr(scip, consdata->expr, &simplified) );

         /* If root expression changed, then we need to take care updating the locks as well (the consdata is the one holding consdata->expr "as a child").
          * If root expression did not change, some subexpression may still have changed, but the locks were taking care of in the corresponding SCIPreplaceConsExprExprChild() call.
          */
         if( simplified != consdata->expr )
         {
            /* remove locks on old expression */
            SCIP_CALL( propagateLocks(scip, consdata->expr, -consdata->nlockspos, -consdata->nlocksneg) );

            /* release old expression */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );

            /* store simplified expression */
            consdata->expr = simplified;

            /* add locks on new expression */
            SCIP_CALL( propagateLocks(scip, consdata->expr, consdata->nlockspos, consdata->nlocksneg) );
         }
         else
         {
            /* The simplify captures simplified in any case, also if nothing has changed.
             * Therefore, we have to release it here.
             */
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplified) );
         }
      }
   }

   return SCIP_OKAY;
}

/** @name Parsing methods
 * @{
 * Here is an attempt at defining the grammar of an expression.
 * We use upper case names for variables (in the grammar sense) and terminals are between "".
 * Loosely speaking, a Base will be any "block", a Factor is a Base to a power, a Term is a product of Factors
 * and an Expression is a sum of terms.
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
 * parse(Expr|Term|Base) returns an SCIP_CONSEXPR_EXPR, while parseFactor returns also the exponent
 *
 * @todo We can change the grammar so that Factor becomes base and we allow a Term to be
 *       <pre> Term       -> Factor { ("*" | "/" | "^") Factor } </pre>
 */

#ifdef PARSE_DEBUG
#define debugParse                      printf
#else
#define debugParse                      while( FALSE ) printf
#endif
static
SCIP_RETCODE parseExpr(SCIP*, SCIP_CONSHDLR*, SCIP_HASHMAP*, const char*, const char**, SCIP_CONSEXPR_EXPR**);

/** Parses base to build a value, variable, sum, or function-like ("func(...)") expression.
 * <pre>
 * Base       -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ")
 * </pre>
 */
static
SCIP_RETCODE parseBase(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between SCIP vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  basetree            /**< buffer to store the expr parsed by Base */
   )
{
   SCIP_VAR* var;

   debugParse("parsing base from %s\n", expr); /*lint !e506 !e681*/

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   if( *expr == '\0' )
   {
      SCIPerrorMessage("Unexpected end of expression string\n");
      return SCIP_READERROR;
   }

   if( *expr == '<' )
   {
      /* parse a variable */
      SCIP_CALL( SCIPparseVarName(scip, expr, &var, (char**)newpos) );

      if( var == NULL )
      {
         SCIPerrorMessage("Could not find variable with name '%s'\n", expr);
         return SCIP_READERROR;
      }
      expr = *newpos;

      /* check if we have already created an expression out of this var */
      if( SCIPhashmapExists(vartoexprvarmap, (void *)var) )
      {
         debugParse("Variable %s has been parsed, capturing its expression\n", SCIPvarGetName(var)); /*lint !e506 !e681*/
         *basetree = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(vartoexprvarmap, (void *)var);
         SCIPcaptureConsExprExpr(*basetree);
      }
      else
      {
         debugParse("First time parsing variable %s, creating varexpr and adding it to hashmap\n", SCIPvarGetName(var)); /*lint !e506 !e681*/
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, conshdlr, basetree, var) );
         SCIP_CALL( SCIPhashmapInsert(vartoexprvarmap, (void*)var, (void*)(*basetree)) );
      }
   }
   else if( *expr == '(' )
   {
      /* parse expression */
      SCIP_CALL( parseExpr(scip, conshdlr, vartoexprvarmap, ++expr, newpos, basetree) );
      expr = *newpos;

      /* expect ')' */
      if( *expr != ')' )
      {
         SCIPerrorMessage("Expected ')', but got <%c> from <%s>\n", *expr, expr);
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, basetree) );
         return SCIP_READERROR;
      }
      ++expr;
      debugParse("Done parsing expression, continue with <%s>\n", expr); /*lint !e506 !e681*/
   }
   else if( isdigit(*expr) )
   {
      /* parse number */
      SCIP_Real value;
      if( !SCIPstrToRealValue(expr, &value, (char**)&expr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", expr);
         return SCIP_READERROR;
      }
      debugParse("Parsed value %g, creating a value-expression.\n", value); /*lint !e506 !e681*/
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, conshdlr, basetree, value) );
   }
   else if( isalpha(*expr) )
   {
      /* a (function) name is coming, should find exprhandler with such name */
      int i;
      char operatorname[SCIP_MAXSTRLEN];
      SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
      SCIP_Bool success;

      /* get name */
      i = 0;
      while( *expr != '(' && !isspace((unsigned char)*expr) && *expr != '\0' )
      {
         operatorname[i] = *expr;
         ++expr;
         ++i;
      }
      operatorname[i] = '\0';

      /* after name we must see a '(' */
      if( *expr != '(' )
      {
         SCIPerrorMessage("Expected '(' after operator name <%s>, but got %s.\n", operatorname, expr);
         return SCIP_READERROR;
      }

      /* search for expression handler */
      exprhdlr = SCIPfindConsExprExprHdlr(conshdlr, operatorname);

      /* check expression handler exists and has a parsing method */
      if( exprhdlr == NULL )
      {
         SCIPerrorMessage("No expression handler with name <%s> found.\n", operatorname);
         return SCIP_READERROR;
      }
      if( exprhdlr->parse == NULL )
      {
         SCIPerrorMessage("Expression handler <%s> has no parsing method.\n", operatorname);
         return SCIP_READERROR;
      }

      /* give control to exprhdlr's parser */
      ++expr;
      SCIP_CALL( exprhdlr->parse(scip, conshdlr, expr, newpos, basetree, &success) );

      if( !success )
      {
         SCIPerrorMessage("Error while expression handler <%s> was parsing %s\n", operatorname, expr);
         assert(*basetree == NULL);
         return SCIP_READERROR;
      }
      expr = *newpos;

      /* we should see the ')' of Op "(" OpExpression ") */
      assert(*expr == ')');

      /* move one character forward */
      ++expr;
   }
   else
   {
      /* Base -> "number" | "<varname>" | "(" Expression ")" | Op "(" OpExpression ") */
      SCIPerrorMessage("Expected a number, (expression), <varname>, Opname(Opexpr), instead got <%c> from %s\n", *expr, expr);
      return SCIP_READERROR;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** Parses a factor and builds a product-expression if there is an exponent, otherwise returns the base expression.
 * <pre>
 * Factor -> Base [ "^" "number" | "^(" "number" ")" ]
 * </pre>
 */
static
SCIP_RETCODE parseFactor(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_Real*            exponent,           /**< buffer to store exponent of Factor */
   SCIP_CONSEXPR_EXPR**  factortree          /**< buffer to store the expr parsed by Factor */
   )
{
   SCIP_CONSEXPR_EXPR*  basetree;

   debugParse("parsing factor from %s\n", expr); /*lint !e506 !e681*/

   if( *expr == '\0' )
   {
      SCIPerrorMessage("Unexpected end of expression string.\n");
      return SCIP_READERROR;
   }

   /* parse Base */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseBase(scip, conshdlr, vartoexprvarmap, expr, newpos, &basetree) );
   *factortree = basetree;
   expr = *newpos;

   /* check if there is an exponent */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '^' )
   {

      ++expr;
      while( isspace((unsigned char)*expr) )
         ++expr;

      if( *expr == '\0' )
      {
         SCIPerrorMessage("Unexpected end of expression string after '^'.\n");
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
         return SCIP_READERROR;
      }

      if( *expr == '(' )
      {
         ++expr;

         /* it is exponent with parenthesis; expect number possibly starting with + or - */
         if( !SCIPstrToRealValue(expr, exponent, (char**)&expr) )
         {
            SCIPerrorMessage("error parsing number from <%s>\n", expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }

         /* expect the ')' */
         while( isspace((unsigned char)*expr) )
            ++expr;
         if( *expr != ')' )
         {
            SCIPerrorMessage("error in parsing exponent: expected ')', received <%c> from <%s>\n", *expr,  expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }
         ++expr;
      }
      else
      {
         /* no parenthesis, we should see just a positive number */

         /* expect a digit */
         if( isdigit(*expr) )
         {
            if( !SCIPstrToRealValue(expr, exponent, (char**)&expr) )
            {
               SCIPerrorMessage("error parsing number from <%s>\n", expr);
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
               return SCIP_READERROR;
            }
         }
         else
         {
            SCIPerrorMessage("error in parsing exponent, expected a digit, received <%c> from <%s>\n", *expr,  expr);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &basetree) );
            return SCIP_READERROR;
         }
      }

      debugParse("parsed the exponent %g\n", *exponent); /*lint !e506 !e681*/
   }
   else
   {
      /* there is no explicit exponent */
      *exponent = 1.0;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** Parses a term and builds a product-expression, where each factor is a child.
 * <pre>
 * Term -> Factor { ("*" | "/" ) Factor }
 * </pre>
 */
static
SCIP_RETCODE parseTerm(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  termtree            /**< buffer to store the expr parsed by Term */
   )
{
   SCIP_CONSEXPR_EXPR* factortree;
   SCIP_Real exponent;

   debugParse("parsing term from %s\n", expr); /*lint !e506 !e681*/

   /* parse Factor */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseFactor(scip, conshdlr, vartoexprvarmap, expr, newpos, &exponent, &factortree) );
   expr = *newpos;

   debugParse("back to parsing Term (we have a Factor with exponent %g), continue parsing from %s\n", exponent, expr); /*lint !e506 !e681*/

   /* check if Terms has another Factor incoming */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '*' || *expr == '/' )
   {
      /* initialize termtree as a product expression with a single term, so we can append the extra Factors */
      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, termtree, 1, &factortree, &exponent, 1.0) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &factortree) );

      /* loop: parse Factor, find next symbol */
      do
      {
         SCIP_RETCODE retcode;
         SCIP_Bool isdivision;

         isdivision = (*expr == '/') ? TRUE : FALSE;

         debugParse("while parsing term, read char %c\n", *expr); /*lint !e506 !e681*/

         ++expr;
         retcode = parseFactor(scip, conshdlr, vartoexprvarmap, expr, newpos, &exponent, &factortree);

         /* release termtree, if parseFactor fails with a read-error */
         if( retcode == SCIP_READERROR )
         {
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, termtree) );
         }
         SCIP_CALL( retcode );

         /* append newly created factor */
         exponent = isdivision ? -exponent : exponent;
         SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, *termtree, factortree, exponent) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &factortree) );

         /* find next symbol */
         expr = *newpos;
         while( isspace((unsigned char)*expr) )
            ++expr;
      } while( *expr == '*' || *expr == '/' );
   }
   else
   {
      /* Term consists of this unique factor^exponent */
      if( exponent != 1.0 )
      {
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, termtree, 1, &factortree, &exponent, 1.0) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &factortree) );
      }
      else
         *termtree = factortree;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** Parses an expression and builds a sum-expression with children.
 * <pre>
 * Expression -> ["+" | "-"] Term { ("+" | "-" | "number *") ] Term }
 * </pre>
 */
static
SCIP_RETCODE parseExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< expression constraint handler */
   SCIP_HASHMAP*         vartoexprvarmap,    /**< hashmap to map between scip vars and var expressions */
   const char*           expr,               /**< expr that we are parsing */
   const char**          newpos,             /**< buffer to store the position of expr where we finished reading */
   SCIP_CONSEXPR_EXPR**  exprtree            /**< buffer to store the expr parsed by Expr */
   )
{
   SCIP_Real sign;
   SCIP_CONSEXPR_EXPR* termtree;

   debugParse("parsing expression %s\n", expr); /*lint !e506 !e681*/

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   /* if '+' or '-', store it */
   sign = 1.0;
   if( *expr == '+' || *expr == '-' )
   {
      debugParse("while parsing expression, read char %c\n", *expr); /*lint !e506 !e681*/
      sign = *expr == '+' ? 1.0 : -1.0;
      ++expr;
   }

   SCIP_CALL( parseTerm(scip, conshdlr, vartoexprvarmap, expr, newpos, &termtree) );
   expr = *newpos;

   debugParse("back to parsing expression (we have the following term), continue parsing from %s\n", expr); /*lint !e506 !e681*/

   /* check if Expr has another Term incoming */
   while( isspace((unsigned char)*expr) )
      ++expr;
   if( *expr == '+' || *expr == '-' )
   {
      if( SCIPgetConsExprExprHdlr(termtree) == SCIPgetConsExprExprHdlrValue(conshdlr) )
      {
         /* initialize exprtree as a sum expression with a constant only, so we can append the following terms */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 0, NULL, NULL, sign * SCIPgetConsExprExprValueValue(termtree)) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );
      }
      else
      {
         /* initialize exprtree as a sum expression with a single term, so we can append the following terms */
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 1, &termtree, &sign, 0.0) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );
      }

      /* loop: parse Term, find next symbol */
      do
      {
         SCIP_RETCODE retcode;
         SCIP_Real coef;

         /* check if we have a "coef * <term>" */
         if( SCIPstrToRealValue(expr, &coef, (char**)newpos) )
         {
            while( isspace((unsigned char)**newpos) )
               ++(*newpos);

            if( **newpos != '*' )
            {
               /* no '*', so fall back to parsing term after sign */
               coef = (*expr == '+') ? 1.0 : -1.0;
               ++expr;
            }
            else
            {
               /* keep coefficient in coef and continue parsing term after coefficient */
               expr = (*newpos)+1;

               while( isspace((unsigned char)*expr) )
                  ++expr;
            }
         }
         else
         {
            coef = (*expr == '+') ? 1.0 : -1.0;
            ++expr;
         }

         debugParse("while parsing expression, read coefficient %g\n", coef); /*lint !e506 !e681*/

         retcode = parseTerm(scip, conshdlr, vartoexprvarmap, expr, newpos, &termtree);

         /* release exprtree if parseTerm fails with an read-error */
         if( retcode == SCIP_READERROR )
         {
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, exprtree) );
         }
         SCIP_CALL( retcode );

         /* append newly created term */
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *exprtree, termtree, coef) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );

         /* find next symbol */
         expr = *newpos;
         while( isspace((unsigned char)*expr) )
            ++expr;
      } while( *expr == '+' || *expr == '-' );
   }
   else
   {
      /* Expr consists of this unique ['+' | '-'] Term */
      if( sign  < 0.0 )
      {
         assert(sign == -1.0);
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, conshdlr, exprtree, 1, &termtree, &sign, 0.0) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &termtree) );
      }
      else
         *exprtree = termtree;
   }

   *newpos = expr;

   return SCIP_OKAY;
}

/** given a cons_expr expression, creates an equivalent classic (nlpi-) expression */
static
SCIP_RETCODE makeClassicExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   sourceexpr,         /**< expression to convert */
   SCIP_EXPR**           targetexpr,         /**< buffer to store pointer to created expression */
   SCIP_CONSEXPR_EXPR**  varexprs,           /**< variable expressions that might occur in expr, their position in this array determines the varidx */
   int                   nvarexprs           /**< number of variable expressions */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   SCIP_EXPR** children = NULL;
   int nchildren;
   int c;

   assert(scip != NULL);
   assert(sourceexpr != NULL);
   assert(targetexpr != NULL);

   exprhdlr = SCIPgetConsExprExprHdlr(sourceexpr);
   nchildren = SCIPgetConsExprExprNChildren(sourceexpr);

   /* collect children expressions from children, if any */
   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );
      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( makeClassicExpr(scip, SCIPgetConsExprExprChildren(sourceexpr)[c], &children[c], varexprs, nvarexprs) );
         assert(children[c] != NULL);
      }
   }

   /* create target expression */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "var") == 0 )
   {
      int varidx;

      /* find variable expression in varexprs array
       * the position in the array determines the index of the variable in the classic expression
       * TODO if varexprs are sorted, then can do this more efficient
       */
      for( varidx = 0; varidx < nvarexprs; ++varidx )
         if( varexprs[varidx] == sourceexpr )
            break;
      assert(varidx < nvarexprs);

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_VARIDX, varidx) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "val") == 0 )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_CONST, SCIPgetConsExprExprValueValue(sourceexpr)) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "sum") == 0 )
   {
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), targetexpr, nchildren, children, SCIPgetConsExprExprSumCoefs(sourceexpr), SCIPgetConsExprExprSumConstant(sourceexpr)) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "prod") == 0 )
   {
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomial, SCIPgetConsExprExprProductCoef(sourceexpr), nchildren, NULL, SCIPgetConsExprExprProductExponents(sourceexpr)) );
      SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), targetexpr, nchildren, children, 1, &monomial, 0.0, FALSE) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "abs") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_ABS, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "exp") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_EXP, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "log") == 0 )
   {
      assert(nchildren == 1);
      assert(children != NULL && children[0] != NULL);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), targetexpr, SCIP_EXPR_LOG, children[0]) );
   }
   else
   {
      SCIPerrorMessage("unsupported expression handler <%s>, cannot convert to classical expression\n", SCIPgetConsExprExprHdlrName(exprhdlr));
      return SCIP_ERROR;
   }

   SCIPfreeBufferArrayNull(scip, &children);

   return SCIP_OKAY;
}

/** given an expression and an array of occurring variable expressions, construct a classic expression tree */
static
SCIP_RETCODE makeClassicExprTree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to convert */
   SCIP_CONSEXPR_EXPR**  varexprs,           /**< variable expressions that occur in expr */
   int                   nvarexprs,          /**< number of variable expressions */
   SCIP_EXPRTREE**       exprtree            /**< buffer to store classic expression tree, or NULL if failed */
)
{
   SCIP_EXPR* classicexpr;
   SCIP_VAR** vars;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(varexprs != NULL);  /* we could also create this here, if NULL; but for now, assume it is given by called */
   assert(exprtree != NULL);

   /* make classic expression */
   SCIP_CALL( makeClassicExpr(scip, expr, &classicexpr, varexprs, nvarexprs) );

   /* make classic expression tree */
   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), exprtree, classicexpr, nvarexprs, 0, NULL) );

   /* set variables in expression tree */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvarexprs) );
   for( i = 0; i < nvarexprs; ++i )
      vars[i] = SCIPgetConsExprExprVarVar(varexprs[i]);
   SCIP_CALL( SCIPexprtreeSetVars(*exprtree, nvarexprs, vars) );
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}

/** create a nonlinear row representation of an expr constraint and stores them in consdata */
static
SCIP_RETCODE createNlRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< expression constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
   }

   if( consdata->expr == NULL )
   {
      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            0, NULL, NULL, 0, NULL, 0, NULL, NULL, consdata->lhs, consdata->rhs) );
   }
   else
   {
      /* get an exprtree representation of the cons-expr-expression */
      SCIP_EXPRTREE* exprtree;

      SCIP_CALL( makeClassicExprTree(scip, consdata->expr, consdata->varexprs, consdata->nvarexprs, &exprtree) );
      if( exprtree == NULL )
      {
         SCIPerrorMessage("could not create classic expression tree from cons_expr expression\n");
         return SCIP_ERROR;
      }

      SCIP_CALL( SCIPcreateNlRow(scip, &consdata->nlrow, SCIPconsGetName(cons), 0.0,
            0, NULL, NULL, 0, NULL, 0, NULL, exprtree, consdata->lhs, consdata->rhs) );
      SCIP_CALL( SCIPexprtreeFree(&exprtree) );
   }

   return SCIP_OKAY;
}

/** registers branching candidates */
static
SCIP_RETCODE registerBranchingCandidates(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< constraints to check */
   int                   nconss,             /**< number of constraints to check */
   int*                  nnotify             /**< counter for number of notifications performed */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR* var;
   int c;
   int i;

   *nnotify = 0;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* violations have been computed during CONSENFOLP */
      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         assert(consdata->varexprs != NULL);

         /* introduce all variables which do not have been fixed */
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            var = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
            assert(var != NULL);

            if( !SCIPisEQ(scip, SCIPcomputeVarLbLocal(scip, var), SCIPcomputeVarUbLocal(scip, var)) )
            {
               SCIP_CALL( SCIPaddExternBranchCand(scip, var, MAX(consdata->lhsviol, consdata->rhsviol), SCIP_INVALID) );
               ++(*nnotify);
            }
         }
      }
   }

   return SCIP_OKAY;
}

/** expression walk callback to create and add auxiliary variables for the outer approximation */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(createAuxVarsEnterExpr)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;
   CREATE_AUXVARS_DATA* createdata;

   assert(expr != NULL);
   assert(result != NULL);
   assert(data != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   createdata = (CREATE_AUXVARS_DATA *)data;
   conshdlr = createdata->conshdlr;
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(conshdlrdata->auxvarid >= 0);

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   /* create and capture auxiliary variable; because of common sub-expressions it might happen that we already added an
    * auxiliary variable to an expression
    */
   if( expr->auxvar == NULL && expr->exprhdlr != SCIPgetConsExprExprHdlrVar(conshdlr)
      && expr->exprhdlr != SCIPgetConsExprExprHdlrValue(conshdlr) )
   {
      char name[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "auxvar_%d", conshdlrdata->auxvarid);

      /* @todo add an unique variable name */
      SCIP_CALL( SCIPcreateVarBasic(scip, &expr->auxvar, name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0,
            SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, expr->auxvar) );
     ++(conshdlrdata->auxvarid);

      SCIPdebugMessage("add auxiliary variable %s for expression %p\n", SCIPvarGetName(expr->auxvar), (void*)expr);

      /* add variable locks in both directions */
      SCIP_CALL( SCIPaddVarLocks(scip, expr->auxvar, 1, 1) );

#ifdef SCIP_DEBUG_SOLUTION
      if( SCIPdebugIsMainscip(scip) )
      {
         /* store debug solution of auxiliary variable */
         assert(createdata->debugsol != NULL);
         SCIP_CALL( SCIPevalConsExprExpr(scip, expr, createdata->debugsol, 0) );
         SCIP_CALL( SCIPdebugAddSolVal(scip, expr->auxvar, SCIPgetConsExprExprValue(expr)) );
      }
#endif
   }
   else
   {
      /* skip nodes which have been already explored */
      *result = SCIP_CONSEXPREXPRWALK_SKIP;
   }

   return SCIP_OKAY;
}

/** expression walk callback to free auxiliary variables created for the outer approximation */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(freeAuxVarsEnterExpr)
{
   assert(expr != NULL);
   assert(result != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   assert((SCIP_CONSHDLR*)data != NULL);
   assert(strcmp(SCIPconshdlrGetName((SCIP_CONSHDLR*)data), CONSHDLR_NAME) == 0);

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   if( expr->auxvar != NULL )
   {
      assert(expr->exprhdlr != SCIPgetConsExprExprHdlrVar((SCIP_CONSHDLR*)data));
      assert(expr->exprhdlr != SCIPgetConsExprExprHdlrValue((SCIP_CONSHDLR*)data));

      SCIPdebugMessage("remove auxiliary variable %s for expression %p\n", SCIPvarGetName(expr->auxvar), (void*)expr);

      /* remove variable locks */
      SCIP_CALL( SCIPaddVarLocks(scip, expr->auxvar, -1, -1) );

      SCIP_CALL( SCIPreleaseVar(scip, &expr->auxvar) );
      assert(expr->auxvar == NULL);
   }
   else
   {
      /* skip nodes which have been already explored */
      *result = SCIP_CONSEXPREXPRWALK_SKIP;
   }

   return SCIP_OKAY;
}

/** creates and adds auxiliary variables for outer approximation; we do not add variables for value and variable
 *  expressions; common sub-expression will share the same auxiliary variable
 */
static
SCIP_RETCODE createAuxVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss              /**< total number of constraints */
   )
{
   CREATE_AUXVARS_DATA createdata;
   SCIP_CONSDATA* consdata;
   int i;

   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);

   createdata.conshdlr = conshdlr;
#ifdef SCIP_DEBUG_SOLUTION
   if( SCIPdebugIsMainscip(scip) )
   {
      createdata.debugsol = NULL;
      SCIP_CALL( SCIPdebugGetSol(scip, &createdata.debugsol) );
      assert(createdata.debugsol != NULL);
   }
#endif
   for( i = 0; i < nconss; ++i )
   {
      assert(conss != NULL && conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( consdata->expr != NULL && consdata->expr->auxvar == NULL )
      {
         SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, createAuxVarsEnterExpr, NULL, NULL, NULL, &createdata) );

         /* set the bounds of the auxiliary variable of the root node to [lhs,rhs] */
         assert(SCIPisInfinity(scip, -SCIPvarGetLbLocal(consdata->expr->auxvar)));
         assert(SCIPisInfinity(scip, SCIPvarGetUbLocal(consdata->expr->auxvar)));

         if( !SCIPisInfinity(scip, -consdata->lhs) )
         {
            SCIP_CALL( SCIPchgVarLb(scip, consdata->expr->auxvar, consdata->lhs) );
         }

         if( !SCIPisInfinity(scip, consdata->rhs) )
         {
            SCIP_CALL( SCIPchgVarUb(scip, consdata->expr->auxvar, consdata->rhs) );
         }
      }
   }

   return SCIP_OKAY;
}

/** frees auxiliary variables which have been added to compute an outer approximation */
static
SCIP_RETCODE freeAuxVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to check for auxiliary variables */
   int                   nconss              /**< total number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int i;

   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);

   for( i = 0; i < nconss; ++i )
   {
      assert(conss != NULL && conss[i] != NULL);

      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      if( consdata->expr != NULL )
      {
         SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, freeAuxVarsEnterExpr, NULL, NULL, NULL, (void*)conshdlr) );
      }
   }

   return SCIP_OKAY;
}

/** expression walk callback for separation initialization */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(initSepaEnterExpr)
{
   INITSEPA_DATA* initsepadata;

   assert(expr != NULL);
   assert(result != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   initsepadata = (INITSEPA_DATA*)data;
   assert(initsepadata != NULL);
   assert(initsepadata->conshdlr != NULL);
   assert(!initsepadata->infeasible);

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   /* skip expression if it has been considered already */
   if( initsepadata->initsepatag == expr->initsepatag )
   {
      *result = SCIP_CONSEXPREXPRWALK_SKIP;
      return SCIP_OKAY;
   }

   if( *(expr->exprhdlr->initsepa) != NULL )
   {
      /* call the separation initialization callback of the expression handler */
      SCIP_CALL( (*expr->exprhdlr->initsepa)(scip, initsepadata->conshdlr, expr, &initsepadata->infeasible) );
   }

   /* stop if we detected infeasibility */
   *result = initsepadata->infeasible ? SCIP_CONSEXPREXPRWALK_ABORT : SCIP_CONSEXPREXPRWALK_CONTINUE;

   /* store the initsepa tag */
   expr->initsepatag = initsepadata->initsepatag;

   return SCIP_OKAY;
}

/** expression walk callback for separation deinitialization */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(exitSepaEnterExpr)
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(result != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   if( *(expr->exprhdlr->exitsepa) != NULL )
   {
      /* call the separation deinitialization callback of the expression handler */
      SCIP_CALL( (*expr->exprhdlr->exitsepa)(scip, expr) );
   }

   return SCIP_OKAY;
}

/** expression walk callback for separating a given solution */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(separateSolEnterExpr)
{
   SEPA_DATA* sepadata;

   assert(expr != NULL);
   assert(result != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_ENTEREXPR);

   sepadata = (SEPA_DATA*)data;
   assert(sepadata != NULL);
   assert(sepadata->result != SCIP_CUTOFF);

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   /* skip expression if it has been considered already */
   if( sepadata->sepatag != 0 && sepadata->sepatag == expr->sepatag )
   {
      *result = SCIP_CONSEXPREXPRWALK_SKIP;
      return SCIP_OKAY;
   }

   /* it only makes sense to call the separation callback if there is a variable attached to the expression */
   if( expr->exprhdlr->sepa != NULL && expr->auxvar != NULL )
   {
      SCIP_RESULT separesult;
      int ncuts;

      separesult = SCIP_DIDNOTFIND;
      ncuts = 0;

      /* call the separation callback of the expression handler */
      SCIP_CALL( (*expr->exprhdlr->sepa)(scip, sepadata->conshdlr, expr, sepadata->sol, sepadata->minviolation, &separesult, &ncuts) );

      assert(ncuts >= 0);
      sepadata->ncuts += ncuts;

      if( separesult == SCIP_CUTOFF )
      {
         assert(ncuts > 0);
         SCIPdebugMessage("found a cutoff -> stop separation\n");
         sepadata->result = SCIP_CUTOFF;
         *result = SCIP_CONSEXPREXPRWALK_ABORT;
      }
      else if( separesult == SCIP_SEPARATED )
      {
         assert(ncuts > 0);
         SCIPdebugMessage("found %d cuts separating the current solution\n", ncuts);
         sepadata->result = SCIP_SEPARATED;
      }
   }

   /* store the separation tag at the expression */
   expr->sepatag = sepadata->sepatag;

   return SCIP_OKAY;
}

/** calls separation initialization callback for each expression */
static
SCIP_RETCODE initSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            infeasible          /**< pointer to store whether the problem is infeasible or not */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   INITSEPA_DATA initsepadata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);
   assert(infeasible != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *infeasible = FALSE;

   initsepadata.infeasible = FALSE;
   initsepadata.conshdlr = conshdlr;
   initsepadata.initsepatag = (++conshdlrdata->lastinitsepatag);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      /* call LP initialization callback for 'initial' constraints only */
      if( SCIPconsIsInitial(conss[c]) )
      {
         consdata = SCIPconsGetData(conss[c]);
         assert(consdata != NULL);

         /* walk through the expression tree and call separation initialization callbacks */
         SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, initSepaEnterExpr, NULL, NULL, NULL, (void*)&initsepadata) );

         if( initsepadata.infeasible )
         {
            SCIPdebugMessage("detect infeasibility for constraint %s during initsepa()\n", SCIPconsGetName(conss[c]));
            *infeasible = TRUE;
            return SCIP_OKAY;
         }
      }
   }

   return SCIP_OKAY;
}

/** calls separation deinitialization callback for each expression */
static
SCIP_RETCODE exitSepa(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_CONSDATA* consdata;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(nconss >= 0);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL);
      assert(conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* walk through the expression tree and call separation deinitialization callbacks */
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, exitSepaEnterExpr, NULL, NULL, NULL, NULL) );
   }

   return SCIP_OKAY;
}

/** tries to separate solution or LP solution by a linear cut
 *
 *  assumes that constraint violations have been computed
 */
static
SCIP_RETCODE separatePoint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   int                   nusefulconss,       /**< number of constraints that seem to be useful */
   SCIP_SOL*             sol,                /**< solution to separate, or NULL if LP solution should be used */
   SCIP_Real             minviolation,       /**< minimal violation of a cut if it should be added to the LP */
   SCIP_RESULT*          result              /**< result of separation */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SEPA_DATA sepadata;
   int c;

   assert(conss != NULL || nconss == 0);
   assert(nconss >= nusefulconss);
   assert(minviolation >= 0.0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      /* skip constraints that are not enabled; skip non-violated constraints */
      if( !SCIPconsIsEnabled(conss[c]) || SCIPconsIsDeleted(conss[c])
         || SCIPisLE(scip, MAX(consdata->lhsviol, consdata->rhsviol), SCIPfeastol(scip)) )
         continue;
      assert(SCIPconsIsActive(conss[c]));

      /* initialize separation data */
      sepadata.conshdlr = conshdlr;
      sepadata.sol = sol;
      sepadata.minviolation = minviolation;
      sepadata.result = SCIP_DIDNOTFIND;
      sepadata.ncuts = 0;
      sepadata.sepatag = ++(conshdlrdata->lastsepatag);

      #ifdef SEPA_DEBUG
      {
         int i;
         printf("separating point\n");
         for( i = 0; i < consdata->nvarexprs; ++i )
         {
            SCIP_VAR* var;
            var = SCIPgetConsExprExprVarVar(consdata->varexprs[i]);
            printf("%s = %g bounds: %g,%g\n", SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
         }
         printf("in constraint\n");
         SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
         SCIPinfoMessage(scip, NULL, ";\n");
      }
      #endif

      /* walk through the expression tree and call separation callback functions */
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, separateSolEnterExpr, NULL, NULL, NULL, (void*)&sepadata) );

      if( sepadata.result == SCIP_CUTOFF || sepadata.result == SCIP_SEPARATED )
      {
         assert(sepadata.ncuts > 0);
         *result = sepadata.result;

         if( *result == SCIP_CUTOFF )
            return SCIP_OKAY;
      }

      /* enforce only useful constraints; others are only checked and enforced if we are still feasible or have not
       * found a separating cut yet
       */
      if( c >= nusefulconss && *result == SCIP_SEPARATED )
         break;
   }

   return SCIP_OKAY;
}

/** removes locally fixed/bound constraints, i.e. constraints for which the root expression is a value or var */
static
SCIP_RETCODE removeFixedAndBoundConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< nonlinear constraints handler */
   SCIP_CONS**           conss,              /**< constraints */
   int                   nconss,             /**< number of constraints */
   SCIP_Bool*            infeasible,         /**< buffer to store whether the node is infeasible */
   int*                  ndelconss           /**< buffer to add the total number of deleted constraints */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_Real value;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(infeasible != NULL);
   assert(ndelconss != NULL);

   *infeasible = FALSE;

   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( consdata->expr->exprhdlr == SCIPgetConsExprExprHdlrValue(conshdlr) )
      {
         value = SCIPgetConsExprExprValueValue(consdata->expr);
         if( (!SCIPisInfinity(scip, -consdata->lhs) && SCIPisLT(scip, value - consdata->lhs, -SCIPfeastol(scip)))
            || (!SCIPisInfinity(scip, consdata->rhs) && SCIPisGT(scip, value - consdata->rhs, SCIPfeastol(scip))) )
         {
            /* we should not stop here since SCIP is calling initSolve() independent if the problem is infeasible or
             * not; having some value expression as the root node of some other constraints trigger asserts at different
             * places, e.g. the creation of a sub-SCIP in the subnlp heuristic
             */
            *infeasible = TRUE;
         }

         /* delete the redundant constraint locally */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
         ++(*ndelconss);
      }

      if( consdata->expr->exprhdlr == SCIPgetConsExprExprHdlrVar(conshdlr) )
      {
         /* backward propagation should have tighthened the bounds of the variable */
         assert(SCIPvarGetLbLocal(SCIPgetConsExprExprVarVar(consdata->expr)) - consdata->lhs >= -SCIPfeastol(scip));
         assert(SCIPvarGetUbLocal(SCIPgetConsExprExprVarVar(consdata->expr)) - consdata->rhs <=  SCIPfeastol(scip));

         /* delete the redundant constraint locally */
         SCIP_CALL( SCIPdelConsLocal(scip, conss[c]) );
         ++(*ndelconss);
      }
   }

   return SCIP_OKAY;
}

/** @} */

/*
 * Callback methods of constraint handler
 */

/* TODO: Implement all necessary constraint handler methods. The methods with #if 0 ... #else #define ... are optional */


/** upgrades quadratic constraint to expr constraint */
static
SCIP_DECL_QUADCONSUPGD(quadconsUpgdExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLR* consexprhdlr;
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_CONSEXPR_EXPR* varexpr;
   SCIP_CONSEXPR_EXPR** varexprs;
   SCIP_CONSEXPR_EXPR* prodexpr;
   SCIP_CONSEXPR_EXPR* twoexprs[2];
   SCIP_QUADVARTERM* quadvarterm;
   SCIP_BILINTERM* bilinterm;
   SCIP_Real two = 2.0;
   int pos;
   int i;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(nupgdconss != NULL);
   assert(upgdconss  != NULL);

   *nupgdconss = 0;

   SCIPdebugMessage("quadconsUpgdExpr called for constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   /* no interest in linear constraints */
   if( SCIPgetNQuadVarTermsQuadratic(scip, cons) == 0 )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* signal that we need more memory */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   if( SCIPgetNBilinTermsQuadratic(scip, cons) > 0 )
   {
      /* we will need SCIPfindQuadVarTermQuadratic later, so ensure now that quad var terms are sorted */
      SCIP_CALL( SCIPsortQuadVarTermsQuadratic(scip, cons) );
   }

   consexprhdlr = SCIPfindConshdlr(scip, "expr");
   assert(consexprhdlr != NULL);

   SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &expr, 0, NULL, NULL, 0.0) );

   /* append linear terms */
   for( i = 0; i < SCIPgetNLinearVarsQuadratic(scip, cons); ++i )
   {
      SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexpr, SCIPgetLinearVarsQuadratic(scip, cons)[i]) );
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, varexpr, SCIPgetCoefsLinearVarsQuadratic(scip, cons)[i]) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexpr) );
   }

   /* array to store variable expression for each quadratic variable */
   SCIP_CALL( SCIPallocBufferArray(scip, &varexprs, SCIPgetNQuadVarTermsQuadratic(scip, cons)) );

   /* create var exprs for quadratic vars; append linear and square part of quadratic terms */
   for( i = 0; i < SCIPgetNQuadVarTermsQuadratic(scip, cons); ++i )
   {
      quadvarterm = &SCIPgetQuadVarTermsQuadratic(scip, cons)[i];

      SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexprs[i], quadvarterm->var) );

      if( quadvarterm->lincoef != 0.0 )
      {
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, varexprs[i], quadvarterm->lincoef) );
      }

      if( quadvarterm->sqrcoef != 0.0 )
      {
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prodexpr, 1, &varexprs[i], &two, 1.0) );
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, prodexpr, quadvarterm->sqrcoef) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
      }
   }

   /* append bilinear terms */
   for( i = 0; i < SCIPgetNBilinTermsQuadratic(scip, cons); ++i)
   {
      bilinterm = &SCIPgetBilinTermsQuadratic(scip, cons)[i];

      SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, cons, bilinterm->var1, &pos) );
      assert(pos >= 0);
      assert(pos < SCIPgetNQuadVarTermsQuadratic(scip, cons));
      assert(SCIPgetQuadVarTermsQuadratic(scip, cons)[pos].var == bilinterm->var1);
      twoexprs[0] = varexprs[pos];

      SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, cons, bilinterm->var2, &pos) );
      assert(pos >= 0);
      assert(pos < SCIPgetNQuadVarTermsQuadratic(scip, cons));
      assert(SCIPgetQuadVarTermsQuadratic(scip, cons)[pos].var == bilinterm->var2);
      twoexprs[1] = varexprs[pos];

      SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prodexpr, 2, twoexprs, NULL, 1.0) );
      SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, prodexpr, bilinterm->coef) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prodexpr) );
   }

   /* release variable expressions */
   for( i = 0; i < SCIPgetNQuadVarTermsQuadratic(scip, cons); ++i )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexprs[i]) );
   }

   SCIPfreeBufferArray(scip, &varexprs);

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsExpr(scip, upgdconss, SCIPconsGetName(cons),
      expr, SCIPgetLhsQuadratic(scip, cons), SCIPgetRhsQuadratic(scip, cons),
      SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
      SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons),  SCIPconsIsLocal(cons),
      SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
      SCIPconsIsStickingAtNode(cons)) );

   SCIPdebugMessage("created expr constraint:\n");
   SCIPdebugPrintCons(scip, *upgdconss, NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** upgrades nonlinear constraint to expr constraint */
static
SCIP_DECL_NONLINCONSUPGD(nonlinconsUpgdExpr)
{
   SCIP_CONSHDLR* consexprhdlr;
   SCIP_EXPRGRAPH* exprgraph;
   SCIP_EXPRGRAPHNODE* node;
   SCIP_CONSEXPR_EXPR* expr;

   assert(nupgdconss != NULL);
   assert(upgdconss != NULL);

   *nupgdconss = 0;

   exprgraph = SCIPgetExprgraphNonlinear(scip, SCIPconsGetHdlr(cons));
   node = SCIPgetExprgraphNodeNonlinear(scip, cons);

   SCIPdebugMessage("nonlinconsUpgdExpr called for constraint <%s>\n", SCIPconsGetName(cons));
   SCIPdebugPrintCons(scip, cons, NULL);

   /* no interest in linear constraints */
   if( node == NULL )
      return SCIP_OKAY;

   consexprhdlr = SCIPfindConshdlr(scip, "expr");
   assert(consexprhdlr != NULL);

   /* try to create a cons_expr expression from an expression graph node */
   SCIP_CALL( SCIPcreateConsExprExpr3(scip, consexprhdlr, &expr, exprgraph, node) );

   /* if that didn't work, then because we do not support a certain expression type yet -> no upgrade */
   if( expr == NULL )
      return SCIP_OKAY;

   if( upgdconsssize < 1 )
   {
      /* request larger upgdconss array */
      *nupgdconss = -1;
      return SCIP_OKAY;
   }

   if( SCIPgetNLinearVarsNonlinear(scip, cons) > 0 )
   {
      /* add linear terms */
      SCIP_CONSEXPR_EXPR* varexpr;
      int i;

      /* ensure expr is a sum expression */
      if( SCIPgetConsExprExprHdlr(expr) != SCIPgetConsExprExprHdlrSum(consexprhdlr) )
      {
         SCIP_CONSEXPR_EXPR* sumexpr;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, &sumexpr, 1, &expr, NULL, 0.0) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

         expr = sumexpr;
      }

      for( i = 0; i < SCIPgetNLinearVarsNonlinear(scip, cons); ++i )
      {
         SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, &varexpr, SCIPgetLinearVarsNonlinear(scip, cons)[i]) );
         SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, expr, varexpr, SCIPgetLinearCoefsNonlinear(scip, cons)[i]) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &varexpr) );
      }
   }

   *nupgdconss = 1;
   SCIP_CALL( SCIPcreateConsExpr(scip, upgdconss, SCIPconsGetName(cons),
      expr, SCIPgetLhsNonlinear(scip, cons), SCIPgetRhsNonlinear(scip, cons),
      SCIPconsIsInitial(cons), SCIPconsIsSeparated(cons), SCIPconsIsEnforced(cons),
      SCIPconsIsChecked(cons), SCIPconsIsPropagated(cons), SCIPconsIsLocal(cons),
      SCIPconsIsModifiable(cons), SCIPconsIsDynamic(cons), SCIPconsIsRemovable(cons),
      SCIPconsIsStickingAtNode(cons)) );

   SCIPdebugMessage("created expr constraint:\n");
   SCIPdebugPrintCons(scip, *upgdconss, NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );

   return SCIP_OKAY;
}

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_CONSHDLRCOPY(conshdlrCopyExpr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(valid != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);

   /* create basic data of constraint handler and include it to scip */
   SCIP_CALL( includeConshdlrExprBasic(scip) );

   /* copy expression handlers */
   SCIP_CALL( copyConshdlrExprExprHdlr(scip, conshdlr, valid) );

   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nexprhdlrs; ++i )
   {
      exprhdlr = conshdlrdata->exprhdlrs[i];
      assert(exprhdlr != NULL);

      if( exprhdlr->freehdlr != NULL )
      {
         SCIP_CALL( (*exprhdlr->freehdlr)(scip, conshdlr, exprhdlr, &exprhdlr->data) );
      }

      SCIPfreeMemory(scip, &exprhdlr->name);
      SCIPfreeMemoryNull(scip, &exprhdlr->desc);

      SCIPfreeMemory(scip, &exprhdlr);
   }

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->exprhdlrs, conshdlrdata->exprhdlrssize);

   SCIPfreeMemory(scip, &conshdlrdata);
   SCIPconshdlrSetData(conshdlr, NULL);

   return SCIP_OKAY;
}


/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->restart )
      return SCIP_OKAY;

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( storeVarExprs(scip, SCIPconsGetData(conss[i])) );
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[i]) );
   }

   return SCIP_OKAY;
}


/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[i]) );
      SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(conss[i])) );
   }

   return SCIP_OKAY;
}


/** presolving initialization method of constraint handler (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_CONSINITPRE(consInitpreExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitpreExpr NULL
#endif


/** presolving deinitialization method of constraint handler (called after presolving has been finished) */
static
SCIP_DECL_CONSEXITPRE(consExitpreExpr)
{  /*lint --e{715}*/

   if( nconss > 0 )
   {
      /* tell SCIP that we have something nonlinear */
      SCIPenableNLP(scip);
   }

   return SCIP_OKAY;
}


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
static
SCIP_DECL_CONSINITSOL(consInitsolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->restart )
      return SCIP_OKAY;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* add nlrow representation to NLP, if NLP had been constructed */
      if( SCIPisNLPConstructed(scip) && SCIPconsIsEnabled(conss[c]) )
      {
         if( consdata->nlrow == NULL )
         {
            SCIP_CALL( createNlRow(scip, conss[c]) );
            assert(consdata->nlrow != NULL);
         }
         SCIP_CALL( SCIPaddNlRow(scip, consdata->nlrow) );
      }
   }

   return SCIP_OKAY;
}

/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
static
SCIP_DECL_CONSEXITSOL(consExitsolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* call separation deinitialization callbacks; after a restart, the rows stored in the expressions are broken */
   SCIP_CALL( exitSepa(scip, conshdlr, conss, nconss) );

   conshdlrdata->restart = restart;
   if( conshdlrdata->restart )
      return SCIP_OKAY;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* free nonlinear row representation */
      if( consdata->nlrow != NULL )
      {
         SCIP_CALL( SCIPreleaseNlRow(scip, &consdata->nlrow) );
      }
   }

   /* remove auxiliary variables from expressions */
   SCIP_CALL( freeAuxVars(scip, conshdlr, conss, nconss) );

   return SCIP_OKAY;
}


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExpr)
{  /*lint --e{715}*/
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->expr != NULL);
   assert((*consdata)->nvarexprs == 0);
   assert((*consdata)->varexprs == NULL);

   /* constraint locks should have been removed */
   assert((*consdata)->nlockspos == 0);
   assert((*consdata)->nlocksneg == 0);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*consdata)->expr) );

   /* free nonlinear row representation */
   if( (*consdata)->nlrow != NULL )
   {
      SCIP_CALL( SCIPreleaseNlRow(scip, &(*consdata)->nlrow) );
   }

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransExpr)
{  /*lint --e{715}*/
   COPY_DATA copydata;
   SCIP_CONSEXPR_EXPR* sourceexpr;
   SCIP_CONSEXPR_EXPR* targetexpr;
   SCIP_CONSDATA* sourcedata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   sourceexpr = sourcedata->expr;

   copydata.targetscip = scip;
   copydata.mapvar = transformVar;
   copydata.mapvardata = NULL;

   /* get a copy of sourceexpr with transformed vars */
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, sourceexpr, copyExpr, NULL, copyExpr, copyExpr, &copydata) );
   targetexpr = copydata.targetexpr;

   if( targetexpr == NULL )
   {
      SCIPerrorMessage("Copying expression in consTransExpr failed.\n");
      return SCIP_ERROR;
   }

   /* create transformed cons (captures targetexpr) */
   SCIP_CALL( SCIPcreateConsExpr(scip, targetcons, SCIPconsGetName(sourcecons),
      targetexpr, sourcedata->lhs, sourcedata->rhs,
      SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
      SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
      SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
      SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   /* release target expr */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &targetexpr) );

   return SCIP_OKAY;
}


/** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved) */
static
SCIP_DECL_CONSINITLP(consInitlpExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* add auxiliary variables to expressions */
   if( !conshdlrdata->restart )
   {
      SCIP_CALL( createAuxVars(scip, conshdlr, conss, nconss) );
   }

   /* call LP initialization callbacks of the expression handlers */
   SCIP_CALL( initSepa(scip, conshdlr, conss, nconss, infeasible) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlr != NULL);

   *result = SCIP_DIDNOTFIND;

   /* compute violations */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], NULL, 0) );
   }

   /* call separation */
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, conshdlrdata->mincutviolationsepa,
         result) );

   return SCIP_OKAY;
}


/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlr != NULL);

   *result = SCIP_DIDNOTFIND;

   /* compute violations */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss[c] != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], NULL, 0) );
   }

   /* call separation */
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, sol, conshdlrdata->mincutviolationsepa,
         result) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Real minviolation;
   SCIP_Real maxviol;
   SCIP_RESULT propresult;
   int nnotify;
   int nchgbds;
   int ndelconss;
   int c;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlr != NULL);

   maxviol = 0.0;

   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], NULL, 0) );
      consdata = SCIPconsGetData(conss[c]);

      /* compute max violation */
      maxviol = MAX3(maxviol, consdata->lhsviol, consdata->rhsviol);
   }
   SCIPdebugMessage("maxviol=%e\n", maxviol);

   *result = SCIPisGT(scip, maxviol, SCIPfeastol(scip)) ? SCIP_INFEASIBLE : SCIP_FEASIBLE;

   if( *result == SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* try to propagate */
   nchgbds = 0;
   ndelconss = 0;
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, &propresult, &nchgbds, &ndelconss) );

   if( propresult == SCIP_CUTOFF || propresult == SCIP_REDUCEDDOM )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* try to separate the LP solution */
   minviolation = MIN(0.75*maxviol, conshdlrdata->mincutviolationenfofac * SCIPfeastol(scip));  /*lint !e666*/
   minviolation = MAX(minviolation, SCIPfeastol(scip));  /*lint !e666*/
   SCIP_CALL( separatePoint(scip, conshdlr, conss, nconss, nusefulconss, NULL, minviolation, result) );

   if( *result == SCIP_CUTOFF || *result == SCIP_SEPARATED )
      return SCIP_OKAY;

   /* find branching candidates */
   SCIP_CALL( registerBranchingCandidates(scip, conss, nconss, &nnotify) );
   SCIPdebugMessage("registered %d external branching candidates\n", nnotify);

   /* all variables have been fixed -> cutoff node */
   /* @todo If we do not branch on linear variables any more this should be changed. We need to introduce linear
    * constraints which are obtained by replacing all fixed non-linear variables as it is done in cons_nonlinear.
    */
   if( nnotify == 0 )
      *result = SCIP_CUTOFF;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_RESULT propresult;
   int nchgbds;
   int ndelconss;
   int nnotify;
   int c;

   *result = SCIP_FEASIBLE;
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( computeViolation(scip, conss[c], NULL, 0) );

      consdata = SCIPconsGetData(conss[c]);
      if( SCIPisGT(scip, MAX(consdata->lhsviol, consdata->rhsviol), SCIPfeastol(scip)) )
      {
         *result = SCIP_INFEASIBLE;
         break;
      }
   }

   if( *result == SCIP_FEASIBLE )
      return SCIP_OKAY;

   /* try to propagate */
   nchgbds = 0;
   ndelconss = 0;
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, &propresult, &nchgbds, &ndelconss) );

   if( (propresult == SCIP_CUTOFF) || (propresult == SCIP_REDUCEDDOM) )
   {
      *result = propresult;
      return SCIP_OKAY;
   }

   /* find branching candidates */
   SCIP_CALL( registerBranchingCandidates(scip, conss, nconss, &nnotify) );
   SCIPdebugMessage("registered %d external branching candidates\n", nnotify);

   *result = nnotify == 0 ? SCIP_SOLVELP : SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   unsigned int soltag;
   int c;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(conss != NULL || nconss == 0);
   assert(result != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_FEASIBLE;
   soltag = ++(conshdlrdata->lastsoltag);

   /* check all nonlinear constraints for feasibility */
   for( c = 0; c < nconss; ++c )
   {
      assert(conss != NULL && conss[c] != NULL);
      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) || SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
      {
         *result = SCIP_INFEASIBLE;

         /* print reason for infeasibility */
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");

            if( SCIPisGT(scip, consdata->lhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhsviol);
            }
            if( SCIPisGT(scip, consdata->rhsviol, SCIPfeastol(scip)) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", consdata->rhsviol);
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropExpr)
{  /*lint --e{715}*/
   int nchgbds;
   int ndelconss;

   nchgbds = 0;
   ndelconss = 0;

   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, result, &nchgbds, &ndelconss) );
   assert(nchgbds >= 0);

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool infeasible;
   int i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   *result = SCIP_DIDNOTRUN;
   if( conshdlrdata->restart )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* simplify constraints */
   SCIP_CALL( simplifyConstraints(scip, conss, nconss) );

   /* replace common subexpressions */
   SCIP_CALL( replaceCommonSubexpressions(scip, conss, nconss) );

   /* FIXME: this is a dirty hack for updating the variable expressions stored inside an expression which might have
    * been changed after simplification; now we completely recollect all variable expression and variable events
    */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, conss[i]) );
      SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(conss[i])) );
   }
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL( storeVarExprs(scip, SCIPconsGetData(conss[i])) );
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, conss[i]) );
   }

   /* propagate constraints */
   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, result, nchgbds, ndelconss) );
   assert(*nchgbds >= 0);
   assert(*ndelconss >= 0);

   /* it might be possible that after simplification only a value expression remains in the root node
    * FIXME: we can't terminate presolve before calling this function, since constraints' expression
    * can't be value nor variable expression in several functions! */
   SCIP_CALL( removeFixedAndBoundConstraints(scip, conshdlr, conss, nconss, &infeasible, ndelconss) );
   assert(*ndelconss >= 0);

   if( infeasible )
   {
      *result = SCIP_CUTOFF;
   }

   if( *result == SCIP_CUTOFF )
      return SCIP_OKAY;

   if( *ndelconss > 0 || *nchgbds > 0 )
      *result = SCIP_SUCCESS;
   else
      *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSRESPROP(consRespropExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consRespropExpr NULL
#endif


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   EXPRLOCK_DATA lockdata;

   assert(conshdlr != NULL);
   assert(cons != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   if( consdata->expr == NULL )
      return SCIP_OKAY;

   lockdata.exprvarhdlr = conshdlrdata->exprvarhdlr;
   lockdata.nlockspos = nlockspos + nlocksneg;
   lockdata.nlocksneg = nlockspos + nlocksneg;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, lockVar, NULL, NULL, NULL, &lockdata) );

   /* remember how the constraint was locked */
   consdata->nlockspos += nlockspos + nlocksneg;
   consdata->nlocksneg += nlockspos + nlocksneg;

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveExpr)
{  /*lint --e{715}*/

   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( storeVarExprs(scip, SCIPconsGetData(cons)) );
   }

   return SCIP_OKAY;
}


/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) > SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
      SCIP_CALL( freeVarExprs(scip, SCIPconsGetData(cons)) );
   }

   return SCIP_OKAY;
}

/** constraint enabling notification method of constraint handler */
static
SCIP_DECL_CONSENABLE(consEnableExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( catchVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** constraint disabling notification method of constraint handler */
static
SCIP_DECL_CONSDISABLE(consDisableExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( dropVarEvents(scip, conshdlrdata->eventhdlr, cons) );
   }

   return SCIP_OKAY;
}

/** variable deletion of constraint handler */
#if 0
static
SCIP_DECL_CONSDELVARS(consDelvarsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDelvarsExpr NULL
#endif


/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintExpr)
{  /*lint --e{715}*/

   SCIP_CONSDATA* consdata;

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   /* print left hand side for ranged constraints */
   if( !SCIPisInfinity(scip, -consdata->lhs)
      && !SCIPisInfinity(scip, consdata->rhs)
      && !SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, "%.15g <= ", consdata->lhs);

   /* print expression */
   if( consdata->expr != NULL )
   {
      SCIP_CALL( SCIPprintConsExprExpr(scip, consdata->expr, file) );
   }
   else
   {
      SCIPinfoMessage(scip, file, "0");
   }

   /* print right hand side */
   if( SCIPisEQ(scip, consdata->lhs, consdata->rhs) )
      SCIPinfoMessage(scip, file, " == %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, consdata->rhs) )
      SCIPinfoMessage(scip, file, " <= %.15g", consdata->rhs);
   else if( !SCIPisInfinity(scip, -consdata->lhs) )
      SCIPinfoMessage(scip, file, " >= %.15g", consdata->lhs);
   else
      SCIPinfoMessage(scip, file, " [free]");

   return SCIP_OKAY;
}


/** constraint copying method of constraint handler */
static
SCIP_DECL_CONSCOPY(consCopyExpr)
{  /*lint --e{715}*/
   COPY_DATA copydata;
   COPY_MAPVAR_DATA mapvardata;
   SCIP_CONSEXPR_EXPR* sourceexpr;
   SCIP_CONSEXPR_EXPR* targetexpr;
   SCIP_CONSDATA* sourcedata;

   assert(cons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   sourceexpr = sourcedata->expr;

   mapvardata.varmap = varmap;
   mapvardata.consmap = consmap;
   mapvardata.global = global;
   mapvardata.valid = TRUE; /* hope the best */

   copydata.targetscip = scip;
   copydata.mapvar = copyVar;
   copydata.mapvardata = &mapvardata;

   /* get a copy of sourceexpr with transformed vars */
   SCIP_CALL( SCIPwalkConsExprExprDF(sourcescip, sourceexpr, copyExpr, NULL, copyExpr, copyExpr, &copydata) );
   targetexpr = copydata.targetexpr;

   if( targetexpr == NULL )
   {
      *cons = NULL;
      *valid = FALSE;

      return SCIP_OKAY;
   }

   /* validity depends only on the SCIPgetVarCopy() returns from copyVar, which are accumulated in mapvardata.valid */
   *valid = mapvardata.valid;

   /* create copy (captures targetexpr) */
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name != NULL ? name : SCIPconsGetName(sourcecons),
      targetexpr, sourcedata->lhs, sourcedata->rhs,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );

   /* release target expr */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &targetexpr) );

   return SCIP_OKAY;
}


/** constraint parsing method of constraint handler */
static
SCIP_DECL_CONSPARSE(consParseExpr)
{  /*lint --e{715}*/
   SCIP_Real  lhs;
   SCIP_Real  rhs;
   const char* endptr;
   SCIP_CONSEXPR_EXPR* consexprtree;

   SCIPdebugMessage("cons_nonlinear::consparse parsing %s\n",str);

   assert(scip != NULL);
   assert(success != NULL);
   assert(str != NULL);
   assert(name != NULL);
   assert(cons != NULL);

   *success = FALSE;

   /* return if string empty */
   if( !*str )
      return SCIP_OKAY;

   endptr = str;

   /* set left and right hand side to their default values */
   lhs = -SCIPinfinity(scip);
   rhs =  SCIPinfinity(scip);

   /* parse constraint to get lhs, rhs, and expression in between (from cons_linear.c::consparse, but parsing whole string first, then getting expression) */

   /* check for left hand side */
   if( isdigit((unsigned char)str[0]) || ((str[0] == '-' || str[0] == '+') && isdigit((unsigned char)str[1])) )
   {
      /* there is a number coming, maybe it is a left-hand-side */
      if( !SCIPstrToRealValue(str, &lhs, (char**)&endptr) )
      {
         SCIPerrorMessage("error parsing number from <%s>\n", str);
         return SCIP_READERROR;
      }

      /* ignore whitespace */
      while( isspace((unsigned char)*endptr) )
         ++endptr;

      if( endptr[0] != '<' || endptr[1] != '=' )
      {
         /* no '<=' coming, so it was the beginning of the expression and not a left-hand-side */
         lhs = -SCIPinfinity(scip);
      }
      else
      {
         /* it was indeed a left-hand-side, so continue parsing after it */
         str = endptr + 2;

         /* ignore whitespace */
         while( isspace((unsigned char)*str) )
            ++str;
      }
   }

   debugParse("str should start at beginning of expr: %s\n", str); /*lint !e506 !e681*/

   /* parse expression: so far we did not allocate memory, so can just return in case of readerror */
   SCIP_CALL( SCIPparseConsExprExpr(scip, conshdlr, str, &str, &consexprtree) );

   /* check for left or right hand side */
   while( isspace((unsigned char)*str) )
      ++str;

   /* check for free constraint */
   if( strncmp(str, "[free]", 6) == 0 )
   {
      if( !SCIPisInfinity(scip, -lhs) )
      {
         SCIPerrorMessage("cannot have left hand side and [free] status \n");
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
         return SCIP_OKAY;
      }
      *success = TRUE;
   }
   else
   {
      switch( *str )
      {
         case '<':
            *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
            break;
         case '=':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have == on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &rhs, (char**)&endptr);
               lhs = rhs;
            }
            break;
         case '>':
            if( !SCIPisInfinity(scip, -lhs) )
            {
               SCIPerrorMessage("cannot have => on rhs if there was a <= on lhs\n");
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
               return SCIP_OKAY;
            }
            else
            {
               *success = SCIPstrToRealValue(str+2, &lhs, (char**)&endptr);
               break;
            }
         case '\0':
            *success = TRUE;
            break;
         default:
            SCIPerrorMessage("unexpected character %c\n", *str);
            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );
            return SCIP_OKAY;
      }
   }

   /* create constraint */
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name,
      consexprtree, lhs, rhs,
      initial, separate, enforce, check, propagate, local, modifiable, dynamic, removable, stickingatnode) );
   assert(*cons != NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consexprtree) );

   debugParse("created expression constraint: <%s>\n", SCIPconsGetName(*cons)); /*lint !e506 !e681*/

   return SCIP_OKAY;
}


/** constraint method of constraint handler which returns the variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETVARS(consGetVarsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetVarsExpr NULL
#endif

/** constraint method of constraint handler which returns the number of variables (if possible) */
#if 0
static
SCIP_DECL_CONSGETNVARS(consGetNVarsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetNVarsExpr NULL
#endif

/** constraint handler method to suggest dive bound changes during the generic diving algorithm */
#if 0
static
SCIP_DECL_CONSGETDIVEBDCHGS(consGetDiveBdChgsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consGetDiveBdChgsExpr NULL
#endif



/** creates the handler for an expression handler and includes it into the expression constraint handler */
SCIP_RETCODE SCIPincludeConsExprExprHdlrBasic(
   SCIP*                       scip,         /**< SCIP data structure */
   SCIP_CONSHDLR*              conshdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR**    exprhdlr,     /**< buffer where to store expression handler */
   const char*                 name,         /**< name of expression handler (must not be NULL) */
   const char*                 desc,         /**< description of expression handler (can be NULL) */
   unsigned int                precedence,   /**< precedence of expression operation (used for printing) */
   SCIP_DECL_CONSEXPR_EXPREVAL((*eval)),     /**< point evaluation callback (can not be NULL) */
   SCIP_CONSEXPR_EXPRHDLRDATA* data          /**< data of expression handler (can be NULL) */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(name != NULL);
   assert(exprhdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocClearMemory(scip, exprhdlr) );

   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*exprhdlr)->name, name, strlen(name)+1) );
   if( desc != NULL )
   {
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*exprhdlr)->desc, desc, strlen(desc)+1) );
   }

   (*exprhdlr)->precedence = precedence;
   (*exprhdlr)->eval = eval;
   (*exprhdlr)->data = data;

   ENSUREBLOCKMEMORYARRAYSIZE(scip, conshdlrdata->exprhdlrs, conshdlrdata->exprhdlrssize, conshdlrdata->nexprhdlrs+1);

   conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs] = *exprhdlr;
   ++conshdlrdata->nexprhdlrs;

   return SCIP_OKAY;
}

/** set the expression handler callbacks to copy and free an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeHdlr(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOPYHDLR((*copyhdlr)), /**< handler copy callback (can be NULL) */
   SCIP_DECL_CONSEXPR_EXPRFREEHDLR((*freehdlr))  /**< handler free callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->copyhdlr = copyhdlr;
   exprhdlr->freehdlr = freehdlr;

   return SCIP_OKAY;
}

/** set the expression handler callbacks to copy and free expression data */
SCIP_RETCODE SCIPsetConsExprExprHdlrCopyFreeData(
   SCIP*                      scip,              /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,          /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,          /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA((*copydata)), /**< expression data copy callback (can be NULL for expressions without data) */
   SCIP_DECL_CONSEXPR_EXPRFREEDATA((*freedata))  /**< expression data free callback (can be NULL if data does not need to be freed) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->copydata = copydata;
   exprhdlr->freedata = freedata;

   return SCIP_OKAY;
}

/** set the simplify callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrSimplify(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRSIMPLIFY((*simplify))  /**< simplify callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->simplify = simplify;

   return SCIP_OKAY;
}

/** set the compare callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrCompare(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRCMP((*compare))    /**< compare callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->compare = compare;

   return SCIP_OKAY;
}

/** set the print callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrPrint(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPRINT((*print))    /**< print callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->print = print;

   return SCIP_OKAY;
}

/** set the parse callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrParse(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPARSE((*parse))    /**< parse callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->parse = parse;

   return SCIP_OKAY;
}

/** set the interval evaluation callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrIntEval(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINTEVAL((*inteval))/**< interval evaluation callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->inteval = inteval;

   return SCIP_OKAY;
}


/** set the separation initialization callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrInitSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRINITSEPA((*initsepa))  /**< separation initialization callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->initsepa = initsepa;

   return SCIP_OKAY;
}

/** set the separation deinitialization callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrExitSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPREXITSEPA((*exitsepa))  /**< separation deinitialization callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->exitsepa = exitsepa;

   return SCIP_OKAY;
}

/** set the separation callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrSepa(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRSEPA((*sepa))      /**< separation callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->sepa = sepa;

   return SCIP_OKAY;
}

/** set the reverse propagation callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrReverseProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_REVERSEPROP((*reverseprop))/**< reverse propagation callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->reverseprop = reverseprop;

   return SCIP_OKAY;
}

/** set the hash callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrHash(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRHASH((*hash))      /**< hash callback (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(exprhdlr != NULL);

   exprhdlr->hash = hash;

   return SCIP_OKAY;
}

/** gives expression handlers */
SCIP_CONSEXPR_EXPRHDLR** SCIPgetConsExprExprHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprhdlrs;
}

/** gives number of expression handlers */
int SCIPgetConsExprExprNHdlrs(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
)
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->nexprhdlrs;
}

/** returns an expression handler of a given name (or NULL if not found) */
SCIP_CONSEXPR_EXPRHDLR* SCIPfindConsExprExprHdlr(
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   const char*                name           /**< name of expression handler */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int h;

   assert(conshdlr != NULL);
   assert(name != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( h = 0; h < conshdlrdata->nexprhdlrs; ++h )
      if( strcmp(SCIPgetConsExprExprHdlrName(conshdlrdata->exprhdlrs[h]), name) == 0 )
         return conshdlrdata->exprhdlrs[h];

   return NULL;
}

/** returns expression handler for variable expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrVar(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprvarhdlr;
}

/** returns expression handler for constant value expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrValue(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprvalhdlr;
}

/** returns expression handler for sum expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrSum(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprsumhdlr;
}

/** returns expression handler for product expressions */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlrProduct(
   SCIP_CONSHDLR*             conshdlr       /**< expression constraint handler */
   )
{
   assert(conshdlr != NULL);

   return SCIPconshdlrGetData(conshdlr)->exprprodhdlr;
}


/** gives the name of an expression handler */
const char* SCIPgetConsExprExprHdlrName(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->name;
}

/** gives the description of an expression handler (can be NULL) */
const char* SCIPgetConsExprExprHdlrDescription(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->desc;
}

/** gives the precedence of an expression handler */
unsigned int SCIPgetConsExprExprHdlrPrecedence(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr       /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->precedence;
}

/** gives the data of an expression handler */
SCIP_CONSEXPR_EXPRHDLRDATA* SCIPgetConsExprExprHdlrData(
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr      /**< expression handler */
)
{
   assert(exprhdlr != NULL);

   return exprhdlr->data;
}

/** creates and captures an expression with given expression data and children */
SCIP_RETCODE SCIPcreateConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data (expression assumes ownership) */
   int                     nchildren,        /**< number of children */
   SCIP_CONSEXPR_EXPR**    children          /**< children (can be NULL if nchildren is 0) */
   )
{
   int c;

   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(children != NULL || nchildren == 0);

   SCIP_CALL( SCIPallocClearBlockMemory(scip, expr) );

   (*expr)->exprhdlr = exprhdlr;
   (*expr)->exprdata = exprdata;

   /* initialize an empty interval for interval evaluation */
   SCIPintervalSetEntire(SCIPinfinity(scip), &(*expr)->interval);

   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*expr)->children, children, nchildren) );
      (*expr)->nchildren = nchildren;
      (*expr)->childrensize = nchildren;

      for( c = 0; c < nchildren; ++c )
         SCIPcaptureConsExprExpr((*expr)->children[c]);
   }

   SCIPcaptureConsExprExpr(*expr);

   return SCIP_OKAY;
}

/** creates and captures an expression with up to two children */
SCIP_RETCODE SCIPcreateConsExprExpr2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr,         /**< expression handler */
   SCIP_CONSEXPR_EXPRDATA* exprdata,         /**< expression data */
   SCIP_CONSEXPR_EXPR*     child1,           /**< first child (can be NULL) */
   SCIP_CONSEXPR_EXPR*     child2            /**< second child (can be NULL) */
   )
{  /*lint --e{715}*/
   assert(expr != NULL);
   assert(exprhdlr != NULL);

   if( child1 != NULL && child2 != NULL )
   {
      SCIP_CONSEXPR_EXPR* pair[2];
      pair[0] = child1;
      pair[1] = child2;

      SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, 2, pair) );
   }
   else if( child2 == NULL )
   {
      SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, child1 == NULL ? 0 : 1, &child1) );
   }
   else
   {
      /* child2 != NULL, child1 == NULL */
      SCIP_CALL( SCIPcreateConsExprExpr(scip, expr, exprhdlr, exprdata, 1, &child2) );
   }

   return SCIP_OKAY;
}

/** creates and captures an expression from a node in an (old-style) expression graph */
SCIP_RETCODE SCIPcreateConsExprExpr3(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSHDLR*          consexprhdlr,     /**< expression constraint handler */
   SCIP_CONSEXPR_EXPR**    expr,             /**< pointer where to store expression */
   SCIP_EXPRGRAPH*         exprgraph,        /**< expression graph */
   SCIP_EXPRGRAPHNODE*     node              /**< expression graph node */
   )
{
   SCIP_CONSEXPR_EXPR** children = NULL;
   int nchildren;
   int c = 0;

   assert(expr != NULL);
   assert(node != NULL);

   *expr = NULL;
   nchildren = SCIPexprgraphGetNodeNChildren(node);

   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

      for( c = 0; c < nchildren; ++c )
      {
         SCIP_CALL( SCIPcreateConsExprExpr3(scip, consexprhdlr, &children[c], exprgraph, SCIPexprgraphGetNodeChildren(node)[c]) );
         if( children[c] == NULL )
            goto TERMINATE;
      }

   }

   switch( SCIPexprgraphGetNodeOperator(node) )
   {
      case SCIP_EXPR_CONST :
         SCIP_CALL( SCIPcreateConsExprExprValue(scip, consexprhdlr, expr, SCIPexprgraphGetNodeOperatorReal(node)) );
         break;

      case SCIP_EXPR_VARIDX :
      {
         int varidx;

         varidx = SCIPexprgraphGetNodeOperatorIndex(node);
         assert(varidx >= 0);
         assert(varidx < SCIPexprgraphGetNVars(exprgraph));

         SCIP_CALL( SCIPcreateConsExprExprVar(scip, consexprhdlr, expr, SCIPexprgraphGetVars(exprgraph)[varidx]) );

         break;
      }

      case SCIP_EXPR_PLUS:
      {
         assert(nchildren == 2);
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 2, children, NULL, 0.0) );

         break;
      }

      case SCIP_EXPR_MINUS:
      {
         SCIP_Real coefs[2] = {1.0, -1.0};

         assert(nchildren == 2);
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 2, children, coefs, 0.0) );

         break;
      }

      case SCIP_EXPR_MUL:
      {
         assert(nchildren == 2);

         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 2, children, NULL, 1.0) );

         break;
      }

      case SCIP_EXPR_DIV:
      {
         SCIP_Real coefs[2] = {1.0, -1.0};

         assert(nchildren == 2);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 2, children, coefs, 1.0) );

         break;
      }

      case SCIP_EXPR_SQUARE:
      {
         SCIP_Real two = 2.0;

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &two, 1.0) );

         break;
      }

      case SCIP_EXPR_SQRT:
      {
         SCIP_Real half = 0.5;

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &half, 1.0) );

         break;
      }

      case SCIP_EXPR_REALPOWER:
      {
         SCIP_Real exponent;

         exponent = SCIPexprgraphGetNodeRealPowerExponent(node);

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &exponent, 1.0) );

         break;
      }

      case SCIP_EXPR_INTPOWER:
      {
         SCIP_Real exponent;

         exponent = (SCIP_Real)SCIPexprgraphGetNodeIntPowerExponent(node);

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &exponent, 1.0) );

         break;
      }

      case SCIP_EXPR_SUM:
      {
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, nchildren, children, NULL, 0.0) );

         break;
      }

      case SCIP_EXPR_PRODUCT:
      {
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, nchildren, children, NULL, 1.0) );

         break;
      }

      case SCIP_EXPR_LINEAR:
      {
         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, nchildren, children, SCIPexprgraphGetNodeLinearCoefs(node), SCIPexprgraphGetNodeLinearConstant(node)) );

         break;
      }

      case SCIP_EXPR_QUADRATIC:
      {
         SCIP_QUADELEM quadelem;
         SCIP_CONSEXPR_EXPR* prod;
         int i;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 0, NULL, NULL, SCIPexprgraphGetNodeQuadraticConstant(node)) );

         /* append linear terms */
         if( SCIPexprgraphGetNodeQuadraticLinearCoefs(node) != NULL )
         {
            for( i = 0; i < nchildren; ++i )
            {
               if( SCIPexprgraphGetNodeQuadraticLinearCoefs(node)[i] != 0.0 )
               {
                  assert(children != NULL);
                  SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, children[i], SCIPexprgraphGetNodeQuadraticLinearCoefs(node)[i]) );
               }
            }
         }

         /* append quadratic terms */
         for( i = 0; i < SCIPexprgraphGetNodeQuadraticNQuadElements(node); ++i )
         {
            quadelem = SCIPexprgraphGetNodeQuadraticQuadElements(node)[i];

            if( quadelem.idx1 == quadelem.idx2 )
            {
               SCIP_Real two = 2.0;
               assert(children != NULL);
               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prod, 1, &children[quadelem.idx1], &two, 1.0) );
            }
            else
            {
               SCIP_CONSEXPR_EXPR* prodchildren[2];

               assert(children != NULL);

               prodchildren[0] = children[quadelem.idx1];
               prodchildren[1] = children[quadelem.idx2];

               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prod, 2, prodchildren, NULL, 1.0) );
            }

            SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, prod, quadelem.coef) );

            SCIP_CALL( SCIPreleaseConsExprExpr(scip, &prod) );
         }

         break;
      }

      case SCIP_EXPR_POLYNOMIAL:
      {
         SCIP_EXPRDATA_MONOMIAL* monom;
         int m;

         SCIP_CALL( SCIPcreateConsExprExprSum(scip, consexprhdlr, expr, 0, NULL, NULL, SCIPexprgraphGetNodePolynomialConstant(node)) );

         /* append monomials */
         for( m = 0; m < SCIPexprgraphGetNodePolynomialNMonomials(node); ++m )
         {
            SCIP_Real* exponents;

            monom = SCIPexprgraphGetNodePolynomialMonomials(node)[m];
            exponents = SCIPexprGetMonomialExponents(monom);

            if( SCIPexprGetMonomialNFactors(monom) == 1 && (exponents == NULL || exponents[0] == 1.0) )
            {
               assert(children != NULL && children[SCIPexprGetMonomialChildIndices(monom)[0]] != NULL);

               /* monom is linear in child -> append child itself */
               SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, children[SCIPexprGetMonomialChildIndices(monom)[0]], SCIPexprGetMonomialCoef(monom)) );
            }
            else
            {
               /* monom is nonlinear -> translate into a product expression */
               SCIP_CONSEXPR_EXPR* monomial;
               int f;

               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &monomial, 0, NULL, NULL, 1.0) );

               for( f = 0; f < SCIPexprGetMonomialNFactors(monom); ++f )
               {
                  assert(children != NULL && children[SCIPexprGetMonomialChildIndices(monom)[f]] != NULL);
                  SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, monomial, children[SCIPexprGetMonomialChildIndices(monom)[f]], exponents ? exponents[f] : 1.0) );
               }

               SCIP_CALL( SCIPappendConsExprExprSumExpr(scip, *expr, monomial, SCIPexprGetMonomialCoef(monom)) );
               SCIP_CALL( SCIPreleaseConsExprExpr(scip, &monomial) );
            }
         }

         break;
      }

      case SCIP_EXPR_EXP:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprExp(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_LOG:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprLog(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_ABS:
      {
         assert(nchildren == 1);
         assert(children != NULL && children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprAbs(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_SIGNPOWER:
      case SCIP_EXPR_SIN:
      case SCIP_EXPR_COS:
      case SCIP_EXPR_TAN:
      case SCIP_EXPR_MIN:
      case SCIP_EXPR_MAX:
      case SCIP_EXPR_SIGN:
      case SCIP_EXPR_USER:
      case SCIP_EXPR_PARAM:
      case SCIP_EXPR_LAST:
      default:
         goto TERMINATE;
   }


TERMINATE:
   /* release all created children expressions (c-1...0) */
   for( --c; c >= 0; --c )
   {
      assert(children != NULL && children[c] != NULL);
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[c]) );
   }

   SCIPfreeBufferArrayNull(scip, &children);

   return SCIP_OKAY;
}

/** gets the number of times the expression is currently captured */
int SCIPgetConsExprExprNUses(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   return expr->nuses;
}

/** captures an expression (increments usage count) */
void SCIPcaptureConsExprExpr(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   ++expr->nuses;
}

/** releases an expression (decrements usage count and possibly frees expression) */
SCIP_RETCODE SCIPreleaseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to expression to be released */
   )
{
   assert(expr != NULL);
   assert(*expr != NULL);

   if( (*expr)->nuses == 1 )
   {
      /* release the auxiliary variable */
      if( (*expr)->auxvar != NULL )
      {
         /* remove variable locks of the auxiliary variable */
         SCIP_CALL( SCIPaddVarLocks(scip, (*expr)->auxvar, -1, -1) );

         SCIP_CALL( SCIPreleaseVar(scip, &(*expr)->auxvar) );
      }
      assert((*expr)->auxvar == NULL);

      /* handle the root expr separately: free its data here */
      if( (*expr)->exprdata != NULL && (*expr)->exprhdlr->freedata != NULL )
      {
         SCIP_CALL( (*expr)->exprhdlr->freedata(scip, *expr) );
      }

      SCIP_CALL( SCIPwalkConsExprExprDF(scip, *expr, NULL, freeExpr, freeExpr, NULL,  NULL) );

      /* handle the root expr separately: free its children and itself here */
      assert((*expr)->nuses == 1);

      /* free children array, if any */
      SCIPfreeBlockMemoryArrayNull(scip, &(*expr)->children, (*expr)->childrensize);

      /* expression should not be locked anymore */
      assert((*expr)->nlockspos == 0);
      assert((*expr)->nlocksneg == 0);

      SCIPfreeBlockMemory(scip, expr);
      assert(*expr == NULL);

      return SCIP_OKAY;
   }

   --(*expr)->nuses;
   assert((*expr)->nuses > 0);
   *expr = NULL;

   return SCIP_OKAY;
}

/** gives the number of children of an expression */
int SCIPgetConsExprExprNChildren(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   return expr->nchildren;
}

/** gives the children of an expression (can be NULL if no children) */
SCIP_CONSEXPR_EXPR** SCIPgetConsExprExprChildren(
   SCIP_CONSEXPR_EXPR*   expr               /**< expression */
   )
{
   assert(expr != NULL);

   return expr->children;
}

/** gets the handler of an expression
 *
 * This identifies the type of the expression (sum, variable, ...).
 */
SCIP_CONSEXPR_EXPRHDLR* SCIPgetConsExprExprHdlr(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->exprhdlr;
}

/** gets the expression data of an expression */
SCIP_CONSEXPR_EXPRDATA* SCIPgetConsExprExprData(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return expr->exprdata;
}

/** returns the variable used for linearizing a given expression (return value might be NULL)
 *
 * @note for variable expression it returns the corresponding variable
 */
SCIP_VAR* SCIPgetConsExprExprLinearizationVar(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression */
   )
{
   assert(expr != NULL);

   return strcmp(expr->exprhdlr->name, "var") == 0 ? SCIPgetConsExprExprVarVar(expr) : expr->auxvar;
}

/** sets the expression data of an expression
 *
 * The pointer to possible old data is overwritten and the
 * freedata-callback is not called before.
 * This function is intended to be used by expression handler.
 */
void SCIPsetConsExprExprData(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_CONSEXPR_EXPRDATA* exprdata          /**< expression data to be set (can be NULL) */
   )
{
   assert(expr != NULL);

   expr->exprdata = exprdata;
}

/** print an expression as info-message */
SCIP_RETCODE SCIPprintConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be printed */
   FILE*                   file              /**< file to print to, or NULL for stdout */
   )
{
   assert(expr != NULL);

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, printExpr, printExpr, printExpr, printExpr, (void*)file) );

   return SCIP_OKAY;
}

/** initializes printing of expressions in dot format */
SCIP_RETCODE SCIPprintConsExprExprDotInit(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata,     /**< buffer to store dot printing data */
   FILE*                   file,             /**< file to print to, or NULL for stdout */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint   /**< info on what to print for each expression */
   )
{
   assert(dotdata != NULL);

   if( file == NULL )
      file = stdout;

   SCIP_CALL( SCIPallocBlockMemory(scip, dotdata) );

   (*dotdata)->file = file;
   (*dotdata)->closefile = FALSE;
   SCIP_CALL( SCIPhashmapCreate(&(*dotdata)->visitedexprs, SCIPblkmem(scip), 1000) );
   (*dotdata)->whattoprint = whattoprint;

   SCIPinfoMessage(scip, file, "strict digraph exprgraph {\n");
   SCIPinfoMessage(scip, file, "node [fontcolor=white, style=filled, rankdir=LR]\n");

   return SCIP_OKAY;
}

/** initializes printing of expressions in dot format to a file with given filename */
SCIP_RETCODE SCIPprintConsExprExprDotInit2(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata,     /**< buffer to store dot printing data */
   const char*             filename,         /**< name of file to print to */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint   /**< info on what to print for each expression */
   )
{
   FILE* f;

   assert(dotdata != NULL);
   assert(filename != NULL);

   f = fopen(filename, "w");
   if( f == NULL )
   {
      SCIPerrorMessage("could not open file <%s> for writing\n", filename);  /* error code would be in errno */
      return SCIP_FILECREATEERROR;
   }

   SCIP_CALL( SCIPprintConsExprExprDotInit(scip, dotdata, f, whattoprint) );
   (*dotdata)->closefile = TRUE;

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPprintConsExprExprDot(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_PRINTDOTDATA* dotdata,      /**< data as initialized by \ref SCIPprintConsExprExprDotInit() */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to be printed */
   )
{
   assert(dotdata != NULL);
   assert(expr != NULL);

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, printExprDot, NULL, NULL, NULL, (void*)dotdata) );

   return SCIP_OKAY;
}

/** finishes printing of expressions in dot format */
SCIP_RETCODE SCIPprintConsExprExprDotFinal(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_PRINTDOTDATA** dotdata      /**< buffer where dot printing data has been stored */
   )
{
   SCIP_HASHMAPLIST* list;
   FILE* file;
   int l;

   assert(dotdata != NULL);
   assert(*dotdata != NULL);

   file = (*dotdata)->file;
   assert(file != NULL);

   /* tell dot that all expressions without children have the same rank */
   SCIPinfoMessage(scip, file, "{rank=same;");
   for( l = 0; l < SCIPhashmapGetNLists((*dotdata)->visitedexprs); ++l )
      for( list = SCIPhashmapGetList((*dotdata)->visitedexprs, l); list != NULL; list = SCIPhashmapListGetNext(list) )
         if( SCIPgetConsExprExprNChildren((SCIP_CONSEXPR_EXPR*)SCIPhashmapListGetOrigin(list)) == 0 )
            SCIPinfoMessage(scip, file, " n%p", SCIPhashmapListGetOrigin(list));
   SCIPinfoMessage(scip, file, "}\n");

   SCIPinfoMessage(scip, file, "}\n");

   SCIPhashmapFree(&(*dotdata)->visitedexprs);

   if( (*dotdata)->closefile )
      fclose((*dotdata)->file);

   SCIPfreeBlockMemory(scip, dotdata);

   return SCIP_OKAY;
}

/** shows a single expression by use of dot and gv
 *
 * This function is meant for debugging purposes.
 * It prints the expression into a temporary file in dot format, then calls dot to create a postscript file, then calls ghostview (gv) to show the file.
 * SCIP will hold until ghostscript is closed.
 */
SCIP_RETCODE SCIPshowConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to be printed */
   )
{
   /* this function is for developers, so don't bother with C variants that don't have popen() */
#if _POSIX_C_SOURCE < 2
   SCIPerrorMessage("No POSIX version 2. Try http://distrowatch.com/.");
   return SCIP_ERROR;
#else
   SCIP_CONSEXPR_PRINTDOTDATA* dotdata;
   FILE* f;

   assert(expr != NULL);

   /* call dot to generate postscript output and show it via ghostview */
   f = popen("dot -Tps | gv -", "w");
   if( f == NULL )
   {
      SCIPerrorMessage("Calling popen() failed");
      return SCIP_FILECREATEERROR;
   }

   /* print all of the expression into the pipe */
   SCIP_CALL( SCIPprintConsExprExprDotInit(scip, &dotdata, f, SCIP_CONSEXPR_PRINTDOT_ALL) );
   SCIP_CALL( SCIPprintConsExprExprDot(scip, dotdata, expr) );
   SCIP_CALL( SCIPprintConsExprExprDotFinal(scip, &dotdata) );

   /* close the pipe */
   (void) pclose(f);

   return SCIP_OKAY;
#endif
}


/** evaluate an expression in a point
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * Value can be received via SCIPgetConsExprExprEvalValue().
 * If an evaluation error (division by zero, ...) occurs, this value will
 * be set to SCIP_INVALID.
 *
 * If a nonzero \p soltag is passed, then only (sub)expressions are
 * reevaluated that have a different solution tag. If a soltag of 0
 * is passed, then subexpressions are always reevaluated.
 * The tag is stored together with the value and can be received via
 * SCIPgetConsExprExprEvalTag().
 */
SCIP_RETCODE SCIPevalConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be evaluated */
   SCIP_SOL*               sol,              /**< solution to be evaluated */
   unsigned int            soltag            /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   EXPREVAL_DATA evaldata;

   /* if value is up-to-date, then nothing to do */
   if( soltag != 0 && expr->evaltag == soltag )
      return SCIP_OKAY;

   evaldata.sol = sol;
   evaldata.soltag = soltag;
   evaldata.aborted = FALSE;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, evalExprVisitChild, NULL, evalExprLeaveExpr, &evaldata) );

   if( evaldata.aborted )
   {
      expr->evalvalue = SCIP_INVALID;
      expr->evaltag = soltag;
   }

   return SCIP_OKAY;
}

/** evaluates an expression over a box
 *
 * Initiates an expression walk to also evaluate children, if necessary.
 * The resulting interval can be received via SCIPgetConsExprExprEvalInterval().
 * If the box does not overlap with the domain of the function behind the expression
 * (e.g., sqrt([-2,-1]) or 1/[0,0]) this interval will be empty.
 *
 * For variables, the local variable bounds, possibly relaxed by the amount given
 * by varboundrelax, are used as interval. In the current implementation, variable
 * bounds are relaxed by varboundrelax if that does not change the sign of the bound,
 * and to 0.0 otherwise.
 *
 * If a nonzero \p boxtag is passed, then only (sub)expressions are
 * reevaluated that have a different tag. If a tag of 0 is passed,
 * then subexpressions are always reevaluated.
 * The tag is stored together with the interval and can be received via
 * SCIPgetConsExprExprEvalIntervalTag().
 */
SCIP_RETCODE SCIPevalConsExprExprInterval(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be evaluated */
   SCIP_Bool               intersect,        /**< should the new expr. bounds be intersected with the previous ones? */
   unsigned int            boxtag,           /**< tag that uniquely identifies the current variable domains (with its values), or 0 */
   SCIP_Real               varboundrelax     /**< amount by which variable bounds should be relaxed (at most) */
   )
{
   EXPRINTEVAL_DATA propdata;

   /* if value is up-to-date, then nothing to do */
   if( boxtag != 0 && expr->intevaltag == boxtag )
      return SCIP_OKAY;

   propdata.aborted = FALSE;
   propdata.boxtag = boxtag;
   propdata.intersect = intersect;
   propdata.varboundrelax = varboundrelax;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, intevalExprVisitChild, NULL, intevalExprLeaveExpr, &propdata) );

   if( propdata.aborted )
   {
      SCIPintervalSetEmpty(&expr->interval);
      expr->intevaltag = boxtag;
   }

   return SCIP_OKAY;
}

/** tightens the bounds of an expression and stores the result in the expression interval; variables in variable
 *  expression will be tightened immediately if SCIP is in a stage above SCIP_STAGE_TRANSFORMED
 */
SCIP_RETCODE SCIPtightenConsExprExprInterval(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be tightened */
   SCIP_INTERVAL           newbounds,        /**< new bounds for the expression */
   SCIP_Bool*              cutoff,           /**< buffer to store whether a node's bounds were propagated to an empty interval */
   int*                    ntightenings      /**< buffer to add the total number of tightenings */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(cutoff != NULL);
   assert(ntightenings != NULL);
   assert(!SCIPintervalIsEmpty(SCIPinfinity(scip), expr->interval));

   *cutoff = FALSE;

   /* check if the new bounds lead to an empty interval */
   if( SCIPintervalGetInf(newbounds) > SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr))
      || SCIPintervalGetSup(newbounds) < SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)) )
   {
      SCIPintervalSetEmpty(&expr->interval);
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* do not consider very small tightenings */
   if( SCIPisLbBetter(scip, SCIPintervalGetInf(newbounds), SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)), SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)))
      || SCIPisUbBetter(scip, SCIPintervalGetSup(newbounds), SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)), SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr))) )
   {
      SCIP_VAR* var;

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "tighten bounds of ");
      SCIP_CALL( SCIPprintConsExprExpr(scip, expr, NULL) );
      SCIPinfoMessage(scip, NULL, "");
      SCIPinfoMessage(scip, NULL, " from [%e,%e] to ", SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)), SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)));
#endif

      /* intersect old bound with the found one */
      SCIPintervalIntersect(&expr->interval, expr->interval, newbounds);

#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip, NULL, "[%e,%e]\n", SCIPintervalGetInf(SCIPgetConsExprExprInterval(expr)), SCIPintervalGetSup(SCIPgetConsExprExprInterval(expr)));
#endif

      /* mark expression as tightened; important for reverse propagation to ignore irrelevant sub-expressions*/
      expr->hastightened = TRUE;

      /* do not tighten variable in problem stage (important for unittests)
       * TODO put some kind of #ifdef UNITTEST around this once the unittest are modified to include the .c file (again)? */
      if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
         return SCIP_OKAY;

      var = SCIPgetConsExprExprLinearizationVar(expr);

      if( var != NULL )
      {
         SCIP_Bool tightened;

#ifdef SCIP_DEBUG
         SCIP_Real oldlb;
         SCIP_Real oldub;
         oldlb = SCIPvarGetLbLocal(SCIPgetConsExprExprVarVar(expr));
         oldub = SCIPvarGetUbLocal(SCIPgetConsExprExprVarVar(expr));
#endif

         /* tighten lower bound */
         SCIP_CALL( SCIPtightenVarLb(scip, var, SCIPintervalGetInf(expr->interval), FALSE, cutoff, &tightened) );

         if( tightened )
            ++(*ntightenings);

         /* tighten upper bound */
         if( !(*cutoff) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var, SCIPintervalGetSup(expr->interval), FALSE, cutoff, &tightened) );

            if( tightened )
               ++(*ntightenings);
         }

#ifdef SCIP_DEBUG
         SCIPdebugMessage("tighten bounds of %s from [%e, %e] -> [%e, %e]\n", SCIPvarGetName(var), oldlb, oldub, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var));
#endif
      }
   }

   return SCIP_OKAY;
}

/** gives the value from the last evaluation of an expression (or SCIP_INVALID if there was an eval error) */
SCIP_Real SCIPgetConsExprExprValue(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evalvalue;
}

/** returns the interval from the last interval evaluation of an expression (interval can be empty) */
SCIP_INTERVAL SCIPgetConsExprExprInterval(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->interval;
}

/** gives the evaluation tag from the last evaluation, or 0 */
unsigned int SCIPgetConsExprExprEvalTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->evaltag;
}

/** gives the box tag from the last interval evaluation, or 0 */
unsigned int SCIPgetConsExprExprEvalIntervalTag(
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   assert(expr != NULL);

   return expr->intevaltag;
}

/** sets the evaluation value */
void SCIPsetConsExprExprEvalValue(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_Real               value,            /**< value to set */
   unsigned int            tag               /**< tag of solution that was evaluated, or 0 */
   )
{
   assert(expr != NULL);

   expr->evalvalue = value;
   expr->evaltag = tag;
}

/** sets the evaluation interval */
void SCIPsetConsExprExprEvalInterval(
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   SCIP_INTERVAL*          interval,         /**< interval to set */
   unsigned int            tag               /**< tag of variable domains that were evaluated, or 0. */
   )
{
   assert(expr != NULL);

   SCIPintervalSetBounds(&expr->interval, SCIPintervalGetInf(*interval), SCIPintervalGetSup(*interval));
   expr->intevaltag = tag;
}

/** returns the hash key of an expression */
SCIP_RETCODE SCIPgetConsExprExprHashkey(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression */
   unsigned int*           hashkey           /**< pointer to store the hash key */
   )
{
   SCIP_HASHMAP* expr2key;

   assert(expr != NULL);

   SCIP_CALL( SCIPhashmapCreate(&expr2key, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip))) );

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, NULL, NULL, hashExprLeaveExpr, (void*)expr2key) );

   assert(SCIPhashmapExists(expr2key, (void*)expr));  /* we just computed the hash, so should be in the map */
   *hashkey = (unsigned int)(size_t)SCIPhashmapGetImage(expr2key, (void*)expr);

   SCIPhashmapFree(&expr2key);

   return SCIP_OKAY;
}



/** walks the expression graph in depth-first manner and executes callbacks at certain places
 *
 * Many algorithms over expression trees need to traverse the tree in depth-first manner and a
 * natural way of implementing this algorithms is using recursion.
 * In general, a function which traverses the tree in depth-first looks like
 * <pre>
 * fun( expr )
 *    enterexpr()
 *    continue skip or abort
 *       for( child in expr->children )
 *          visitingchild()
 *          continue skip or abort
 *          fun(child, data, proceed)
 *          visitedchild()
 *          continue skip or abort
 *    leaveexpr()
 * </pre>
 * Given that some expressions might be quite deep we provide this functionality in an iterative fashion.
 *
 * Consider an expression (x*y) + z + log(x-y).
 * The corresponding expression graph is
 * <pre>
 *           [+]
 *       /    |   \
 *    [*]     |    [log]
 *    / \     |      |
 *   /   \    |     [-]
 *   |   |    |     / \
 *  [x] [y]  [z]  [x] [y]
 * </pre>
 * (where [x] and [y] are actually the same expression).
 *
 * If given a pointer to the [+] expression is given as root to this expression, it will walk
 * the graph in a depth-first manner and call the given callback methods at various stages.
 * - When entering an expression, it calls the enterexpr callback.
 *   The SCIPgetConsExprExprWalkParent() function indicates from where the expression has been entered (NULL for the root expression).
 * - Before visiting a child of an expression, it calls the visitingchild callback.
 *   The SCIPgetConsExprExprWalkCurrentChild() function returns which child will be visited (as an index in the current expr's children array).
 * - When returning from visiting a child of an expression, the visitedchild callback is called.
 *   Again the SCIPgetConsExprExprWalkCurrentChild() function returns which child has been visited.
 * - When leaving an expression, it calls the leaveexpr callback.
 *
 * Thus, for the above expression, the callbacks are called in the following order:
 * - enterexpr([+])
 * - visitingchild([+])  currentchild == 0
 * - enterexpr([*])
 * - visitingchild([*])  currentchild == 0
 * - enterexpr([x])
 * - leaveexpr([x])
 * - visitedchild([*])   currentchild == 0
 * - visitingchild([*])  currentchild == 1
 * - enterexpr([y])
 * - leaveexpr([y])
 * - visitedchild([*])   currentchild == 1
 * - leaveexpr([*])
 * - visitedchild([+])   currentchild == 0
 * - visitingchild([+])  currentchild == 1
 * - enterexpr([z])
 * - leaveexpr([z])
 * - visitedchild([+])   currentchild == 1
 * - visitingchild([+])  currentchild == 2
 * - enterexpr([log])
 * - visitingchild([log]) currentchild == 0
 * - enterexpr([-])
 * - visitingchild([-])  currentchild == 0
 * - enterexpr([x])
 * - leaveexpr([x])
 * - visitedchild([-])   currentchild == 0
 * - visitingchild([-])  currentchild == 1
 * - enterexpr([y])
 * - leaveexpr([y])
 * - visitedchild([-])   currentchild == 1
 * - leaveexpr([-])
 * - visitedchild([log]) currentchild == 0
 * - leaveexpr([log])
 * - visitedchild([+])   currentchild == 2
 * - leaveexpr([+])
 *
 * The callbacks can direct the walking method to skip parts of the tree or abort.
 * If returning SCIP_CONSEXPREXPRWALK_SKIP as result of an enterexpr callback, all children of that expression will be skipped. The leaveexpr callback will still be called.
 * If returning SCIP_CONSEXPREXPRWALK_SKIP as result of an visitingchild callback, visiting the current child will be skipped.
 * If returning SCIP_CONSEXPREXPRWALK_SKIP as result of an visitedchild callback, visiting the remaining children will be skipped.
 * If returning SCIP_CONSEXPREXPRWALK_ABORT in any of the callbacks, the walk will be aborted immediately.
 *
 * @note The walkio member of the root expression is reset to its previous value when the walk finishes.
 */
SCIP_RETCODE SCIPwalkConsExprExprDF(
   SCIP*                 scip,                         /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   root,                         /**< the root expression from where to start the walk */
   SCIP_DECL_CONSEXPREXPRWALK_VISIT((*enterexpr)),     /**< callback to be called when entering an expression, or NULL */
   SCIP_DECL_CONSEXPREXPRWALK_VISIT((*visitingchild)), /**< callback to be called before visiting a child, or NULL */
   SCIP_DECL_CONSEXPREXPRWALK_VISIT((*visitedchild)),  /**< callback to be called when returning from a child, or NULL */
   SCIP_DECL_CONSEXPREXPRWALK_VISIT((*leaveexpr)),     /**< callback to be called when leaving an expression, or NULL */
   void*                 data                          /**< data to be passed on to callback methods, or NULL */
   )
{
   SCIP_CONSEXPREXPRWALK_STAGE  stage;
   SCIP_CONSEXPREXPRWALK_RESULT result;
   SCIP_CONSEXPR_EXPR*          child;
   SCIP_CONSEXPR_EXPR*          oldroot;
   SCIP_CONSEXPR_EXPR*          oldparent;
   SCIP_CONSEXPREXPRWALK_IO     oldwalkio;
   int                          oldcurrentchild;
   SCIP_Bool                    aborted;

   assert(root != NULL);

   /* remember the current root, data, child and parent of incoming root, in case we are called from within another walk
    * furthermore, we need to capture the root, because we don't want nobody somebody to invalidate it while we have it
    * NOTE: no callback should touch walkparent, nor walkcurrentchild: these are internal fields of the walker!
    */
   SCIPcaptureConsExprExpr(root);
   oldroot         = root;
   oldcurrentchild = root->walkcurrentchild;
   oldparent       = root->walkparent;
   oldwalkio       = root->walkio;

   /* traverse the tree */
   root->walkcurrentchild = 0;
   root->walkparent = NULL;
   result = SCIP_CONSEXPREXPRWALK_CONTINUE;
   stage = SCIP_CONSEXPREXPRWALK_ENTEREXPR;
   aborted = FALSE;
   while( !aborted )
   {
      switch( stage )
      {
         case SCIP_CONSEXPREXPRWALK_ENTEREXPR:
            assert(root->walkcurrentchild == 0);
            if( enterexpr != NULL )
            {
               SCIP_CALL( (*enterexpr)(scip, root, stage, data, &result) );
               switch( result )
               {
                  case SCIP_CONSEXPREXPRWALK_CONTINUE :
                     break;
                  case SCIP_CONSEXPREXPRWALK_SKIP :
                     root->walkcurrentchild = root->nchildren;
                     break;
                  case SCIP_CONSEXPREXPRWALK_ABORT :
                     aborted = TRUE;
                     break;
               }
            }
            /* goto start visiting children */
            stage = SCIP_CONSEXPREXPRWALK_VISITINGCHILD;
            break;

         case SCIP_CONSEXPREXPRWALK_VISITINGCHILD:
            /* if there are no more children to visit, goto leave */
            if( root->walkcurrentchild >= root->nchildren )
            {
               assert(root->walkcurrentchild == root->nchildren);
               stage = SCIP_CONSEXPREXPRWALK_LEAVEEXPR;
               break;
            }
            /* prepare visit */
            if( visitingchild != NULL )
            {
               SCIP_CALL( (*visitingchild)(scip, root, stage, data, &result) );
               if( result == SCIP_CONSEXPREXPRWALK_SKIP )
               {
                  /* ok, we don't go down, but skip the child: continue and try again with next child (if any) */
                  ++root->walkcurrentchild;
                  continue;
               }
               else if( result == SCIP_CONSEXPREXPRWALK_ABORT )
               {
                  aborted = TRUE;
                  break;
               }
            }
            /* remember the parent and set the first child that should be visited of the new root */
            child = root->children[root->walkcurrentchild];
            child->walkparent = root;
            child->walkcurrentchild = 0;
            root = child;
            /* visit child */
            stage = SCIP_CONSEXPREXPRWALK_ENTEREXPR;
            break;

         case SCIP_CONSEXPREXPRWALK_VISITEDCHILD:
            if( visitedchild != NULL )
            {
               SCIP_CALL( (*visitedchild)(scip, root, stage, data, &result) );
               switch( result )
               {
                  case SCIP_CONSEXPREXPRWALK_CONTINUE :
                     /* visit next (if any) */
                     ++root->walkcurrentchild;
                     break;
                  case SCIP_CONSEXPREXPRWALK_SKIP :
                     /* skip visiting the rest of the children */
                     root->walkcurrentchild = root->nchildren;
                     break;
                  case SCIP_CONSEXPREXPRWALK_ABORT :
                     aborted = TRUE;
                     break;
               }
            }
            else
            {
               /* visit next child (if any) */
               ++root->walkcurrentchild;
            }
            /* goto visiting (next) */
            stage = SCIP_CONSEXPREXPRWALK_VISITINGCHILD;
            break;

         case SCIP_CONSEXPREXPRWALK_LEAVEEXPR:
            if( leaveexpr != NULL )
            {
               SCIP_CONSEXPR_EXPR* parent;

               /* store parent in case the callback frees root */
               parent = root->walkparent;

               SCIP_CALL( (*leaveexpr)(scip, root, stage, data, &result) );
               switch( result )
               {
                  case SCIP_CONSEXPREXPRWALK_CONTINUE :
                     break;
                  case SCIP_CONSEXPREXPRWALK_SKIP :
                     /* this result is not allowed here */
                     SCIPABORT();
                     break;
                  case SCIP_CONSEXPREXPRWALK_ABORT :
                     aborted = TRUE;
                     break;
               }

               /* go up */
               root = parent;
            }
            else
            {
               /* go up */
               root = root->walkparent;
            }
            /* if we finished with the real root (walkparent == NULL), we are done */
            if( root == NULL )
               aborted = TRUE;

            /* goto visited */
            stage = SCIP_CONSEXPREXPRWALK_VISITEDCHILD;
            break;
         default:
            /* unknown stage */
            SCIPABORT();
      }
   }

   /* recover previous information */
   root                   = oldroot;
   root->walkcurrentchild = oldcurrentchild;
   root->walkparent       = oldparent;
   root->walkio           = oldwalkio;

   /* release root captured by walker */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &root) );

   return SCIP_OKAY;
}

/** Gives the parent of an expression in an expression graph walk.
 *
 * During an expression walk, this function returns the expression from which the given expression has been accessed.
 * If not in an expression walk, the returned pointer is undefined.
 */
SCIP_CONSEXPR_EXPR* SCIPgetConsExprExprWalkParent(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression which parent to get */
   )
{
   assert(expr != NULL);

   return expr->walkparent;
}

/** Gives the index of the child that will be visited next (or is currently visited) by an expression walk. */
int SCIPgetConsExprExprWalkCurrentChild(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression which nextchild to get */
   )
{
   assert(expr != NULL);

   return expr->walkcurrentchild;
}

/** Gives the precedence of the expression handler of the parent expression in an expression graph walk.
 *
 * If there is no parent, then 0 is returned.
 */
unsigned int SCIPgetConsExprExprWalkParentPrecedence(
   SCIP_CONSEXPR_EXPR*   expr                /**< expression which parent to get */
   )
{
   assert(expr != NULL);

   if( expr->walkparent == NULL )
      return 0;

   return expr->walkparent->exprhdlr->precedence;
}

/*
 * constraint specific interface methods
 */

/** create and include conshdlr to SCIP and set everything except for expression handlers */
static
SCIP_RETCODE includeConshdlrExprBasic(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create expr constraint handler data */
   SCIP_CALL( SCIPallocClearMemory(scip, &conshdlrdata) );
   conshdlrdata->lastsoltag = 1;

   conshdlr = NULL;

   /* include constraint handler */
#if 0
   /* use SCIPincludeConshdlr() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlr(scip, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY,
         CONSHDLR_SEPAFREQ, CONSHDLR_PROPFREQ, CONSHDLR_EAGERFREQ, CONSHDLR_MAXPREROUNDS,
         CONSHDLR_DELAYSEPA, CONSHDLR_DELAYPROP, CONSHDLR_NEEDSCONS,
         CONSHDLR_PROP_TIMING, CONSHDLR_PRESOLTIMING,
         conshdlrCopyExpr,
         consFreeExpr, consInitExpr, consExitExpr,
         consInitpreExpr, consExitpreExpr, consInitsolExpr, consExitsolExpr,
         consDeleteExpr, consTransExpr, consInitlpExpr,
         consSepalpExpr, consSepasolExpr, consEnfolpExpr, consEnfopsExpr, consCheckExpr,
         consPropExpr, consPresolExpr, consRespropExpr, consLockExpr,
         consActiveExpr, consDeactiveExpr,
         consEnableExpr, consDisableExpr, consDelvarsExpr,
         consPrintExpr, consCopyExpr, consParseExpr,
         consGetVarsExpr, consGetNVarsExpr, consGetDiveBdChgsExpr, conshdlrdata) );
#else
   /* use SCIPincludeConshdlrBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpExpr, consEnfopsExpr, consCheckExpr, consLockExpr,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveExpr) );
   SCIP_CALL( SCIPsetConshdlrCopy(scip, conshdlr, conshdlrCopyExpr, consCopyExpr) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveExpr) );
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteExpr) );
   SCIP_CALL( SCIPsetConshdlrDelvars(scip, conshdlr, consDelvarsExpr) );
   SCIP_CALL( SCIPsetConshdlrDisable(scip, conshdlr, consDisableExpr) );
   SCIP_CALL( SCIPsetConshdlrEnable(scip, conshdlr, consEnableExpr) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitExpr) );
   SCIP_CALL( SCIPsetConshdlrExitpre(scip, conshdlr, consExitpreExpr) );
   SCIP_CALL( SCIPsetConshdlrExitsol(scip, conshdlr, consExitsolExpr) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeExpr) );
   SCIP_CALL( SCIPsetConshdlrGetDiveBdChgs(scip, conshdlr, consGetDiveBdChgsExpr) );
   SCIP_CALL( SCIPsetConshdlrGetVars(scip, conshdlr, consGetVarsExpr) );
   SCIP_CALL( SCIPsetConshdlrGetNVars(scip, conshdlr, consGetNVarsExpr) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitExpr) );
   SCIP_CALL( SCIPsetConshdlrInitpre(scip, conshdlr, consInitpreExpr) );
   SCIP_CALL( SCIPsetConshdlrInitsol(scip, conshdlr, consInitsolExpr) );
   SCIP_CALL( SCIPsetConshdlrInitlp(scip, conshdlr, consInitlpExpr) );
   SCIP_CALL( SCIPsetConshdlrParse(scip, conshdlr, consParseExpr) );
   SCIP_CALL( SCIPsetConshdlrPresol(scip, conshdlr, consPresolExpr, CONSHDLR_MAXPREROUNDS, CONSHDLR_PRESOLTIMING) );
   SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintExpr) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropExpr, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrResprop(scip, conshdlr, consRespropExpr) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpExpr, consSepasolExpr, CONSHDLR_SEPAFREQ, CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransExpr) );
#endif

   if( SCIPfindConshdlr(scip, "quadratic") != NULL )
   {
      /* include function that upgrades quadratic constraint to expr constraints */
      SCIP_CALL( SCIPincludeQuadconsUpgrade(scip, quadconsUpgdExpr, QUADCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   if( SCIPfindConshdlr(scip, "nonlinear") != NULL )
   {
      /* include the linear constraint upgrade in the linear constraint handler */
      SCIP_CALL( SCIPincludeNonlinconsUpgrade(scip, nonlinconsUpgdExpr, NULL, NONLINCONSUPGD_PRIORITY, TRUE, CONSHDLR_NAME) );
   }

   /* add expr constraint handler parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a set of constraints within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 10, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/minviolationsepa",
         "minimal violation for a cut to be added to the LP during separation; overwrites separating/efficacy",
         &conshdlrdata->mincutviolationsepa, TRUE, 0.0001, 0.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "constraints/" CONSHDLR_NAME "/minviolationenfofac",
         "minimal target violation of a cut in order to add it to relaxation during enforcement as a factor of the feasibility tolerance (may be ignored)",
         &conshdlrdata->mincutviolationenfofac, TRUE, 2.0, 1.0, SCIPinfinity(scip), NULL, NULL) );

   /* include handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &conshdlrdata->eventhdlr, CONSHDLR_NAME "_boundchange",
         "signals a bound change to an expression constraint", processVarEvent, NULL) );
   assert(conshdlrdata->eventhdlr != NULL);

   return SCIP_OKAY;
}


/** creates the handler for expr constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrExpr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   SCIP_CALL( includeConshdlrExprBasic(scip) );

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* include and remember handler for variable expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrVar(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "var") == 0);
   conshdlrdata->exprvarhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for constant value expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrValue(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "val") == 0);
   conshdlrdata->exprvalhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for sum expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrSum(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "sum") == 0);
   conshdlrdata->exprsumhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include and remember handler for product expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrProduct(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "prod") == 0);
   conshdlrdata->exprprodhdlr = conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1];

   /* include handler for exponential expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrExp(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "exp") == 0);

   /* include handler for logarithmic expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrLog(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "log") == 0);

   /* include handler for absolute expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrAbs(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "abs") == 0);

   /* include handler for power expression */
   SCIP_CALL( SCIPincludeConsExprExprHdlrPow(scip, conshdlr) );
   assert(conshdlrdata->nexprhdlrs > 0 && strcmp(conshdlrdata->exprhdlrs[conshdlrdata->nexprhdlrs-1]->name, "pow") == 0);

   return SCIP_OKAY;
}

/** creates and captures a expr constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs,                /**< right hand side of constraint */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   /* TODO: (optional) modify the definition of the SCIPcreateConsExpr() call, if you don't need all the information */

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   assert(expr != NULL);

   /* find the expr constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("expr constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   SCIP_CALL( SCIPallocClearBlockMemory(scip, &consdata) );
   consdata->expr = expr;
   consdata->lhs = lhs;
   consdata->rhs = rhs;

   /* capture expression */
   SCIPcaptureConsExprExpr(consdata->expr);

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a expr constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsExprBasic(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression of constraint (must not be NULL) */
   SCIP_Real             lhs,                /**< left hand side of constraint */
   SCIP_Real             rhs                 /**< right hand side of constraint */
   )
{
   SCIP_CALL( SCIPcreateConsExpr(scip, cons, name, expr, lhs, rhs,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** returns the expression of the given expression constraint */
SCIP_CONSEXPR_EXPR* SCIPgetExprConsExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< constraint data */
   )
{
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);

   if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(cons)), CONSHDLR_NAME) != 0 )
   {
      SCIPerrorMessage("constraint is not expression\n");
      SCIPABORT();
      return NULL;  /*lint !e527*/
   }

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->expr;
}

/** Creates an expression from a string.
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
 * See also @ref parseExpr.
 */
SCIP_RETCODE SCIPparseConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        consexprhdlr,       /**< expression constraint handler */
   const char*           exprstr,            /**< string with the expr to parse */
   const char**          finalpos,           /**< buffer to store the position of exprstr where we finished reading, or NULL if not of interest */
   SCIP_CONSEXPR_EXPR**  expr                /**< pointer to store the expr parsed */
   )
{
   const char* finalpos_;
   SCIP_RETCODE retcode;
   SCIP_HASHMAP* vartoexprvarmap;

   SCIP_CALL( SCIPhashmapCreate(&vartoexprvarmap, SCIPblkmem(scip), SCIPcalcHashtableSize(5 * SCIPgetNVars(scip))) );

   /* if parseExpr fails, we still want to free hashmap */
   retcode = parseExpr(scip, consexprhdlr, vartoexprvarmap, exprstr, &finalpos_, expr);

   SCIPhashmapFree(&vartoexprvarmap);

   if( finalpos != NULL )
      *finalpos = finalpos_;

   return retcode;
}

/** appends child to the children list of expr */
SCIP_RETCODE SCIPappendConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression */
   SCIP_CONSEXPR_EXPR*   child               /**< expression to be appended */
   )
{
   ENSUREBLOCKMEMORYARRAYSIZE(scip, expr->children, expr->childrensize, expr->nchildren + 1);

   expr->children[expr->nchildren] = child;
   ++expr->nchildren;

   /* capture child */
   SCIPcaptureConsExprExpr(child);

   /* update locks in child */
   SCIP_CALL( propagateLocks(scip, child, expr->nlockspos, expr->nlocksneg) );

   return SCIP_OKAY;
}

/** duplicates the given expression
 *
 * If a copy could not be created (e.g., due to missing copy callbacks in expression handlers), *copyexpr will be set to NULL.
 */
SCIP_RETCODE SCIPduplicateConsExprExpr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< original expression */
   SCIP_CONSEXPR_EXPR**  copyexpr            /**< buffer to store duplicate of expr */
   )
{
   COPY_DATA copydata;

   copydata.targetscip = scip;
   copydata.mapvar = NULL;
   copydata.mapvardata = NULL;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, copyExpr, NULL, copyExpr, copyExpr, &copydata) );
   *copyexpr = copydata.targetexpr;

   return SCIP_OKAY;
}

/** simplifies an expression
 * The given expression will be released and overwritten with the simplified expression.
 * To keep the expression, duplicate it via SCIPduplicateConsExprExpr before calling this method.
 */
SCIP_RETCODE SCIPsimplifyConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**    simplified        /**< buffer to store simplified expression */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(simplified != NULL);

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, NULL, simplifyExpr, simplifyExpr, (void*)simplified) );
   assert(*simplified != NULL);

   return SCIP_OKAY;
}

/** prints structure of an expression a la Maple's dismantle */
SCIP_RETCODE SCIPdismantleConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression to dismantle */
   )
{
   int depth;

   depth = -1;
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, dismantleExpr, dismantleExpr, NULL, dismantleExpr, &depth) );
   assert(depth == -1);

   return SCIP_OKAY;
}

/** overwrites/replaces a child of an expressions
 *
 * @note the old child is released and the newchild is captured
 */
SCIP_RETCODE SCIPreplaceConsExprExprChild(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr,             /**< expression which is going to replace a child */
   int                     childidx,         /**< index of child being replaced */
   SCIP_CONSEXPR_EXPR*     newchild          /**< the new child */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(newchild != NULL);
   assert(childidx < SCIPgetConsExprExprNChildren(expr));

   /* capture new child (do this before releasing the old child in case there are equal */
   SCIPcaptureConsExprExpr(newchild);

   /* update locks in old child */
   SCIP_CALL( propagateLocks(scip, expr->children[childidx], -expr->nlockspos, -expr->nlocksneg) );

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(expr->children[childidx])) );
   expr->children[childidx] = newchild;

   /* update locks in new child */
   SCIP_CALL( propagateLocks(scip, expr->children[childidx], expr->nlockspos, expr->nlocksneg) );

   return SCIP_OKAY;
}

/* maybe should make this a parameter (was cutmaxrange in other conshdlr)
 * maybe should derive this from the current feastol (e.g., 10/feastol)
 */
#define SCIP_CONSEXPR_CUTMAXRANGE 1.0e7

/** checks a cut for violation and numerical stability and possibly tries to improve it
 *
 * If the numerical properties of the cut are too bad, the routines tries to improve this.
 * If the violation of the cut in the given solution will end up to be below the given minviolation,
 * the cut will be released.
 * Passing -SCIPinfinity(scip) as minviolation will disable the violation check.
 */
SCIP_RETCODE SCIPmassageConsExprExprCut(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_ROW**              cut,              /**< cut to be checked and maybe modified */
   SCIP_SOL*               sol,              /**< solution that we try to cut off */
   SCIP_Real               minviolation      /**< minimal violation requirement (need to be nonnegative or -SCIPinfinity(scip)) */
   )
{
   SCIP_Real violation = SCIP_INVALID;
   SCIP_Real mincoef;
   SCIP_Real maxcoef;
   SCIP_SIDETYPE side;

   assert(scip != NULL);
   assert(cut != NULL);
   assert(*cut != NULL);
   assert(minviolation >= 0.0 || minviolation == -SCIPinfinity(scip));

   if( minviolation != -SCIPinfinity(scip) )
   {
      /* get current violation */
      violation = -SCIPgetRowSolFeasibility(scip, *cut, sol);

      /* release cut if its violation is too low */
      if( violation < minviolation )
      {
         SCIP_CALL( SCIPreleaseRow(scip, cut) );
         return SCIP_OKAY;
      }
   }

   /* check that there is either a lhs or a rhs (if not, then probably we were overflowing SCIPinfinity for one side) */
   if( SCIPisInfinity(scip, -SCIProwGetLhs(*cut)) && SCIPisInfinity(scip, SCIProwGetRhs(*cut)) )
   {
      SCIPdebugMessage("cut <%s> sides are both infinite: left %g right %g\n", SCIProwGetName(*cut), SCIProwGetLhs(*cut), SCIProwGetRhs(*cut));
      SCIP_CALL( SCIPreleaseRow(scip, cut) );
      return SCIP_OKAY;
   }

   /* assert that cut is an inequality */
   assert( SCIPisInfinity(scip, -SCIProwGetLhs(*cut)) ||  SCIPisInfinity(scip, SCIProwGetRhs(*cut)));
   /* check which side of the cut is not infinity */
   side = SCIPisInfinity(scip, SCIProwGetRhs(*cut)) ? SCIP_SIDETYPE_LEFT : SCIP_SIDETYPE_RIGHT;

   mincoef = SCIPgetRowMinCoef(scip, *cut);
   maxcoef = SCIPgetRowMaxCoef(scip, *cut);

   /* SCIPdebugMessage("cut <%s> violation %g mincoef %g maxcoef %g range %.2e side %g \n", SCIProwGetName(*cut), violation, mincoef, maxcoef, maxcoef/mincoef, side == SCIP_SIDETYPE_LEFT ? SCIProwGetLhs(*cut) : SCIProwGetRhs(*cut)); */

   /* if maximal coefficient is infinity, give up on the cut */
   if( SCIPisInfinity(scip, maxcoef) )
   {
      SCIPdebugMessage("maximal coefficients of cut <%s> is too large: %g\n", SCIProwGetName(*cut), maxcoef);
      SCIP_CALL( SCIPreleaseRow(scip, cut) );
      return SCIP_OKAY;
   }

   /* check and possibly try to improve cut range */
   while( maxcoef / mincoef > SCIP_CONSEXPR_CUTMAXRANGE )
   {
      SCIP_VAR* var;
      SCIP_Real coef;
      SCIP_Real constant;
      int j;

      /* if range of coefficients is bad, find very small coefficients and make them zero */
      SCIPdebugMessage("cut coefficients for cut <%s> have very large range: mincoef = %g maxcoef = %g\n", SCIProwGetName(*cut), mincoef, maxcoef);

      /* TODO in the original code, we were not looking into cases where the minimal coefficient is due to a linear variable.
       * This would correspond to the auxiliary variable here, which we might not want to eliminate from the cut.
       * Maybe we pass in a variable to this function which should not be removed from the cut?
       */

      /* eliminate variables with minimal coefficient from cut */
      constant = 0.0;
      for( j = 0; j < SCIProwGetNNonz(*cut); ++j )
      {
         coef = SCIProwGetVals(*cut)[j];
         if( !SCIPisEQ(scip, REALABS(coef), mincoef) )
            continue;

         var = SCIPcolGetVar(SCIProwGetCols(*cut)[j]);
         assert(var != NULL);

         /* try to eliminate coefficient with minimal absolute value by weakening cut and try again */
         if( (coef > 0.0 && side == SCIP_SIDETYPE_RIGHT) || (coef < 0.0 && side == SCIP_SIDETYPE_LEFT) )
         {
            SCIP_Real lb;

            lb = SCIProwIsLocal(*cut) ? SCIPvarGetLbLocal(var) : SCIPvarGetLbGlobal(var);
            if( !SCIPisInfinity(scip, -lb) )
            {
               SCIPdebugMessage("eliminate coefficient %g for <%s> = %g [>= %g]\n", coef, SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), lb);

               constant += coef * lb;
               SCIP_CALL( SCIPaddVarToRow(scip, *cut, var, -coef) );
               continue;
            }
         }

         if( (coef < 0.0 && side == SCIP_SIDETYPE_RIGHT) || (coef > 0.0 && side == SCIP_SIDETYPE_LEFT) )
         {
            SCIP_Real ub;

            ub = SCIProwIsLocal(*cut) ? SCIPvarGetUbLocal(var) : SCIPvarGetUbGlobal(var);
            if( !SCIPisInfinity(scip, ub) )
            {
               SCIPdebugMessage("eliminate coefficient %g for <%s> = %g [<= %g]\n", coef, SCIPvarGetName(var), SCIPgetSolVal(scip, sol, var), ub);

               constant += coef * ub;
               SCIP_CALL( SCIPaddVarToRow(scip, *cut, var, -coef) );
               continue;
            }
         }

         /* if variable could not be eliminated, then give up */
         SCIPdebugMessage("could not eliminate small coefficient %g of variable <%s>\n", coef, SCIPvarGetName(var));
         SCIP_CALL( SCIPreleaseRow(scip, cut) );
         return SCIP_OKAY;
      }

      /* adapt lhs/rhs */
      if( side == SCIP_SIDETYPE_LEFT )
      {
         SCIP_CALL( SCIPchgRowLhs(scip, *cut, SCIProwGetLhs(*cut) - constant) );
      }
      else
      {
         SCIP_CALL( SCIPchgRowRhs(scip, *cut, SCIProwGetRhs(*cut) - constant) );
      }

      /* update min/max coefficient */
      mincoef = SCIPgetRowMinCoef(scip, *cut);
      maxcoef = SCIPgetRowMaxCoef(scip, *cut);

      /* remember that we changed the cut */
      violation = SCIP_INVALID;
   }
   assert(maxcoef / mincoef < SCIP_CONSEXPR_CUTMAXRANGE);

   /* check that left/right hand side are finite */
   if( (side == SCIP_SIDETYPE_LEFT  && SCIPisInfinity(scip, -SCIProwGetLhs(*cut))) ||
       (side == SCIP_SIDETYPE_RIGHT && SCIPisInfinity(scip,  SCIProwGetRhs(*cut))) )
   {
      SCIPdebugMessage("cut <%s> has very large side: %g\n", SCIProwGetName(*cut), side == SCIP_SIDETYPE_LEFT ? -SCIProwGetLhs(*cut) : SCIProwGetRhs(*cut));
      SCIP_CALL( SCIPreleaseRow(scip, cut) );
      return SCIP_OKAY;
   }

   /* check violation again if cut was changed */
   if( violation == SCIP_INVALID && minviolation != -SCIPinfinity(scip) )
   {
      violation = -SCIPgetRowSolFeasibility(scip, *cut, sol);

      /* release cut if its violation is too low */
      if( violation < minviolation )
      {
         SCIP_CALL( SCIPreleaseRow(scip, cut) );
         return SCIP_OKAY;
      }
   }

   /* TODO we could scale up the cut (issue #7) if violation is >0 and <minviolation */

   return SCIP_OKAY;
}
