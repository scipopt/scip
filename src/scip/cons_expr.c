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

/** constraint data for expr constraints */
struct SCIP_ConsData
{
   SCIP_CONSEXPR_EXPR*   expr;               /**< expression that represents this constraint (must evaluate to 0 (FALSE) or 1 (TRUE)) */
   SCIP_Real             lhs;                /**< left-hand side */
   SCIP_Real             rhs;                /**< right-hand side */

   SCIP_Real             lhsviol;            /**< violation of left-hand side by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of right-hand side by current solution (used temporarily inside constraint handler) */
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

   unsigned int             lastsoltag;      /**< last solution tag used to evaluate current solution */
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
} EXPRINTEVAL_DATA;

/** data passed on during variable locking  */
typedef struct
{
   SCIP_CONSEXPR_EXPRHDLR* exprvarhdlr;      /**< handler for variable expressions (to recognize variable expressions) */
   int                     nlockspos;        /**< number of positive locks */
   int                     nlocksneg;        /**< number of negative locks */
} EXPRLOCK_DATA;

/** data passed on during copying expressions */
typedef struct
{
   SCIP*                   targetscip;                 /**< target SCIP pointer */
   SCIP_DECL_CONSEXPR_EXPRCOPYDATA_MAPVAR((*mapvar));  /**< variable mapping function, or NULL for identity mapping (used in handler for var-expressions) */
   void*                   mapvardata;                 /**< data of variable mapping function */
} COPY_DATA;

/** variable mapping data passed on during copying expressions when copying SCIP instances */
typedef struct
{
   SCIP_HASHMAP*           varmap;           /**< SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP */
   SCIP_HASHMAP*           consmap;          /**< SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP */
   SCIP_Bool               global;           /**< should a global or a local copy be created */
   SCIP_Bool               valid;            /**< indicates whether every variable copy was valid */
} COPY_MAPVAR_DATA;

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
            assert(targetexprdata != NULL);
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
      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      default:
      {
         SCIPABORT(); /* we should never be called in this stage */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }
   }
}

static
SCIP_DECL_CONSEXPR_EXPRCOPYDATA_MAPVAR(transformVar)
{
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
{
   assert(expr != NULL);

   /* the expression should not be used by more than the parent that is also freed */
   assert(expr->nuses == 1);

   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      {
         /* free expression data, if any, when entering expression */

         if( expr->exprdata != NULL && expr->exprhdlr->freedata != NULL )
         {
            SCIP_CALL( expr->exprhdlr->freedata(scip, expr) );
            assert(expr->exprdata == NULL);
         }

         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }

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
            *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         }

         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR :
      {
         /* free expression when leaving it */

         /* free children array, if any */
         SCIPfreeBlockMemoryArrayNull(scip, &expr->children, expr->childrensize);

         SCIPfreeBlockMemory(scip, &expr);
         assert(expr == NULL);

         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD :
      default:
      {
         SCIPABORT(); /* we should never be called in this stage */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
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

/** expression walk callback when evaluating expression, called before child is visited */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(evalExprVisitChild)
{
   EXPREVAL_DATA* evaldata;

   assert(expr != NULL);
   assert(data != NULL);

   evaldata = (EXPREVAL_DATA*)data;

   /* skip child if it has been evaluated for that solution already */
   if( evaldata->soltag != 0 && evaldata->soltag == expr->children[expr->walkcurrentchild]->evaltag )
   {
      if( expr->children[expr->walkcurrentchild]->evalvalue == SCIP_INVALID )
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
{
   EXPREVAL_DATA* evaldata;

   assert(expr != NULL);
   assert(data != NULL);
   assert(expr->exprhdlr->eval != NULL);

   evaldata = (EXPREVAL_DATA*)data;

   SCIP_CALL( (*expr->exprhdlr->eval)(scip, expr, &expr->evalvalue, evaldata->sol) );
   expr->evaltag = evaldata->soltag;

   if( expr->evalvalue == SCIP_INVALID )
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
{
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
{
   EXPRINTEVAL_DATA* propdata;

   assert(expr != NULL);
   assert(data != NULL);

   propdata = (EXPRINTEVAL_DATA*)data;

   /* set tag in any case */
   expr->intevaltag = propdata->boxtag;

   /* set interval to [-inf,+inf] if interval evaluation callback is not implemented */
   if( expr->exprhdlr->inteval == NULL )
   {
      SCIPintervalSetEntire(SCIPinfinity(scip), &expr->interval);
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

      return SCIP_OKAY;
   }

   /* evaluate current expression and move on */
   SCIP_CALL( (*expr->exprhdlr->inteval)(scip, expr, &expr->interval) );

   /* stop if callback returned an empty interval */
   if( SCIPintervalIsEmpty(SCIPinfinity(scip), expr->interval) )
   {
      propdata->aborted = TRUE;
      *result = SCIP_CONSEXPREXPRWALK_ABORT;
   }
   else
   {
      *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
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

   if( SCIPgetConsExprExprHdlr(expr) == lockdata->exprvarhdlr )
   {
      /* if a variable, lock in both directions */
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
      default:
      {
         /* shouldn't be here */
         SCIPABORT();
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }
   }
   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

   return SCIP_OKAY;
}

/** TODO: FIXME: DELETE expression walk callback to simplify an expression
 * simplifies bottom up */
static
SCIP_RETCODE SCIPsimplifyConsExprExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**  simplifiedexpr      /**< pointer to store the simplified expression */
   );
static
SCIP_RETCODE SCIPsimplifyConsExprExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**  simplifiedexpr      /**< pointer to store the simplified expression */
   );

/** expression walk callback to simplify an expression
 * simplifies bottom up; when leaving an expression it simplifies it and stores the simplified expr in its walkio ptr;
 * after the child was visited, it is replaced with the simpified expr
 */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(simplifyExpr)
{
   assert(expr != NULL);
   switch( stage )
   {
      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD:
      {
         int currentchild;

         currentchild = SCIPgetConsExprExprWalkCurrentChild(expr);

         SCIP_CALL( SCIPreplaceConsExprExprChild(scip, expr, currentchild,
                     (SCIP_CONSEXPR_EXPR*)expr->children[currentchild]->walkio.ptrval) );

         /* continue */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }
      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR:
      {
         SCIP_CONSEXPR_EXPR* simplifiedexpr;
         /* TODO: this should probably become SCIP_CALL( (*expr->exprhdlr->simplify)(scip, expr, &simplifiedexpr) ); */
         if( strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "sum") == 0 )
         {
            SCIPsimplifyConsExprExprSum(scip, expr, &simplifiedexpr);
         }
         else if( strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(expr)), "prod") == 0 )
         {
            SCIPsimplifyConsExprExprProduct(scip, expr, &simplifiedexpr);
         }
         else
         {
            simplifiedexpr = expr;
            SCIPcaptureConsExprExpr(simplifiedexpr);
         }
         expr->walkio.ptrval = (void *)simplifiedexpr;
         /* continue */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR :
      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD :
      default:
      {
         SCIPABORT(); /* we should never be called in this stage */
         *result = SCIP_CONSEXPREXPRWALK_CONTINUE;
         return SCIP_OKAY;
      }
   }
}

/**@} */  /* end of walking methods */

/** @name Simplifying methods
 *
 * This is largely inspired in Joel Cohen's
 * Computer algebra and symbolic computation: Mathematical methods
 * In particular Chapter 3
 * The other fountain of inspiration is the current simplifying methods in expr.c.
 *
 * Note: 1) some parts might not apply to what we want to do. This is just for recording
 *          the information and discussion
 *       2) one can think that a child of a product with a non integer exponent is another
 *          operator
 *
 * ============================================
 * Definition of simplified expressions
 * ============================================
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
 * - it is a function with simplified arguments
 * ? a logarithm doesn't have a product as a child
 * ? the exponent of an exponential is always 1
 *
 * ============================================
 * ORDER
 * ============================================
 * This is a partial order for *simplified* expressions. Just a copy of the order from
 * the book so feel very free to modify.
 *
 * Equal type expressions:
 * - u,v value expressions: u < v <=> val(u) < val(v)
 * - u,v var expressions: u < v <=> SCIPvarGetIndex(var(u)) < SCIPvarGetIndex(var(v)) <=> SCIPvarCompare(var(u),var(v))
 * - u,v are both sum or product expression: < is a lexicographical order on the terms,
 *    starting from the _last_ finds the first index i where they differ and u < v <=> u_i < v_i
 *    If they are the same in all indices, then u < v <=> nchildren(u) < nchildren(v)
 *       Note: we are assuming expression are simplified, so within u, we have u_1 < u_2, etc
 *       Example: y + z < x + y + z, 2*x + 3*y < 3*x + 3*y
 *       Question: Quadratics are one of the most important cases, does this make sense for quadratics?
 * - u, v expressions p,q numbers: u^q < v^p <=> u < v, and in case they are equal, q < p
 * - u, v are functional expressions (exp, log, etc): u < v <=> Kind(u) < Kind(v), or if they are equal, args(u) < args(v)
 *    (the first argument and so on), if all common arguments are equal, then the one with less arguments < other one
 *       Example: f(x) < f(y), g(x) < g(x,y)
 *       Note: Kind is the type of operator
 *
 * Different type expressions:
 * - u value expr, v other: u < v always
 * - u product, v sum, var or func: u < v <=> u < 1*v
 *       Note: This means we compare u with the 1*v product. Though 1*v is unsimplified, the rule applies.
 *       Example: 2*x^0.5 < x [Note that x is a var expression]
 * - u^p, v sum, var or func: u^p < v <=> u^p < v^1 (I think this is the same as the previous one)
 * - u sum, v var or func: u < v <=> u < 0+v
 * - u sum, v var or func: u < v <=> u < 0+v
 * - u var, v func: u < v always
 * - u, v and none of the rules apply: u < v <=> ! v < u
 *    Example: is x < x^2 ? x is var and x^2 product, so none applies, then
 *    we try to answer if x^2 < x <=> x^2 < x^1 <=> 2 < 1 <=> False, so x < x^2 is True
 *
 * ============================================
 * RULES
 * ============================================
 * - Distributive: sums up identical terms (a*u + b*u -> (a+b)*u)
 *                 multiplies up identical factors ( u^a * u^b -> u^(a+b) )
 *    Note: 1 + x + (1 + x) is not transformed into 2(1 + x), since the tree is
 *          +--------
 *          |   |   |
 *          1   x   +----
 *                  |   |
 *                  1   x
 *          and no two children are identical
 *
 * - Simple associative: sums/products don't have sums/products as children with coefficient/exponent 1
 *    Example: 1) x + (y + z) = x + y + z
 *             2) x * (y * z) = x * y * z
 *
 * - Associative: sums don't have sums as children
 *                products don't have products with integer exponent as children
 *    Example: 1) x + 2*(y + z) = x + 2*y + 2*z (if represented correctly)
 *             2) x * (y * z)^n = x * y^n * z^n
 *
 * - Commutative: sorts operands in sums and products in a specified order
 *    Example: 1) x*z*3*y = 3*x*y*z
 *
 * - Basic identities: u * 0 -> 0
 *                     u + 0 -> u
 *                     u * 1 -> u
 *                     0^p   -> 0 where p positive otherwise invalid expr (= infeasible?)
 *                     1^w   -> 1 with w whatever (!= infinity)
 *                     v^0   -> 1 (this is actually tricky, because v has to be != 0)
 *                     v^1   -> v
 *
 * - Numerical: no sum/product has more than one constant operand
 *              no function has only constant arguments
 *
 * - Unary: no sum/product has only one child with 1 as coefficient/exponent
 *
 * ============================================
 * Algorithm
 * ============================================
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

#ifdef SIMPLIFY_DEBUG
#define debugSimplify                   printf
#else
#define debugSimplify                   while( FALSE ) printf
#endif

/** node for linked list of expressions */
typedef struct exprnode
{
   SCIP_CONSEXPR_EXPR*   expr;               /**< expression in node */
   SCIP_Real             coef;               /**< coefficient or exponent of expr*/
   struct exprnode*      next;               /**< next node */
} EXPRNODE;

static
void insertFirstList(EXPRNODE* newnode, EXPRNODE** list)
{
   newnode->next = *list;
   *list = newnode;
}

static
EXPRNODE* listPopFirst(EXPRNODE** list)
{
   EXPRNODE* first;

   if( *list == NULL )
      return NULL;

   first = *list;
   *list = (*list)->next;
   first->next = NULL;

   return first;
}

static
int listLength(EXPRNODE* list)
{
   int length;

   if( list == NULL )
      return 0;

   length = 1;
   while( (list=list->next) != NULL )
      ++length;

   return length;
}

/** creates node and captures expression */
static
SCIP_RETCODE createExprNode(SCIP* scip, SCIP_CONSEXPR_EXPR* expr, SCIP_Real coef, EXPRNODE** newnode)
{
   SCIP_CALL( SCIPallocBuffer(scip, newnode) );
   (*newnode)->expr = expr;
   (*newnode)->coef = coef;
   (*newnode)->next = NULL;
   SCIPcaptureConsExprExpr(expr);

   return SCIP_OKAY;
}

/** frees node and releases expression */
static
SCIP_RETCODE freeExprNode(SCIP* scip, EXPRNODE** newnode)
{
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*newnode)->expr) );
   SCIPfreeBuffer(scip, newnode);

   return SCIP_OKAY;
}

/** frees list and releases expressions */
static
SCIP_RETCODE freeExprlist(SCIP* scip, EXPRNODE** exprlist)
{
   EXPRNODE* current;

   if( *exprlist == NULL )
      return SCIP_OKAY;

   current = *exprlist;
   while( current != NULL )
   {
      EXPRNODE* aux;

      aux = current->next;
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(current->expr)) );
      SCIPfreeBuffer(scip, &current);
      current = aux;
   }

   assert(current == NULL);
   *exprlist = NULL;

   return SCIP_OKAY;
}

static
SCIP_RETCODE buildExprlistFromExprs(SCIP* scip, SCIP_CONSEXPR_EXPR** exprs, SCIP_Real* coefs, int nexprs, EXPRNODE** list)
{
   int i;

   assert(nexprs > 0);
   assert(*list == NULL);

   debugSimplify("building expr list from %d expressions\n", nexprs);
   for( i = nexprs - 1; i >= 0; --i )
   {
      EXPRNODE* newnode;

      SCIP_CALL( createExprNode(scip, exprs[i], coefs[i], &newnode) );
      insertFirstList(newnode, list);
   }

   assert(nexprs > 1 || (*list)->next == NULL);

   return SCIP_OKAY;
}

/** merges tomerge into finalchildren
 *
 * Both, tomerge and finalchildren contain expressions that could be the children of a simplified sum
 * (except for SS6 and SS7 which are enforced later).
 * However, the concatenation of both lists, will not in general yield a simplified sum expression,
 * because both SS4 and SS5 could be violated. So the purpose of this method is to enforce SS4 and SS5.
 * In the process of enforcing SS4, it could happen that SS8 is violated, but this is easy to fix.
 * note: if list has more than one element, then they are the children of a simplified sum expression
 * (ie, no val, nor sum), otherwise is a product, a variable or a function node
 */
static
SCIP_RETCODE mergeSumExprlist(SCIP* scip, EXPRNODE* tomerge, EXPRNODE** finalchildren)
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

      debugSimplify("comparing exprs:\n");
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
            debugSimplify("GOT 0 WHILE ADDING UP\n");
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

   /* there are still nodes of tomerge unmerged; these nodes are larger than finalchildren, so append at end */
   assert(current == NULL);
   assert(previous != NULL && previous->next == NULL);
   previous->next = tomergenode;

   return SCIP_OKAY;
}


/** builds a sum expression with the elements of exprlist as its children */
static
SCIP_RETCODE buildExprSumFromExprlist(SCIP* scip, EXPRNODE* exprlist, SCIP_Real constant, SCIP_CONSEXPR_EXPR** expr)
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

/* TODO: both functions are the same (buildExprSum|ProductFromExprlist), do something about it */
static
SCIP_RETCODE buildExprProductFromExprlist(SCIP* scip, EXPRNODE* exprlist, SCIP_Real coef, SCIP_CONSEXPR_EXPR** expr)
{
   int i;
   int nchildren;
   SCIP_CONSEXPR_EXPR** children;
   SCIP_Real* exponents;

   nchildren = listLength(exprlist);

   SCIP_CALL( SCIPallocBufferArray(scip, &exponents, nchildren) );
   SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );

   for( i = 0; i < nchildren; ++i )
   {
      children[i] = exprlist->expr;
      exponents[i] = exprlist->coef;
      exprlist = exprlist->next;
   }

   assert(exprlist == NULL);

   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, SCIPfindConshdlr(scip, "expr"), expr, nchildren, children, exponents, coef) );

   SCIPfreeBufferArray(scip, &children);
   SCIPfreeBufferArray(scip, &exponents);

   return SCIP_OKAY;
}

/** simplifies a term expression constant * expr
 * @note: in contrast to other simplify methods, this does *not* return a simplified expression.
 * Instead, the method is intended to be called only when simplifying a sum expression,
 * Since in general, constant*expr is not a simplified child of a sum expression, this method returns
 * a list of expressions L, such that (sum L) = constant * expr *and* each expression in L
 * is a valid child of a simplified sum expression.
 */
static
SCIP_RETCODE simplifyTerm(
   SCIP*                 scip,
   SCIP_CONSEXPR_EXPR*   expr,
   SCIP_Real             coef,
   SCIP_Real*            simplifiedconstant,
   EXPRNODE**            simplifiedterm
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

   /* enforces SS2 */
   if( strcmp(exprtype, "sum") == 0 )
   {
      *simplifiedconstant += coef * SCIPgetConsExprExprSumConstant(expr);
      SCIPsetConsExprExprSumConstant(expr, 0.0);

      /* if the coef of this expr is not 1.0, we have to multiply (distribute) it with coef;
       * this can render expr unsimplified (e.g., 2 * (1/2 * x) -> x)
       * we simplify it and re-call the method with the simplified expr since its type could have changed
       * Otherwise, just build a list with the children to merge it to the finalchildren
       */
      if( coef != 1.0 )
      {
         /* at this point, expr is a simplified sum with constant 0: expr = (sum coef1 expr1 coef2 expr2...). We have to
          * multiply expr by coef before grabbing its children for merging (SS2). After multiplying we obtain
          * expr' = (sum coef1' expr1 coef2' expr2...) which will clearly satisfy SS1-SS4, SS6 and SS8.
          * SS5 is satisfied, because if coef1 expr1 < coef2 expr2 are children in a simplified sum, then expr1 != expr2.
          * Therefore expr1 < expr2, which implies that C1 * expr1 < C2 * expr2 for any C1, C2 different from 0
          * So the only condition that can fail is SS7. In that case, expr = (sum coef1 expr1) and expr' = (sum 1 expr1)
          * and so simplifying expr' gives  expr1
          */
         SCIPmultiplyConsExprExprSumByConstant(expr, coef);
         if( SCIPgetConsExprExprNChildren(expr) == 1 && SCIPgetConsExprExprSumCoefs(expr)[0] == 1.0 )
         {
            /* child of expr (a simplified sum expression) cannot be another sum */
            assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(SCIPgetConsExprExprChildren(expr)[0])), "sum") != 0);
            /* this is not really a recursive call */
            SCIP_CALL( simplifyTerm(scip, SCIPgetConsExprExprChildren(expr)[0], 1.0, simplifiedconstant, simplifiedterm) );
            return SCIP_OKAY;
         }
      }
      SCIP_CALL( buildExprlistFromExprs(scip, SCIPgetConsExprExprChildren(expr),
               SCIPgetConsExprExprSumCoefs(expr), SCIPgetConsExprExprNChildren(expr), simplifiedterm) );

      return SCIP_OKAY;
   }
   else
   {
      /* other types of (simplified) expressions can be children of a simplified sum */
      assert(strcmp(exprtype, "sum") != 0);
      assert(strcmp(exprtype, "val") != 0);

      SCIP_CALL( createExprNode(scip, expr, coef, simplifiedterm) );
      return SCIP_OKAY;
   }

   return SCIP_OKAY;
}

/** simplifies a product expression (base)^exponent
 * @note: in contrast to other simplify methods, this does *not* return a simplified expression.
 * Instead, the method is intended to be called only when simplifying a product expression,
 * Since in general, (base)^exponent is not a simplified child of a product expression, this method returns
 * a list of expressions (with exponents) L, such that (prod L) = (base)^exponent *and* each expression in L
 * is a valid child of a simplified product expression.
 * TODO: handle more cases when exponent is non-integer and base >= 0 (of course one has to adapt the definition
 * of a simplified product)
 */
static
SCIP_RETCODE simplifyPower(
   SCIP*                 scip,
   SCIP_CONSEXPR_EXPR*   base,
   SCIP_Real             exponent,
   SCIP_Real*            simplifiedcoef,
   EXPRNODE**            simplifiedpower
   )
{
   const char* basetype;
   SCIP_CONSEXPR_EXPR* simplifiedbase;

   assert(simplifiedpower != NULL);
   assert(*simplifiedpower == NULL);
   assert(base != NULL);

   basetype = SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(base));

   /* enforces SP7 */
   if( strcmp(basetype, "val") == 0 )
   {
      /* TODO: if val < 0 and exponent non integer -> domain error/undefined etc */
      debugSimplify("[simplifyPower] seeing value %g, exponent %g -> include in coef\n", SCIPgetConsExprExprValueValue(base), exponent);
      *simplifiedcoef *= pow(SCIPgetConsExprExprValueValue(base), exponent);
      return SCIP_OKAY;
   }

   /* currently, for a non-value base, we only try to do something for integer exponents
    * TODO: maybe put this as an extra condition in each sub-case so that later is easier to extend */
   if( !SCIPisIntegral(scip, exponent) )
   {
      debugSimplify("[simplifyPower] seeing some expr with non-integer exponent %g -> potential child\n", exponent);
      SCIP_CALL( createExprNode(scip, base, exponent, simplifiedpower) );
      return SCIP_OKAY;
   }
   /* round exponent so that is actually an integer */
   exponent = SCIPround(scip, exponent);

   /* enforces SP6
    * (base)^0 return empty list, which is the same as value 1
    */
   if( exponent == 0.0 )
   {
      debugSimplify("[simplifyPower] exponent %g (zero), ignore child\n", exponent);
      return SCIP_OKAY;
   }

   /* enforces SP2 */
   if( strcmp(basetype, "prod") == 0 )
   {
      assert(SCIPgetConsExprExprProductCoef(base) == 1.0);
      debugSimplify("[simplifyPower] seing a producut with exponent %g: include its children\n", exponent);

      /* if base is a product and exponent is not 1, we distribute the exponent among the base children.
       * this operation can render base unsimplified (e.g., ((x^0.5 * y^0.5)^0.5)^4 -> (x^0.5 * y^0.5)^2).
       * we simplify it and re-call the method with the simplified base since its type could have changed
       * e.g., (<any_expr>^0.5)^2 is a product-> <any_expr> which is any expression
       */
      if( exponent != 1.0 )
      {
         SCIPexponentiateConsExprExprProductByConstant(base, exponent);
         SCIP_CALL( SCIPsimplifyConsExprExprProduct(scip, base, &simplifiedbase) );
         SCIP_CALL( simplifyPower(scip, simplifiedbase, 1.0, simplifiedcoef, simplifiedpower) );
         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &simplifiedbase) );
      }
      else
      {
         SCIP_CALL( buildExprlistFromExprs(scip, SCIPgetConsExprExprChildren(base),
                  SCIPgetConsExprExprProductExponents(base), SCIPgetConsExprExprNChildren(base), simplifiedpower) );
      }

      return SCIP_OKAY;
   }

   /* enforces SP3
    * given (prod C ... n (sum 0.0 coef expr) ...) we can take coef out of the sum:
    * (prod C*coef^n ... n (sum 0.0 1 expr) ...) -> (prod C*coef^n ... n expr ...)
    * se we have to update simplifiedcoef and base = (sum 0.0 coef expr) changes to expr
    * notes: - since base is simplified and its constant is 0, then coef != 1.0 (SS7)
    *        - n is an integer (including 1, but not 0; see SP6 above)
    */
   if( strcmp(basetype, "sum") == 0 && SCIPgetConsExprExprNChildren(base) == 1 && SCIPgetConsExprExprSumConstant(base) == 0.0 )
   {
      debugSimplify("[simplifyPower] seing a sum with one term, exponent %g: base should be its child\n", exponent);
      /* assert SS7 holds */
      assert(SCIPgetConsExprExprSumCoefs(base)[0] != 1.0);

      /* update simplifiedcoef and simplify new base */
      *simplifiedcoef *= pow(SCIPgetConsExprExprSumCoefs(base)[0], exponent);
      SCIP_CALL( simplifyPower(scip, SCIPgetConsExprExprChildren(base)[0], exponent, simplifiedcoef, simplifiedpower) );

      return SCIP_OKAY;
   }
   else
   {
      /* other types of (simplified) expressions can be children of a simplified sum */
      assert(strcmp(basetype, "prod") != 0);
      assert(strcmp(basetype, "val") != 0);
      SCIP_CALL( createExprNode(scip, base, exponent, simplifiedpower) );
      return SCIP_OKAY;
   }

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
 * it can contain products, but only because they are acting as powers!
 * TODO: this function and mergeSumExprlist are very similar... merge them?
 */
static
SCIP_RETCODE mergeProductExprlist(SCIP* scip, EXPRNODE* tomerge, EXPRNODE** finalchildren, EXPRNODE** unsimplifiedchildren)
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

      /* assert invariants */
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(tomergenode->expr)), "val") != 0);
      assert(strcmp(SCIPgetConsExprExprHdlrName(SCIPgetConsExprExprHdlr(current->expr)), "val") != 0);
      assert(previous == NULL || previous->next == current);

      compareres = SCIPcompareConsExprExprs(current->expr, tomergenode->expr);
      if( compareres == 0 )
      {
         /* enforces SP4 */
         current->coef += tomergenode->coef;

         /* destroy tomergenode (since will not use the node again) */
         aux = tomergenode;
         tomergenode = tomergenode->next;
         SCIP_CALL( freeExprNode(scip, &aux) );

         /* the product might have render current unsimplified, so remove it from finalchildren list and put it back
          * to unsimplified children
          * TODO: we can say if current became unsimplified and don't touch it if it didn't
          */
         /* remove: if current is the first node, then pop; otherwise, use previous and current to remove */
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

         insertFirstList(aux, unsimplifiedchildren);
      }
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

   /* there are still nodes of tomerge unmerged; these nodes are larger than finalchildren, so append at end */
   assert(current == NULL);
   assert(previous != NULL && previous->next == NULL);
   previous->next = tomergenode;

   return SCIP_OKAY;
}

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
SCIP_RETCODE SCIPsimplifyConsExprExprSum(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**  simplifiedexpr      /**< pointer to store the simplified expression */
   )
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
   debugSimplify("what to do? finalchildren = %p has length %d\n", (void *)finalchildren, listLength(finalchildren));

   /* enforces SS6: if list is empty, return value */
   if( finalchildren == NULL )
   {
      debugSimplify("got empty list, return value\n");
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
      SCIP_CALL( buildExprSumFromExprlist(scip, finalchildren, simplifiedconstant, simplifiedexpr) );
   }

   /* free memory */
   freeExprlist(scip, &finalchildren);
   assert(finalchildren == NULL);

   return SCIP_OKAY;
}

/** simplifies a product expression
 *
 * Summary: we first build a list of expressions (called finalchildren) which will be the children of the simplified product
 * and then we process this list in order to enforce SP8-10
 * Description: In order to build finalchildren, we first build list of unsimplified children (called unsimplifiedchildren)
 * with the children of the product. Each node of the list is manipulated (see simplifyPower) in order to satisfy
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
SCIP_RETCODE SCIPsimplifyConsExprExprProduct(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*   expr,               /**< expression to be simplified */
   SCIP_CONSEXPR_EXPR**  simplifiedexpr      /**< pointer to store the simplified expression */
   )
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
   SCIP_CALL( buildExprlistFromExprs(scip, SCIPgetConsExprExprChildren(expr),
            SCIPgetConsExprExprProductExponents(expr), SCIPgetConsExprExprNChildren(expr), &unsimplifiedchildren) );

   /* while there are still children to process */
   finalchildren  = NULL;
   simplifiedcoef = SCIPgetConsExprExprProductCoef(expr);
   while( unsimplifiedchildren != NULL )
   {
      EXPRNODE* tomerge;
      EXPRNODE* first;

      /* if the simplified coefficient is 0, we can return value 0 */
      if( simplifiedcoef == 0.0 )
      {
         freeExprlist(scip, &finalchildren);
         freeExprlist(scip, &unsimplifiedchildren);
         assert(finalchildren == NULL);
         break;
      }

      first = listPopFirst(&unsimplifiedchildren);
      assert(first != NULL);

      /* enforces SP2, SP3, SP6 and SP7 */
      tomerge = NULL;
      SCIP_CALL( simplifyPower(scip, first->expr, first->coef, &simplifiedcoef, &tomerge) );

      /* enforces SP4 and SP5
       * note: merge frees (or uses) the nodes of the list tomerge */
      SCIP_CALL( mergeProductExprlist(scip, tomerge, &finalchildren, &unsimplifiedchildren) );

      /* free first */
      SCIP_CALL( freeExprlist(scip, &first) );
   }

   /* build product expression from finalchildren and post-simplify */
   debugSimplify("what to do? finalchildren = %p has length %d\n", (void *)finalchildren, listLength(finalchildren));

   /* enforces SP10: if list is empty, return value */
   if( finalchildren == NULL )
   {
      debugSimplify("got empty list, return value\n");
      SCIP_CALL( SCIPcreateConsExprExprValue(scip, SCIPfindConshdlr(scip, "expr"), simplifiedexpr, simplifiedcoef) );
   }
   /* enforces SP9
    * if finalchildren has only one expr and its exponent is 1.0 and coef 1.0, return that expr */
   else if( finalchildren->next == NULL && finalchildren->coef == 1.0 && simplifiedcoef == 1.0 )
   {
      *simplifiedexpr = finalchildren->expr;
      SCIPcaptureConsExprExpr(*simplifiedexpr);
   }
   /* enforces SP8 and SP9
    * if finalchildren has only one expr and its exponent is 1.0, but coef != 1.0, return sum with the expr as unique child */
   else if( finalchildren->next == NULL && finalchildren->coef == 1.0 )
   {
      SCIP_CONSEXPR_EXPR* aux;

      SCIP_CALL( SCIPcreateConsExprExprSum(scip, SCIPfindConshdlr(scip, "expr"), &aux,
               1, &(finalchildren->expr), &simplifiedcoef, 0.0) );

      /* simplifying here is necessary, the product could have sums as children
       * e.g., (prod 2 (sum 1 <x>)) -> (sum 2 (sum 1 <x>)) and that needs to be simplified
       */
      SCIP_CALL( SCIPsimplifyConsExprExprSum(scip, aux, simplifiedexpr) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );
   }
   /* enforces SP8
    * if simplifiedcoef != 1.0, transform it into a sum with the (simplified) product as child */
   else if( simplifiedcoef != 1.0 )
   {
      SCIP_CONSEXPR_EXPR* aux;

      SCIP_CALL( buildExprProductFromExprlist(scip, finalchildren, 1.0, &aux) );
      SCIP_CALL( SCIPcreateConsExprExprSum(scip, SCIPfindConshdlr(scip, "expr"), simplifiedexpr,
               1, &aux, &simplifiedcoef, 0.0) );
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );
   }
   /* build product expression from list */
   else
   {
      SCIP_CALL( buildExprProductFromExprlist(scip, finalchildren, simplifiedcoef, simplifiedexpr) );
   }

   /* free memory */
   freeExprlist(scip, &finalchildren);
   assert(finalchildren == NULL);

   assert(*simplifiedexpr != NULL);
   return SCIP_OKAY;
}

/*
 * ============================================
 * ORDER
 * ============================================
 * This is a partial order for *simplified* expressions. Just a copy of the order from
 * the book so feel very free to modify.
 * Comparing equal type expressions:
 * - u,v value expressions: u < v <=> val(u) < val(v) [DONE]
 * - u,v var expressions: u < v <=> SCIPvarGetIndex(var(u)) < SCIPvarGetIndex(var(v)) <=> SCIPvarCompare(var(u),var(v)) [DONE]
 * - u,v are both sum or product expression: < is a lexicographical order on the terms, [DONE]
 *    starting from the _last_ finds the first index i where they differ and u < v <=> u_i < v_i
 *    If they are the same in all indices, then u < v <=> nchildren(u) < nchildren(v)
 *       Note: we are assuming expression are simplified, so within u, we have u_1 < u_2, etc
 *       Example: y + z < x + y + z, 2*x + 3*y < 3*x + 3*y
 *       Question: Quadratics are one of the most important cases, does this make sense for quadratics?
 * - u, v expressions p,q numbers: u^q < v^p <=> u < v, and in case they are equal, q < p [DONE:considered in the above case]
 * - u, v are functional expressions (exp, log, etc): u < v <=> Kind(u) < Kind(v), or if they are equal, args(u) < args(v) [not quite]
 *    (the first argument and so on), if all common arguments are equal, then the one with less arguments < other one
 *       Example: f(x) < f(y), g(x) < g(x,y)
 *       Note: Kind is the type of operator
 *
 * Different type expressions:
 * - u value expr, v other: u < v always [DONE]
 * - u product (this is a proper product in the book, not a power), v sum, var or func: u < v <=> u < 1*v <=> u_n < v [Done]
 *       Note: This means we compare u with the 1*v product. Though 1*v is unsimplified, the rule applies.
 *       Example: 2*x^0.5 < x [Note that x is a var expression]
 * - u^p, v sum, var or func: u^p < v <=> u^p < v^1 (I think this is the same as the previous one) [Done]
 *       This means that u^p < v <=> u < v and if they are equal, if p < 1
 * - u sum, v var or func: u < v <=> u < 0+v [Done]
 * - u, v and none of the rules apply: u < v <=> ! v < u [Done]
 *    Example: is x < x^2 ? x is var and x^2 product, so none applies, then
 *    we try to answer if x^2 < x <=> x^2 < x^1 <=> 2 < 1 <=> False, so x < x^2 is True
 */

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

   exprhdlr1 = SCIPgetConsExprExprHdlr(expr1);
   exprhdlr2 = SCIPgetConsExprExprHdlr(expr2);

   if( exprhdlr1 == exprhdlr2 )
   { /* expressions are of the same kind/type */
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "val") == 0 )
      {
         SCIP_Real val1;
         SCIP_Real val2;

         val1 = SCIPgetConsExprExprValueValue(expr1);
         val2 = SCIPgetConsExprExprValueValue(expr2);

         return val1 < val2 ? -1 : val1 == val2 ? 0 : 1;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "var") == 0 )
      {
         int index1;
         int index2;

         index1 = SCIPvarGetIndex(SCIPgetConsExprExprVarVar(expr1));
         index2 = SCIPvarGetIndex(SCIPgetConsExprExprVarVar(expr2));

         return index1 < index2 ? -1 : index1 == index2 ? 0 : 1;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "sum") == 0 )
      {
         int nchildren1;
         int nchildren2;
         int compareresult;
         int i;
         int j;
         SCIP_Real* coefs1;
         SCIP_Real* coefs2;
         SCIP_Real  const1;
         SCIP_Real  const2;

         nchildren1 = SCIPgetConsExprExprNChildren(expr1);
         nchildren2 = SCIPgetConsExprExprNChildren(expr2);
         coefs1 = SCIPgetConsExprExprSumCoefs(expr1);
         coefs2 = SCIPgetConsExprExprSumCoefs(expr2);

         for( i = nchildren1 - 1, j = nchildren2 -1; i >= 0 && j >= 0; --i, --j )
         {
            compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[i], SCIPgetConsExprExprChildren(expr2)[j]);
            if( compareresult != 0 )
               return compareresult;
            else
            {
               /* expressions are equal, compare coefficient */
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

         /* everything is equal, use constant as tie-breaker */
         assert(i == -1 && j == -1);
         const1 = SCIPgetConsExprExprSumConstant(expr1);
         const2 = SCIPgetConsExprExprSumConstant(expr2);
         if( const1 < const2 )
            return -1;
         if( const1 > const2 )
            return 1;

         /* they are equal */
         return 0;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "prod") == 0 )
      {
         int nchildren1;
         int nchildren2;
         int compareresult;
         int i;
         int j;
         SCIP_Real* exponents1;
         SCIP_Real* exponents2;
         SCIP_Real  coef1;
         SCIP_Real  coef2;

         nchildren1 = SCIPgetConsExprExprNChildren(expr1);
         nchildren2 = SCIPgetConsExprExprNChildren(expr2);
         exponents1 = SCIPgetConsExprExprProductExponents(expr1);
         exponents2 = SCIPgetConsExprExprProductExponents(expr2);

         for( i = nchildren1 - 1, j = nchildren2 -1; i >= 0 && j >= 0; --i, --j )
         {
            compareresult = SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[i], SCIPgetConsExprExprChildren(expr2)[j]);
            if( compareresult != 0 )
               return compareresult;
            else
            {
               /* expressions are equal, compare exponents */
               if( exponents1[i] < exponents2[j] )
                  return -1;
               if( exponents1[i] > exponents2[j] )
                  return 1;

               /* exponents are equal, continue */
            }
         }

         /* all children of one expression are children of the other expression, use number of children as a tie-breaker */
         if( i < j )
         {
            assert(i == -1);
            return -1;
         }
         if( i > j )
         {
            assert(j == -1);
            return 1;
         }

         /* everything is equal, use coefficient as tie-breaker */
         assert(i == -1 && j == -1);
         coef1 = SCIPgetConsExprExprProductCoef(expr1);
         coef2 = SCIPgetConsExprExprProductCoef(expr2);
         if( coef1 < coef2 )
            return -1;
         if( coef1 > coef2 )
            return 1;

         /* they are equal */
         return 0;
      }
      /* TODO: this should be a callback (as well as the others), for testing purposes we handle it here
       * only for exp and log */
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "abs") == 0 )
      {
         SCIPerrorMessage("Not implemented yet\n");
         SCIPABORT();
      }
      else
      {
         return SCIPcompareConsExprExprs(SCIPgetConsExprExprChildren(expr1)[0], SCIPgetConsExprExprChildren(expr2)[0]);
      }
   }
   else
   { /* expressions are of different kind/type */
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "val") == 0 )
      {
         return -1;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "val") == 0 )
         return -SCIPcompareConsExprExprs(expr2, expr1);

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

         /* largest expression of product is larger or equal than expr2 => expr1 >= expr2 */
         return 1;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "prod") == 0 )
         return -SCIPcompareConsExprExprs(expr2, expr1);

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

         /* largest expression of sum is larger or equal than expr2 => expr1 >= expr2 */
         return 1;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "sum") == 0 )
         return -SCIPcompareConsExprExprs(expr2, expr1);

      /* at this point we know type(expr1) != type(expr2) and neither is value, product nor sum;
       * if type(expr2) is var, then exp1 is some function */
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "var") == 0 )
         return 1;
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "var") == 0 )
         return -1;

      /* both are a different function */
      return strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), SCIPgetConsExprExprHdlrName(exprhdlr2));
   }

   /* should not get here */
   SCIPerrorMessage("Unexpected behavior in comparison\n");
   assert(0);
   return -9;
}

/**@} */  /* end of simplifying methods */

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
   if( activity == SCIP_INVALID )
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

   debugParse("parsing base from %s\n", expr);

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
         debugParse("Variable %s has been parsed, capturing its expression\n", SCIPvarGetName(var));
         *basetree = (SCIP_CONSEXPR_EXPR*)SCIPhashmapGetImage(vartoexprvarmap, (void *)var);
         SCIPcaptureConsExprExpr(*basetree);
      }
      else
      {
         debugParse("First time parsing variable %s, creating varexpr and adding it to hashmap\n", SCIPvarGetName(var));
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
      debugParse("Done parsing expression, continue with <%s>\n", expr);
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
      debugParse("Parsed value %g, creating a value-expression.\n", value);
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

   debugParse("parsing factor from %s\n", expr);

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

      debugParse("parsed the exponent %g\n", *exponent);
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

   debugParse("parsing term from %s\n", expr);

   /* parse Factor */
   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   SCIP_CALL( parseFactor(scip, conshdlr, vartoexprvarmap, expr, newpos, &exponent, &factortree) );
   expr = *newpos;

   debugParse("back to parsing Term (we have a Factor with exponent %g), continue parsing from %s\n", exponent, expr);

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

         debugParse("while parsing term, read char %c\n", *expr);

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

   debugParse("parsing expression %s\n", expr);

   /* ignore whitespace */
   while( isspace((unsigned char)*expr) )
      ++expr;

   /* if '+' or '-', store it */
   sign = 1.0;
   if( *expr == '+' || *expr == '-' )
   {
      debugParse("while parsing expression, read char %c\n", *expr);
      sign = *expr == '+' ? 1.0 : -1.0;
      ++expr;
   }

   SCIP_CALL( parseTerm(scip, conshdlr, vartoexprvarmap, expr, newpos, &termtree) );
   expr = *newpos;

   debugParse("back to parsing expression (we have the following term), continue parsing from %s\n", expr);

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

         debugParse("while parsing expression, read coefficient %g\n", coef);

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
#if 0
static
SCIP_DECL_CONSINIT(consInitExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitExpr NULL
#endif


/** deinitialization method of constraint handler (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_CONSEXIT(consExitExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitExpr NULL
#endif


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
#if 0
static
SCIP_DECL_CONSEXITPRE(consExitpreExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitpreExpr NULL
#endif


/** solving process initialization method of constraint handler (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_CONSINITSOL(consInitsolExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitsolExpr NULL
#endif


/** solving process deinitialization method of constraint handler (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_CONSEXITSOL(consExitsolExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consExitsolExpr NULL
#endif


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteExpr)
{  /*lint --e{715}*/
   assert(consdata != NULL);
   assert(*consdata != NULL);
   assert((*consdata)->expr != NULL);

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(*consdata)->expr) );

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

   /* get a copy of sourceexpr with transformed vars */
   SCIP_CALL( SCIPwalkConsExprExprDF(scip, sourceexpr, copyExpr, NULL, copyExpr, NULL, &copydata) );
   targetexpr = (SCIP_CONSEXPR_EXPR*)sourceexpr->walkio.ptrval;

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
#if 0
static
SCIP_DECL_CONSINITLP(consInitlpExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consInitlpExpr NULL
#endif


/** separation method of constraint handler for LP solutions */
#if 1
static
SCIP_DECL_CONSSEPALP(consSepalpExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepalpExpr NULL
#endif


/** separation method of constraint handler for arbitrary primal solutions */
#if 0
static
SCIP_DECL_CONSSEPASOL(consSepasolExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consSepasolExpr NULL
#endif


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckExpr)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA*     consdata;
   int soltag;
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
      SCIP_CALL( computeViolation(scip, conss[c], sol, soltag) );

      consdata = SCIPconsGetData(conss[c]);
      assert(consdata != NULL);

      if( SCIPisFeasGT(scip, consdata->lhsviol, 0.0) || SCIPisFeasGT(scip, consdata->rhsviol, 0.0) )
      {
         *result = SCIP_INFEASIBLE;

         /* print reason for infeasibility */
         if( printreason )
         {
            SCIP_CALL( SCIPprintCons(scip, conss[c], NULL) );
            SCIPinfoMessage(scip, NULL, ";\n");

            if( SCIPisFeasGE(scip, consdata->lhsviol, 0.0) )
            {
               SCIPinfoMessage(scip, NULL, "violation: left hand side is violated by %.15g\n", consdata->lhsviol);
            }
            if( SCIPisFeasGE(scip, consdata->rhsviol, 0.0) )
            {
               SCIPinfoMessage(scip, NULL, "violation: right hand side is violated by %.15g\n", consdata->rhsviol);
            }
         }
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
#if 1
static
SCIP_DECL_CONSPROP(consPropExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPropExpr NULL
#endif


/** presolving method of constraint handler */
#if 0
static
SCIP_DECL_CONSPRESOL(consPresolExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consPresolExpr NULL
#endif


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
   lockdata.nlockspos = nlockspos;
   lockdata.nlocksneg = nlocksneg;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, lockVar, NULL, NULL, NULL, &lockdata) );

   return SCIP_OKAY;
}


/** constraint activation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSACTIVE(consActiveExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consActiveExpr NULL
#endif


/** constraint deactivation notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDEACTIVE(consDeactiveExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDeactiveExpr NULL
#endif


/** constraint enabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSENABLE(consEnableExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consEnableExpr NULL
#endif


/** constraint disabling notification method of constraint handler */
#if 0
static
SCIP_DECL_CONSDISABLE(consDisableExpr)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of expr constraint handler not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define consDisableExpr NULL
#endif

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
   SCIP_CALL( SCIPwalkConsExprExprDF(sourcescip, sourceexpr, copyExpr, NULL, copyExpr, NULL, &copydata) );
   targetexpr = (SCIP_CONSEXPR_EXPR*)sourceexpr->walkio.ptrval;

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

   debugParse("str should start at beginning of expr: %s\n", str);

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

   debugParse("created expression constraint: <%s>\n", SCIPconsGetName(*cons));

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
{
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
{
   assert(exprhdlr != NULL);

   exprhdlr->copydata = copydata;
   exprhdlr->freedata = freedata;

   return SCIP_OKAY;
}

/** set the print callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrPrint(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_EXPRPRINT((*print))    /**< print callback (can be NULL) */
)
{
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
{
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
{
   assert(exprhdlr != NULL);

   exprhdlr->inteval = inteval;

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
{
   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(child1 != NULL);
   assert(child2 != NULL);

   if( child1 != NULL && child2 != NULL )
   {
      SCIP_CONSEXPR_EXPR* pair[2] = {child1, child2};

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

         exponent = SCIPexprgraphGetNodeOperatorReal(node);

         assert(nchildren == 1);
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, 1, children, &exponent, 1.0) );

         break;
      }

      case SCIP_EXPR_INTPOWER:
      {
         SCIP_Real exponent;

         exponent = (SCIP_Real)SCIPexprgraphGetNodeOperatorIndex(node);

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
         SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, expr, nchildren, children, NULL, 0.0) );

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
               SCIP_CALL( SCIPcreateConsExprExprProduct(scip, consexprhdlr, &prod, 1, &children[quadelem.idx1], &two, 1.0) );
            }
            else
            {
               SCIP_CONSEXPR_EXPR* prodchildren[2] = {children[quadelem.idx1], children[quadelem.idx2]};
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
         assert(children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprExp(scip, consexprhdlr, expr, children[0]) );

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
      default:
         goto TERMINATE;
   }


TERMINATE:
   /* release all created children expressions (c-1...0) */
   for( --c; c >= 0; --c )
   {
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &children[c]) );
   }

   SCIPfreeBufferArrayNull(scip, &children);

   return SCIP_OKAY;
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
      SCIP_CALL( SCIPwalkConsExprExprDF(scip, *expr, freeExpr, freeExpr, NULL, freeExpr, NULL) );
      *expr = NULL;

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
 * For variables, the local variable bounds are used as interval.
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
   unsigned int            boxtag            /**< tag that uniquely identifies the current variable domains (with its values), or 0 */
   )
{
   EXPRINTEVAL_DATA propdata;

   /* if value is up-to-date, then nothing to do */
   if( boxtag != 0 && expr->intevaltag == boxtag )
      return SCIP_OKAY;

   propdata.aborted = FALSE;
   propdata.boxtag = boxtag;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, intevalExprVisitChild, NULL, intevalExprLeaveExpr, &propdata) );

   if( propdata.aborted )
   {
      SCIPintervalSetEmpty(&expr->interval);
      expr->intevaltag = boxtag;
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
   SCIP_CONSEXPR_EXPR*          oldparent;
   SCIP_CONSEXPREXPRWALK_IO     oldwalkio;
   int                          oldcurrentchild;

   assert(root != NULL);

   /* remember the current data, child and parent of incoming root, in case we are called from within another walk */
   oldcurrentchild = root->walkcurrentchild;
   oldparent       = root->walkparent;
   oldwalkio       = root->walkio;

   root->walkcurrentchild = 0;
   root->walkparent = NULL;

   /* traverse the tree */
   stage = SCIP_CONSEXPREXPRWALK_ENTEREXPR;
   while( TRUE )
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
                     return SCIP_OKAY;
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
                  return SCIP_OKAY;
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
                     return SCIP_OKAY;
               }
            }
            else
            {
               /* visit next child (if any) */
               ++root->walkcurrentchild;
            }
            /* goto visitng (next) */
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
                     return SCIP_OKAY;
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
               return SCIP_OKAY;

            /* goto visited */
            stage = SCIP_CONSEXPREXPRWALK_VISITEDCHILD;
            break;
         default:
            /* unknown stage */
            SCIPABORT();
      }
   }

   /* recover previous information */
   root->walkcurrentchild = oldcurrentchild;
   root->walkparent       = oldparent;
   root->walkio           = oldwalkio;

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
   /* TODO: (optional) add constraint handler specific parameters with SCIPaddTypeParam() here */

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

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, copyExpr, NULL, copyExpr, NULL, &copydata) );
   *copyexpr = (SCIP_CONSEXPR_EXPR*)expr->walkio.ptrval;

   return SCIP_OKAY;
}

/** simplifies an expression */
SCIP_RETCODE SCIPsimplifyConsExprExpr(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR**    expr              /**< expression to be simplified */
   )
{
   SCIP_CONSEXPR_EXPR* aux;

   assert(scip != NULL);
   assert(*expr != NULL);

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, *expr, NULL, NULL, simplifyExpr, simplifyExpr, NULL) );

   aux = *expr;
   *expr = (SCIP_CONSEXPR_EXPR*)(*expr)->walkio.ptrval;

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &aux) );

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

/** overwrites/replaces a child of an expressions */
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

   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &(expr->children[childidx])) );
   expr->children[childidx] = newchild;

   return SCIP_OKAY;
}
