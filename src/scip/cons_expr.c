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
#include "scip/cons_expr_log.h"
#include "scip/cons_expr_abs.h"

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

   SCIP_CONSEXPR_EXPR*   expr;               /**< expression that represents this constraint (must evaluate to 0 (FALSE) or 1 (TRUE)) */
   SCIP_Real             lhs;                /**< left-hand side */
   SCIP_Real             rhs;                /**< right-hand side */

   SCIP_Real             lhsviol;            /**< violation of left-hand side by current solution (used temporarily inside constraint handler) */
   SCIP_Real             rhsviol;            /**< violation of right-hand side by current solution (used temporarily inside constraint handler) */

   unsigned int          ispropagated:1;     /**< did we propagate the current bounds already? */

   SCIP_NLROW*           nlrow;              /**< a nonlinear row representation of this constraint */
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

   SCIP_EVENTHDLR*          eventhdlr;       /**< handler for variable bound change events */

   int                      maxproprounds;   /**< limit on number of propagation rounds for a set of constraints within one round of SCIP propagation */
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
} COPY_DATA;

/** variable mapping data passed on during copying expressions when copying SCIP instances */
typedef struct
{
   SCIP_HASHMAP*           varmap;           /**< SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP */
   SCIP_HASHMAP*           consmap;          /**< SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP */
   SCIP_Bool               global;           /**< should a global or a local copy be created */
   SCIP_Bool               valid;            /**< indicates whether every variable copy was valid */
} COPY_MAPVAR_DATA;

/** data passed on during conversion to classic expression */
typedef struct
{
   SCIP_CONSEXPR_EXPR**    varexprs;         /**< variable expressions */
   int                     nvarexprs;        /**< number of variable expressions */
} MAKECLASSICEXPR_DATA;

struct SCIP_ConsExpr_PrintDotData
{
   FILE*                   file;             /**< file to print to */
   SCIP_Bool               closefile;        /**< whether file need to be closed when finished printing */
   SCIP_HASHMAP*           visitedexprs;     /**< hashmap storing expressions that have been printed already */
   SCIP_CONSEXPR_PRINTDOT_WHAT whattoprint;  /**< flags that indicate what to print for each expression */
};

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
{
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
         assert(SCIPcompareExprs(expr, *newexpr) == 0);
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
   while( TRUE );

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

/** expression walk callback to print an expression in dot format */
static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(printExprDot)
{
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
   SCIP_CALL( (*expr->exprhdlr->inteval)(scip, expr, &interval) );

   /* stop if callback returned an empty interval */
   if( SCIPintervalIsEmpty(SCIPinfinity(scip), interval) )
   {
      SCIPintervalSetEmpty(&expr->interval);
      propdata->aborted = TRUE;
      *result = SCIP_CONSEXPREXPRWALK_ABORT;
   }
   else
   {
      /* intersect new interval with the previous one */
      if( propdata->intersect )
         SCIPintervalIntersect(&expr->interval, expr->interval, interval);
      else
         SCIPintervalSetBounds(&expr->interval, interval.inf, interval.sup);

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
         hashkey += (unsigned int) expr->exprhdlr->name[i];

      hashkey = SCIPcalcFibHash(hashkey);
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
      assert(SCIPcompareExprs(child, newchild) == 0);

      /** @todo use SCIPsetConsExprExprChild() (simplify-branch) to replace child */
      SCIP_CALL( SCIPreleaseConsExprExpr(scip, &child) );

      expr->children[expr->walkcurrentchild] = newchild;
      SCIPcaptureConsExprExpr(newchild);

      SCIPdebugMessage("replace common inner node expression\n");

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

   /* add variable expression if not seen so far; there is only one variable expression representing a variable */
   if( strcmp(expr->exprhdlr->name, "var") == 0 && !SCIPhashmapExists(getvarsdata->varexprsmap, (void*) expr) )
   {
      assert(SCIPgetNVars(scip) >= getvarsdata->nvarexprs + 1);

      getvarsdata->varexprs[ getvarsdata->nvarexprs ] = expr;
      ++(getvarsdata->nvarexprs);
      SCIP_CALL( SCIPhashmapInsert(getvarsdata->varexprsmap, (void*) expr, NULL) );
   }

   return SCIP_OKAY;
}

/**@} */  /* end of walking methods */

/* export this function here, so it can be used by unittests but is not really part of the API */
/** propagates bounds for each sub-expression in the constraint by using variable bounds; the resulting bounds for the
 *  root expression will be intersected with the [lhs,rhs] which might lead to an empty interval
 */
SCIP_RETCODE forwardPropCons(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONS*              cons,             /**< constraint to propagate */
   SCIP_Bool               intersect,        /**< should the new expr. bounds be intersected with the previous ones? */
   SCIP_Bool*              infeasible        /**< buffer to store whether an expression's bounds were propagated to an empty interval */
   );
SCIP_RETCODE forwardPropCons(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONS*              cons,             /**< constraint to propagate */
   SCIP_Bool               intersect,        /**< should the new expr. bounds be intersected with the previous ones? */
   SCIP_Bool*              infeasible        /**< buffer to store whether an expression's bounds were propagated to an empty interval */
   )
{
   SCIP_INTERVAL interval;
   SCIP_CONSDATA* consdata;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(infeasible != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   *infeasible = FALSE;

   /* propagate active constraints only */
   if( !SCIPconsIsActive(cons) && SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
      return SCIP_OKAY;

   /* use 0 tag to recompute intervals */
   SCIP_CALL( SCIPevalConsExprExprInterval(scip, consdata->expr, intersect, 0) );

   /* compare root expression interval with constraint sides; store the result in the root expression */
   SCIPintervalSetBounds(&interval, consdata->lhs, consdata->rhs);
   SCIPtightenConsExprExprInterval(scip, consdata->expr, interval, infeasible, NULL);

#ifdef SCIP_DEBUG
   if( *infeasible )
   {
      SCIPdebugMessage(" -> found empty bound for an expression during forward propagation of constraint %s\n",
         SCIPconsGetName(cons));
   }
#endif

   /* TODO if the root expression interval could not be tightened by constraint sides,
    * then the constraint is redundant and should be deleted (locally)
    */

   return SCIP_OKAY;
}

/* export this function here, so it can be used by unittests but is not really part of the API */
/** propagates bounds for each sub-expression of a given set of constraints by starting from the root expressions; the
 *  expression will be traversed in breadth first search by using a queue
 *
 *  @note calling this function requires feasible intervals for each sub-expression; this is guaranteed by calling
 *  forwardPropCons() before calling this function
 */
SCIP_RETCODE reversePropConss(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONS**             conss,            /**< constraints to propagate */
   int                     nconss,           /**< total number of constraints to propagate */
   SCIP_Bool*              infeasible,       /**< buffer to store whether an expression's bounds were propagated to an empty interval */
   int*                    ntightenings      /**< buffer to store the number of (variable) tightenings */
   );
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
   SCIP_CALL( SCIPqueueCreate(&queue, SCIPgetNVars(scip), 2) );

   /* add root expressions to the queue */
   for( i = 0; i < nconss; ++i )
   {
      assert(conss[i] != NULL);
      consdata = SCIPconsGetData(conss[i]);
      assert(consdata != NULL);

      /* propagate active constraints only */
      if( !SCIPconsIsActive(conss[i]) && SCIPgetStage(scip) >= SCIP_STAGE_TRANSFORMED )
         return SCIP_OKAY;

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
SCIP_RETCODE propConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  nchgbds             /**< buffer to add the number of changed bounds */
   );
SCIP_RETCODE propConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to propagate */
   int                   nconss,             /**< total number of constraints */
   SCIP_RESULT*          result,             /**< pointer to store the result */
   int*                  nchgbds             /**< buffer to add the number of changed bounds */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSDATA* consdata;
   SCIP_Bool cutoff;
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
   ntightenings = 0;
   cutoff = FALSE;

   /* main propagation loop */
   do
   {
      SCIPdebugMessage("start propagation round %d\n", roundnr);
      success = FALSE;

      /* apply forward propagation; recompute expression intervals if it is called for the first time (this also marks
       * all expressions as non-tightened)
       */
      for( i = 0; i < nconss; ++i )
      {
         consdata = SCIPconsGetData(conss[i]);
         assert(consdata != NULL);

         if( SCIPconsIsActive(conss[i]) && !consdata->ispropagated )
         {
            SCIP_CALL( forwardPropCons(scip, conss[i], (roundnr != 0), &cutoff) );

            if( cutoff )
            {
               SCIPdebugMessage(" -> cutoff\n");
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }

            /* mark constraint as propagated; this will be reset via the event system when we find a variable tightening */
            consdata->ispropagated = TRUE;
         }
      }

      /* apply backward propagation; mark constraint as propagated */
      SCIP_CALL( reversePropConss(scip, conss, nconss, &cutoff, &ntightenings) );

      /* @todo add parameter for the minimum number of tightenings to trigger a new propagation round */
      success = ntightenings > 0;

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
   SCIP_CALL( SCIPhashmapCreate(&getvarsdata.varexprsmap, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip))) );

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
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &consdata->varexprs, SCIPgetNVars(scip)) );

   SCIP_CALL( getVarExprs(scip, consdata->expr, consdata->varexprs, &(consdata->nvarexprs)) );
   assert(SCIPgetNVars(scip) >= consdata->nvarexprs);

   /* realloc array if there are less variable expression than variables */
   if( SCIPgetNVars(scip) > consdata->nvarexprs )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &consdata->varexprs, SCIPgetNVars(scip), consdata->nvarexprs) );
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
   assert(consdata != NULL);

   /* skip if we have stored the variable expressions already*/
   if( consdata->varexprs == NULL )
      return SCIP_OKAY;

   assert(consdata->varexprs != NULL);
   assert(consdata->nvarexprs >= 0);

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

      SCIP_CALL( SCIPallocBlockMemory(scip, &(consdata->vareventdata[i])) );
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

      SCIPfreeBlockMemory(scip, &consdata->vareventdata[i]);
      consdata->vareventdata[i] = NULL;
   }

   SCIPfreeBlockMemoryArray(scip, &consdata->vareventdata, consdata->nvarexprs);
   consdata->vareventdata = NULL;

   return SCIP_OKAY;
}

/** processes variable fixing or bound change event */
static
SCIP_DECL_EVENTEXEC(processVarEvent)
{
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

   SCIPdebugMessage("  exec event %d for %s in %s\n", eventtype, SCIPvarGetName(var), SCIPconsGetName(cons));

   /* mark constraint to be propagated again */
   if( (eventtype & SCIP_EVENTTYPE_BOUNDTIGHTENED) != 0 )
   {
      SCIPdebugMessage("  propagate %s again\n", SCIPconsGetName(cons));
      consdata->ispropagated = FALSE;
   }

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

   return expr1 == expr2 || SCIPcompareExprs(expr1, expr2) == 0;
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
 *  (with the same hash) are structurally the same we use the function SCIPcompareExprs()
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
         assert(SCIPcompareExprs(consdata->expr, newroot) == 0);

         SCIP_CALL( SCIPreleaseConsExprExpr(scip, &consdata->expr) );

         consdata->expr = newroot;
         SCIPcaptureConsExprExpr(newroot);

         SCIPdebugMessage("replace common root expression\n");
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

static
SCIP_DECL_CONSEXPREXPRWALK_VISIT(makeClassicExpr)
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr;
   MAKECLASSICEXPR_DATA* walkdata;
   SCIP_EXPR** children = NULL;
   int nchildren;
   int c;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(data != NULL);
   assert(stage == SCIP_CONSEXPREXPRWALK_LEAVEEXPR);

   exprhdlr = SCIPgetConsExprExprHdlr(expr);
   walkdata = (MAKECLASSICEXPR_DATA*) data;

   nchildren = SCIPgetConsExprExprNChildren(expr);

   /* collect children expressions from children, if any */
   if( nchildren > 0 )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &children, nchildren) );
      for( c = 0; c < nchildren; ++c )
      {
         children[c] = (SCIP_EXPR*)SCIPgetConsExprExprChildren(expr)[c]->walkio.ptrval;
         assert(children[c] != NULL);
      }
   }

   /* create expression and store in walkio */
   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "var") == 0 )
   {
      int varidx;

      /* find variable expression in varexprs array
       * the position in the array determines the index of the variable in the classic expression
       * TODO if varexprs are sorted, then can do this more efficient
       */
      for( varidx = 0; varidx < walkdata->nvarexprs; ++varidx )
         if( walkdata->varexprs[varidx] == expr )
            break;
      assert(varidx < walkdata->nvarexprs);

      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), (SCIP_EXPR**)&expr->walkio.ptrval, SCIP_EXPR_VARIDX, varidx) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "val") == 0 )
   {
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), (SCIP_EXPR**)&expr->walkio.ptrval, SCIP_EXPR_CONST, SCIPgetConsExprExprValueValue(expr)) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "sum") == 0 )
   {
      SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), (SCIP_EXPR**)&expr->walkio.ptrval, nchildren, children, SCIPgetConsExprExprSumCoefs(expr), SCIPgetConsExprExprSumConstant(expr)) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "prod") == 0 )
   {
      SCIP_EXPRDATA_MONOMIAL* monomial;
      SCIP_CALL( SCIPexprCreateMonomial(SCIPblkmem(scip), &monomial, SCIPgetConsExprExprProductCoef(expr), nchildren, NULL, SCIPgetConsExprExprProductExponents(expr)) );
      SCIP_CALL( SCIPexprCreatePolynomial(SCIPblkmem(scip), (SCIP_EXPR**)&expr->walkio.ptrval, nchildren, children, 1, &monomial, 0.0, FALSE) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "abs") == 0 )
   {
      assert(nchildren == 1);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), (SCIP_EXPR**)&expr->walkio.ptrval, SCIP_EXPR_ABS, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "exp") == 0 )
   {
      assert(nchildren == 1);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), (SCIP_EXPR**)&expr->walkio.ptrval, SCIP_EXPR_EXP, children[0]) );
   }
   else if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr), "log") == 0 )
   {
      assert(nchildren == 1);
      SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), (SCIP_EXPR**)&expr->walkio.ptrval, SCIP_EXPR_LOG, children[0]) );
   }
   else
   {
      SCIPerrorMessage("unsupported expression handler <%s>, cannot convert to classical expression\n", SCIPgetConsExprExprHdlrName(exprhdlr));
      return SCIP_ERROR;
   }

   SCIPfreeBufferArrayNull(scip, &children);

   *result = SCIP_CONSEXPREXPRWALK_CONTINUE;

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
   MAKECLASSICEXPR_DATA walkdata = { .varexprs = varexprs, .nvarexprs = nvarexprs };
   SCIP_VAR** vars;
   int i;

   assert(scip != NULL);
   assert(expr != NULL);
   assert(varexprs != NULL);  /* we could also create this here, if NULL; but for now, assume it is given by called */
   assert(exprtree != NULL);

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvarexprs) );
   for( i = 0; i < nvarexprs; ++i )
      vars[i] = SCIPgetConsExprExprVarVar(varexprs[i]);

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, NULL, NULL, makeClassicExpr, &walkdata) );
   assert(expr->walkio.ptrval != NULL);

   SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), exprtree, (SCIP_EXPR*)expr->walkio.ptrval, nvarexprs, 0, NULL) );
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
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_EXPRTREE* exprtree;

      conshdlrdata = SCIPconshdlrGetData(SCIPconsGetHdlr(cons));
      assert(conshdlrdata != NULL);

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
   SCIP_CONSDATA* consdata;
   int c;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);  /*lint !e613*/
      assert(consdata != NULL);

      /* add nlrow respresentation to NLP, if NLP had been constructed */
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
   SCIP_CONSDATA* consdata;
   int c;

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
static
SCIP_DECL_CONSSEPALP(consSepalpExpr)
{  /*lint --e{715}*/

   /* TODO separate */
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}


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
   SCIP_CONSDATA* consdata;
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

   /* TODO do domain propagation */
   /* TODO try to separate */
   /* TODO register branching candidates */

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsExpr)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
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

   /* TODO do domain propagation */
   /* TODO ask for solving LP */

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
static
SCIP_DECL_CONSPROP(consPropExpr)
{  /*lint --e{715}*/
   int nchgbds;

   nchgbds = 0;

   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, result, &nchgbds) );
   assert(nchgbds >= 0);

   return SCIP_OKAY;
}


/** presolving method of constraint handler */
static
SCIP_DECL_CONSPRESOL(consPresolExpr)
{  /*lint --e{715}*/
   SCIP_CALL( replaceCommonSubexpressions(scip, conss, nconss) );

   SCIP_CALL( propConss(scip, conshdlr, conss, nconss, result, nchgbds) );
   assert(*nchgbds >= 0);

   if( *nchgbds > 0 )
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
   lockdata.nlockspos = nlockspos;
   lockdata.nlocksneg = nlocksneg;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, consdata->expr, lockVar, NULL, NULL, NULL, &lockdata) );

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

/** set the reverse propagation callback of an expression handler */
SCIP_RETCODE SCIPsetConsExprExprHdlrReverseProp(
   SCIP*                      scip,          /**< SCIP data structure */
   SCIP_CONSHDLR*             conshdlr,      /**< expression constraint handler */
   SCIP_CONSEXPR_EXPRHDLR*    exprhdlr,      /**< expression handler */
   SCIP_DECL_CONSEXPR_REVERSEPROP((*reverseprop))/**< reverse propagation callback (can be NULL) */
   )
{
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
{
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
{
   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(child1 != NULL);
   assert(child2 != NULL);

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
               SCIP_CONSEXPR_EXPR* prodchildren[2];
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
      case SCIP_EXPR_LOG:
      {
         assert(nchildren == 1);
         assert(children[0] != NULL);

         SCIP_CALL( SCIPcreateConsExprExprLog(scip, consexprhdlr, expr, children[0]) );

         break;
      }
      case SCIP_EXPR_ABS:
      {
         assert(nchildren == 1);
         assert(children[0] != NULL);

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
   pclose(f);

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
   SCIP_Bool               intersect,        /**< should the new expr. bounds be intersected with the previous ones? */
   unsigned int            boxtag            /**< tag that uniquely identifies the current variable domains (with its values), or 0 */
   )
{
   EXPRINTEVAL_DATA propdata;

   /* if value is up-to-date, then nothing to do */
   if( boxtag != 0 && expr->intevaltag == boxtag )
      return SCIP_OKAY;

   propdata.aborted = FALSE;
   propdata.boxtag = boxtag;
   propdata.intersect = intersect;

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
   int*                    ntightenings      /**< buffer to add the total number of tightenings (NULL if not needed) */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(cutoff != NULL);
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

      /* apply bound tightening directly for variable expressions */
      if( strcmp(expr->exprhdlr->name, "var") == 0 )
      {
         SCIP_VAR* var;
         SCIP_Bool tightened;

#ifdef SCIP_DEBUG
         SCIP_Real oldlb;
         SCIP_Real oldub;
         oldlb = SCIPvarGetLbLocal(SCIPgetConsExprExprVarVar(expr));
         oldub = SCIPvarGetUbLocal(SCIPgetConsExprExprVarVar(expr));
#endif

         var = SCIPgetConsExprExprVarVar(expr);
         assert(var != NULL);

         /* tighten lower bound */
         SCIP_CALL( SCIPtightenVarLb(scip, var, SCIPintervalGetInf(expr->interval), FALSE, cutoff, &tightened) );

         if( ntightenings != NULL && tightened )
            ++*ntightenings;

         /* tighten upper bound */
         if( !(*cutoff) )
         {
            SCIP_CALL( SCIPtightenVarUb(scip, var, SCIPintervalGetSup(expr->interval), FALSE, cutoff, &tightened) );

            if( ntightenings != NULL && tightened )
               ++*ntightenings;
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
unsigned int SCIPgetConsExprExprHashkey(
   SCIP*                   scip,             /**< SCIP data structure */
   SCIP_CONSEXPR_EXPR*     expr              /**< expression */
   )
{
   SCIP_HASHMAP* expr2key;
   unsigned int hashkey;

   assert(expr != NULL);

   SCIP_CALL( SCIPhashmapCreate(&expr2key, SCIPblkmem(scip), SCIPcalcHashtableSize(SCIPgetNVars(scip))) );

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, NULL, NULL, NULL, hashExprLeaveExpr, (void*)expr2key) );

   assert(SCIPhashmapExists(expr2key, (void*)expr));  /* we just computed the hash, so should be in the map */
   hashkey = (unsigned int)(size_t)SCIPhashmapGetImage(expr2key, (void*)expr);

   SCIPhashmapFree(&expr2key);

   return hashkey;
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
   result = SCIP_CONSEXPREXPRWALK_CONTINUE;
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
   SCIP_CALL( SCIPaddIntParam(scip, "constraints/" CONSHDLR_NAME "/maxproprounds",
         "limit on number of propagation rounds for a set of constraints within one round of SCIP propagation",
         &conshdlrdata->maxproprounds, FALSE, 10, 0, INT_MAX, NULL, NULL) );

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
 * - u, v are functional expressions (exp, log, etc): u < v <=> Kind(u) < Kind(v), or if they are equal, args(u) < args(v)
 *    (the first argument and so on), if all common arguments are equal, then the one with less arguments < other one
 *       Example: f(x) < f(y), g(x) < g(x,y)
 *       Note: Kind is the type of operator
 *
 * Different type expressions:
 * - u value expr, v other: u < v always [DONE]
 * - u product (this is a proper product in the book, not a power), v sum, var or func: u < v <=> u < 1*v <=> u_n < v
 *       Note: This means we compare u with the 1*v product. Though 1*v is unsimplified, the rule applies.
 *       Example: 2*x^0.5 < x [Note that x is a var expression]
 * - u^p, v sum, var or func: u^p < v <=> u^p < v^1 (I think this is the same as the previous one)
 *       This means that u^p < v <=> u < v and if they are equal, if p < 1
 * - u sum, v var or func: u < v <=> u < 0+v
 * - u sum, v var or func: u < v <=> u < 0+v
 * - u, v and none of the rules apply: u < v <=> ! v < u
 *    Example: is x < x^2 ? x is var and x^2 product, so none applies, then
 *    we try to answer if x^2 < x <=> x^2 < x^1 <=> 2 < 1 <=> False, so x < x^2 is True
 */

/** compare expressions
 * The given expressions are assumed to be simplified */
int SCIPcompareExprs(
   SCIP_CONSEXPR_EXPR*   expr1,              /**< first expression */
   SCIP_CONSEXPR_EXPR*   expr2               /**< second expression */
   )
{
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr1;
   SCIP_CONSEXPR_EXPRHDLR* exprhdlr2;

   exprhdlr1 = SCIPgetConsExprExprHdlr(expr1);
   exprhdlr2 = SCIPgetConsExprExprHdlr(expr2);

   if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), SCIPgetConsExprExprHdlrName(exprhdlr2)) == 0 )
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
            compareresult = SCIPcompareExprs(SCIPgetConsExprExprChildren(expr1)[i], SCIPgetConsExprExprChildren(expr2)[j]);
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

         /* all children of one expression are children of the other expression, use amount of children as a tie-breaker */
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
            compareresult = SCIPcompareExprs(SCIPgetConsExprExprChildren(expr1)[i], SCIPgetConsExprExprChildren(expr2)[j]);
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

         /* all children of one expression are children of the other expression, use amount of children as a tie-breaker */
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
      /* TODO: the same function! */
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "abs") == 0 )
      {
         return SCIPcompareExprs(SCIPgetConsExprExprChildren(expr1)[0], SCIPgetConsExprExprChildren(expr2)[0]);
      }

      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "exp") == 0 )
      {
         return SCIPcompareExprs(SCIPgetConsExprExprChildren(expr1)[0], SCIPgetConsExprExprChildren(expr2)[0]);
      }

      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "log") == 0 )
      {
         return SCIPcompareExprs(SCIPgetConsExprExprChildren(expr1)[0], SCIPgetConsExprExprChildren(expr2)[0]);
      }
   }
   else
   { /* expressions are of different kind/type */
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "val") == 0 )
      {
         return -1;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "val") == 0 )
         return -1 * SCIPcompareExprs(expr2, expr1);

      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "prod") == 0 )
      {
         int compareresult;
         int nchildren;

         nchildren = SCIPgetConsExprExprNChildren(expr1);
         compareresult = SCIPcompareExprs(SCIPgetConsExprExprChildren(expr1)[nchildren-1], expr2);

         if( compareresult != 0 )
            return compareresult;

         /* base of the largest expression of the product is equal to expr2, exponent might tell us that expr2 is larger */
         if( SCIPgetConsExprExprProductExponents(expr1)[nchildren-1] < 1.0 )
            return -1;

         /* largest expression of product is larger or equal than expr2 => expr1 >= expr2 */
         return 1;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "prod") == 0 )
         return -1 * SCIPcompareExprs(expr2, expr1);

      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr1), "sum") == 0 )
      {
         int compareresult;
         int nchildren;

         nchildren = SCIPgetConsExprExprNChildren(expr1);
         compareresult = SCIPcompareExprs(SCIPgetConsExprExprChildren(expr1)[nchildren-1], expr2);

         if( compareresult != 0 )
            return compareresult;

         /* "base" of the largest expression of the sum is equal to expr2, coefficient might tell us that expr2 is larger */
         if( SCIPgetConsExprExprSumCoefs(expr1)[nchildren-1] < 1.0 )
            return -1;

         /* largest expression of sum is larger or equal than expr2 => expr1 >= expr2 */
         return 1;
      }
      if( strcmp(SCIPgetConsExprExprHdlrName(exprhdlr2), "sum") == 0 )
         return -1 * SCIPcompareExprs(expr2, expr1);

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
   return SCIP_ERROR;
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
   copydata.mapvardata = NULL;

   SCIP_CALL( SCIPwalkConsExprExprDF(scip, expr, copyExpr, NULL, copyExpr, NULL, &copydata) );
   *copyexpr = (SCIP_CONSEXPR_EXPR*)expr->walkio.ptrval;

   return SCIP_OKAY;
}
