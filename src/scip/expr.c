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

/**@file   expr.c
 * @brief  functions for algebraic expressions
 * @author Ksenia Bestuzheva
 * @author Benjamin Mueller
 * @author Felipe Serrano
 * @author Stefan Vigerske
 */

#include "scip/scip_expr.h"
#include "scip/expr.h"
#include "scip/struct_expr.h"
#include "nlpi/nlpi_ipopt.h" /* for LAPACK */

/*
 * Data structures
 */

/** variable mapping data passed on during copying expressions when copying SCIP instances */
typedef struct
{
   SCIP_HASHMAP*         varmap;             /**< SCIP_HASHMAP mapping variables of the source SCIP to corresponding variables of the target SCIP */
   SCIP_HASHMAP*         consmap;            /**< SCIP_HASHMAP mapping constraints of the source SCIP to corresponding constraints of the target SCIP */
   SCIP_Bool             global;             /**< should a global or a local copy be created */
   SCIP_Bool             valid;              /**< indicates whether every variable copy was valid */
} COPY_MAPVAR_DATA;

/** printing to file data */
struct SCIP_ExprPrintData
{
   FILE*                 file;               /**< file to print to */
   SCIP_EXPRITER*        iterator;           /**< iterator to use */
   SCIP_Bool             closefile;          /**< whether file need to be closed when finished printing */
   SCIP_HASHMAP*         leaveexprs;         /**< hashmap storing leave (no children) expressions */
   SCIP_EXPRPRINT_WHAT   whattoprint;        /**< flags that indicate what to print for each expression */
};

/*
 * Local methods
 */

/** creates an expression */
static
SCIP_RETCODE createExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
//   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr,               /**< pointer where to store expression */
   SCIP_EXPRHDLR*        exprhdlr,           /**< expression handler */
   SCIP_EXPRDATA*        exprdata,           /**< expression data (expression assumes ownership) */
   int                   nchildren,          /**< number of children */
   SCIP_EXPR**           children,           /**< children (can be NULL if nchildren is 0) */
   SCIP_DECL_EXPR_OWNERDATACREATE((*ownerdatacreate)), /**< function to call to create ownerdata */
   SCIP_EXPR_OWNERDATACREATEDATA* ownerdatacreatedata, /**< data to pass to ownerdatacreate */
   SCIP_DECL_EXPR_OWNERDATAFREE((*ownerdatafree)),     /**< function to call when freeing expression, e.g., to free ownerdata */
   )
{
   int c;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);
   assert(exprhdlr != NULL);
   assert(children != NULL || nchildren == 0);
   assert(exprdata == NULL || exprhdlr->copydata != NULL); /* copydata must be available if there is expression data */
   assert(exprdata == NULL || exprhdlr->freedata != NULL); /* freedata must be available if there is expression data */

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, expr) );

   (*expr)->exprhdlr = exprhdlr;
   (*expr)->exprdata = exprdata;
   (*expr)->curvature = SCIP_EXPRCURV_UNKNOWN;

   /* initialize activity to entire interval */
   SCIPintervalSetEntire(SCIP_INTERVAL_INFINITY, &(*expr)->activity);

   if( nchildren > 0 )
   {
      SCIP_ALLOC( BMSduplicateBlockMemoryArray(blkmem, &(*expr)->children, children, nchildren) );
      (*expr)->nchildren = nchildren;
      (*expr)->childrensize = nchildren;

      for( c = 0; c < nchildren; ++c )
         SCIPcaptureExpr((*expr)->children[c]);
   }

   SCIPcaptureExpr(*expr);

   if( ownerdatacreate != NULL )
   {
      SCIP_CALL( ownerdatacreate(set->scip, *expr, &(*expr)->ownerdata, ownerdatacreatedata) );
   }
   (*expr)->ownerdatafree = ownerdatafree;

   return SCIP_OKAY;
}

/** frees an expression */
static
SCIP_RETCODE freeExpr(
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR**           expr                /**< pointer to free the expression */
   )
{
   assert(scip != NULL);
   assert(expr != NULL);
   assert(*expr != NULL);
   assert((*expr)->nuses == 1);

   /* call ownerdatafree callback, if given
    * we intentially call this also if ownerdata is NULL, so owner can be notified without storing data
    */
   if( (*expr)->ownerdatafree != NULL )
   {
      SCIP_CALL( (*expr)->ownerdatafree(set->scip, *expr, &(*expr)->ownerdata) );
   }
   assert((*expr)->ownerdata == NULL);

   /* free children array, if any */
   BMSfreeBlockMemoryArrayNull(blkmem, &(*expr)->children, (*expr)->childrensize);

   BMSfreeBlockMemory(blkmem, expr);
   assert(*expr == NULL);

   return SCIP_OKAY;
}

static
SCIP_RETCODE evalAndDiff(
   SCIP_SET*             set,                /**< global SCIP settings */
//   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPR*            expr,               /**< expression to be evaluated */
   SCIP_SOL*             sol,                /**< solution to be evaluated */
   SCIP_Longint          soltag              /**< tag that uniquely identifies the solution (with its values), or 0. */
   )
{
   SCIP_EXPRITER* it;

   assert(set != NULL);
   assert(blkmem != NULL);
   assert(expr != NULL);

   /* assume we'll get a domain error, so we don't have to get this expr back if we abort the iteration
    * if there is no domain error, then we will overwrite the evalvalue in the last leaveexpr stage
    */
   expr->evalvalue = SCIP_INVALID;
   expr->evaltag = soltag;
   expr->dot = SCIP_INVALID;

   SCIP_CALL( SCIPexpriteratorCreate(set, &it, blkmem) );
   SCIP_CALL( SCIPexpriteratorInit(it, expr, SCIP_EXPRITER_DFS, TRUE) );
   SCIPexpriteratorSetStagesDFS(it, SCIP_EXPRITER_LEAVEEXPR);

   for( expr = SCIPexpriteratorGetCurrent(it); !SCIPexpriteratorIsEnd(it); expr = SCIPexpriteratorGetNext(it) )  /*lint !e441*/
   {
      /* evaluate expression only if necessary */
      if( soltag == 0 || expr->evaltag != soltag )
      {
         SCIP_CALL( SCIPcallExprhdlrEval(set->scip, expr, &expr->evalvalue, NULL, sol) );

         expr->evaltag = soltag;
      }

      if( expr->evalvalue == SCIP_INVALID ) /*lint !e777*/
         break;

      /* compute forward diff */
      SCIP_CALL( SCIPcallExprhdlrFwdiff(set->scip, expr, &expr->dot) );

      if( expr->dot == SCIP_INVALID ) /*lint !e777*/
         break;
   }

   SCIPexpriteratorFree(&it);

   return SCIP_OKAY;
}
