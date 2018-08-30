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

/**@file   cons_expr_iterator.h
 * @brief  expression tree iterators
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_EXPR_ITERATOR_H__
#define __SCIP_CONS_EXPR_ITERATOR_H__


#include "scip/scip.h"
#include "scip/type_cons_expr.h"


#ifdef __cplusplus
extern "C" {
#endif


/**@name Expression Iterator Methods */
/**@{ */

/** creates an expression iterator */
EXTERN
SCIP_RETCODE SCIPexpriteratorCreate(
   SCIP_CONSEXPR_ITERATOR**    iterator,    /**< buffer to store expression iterator */
   SCIP_CONSHDLR*              consexprhdlr,/**< expr constraint handler, might be NULL */
   BMS_BLKMEM*                 blkmem       /**< block memory used to store hash map entries */
   );

/** frees an expression iterator */
EXTERN
void SCIPexpriteratorFree(
   SCIP_CONSEXPR_ITERATOR**    iterator     /**< pointer to the expression iterator */
   );

/** initializes an expression iterator
 *
 * \note If no conshdlr has been given when creating the iterator, then allowrevisit must be TRUE and type must not be DFS.
 */
EXTERN
SCIP_RETCODE SCIPexpriteratorInit(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*         expr,        /**< expression of the iterator */
   SCIP_CONSEXPRITERATOR_TYPE  type,        /**< type of expression iterator */
   SCIP_Bool                   allowrevisit /**< whether expression are allowed to be visited more than once */
   );

/** gets the current expression that the expression iterator points to */
EXTERN
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetCurrent(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gets the next expression according to the mode of the expression iterator */
EXTERN
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetNext(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** returns whether the iterator visited all expressions already */
EXTERN
SCIP_Bool SCIPexpriteratorIsEnd(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
