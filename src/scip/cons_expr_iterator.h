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
 *
 * If type is DFS, then stopstages will be set to ENTEREXPR. Use SCIPexpriteratorSetStagesDFS to change this.
 */
EXTERN
SCIP_RETCODE SCIPexpriteratorInit(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*         expr,        /**< expression of the iterator */
   SCIP_CONSEXPRITERATOR_TYPE  type,        /**< type of expression iterator */
   SCIP_Bool                   allowrevisit /**< whether expression are allowed to be visited more than once */
   );

/** specifies in which stages to stop a DFS iterator
 *
 * @param stopstages should be a bitwise OR of different SCIP_CONSEXPREXPRWALK_STAGE values
 *
 * If the current stage is not one of the stopstages, then the iterator will be moved on.
 */
EXTERN
void SCIPexpriteratorSetStagesDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   unsigned int                stopstages   /**< the stages in which to stop when iterating via DFS */
   );

/** gets the current expression that the expression iterator points to */
EXTERN
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetCurrent(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gets the current stage that the expression iterator is in when using DFS
 *
 * If the iterator has finished (IsEnd() is TRUE), then the stage is undefined.
 */
EXTERN
SCIP_CONSEXPREXPRWALK_STAGE SCIPexpriteratorGetStageDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gets the child index that the expression iterator considers when in DFS mode and stage visitingchild or visitedchild */
EXTERN
int SCIPexpriteratorGetChildIdxDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gets the child expression that the expression iterator considers when in DFS mode and stage visitingchild or visitedchild */
EXTERN
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetChildExprDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gives the parent of the current expression of an expression iteration if in DFS mode
 *
 * @return the expression from which the current expression has been accessed
 */
EXTERN
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetParentDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gives the iterator specific user data of the current expression
 *
 * \note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
EXTERN
SCIP_CONSEXPREXPRWALK_IO SCIPexpriteratorGetUserData(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** sets the iterator specific user data of the current expression for an expression iteration if in DFS mode
 *
 * \note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
EXTERN
void SCIPexpriteratorSetUserData(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   SCIP_CONSEXPREXPRWALK_IO    userdata     /**< data to be stored */
   );

/** moves the iterator to the next expression according to the mode of the expression iterator
 *
 * @return the next expression, if any, and NULL otherwise
 */
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
