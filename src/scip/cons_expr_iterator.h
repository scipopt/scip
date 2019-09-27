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
SCIP_EXPORT
SCIP_RETCODE SCIPexpriteratorCreate(
   SCIP_CONSEXPR_ITERATOR**    iterator,    /**< buffer to store expression iterator */
   SCIP_CONSHDLR*              consexprhdlr,/**< expr constraint handler, might be NULL */
   BMS_BLKMEM*                 blkmem       /**< block memory used to store hash map entries */
   );

/** frees an expression iterator */
SCIP_EXPORT
void SCIPexpriteratorFree(
   SCIP_CONSEXPR_ITERATOR**    iterator     /**< pointer to the expression iterator */
   );

/** returns whether expression iterator is current initialized */
SCIP_EXPORT
SCIP_Bool SCIPexpriteratorIsInit(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** initializes an expression iterator
 *
 * @note If no conshdlr has been given when creating the iterator, then allowrevisit must be TRUE and type must not be DFS.
 *
 * @note If expr is NULL, then iterator will be ended (SCIPexpriteratorIsEnd() is TRUE). Useful if following with SCIPexpriteratorRestartDFS().
 *
 * If type is DFS, then stopstages will be set to ENTEREXPR. Use SCIPexpriteratorSetStagesDFS to change this.
 *
 * More details on the DFS mode:
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
 * If given a pointer to the [+] expression is given as root to this expression, it will iterate
 * the graph in a depth-first manner and stop at various stages.
 * - When entering an expression, it stops in the enterexpr stage.
 *   The SCIPexpriteratorGetParentDFS() function indicates from where the expression has been entered (NULL for the root expression).
 * - Before visiting a child of an expression, it stops in the visitingchild stage.
 *   The SCIPexpriteratorGetChildIdxDFS() function returns which child will be visited (as an index in the current expr's children array).
 *   Use SCIPexpriteratorGetChildExprDFS() to obtain the corresponding expression.
 * - When returning from visiting a child of an expression, it stops in the visitedchild stage.
 *   Again the SCIPexpriteratorGetChildExprDFS() function returns which child has been visited.
 * - When leaving an expression, it stops in the leaveexpr stage.
 *
 * Thus, for the above expression, the expression are visited in the following order and stages:
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
 * The caller can direct the iterator to skip parts of the tree:
 * If calling SCIPexpriteratorSkipDFS() in enterexpr stage, all children of that expression will be skipped. The leaveexpr stage will still be next.
 * If calling SCIPexpriteratorSkipDFS() in visitingchild stage, visiting the current child will be skipped.
 * If calling SCIPexpriteratorSkipDFS() in visitedchild child, visiting the remaining children will be skipped.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPexpriteratorInit(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*         expr,        /**< expression of the iterator, can be NULL */
   SCIP_CONSEXPRITERATOR_TYPE  type,        /**< type of expression iterator */
   SCIP_Bool                   allowrevisit /**< whether expression are allowed to be visited more than once */
   );

/** restarts an already initialized expression iterator in DFS mode
 *
 * The expression iterator will continue from the given expression, not revisiting expressions that
 * this iterator has already been visited (if initialized with allowrevisit==FALSE) and giving access
 * to the same iterator specified expression data that may have been set already.
 * Also the stop-stages are not reset.
 *
 * If revisiting is forbidden and given expr has already been visited, then the iterator will behave
 * as on the end of iteration (IsEnd() is TRUE).
 * If the enterexpr stage is not one of the stop stages, then the iterator will be moved forward
 * (GetNext() is called).
 *
 * @return The current expression.
 */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR* SCIPexpriteratorRestartDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*         expr         /**< expression of the iterator */
   );

/** specifies in which stages to stop a DFS iterator
 *
 * @param stopstages should be a bitwise OR of different SCIP_CONSEXPRITERATOR_STAGE values
 *
 * If the current stage is not one of the stopstages, then the iterator will be moved on.
 */
SCIP_EXPORT
void SCIPexpriteratorSetStagesDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   SCIP_CONSEXPRITERATOR_STAGE stopstages   /**< the stages in which to stop when iterating via DFS */
   );

/** gets the current expression that the expression iterator points to */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetCurrent(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gets the current stage that the expression iterator is in when using DFS
 *
 * If the iterator has finished (IsEnd() is TRUE), then the stage is undefined.
 */
SCIP_EXPORT
SCIP_CONSEXPRITERATOR_STAGE SCIPexpriteratorGetStageDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gets the child index that the expression iterator considers when in DFS mode and stage visitingchild or visitedchild */
SCIP_EXPORT
int SCIPexpriteratorGetChildIdxDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gets the child expression that the expression iterator considers when in DFS mode and stage visitingchild or visitedchild */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetChildExprDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gives the parent of the current expression of an expression iteration if in DFS mode
 *
 * @return the expression from which the current expression has been accessed
 */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetParentDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gives the iterator specific user data of the current expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPORT
SCIP_CONSEXPRITERATOR_USERDATA SCIPexpriteratorGetCurrentUserData(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gives the iterator specific user data of the current expressions current child
 *
 * @note The expression iterator mode must be in DFS mode and stage visitingchild or visitedchild
 */
SCIP_EXPORT
SCIP_CONSEXPRITERATOR_USERDATA SCIPexpriteratorGetChildUserDataDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** gives the iterator specific user data of a given expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPORT
SCIP_CONSEXPRITERATOR_USERDATA SCIPexpriteratorGetExprUserData(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*         expr         /**< expression for which to get the userdata of this iterator */
   );

/** sets the iterator specific user data of the current expression for an expression iteration if in DFS mode
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPORT
void SCIPexpriteratorSetCurrentUserData(
   SCIP_CONSEXPR_ITERATOR*         iterator, /**< expression iterator */
   SCIP_CONSEXPRITERATOR_USERDATA  userdata  /**< data to be stored */
   );

/** sets the iterator specific user data of a given expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPORT
void SCIPexpriteratorSetExprUserData(
   SCIP_CONSEXPR_ITERATOR*         iterator, /**< expression iterator */
   SCIP_CONSEXPR_EXPR*             expr,     /**< expression where to set iterator data */
   SCIP_CONSEXPRITERATOR_USERDATA  userdata  /**< data to be stored in current child */
   );


/** sets the iterator specific user data of the current expressions current child
 *
 * @note The expression iterator mode must be in DFS mode and stage visitingchild or visitedchild
 */
SCIP_EXPORT
void SCIPexpriteratorSetChildUserData(
   SCIP_CONSEXPR_ITERATOR*         iterator, /**< expression iterator */
   SCIP_CONSEXPRITERATOR_USERDATA  userdata  /**< data to be stored in current child */
   );

/** moves the iterator to the next expression according to the mode of the expression iterator
 *
 * @return the next expression, if any, and NULL otherwise
 */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetNext(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** moves a DFS iterator to one of the next expressions
 *
 * If in ENTEREXPR stage, then all children of that expression will be skipped.
 *   If LEAVEEXPR is one of the stopstages, then it will be the next stage. Otherwise, the iterator will move further on (go the parent, etc).
 * If in VISITINGCHILD stage, then the child that was going to be visited next will be skipped and the iterator will be moved on to the next child (if any).
 * If in VISITEDCHILD stage, then all remaining children will be skipped and we move on to the LEAVEEXPR stage (if a stop stage, otherwise further on).
 * It is not allowed to call this function when in LEAVEEXPR stage.
 *
 * @return the next expression, if any, and NULL otherwise
 */
SCIP_EXPORT
SCIP_CONSEXPR_EXPR* SCIPexpriteratorSkipDFS(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** returns whether the iterator visited all expressions already */
SCIP_EXPORT
SCIP_Bool SCIPexpriteratorIsEnd(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   );

/** @} */

#ifdef __cplusplus
}
#endif

#endif
