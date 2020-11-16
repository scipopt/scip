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

/**@file   expriter.c
 * @brief  functions for iterating over algebraic expressions
 * @author Benjamin Mueller
 * @author Stefan Vigerske
 */

#include <assert.h>

#include "scip/expr.h"
#include "scip/pub_misc.h"
#include "scip/struct_expr.h"
#include "scip/struct_stat.h"

#define MINDFSSIZE  16 /**< minimum stack size for DFS*/
#define MINBFSSIZE  16 /**< minimum queue size for BFS */

/*
 * local functions
 */

static
void deinit(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL );

   if( !iterator->initialized )
      return;

   if( iterator->iterindex >= 0 )
   {
      /* the iterindex must be the one of the last initialized iterator */
      assert(iterator->iterindex == iterator->stat->nactiveexpriter-1);

      /* tell core that this iterator is no longer active */
      --iterator->stat->nactiveexpriter;

      iterator->iterindex = -1;
   }

   switch( iterator->itertype )
   {
      case SCIP_EXPRITER_BFS :
      {
         assert(iterator->queue != NULL);

         SCIPqueueFree(&iterator->queue);

         break;
      }

      case SCIP_EXPRITER_RTOPOLOGIC :
      {
         assert(iterator->dfsnvisited != NULL);
         assert(iterator->dfsexprs != NULL);

         /* free dfs arrays */
         BMSfreeBlockMemoryArray(iterator->blkmem, &iterator->dfsnvisited, iterator->dfssize);
         BMSfreeBlockMemoryArray(iterator->blkmem, &iterator->dfsexprs, iterator->dfssize);
         iterator->dfssize = 0;

         break;
      }

      case SCIP_EXPRITER_DFS :
      default: break;
   }
}

/** ensures minimum stack size of iterator's data */
static
SCIP_RETCODE ensureStackSize(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   int                   size                /**< minimum requires size */
   )
{
   assert(iterator != NULL);
   assert(iterator->blkmem != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_RTOPOLOGIC);
   assert(size >= 0);

   if( size > iterator->dfssize )
   {
      int newsize = size * 2;

      SCIP_ALLOC( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsexprs, iterator->dfssize, newsize) );
      SCIP_ALLOC( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsnvisited, iterator->dfssize, newsize) );
      iterator->dfssize = newsize;
   }

   return SCIP_OKAY;
}

/** adds an expression to the DFS stack */
static
void reverseTopologicalInsert(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr                /**< expression */
   )
{
   assert(iterator != NULL);
   assert(expr != NULL);

   SCIP_CALL_ABORT( ensureStackSize(iterator, iterator->dfsnexprs + 1) );
   iterator->dfsexprs[iterator->dfsnexprs] = expr;
   iterator->dfsnvisited[iterator->dfsnexprs] = 0;
   ++(iterator->dfsnexprs);
}

/** moves to the next expression according to a reverse topological order */
static
SCIP_EXPR* doReverseTopologicalNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   SCIP_EXPR* expr;
   int childidx;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_RTOPOLOGIC);

   /* no expression left */
   if( iterator->dfsnexprs == 0 )
      return NULL;

   /* get expression on the top of the stack */
   expr = iterator->dfsexprs[iterator->dfsnexprs - 1];
   childidx = iterator->dfsnvisited[iterator->dfsnexprs - 1];

   /* remove the expression if all children have been visited */
   if( childidx >= SCIPexprGetNChildren(expr) )
   {
      --(iterator->dfsnexprs);
      return expr;
   }
   /* go to the next child */
   else
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[childidx];
      assert(child != NULL);

      /* mark that the child has been visited */
      ++(iterator->dfsnvisited[iterator->dfsnexprs-1]);

      /* do left-most step */
      while( SCIPexprGetNChildren(child) > 0 )
      {
         /* add child to the DFS stack */
         reverseTopologicalInsert(iterator, child);

         /* mark that the child has been visited; note that child is on top of the DFS stack */
         ++(iterator->dfsnvisited[iterator->dfsnexprs-1]);

         child = SCIPexprGetChildren(child)[0];
      }

      /* return last child; NOTE this child is not been added to the stack */
      return child;
   }
}

/** moves to the next expression according to the BFS rule */
static
SCIP_EXPR* doBfsNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   SCIP_EXPR* expr;
   int i;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_BFS);
   assert(iterator->queue != NULL);

   /* no expression left */
   if( SCIPqueueIsEmpty(iterator->queue) )
      return NULL;

   expr = (SCIP_EXPR*) SCIPqueueRemove(iterator->queue);
   assert(expr != NULL);

   assert(iterator->visitedtag == 0 || iterator->iterindex >= 0);
   assert(iterator->visitedtag == 0 || iterator->iterindex < SCIP_EXPRITER_MAXNACTIVE);
   /* we should have set the visitedtag when adding the expression to the queue */
   assert(iterator->visitedtag == 0 || expr->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag);

   /* add all (possibly non-visited) children to the queue */
   for( i = 0; i < SCIPexprGetNChildren(expr); ++i )
   {
      SCIP_EXPR* child = SCIPexprGetChildren(expr)[i];
      assert(child != NULL);

      if( iterator->visitedtag != 0 )
      {
         assert(iterator->iterindex >= 0);
         assert(iterator->iterindex < SCIP_EXPRITER_MAXNACTIVE);

         /* skip children that have already been visited or have already been added to the queue */
         if( child->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag )
            continue;

         /* mark child as being in the queue (will be inserted next) */
         child->iterdata[iterator->iterindex].visitedtag = iterator->visitedtag;
      }

      /* add child to the queue */
      SCIP_CALL_ABORT( SCIPqueueInsert(iterator->queue, child) );
   }

   return expr;
}

static
SCIP_EXPR* doDfsNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   SCIP_EXPRITERDATA* iterdata;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->iterindex >= 0);

   if( iterator->curr == NULL )
      return NULL;

   iterdata = &iterator->curr->iterdata[iterator->iterindex];

   switch( iterator->dfsstage )
   {
      case SCIP_EXPRITER_VISITEDCHILD:
         /* consider next child */
         ++iterdata->currentchild;
         /* fall through */ /* no break */ /*lint -fallthrough*/

      case SCIP_EXPRITER_ENTEREXPR:
      {
         /* if there is an unvisited child (left), then go into visitingchild stage, otherwise go to leave stage */
         iterator->dfsstage = SCIP_EXPRITER_LEAVEEXPR;  /* expect that we will leave expr, and change mind to visitingchild below */
         while( iterdata->currentchild < iterator->curr->nchildren )
         {
            if( iterator->visitedtag == 0 || iterator->visitedtag != iterator->curr->children[iterdata->currentchild]->iterdata[iterator->iterindex].visitedtag )
            {
               /* if visitedtag is not used or child "currentchild" has not been visited yet, then go into visitingchild stage for this child */
               iterator->dfsstage = SCIP_EXPRITER_VISITINGCHILD;
               break;
            }
            ++iterdata->currentchild;
         }
         assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterdata->currentchild == iterator->curr->nchildren); /* if leaving expr, then currentchild should be at nchildren */
         assert(iterator->dfsstage == SCIP_EXPRITER_LEAVEEXPR || iterdata->currentchild < iterator->curr->nchildren); /* if visiting child, then currentchild should be a valid index */
         assert(iterator->dfsstage == SCIP_EXPRITER_LEAVEEXPR || iterator->visitedtag == 0 || iterator->visitedtag != iterator->curr->children[iterdata->currentchild]->iterdata[iterator->iterindex].visitedtag); /* if visiting child, then either we don't care whether we visited it already or it has not been visited yet */

         return iterator->curr;
      }

      case SCIP_EXPRITER_VISITINGCHILD:
      {
         SCIP_EXPR* child;

         assert(iterdata->currentchild < iterator->curr->nchildren);

         /* remember the parent and set the first child that should be visited of the new root */
         child = iterator->curr->children[iterdata->currentchild];
         child->iterdata[iterator->iterindex].parent = iterator->curr;
         child->iterdata[iterator->iterindex].currentchild = 0;

         /* visit child */
         iterator->dfsstage = SCIP_EXPRITER_ENTEREXPR;

         return child;
      }

      case SCIP_EXPRITER_LEAVEEXPR:
      {
         /* go back to parent expression */

         /* remember that this expression has been visited */
         iterdata->visitedtag = iterator->visitedtag;

         /* be in visitedchild stage for the parent */
         iterator->dfsstage = SCIP_EXPRITER_VISITEDCHILD;

         return iterdata->parent;
      }

      default:
         /* unknown stage */
         SCIPABORT();
         return NULL;
   }
}

/*
 * private functions (expr.h)
 */

/** creates an expression iterator */
SCIP_RETCODE SCIPexpriterCreate(
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_EXPRITER**       iterator            /**< buffer to store expression iterator */
   )
{
   assert(iterator != NULL);
   assert(blkmem  != NULL);

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, iterator) );

   (*iterator)->stat = stat;
   (*iterator)->blkmem = blkmem;

   return SCIP_OKAY;
}

/** frees an expression iterator */
void SCIPexpriterFree(
   SCIP_EXPRITER**       iterator            /**< pointer to the expression iterator */
   )
{
   assert(iterator != NULL);
   assert(*iterator != NULL);
   assert((*iterator)->blkmem != NULL);

   deinit(*iterator);

   assert((*iterator)->queue == NULL);
   assert((*iterator)->dfsnvisited == NULL);
   assert((*iterator)->dfsexprs == NULL);

   /* free iterator */
   BMSfreeBlockMemory((*iterator)->blkmem, iterator);
}

/*
 * public functions (pub_expr.h)
 */

/** returns whether expression iterator is current initialized */
SCIP_Bool SCIPexpriterIsInit(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);

   return iterator->initialized;
}

/** initializes an expression iterator
 *
 * @note If no conshdlr has been given when creating the iterator, then allowrevisit must be TRUE and type must not be DFS.
 *
 * @note If expr is NULL, then iterator will be ended (SCIPexpriterIsEnd() is TRUE). Useful if following with SCIPexpriterRestartDFS().
 *
 * If type is DFS, then stopstages will be set to ENTEREXPR. Use SCIPexpriterSetStagesDFS to change this.
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
 *   The SCIPexpriterGetParentDFS() function indicates from where the expression has been entered (NULL for the root expression).
 * - Before visiting a child of an expression, it stops in the visitingchild stage.
 *   The SCIPexpriterGetChildIdxDFS() function returns which child will be visited (as an index in the current expr's children array).
 *   Use SCIPexpriterGetChildExprDFS() to obtain the corresponding expression.
 * - When returning from visiting a child of an expression, it stops in the visitedchild stage.
 *   Again the SCIPexpriterGetChildExprDFS() function returns which child has been visited.
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
 * If calling SCIPexpriterSkipDFS() in enterexpr stage, all children of that expression will be skipped. The leaveexpr stage will still be next.
 * If calling SCIPexpriterSkipDFS() in visitingchild stage, visiting the current child will be skipped.
 * If calling SCIPexpriterSkipDFS() in visitedchild child, visiting the remaining children will be skipped.
 */
SCIP_RETCODE SCIPexpriterInit(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr,               /**< expression of the iterator, can be NULL */
   SCIP_EXPRITER_TYPE    type,               /**< type of expression iterator */
   SCIP_Bool             allowrevisit        /**< whether expression are allowed to be visited more than once */
   )
{
   assert(iterator != NULL);

   deinit(iterator);

   /* store the new type of the iterator */
   iterator->itertype = type;

   /* get iterindex, if necessary */
   if( !allowrevisit || type == SCIP_EXPRITER_DFS )
   {
      if( iterator->stat->nactiveexpriter + 1 >= SCIP_EXPRITER_MAXNACTIVE )
      {
         SCIPerrorMessage("Maximal number of active expression iterators reached.\n");
         return SCIP_MAXDEPTHLEVEL;
      }

      iterator->iterindex = iterator->stat->nactiveexpriter++;
   }
   else
   {
      iterator->iterindex = -1;
   }

   /* get new tag to recognize visited expressions */
   if( !allowrevisit )
   {
      iterator->visitedtag = ++iterator->stat->exprlastvisitedtag;
   }
   else
   {
      iterator->visitedtag = 0;
   }

   switch( iterator->itertype )
   {
      case SCIP_EXPRITER_BFS:
      {
         SCIP_CALL( SCIPqueueCreate(&iterator->queue, MINBFSSIZE, 2.0) );

         assert(iterator->queue != NULL);
         SCIPqueueClear(iterator->queue);

         if( expr == NULL )
         {
            iterator->curr = NULL;
            break;
         }

         SCIP_CALL( SCIPqueueInsert(iterator->queue, expr) );

         if( iterator->visitedtag != 0 )
         {
            assert(iterator->iterindex >= 0);
            assert(iterator->iterindex < SCIP_EXPRITER_MAXNACTIVE);
            assert(expr->iterdata[iterator->iterindex].visitedtag != iterator->visitedtag);

            /* mark expression as being in the queue */
            expr->iterdata[iterator->iterindex].visitedtag = iterator->visitedtag;
         }

         iterator->curr = SCIPexpriterGetNext(iterator);
         break;
      }

      case SCIP_EXPRITER_RTOPOLOGIC :
      {
         SCIP_CALL( ensureStackSize(iterator, MINDFSSIZE) );

         if( expr != NULL )
         {
            reverseTopologicalInsert(iterator, expr);
            iterator->curr = SCIPexpriterGetNext(iterator);
         }
         else
         {
            iterator->curr = NULL;
         }

         break;
      }

      case SCIP_EXPRITER_DFS :
      {
         assert(iterator->iterindex >= 0);

         iterator->stopstages = SCIP_EXPRITER_ENTEREXPR;
         iterator->curr = expr;

         if( expr == NULL )
            break;

         expr->iterdata[iterator->iterindex].currentchild = 0;
         expr->iterdata[iterator->iterindex].parent = NULL;
         iterator->dfsstage = SCIP_EXPRITER_ENTEREXPR;

         break;
      }
   }

   iterator->initialized = TRUE;

   return SCIP_OKAY;
}

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
SCIP_EXPR* SCIPexpriterRestartDFS(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr                /**< expression of the iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->initialized);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);

   /* if we forbid revisiting and root expr has already been visited, then set curr to NULL, that is, be at end of iterator */
   if( iterator->visitedtag > 0 && expr->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag )
   {
      iterator->curr = NULL;
      return NULL;
   }

   /* set current to given expr, make it the root, and set stage to enterexpr */
   iterator->curr = expr;
   expr->iterdata[iterator->iterindex].currentchild = 0;
   expr->iterdata[iterator->iterindex].parent = NULL;
   iterator->dfsstage = SCIP_EXPRITER_ENTEREXPR;

   if( (iterator->stopstages & SCIP_EXPRITER_ENTEREXPR) == 0 )
      return SCIPexpriterGetNext(iterator);

   return iterator->curr;
}

/** specifies in which stages to stop a DFS iterator
 *
 * @param stopstages should be a bitwise OR of different SCIP_EXPRITER_STAGE values
 *
 * If the current stage is not one of the stopstages, then the iterator will be moved on.
 */
void SCIPexpriterSetStagesDFS(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPRITER_STAGE   stopstages          /**< the stages in which to stop when iterating via DFS */
   )
{
   assert(iterator != NULL);

   if( (iterator->dfsstage & stopstages) == 0 )
   {
      iterator->stopstages = stopstages;
      (void) SCIPexpriterGetNext(iterator);
   }
   else
   {
      iterator->stopstages = stopstages;
   }
}

/** gets the current expression that the expression iterator points to */
SCIP_EXPR* SCIPexpriterGetCurrent(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);

   return iterator->curr;
}

/** gets the current stage that the expression iterator is in when using DFS
 *
 * If the iterator has finished (IsEnd() is TRUE), then the stage is undefined.
 */
SCIP_EXPRITER_STAGE SCIPexpriterGetStageDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);

   return iterator->dfsstage;
}

/** gets the child index that the expression iterator considers when in DFS mode and stage visitingchild or visitedchild */
int SCIPexpriterGetChildIdxDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterator->dfsstage == SCIP_EXPRITER_VISITEDCHILD);

   return iterator->curr->iterdata[iterator->iterindex].currentchild;
}

/** gets the child expression that the expression iterator considers when in DFS mode and stage visitingchild or visitedchild */
SCIP_EXPR* SCIPexpriterGetChildExprDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterator->dfsstage == SCIP_EXPRITER_VISITEDCHILD);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild >= 0);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild < iterator->curr->nchildren);

   return iterator->curr->children[iterator->curr->iterdata[iterator->iterindex].currentchild];
}

/** gives the parent of the current expression of an expression iteration if in DFS mode
 *
 * @return the expression from which the current expression has been accessed
 */
SCIP_EXPR* SCIPexpriterGetParentDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);

   return iterator->curr->iterdata[iterator->iterindex].parent;
}

/** gives the iterator specific user data of the current expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPRITER_USERDATA SCIPexpriterGetCurrentUserData(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);

   return iterator->curr->iterdata[iterator->iterindex].userdata;
}

/** gives the iterator specific user data of the current expressions current child
 *
 * @note The expression iterator mode must be in DFS mode and stage visitingchild or visitedchild
 */
SCIP_EXPRITER_USERDATA SCIPexpriterGetChildUserDataDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterator->dfsstage == SCIP_EXPRITER_VISITEDCHILD);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild >= 0);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild < iterator->curr->nchildren);

   return iterator->curr->children[iterator->curr->iterdata[iterator->iterindex].currentchild]->iterdata[iterator->iterindex].userdata;
}

/** gives the iterator specific user data of a given expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
SCIP_EXPRITER_USERDATA SCIPexpriterGetExprUserData(
   SCIP_EXPRITER*        iterator,           /**< expression iterator */
   SCIP_EXPR*            expr                /**< expression for which to get the userdata of this iterator */
   )
{
   assert(iterator != NULL);
   assert(expr != NULL);
   assert(iterator->iterindex >= 0);

   return expr->iterdata[iterator->iterindex].userdata;
}

/** sets the iterator specific user data of the current expression for an expression iteration if in DFS mode
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
void SCIPexpriterSetCurrentUserData(
   SCIP_EXPRITER*         iterator,          /**< expression iterator */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);

   iterator->curr->iterdata[iterator->iterindex].userdata = userdata;
}

/** sets the iterator specific user data of a given expression
 *
 * @note The expression iterator mode must be DFS or another mode with allowrevisit=FALSE
 */
void SCIPexpriterSetExprUserData(
   SCIP_EXPRITER*         iterator,          /**< expression iterator */
   SCIP_EXPR*             expr,              /**< expression where to set iterator data */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored in current child */
   )
{
   assert(iterator != NULL);
   assert(iterator->iterindex >= 0);

   expr->iterdata[iterator->iterindex].userdata = userdata;
}

/** sets the iterator specific user data of the current expressions current child
 *
 * @note The expression iterator mode must be in DFS mode and stage visitingchild or visitedchild
 */
void SCIPexpriterSetChildUserData(
   SCIP_EXPRITER*         iterator,          /**< expression iterator */
   SCIP_EXPRITER_USERDATA userdata           /**< data to be stored in current child */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->iterindex >= 0);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->dfsstage == SCIP_EXPRITER_VISITINGCHILD || iterator->dfsstage == SCIP_EXPRITER_VISITEDCHILD);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild >= 0);
   assert(iterator->curr->iterdata[iterator->iterindex].currentchild < iterator->curr->nchildren);

   iterator->curr->children[iterator->curr->iterdata[iterator->iterindex].currentchild]->iterdata[iterator->iterindex].userdata = userdata;
}

/** moves the iterator to the next expression according to the mode of the expression iterator
 *
 * @return the next expression, if any, and NULL otherwise
 */
SCIP_EXPR* SCIPexpriterGetNext(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   /* move to the next expression according to iterator type */
   switch( iterator->itertype )
   {
      case SCIP_EXPRITER_BFS:
      {
         iterator->curr = doBfsNext(iterator);
         break;
      }

      case SCIP_EXPRITER_RTOPOLOGIC :
      {
         iterator->curr = doReverseTopologicalNext(iterator);
         if( iterator->visitedtag != 0 )
         {
            assert(iterator->iterindex >= 0);
            assert(iterator->iterindex < SCIP_EXPRITER_MAXNACTIVE);

            /* skip already visited expressions */
            while( iterator->curr != NULL )
            {
               if( iterator->curr->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag )
               {
                  /* if curr has already been visited, get next one
                   * TODO this isn't really efficient, since we still walk through already visited expressions
                   */
                  iterator->curr = doReverseTopologicalNext(iterator);
               }
               else
               {
                  /* curr has not been visited yet, so mark it as visited and interrupt loop */
                  iterator->curr->iterdata[iterator->iterindex].visitedtag = iterator->visitedtag;
                  break;
               }
            }
         }
         break;
      }

      case SCIP_EXPRITER_DFS :
      {
         assert(iterator->iterindex >= 0);

         /* get next until we are in a stopstage again
          * this might give expressions more than once, depending on what the stopstages are
          */
         do
         {
            iterator->curr = doDfsNext(iterator);
         }
         while( iterator->curr != NULL && (iterator->dfsstage & iterator->stopstages) == 0 );

         break;
      }
   }

   return iterator->curr;
}

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
SCIP_EXPR* SCIPexpriterSkipDFS(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->curr != NULL);
   assert(iterator->itertype == SCIP_EXPRITER_DFS);
   assert(iterator->iterindex >= 0);

   switch( iterator->dfsstage )
   {
      case SCIP_EXPRITER_ENTEREXPR :
      case SCIP_EXPRITER_VISITEDCHILD :
      {
         /* move directly to leaveexpr */
         iterator->dfsstage = SCIP_EXPRITER_LEAVEEXPR;
         /* if leaveexpr is not a stopstage, then move on */
         while( iterator->curr != NULL && (iterator->dfsstage & iterator->stopstages) == 0 )
            iterator->curr = doDfsNext(iterator);
         return iterator->curr;
      }

      case SCIP_EXPRITER_VISITINGCHILD :
      {
         /* skip the child to be visited */
         /* pretend we just visited this child and get next */
         iterator->dfsstage = SCIP_EXPRITER_VISITEDCHILD;
         return SCIPexpriterGetNext(iterator);
      }

      case SCIP_EXPRITER_LEAVEEXPR :
      default :
         SCIPerrorMessage("SCIPexpriterSkipDFS called in invalid stage %d", iterator->dfsstage);
         SCIPABORT();
         return iterator->curr;
   }
}

/** returns whether the iterator visited all expressions already */
SCIP_Bool SCIPexpriterIsEnd(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL);

   return iterator->curr == NULL;
}
