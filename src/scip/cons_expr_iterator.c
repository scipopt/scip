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

/**@file   cons_expr_iterator.c
 * @brief  expression tree iterators
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "scip/cons_expr_iterator.h"
#include "scip/cons_expr.h"
#include "scip/struct_cons_expr.h"

#define MINDFSSIZE                       16  /**< minimum stack size for DFS*/
#define MINBFSSIZE                       16  /**< minimum queue size for BFS */

/*
 * Local methods
 */

/** ensures minimum stack size of iterator's data */
static
SCIP_RETCODE ensureStackSize(
   SCIP_CONSEXPR_ITERATOR*    iterator,     /**< expression iterator */
   int                        size          /**< minimum requires size */
   )
{
   assert(iterator != NULL);
   assert(iterator->blkmem != NULL);
   assert(iterator->itertype == SCIP_CONSEXPRITERATOR_RTOPOLOGIC);
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

static
void deinit(
   SCIP_CONSEXPR_ITERATOR*    iterator      /**< expression iterator */
   )
{
   assert(iterator != NULL );

   if( !iterator->initialized )
      return;

   if( iterator->iterindex >= 0 )
   {
      /* tell conshdlr that this iterator is no longer active */
      SCIPdeactivateConsExprExprHdlrIterator(iterator->consexprhdlr, iterator->iterindex);
      iterator->iterindex = -1;
   }

   switch( iterator->itertype )
   {
      case SCIP_CONSEXPRITERATOR_BFS :
      {
         assert(iterator->queue != NULL);

         SCIPqueueFree(&iterator->queue);

         break;
      }

      case SCIP_CONSEXPRITERATOR_RTOPOLOGIC :
      {
         assert(iterator->dfsnvisited != NULL);
         assert(iterator->dfsexprs != NULL);

         /* free dfs arrays */
         BMSfreeBlockMemoryArray(iterator->blkmem, &iterator->dfsnvisited, iterator->dfssize);
         BMSfreeBlockMemoryArray(iterator->blkmem, &iterator->dfsexprs, iterator->dfssize);
         iterator->dfssize = 0;

         break;
      }

      default: break;
   }
}

/** adds an expression to the DFS stack */
static
void reverseTopologicalInsert(
   SCIP_CONSEXPR_ITERATOR*    iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*        expr         /**< expression */
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
SCIP_CONSEXPR_EXPR* doReverseTopologicalNext(
   SCIP_CONSEXPR_ITERATOR*    iterator     /**< expression iterator */
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   int childidx;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_CONSEXPRITERATOR_RTOPOLOGIC);

   /* no expression left */
   if( iterator->dfsnexprs == 0 )
      return NULL;

   /* get expression on the top of the stack */
   expr = iterator->dfsexprs[iterator->dfsnexprs - 1];
   childidx = iterator->dfsnvisited[iterator->dfsnexprs - 1];

   /* remove the expression if all children have been visited */
   if( childidx >= SCIPgetConsExprExprNChildren(expr) )
   {
      --(iterator->dfsnexprs);
      return expr;
   }
   /* go to the next child */
   else
   {
      SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[childidx];
      assert(child != NULL);

      /* mark that the child has been visited */
      ++(iterator->dfsnvisited[iterator->dfsnexprs-1]);

      /* do left-most step */
      while( SCIPgetConsExprExprNChildren(child) > 0 )
      {
         /* add child to the DFS stack */
         reverseTopologicalInsert(iterator, child);

         /* mark that the child has been visited; note that child is on top of the DFS stack */
         ++(iterator->dfsnvisited[iterator->dfsnexprs-1]);

         child = SCIPgetConsExprExprChildren(child)[0];
      }

      /* return last child; NOTE this child is not been added to the stack */
      return child;
   }
}

/** moves to the next expression according to the BFS rule */
static
SCIP_CONSEXPR_EXPR* doBfsNext(
   SCIP_CONSEXPR_ITERATOR*    iterator     /**< expression iterator */
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   int i;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_CONSEXPRITERATOR_BFS);
   assert(iterator->queue != NULL);

   /* no expression left */
   if( SCIPqueueIsEmpty(iterator->queue) )
      return NULL;

   expr = (SCIP_CONSEXPR_EXPR*) SCIPqueueRemove(iterator->queue);
   assert(expr != NULL);

   assert(iterator->visitedtag == 0 || iterator->iterindex >= 0);
   assert(iterator->visitedtag == 0 || iterator->iterindex < SCIP_CONSEXPR_MAXNITER);
   /* we should have set the visitedtag when adding the expression to the queue */
   assert(iterator->visitedtag == 0 || expr->iterdata[iterator->iterindex].visitedtag == iterator->visitedtag);

   /* add all (possibly non-visited) children to the queue */
   for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
   {
      SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[i];
      assert(child != NULL);

      if( iterator->visitedtag != 0 )
      {
         assert(iterator->iterindex >= 0);
         assert(iterator->iterindex < SCIP_CONSEXPR_MAXNITER);

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
SCIP_CONSEXPR_EXPR* doDfsNext(
   SCIP_CONSEXPR_ITERATOR*    iterator     /**< expression iterator */
   )
{
   SCIP_CONSEXPR_EXPR_ITERDATA* iterdata;

   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_CONSEXPRITERATOR_DFS);
   assert(iterator->iterindex >= 0);

   if( iterator->curr == NULL )
      return NULL;

   iterdata = &iterator->curr->iterdata[iterator->iterindex];

   switch( iterator->dfsstage )
   {
      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD:
         /* consider next child */
         ++iterdata->currentchild;
         /* fall through */
         /* no break */

      case SCIP_CONSEXPREXPRWALK_ENTEREXPR:
      {
         /* if there is an unvisited child (left), then go into visitingchild stage, otherwise go to leave stage */
         iterator->dfsstage = SCIP_CONSEXPREXPRWALK_LEAVEEXPR;  /* expect that we will leave expr, and change mind to visitingchild below */
         while( iterdata->currentchild < iterator->curr->nchildren )
         {
            if( iterator->visitedtag == 0 || iterator->visitedtag != iterator->curr->children[iterdata->currentchild]->iterdata[iterator->iterindex].visitedtag )
            {
               /* if visitedtag is not used or child "currentchild" has not been visited yet, then go into visitingchild stage for this child */
               iterator->dfsstage = SCIP_CONSEXPREXPRWALK_VISITINGCHILD;
               break;
            }
            ++iterdata->currentchild;
         }
         assert(iterator->dfsstage == SCIP_CONSEXPREXPRWALK_VISITINGCHILD || iterdata->currentchild == iterator->curr->nchildren); /* if leaving expr, then currentchild should be at nchildren */
         assert(iterator->dfsstage == SCIP_CONSEXPREXPRWALK_LEAVEEXPR || iterdata->currentchild < iterator->curr->nchildren); /* if visiting child, then currentchild should be a valid index */
         assert(iterator->dfsstage == SCIP_CONSEXPREXPRWALK_LEAVEEXPR || iterator->visitedtag == 0 || iterator->visitedtag != iterator->curr->children[iterdata->currentchild]->iterdata[iterator->iterindex].visitedtag); /* if visiting child, then either we don't care whether we visited it already or it has not been visited yet */

         return iterator->curr;
      }

      case SCIP_CONSEXPREXPRWALK_VISITINGCHILD:
      {
         SCIP_CONSEXPR_EXPR* child;

         assert(iterdata->currentchild < iterator->curr->nchildren);

         /* remember the parent and set the first child that should be visited of the new root */
         child = iterator->curr->children[iterdata->currentchild];
         child->iterdata[iterator->iterindex].parent = iterator->curr;
         child->iterdata[iterator->iterindex].currentchild = 0;

         /* visit child */
         iterator->dfsstage = SCIP_CONSEXPREXPRWALK_ENTEREXPR;

         return child;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR:
      {
         /* go back to parent expression */

         /* remember that this expression has been visited */
         iterdata->visitedtag = iterator->visitedtag;

         /* be in visitedchild stage for the parent */
         iterator->dfsstage = SCIP_CONSEXPREXPRWALK_VISITEDCHILD;

         return iterdata->parent;
      }

      default:
         /* unknown stage */
         SCIPABORT();
         return NULL;
   }
}

/*
 * Interface methods
 */

/** creates an expression iterator */
SCIP_RETCODE SCIPexpriteratorCreate(
   SCIP_CONSEXPR_ITERATOR**    iterator,    /**< buffer to store expression iterator */
   SCIP_CONSHDLR*              consexprhdlr,/**< expr constraint handler, might be NULL */
   BMS_BLKMEM*                 blkmem       /**< block memory used to store hash map entries */
   )
{
   assert(iterator != NULL);
   assert(blkmem  != NULL);

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, iterator) );

   (*iterator)->blkmem = blkmem;
   (*iterator)->consexprhdlr = consexprhdlr;

   return SCIP_OKAY;
}

/** frees an expression iterator */
void SCIPexpriteratorFree(
   SCIP_CONSEXPR_ITERATOR**    iterator     /**< pointer to the expression iterator */
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

/** initializes an expression iterator
 *
 * \note If no conshdlr has been given when creating the iterator, then allowrevisit must be TRUE and type must not be DFS.
 */
SCIP_RETCODE SCIPexpriteratorInit(
   SCIP_CONSEXPR_ITERATOR*     iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*         expr,        /**< expression of the iterator */
   SCIP_CONSEXPRITERATOR_TYPE  type,        /**< type of expression iterator */
   SCIP_Bool                   allowrevisit /**< whether expression are allowed to be visited more than once */
   )
{
   assert(iterator != NULL);
   assert(expr != NULL);

   deinit(iterator);

   /* store the new type of the iterator */
   iterator->itertype = type;

   /* get iterindex, if necessary */
   if( !allowrevisit || type == SCIP_CONSEXPRITERATOR_DFS )
   {
      assert(iterator->consexprhdlr != NULL);

      SCIP_CALL( SCIPactivateConsExprExprHdlrIterator(iterator->consexprhdlr, &iterator->iterindex) );
   }
   else
   {
      iterator->iterindex = -1;
   }

   /* get new tag to recognize visited expressions */
   if( !allowrevisit )
   {
      assert(iterator->consexprhdlr != NULL);
      iterator->visitedtag = SCIPgetConsExprExprHdlrNewVisitedTag(iterator->consexprhdlr);
   }
   else
   {
      iterator->visitedtag = 0;
   }

   switch( iterator->itertype )
   {
      case SCIP_CONSEXPRITERATOR_BFS:
      {
         SCIP_CALL( SCIPqueueCreate(&iterator->queue, MINBFSSIZE, 2.0) );

         assert(iterator->queue != NULL);
         SCIPqueueClear(iterator->queue);
         SCIP_CALL( SCIPqueueInsert(iterator->queue, expr) );

         if( iterator->visitedtag != 0 )
         {
            assert(iterator->iterindex >= 0);
            assert(iterator->iterindex < SCIP_CONSEXPR_MAXNITER);
            assert(expr->iterdata[iterator->iterindex].visitedtag != iterator->visitedtag);

            /* mark expression as being in the queue */
            expr->iterdata[iterator->iterindex].visitedtag = iterator->visitedtag;
         }

         iterator->curr = SCIPexpriteratorGetNext(iterator);
         break;
      }

      case SCIP_CONSEXPRITERATOR_RTOPOLOGIC :
      {
         SCIP_CALL( ensureStackSize(iterator, MINDFSSIZE) );

         reverseTopologicalInsert(iterator, expr);

         iterator->curr = SCIPexpriteratorGetNext(iterator);
         break;
      }

      case SCIP_CONSEXPRITERATOR_DFS :
      {
         assert(iterator->iterindex >= 0);

         iterator->curr = expr;
         expr->iterdata[iterator->iterindex].currentchild = 0;
         expr->iterdata[iterator->iterindex].parent = NULL;
         iterator->dfsstage = SCIP_CONSEXPREXPRWALK_ENTEREXPR;

         break;
      }
   }

   iterator->initialized = TRUE;

   return SCIP_OKAY;
}

/** gets the current expression that the expression iterator points to */
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetCurrent(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   )
{
   return iterator->curr;
}

/** moves the iterator to the next expression according to the mode of the expression iterator
 *
 * @return the next expression, if any, and NULL otherwise
 */
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetNext(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   )
{
   /* move to the next expression according to iterator type */
   switch( iterator->itertype )
   {
      case SCIP_CONSEXPRITERATOR_BFS:
      {
         iterator->curr = doBfsNext(iterator);
         break;
      }

      case SCIP_CONSEXPRITERATOR_RTOPOLOGIC :
      {
         iterator->curr = doReverseTopologicalNext(iterator);
         if( iterator->visitedtag != 0 )
         {
            assert(iterator->iterindex >= 0);
            assert(iterator->iterindex < SCIP_CONSEXPR_MAXNITER);

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

      case SCIP_CONSEXPRITERATOR_DFS :
      {
         assert(iterator->iterindex >= 0);

         /* get next until we are in enterexpr state again
          * this will give every expression only once at time it is first "entered"
          * TODO the user should customize which stages it want's to see
          */
         do
         {
            iterator->curr = doDfsNext(iterator);
         }
         while( iterator->curr != NULL && iterator->dfsstage != SCIP_CONSEXPREXPRWALK_ENTEREXPR );

         break;
      }
   }

   return iterator->curr;
}

/** returns whether the iterator visited all expressions already */
SCIP_Bool SCIPexpriteratorIsEnd(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   )
{
   assert(iterator != NULL);

   return iterator->curr == NULL;
}
