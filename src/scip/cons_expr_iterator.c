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
void ensureStackSize(
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

      SCIP_ALLOC_ABORT( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsexprs, iterator->dfssize, newsize) );
      SCIP_ALLOC_ABORT( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsnvisited, iterator->dfssize, newsize) );
      iterator->dfssize = newsize;
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

   ensureStackSize(iterator, iterator->dfsnexprs + 1);
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
      case SCIP_CONSEXPREXPRWALK_ENTEREXPR:
      {
         /* if there is a child, then go into visitingchild stage, otherwise go to leave stage */
         if( iterator->curr->nchildren > 0 )
         {
            iterator->dfsstage = SCIP_CONSEXPREXPRWALK_VISITINGCHILD;
            assert(iterdata->currentchild == 0);
         }
         else
         {
            iterator->dfsstage = SCIP_CONSEXPREXPRWALK_LEAVEEXPR;
         }

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

      case SCIP_CONSEXPREXPRWALK_VISITEDCHILD:
      {
         /* visit next child, if any, otherwise leave expr */
         if( ++iterdata->currentchild < iterator->curr->nchildren )
         {
            iterator->dfsstage = SCIP_CONSEXPREXPRWALK_VISITINGCHILD;
         }
         else
         {
            iterator->dfsstage = SCIP_CONSEXPREXPRWALK_LEAVEEXPR;
         }

         return iterator->curr;
      }

      case SCIP_CONSEXPREXPRWALK_LEAVEEXPR:
      {
         /* go back to parent expression, be in visitedchild stage */
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
   BMS_BLKMEM*                 blkmem,      /**< block memory used to store hash map entries */
   SCIP_CONSEXPRITERATOR_TYPE  type         /**< type of expression iterator */
   )
{
   assert(iterator != NULL);
   assert(blkmem  != NULL);

   SCIP_ALLOC( BMSallocClearBlockMemory(blkmem, iterator) );

   (*iterator)->itertype = type;
   (*iterator)->blkmem = blkmem;
   (*iterator)->iterindex = -1;
   (*iterator)->visitedtag = 0;

   /* allocate memory for type-specific data structure */
   switch( type )
   {
      case SCIP_CONSEXPRITERATOR_RTOPOLOGIC:
      {
         ensureStackSize(*iterator, MINDFSSIZE);
         break;
      }

      case SCIP_CONSEXPRITERATOR_BFS:
      {
         SCIP_CALL( SCIPqueueCreate(&(*iterator)->queue, MINBFSSIZE, 2.0) );
         break;
      }

      case SCIP_CONSEXPRITERATOR_DFS:
      {
         break;
      }
   }

   return SCIP_OKAY;
}

/** creates a more powerful expression iterator */
SCIP_RETCODE SCIPexpriteratorCreate2(
   SCIP_CONSEXPR_ITERATOR**    iterator,    /**< buffer to store expression iterator */
   BMS_BLKMEM*                 blkmem,      /**< block memory used to store hash map entries */
   SCIP_CONSEXPRITERATOR_TYPE  type,        /**< type of expression iterator */
   int                         iterindex,   /**< index of iteration data in expressions */
   unsigned int                visitedtag   /**< tag to mark or recognize visited expressions, or 0 if allow revisiting */
   )
{
   assert(iterator != NULL);
   assert(blkmem  != NULL);

   SCIP_CALL( SCIPexpriteratorCreate(iterator, blkmem, type) );
   assert(*iterator != NULL);

   /* remember where in the expressions we find our data and what marker to use to mark expr as visisted */
   (*iterator)->iterindex = iterindex;
   (*iterator)->visitedtag = visitedtag;

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

   if( (*iterator)->queue != NULL )
   {
      SCIPqueueFree(&(*iterator)->queue);
   }

   /* free iterator arrays */
   BMSfreeBlockMemoryArrayNull((*iterator)->blkmem, &(*iterator)->dfsnvisited, (*iterator)->dfssize);
   BMSfreeBlockMemoryArrayNull((*iterator)->blkmem, &(*iterator)->dfsexprs, (*iterator)->dfssize);

   /* free iterator */
   BMSfreeBlockMemory((*iterator)->blkmem, iterator);
}

/** initializes an expression iterator */
SCIP_CONSEXPR_EXPR* SCIPexpriteratorInit(
   SCIP_CONSEXPR_ITERATOR*    iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*        expr         /**< expression of the iterator */
   )
{
   assert(iterator != NULL);
   assert(expr != NULL);

   switch( iterator->itertype )
   {
      case SCIP_CONSEXPRITERATOR_BFS:
      {
         assert(iterator->queue != NULL);
         SCIPqueueClear(iterator->queue);
         SCIP_CALL_ABORT( SCIPqueueInsert(iterator->queue, expr) );

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

   /* return current expression */
   return iterator->curr;
}

/** gets the current expression that the expression iterator points to */
SCIP_CONSEXPR_EXPR* SCIPexpriteratorGetCurrent(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   )
{
   return iterator->curr;
}

/** gets the next expression according to the mode of the expression iterator */
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
         else
         {
            iterator->curr = doReverseTopologicalNext(iterator);
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
