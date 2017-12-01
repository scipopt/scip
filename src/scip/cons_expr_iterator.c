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

/** ensures minimum size of iterator's data */
static
void ensureStackSize(
   SCIP_CONSEXPR_ITERATOR*    iterator,     /**< expression iterator */
   int                        size          /**< minimum requires size */
   )
{
   assert(iterator != NULL);
   assert(iterator->blkmem != NULL);
   assert(size >= 0);

   if( size < iterator->dfssize )
   {
      int newsize = size * 2;

      SCIP_ALLOC_ABORT( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsexprs, iterator->dfssize, newsize) );
      SCIP_ALLOC_ABORT( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsparents, iterator->dfssize, newsize) );
      SCIP_ALLOC_ABORT( BMSreallocBlockMemoryArray(iterator->blkmem, &iterator->dfsnvisited, iterator->dfssize, newsize) );
      iterator->dfssize = newsize;
   }
}

/** moves to the next expression according to the DFS rule */
static
SCIP_RETCODE doDfsNext(
   SCIP_CONSEXPR_ITERATOR*    iterator     /**< expression iterator */
   )
{
   assert(iterator != NULL);
   assert(iterator->itertype == SCIP_CONSEXPRITERATOR_DFS);

   /* no expression left */
   if( iterator->dfsnexprs == 0 )
      return SCIP_OKAY;

   return SCIP_OKAY;
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

   /* add all children to the queue */
   for( i = 0; i < SCIPgetConsExprExprNChildren(expr); ++i )
   {
      SCIP_CONSEXPR_EXPR* child = SCIPgetConsExprExprChildren(expr)[i];
      assert(child != NULL);

      /* add child to the queue */
      SCIP_CALL_ABORT( SCIPqueueInsert(iterator->queue, child) );
   }

   return expr;
}

/*
 * Interface methods
 */

/**< creates an expression iterator */
SCIP_RETCODE SCIPexpriteratorCreate(
   SCIP_CONSEXPR_ITERATOR**    iterator,    /**< buffer to store expression iterator */
   BMS_BLKMEM*                 blkmem,      /**< block memory used to store hash map entries */
   SCIP_CONSEXPRITERATOR_TYPE  itertype     /**< type of expression iterator */
   )
{
   assert(iterator != NULL);
   assert(blkmem  != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, iterator) );
   BMSclearMemory(*iterator);

   (*iterator)->itertype = itertype;
   (*iterator)->blkmem = blkmem;

   /* allocate memory for DFS or BFS data structure */
   if( itertype == SCIP_CONSEXPRITERATOR_BFS )
   {
      SCIP_CALL( SCIPqueueCreate(&(*iterator)->queue, MINBFSSIZE, 2.0) );
   }
   else
   {
      assert(itertype == SCIP_CONSEXPRITERATOR_DFS);
      ensureStackSize(*iterator, MINDFSSIZE);
   }

   return SCIP_OKAY;
}

/**< frees an expression iterator */
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
   BMSfreeBlockMemoryArrayNull((*iterator)->blkmem, &(*iterator)->dfsparents, (*iterator)->dfssize);
   BMSfreeBlockMemoryArrayNull((*iterator)->blkmem, &(*iterator)->dfsexprs, (*iterator)->dfssize);

   /* free iterator */
   BMSfreeBlockMemory((*iterator)->blkmem, iterator);
}

/**< initializes an expression iterator */
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
         assert(iterator->queue != NULL);
         SCIPqueueClear(iterator->queue);
         SCIP_CALL_ABORT( SCIPqueueInsert(iterator->queue, expr) );
         break;

      case SCIP_CONSEXPRITERATOR_DFS:
         ensureStackSize(iterator, 1);
         iterator->dfsexprs[0] = expr;
         iterator->dfsparents[0] = NULL;
         iterator->dfsnvisited[0] = 0;
         iterator->dfsnexprs = 0;
         break;

      default:
         SCIPABORT();
   }

   /* return next expression */
   return SCIPexpriteratorGetNext(iterator);
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
         return doBfsNext(iterator);

      case SCIP_CONSEXPRITERATOR_DFS:
         SCIP_CALL_ABORT( doDfsNext(iterator) );
         break;

      default:
         SCIPABORT();
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
