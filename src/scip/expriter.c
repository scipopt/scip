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
#include "scip/struct_expr.h"
#include "scip/struct_stat.h"

/*
 * local functions
 */

/** gets the index an expression iterator can use to store iterator specific data in an expression */
SCIP_RETCODE activate(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL );

   if( iterator->stat->nactiveexpriter + 1 >= SCIP_EXPRITER_MAXNACTIVE )
   {
      SCIPerrorMessage("Maximal number of active expression iterators reached.\n");
      return SCIP_MAXDEPTHLEVEL;
   }

   iterator->iterindex = iterator->stat->nactiveexpriter++;

   return SCIP_OKAY;
}

/** returns the index that an expression iterator used to store iterator specific data in an expression for other iterators to use */
static
void deactivate(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL );
   assert(iterator->iterindex >= 0);

   /* the iterindex must be the one of the last initialized iterator */
   assert(iterator->iterindex == iterator->stat->nactiveexpriter-1);

   --iterator->stat->nactiveexpriter;
}

/** get a new tag that can be used to mark an expression as visited */
static
SCIP_Longint newVisitedTag(
   SCIP_EXPRITER*        iterator            /**< expression iterator */
   )
{
   assert(iterator != NULL );

   return ++iterator->stat->exprlastvisitedtag;
}

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
      /* tell core that this iterator is no longer active */
      deactivate(iterator);
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

/* private functions (expr.h) */

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

/* public functions (pub_expr.h) */
