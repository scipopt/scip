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
#include "scip/struct_cons_expr.h"

/*
 * Local methods
 */

/** moves to the next expression according to the DFS rule */
static
void doDfsNext(
   SCIP_CONSEXPR_ITERATOR*    iterator     /**< expression iterator */
   )
{
   /* TODO */
}

/** moves to the next expression according to the BFS rule */
static
void doBfsNext(
   SCIP_CONSEXPR_ITERATOR*    iterator     /**< expression iterator */
   )
{
   /* TODO */
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

   /* initialize variables */
   (*iterator)->itertype = itertype;
   (*iterator)->blkmem = blkmem;

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

   BMSfreeBlockMemory((*iterator)->blkmem, iterator);
}

/**< initializes an expression iterator */
void SCIPexpriteratorInit(
   SCIP_CONSEXPR_ITERATOR*    iterator,    /**< expression iterator */
   SCIP_CONSEXPR_EXPR*        expr         /**< expression of the iterator */
   )
{
   assert(iterator != NULL);
   assert(expr != NULL);

   iterator->expr = expr;
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
         doBfsNext(iterator);
         break;

      case SCIP_CONSEXPRITERATOR_DFS:
         doDfsNext(iterator);
         break;

      default:
         SCIPABORT();
   }

   return iterator->expr;
}

/** returns whether the iterator visited all expressions already */
SCIP_Bool SCIPexpriteratorIsEnd(
   SCIP_CONSEXPR_ITERATOR*     iterator     /**< expression iterator */
   )
{
   assert(iterator != NULL);
   return iterator->expr == NULL;
}
