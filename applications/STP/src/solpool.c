/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   solpool.c
 * @brief  includes solution pool for Steiner tree problems
 * @author Daniel Rehfeldt
 *
 * This file contains several basic methods for a Steiner tree solution pool.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/*lint -esym(766,stdlib.h) -esym(766,malloc.h)         */
/*lint -esym(766,string.h)                             */
//#define SCIP_DEBUG

#include "solpool.h"


/** is given solution in pool? */
SCIP_Bool solpool_isContained(
   const int*            soledges,           /**< edge array of solution to be checked */
   const STPSOLPOOL*     pool                /**< the pool */
)
{
   STPSOL** poolsols = pool->sols;
   const int poolsize = pool->size;
   const int nedges = pool->nedges;

   for( int i = 0; i < poolsize; i++ )
   {
      int j;
      const int* pooledges = poolsols[i]->soledges;
      assert(pooledges != NULL);

      for( j = 0; j < nedges; j++ )
         if( pooledges[j] != soledges[j] )
            break;

      /* pooledges == soledges? */
      if( j == nedges )
      {
         SCIPdebugMessage("Pool: solution is already contained \n");
         return TRUE;
      }
   }
   return FALSE;
}


/** get solution from index */
STPSOL* solpool_solFromIndex(
   STPSOLPOOL*           pool,               /**< the pool */
   const int             soindex             /**< the index */
)
{
   int i;
   int size;

   assert(pool != NULL);
   assert(soindex >= 0 && soindex <= pool->maxindex);

   size = pool->size;

   for( i = 0; i < size; i++ )
      if( pool->sols[i]->index == soindex )
         break;

   if( i == size )
      return NULL;
   else
      return pool->sols[i];
}


/** initializes STPSOL pool */
SCIP_RETCODE solpool_init(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL**          pool,               /**< the pool */
   const int             nedges,             /**< number of edges of solutions to be stored in the pool */
   const int             maxsize             /**< capacity of pool */
)
{
   STPSOLPOOL* dpool;

   assert(pool != NULL);
   assert(nedges > 0);
   assert(maxsize > 0);

   SCIP_CALL( SCIPallocBlockMemory(scip, pool) );

   dpool = *pool;
   SCIP_CALL( SCIPallocMemoryArray(scip, &(dpool->sols), maxsize) );

   for( int i = 0; i < maxsize; i++ )
      dpool->sols[i] = NULL;

   dpool->size = 0;
   dpool->nedges = nedges;
   dpool->maxsize = maxsize;
   dpool->maxindex = -1;

   return SCIP_OKAY;
}


/** frees STPSOL pool */
void solpool_free(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL**          pool                /**< the pool */
   )
{
   STPSOLPOOL* dpool = *pool;
   const int poolsize = dpool->size;

   assert(pool != NULL);
   assert(dpool != NULL);
   assert(poolsize == dpool->maxsize || dpool->sols[poolsize] == NULL);

   for( int i = poolsize - 1; i >= 0; i-- )
   {
      STPSOL* sol = dpool->sols[i];

      assert(sol != NULL);

      SCIPfreeMemoryArray(scip, &(sol->soledges));
      SCIPfreeBlockMemory(scip, &sol);
   }

   SCIPfreeMemoryArray(scip, &(dpool->sols));
   SCIPfreeBlockMemory(scip, pool);
}


/** tries to add STPSOL to pool */
SCIP_RETCODE solpool_addSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real       obj,                /**< objective of solution to be added */
   const int*            soledges,           /**< edge array of solution to be added */
   STPSOLPOOL*           pool,               /**< the pool */
   SCIP_Bool*            success             /**< has solution been added? */
   )
{
   STPSOL** poolsols = pool->sols;
   STPSOL* sol;
   int i;
   int poolsize = pool->size;
   const int nedges = pool->nedges;
   const int poolmaxsize = pool->maxsize;

   assert(scip != NULL);
   assert(pool != NULL);
   assert(poolsols != NULL);
   assert(poolsize >= 0);
   assert(poolmaxsize >= 0);
   assert(poolsize <= poolmaxsize);

   *success = FALSE;

   /* is solution in pool? */
   if( solpool_isContained(soledges, pool) )
      return SCIP_OKAY;

   SCIPdebugMessage("Pool: add to pool (current size: %d, max: %d) \n", poolsize, poolmaxsize);

   /* enlarge pool if possible */
   if( poolsize < poolmaxsize )
   {
      SCIPdebugMessage("Pool: alloc memory at position %d \n", poolsize);

      SCIP_CALL( SCIPallocBlockMemory(scip, &(poolsols[poolsize])) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(poolsols[poolsize]->soledges), nedges) );

      poolsize++;
      pool->size++;
   }
   /* pool is full; new solution worse than worst solution in pool? */
   else if( SCIPisGT(scip, obj, poolsols[poolsize - 1]->obj) )
   {
      return SCIP_OKAY;
   }

   /* overwrite last element of pool (either empty or inferior to current solution) */
   sol = poolsols[poolsize - 1];
   assert(sol != NULL);
   sol->obj = obj;
   sol->index = ++(pool->maxindex);
   BMScopyMemoryArray(sol->soledges, soledges, nedges);

   /* shift solution up */
   for( i = poolsize - 1; i >= 1; i-- )
   {
      if( SCIPisGT(scip, obj, poolsols[i - 1]->obj) )
         break;

      poolsols[i] = poolsols[i - 1];
   }

   poolsols[i] = sol;
   SCIPdebugMessage("Pool: put new solution to position %d \n", i);
   *success = TRUE;

   return SCIP_OKAY;
}
