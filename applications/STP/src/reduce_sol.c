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

/**@file   reduce_sol.c
 * @brief  Reduction solution storage methods for Steiner problems
 * @author Daniel Rehfeldt
 *
 * This file includes methods to save and retain solutions and primal bounds during the reduction
 * process
 *
 * A list of all interface methods can be found in reduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/*lint -esym(750,REDUCE_C) -esym(766,stdlib.h) -esym(766,string.h) */
#define SCIP_DEBUG
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "graph.h"
#include "reduce.h"
#include "scip/scip.h"
#include "portab.h"



struct reduction_primal_bound_storage
{
   SCIP_Real             offset;             /**< offset */
   SCIP_Real             primalbound;        /**< best primal bound */
   int                   decomplevel;        /**< decomposition level */
   // todo: also save level and bound per level
};


/*
 * Interface methods
 */


/** initializes */
SCIP_RETCODE reduce_primalInit(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          g,                  /**< graph data structure */
   REDPRIMAL**           primal              /**< to initialize */
   )
{
   REDPRIMAL* rp;
   assert(scip);

   SCIP_CALL( SCIPallocMemory(scip, primal) );
   rp = *primal;

   rp->offset = 0.0;
   rp->primalbound = FARAWAY;
   rp->decomplevel = 0;

   return SCIP_OKAY;
}


/** frees */
void reduce_primalFree(
   SCIP*                 scip,               /**< SCIP data structure */
   REDPRIMAL**           primal              /**< to free */
   )
{
   assert(scip && primal);
   assert(*primal);


   SCIPfreeMemory(scip, primal);
}

/** sets offset */
void reduce_primalSetOffset(
   SCIP_Real             offsetnew,          /**< new offset */
   REDPRIMAL*            primal              /**< primal */
   )
{
   assert(primal);
   assert(GE(offsetnew, 0.0));

   primal->offset = offsetnew;
}


/** sets level */
void reduce_primalSetLevel(
   int                   level,              /**< level */
   REDPRIMAL*            primal              /**< primal */
   )
{
   assert(primal);
   assert(level >= 0);

   primal->decomplevel = level;
}


/** sets new primal bound if better than old one */
void reduce_primalUpdateUpperBound(
   SCIP_Real             ubnew,              /**< new upper bound, not including offset! */
   REDPRIMAL*            primal              /**< primal */
   )
{
   SCIP_Real ubnew_scaled;

   assert(primal);

   ubnew_scaled = ubnew + primal->offset;

   assert(GE(ubnew_scaled, 0.0));
   assert(LE(ubnew_scaled, FARAWAY));

   if( ubnew_scaled < primal->primalbound )
   {
      primal->primalbound = ubnew_scaled;
   }
}


/** gets */
SCIP_Real reduce_primalGetUpperBound(
   const REDPRIMAL*      primal              /**< primal */
   )
{
   assert(primal);

   if( EQ(primal->primalbound, FARAWAY) )
      return FARAWAY;

   assert(GE(primal->primalbound - primal->offset, 0.0));

   SCIPdebugMessage("returning best bound: %f \n", primal->primalbound - primal->offset);

   return (primal->primalbound - primal->offset);
}


/** gets */
SCIP_Real reduce_primalGetOffset(
   const REDPRIMAL*      primal              /**< primal */
   )
{
   assert(primal);


   return (primal->offset);
}

/** gets */
SCIP_Real* reduce_primalGetOffsetPointer(
   const REDPRIMAL*      primal              /**< primal */
   )
{
   assert(primal);


   return &(primal->offset);
}
