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
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_SOLPOOL_H_
#define APPLICATIONS_STP_SRC_SOLPOOL_H_

#include "scip/scip.h"
#include "graph.h"


/** element of Steiner tree solution pool */
typedef struct stp_solution
{
   SCIP_Real obj;
   int* soledges;
   int index;
} STPSOL;

/** edge based solution pool for Steiner tree problems (in presolving) */
typedef struct stp_solution_pool
{
   STPSOL** sols;
   int size;
   int nedges;
   int maxsize;
   int maxindex;
} STPSOLPOOL;


/** is given solution in pool? */
SCIP_Bool solpool_isContained(
   const int*            soledges,           /**< edge array of solution to be checked */
   const STPSOLPOOL*     pool                /**< the pool */
);

/** get solution from index */
STPSOL* solpool_solFromIndex(
   STPSOLPOOL*           pool,               /**< the pool */
   const int             soindex             /**< the index */
    );

/** initializes STPSOL pool */
SCIP_RETCODE solpool_init(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL**          pool,               /**< the pool */
   const int             nedges,             /**< number of edges of solutions to be stored in the pool */
   const int             maxsize             /**< capacity of pool */
   );

/** tries to add STPSOL to pool */
SCIP_RETCODE solpool_addSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const SCIP_Real       obj,                /**< objective of solution to be added */
   const int*            soledges,           /**< edge array of solution to be added */
   STPSOLPOOL*           pool,               /**< the pool */
   SCIP_Bool*            success             /**< has solution been added? */
   );

/** frees STPSOL pool */
void solpool_free(
   SCIP*                 scip,               /**< SCIP data structure */
   STPSOLPOOL**          pool                /**< the pool */
   );


#endif /* APPLICATIONS_STP_SRC_STPSOLPOOL_H_ */
