/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip.c
 * @brief  SCIP callable library
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "set.h"
#include "mem.h"
#include "prob.h"
#include "tree.h"
#include "lp.h"
#include "stat.h"
#include "scip.h"


/** SCIP main data structure */
struct Scip
{
   MEM*             mem;                /**< block memory buffers */
   SET*             set;                /**< global SCIP settings */
   PROB*            prob;               /**< original problem data */
   TREE*            tree;               /**< branch and bound tree */
   LP*              lp;                 /**< LP data */
   STAT*            stat;               /**< dynamic problem statistics */
};




RETCODE SCIPcreate(                     /**< creates and initializes SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);

   ALLOC_OKAY( allocMemory(*scip) );

   CHECK_OKAY( SCIPmemCreate(&(*scip)->mem) );
   CHECK_OKAY( SCIPsetCreate(&(*scip)->set) );
   CHECK_OKAY( SCIPprobCreate(&(*scip)->prob) );
   CHECK_OKAY( SCIPtreeCreate(&(*scip)->tree, (*scip)->set) );
   CHECK_OKAY( SCIPlpCreate(&(*scip)->lp) );
   CHECK_OKAY( SCIPstatCreate(&(*scip)->stat) );

   return SCIP_OKAY;
}

RETCODE SCIPfree(                       /**< frees SCIP data structures */
   SCIP**           scip                /**< pointer to SCIP data structure */
   )
{
   assert(scip != NULL);
   assert(*scip != NULL);

   CHECK_OKAY( SCIPstatFree(&(*scip)->stat) );
   CHECK_OKAY( SCIPlpFree(&(*scip)->lp) );
   CHECK_OKAY( SCIPtreeCreate(&(*scip)->tree, (*scip)->set) );
   CHECK_OKAY( SCIPprobCreate(&(*scip)->prob) );
   CHECK_OKAY( SCIPsetCreate(&(*scip)->set) );
   CHECK_OKAY( SCIPmemCreate(&(*scip)->mem) );

   ALLOC_OKAY( allocMemory(*scip) );

   return SCIP_OKAY;
}

