/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   decomp.c
 * @brief  methods for working with decompositions
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "scip/struct_decomp.h"
#include "scip/decomp.h"
#include "scip/mem.h"
#include "scip/pub_misc.h"

/* author bzfhende
 *
 * TODO create and free a decomposition
 */

#define INIT_MAP_SIZE 2000

/** create a decomposition */
SCIP_RETCODE SCIPdecompCreate(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(decomp != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, decomp) );
   SCIP_CALL( SCIPhashmapCreate(&(*decomp)->var2block, blkmem, INIT_MAP_SIZE) );
   SCIP_CALL( SCIPhashmapCreate(&(*decomp)->cons2block, blkmem, INIT_MAP_SIZE) );

   (*decomp)->nblocks = 0;
   (*decomp)->score = -1.0;
   (*decomp)->haschanges = FALSE;

   return SCIP_OKAY;
}

/** free a decomposition */
void SCIPdecompFree(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   assert(decomp != NULL);
   assert(blkmem != NULL);

   if( *decomp == NULL )
      return;

   assert((*decomp)->var2block != NULL);
   SCIPhashmapFree(&(*decomp)->var2block);
   SCIPhashmapFree(&(*decomp)->cons2block);

   BMSfreeBlockMemory(blkmem, decomp);
}

/* author bzfhende
 *
 * TODO getter and setter for variable labelling
 */

/** set labels for an array of variables */
SCIP_RETCODE SCIPdecompSetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< array of labels, one per variable */
   int                   nvars               /**< length of variables array */
   )
{
   int i;

   assert(decomp != NULL);
   assert(vars != NULL);
   assert(labels != NULL);

   /* store each label in hash map */
   for( i = 0; i < nvars; ++i )
   {
      /* author bzfhende
       *
       * TODO ensure that labels are okay, e.g., nonnegative integers
       */

      SCIP_CALL( SCIPhashmapInsertInt(decomp->var2block, (void *)vars[i], labels[i]) );
   }

   return SCIP_OKAY;
}

/** query labels for an array of variables */
void SCIPdecompGetVarsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_VAR**            vars,               /**< array of variables */
   int*                  labels,             /**< buffer to store labels, one per variable */
   int                   nvars               /**< length of variables array */
   )
{
   int i;

   assert(decomp != NULL);
   assert(vars != NULL);
   assert(labels != NULL);

   /* store variable labels in buffer array */
   for( i = 0; i < nvars; ++i )
   {
      if( SCIPhashmapExists(decomp->var2block, (void *)vars[i]) )
         labels[i] = SCIPhashmapGetImageInt(decomp->var2block, (void *)vars[i]);
      else
         labels[i] = SCIP_DECOMP_LINKVAR;

   }
}
