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
#include "scip/scip_cons.h"
#include "scip/scip_mem.h"

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

/** set labels for an array of constraints */
SCIP_RETCODE SCIPdecompSetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   )
{
   int i;

   assert(decomp != NULL);
   assert(conss != NULL);
   assert(labels != NULL);

   /* store each label in hash map */
   for( i = 0; i < nconss; ++i )
   {
      /* author bzfhende
       *
       * TODO ensure that labels are okay, e.g., nonnegative integers
       */

      SCIP_CALL( SCIPhashmapInsertInt(decomp->cons2block, (void *)conss[i], labels[i]) );
   }

   return SCIP_OKAY;
}

/** query labels for an array of constraints */
void SCIPdecompGetConsLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int*                  labels,             /**< array of labels, one per constraint */
   int                   nconss              /**< length of constraints array */
   )
{
   int i;

   assert(decomp != NULL);
   assert(conss != NULL);
   assert(labels != NULL);

   /* store variable labels in buffer array */
   for( i = 0; i < nconss; ++i )
   {
      if( SCIPhashmapExists(decomp->cons2block, (void *)conss[i]) )
      {
         labels[i] = SCIPhashmapGetImageInt(decomp->cons2block, (void *)conss[i]);
         printf("Returning %d\n", labels[i]);
      }
      else
         labels[i] = SCIP_DECOMP_LINKCONS;
   }
}

/** clears the corresponding labeling (constraints, variables, or both) of this decomposition */
SCIP_RETCODE SCIPdecompClear(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             clearvarlabels,     /**< should the variable labels be cleared? */
   SCIP_Bool             clearconslabels     /**< should the constraint labels be cleared? */
   )
{
   assert(decomp != NULL);

   if( clearvarlabels )
   {
      SCIP_CALL( SCIPhashmapRemoveAll(decomp->var2block) );
   }

   if( clearconslabels )
   {
      SCIP_CALL( SCIPhashmapRemoveAll(decomp->cons2block) );
   }

   return SCIP_OKAY;
}

/** raises an error if the condition is not TRUE */
static
SCIP_RETCODE ensureCondition(
   SCIP_Bool             condition           /**< some condition that must hold */
   )
{
   return condition ? SCIP_OKAY : SCIP_ERROR;
}

#define LABEL_UNASSIGNED INT_MIN
/** computes constraint labels from variable labels. Existing labels for the constraints are simply overridden */
SCIP_RETCODE SCIPdecompComputeConsLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   )
{
   int c;

   assert(decomp != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   /* assign label to each individual constraint */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_VAR** consvarsbuff;
      int nconsvars;
      int v;
      int* varlabels;
      int nlinkingvars;
      int conslabel;

      SCIP_Bool success;

      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, &success) );

      SCIP_CALL( ensureCondition(success) );

      SCIP_CALL( SCIPallocBufferArray(scip, &consvarsbuff, nconsvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, nconsvars) );


      SCIP_CALL( SCIPgetConsVars(scip, conss[c], consvarsbuff, nconsvars, &success) );

      SCIP_CALL( ensureCondition(success) );

      SCIPdecompGetVarsLabels(decomp, consvarsbuff, varlabels, nconsvars);

      /* loop over variable labels to compute the constraint label */
      conslabel = LABEL_UNASSIGNED;
      for( v = 0; v < nconsvars; ++v )
      {
         int varlabel = varlabels[v];

         if( varlabels[v] == SCIP_DECOMP_LINKVAR )
            ++nlinkingvars;
         else if( conslabel == LABEL_UNASSIGNED )
            conslabel = varlabel;
         else if( conslabel != varlabel )
         {
            conslabel = SCIP_DECOMP_LINKCONS;
            break;
         }
      }

      assert(nlinkingvars == nconsvars || conslabel != LABEL_UNASSIGNED);
      assert(v == nconsvars || conslabel == SCIP_DECOMP_LINKCONS);

      if( conslabel == LABEL_UNASSIGNED )
         conslabel = SCIP_DECOMP_LINKCONS;

      SCIP_CALL( SCIPdecompSetConsLabels(decomp, &(conss[c]), &conslabel, 1) );

      SCIPfreeBufferArray(scip, &varlabels);
      SCIPfreeBufferArray(scip, &consvarsbuff);
   }

   return SCIP_OKAY;
}
