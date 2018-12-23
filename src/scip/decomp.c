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
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_prob.h"
#include "scip/scip_var.h"
#include "scip/scip_mem.h"
#include "scip/struct_scip.h"
#include "scip/pub_cons.h"
#include "scip/scip_general.h"
#include "scip/scip_var.h"

/* author bzfhende
 *
 * TODO create and free a decomposition
 */

#define INIT_MAP_SIZE 2000

/** create a decomposition */
SCIP_RETCODE SCIPdecompCreate(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_Bool             original            /**< is this a decomposition in the original (TRUE) or transformed space? */
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
   (*decomp)->original = original;

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
      SCIP_CALL( SCIPhashmapSetImageInt(decomp->var2block, (void *)vars[i], labels[i]) );

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

      SCIP_CALL( SCIPhashmapSetImageInt(decomp->cons2block, (void *)conss[i], labels[i]) );
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

/** returns TRUE if decomposition is in the original space */
SCIP_Bool SCIPdecompIsOriginal(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->original;
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
   SCIP_VAR** varbuffer;
   int c;
   int twicenvars;
   int* varlabels;
   int* conslabels;

   assert(decomp != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   twicenvars = 2 * SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &varbuffer, twicenvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, twicenvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );


   /* assign label to each individual constraint */
   for( c = 0; c < nconss; ++c )
   {
      int nconsvars;
      int v;
      int nlinkingvars;
      int conslabel;
      int requiredsize;

      SCIP_Bool success;

      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, &success) );
      SCIP_CALL( ensureCondition(success) );
      SCIP_CALL( SCIPgetConsVars(scip, conss[c], varbuffer, twicenvars, &success) );
      SCIP_CALL( ensureCondition(success) );

      if( ! SCIPdecompIsOriginal(decomp) )
      {
         SCIP_CALL( SCIPgetActiveVars(scip, varbuffer, &nconsvars, twicenvars, &requiredsize) );
         assert(nconsvars <= twicenvars);
      }

      SCIPdecompGetVarsLabels(decomp, varbuffer, varlabels, nconsvars);

      /* loop over variable labels to compute the constraint label */
      conslabel = LABEL_UNASSIGNED;
      for( v = 0; v < nconsvars; ++v )
      {
         int varlabel = varlabels[v];

         /* count the number of linking variables, and keep track if there are two variables with different labels */
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

      /* if there are only linking variables, the constraint is unassigned */
      if( conslabel == LABEL_UNASSIGNED )
         conslabel = SCIP_DECOMP_LINKCONS;
      conslabels[c] = conslabel;
   }

   SCIP_CALL( SCIPdecompSetConsLabels(decomp, conss, conslabels, nconss) );

   SCIPfreeBufferArray(scip, &conslabels);
   SCIPfreeBufferArray(scip, &varlabels);
   SCIPfreeBufferArray(scip, &varbuffer);

   return SCIP_OKAY;
}

/** create a decomposition of the variables from a labeling of the constraints */
SCIP_RETCODE SCIPdecompComputeVarsLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of constraints */
   int                   nconss              /**< number of constraints */
   )
{
   int c;
   int* conslabels;
   SCIP_VAR** varbuffer;
   int twicenvars;


   assert(scip != NULL);
   assert(decomp != NULL);
   assert(conss != NULL);

   /* make the buffer array larger than necessary */
   twicenvars = 2 * SCIPgetNOrigVars(scip);

   /* allocate buffer to store constraint variables and labels */
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varbuffer, twicenvars) );

   /* query constraint labels */
   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

   /* iterate over constraints and query the corresponding constraint labels */
   for( c = 0; c < nconss; ++c )
   {
      int conslabel;
      int v;
      int nconsvars;
      SCIP_Bool success;
      int requiredsize;

      /* skip linking constraints */
      conslabel = conslabels[c];

      if( conslabel == SCIP_DECOMP_LINKCONS )
         continue;

      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, &success) );
      SCIP_CALL( ensureCondition(success) );
      SCIP_CALL( SCIPgetConsVars(scip, conss[c], varbuffer, twicenvars, &success) );
      SCIP_CALL( ensureCondition(success) );

      if( ! SCIPdecompIsOriginal(decomp) )
      {
         SCIP_CALL( SCIPgetActiveVars(scip, varbuffer, &nconsvars, twicenvars, &requiredsize) );
         assert(nconsvars <= twicenvars);
      }

      /** each variable in this constraint gets the constraint label unless it already has a different label -> make it a linking variable */
      for( v = 0; v < nconsvars; ++v )
      {
         SCIP_VAR* var = varbuffer[v];

         assert(SCIPvarIsActive(var));

         if( SCIPhashmapExists(decomp->var2block, (void *)var) )
         {
            int varlabel = SCIPhashmapGetImageInt(decomp->var2block, (void *)var);

            /* store the label linking variable explicitly to distinguish it from the default */
            if( varlabel != SCIP_DECOMP_LINKVAR && varlabel != conslabel )
               SCIP_CALL( SCIPhashmapSetImageInt(decomp->var2block, (void *)var, SCIP_DECOMP_LINKVAR) );
         }
         else
         {
            SCIP_CALL( SCIPhashmapInsertInt(decomp->var2block, (void *)var, conslabel) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &varbuffer);
   SCIPfreeBufferArray(scip, &conslabels);

   return SCIP_OKAY;
}

/** create a decomposition storage */
SCIP_RETCODE SCIPdecompstoreCreate(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   int                   nslots              /**< maximum number of decomposition slots in storage */
   )
{
   assert(decompstore != NULL);
   assert(blkmem != NULL);
   assert(nslots > 0);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, decompstore) );

   (*decompstore)->ndecomps = 0;
   (*decompstore)->norigdecomps = 0;
   (*decompstore)->decompssize = nslots;

   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decompstore)->decomps, nslots) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decompstore)->origdecomps, nslots) );

   return SCIP_OKAY;
}

/** free array of decompositions */
static
void freeDecompositions(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_DECOMP**         decomps,            /**< decomposition array */
   int*                  ndecomps            /**< pointer for initial number of decompositions, will be set to 0 */
   )
{
   int d;

   assert(decomps != NULL);
   assert(ndecomps != NULL);

   /* delete all remaining decompositions from this store */
   for( d = 0; d < *ndecomps; ++d )
      SCIPdecompFree(&decomps[d], blkmem);

   *ndecomps = 0;
}

/** free all decompositions in transformed space */
void SCIPexitSolveDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_DECOMPSTORE* decompstore = scip->decompstore;

   assert(decompstore != NULL);

   freeDecompositions(SCIPblkmem(scip), decompstore->decomps, &decompstore->ndecomps);
}

/** free a decomposition storage */
void SCIPdecompstoreFree(
   SCIP_DECOMPSTORE**    decompstore,        /**< pointer to store decomposition storage */
   BMS_BLKMEM*           blkmem              /**< block memory data structure */
   )
{
   assert(decompstore != NULL);

   if( *decompstore == NULL )
      return;

   freeDecompositions(blkmem, (*decompstore)->origdecomps, &(*decompstore)->norigdecomps);
   freeDecompositions(blkmem, (*decompstore)->decomps, &(*decompstore)->ndecomps);


   BMSfreeBlockMemoryArray(blkmem, &(*decompstore)->decomps, (*decompstore)->decompssize);
   BMSfreeBlockMemoryArray(blkmem, &(*decompstore)->origdecomps, (*decompstore)->decompssize);

   BMSfreeBlockMemory(blkmem, decompstore);
}

/** add decomposition to storage */
SCIP_RETCODE SCIPdecompstoreAdd(
   SCIP_DECOMPSTORE*     decompstore,        /**< decomposition storage */
   SCIP_DECOMP*          decomp              /**< decomposition to add */
   )
{
   SCIP_DECOMP** decomps;
   int* ndecompsptr;

   assert(decompstore != NULL);
   assert(decomp != NULL);

   /* distinguish between storage for original or transformed decompositions */
   if( SCIPdecompIsOriginal(decomp) )
   {
      decomps = decompstore->origdecomps;
      ndecompsptr = &decompstore->norigdecomps;
   }
   else
   {
      decomps = decompstore->decomps;
      ndecompsptr = &decompstore->ndecomps;
   }

   /* ensure that storage capacity is not exceeded */
   if( *ndecompsptr == decompstore->decompssize )
   {
      SCIPerrorMessage("Error: Decomposition storage size exceeded, maximum is %d decompositions\n", decompstore->decompssize);
      return SCIP_ERROR;
   }

   decomps[(*ndecompsptr)++] = decomp;

   return SCIP_OKAY;
}

/** get decompositions from this storage */
SCIP_DECOMP** SCIPdecompstoreGetDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   )
{
   assert(decompstore != NULL);

   return decompstore->decomps;
}

/** get number of decompositions in this storage */
int SCIPdecompstoreGetNDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   )
{
   assert(decompstore != NULL);
   return decompstore->ndecomps;
}

/** get decompositions from this storage */
SCIP_DECOMP** SCIPdecompstoreGetOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   )
{
   assert(decompstore != NULL);

   return decompstore->origdecomps;
}

/** get number of decompositions in this storage */
int SCIPdecompstoreGetNOrigDecomps(
   SCIP_DECOMPSTORE*     decompstore         /**< decomposition storage */
   )
{
   assert(decompstore != NULL);
   return decompstore->norigdecomps;
}

/** get decomposition store from SCIP */
SCIP_DECOMPSTORE* SCIPgetDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   return scip->decompstore;
}

/** transform all available original decompositions into transformed space */
SCIP_RETCODE SCIPtransformDecompstore(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   int d;
   int v;
   SCIP_DECOMPSTORE* decompstore;
   SCIP_VAR** vars;
   SCIP_VAR** origvars;
   SCIP_VAR** varssorted;
   SCIP_CONS** conss;
   int nconss;
   int nvars;
   int nvarsoriginal;
   int nvarsintroduced;
   int* varslabels;
   SCIP_Bool original = FALSE;



   assert(scip != NULL);
   assert(scip->decompstore != NULL);

   decompstore = scip->decompstore;
   assert(decompstore->ndecomps == 0);

   assert(SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED);

   nvars = SCIPgetNVars(scip);
   vars = SCIPgetVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &varssorted, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &origvars, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );


   /* determine if variable has an original counterpart or not, and put it into varssorted array at the front or back */
   nvarsoriginal = nvarsintroduced = 0;
   for( v = 0; v < nvars; ++v )
   {
      SCIP_Real scalar;
      SCIP_Real constant;
      SCIP_VAR* origvar;

      origvar = vars[v];
      scalar = 1.0;
      constant = 0.0;
      SCIP_CALL( SCIPvarGetOrigvarSum(&origvar, &scalar, &constant) );

      /* the variable has no original counterpart and is therefore put at the end of the array */
      if( origvar == NULL )
      {
         varssorted[nvars - 1 - nvarsintroduced] = vars[v];
         ++nvarsintroduced;
      }
      else
      {
         varssorted[nvarsoriginal] = vars[v];
         origvars[nvarsoriginal] = origvar;
         ++nvarsoriginal;
      }

      assert(nvarsoriginal + nvarsintroduced <= nvars);
   }

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   /* loop over available, original decompositions, transform and add them to the storage */
   for( d = 0; d < decompstore->norigdecomps; ++d )
   {
      SCIP_DECOMP* origdecomp = decompstore->origdecomps[d];
      SCIP_DECOMP* decomp;

      /* 1. query the decomposition labels of the original variables and set them for the transformed variables
       * that have original counterparts
       */

      SCIP_CALL( SCIPdecompCreate(&decomp, SCIPblkmem(scip), original) );

      SCIPdecompGetVarsLabels(origdecomp, origvars, varslabels, nvarsoriginal);

      SCIP_CALL( SCIPdecompSetVarsLabels(decomp, varssorted, varslabels, nvarsoriginal) );

      /* 2. compute the constraint labels based on the preliminary variable labels */
      SCIP_CALL( SCIPdecompComputeConsLabels(scip, decomp, conss, nconss) );

      /* 3. remove the variable labels now that we have constraint labels */
      SCIP_CALL( SCIPdecompClear(decomp, TRUE, FALSE) );

      /* 4. use the constraint labels for the final variable labeling */
      SCIP_CALL( SCIPdecompComputeVarsLabels(scip, decomp, conss, nconss) );

      SCIP_CALL( SCIPdecompstoreAdd(decompstore, decomp) );
   }

   SCIPfreeBufferArray(scip, &varslabels);
   SCIPfreeBufferArray(scip, &origvars);
   SCIPfreeBufferArray(scip, &varssorted);

   return SCIP_OKAY;
}
