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
#include "scip/scip_datastructures.h"

/* create and free a decomposition */
#define INIT_MAP_SIZE 2000

/** create a decomposition */
SCIP_RETCODE SCIPdecompCreate(
   SCIP_DECOMP**         decomp,             /**< pointer to store the decomposition data structure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   int                   nblocks,            /**< the number of blocks (without the linking block) */
   SCIP_Bool             original,           /**< is this a decomposition in the original (TRUE) or transformed space? */
   SCIP_Bool             benderslabels       /**< should the variables be labeled for the application of Benders' decomposition */
   )
{
   assert(decomp != NULL);
   assert(blkmem != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, decomp) );
   SCIP_CALL( SCIPhashmapCreate(&(*decomp)->var2block, blkmem, INIT_MAP_SIZE) );
   SCIP_CALL( SCIPhashmapCreate(&(*decomp)->cons2block, blkmem, INIT_MAP_SIZE) );

   /* we allocate one extra slot for the linking block */
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decomp)->varssize, nblocks + 1) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decomp)->consssize, nblocks + 1) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*decomp)->labels, nblocks + 1) );

   (*decomp)->nblocks = nblocks;
   (*decomp)->score = -1.0;
   (*decomp)->modularity = -1.0;
   (*decomp)->idxsmallestblock = 0;
   (*decomp)->idxlargestblock = 0;
   (*decomp)->haschanges = FALSE;
   (*decomp)->original = original;
   (*decomp)->benderslabels = benderslabels;
   (*decomp)->areascore = -1.0;
   (*decomp)->nedges = 0;
   (*decomp)->mindegree = 0;
   (*decomp)->maxdegree = 0;
   (*decomp)->ncomponents = 0;

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

   BMSfreeBlockMemoryArray(blkmem, &(*decomp)->varssize, (*decomp)->nblocks + 1);
   BMSfreeBlockMemoryArray(blkmem, &(*decomp)->consssize, (*decomp)->nblocks + 1);
   BMSfreeBlockMemoryArray(blkmem, &(*decomp)->labels, (*decomp)->nblocks + 1);

   BMSfreeBlockMemory(blkmem, decomp);
}

/* getter and setter for variable labels */

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

/** sets the parameter that indicates whether the variables must be labeled for the application of Benders'
 * decomposition
 */
void SCIPdecompSetUseBendersLabels(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Bool             benderslabels       /**< whether Benders' variable labels should be used */
   )
{
   assert(decomp != NULL);

   decomp->benderslabels = benderslabels;
}

/** returns TRUE if the variables must be labeled for the application of Benders' decomposition */
SCIP_Bool SCIPdecompUseBendersLabels(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->benderslabels;
}

/** gets number of blocks of this decomposition */
int SCIPdecompGetNBlocks(
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   assert(decomp != NULL);

   return decomp->nblocks;
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
/** computes constraint labels from variable labels.
 *
 *  Existing labels for the constraints are simply overridden
 *
 *  The computed labels depend on the flag SCIPdecompUseBendersLabels() of the decomposition.
 *
 *  If the flag is set to FALSE, the labeling assigns
 *
 *  - label i, if only variables labeled i are present in the constraint (and optionally linking variables)
 *  - SCIP_DECOMP_LINKCONS, if there are either only variables labeled with SCIP_DECOMP_LINKVAR present, or
 *    if there are variables with more than one block label.
 *
 *  If the flag is set to TRUE, the assignment is the same, unless variables from 2 named blocks occur in the same
 *  constraint, which is an invalid labeling for the Benders case.
 *   */
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
   SCIP_Bool benderserror;
   SCIP_Bool benderslabels;

   assert(decomp != NULL);

   if( nconss == 0 )
      return SCIP_OKAY;

   twicenvars = 2 * SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &varbuffer, twicenvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, twicenvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );

   benderslabels = SCIPdecompUseBendersLabels(decomp);
   benderserror = FALSE;

   /* assign label to each individual constraint */
   for( c = 0; c < nconss && ! benderserror; ++c )
   {
      int nconsvars;
      int v;
      int nlinkingvars = 0;
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
      else
      {
         for( v = 0; v < nconsvars; ++v )
         {
            assert(SCIPvarIsActive(varbuffer[v]) || SCIPvarIsNegated(varbuffer[v]));

         /* some constraint handlers such as indicator may already return inactive variables */
         if( SCIPvarIsNegated(varbuffer[v]) )
            varbuffer[v] = SCIPvarGetNegatedVar(varbuffer[v]);
         }
      }

      SCIPdecompGetVarsLabels(decomp, varbuffer, varlabels, nconsvars);

      /* loop over variable labels to compute the constraint label */
      conslabel = LABEL_UNASSIGNED;
      for( v = 0; v < nconsvars; ++v )
      {
         int varlabel = varlabels[v];

         /* count the number of linking variables, and keep track if there are two variables with different labels */
         if( varlabel == SCIP_DECOMP_LINKVAR )
            ++nlinkingvars;
         else if( conslabel == LABEL_UNASSIGNED )
            conslabel = varlabel;
         else if( conslabel != varlabel )
         {
            /* there must not be two variables from different named blocks in a single constraint, since the presence
             * of named block variables forbids this constraint from the master (linking) block
             */
            if( benderslabels )
               benderserror = TRUE;

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

   /* throw an error and inform the user if the variable block decomposition does not allow a benders constraint labeling */
   if( benderserror )
   {
      SCIPerrorMessage("Error in constraint label computation; variables from multiple named blocks in a single constraint\n");

      return SCIP_INVALIDDATA;
   }


   return SCIP_OKAY;
}

/** create a decomposition of the variables from a labeling of the constraints.
 *
 *  NOTE: by default, the variable labeling is based on a Dantzig-Wolfe decomposition. This means that constraints in named
 *  blocks have have precedence over linking constraints. If a variable exists in constraints from
 *  two or more named blocks, then this variable is marked as a linking variable.
 *  If a variable occurs in exactly one named block i>=0, it is assigned label i.
 *  Variables which are only in linking constraints are unlabeled. However, SCIPdecompGetVarsLabels() will
 *  label them as linking variables.
 *
 *  If the variables should be labeled for the application of Benders' decomposition, the decomposition must be
 *  flagged explicitly via SCIPdecompSetUseBendersLabels().
 *  With this setting, the presence in linking constraints takes precedence over the presence in named blocks.
 *  Now, a variable is considered linking if it is present in at least one linking constraint and an arbitrary
 *  number of constraints from named blocks.
 */
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
   SCIP_Bool benderslabels;

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

   benderslabels = SCIPdecompUseBendersLabels(decomp);
   /* iterate over constraints and query the corresponding constraint labels */
   for( c = 0; c < nconss; ++c )
   {
      int conslabel;
      int v;
      int nconsvars;
      SCIP_Bool success;
      int requiredsize;
      int newvarlabel;

      conslabel = conslabels[c];

      if( conslabel == SCIP_DECOMP_LINKCONS )
      {
         /* skip linking constraints unless Benders labeling is used */
         if( ! benderslabels )
            continue;
         else
            newvarlabel = SCIP_DECOMP_LINKVAR;
      }
      else
         newvarlabel = conslabel;

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

         assert(SCIPvarIsActive(var) || (SCIPdecompIsOriginal(decomp) && SCIPvarIsNegated(var)));

         /* some constraint handlers such as indicator may already return inactive variables */
         if( SCIPvarIsNegated(var) )
            var = SCIPvarGetNegatedVar(var);

         if( SCIPhashmapExists(decomp->var2block, (void *)var) )
         {
            int varlabel = SCIPhashmapGetImageInt(decomp->var2block, (void *)var);

            /* store the label linking variable explicitly to distinguish it from the default */
            if( varlabel != SCIP_DECOMP_LINKVAR && varlabel != newvarlabel )
               SCIP_CALL( SCIPhashmapSetImageInt(decomp->var2block, (void *)var, SCIP_DECOMP_LINKVAR) );
         }
         else
         {
            SCIP_CALL( SCIPhashmapInsertInt(decomp->var2block, (void *)var, newvarlabel) );
         }
      }
   }

   SCIPfreeBufferArray(scip, &varbuffer);
   SCIPfreeBufferArray(scip, &conslabels);

   return SCIP_OKAY;
}

/** count occurrences of label in array, starting from pos */
static
int countLabelFromPos(
   int*                  labels,             /**< array of labels */
   int                   pos,                /**< position to start counting from */
   int                   nlabels             /**< the number of labels */
   )
{
   int endpos = pos;
   int currlabel;

   assert(labels != NULL);
   assert(pos < nlabels);

   currlabel = labels[pos];

   do
   {
      endpos++;
   }
   while( endpos < nlabels && labels[endpos] == currlabel );

   return endpos - pos;
}

/** compute decomposition modularity */
static
SCIP_RETCODE computeModularity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_Real*            modularity          /**< pointer to store modularity value */
   )
{
   SCIP_CONS** conss;
   SCIP_VAR** varbuf;
   int* varslabels;
   int* conslabels;
   int* totaldegrees; /* the total degree for every block */
   int* withinedges; /* the number of edges within each block */
   int nnonzeroes = 0;
   int nvars;
   int nconss;
   int c;
   int b;

   /* allocate buffer arrays to hold constraint and variable labels,  and store within-edges and total community degrees */
   nvars = SCIPgetNVars(scip);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varbuf, nvars) );

   SCIP_CALL( SCIPallocClearBufferArray(scip, &totaldegrees, decomp->nblocks + 1) );
   SCIP_CALL( SCIPallocClearBufferArray(scip, &withinedges, decomp->nblocks + 1) );

   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

   /*
    * loop over all nonzeroes, consider the labels of the incident nodes (cons and variable)
    * and increase the corresponding counters
    */
   for( c = 0; c < nconss; ++c )
   {
      int nconsvars;
      int conslabel = conslabels[c];
      int blockpos;
      int varblockstart;
      SCIP_Bool success;
      SCIP_Bool found;

      /* linking constraints do not contribute to the modularity */
      if( conslabel == SCIP_DECOMP_LINKCONS )
         continue;

      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, &success) );
      SCIP_CALL( ensureCondition(success) );

      if( nvars == 0 )
         continue;

      SCIP_CALL( SCIPgetConsVars(scip, conss[c], varbuf, nvars, &success) );
      SCIP_CALL( ensureCondition(success) );

      /* author bzfhende
       *
       * TODO replace the variables by their active representatives for transformed decompositions
       */
      SCIPdecompGetVarsLabels(decomp, varbuf, varslabels, nconsvars);

      /* find the position of the constraint label. Constraints of the border always belong to the first block at index 0 */
      found = SCIPsortedvecFindInt(decomp->labels, conslabel, decomp->nblocks + 1, &blockpos);
      assert(found);

      SCIPsortInt(varslabels, nconsvars);

      /* count occurences of labels (blocks) in the sorted labels array */
      varblockstart = 0;
      while( varblockstart < nconsvars )
      {
         int varblockpos;
         int nblockvars = countLabelFromPos(varslabels, varblockstart, nconsvars);

         found = SCIPsortedvecFindInt(decomp->labels, varslabels[varblockstart], decomp->nblocks + 1, &varblockpos);
         assert(found);

         /* don't consider linking variables for modularity statistics */
         if( varslabels[varblockstart] != SCIP_DECOMP_LINKVAR )
         {
            /* increase the number of within edges for variable and constraints from the same block */
            if( varblockpos == blockpos )
               withinedges[varblockpos] += nblockvars;

            /* increase the total degrees and nonzero (edge) counts; it is intended that the total degrees sum up
             * to twice the number of edges
             */
            totaldegrees[blockpos] += nblockvars;
            totaldegrees[varblockpos] += nblockvars;
            nnonzeroes += nblockvars;
         }

         varblockstart += nblockvars;
      }
   }

/* ensure that total degrees sum up to twice the number of edges */
#ifndef NDEBUG
   {
      int totaldegreesum = 0;
      for( b = 1; b < decomp->nblocks + 1; ++b )
         totaldegreesum += totaldegrees[b];

      assert(totaldegreesum == 2 * nnonzeroes);
   }
#endif

   /* compute modularity */
   *modularity = 0.0;
   for( b = 1; b < decomp->nblocks + 1; ++b )
   {
      SCIP_Real expectedval;
      expectedval = totaldegrees[b] / (2.0 * nnonzeroes);
      expectedval = SQR(expectedval);
      *modularity += (withinedges[b] / (SCIP_Real)nnonzeroes) - expectedval;
   }

   SCIPfreeBufferArray(scip, &withinedges);
   SCIPfreeBufferArray(scip, &totaldegrees);
   SCIPfreeBufferArray(scip, &varbuf);
   SCIPfreeBufferArray(scip, &varslabels);
   SCIPfreeBufferArray(scip, &conslabels);

   return SCIP_OKAY;
}

/** compute the area score */
static
void computeAreaScore(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   SCIP_Real score = 1.0;
   int nvars = SCIPgetNVars(scip);
   int nconss = SCIPgetNConss(scip);

   if( nvars > 0 && nconss > 0 )
   {
      int nlinkconss = decomp->consssize[0];
      int nlinkvars = decomp->varssize[0];
      SCIP_Real factor = 1.0 / ((SCIP_Real)nvars * nconss);

      int i;

      /* compute diagonal contributions to the area score */
      for( i = 1; i < decomp->nblocks + 1; ++i )
      {
         score -= (factor * decomp->consssize[i]) * decomp->varssize[i];
      }

      score -= ((SCIP_Real)nlinkconss * nvars + (SCIP_Real)nconss * nlinkvars - (SCIP_Real)nlinkconss * nlinkvars) * factor;
   }

   decomp->areascore = score;
}

/** build the block decomposition graph */
static
void buildBlockGraph(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_VAR** consvars;
   SCIP_DIGRAPH* blocklinkingvargraph;
   SCIP_DIGRAPH* blockgraph;
   int* varlabels;
   int* conslabels;
   int* linkvaridx;
   int* succnodes;
   SCIP_Bool success, found;
   int nvars;
   int nconss;
   int nblocks;
   int nlinkingvars = 0;
   int nconsvars;
   int nsucc, succ1, succ2;
   int tempmin, tempmax;
   int i, j, v, n;

   assert(scip != NULL);
   assert(decomp != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   nblocks = SCIPdecompGetNBlocks(decomp);

   SCIPallocBufferArray(scip, &conslabels, nconss);
   SCIPallocBufferArray(scip, &varlabels, nvars);
   SCIPallocBufferArray(scip, &linkvaridx, nvars);

   /* store variable and constraint labels in buffer arrays */
   SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);
   SCIPdecompGetVarsLabels(decomp, vars, varlabels, nvars);

   /* create a mapping of all linking variables to 0,..,nlinkingvars -1 and store it in array linkvaridx */
   for( v = 0; v < nvars; ++v )
   {
      if( varlabels[v] == SCIP_DECOMP_LINKVAR )
      {
         linkvaridx[v] = nlinkingvars;
         assert(SCIPvarGetProbindex(vars[v]) == v);
         ++nlinkingvars;
      }
   }

   /* create a bi-partite graph composed of block and linking var nodes */
   SCIPcreateDigraph(scip, &blocklinkingvargraph, nblocks + nlinkingvars); /* nblocks does not include the linking constraints block */

   for( i = 0; i < nconss; ++i )
   {
      int requiredsize;
      /* linking constraints are skipped in this checking step */
      if( conslabels[i] == SCIP_DECOMP_LINKCONS )
         continue;

      SCIPgetConsNVars(scip, conss[i], &nconsvars, &success);
      if ( !success )
         continue;
      SCIPallocBufferArray(scip, &consvars, nconsvars);
      SCIPgetConsVars(scip, conss[i], consvars, nconsvars, &success);
      if( !success )
         continue;

      /* re-transform given variables to active variables*/
      SCIPgetActiveVars(scip, consvars, &nconsvars, nconsvars, &requiredsize);
      assert(requiredsize <= nconsvars);
      SCIPdecompGetVarsLabels(decomp, consvars, varlabels, nconsvars);

      /* adding doube-direction arcs between blocks and corresponding linkingvars */
      for( j = 0; j < nconsvars; ++j )
      {
         assert(consvars[j] != NULL);
         if( varlabels[j] == SCIP_DECOMP_LINKVAR )
         {
            int linkingvarnodeidx = linkvaridx[SCIPvarGetProbindex(consvars[j])];
            int blocknodeidx;

            /* find the position of the constraint label. Subtract by 1 to get the node index as the 1st block is reserved for linking constraints */
            found = SCIPsortedvecFindInt(decomp->labels, conslabels[i], decomp->nblocks + 1, &blocknodeidx); /* assuming labels is sorted */
            assert(found);

            SCIPdigraphAddArcSafe(blocklinkingvargraph, SCIPdecompGetNBlocks(decomp) + linkingvarnodeidx, blocknodeidx - 1, NULL);
            SCIPdigraphAddArcSafe(blocklinkingvargraph, blocknodeidx - 1, SCIPdecompGetNBlocks(decomp) + linkingvarnodeidx, NULL);
         }
      }
   }
   assert(SCIPdigraphGetNNodes(blocklinkingvargraph) > 0);

   /* From the information of the above bi-partite graph, build the block-decomposition graph: nodes -> blocks and double-direction arcs -> linking variables */
   SCIPcreateDigraph(scip, &blockgraph, nblocks);

   for( n = nblocks; n < SCIPdigraphGetNNodes(blocklinkingvargraph); ++n )
   {
      nsucc = (int) SCIPdigraphGetNSuccessors(blocklinkingvargraph, n);
      succnodes = (int*) SCIPdigraphGetSuccessors(blocklinkingvargraph, n);
      for( succ1 = 0; succ1 < nsucc; ++succ1 )
      {
         for( succ2 = 0; succ2 < nsucc; ++succ2 )
         {
            if( succnodes[succ1] != succnodes[succ2] ) /* no self-loops */
               SCIPdigraphAddArcSafe(blockgraph, succnodes[succ1], succnodes[succ2], NULL);
         }
      }
   }

   assert(SCIPdigraphGetNNodes(blockgraph) > 0);

   /* Get the number of edges in the block-decomposition graph */
   decomp->nedges = SCIPdigraphGetNArcs(blockgraph) / 2;

   /* Get the minimum and maximum degree of the block-decomposition graph */
   tempmin = (int) SCIPdigraphGetNSuccessors(blockgraph, 0);
   tempmax = (int) SCIPdigraphGetNSuccessors(blockgraph, 0);
   for( n = 1; n < SCIPdigraphGetNNodes(blockgraph); ++n )
   {
      nsucc = (int) SCIPdigraphGetNSuccessors(blockgraph, n);
      if( nsucc < tempmin )
         tempmin = nsucc;
      else
         if( nsucc > tempmax )
            tempmax = nsucc;
   }

   decomp->mindegree = tempmin;
   decomp->maxdegree = tempmax;

   /* Calculate the number of connected components in the block-decomposition graph */
   SCIPdigraphComputeUndirectedComponents(blockgraph, -1, NULL, NULL);
   decomp->ncomponents = SCIPdigraphGetNComponents(blockgraph);

   /* TODO: Calculate the number of articulation nodes in the block-decomposition graph using DFS*/

   SCIPfreeBufferArray(scip, &linkvaridx);
   SCIPfreeBufferArray(scip, &varlabels);
   SCIPfreeBufferArray(scip, &conslabels);
   SCIPdigraphFree(&blocklinkingvargraph);
   SCIPdigraphFree(&blockgraph);
}

/** compute decomposition statistics and store them in the decomp object */
SCIP_RETCODE SCIPcomputeDecompStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   int* varslabels;
   int* conslabels;
   int nvars;
   int nconss;
   int varblockstart;
   int consblockstart;
   int currlabelidx;
   int varidx;
   int considx;

   assert(scip != NULL);
   assert(decomp != NULL);

   /* store variable and constraint labels in buffer arrays */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

  /* return if problem is empty
   *
   * TODO ensure that statistics reflect this correctly
   */
  if( nvars == 0 || nconss == 0 )
  {
     return SCIP_OKAY;
  }

  SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );
  SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );

  SCIPdecompGetVarsLabels(decomp, vars, varslabels, nvars);
  SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

  /* sort both buffer arrays for quick counting */
  SCIPsortInt(varslabels, nvars);
  SCIPsortInt(conslabels, nconss);


  /* the first label is always LINKVAR, even if Benders' variable labels are used. We can ignore the variables
   * labelled as LINKCONS since this label is only required when computing the variable labels for Benders'
   * decomposition.
   */
  decomp->labels[0] = SCIP_DECOMP_LINKVAR;

  /* treating the linking variables first */
  if( varslabels[0] == SCIP_DECOMP_LINKVAR )
     decomp->varssize[0] = countLabelFromPos(varslabels, 0, nvars);
  else
     decomp->varssize[0] = 0;

  /* count border constraints and store their number */
  if( conslabels[0] == SCIP_DECOMP_LINKCONS )
     decomp->consssize[0] = countLabelFromPos(conslabels, 0, nconss);
  else
     decomp->consssize[0] = 0;


  /* merge labels (except for border at position 0) since neither variable nor constraint labels by themselves need to be complete */
  currlabelidx = 1;
  varidx = decomp->varssize[0];
  considx = decomp->consssize[0];

  while( varidx < nvars || considx < nconss )
  {
      int varlabel;
      int conslabel;

      varlabel = varidx < nvars ? varslabels[varidx] : INT_MAX;
      conslabel = considx < nconss ? conslabels[considx] : INT_MAX;

      /* store the smaller of the two current labels */
      decomp->labels[currlabelidx++] = MIN(varlabel, conslabel);

      /* increment the variable or constraint indices or both, depending on one label being strictly smaller */
      if( varlabel <= conslabel )
         varidx += countLabelFromPos(varslabels, varidx, nvars);

      if( conslabel <= varlabel )
         considx += countLabelFromPos(conslabels, considx, nconss);
  }


  currlabelidx = 1;
  varblockstart = decomp->varssize[0];
  assert(varblockstart == nvars || varslabels[varblockstart] > SCIP_DECOMP_LINKVAR);
   /* loop over the variables count the occurrences, storing also the integer labels */
  while( varblockstart < nvars )
  {
     assert(currlabelidx < decomp->nblocks + 1);
     assert(decomp->labels[currlabelidx] <= varslabels[varblockstart]);

     /* check if the current label is present in the variable labels */
     if( decomp->labels[currlabelidx] < varslabels[varblockstart] )
     {
        decomp->varssize[currlabelidx] = 0;
     }
     else
     {
        /* increment variable block start */
        decomp->varssize[currlabelidx] = countLabelFromPos(varslabels, varblockstart, nvars);
        varblockstart += decomp->varssize[currlabelidx];
     }
     currlabelidx++;
  }




  consblockstart = decomp->consssize[0];
  assert(consblockstart == nconss || conslabels[consblockstart] >= 0);

  /* loop over remaining, non-border constraints */
  decomp->idxsmallestblock = decomp->idxlargestblock = -1;
  currlabelidx = 1;
  while( consblockstart < nconss )
  {
     assert(currlabelidx < decomp->nblocks + 1);
     assert(decomp->labels[currlabelidx] <= conslabels[consblockstart]);

     /* the current label may not be present in the conslabels */
     if( decomp->labels[currlabelidx] < conslabels[consblockstart] )
     {
        decomp->consssize[currlabelidx] = 0;
     }
     else
     {
        /* count the number of occurrences and store it */
        decomp->consssize[currlabelidx] = countLabelFromPos(conslabels, consblockstart, nconss);
     }

     /* update index to largest and smallest constraint blocks (don't consider border for this statistic) */
     if( decomp->idxlargestblock == -1 )
     {
        decomp->idxlargestblock = decomp->idxsmallestblock = currlabelidx;
     }
     else if( decomp->consssize[currlabelidx] > decomp->consssize[decomp->idxlargestblock] )
        decomp->idxlargestblock = currlabelidx;
     else if( decomp->consssize[currlabelidx] < decomp->consssize[decomp->idxsmallestblock] )
        decomp->idxsmallestblock = currlabelidx;

     consblockstart += decomp->consssize[currlabelidx];
     currlabelidx++;
  }

  SCIP_CALL( computeModularity(scip, decomp, &decomp->modularity) );

  computeAreaScore(scip, decomp);
  buildBlockGraph(scip, decomp);

  SCIPfreeBufferArray(scip, &conslabels);
  SCIPfreeBufferArray(scip, &varslabels);

   return SCIP_OKAY;
}

/** print decomposition statistics into string buffer */
char* SCIPdecompPrintStats(
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   char*                 strbuf              /**< string buffer storage */
   )
{
   char* ptr;

   assert(decomp != NULL);
   assert(strbuf != NULL);

   ptr = strbuf;

   ptr += snprintf(ptr, SCIP_MAXSTRLEN,
            "Decomposition with %d blocks.\n",
            decomp->nblocks);
   ptr += snprintf(ptr, SCIP_MAXSTRLEN,
            "Largest block: Block %d with %d constraints and %d variables\n",
            decomp->nblocks == 0 ? -1 : decomp->labels[decomp->idxlargestblock],
            decomp->nblocks == 0 ? 0 : decomp->consssize[decomp->idxlargestblock],
            decomp->nblocks == 0 ? 0 : decomp->varssize[decomp->idxlargestblock]);
   ptr += snprintf(ptr, SCIP_MAXSTRLEN,
            "Smallest block: Block %d with %d constraints and %d variables\n",
            decomp->nblocks == 0 ? 0 : decomp->labels[decomp->idxsmallestblock],
            decomp->nblocks == 0 ? 0 : decomp->consssize[decomp->idxsmallestblock],
            decomp->nblocks == 0 ? 0 : decomp->varssize[decomp->idxsmallestblock]);
   ptr += snprintf(ptr, SCIP_MAXSTRLEN,
            "Border has %d constraints and %d variables\n",
            decomp->labels[0] == SCIP_DECOMP_LINKVAR ? decomp->consssize[0] : 0,
            decomp->labels[0] == SCIP_DECOMP_LINKVAR ? decomp->varssize[0] : 0
            );

   ptr += snprintf(ptr, SCIP_MAXSTRLEN,
            "Modularity: %.3f, Area Score: %.3f\n",
            decomp->modularity, decomp->areascore);
   ptr += snprintf(ptr, SCIP_MAXSTRLEN,
            "Constraint Block Graph: %d edges, %d articulation nodes, %d connected components, %d min., %d max. degree\n",
            decomp->nedges, -1, decomp->ncomponents, decomp->mindegree, decomp->maxdegree);

   return strbuf;
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
      SCIP_CALL( SCIPdecompCreate(&decomp, SCIPblkmem(scip), SCIPdecompGetNBlocks(origdecomp), original, SCIPdecompUseBendersLabels(origdecomp)) );

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
