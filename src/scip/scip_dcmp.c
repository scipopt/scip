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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   decomp.c
 * @brief  methods for working with decompositions
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include "scip/struct_dcmp.h"
#include "scip/debug.h"
#include "scip/dcmp.h"
#include "scip/mem.h"
#include "scip/pub_misc.h"
#include "scip/pub_var.h"
#include "scip/scip_cons.h"
#include "scip/scip_prob.h"
#include "scip/scip_var.h"
#include "scip/scip_mem.h"
#include "scip/struct_scip.h"
#include "scip/pub_cons.h"
#include "scip/pub_dcmp.h"
#include "scip/scip_dcmp.h"
#include "scip/scip_general.h"
#include "scip/scip_var.h"
#include "scip/scip_datastructures.h"
#include "scip/scip_message.h"


#define LABEL_UNASSIGNED INT_MIN /* label constraints or variables as unassigned. Only for internal use */

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

/** raises an error if the condition is not TRUE */
static
SCIP_RETCODE ensureCondition(
   SCIP_Bool             condition           /**< some condition that must hold */
   )
{
   return condition ? SCIP_OKAY : SCIP_ERROR;
}

/** query the constraints active variables and their labels */
static
SCIP_RETCODE decompGetConsVarsAndLabels(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS*            cons,               /**< the constraint */
   SCIP_VAR**            varbuf,             /**< variable buffer array */
   int*                  labelbuf,           /**< buffer to store labels, or NULL if not needed */
   int                   bufsize,            /**< size of buffer arrays */
   int*                  nvars,              /**< pointer to store number of variables */
   int*                  requiredsize,       /**< pointer to store required size */
   SCIP_Bool*            success             /**< pointer to store whether variables and labels were successfully inserted */
   )
{
   SCIP_Bool success2;

   assert(scip != NULL);
   assert(decomp != NULL);
   assert(cons != NULL);
   assert(varbuf != NULL);
   assert(nvars != NULL);
   assert(requiredsize != NULL);
   assert(success != NULL);

   *success = FALSE;
   *requiredsize = 0;
   *nvars = 0;
   SCIP_CALL( SCIPgetConsNVars(scip, cons, nvars, &success2) );

   /* the constraint does not have the corresponding callback */
   if( ! success2 )
   {
      return SCIP_OKAY;
   }

   if( bufsize < *nvars )
   {
      *requiredsize = *nvars;

      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPgetConsVars(scip, cons, varbuf, bufsize, &success2) );
   /* the constraint does not have the corresponding callback */
   if( ! success2 )
   {
      return SCIP_OKAY;
   }

   if( ! SCIPdecompIsOriginal(decomp) )
   {
      SCIP_CALL( SCIPgetActiveVars(scip, varbuf, nvars, bufsize, requiredsize) );

      if( *requiredsize > bufsize )
         return SCIP_OKAY;
   }
   else
   {
      int v;
      for( v = 0; v < *nvars; ++v )
      {
         assert(SCIPvarIsActive(varbuf[v]) || SCIPvarIsNegated(varbuf[v]));

         /* some constraint handlers such as indicator may already return inactive variables */
         if( SCIPvarIsNegated(varbuf[v]) )
            varbuf[v] = SCIPvarGetNegatedVar(varbuf[v]);
      }
   }

   /* get variables labels, if requested */
   if( labelbuf != NULL )
   {
      SCIPdecompGetVarsLabels(decomp, varbuf, labelbuf, *nvars);
   }

   *success = TRUE;

   return SCIP_OKAY;
}

/** add decomposition to SCIP */
SCIP_RETCODE SCIPaddDecomposition(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition to add */
   )
{
   SCIP_Bool isoriginal;

   assert(scip != NULL);
   assert(decomp != NULL);

   isoriginal = SCIPdecompIsOriginal(decomp);

   SCIP_CALL( SCIPcheckStage(scip, "SCIPaddDecomposition", FALSE, isoriginal, isoriginal, isoriginal, isoriginal, TRUE, TRUE, TRUE, !isoriginal, !isoriginal, !isoriginal, FALSE, FALSE, FALSE) );

   SCIP_CALL( SCIPdecompstoreAdd(scip->decompstore, decomp) );

   return SCIP_OKAY;
}

/** get available user decompositions for either the original or transformed problem */
void SCIPgetDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP***        decomps,            /**< pointer to store decompositions array */
   int*                  ndecomps,           /**< pointer to store number of decompositions */
   SCIP_Bool             original            /**< should the decompositions for the original problem be returned? */
   )
{
   assert(scip != NULL);

   SCIP_CALL_ABORT( SCIPcheckStage(scip, "SCIPaddDecomposition", FALSE, original, original, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE) );

   if( decomps != NULL )
      *decomps = original ? SCIPdecompstoreGetOrigDecomps(scip->decompstore) : SCIPdecompstoreGetDecomps(scip->decompstore);

   if( ndecomps != NULL )
      *ndecomps = original ? SCIPdecompstoreGetNOrigDecomps(scip->decompstore) : SCIPdecompstoreGetNDecomps(scip->decompstore);
}

/** returns TRUE if this constraint contains only linking variables */
SCIP_RETCODE SCIPhasConsOnlyLinkVars(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS*            cons,               /**< the constraint */
   SCIP_Bool*            hasonlylinkvars     /**< will be set to TRUE if this constraint has only linking variables */
   )
{
   SCIP_VAR** consvars;
   int nvars;
   int i;
   SCIP_Bool success;

   assert(scip != NULL);
   assert(cons != NULL);
   assert(decomp != NULL);
   assert(hasonlylinkvars != NULL);

   SCIP_CALL( SCIPgetConsNVars(scip, cons, &nvars, &success) );
   SCIP_CALL( ensureCondition(success) );

   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );

   SCIP_CALL( SCIPgetConsVars(scip, cons, consvars, nvars, &success) );
   SCIP_CALL( ensureCondition(success) );

   if( ! SCIPdecompIsOriginal(decomp) )
   {
      int requiredsize;
      SCIP_CALL( SCIPgetActiveVars(scip, consvars, &nvars, nvars, &requiredsize) );
      assert(requiredsize <= nvars);
   }

   *hasonlylinkvars = TRUE;
   /* check if variables are all linking variables */
   for( i = 0; i < nvars && *hasonlylinkvars; ++i )
   {
      int label;

      SCIPdecompGetVarsLabels(decomp, &consvars[i], &label, 1);

      *hasonlylinkvars = (label == SCIP_DECOMP_LINKVAR);
   }

   SCIPfreeBufferArray(scip, &consvars);

   return SCIP_OKAY;
}


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

      SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, conss[c], varbuffer, varlabels,
            twicenvars, &nconsvars, &requiredsize, &success) );
      SCIP_CALL( ensureCondition(success) );

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

      SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, conss[c], varbuffer, NULL,
            twicenvars, &nconsvars, &requiredsize, &success) );
      SCIP_CALL( ensureCondition(success) );

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

/** assign linking constraints to blocks
 *
 * Each linking constraint is assigned to the most frequent block among its variables.
 * Variables of other blocks are relabeled as linking variables.
 *
 * @note: In contrast to SCIPdecompComputeConsLabels(), this method potentially relabels variables.
 */
SCIP_RETCODE SCIPdecompAssignLinkConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp,             /**< decomposition data structure */
   SCIP_CONS**           conss,              /**< array of linking constraints that should be reassigned */
   int                   nconss              /**< number of constraints */
   )
{
   SCIP_VAR** vars;
   int* varslabels;
   int requiredsize;
   int nconsvars;
   int nvars;
   int c;

   assert(scip != NULL);
   assert(decomp != NULL);

   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );

   for( c = 0; c < nconss; c++ )
   {
      SCIP_Bool success;

      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, &success) );
      SCIP_CALL( SCIPgetConsVars(scip, conss[c], vars, nvars, &success) );

      if( ! SCIPdecompIsOriginal(decomp) )
      {
         SCIP_CALL( SCIPgetActiveVars(scip, vars, &nconsvars, nconsvars, &requiredsize) );
         assert(requiredsize <= nvars);
      }

      SCIPdecompGetVarsLabels(decomp, vars, varslabels, nconsvars);

      SCIPsortIntPtr(varslabels, (void **)vars, nconsvars);
      /* constraint contains only linking variables */
      if( varslabels[nconsvars - 1] == SCIP_DECOMP_LINKVAR )
      {
         /* todo */
         continue;
      }
      else
      {
         int startposs[2];
         int endposs[2];
         int nlinkvars;
         int block;
         int maxnblockvars;
         int curr;
         int v;
         int p;

         /* count linking variables */
         if( varslabels[0] == SCIP_DECOMP_LINKVAR )
         {
            nlinkvars = countLabelFromPos(varslabels, 0, nconsvars);
         }
         else
         {
            nlinkvars = 0;
         }

         assert(nlinkvars < nconsvars);

         curr = nlinkvars;
         /* find the most frequent block label among the nonlinking variables */
         maxnblockvars = 0;
         block = -1;
         do
         {
            int nblockvars = countLabelFromPos(varslabels, curr, nconsvars);
            if (nblockvars > maxnblockvars)
            {
               maxnblockvars = nblockvars;
               block = curr;
            }
            curr += nblockvars;
         }
         while (curr < nconsvars);

         /* reassign all variables from other blocks as linking variables */
         startposs[0] = nlinkvars;
         endposs[0] = block;
         startposs[1] = block + maxnblockvars;
         endposs[1] = nconsvars;

         p = 0;
         /* loop over all variables before (p==0) and after (p==1) the most frequent block label */
         for( p = 0; p < 2; ++p )
         {
            /* relabel */
            for( v = startposs[p]; v < endposs[p]; ++v)
               varslabels[v] = SCIP_DECOMP_LINKVAR;

            /* set labels in the decomposition */
            SCIP_CALL( SCIPdecompSetVarsLabels(decomp, &vars[startposs[p]], &varslabels[startposs[p]], endposs[p] - startposs[p]) );
         }

         SCIP_CALL( SCIPdecompSetConsLabels(decomp, &conss[c], &varslabels[block], 1) );
      }
   }

   SCIPfreeBufferArray(scip, &vars);
   SCIPfreeBufferArray(scip, &varslabels);

   return SCIP_OKAY;
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
      int requiredsize;
      SCIP_Bool success;
      SCIP_Bool found;

      /* linking constraints do not contribute to the modularity */
      if( conslabel == SCIP_DECOMP_LINKCONS )
         continue;

      /* find the position of the constraint label. Constraints of the border always belong to the first block at index 0 */
      found = SCIPsortedvecFindInt(decomp->labels, conslabel, decomp->nblocks + 1, &blockpos);
      assert(found);

      SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, conss[c], varbuf, varslabels,
               nvars, &nconsvars, &requiredsize, &success) );
      SCIP_CALL( ensureCondition(success) );

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
   nnonzeroes = MAX(nnonzeroes, 1);
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
   SCIP_Real areascore = 1.0;
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
         areascore -= (factor * decomp->consssize[i]) * decomp->varssize[i];
      }

      areascore -= ((SCIP_Real)nlinkconss * nvars + (SCIP_Real)nconss * nlinkvars - (SCIP_Real)nlinkconss * nlinkvars) * factor;
   }

   decomp->areascore = areascore;
}

/** build the block decomposition graph */
static
SCIP_RETCODE buildBlockGraph(
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

   /* capture the trivial case that no linking variables are present */
   if( decomp->varssize[0] == 0 || decomp->nblocks == 0 )
   {
      decomp->mindegree = 0;
      decomp->maxdegree = 0;
      decomp->nedges = 0;
      decomp->ncomponents = SCIPdecompGetNBlocks(decomp);
      decomp->narticulations = 0;

      return SCIP_OKAY;
   }

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   nblocks = SCIPdecompGetNBlocks(decomp);

   SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
   SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &linkvaridx, nvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nvars) );

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
      else
         linkvaridx[v] = -1;
   }

   /* create a bipartite graph composed of block and linking var nodes */
   SCIP_CALL( SCIPcreateDigraph(scip, &blocklinkingvargraph, nblocks + nlinkingvars) );/* nblocks does not include the linking constraints block */

   for( i = 0; i < nconss; ++i )
   {
      int requiredsize;
      /* linking constraints are skipped in this checking step */
      if( conslabels[i] == SCIP_DECOMP_LINKCONS )
         continue;

      SCIP_CALL( decompGetConsVarsAndLabels(scip, decomp, conss[i], consvars, varlabels,
            nvars, &nconsvars, &requiredsize, &success) );
      SCIP_CALL( ensureCondition(success) );

      /* adding double-direction arcs between blocks and corresponding linking variables */
      for( j = 0; j < nconsvars; ++j )
      {
         assert(consvars[j] != NULL);
         if( varlabels[j] == SCIP_DECOMP_LINKVAR )
         {
            int linkingvarnodeidx = linkvaridx[SCIPvarGetProbindex(consvars[j])];
            int blocknodeidx;

            assert(linkingvarnodeidx >= 0);
            /* find the position of the constraint label. Subtract later by 1 to get the node index as the 1st block is reserved for linking constraints */
            found = SCIPsortedvecFindInt(decomp->labels, conslabels[i], decomp->nblocks + 1, &blocknodeidx); /* assuming labels are sorted */
            assert(found);

            SCIP_CALL( SCIPdigraphAddArcSafe(blocklinkingvargraph, nblocks + linkingvarnodeidx, blocknodeidx - 1, NULL) );
            SCIP_CALL( SCIPdigraphAddArcSafe(blocklinkingvargraph, blocknodeidx - 1, nblocks + linkingvarnodeidx, NULL) );
         }
      }
   }
   assert(SCIPdigraphGetNNodes(blocklinkingvargraph) > 0);

   /* From the information of the above bipartite graph, build the block-decomposition graph: nodes -> blocks and double-direction arcs -> linking variables */
   SCIP_CALL( SCIPcreateDigraph(scip, &blockgraph, nblocks) );

   for( n = nblocks; n < SCIPdigraphGetNNodes(blocklinkingvargraph); ++n )
   {
      nsucc = (int) SCIPdigraphGetNSuccessors(blocklinkingvargraph, n);
      succnodes = (int*) SCIPdigraphGetSuccessors(blocklinkingvargraph, n);
      for( succ1 = 0; succ1 < nsucc; ++succ1 )
      {
         for( succ2 = 0; succ2 < nsucc; ++succ2 )
         {
            if( succnodes[succ1] != succnodes[succ2] ) /* no self-loops */
               SCIP_CALL( SCIPdigraphAddArcSafe(blockgraph, succnodes[succ1], succnodes[succ2], NULL) );
         }
      }
   }

   assert(SCIPdigraphGetNNodes(blockgraph) > 0);

   /* Get the number of edges in the block-decomposition graph.*/
   decomp->nedges = SCIPdigraphGetNArcs(blockgraph) / 2;

   /* Get the minimum and maximum degree of the block-decomposition graph */
   tempmin = (int) SCIPdigraphGetNSuccessors(blockgraph, 0);
   tempmax = (int) SCIPdigraphGetNSuccessors(blockgraph, 0);
   for( n = 1; n < SCIPdigraphGetNNodes(blockgraph); ++n )
   {
      nsucc = (int) SCIPdigraphGetNSuccessors(blockgraph, n);
      if( nsucc < tempmin )
         tempmin = nsucc;
      else if( nsucc > tempmax )
         tempmax = nsucc;
   }

   decomp->mindegree = tempmin;
   decomp->maxdegree = tempmax;

   /* Get the number of connected components in the block-decomposition graph.*/
   SCIP_CALL( SCIPdigraphComputeUndirectedComponents(blockgraph, -1, NULL, NULL) );
   decomp->ncomponents = SCIPdigraphGetNComponents(blockgraph);

   /* Get the number of articulation nodes in the block-decomposition graph using DFS.*/
   SCIP_CALL( SCIPdigraphGetArticulationPoints(blockgraph, NULL) );
   decomp->narticulations = SCIPdigraphGetNArticulationPoints(blockgraph);

   SCIPfreeBufferArray(scip, &consvars);
   SCIPfreeBufferArray(scip, &linkvaridx);
   SCIPfreeBufferArray(scip, &varlabels);
   SCIPfreeBufferArray(scip, &conslabels);
   SCIPdigraphFree(&blockgraph);
   SCIPdigraphFree(&blocklinkingvargraph);

   return SCIP_OKAY;
}

/** compute decomposition statistics and store them in the decomp object */
SCIP_RETCODE SCIPcomputeDecompStats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_DECOMP*          decomp              /**< decomposition data structure */
   )
{
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_CONS** conssarray;
   SCIP_VAR** varsarray;
   int* varslabels;
   int* conslabels;
   int nvars;
   int nconss;
   int varblockstart;
   int consblockstart;
   int currlabelidx;
   int varidx;
   int considx;
   int i;

   assert(scip != NULL);
   assert(decomp != NULL);

   /* store variable and constraint labels in buffer arrays */
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

  /* return if problem is empty
   *
   * TODO ensure that statistics reflect this correctly
   */
  if( nvars == 0 || nconss == 0 )
  {
     return SCIP_OKAY;
  }
  conss = SCIPgetConss(scip);
  vars = SCIPgetVars(scip);

  SCIP_CALL( SCIPduplicateBufferArray(scip, &conssarray, conss, nconss) );
  SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
  SCIP_CALL( SCIPduplicateBufferArray(scip, &varsarray, vars, nvars) );
  SCIP_CALL( SCIPallocBufferArray(scip, &varslabels, nvars) );

  SCIPdecompGetVarsLabels(decomp, varsarray, varslabels, nvars);
  SCIPdecompGetConsLabels(decomp, conssarray, conslabels, nconss);

  /* sort both buffer arrays for quick counting */
  SCIPsortIntPtr(varslabels, (void**)varsarray, nvars);
  SCIPsortIntPtr(conslabels, (void**)conssarray, nconss);

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

      assert(currlabelidx < decomp->memsize);
      /* store the smaller of the two current labels */
      decomp->labels[currlabelidx] = MIN(varlabel, conslabel);

      /* a strictly larger variable label means that there are no variables for the current label */
      if( varlabel <= conslabel )
         decomp->varssize[currlabelidx] = countLabelFromPos(varslabels, varidx, nvars);
      else
         decomp->varssize[currlabelidx] = 0;

      /* the same for constraint labels */
      if( conslabel <= varlabel )
         decomp->consssize[currlabelidx] = countLabelFromPos(conslabels, considx, nconss);
      else
         decomp->consssize[currlabelidx] = 0;

      /* increase indices appropriately */
      varidx += decomp->varssize[currlabelidx];
      considx += decomp->consssize[currlabelidx];

      currlabelidx++;
  }

  SCIPdebugMsg(scip, "Counted %d different labels (should be %d)\n", currlabelidx, decomp->nblocks + 1);

  /* strip the remaining, unused blocks */
  if( currlabelidx < decomp->nblocks + 1 )
     decomp->nblocks = currlabelidx - 1;

  /* delete empty blocks from statistics, relabel the corresponding constraints/variables as linking */
  varblockstart = decomp->varssize[0];
  consblockstart = decomp->consssize[0];

  for( i = 1; i < decomp->nblocks + 1; ++i )
  {
     assert(MAX(decomp->varssize[i], decomp->consssize[i]) > 0);
     /* relabel constraint blocks as linking, if there are no corresponding variables */
     if( decomp->varssize[i] == 0 )
     {
        int nblockconss = decomp->consssize[i];
        int c;
        /* relabel these constraints as linking */
        for( c = consblockstart; c < consblockstart + nblockconss; ++c )
           conslabels[c] = SCIP_DECOMP_LINKCONS;

        SCIP_CALL( SCIPdecompSetConsLabels(decomp, &conssarray[consblockstart], &conslabels[consblockstart], nblockconss) );

        /* increase number of linking constraints */
        decomp->consssize[0] += nblockconss;
     }

     /* same for constraints */
     if( decomp->consssize[i] == 0 )
     {
        int nblockvars = decomp->varssize[i];
        int v;

        /* relabel the variables as linking variables */
        for( v = varblockstart; v < varblockstart + nblockvars; ++v )
           varslabels[v] = SCIP_DECOMP_LINKVAR;

        SCIP_CALL( SCIPdecompSetVarsLabels(decomp, &varsarray[varblockstart], &varslabels[varblockstart], nblockvars) );

        /* increase number of linking variables */
        decomp->varssize[0] += nblockvars;
     }

     varblockstart += decomp->varssize[i];
     consblockstart += decomp->consssize[i];
  }

  currlabelidx = 1;

  /* delete empty blocks; they are no longer present */
  for( i = 1; i < decomp->nblocks + 1; ++i )
  {
     /* keep only nonempty blocks */
     if( decomp->varssize[i] > 0 && decomp->consssize[i] > 0 )
     {
        decomp->labels[currlabelidx] = decomp->labels[i];
        decomp->varssize[currlabelidx] = decomp->varssize[i];
        decomp->consssize[currlabelidx] = decomp->consssize[i];

        currlabelidx++;
     }
  }

  decomp->nblocks = currlabelidx - 1;

  decomp->idxsmallestblock = decomp->idxlargestblock = -1;
  /* now that indices are fixed, store indices with largest and smallest number of constraints */
  for( i = 1; i < decomp->nblocks + 1; ++i )
  {
     if( decomp->idxsmallestblock == -1 )
        decomp->idxsmallestblock = decomp->idxlargestblock = i;
     else if( decomp->consssize[decomp->idxsmallestblock] > decomp->consssize[i] )
        decomp->idxsmallestblock = i;
     else if( decomp->consssize[decomp->idxlargestblock] < decomp->consssize[i] )
        decomp->idxlargestblock = i;
  }

  /* compute more involved statistics such as the area score, the modularity, and the block graph statistics */
  SCIP_CALL( computeModularity(scip, decomp, &decomp->modularity) );

  computeAreaScore(scip, decomp);

  SCIP_CALL( buildBlockGraph(scip, decomp) );

  SCIPfreeBufferArray(scip, &varslabels);
  SCIPfreeBufferArray(scip, &varsarray);
  SCIPfreeBufferArray(scip, &conslabels);
  SCIPfreeBufferArray(scip, &conssarray);

   return SCIP_OKAY;
}
