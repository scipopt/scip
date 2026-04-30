/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   iisfinder_greedy.c
 * @brief  greedy deletion and addition filter heuristic to compute (I)ISs
 * @author Marc Pfetsch
 * @author Mark Turner
 * @author Paul Meinhold
 */

#include "scip/iisfinder_greedy.h"
#include "scip/struct_iisfinder.h"
#include "scip/scip_datastructures.h"
#include "scip/pub_misc_sort.h"

#define IISFINDER_NAME           "greedy"
#define IISFINDER_DESC           "greedy deletion or addition constraint deletion"
#define IISFINDER_PRIORITY        8000
#define IISFINDER_ENABLE          TRUE

#define DEFAULT_TIMELIMPERITER   1e+20 /**< time limit of optimization process for each individual subproblem */
#define DEFAULT_NODELIMPERITER   -1L   /**< node limit of optimization process for each individual subproblem */

#define DEFAULT_ADDITIVE         TRUE  /**< should an additive constraint approach be used instead of deletion */
#define DEFAULT_CONSERVATIVE     TRUE  /**< should an unsolved problem (by e.g. user interrupt, node limit, time limit) be considered feasible when deleting constraints */
#define DEFAULT_DELAFTERADD      TRUE  /**< should the deletion routine be performed after the addition routine (in the case of additive) */
#define DEFAULT_DYNAMICREORDERING TRUE /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */
#define DEFAULT_DETECTCOMPONENTS TRUE  /**< should the deletion filter detect and delete feasible disconnected components */
#define DEFAULT_COMPONENTMINSIZE 8     /**< number of constraints a component must have at least to be detected */

#define DEFAULT_INITBATCHSIZE    16    /**< the initial batchsize for the first iteration, ignored if initrelbatchsize is positive */
#define DEFAULT_INITRELBATCHSIZE 0.03125 /**< the initial batchsize relative to the original problem for the first iteration (0.0: use initbatchsize) */
#define DEFAULT_MAXBATCHSIZE     INT_MAX /**< the maximum batchsize per iteration */
#define DEFAULT_MAXRELBATCHSIZE  0.5   /**< the maximum batchsize relative to the original problem per iteration */
#define DEFAULT_BATCHINGFACTOR   2.0   /**< the factor with which the batchsize is multiplied in every update */
#define DEFAULT_BATCHINGOFFSET   0.0   /**< the offset which is added to the multiplied batchsize in every update */
#define DEFAULT_BATCHUPDATEINTERVAL 1  /**< the number of iterations to run with a constant batchsize before updating (1: always update) */


/*
 * Data structures
 */

/** IIS finder data */
struct SCIP_IISfinderData
{
   SCIP_Real             timelimperiter;     /**< time limit of optimization process for each individual subproblem */
   SCIP_Longint          nodelimperiter;     /**< node limit of optimization process for each individual subproblem */

   SCIP_Bool             additive;           /**< should an additive constraint approach be used instead of deletion */
   SCIP_Bool             conservative;       /**< should an unsolved problem (by e.g. user interrupt, node limit, time limit) be considered feasible when deleting constraints */
   SCIP_Bool             delafteradd;        /**< should the deletion routine be performed after the addition routine (in the case of additive) */
   SCIP_Bool             dynamicreordering;  /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */
   SCIP_Bool             detectcomponents;   /**< should the deletion filter detect and delete feasible disconnected components */
   int                   componentminsize;   /**< number of constraints a component must have at least to be detected */

   int                   initbatchsize;      /**< the initial batchsize for the first iteration, ignored if initrelbatchsize is positive */
   SCIP_Real             initrelbatchsize;   /**< the initial batchsize relative to the original problem for the first iteration (0.0: use initbatchsize) */
   int                   maxbatchsize;       /**< the maximum batchsize per iteration */
   SCIP_Real             maxrelbatchsize;    /**< the maximum batchsize relative to the original problem per iteration */
   SCIP_Real             batchingfactor;     /**< the factor with which the batchsize is multiplied in every update */
   SCIP_Real             batchingoffset;     /**< the offset which is added to the multiplied batchsize in every update */
   int                   batchupdateinterval; /**< the number of iterations to run with a constant batchsize before updating (1: always update) */
};

/*
 * Local methods
 */

/* Set time and node limits on the subproblem */
static
SCIP_RETCODE setLimits(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip  */
   SCIP_Real             timelim,            /**< total time limit allowed of the whole call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelim,            /**< node limit allowed for the whole call */
   SCIP_Longint          nodelimperiter      /**< maximum number of nodes per individual solve call */
   )
{
   SCIP_Real currtime;
   SCIP_Real mintimelim;
   SCIP_Longint globalnodelim;

   /* Set the time limit for the solve call. Take into account the global time limit, the current time used, and the time lim on each individual call */
   currtime = SCIPiisGetTime(iis);
   mintimelim = MIN(timelim - currtime, timelimperiter);
   mintimelim = MAX(mintimelim, 0);
   SCIP_CALL( SCIPsetRealParam(scip, "limits/time", mintimelim) );

   /* Set the node limit for the solve call. Take into account the global node limit, the current nodes processed, and the node lim on each individual call */
   if( nodelim == -1 )
      globalnodelim = -1;
   else
   {
      globalnodelim = nodelim - SCIPiisGetNNodes(iis);
      assert( globalnodelim >= 0 );
   }
   if( globalnodelim == -1 && nodelimperiter == -1 )
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", nodelim) );
   }
   else if( globalnodelim == -1 || nodelimperiter == -1 )
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", MAX(globalnodelim, nodelimperiter)) );
   }
   else
   {
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", MIN(globalnodelim, nodelimperiter)) );
   }
   return SCIP_OKAY;
}

/* Revert bound changes made in the subproblem */
static
SCIP_RETCODE revertBndChgs(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_VAR**            vars,               /**< the array of variables whose bounds are changed */
   SCIP_Real*            bounds,             /**< the array of original bounds for the variables */
   int*                  idxs,               /**< the indices of the vars (in the vars array) that have been deleted */
   int                   ndelbounds,         /**< the number of bounds that will be deleted */
   SCIP_Bool             islb                /**< are the bounds that are being deleted LBs? */
   )
{
   int i;

   /* Iterate over the affected variables and restore the bounds back to their original values */
   for (i = 0; i < ndelbounds; ++i)
   {
      if( islb )
      {
         if( !SCIPisInfinity(scip, -bounds[i]) )
            SCIP_CALL( SCIPchgVarLb(scip, vars[idxs[i]], bounds[i]) );
      }
      else
      {
         if( !SCIPisInfinity(scip, bounds[i]) )
            SCIP_CALL( SCIPchgVarUb(scip, vars[idxs[i]], bounds[i]) );
      }
   }
   return SCIP_OKAY;
}

/* Revert deleted constraint changes made in the subproblem */
static
SCIP_RETCODE revertConssDeletions(
   SCIP*                 scip,               /**< SCIP instance to analyze */
   SCIP_CONS**           conss,              /**< the array of constraints where some have been deleted */
   int*                  idxs,               /**< the indices of the cons (in the conss array) that have been deleted */
   int                   ndelconss,          /**< the number of constraints that have been deleted */
   SCIP_Bool             releaseonly         /**< Should the constraints just be released instead of added back */
   )
{
   int i;
   SCIP_CONS* copycons;

   for( i = 0; i < ndelconss; ++i )
   {
      if( releaseonly )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &conss[idxs[i]]) );
      }
      else
      {
         SCIP_CALL( SCIPaddCons(scip, conss[idxs[i]]) );
         copycons = conss[idxs[i]];
         assert(SCIPconsGetNUses(copycons) > 1);
         SCIP_CALL( SCIPreleaseCons(scip, &copycons) );
      }
   }

   return SCIP_OKAY;
}

/* Set initial and maximum batchsize for given parameters */
static
void setInitAndMaxBatchsize(
   int                   probsize,           /**< the size of the problem (e.g., nvars or nconss) to compute meaningful ratios */
   SCIP_Real             initrelbatchsize,   /**< the initial batchsize relative to the original problem for the first iteration (0.0: use initbatchsize) */
   SCIP_Real             maxrelbatchsize,    /**< the maximum batchsize relative to the original problem for the first iteration */
   int*                  initbatchsize,      /**< the initial batchsize */
   int*                  maxbatchsize        /**< the maximum batchsize per iteration */
   )
{
   int maxrel;

   assert(initbatchsize != NULL);
   assert(*initbatchsize >= 1);
   assert(maxbatchsize != NULL);
   assert(*maxbatchsize >= 1);

   /* compute the maximum batchsize */
   maxrel = (int)ceil(maxrelbatchsize * probsize);
   if( *maxbatchsize > maxrel )
      *maxbatchsize = maxrel;
   if( *maxbatchsize < 1 )
      *maxbatchsize = 1;

   /* compute the initial batchsize */
   if( initrelbatchsize > 0.0 )
      *initbatchsize = (int)ceil(initrelbatchsize * probsize);
   if( *initbatchsize > *maxbatchsize )
      *initbatchsize = *maxbatchsize;
   if( *initbatchsize < 1 )
      *initbatchsize = 1;
}

/* Update the batchsize according to the update formula and keep it in its limits */
static
SCIP_RETCODE updateBatchsize(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   initbatchsize,      /**< the initial batchsize */
   int                   maxbatchsize,       /**< the maximum batchsize per iteration */
   int                   iteration,          /**< the current iteration */
   SCIP_Bool             resettoinit,        /**< should the batchsize be reset to the initial batchsize? */
   SCIP_Real             batchingfactor,     /**< the factor with which the batchsize is multiplied in every update */
   SCIP_Real             batchingoffset,     /**< the offset which is added to the multiplied batchsize in every update */
   int                   batchupdateinterval, /**< the number of iterations to run with a constant batchsize before updating (1: always update) */
   int*                  batchsize           /**< the batchsize to be updated */
   )
{
   if( resettoinit )
      *batchsize = initbatchsize;
   else if( iteration % batchupdateinterval == 0 )
      *batchsize = (int)ceil(batchingfactor * (*batchsize) + batchingoffset);

   /* respect limits and maximum */
   *batchsize = MIN(*batchsize, maxbatchsize);
   *batchsize = MAX(*batchsize, 1);
   SCIPdebugMsg(scip, "Updated batchsize to %d\n", *batchsize);

   return SCIP_OKAY;
}

/* Detect disconnected components and sort them by size */
static
SCIP_RETCODE detectComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< the array of constraints */
   int**                 varconsidxs,        /**< the 2d-array containing the problem indices of the constraints of a variable */
   int*                  nvarconsidxs,       /**< the array containing the numbers of constraint indices */
   int                   nconss,             /**< the number of constraints */
   int                   nvars,              /**< the number of all variables */
   int                   ndeleted,           /**< the number of deleted constraints */
   int                   componentminsize,   /**< number of constraints a component must have at least to be detected */
   SCIP_DIGRAPH**        digraph,            /**< digraph data structure */
   int**                 components,         /**< array of indices of sorted valid components, will be NULL if only one component exists */
   int*                  ncomponents         /**< the number of valid components, will be 0 if only one component exists */
   )
{
   int* compsizes;
   int* nodes;
   int* sizes;
   int* firstconss;
   int nallcomponents;
   int nnodes;
   int v;
   int c;
   int i;
   int j;

   assert(ndeleted <= nconss);

   /* create digraph with nconss nodes, which will be built as a sparse connection digraph,
    * i.e., only the first constraint of a variable will be connected with all its other constraints
    * (note: indicator and referenced linear constraints will always be connected by their slack variable) */
   SCIP_CALL( SCIPcreateDigraph(scip, digraph, nconss) );

   /* allocate the array to hold number of successors of each node in the graph */
   SCIP_CALL( SCIPallocBufferArray(scip, &sizes, nconss) );
   for( c = 0; c < nconss; ++c )
      sizes[c] = 0;

   /* for each variable, find the index of its first constraint and set node successor sizes */
   SCIP_CALL( SCIPallocBufferArray(scip, &firstconss, nvars) );
   for( v = 0; v < nvars; ++v )
   {
      /* default value of -1 if there is no first constraint */
      firstconss[v] = -1;

      for( i = 0; i < nvarconsidxs[v]; ++i )
      {
         c = varconsidxs[v][i];

         /* skip deleted constraints */
         if( conss[c] == NULL )
            continue;

         if( firstconss[v] == -1 )
            firstconss[v] = c;
         else
         {
            /* the first constraint will have the current variable constraint as a successor */
            ++sizes[firstconss[v]];

            /* due to future undirectedness, the current variable constraint will have the first constraint as a successor */
            ++sizes[c];
         }
      }
   }
   SCIP_CALL( SCIPdigraphSetSizes(*digraph, sizes) );
   SCIPfreeBufferArray(scip, &sizes);

   /* fill sparse connection digraph */
   for( v = 0; v < nvars; ++v )
   {
      /* skip variables without a first constraint */
      if( firstconss[v] == -1 )
         continue;

      /* connect the first constraint with all subsequent constraints of this variable */
      for( i = 1; i < nvarconsidxs[v]; ++i )
      {
         c = varconsidxs[v][i];

         /* skip deleted constraints and the first constraint itself */
         if( conss[c] == NULL || c == firstconss[v] )
            continue;

         SCIP_CALL( SCIPdigraphAddArc(*digraph, firstconss[v], c, NULL) );
      }
   }

   /* compute components */
   SCIP_CALL( SCIPdigraphComputeUndirectedComponents(*digraph, componentminsize, NULL, NULL) );
   *ncomponents = nallcomponents = SCIPdigraphGetNComponents(*digraph);

   /* special case if componentminsize == 1, since deleted constraints will be returned as singleton components, which we need to exclude */
   if( componentminsize == 1 )
   {
      *ncomponents -= ndeleted;
      assert(*ncomponents >= 1);
   }

   if( *ncomponents <= 1 )
   {
      *components = NULL;
      *ncomponents = 0;
   }
   else
   {
      /* build the valid components array */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, components, *ncomponents) );
      SCIP_CALL( SCIPallocBufferArray(scip, &compsizes, *ncomponents) );
      j = 0;
      for( i = 0; i < nallcomponents; ++i )
      {
         SCIPdigraphGetComponent(*digraph, i, &nodes, &nnodes);

         /* skip deleted constraints that make up a component */
         if( nnodes == 1 && conss[nodes[0]] == NULL )
            continue;

         (*components)[j] = i;
         compsizes[j] = nnodes;
         ++j;
      }
      assert(j == *ncomponents);

      /* sort by size */
      SCIPsortIntInt(compsizes, *components, *ncomponents);

      SCIPfreeBufferArray(scip, &compsizes);
   }

   SCIPfreeBufferArray(scip, &firstconss);

   return SCIP_OKAY;
}

/* for a given component, get the complement constraints of that component */
static
SCIP_RETCODE getComplementConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< the array of all constraints */
   int                   nconss,             /**< the number of all constraints */
   int                   ndeleted,           /**< the number of deleted constraints */
   int*                  nodes,              /**< array of component constraint indices */
   int                   nnodes,             /**< number of component constraint indices */
   int*                  cplidxs,            /**< array of complement constraint indices */
   int                   ncplidxs            /**< number of complement constraint indices */
   )
{
   int nodeidx;
   int c;
   int k;

   assert(nconss == ndeleted + nnodes + ncplidxs);

   /* sort nodes */
   SCIPsortInt(nodes, nnodes);

   /* fill index array */
   k = 0;
   nodeidx = 0;
   for( c = 0; c < nconss; ++c )
   {
      /* skip deleted constraints */
      if( conss[c] == NULL )
         continue;

      /* skip component constraints */
      if( nodeidx < nnodes && c == nodes[nodeidx] )
      {
         ++nodeidx;
         continue;
      }

      /* add that constraint index to the complement array */
      cplidxs[k] = c;
      ++k;
   }
   assert(nodeidx == nnodes);
   assert(k == ncplidxs);

   return SCIP_OKAY;
}

/* build the varconsidxs array that maps a variable index to an array of associated constraint indices */
static
SCIP_RETCODE buildVarConsIdxs(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< the array of constraints */
   int                   nconss,             /**< the number of constraints */
   int                   nvars,              /**< the number of original variables */
   int***                varconsidxs,        /**< the 2d-array containing the indices of the constraints of a variable */
   int**                 nvarconsidxs,       /**< the array of the number of constraints of a variable */
   SCIP_Bool*            success             /**< whether the varconsidxs array was built correctly */
   )
{
   SCIP_VAR** consvars;
   SCIP_Bool* varoccurred;
   int nconsvars;
   int c;
   int v;
   int i;

   assert(scip != NULL);
   assert(conss != NULL || nconss == 0);
   assert(varconsidxs != NULL);
   assert(nvarconsidxs != NULL);
   assert(success != NULL);

   *success = TRUE;

   /* allocate and initialize count array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, nvarconsidxs, nvars) );
   for( v = 0; v < nvars; ++v )
      (*nvarconsidxs)[v] = 0;

   /* allocate array of pointers */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, varconsidxs, nvars) );
   for( v = 0; v < nvars; ++v )
      (*varconsidxs)[v] = NULL;

   /* allocate temporary arrays to track unique occurrences of variables and constraint variables */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &varoccurred, nvars) );

   /* first pass: count in how many constraints each variable appears in */
   for( c = 0; c < nconss; ++c )
   {
      /* get the (possibly duplicate) variables in this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, success) );
      if( !(*success) )
         break;
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
      SCIP_CALL( SCIPgetConsVars(scip, conss[c], consvars, nconsvars, success) );
      if( !(*success) )
      {
         SCIPfreeBufferArray(scip, &consvars);
         break;
      }

      /* for each unique constraint variable, increment the constraint counter (variable corresponds to this constraint) */
      for( i = 0; i < nconsvars; ++i )
      {
         /* make sure to get non-negated variables */
         v = SCIPvarGetProbindex((SCIPvarIsNegated(consvars[i])) ? SCIPvarGetNegatedVar(consvars[i]) : consvars[i]);
         assert(v >= 0);
         assert(v < nvars);

         /* skip multiple occurrences of variables for the current constraint */
         if( varoccurred[v] )
            continue;

         varoccurred[v] = TRUE;
         ++(*nvarconsidxs)[v];
      }

      /* clear occurrence array for the next constraint */
      for( i = 0; i < nconsvars; ++i )
      {
         /* make sure to get non-negated variables */
         v = SCIPvarGetProbindex((SCIPvarIsNegated(consvars[i])) ? SCIPvarGetNegatedVar(consvars[i]) : consvars[i]);

         varoccurred[v] = FALSE;
      }

      SCIPfreeBufferArray(scip, &consvars);
   }

   /* clean up arrays after an unsuccessful SCIPgetCons[N]Vars() */
   if( !(*success) )
   {
      for( v = 0; v < nvars; ++v )
         varoccurred[v] = FALSE;

      SCIPfreeCleanBufferArray(scip, &varoccurred);
      SCIPfreeBlockMemoryArray(scip, varconsidxs, nvars);
      SCIPfreeBlockMemoryArray(scip, nvarconsidxs, nvars);

      return SCIP_OKAY;
   }

   /* allocate arrays for each variable based on counts */
   for( v = 0; v < nvars; ++v )
   {
      if( (*nvarconsidxs)[v] > 0 )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*varconsidxs)[v], (*nvarconsidxs)[v]) );

         /* reset for second pass (filling phase) to use as index array */
         (*nvarconsidxs)[v] = 0;
      }
   }

   /* second pass: fill the arrays with constraint indices */
   for( c = 0; c < nconss; ++c )
   {
      /* get the (possibly duplicate) variables in this constraint */
      SCIP_CALL( SCIPgetConsNVars(scip, conss[c], &nconsvars, success) );
      assert(*success);
      SCIP_CALL( SCIPallocBufferArray(scip, &consvars, nconsvars) );
      SCIP_CALL( SCIPgetConsVars(scip, conss[c], consvars, nconsvars, success) );
      assert(*success);

      /* for each unique constraint variable, write the constraint index to varconsidxs */
      for( i = 0; i < nconsvars; ++i )
      {
         /* make sure to get non-negated variables */
         v = SCIPvarGetProbindex((SCIPvarIsNegated(consvars[i])) ? SCIPvarGetNegatedVar(consvars[i]) : consvars[i]);

         /* skip multiple occurrences of variables for the current constraint */
         if( varoccurred[v] )
            continue;

         varoccurred[v] = TRUE;
         (*varconsidxs)[v][(*nvarconsidxs)[v]] = c;
         ++(*nvarconsidxs)[v];
      }

      /* clear occurrence array for the next constraint */
      for( i = 0; i < nconsvars; ++i )
      {
         /* make sure to get non-negated variables */
         v = SCIPvarGetProbindex((SCIPvarIsNegated(consvars[i])) ? SCIPvarGetNegatedVar(consvars[i]) : consvars[i]);

         varoccurred[v] = FALSE;
      }

      SCIPfreeBufferArray(scip, &consvars);
   }

   SCIPfreeCleanBufferArray(scip, &varoccurred);

   return SCIP_OKAY;
}

/* free the varconsidxs array */
static
void freeVarConsIdxs(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nvars,              /**< the number of all variables */
   int***                varconsidxs,        /**< the 2d-array containing the indices of the constraints of a variable */
   int**                 nvarconsidxs        /**< the array of the number of constraints of a variable */
   )
{
   int i;

   assert(scip != NULL);
   assert(varconsidxs != NULL);
   assert(nvarconsidxs != NULL);

   /* free individual arrays */
   for( i = 0; i < nvars; ++i )
   {
      if( (*varconsidxs)[i] != NULL )
      {
         SCIPfreeBlockMemoryArray(scip, &(*varconsidxs)[i], (*nvarconsidxs)[i]);
      }
   }

   SCIPfreeBlockMemoryArray(scip, varconsidxs, nvars);
   SCIPfreeBlockMemoryArray(scip, nvarconsidxs, nvars);
}

/** solve subproblem for deletionFilter */
static
SCIP_RETCODE deletionSubproblem(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip */
   SCIP_CONS**           conss,              /**< The array of constraints (may be a superset of the current constraints) */
   SCIP_VAR**            vars,               /**< the array of vars */
   int*                  idxs,               /**< the indices of the constraints / vars that will be deleted / bounds removed */
   int                   ndels,              /**< the number of bounds that will be deleted */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reaches some status that could result in infeasible, e.g. node limit */
   SCIP_Bool             delbounds,          /**< whether bounds should be deleted instead of constraints */
   SCIP_Bool             islb,               /**< are the bounds that are being deleted LBs? */
   SCIP_Bool             forcetest,          /**< do an infeasibility test even though no changes are made to the problem? */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */
   SCIP_Bool*            deleted,            /**< have the deleted bounds or constraints stayed deleted */
   SCIP_Bool*            feasible,           /**< whether the subproblem is proven feasible */
   SCIP_Bool*            stop,               /**< pointer to store whether we have to stop */
   SCIP_Bool*            alldeletionssolved  /**< pointer to store whether all the subscips solved or NULL if not used */
   )
{
   SCIP* scip;
   SCIP_Real* bounds = NULL;
   SCIP_RETCODE retcode;
   SCIP_STATUS status;
   SCIP_Bool chgmade = FALSE;
   int i;

   assert(deleted != NULL);
   assert(feasible != NULL);
   assert(stop != NULL);
   assert(alldeletionssolved != NULL);

   *deleted = FALSE;
   *feasible = FALSE;
   *stop = FALSE;
   scip = SCIPiisGetSubscip(iis);

   /* remove bounds or constraints */
   if( delbounds )
   {
      assert(vars != NULL);
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &bounds, ndels) );
      for (i = 0; i < ndels; ++i)
      {
         if( islb )
         {
            bounds[i] = SCIPvarGetLbOriginal(vars[idxs[i]]);
            if( !SCIPisInfinity(scip, -bounds[i]) )
            {
               SCIP_CALL( SCIPchgVarLb(scip, vars[idxs[i]], -SCIPinfinity(scip)) );
               chgmade = TRUE;
            }
         }
         else
         {
            bounds[i] = SCIPvarGetUbOriginal(vars[idxs[i]]);
            if( !SCIPisInfinity(scip, bounds[i]) )
            {
               SCIP_CALL( SCIPchgVarUb(scip, vars[idxs[i]], SCIPinfinity(scip)) );
               chgmade = TRUE;
            }
         }
      }
   }
   else
   {
      assert(conss != NULL);

      if( ndels > 0 )
         chgmade = TRUE;
      for (i = 0; i < ndels; ++i)
      {
         assert( SCIPconsIsInProb(conss[idxs[i]]) );
         SCIP_CALL( SCIPcaptureCons(scip, conss[idxs[i]]) );
         SCIP_CALL( SCIPdelCons(scip, conss[idxs[i]]) );
      }
   }

   if( !chgmade && !forcetest )
   {
      if( delbounds )
         SCIPfreeBlockMemoryArray(scip, &bounds, ndels);
      return SCIP_OKAY;
   }

   /* solve problem until first solution is found or infeasibility has been proven */
   SCIP_CALL( setLimits(scip, iis, timelim, timelimperiter, nodelim, nodelimperiter) );
   retcode = SCIPsolve(scip);
   SCIPiisAddNNodes(iis, SCIPgetNTotalNodes(scip));

   if( retcode != SCIP_OKAY )
   {
      SCIP_CALL( SCIPfreeTransform(scip) );
      SCIPdebugMsg(scip, "Error in sub-scip with deleted constraints / bounds. Re-adding them.\n");
      if( delbounds )
      {
         SCIP_CALL( revertBndChgs(scip, vars, bounds, idxs, ndels, islb) );
         SCIPfreeBlockMemoryArray(scip, &bounds, ndels);
      }
      else
      {
         SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, FALSE) );
      }
      *alldeletionssolved = FALSE;
      return SCIP_OKAY;
   }

   status = SCIPgetStatus(scip);
   SCIP_CALL( SCIPfreeTransform(scip) );

   /* check status and handle accordingly */
   switch ( status )
   {
      case SCIP_STATUS_USERINTERRUPT:    /* if a user interrupt occurred, just stop */
      case SCIP_STATUS_TERMINATE:
         SCIPdebugMsg(scip, "User interrupt. Stopping.\n");
         if( delbounds )
         {
            SCIP_CALL( revertBndChgs(scip, vars, bounds, idxs, ndels, islb) );
         }
         else
         {
            SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, !conservative) );
         }
         *stop = TRUE;
         *alldeletionssolved = FALSE;
         break;

      case SCIP_STATUS_TIMELIMIT:        /* if we reached some status that may have ended up in an infeasible problem */
      case SCIP_STATUS_NODELIMIT:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_DUALLIMIT:
         *alldeletionssolved = FALSE;
         SCIPdebugMsg(scip, "Some limit reached. Keeping bounds / constraints removed if non-conservative.\n");
         if( !conservative )
         {
            SCIPiisSetSubscipInfeasible(iis, FALSE);
            *deleted = TRUE;
         }
         if( conservative && delbounds )
         {
            SCIP_CALL( revertBndChgs(scip, vars, bounds, idxs, ndels, islb) );
         }
         if( !delbounds )
         {
            SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, !conservative) );
         }
         break;

      case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
         SCIPdebugMsg(scip, "Subproblem with bounds / constraints removed infeasible. Keep them removed.\n");
         SCIPiisSetSubscipInfeasible(iis, TRUE);
         if( !delbounds )
         {
            SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, TRUE) );
         }
         *deleted = TRUE;
         break;

      case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_PRIMALLIMIT:
         SCIPdebugMsg(scip, "Found solution to subproblem with bounds / constraints removed. Add them back.\n");
         if( delbounds )
         {
            SCIP_CALL( revertBndChgs(scip, vars, bounds, idxs, ndels, islb) );
         }
         else
         {
            SCIP_CALL( revertConssDeletions(scip, conss, idxs, ndels, FALSE) );
         }
         *feasible = TRUE;
         break;

      case SCIP_STATUS_UNKNOWN:
      default:
         *alldeletionssolved = FALSE;
         SCIPerrorMessage("Unexpected return status %d in removed bounds subproblem. Exiting...\n", status);
         if( delbounds )
            SCIPfreeBlockMemoryArray(scip, &bounds, ndels);
         return SCIP_ERROR;
   }

   if( !silent && *deleted )
      SCIPiisfinderInfoMessage(iis, FALSE);

   if( delbounds )
      SCIPfreeBlockMemoryArray(scip, &bounds, ndels);

   assert(!(*deleted && *feasible));

   return SCIP_OKAY;
}

/* delete all but the smallest infeasible component */
static
SCIP_RETCODE deleteComponents(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip */
   SCIP_DIGRAPH*         digraph,            /**< digraph data structure */
   SCIP_CONS**           conss,              /**< the array of constraints */
   int                   nconss,             /**< the number of constraints */
   int*                  components,         /**< array of indices of sorted valid components, will be NULL if only one component exists */
   int                   ncomponents,        /**< the number of valid components, will be 0 if only one component exists */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool             conservative,       /**< whether we treat a subproblem to be feasible, if it reaches some status that could result in infeasible, e.g. node limit */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */
   SCIP_Bool*            stopiter,           /**< pointer to store whether we have to stop */
   int*                  ndeleted            /**< the number of deleted constraints */
   )
{
   SCIP* scip;
   int* cplidxs;
   int* nodes;
   SCIP_Bool deleted;
   SCIP_Bool feasible;
   SCIP_Bool allfeasibilitychecked;
   SCIP_Bool forcetestlastcomponent;
   int nnodestotal;
   int ncplidxs;
   int nnodes;
   int i;
   int j;

   assert(iis != NULL);
   assert(digraph != NULL);
   assert(conss != NULL);
   assert(stopiter != NULL);

   scip = SCIPiisGetSubscip(iis);

   if( ncomponents >= 2 )
   {
      assert(components != NULL);

      nnodestotal = *ndeleted;
      allfeasibilitychecked = TRUE;
      forcetestlastcomponent = FALSE;
      for( j = 0; j < ncomponents; ++j )
      {
         /* get the current component */
         SCIPdigraphGetComponent(digraph, components[j], &nodes, &nnodes);

         /* get the complement of the current component */
         ncplidxs = nconss - nnodes - *ndeleted;
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &cplidxs, ncplidxs) );
         SCIP_CALL( getComplementConstraints(scip, conss, nconss, *ndeleted, nodes, nnodes, cplidxs, ncplidxs) );

         /* if all components but the last one are proven feasible, test the last component for feasibility
          * (unless any component was ignored in component detection due to being smaller that componentminsize) */
         nnodestotal += nnodes;
         if( nnodestotal == nconss && allfeasibilitychecked )
            forcetestlastcomponent = TRUE;

         /* try to delete the complement of current component */
         SCIP_CALL( deletionSubproblem(iis, conss, NULL, cplidxs, ncplidxs, timelim, timelimperiter, nodelim, nodelimperiter,
               conservative, FALSE, FALSE, forcetestlastcomponent, silent, &deleted, &feasible, stopiter, &allfeasibilitychecked) );

         SCIPfreeBlockMemoryArray(scip, &cplidxs, ncplidxs);

         /* if the complement has been deleted, the component is proven infeasible */
         if( deleted )
         {
            *ndeleted += ncplidxs;
            break;
         }

         if( feasible )
         {
            /* if the last component is feasible as well, exit with error */
            if( forcetestlastcomponent )
            {
               SCIPinfoMessage(scip, NULL, "Error during component detection. All components are feasible. Abort.\n");
               *stopiter = TRUE;
               break;
            }

            /* if the current component is feasible, delete it */
            for( i = 0; i < nnodes; ++i )
            {
               SCIP_CALL( SCIPdelCons(scip, conss[nodes[i]]) );
               conss[nodes[i]] = NULL;
            }
            *ndeleted += nnodes;
         }

         if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || *stopiter )
            break;
      }

      SCIPfreeBlockMemoryArray(scip, &components, ncomponents);
   }

   return SCIP_OKAY;
}

/** solve subproblem for additionFilter */
static
SCIP_RETCODE additionSubproblem(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip  */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool*            feasible,           /**< pointer to store whether the problem is feasible */
   SCIP_Bool*            stop                /**< pointer to store whether we have to stop */
   )
{
   SCIP* scip;
   SCIP_RETCODE retcode;
   SCIP_STATUS status;

   assert( stop != NULL );
   scip = SCIPiisGetSubscip(iis);

   *stop = FALSE;

   /* solve problem until first solution is found or infeasibility has been proven */
   SCIP_CALL( setLimits(scip, iis, timelim, timelimperiter, nodelim, nodelimperiter) );
   retcode = SCIPsolve(scip);

   if( retcode != SCIP_OKAY )
   {
      SCIPdebugMsg(scip, "Error in sub-scip with added constraints. Keep added constraints.\n");
      return SCIP_ERROR;
   }

   SCIPiisAddNNodes(iis, SCIPgetNTotalNodes(scip));
   status = SCIPgetStatus(scip);

   /* check status */
   switch ( status )
   {
      case SCIP_STATUS_USERINTERRUPT:    /* if an user interrupt occurred, just stop */
      case SCIP_STATUS_TERMINATE:
         SCIPdebugMsg(scip, "User interrupt. Stopping.\n");
         *stop = TRUE;
         break;

      case SCIP_STATUS_NODELIMIT:        /* if we reached some limit */
      case SCIP_STATUS_TIMELIMIT:
      case SCIP_STATUS_TOTALNODELIMIT:
      case SCIP_STATUS_STALLNODELIMIT:
      case SCIP_STATUS_MEMLIMIT:
      case SCIP_STATUS_GAPLIMIT:
      case SCIP_STATUS_SOLLIMIT:
      case SCIP_STATUS_RESTARTLIMIT:
      case SCIP_STATUS_INFORUNBD:
      case SCIP_STATUS_DUALLIMIT:
         SCIPdebugMsg(scip, "Some limit reached. Added constraint batch failed to induce infeasibility. Continue adding.\n");
         break;

      case SCIP_STATUS_INFEASIBLE:       /* if the problem is infeasible */
         SCIPdebugMsg(scip, "Subproblem with added constraints infeasible. Final batch of constraints added.\n");
         *feasible = FALSE;
         break;

      case SCIP_STATUS_BESTSOLLIMIT:     /* we found a solution */
      case SCIP_STATUS_OPTIMAL:
      case SCIP_STATUS_UNBOUNDED:
      case SCIP_STATUS_PRIMALLIMIT:
         SCIPdebugMsg(scip, "Found solution of subproblem with added constraints. Keep adding constraint batches.\n");
         *feasible = TRUE;
         break;

      case SCIP_STATUS_UNKNOWN:
      default:
         SCIPerrorMessage("Unexpected return status %d. Exiting ...\n", status);
         return SCIP_ERROR;
   }

   return SCIP_OKAY;
}

/** Deletion filter to greedily remove constraints to obtain an (I)IS */
static
SCIP_RETCODE deletionFilterBatch(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip  */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Bool             removebounds,       /**< Whether the algorithm should remove bounds as well as constraints */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */

   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool             conservative,       /**< should a node or time limit solve be counted as feasible when deleting constraints */

   int                   initbatchsize,      /**< the initial batchsize for the first iteration */
   SCIP_Real             initrelbatchsize,   /**< the initial batchsize relative to the original problem for the first iteration (0.0: use initbatchsize) */
   int                   maxbatchsize,       /**< the maximum batchsize per iteration */
   SCIP_Real             maxrelbatchsize,    /**< the maximum batchsize relative to the original problem for the first iteration */
   SCIP_Real             batchingfactor,     /**< the factor with which the batchsize is multiplied in every update */
   SCIP_Real             batchingoffset,     /**< the offset which is added to the multiplied batchsize in every update */
   int                   batchupdateinterval, /**< the number of iterations to run with a constant batchsize before updating (1: always update) */

   SCIP_Bool             detectcomponents,   /**< should the deletion filter detect and delete disconnected components */
   int                   componentminsize,   /**< number of constraints a component must have at least to be detected */
   SCIP_Bool*            alldeletionssolved  /**< pointer to store whether all the subscips solved */
   )
{
   SCIP* scip;
   SCIP_CONS** origconss;
   SCIP_CONS** conss;
   SCIP_VAR** origvars;
   SCIP_VAR** vars;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_DIGRAPH* digraph;
   int** varconsidxs;
   int* order;
   int* idxs;
   int* nvarconsidxs;
   int* components;
   SCIP_Bool success;
   SCIP_Bool stopiter;
   SCIP_Bool deleted;
   SCIP_Bool feasible;
   int nconss;
   int nvars;
   int ndeleted;
   int batchindex;
   int batchsize;
   int appliedinitbatchsize;
   int appliedmaxbatchsize;
   int iteration;
   int ncomponents;
   int i;
   int j;

   /* get current subscip */
   scip = SCIPiisGetSubscip(iis);
   assert( scip != NULL );
   assert( SCIPiisIsSubscipInfeasible(iis) );

   /* get random generator */
   randnumgen = SCIPiisGetRandnumgen(iis);
   assert( randnumgen != NULL );

   /* get variable and constraint information */
   nvars = SCIPgetNOrigVars(scip);
   nconss = SCIPgetNOrigConss(scip);
   origconss = SCIPgetOrigConss(scip);
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, origconss, nconss) );

   /* get initial and maximum batchsize */
   appliedinitbatchsize = initbatchsize;
   appliedmaxbatchsize = maxbatchsize;
   setInitAndMaxBatchsize(nconss, initrelbatchsize, maxrelbatchsize, &appliedinitbatchsize, &appliedmaxbatchsize);

   /* allocate indices array */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idxs, appliedmaxbatchsize) );

   /* reset problem */
   SCIP_CALL( SCIPfreeTransform(scip) );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

   /* prepare random order for constraints */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);

   /* create the varconsidxs array for quicker component graph creation */
   if( detectcomponents )
   {
      SCIP_CALL( buildVarConsIdxs(scip, conss, nconss, nvars, &varconsidxs, &nvarconsidxs, &success) );
      if( !success )
      {
         SCIPinfoMessage(scip, NULL, "Failed to successfully setup component detection. Skip it.");
         detectcomponents = FALSE;
      }
   }
   else
   {
      varconsidxs = NULL;
      nvarconsidxs = NULL;
   }

   /* Loop through all batches of constraints in random order */
   i = 0;
   ndeleted = 0;
   iteration = 0;
   ncomponents = 0;
   deleted = TRUE;
   stopiter = FALSE;
   batchsize = appliedinitbatchsize;
   while( i < nconss )
   {
      /* do component detection and deletion */
      if( detectcomponents && deleted )
      {
         assert(nvarconsidxs != NULL);
         assert(varconsidxs != NULL);

         SCIP_CALL( detectComponents(scip, conss, varconsidxs, nvarconsidxs, nconss, nvars, ndeleted, componentminsize, &digraph, &components, &ncomponents) );
         SCIP_CALL( deleteComponents(iis, digraph, conss, nconss, components, ncomponents, timelim, timelimperiter, nodelim, nodelimperiter, conservative, silent, &stopiter, &ndeleted) );

         SCIPdebugMsg(scip, "Detected %d components.\n", ncomponents);

         SCIPdigraphFree(&digraph);

         if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
            break;
      }

      j = 0;
      batchindex = i;
      while( i < nconss && j < batchsize )
      {
         /* only add non-deleted independent constraints to the batch, e.g., exclude indicated constraints */
         if( conss[order[i]] != NULL && SCIPconsGetNUses(conss[order[i]]) == 1 )
         {
            idxs[j] = order[i];
            ++j;
         }
         ++i;
      }
      if( j == 0 )
         break;

      /* treat subproblem */
      SCIP_CALL( deletionSubproblem(iis, conss, NULL, idxs, j, timelim, timelimperiter, nodelim, nodelimperiter,
            conservative, FALSE, FALSE, FALSE, silent, &deleted, &feasible, &stopiter, alldeletionssolved) );
      if( deleted )
         ndeleted += j;
      if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
         break;

      /* reset i to beginning of current batch if batch has not been deleted and j was large */
      if( !deleted && j > appliedinitbatchsize )
         i = batchindex;

      ++iteration;

      /* update batchsize */
      SCIP_CALL( updateBatchsize(scip, appliedinitbatchsize, appliedmaxbatchsize, iteration, !deleted, batchingfactor, batchingoffset, batchupdateinterval, &batchsize) );

      assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
   }

   if( detectcomponents )
      freeVarConsIdxs(scip, nvars, &varconsidxs, &nvarconsidxs);

   SCIPfreeBlockMemoryArray(scip, &order, nconss);
   SCIPfreeBlockMemoryArray(scip, &idxs, appliedmaxbatchsize);
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);

   if( *alldeletionssolved && appliedinitbatchsize == 1 )
      SCIPiisSetSubscipIrreducible(iis, TRUE);

   if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
      return SCIP_OKAY;

   /* Repeat the above procedure but for bounds instead of constraints */
   if( removebounds )
   {
      /* get variables */
      origvars = SCIPgetOrigVars(scip);
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &vars, origvars, nvars) );

      /* get initial and maximum batchsize */
      appliedinitbatchsize = initbatchsize;
      appliedmaxbatchsize = maxbatchsize;
      setInitAndMaxBatchsize(nvars, initrelbatchsize, maxrelbatchsize, &appliedinitbatchsize, &appliedmaxbatchsize);
      batchsize = appliedinitbatchsize;

      /* allocate indices array */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &idxs, appliedmaxbatchsize) );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nvars) );
      for (i = 0; i < nvars; ++i)
         order[i] = i;
      SCIPrandomPermuteIntArray(randnumgen, order, 0, nvars);

      i = 0;
      iteration = 0;
      deleted = FALSE;
      while( i < nvars )
      {
         j = 0;
         batchindex = i;
         /* Do not delete bounds of binary variables or bother with calculations of free variables */
         while( i < nvars && j < batchsize )
         {
            if( (SCIPvarGetType(vars[order[i]]) != SCIP_VARTYPE_BINARY) && (!SCIPisInfinity(scip, -SCIPvarGetLbOriginal(vars[order[i]])) || !SCIPisInfinity(scip, SCIPvarGetUbOriginal(vars[order[i]]))) )
            {
               idxs[j] = order[i];
               ++j;
            }
            ++i;
         }
         if( j == 0 )
            break;

         /* treat subproblem with LB deletions */
         SCIP_CALL( deletionSubproblem(iis, NULL, vars, idxs, j, timelim, timelimperiter, nodelim, nodelimperiter,
               conservative, TRUE, TRUE, FALSE, silent, &deleted, &feasible, &stopiter, alldeletionssolved) );
         if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
            break;

         /* treat subproblem with UB deletions */
         SCIP_CALL( deletionSubproblem(iis, NULL, vars, idxs, j, timelim, timelimperiter, nodelim, nodelimperiter,
               conservative, TRUE, FALSE, FALSE, silent, &deleted, &feasible, &stopiter, alldeletionssolved) );
         if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) || stopiter )
            break;

         /* reset i to beginning of current batch if batch has not been deleted and j was large */
         if( !deleted && j > appliedinitbatchsize )
            i = batchindex;

         ++iteration;

         /* update batchsize */
         SCIP_CALL( updateBatchsize(scip, appliedinitbatchsize, appliedmaxbatchsize, iteration, !deleted, batchingfactor, batchingoffset, batchupdateinterval, &batchsize) );

         assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
      }

      SCIPfreeBlockMemoryArray(scip, &order, nvars);
      SCIPfreeBlockMemoryArray(scip, &idxs, appliedmaxbatchsize);
      SCIPfreeBlockMemoryArray(scip, &vars, nvars);
   }

   return SCIP_OKAY;
}

/** Addition filter to greedily add constraints to obtain an (I)IS */
static
SCIP_RETCODE additionFilterBatch(
   SCIP_IIS*             iis,                /**< IIS data structure containing subscip  */
   SCIP_Real             timelim,            /**< The global time limit on the IIS finder call */
   SCIP_Longint          nodelim,            /**< The global node limit on the IIS finder call */
   SCIP_Bool             silent,             /**< should the run be performed silently without printing progress information */

   SCIP_Real             timelimperiter,     /**< time limit per individual solve call */
   SCIP_Longint          nodelimperiter,     /**< maximum number of nodes per individual solve call */
   SCIP_Bool             dynamicreordering,  /**< should satisfied constraints outside the batch of an intermediate solve be added during the additive method */

   int                   initbatchsize,      /**< the initial batchsize for the first iteration */
   SCIP_Real             initrelbatchsize,   /**< the initial batchsize relative to the original problem for the first iteration (0.0: use initbatchsize) */
   int                   maxbatchsize,       /**< the maximum batchsize per iteration */
   SCIP_Real             maxrelbatchsize,    /**< the maximum batchsize relative to the original problem for the first iteration */
   SCIP_Real             batchingfactor,     /**< the factor with which the batchsize is multiplied in every update */
   SCIP_Real             batchingoffset,     /**< the offset which is added to the multiplied batchsize in every update */
   int                   batchupdateinterval /**< the number of iterations to run with a constant batchsize before updating (1: always update) */
   )
{
   SCIP* scip;
   SCIP_CONS** origconss;
   SCIP_CONS** conss;
   SCIP_RANDNUMGEN* randnumgen;
   SCIP_SOL* sol;
   SCIP_SOL* copysol;
   SCIP_Bool* inIS;
   int* order;
   SCIP_Bool feasible;
   SCIP_Bool stopiter;
   SCIP_RETCODE retcode;
   SCIP_RESULT result;
   int nconss;
   int batchsize;
   int iteration;
   int i;
   int j;
   int k;

   /* get current subscip */
   scip = SCIPiisGetSubscip(iis);
   assert( scip != NULL );
   assert( SCIPiisIsSubscipInfeasible(iis) );

   /* get random generator */
   randnumgen = SCIPiisGetRandnumgen(iis);
   assert( randnumgen != NULL );

   /* get constraint information */
   nconss = SCIPgetNOrigConss(scip);
   origconss = SCIPgetOrigConss(scip);
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &conss, origconss, nconss) );

   /* get initial and maximum batchsize */
   setInitAndMaxBatchsize(nconss, initrelbatchsize, maxrelbatchsize, &initbatchsize, &maxbatchsize);
   batchsize = initbatchsize;

   /* Initialise information for whether a constraint is in the final infeasible system */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &inIS, nconss) );

   /* First capture and then delete all constraints */
   assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
   for( i = 0; i < nconss; ++i )
   {
      assert( SCIPconsIsInProb(conss[i]) );
      if( SCIPconsGetNUses(conss[i]) > 1 )
      {
         inIS[i] = TRUE;
      }
      else
      {
         SCIP_CALL( SCIPcaptureCons(scip, conss[i]) );
         SCIP_CALL( SCIPdelCons(scip, conss[i]) );
         inIS[i] = FALSE;
      }
   }
   SCIPiisSetSubscipInfeasible(iis, FALSE);

   /* Prepare random order in which the constraints will be added back */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &order, nconss) );
   for (i = 0; i < nconss; ++i)
      order[i] = i;
   SCIPrandomPermuteIntArray(randnumgen, order, 0, nconss);

   /* Continue to add constraints until an infeasible status is reached */
   i = 0;
   iteration = 0;
   feasible = TRUE;
   stopiter = FALSE;
   while( i < nconss )
   {
      /* Add the next batch of constraints */
      k = 0;
      while( i < nconss && k < batchsize )
      {
         if( !inIS[order[i]] )
         {
            SCIP_CALL( SCIPaddCons(scip, conss[order[i]]) );
            SCIP_CALL( SCIPreleaseCons(scip, &conss[order[i]]) );
            inIS[order[i]] = TRUE;
            ++k;
         }
         i++;
      }

      /* We have the full infeasible problem again */
      if( i == nconss )
      {
         feasible = FALSE;
         break;
      }

      /* Solve the reduced problem */
      retcode = additionSubproblem(iis, timelim, timelimperiter, nodelim, nodelimperiter, &feasible, &stopiter);
      if( !feasible || stopiter || SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) )
      {
         SCIP_CALL( SCIPfreeTransform(scip) );
         break;
      }

      if( !silent )
         SCIPiisfinderInfoMessage(iis, FALSE);

      if( dynamicreordering && retcode == SCIP_OKAY )
      {
         /* free transform and copy solution if there is one */
         copysol = NULL;
         sol = SCIPgetBestSol(scip);
         if( sol != NULL )
         {
            SCIP_CALL( SCIPcreateSolCopyOrig(scip, &copysol, sol) );
            SCIP_CALL( SCIPunlinkSol(scip, copysol) );
         }
         SCIP_CALL( SCIPfreeTransform(scip) );
         assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );

         /* Add any other constraints that are also feasible for the current solution */
         if( copysol != NULL )
         {
            k = 0;
            for( j = i; j < nconss; ++j )
            {
               /* Don't dynamically add indicator constraints */
               if( !inIS[order[j]]
                  && strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[order[j]])), "indicator") != 0 )
               {
                  SCIP_CALL( SCIPcheckCons(scip, conss[order[j]], copysol, FALSE, FALSE, FALSE, &result) );
                  if( result == SCIP_FEASIBLE )
                  {
                     SCIP_CALL( SCIPaddCons(scip, conss[order[j]]) );
                     SCIP_CALL( SCIPreleaseCons(scip, &conss[order[j]]) );
                     inIS[order[j]] = TRUE;
                     k++;
                  }
               }
            }
            if( k > 0 )
            {
               SCIPdebugMsg(scip, "Added %d constraints by reordering dynamically.\n", k);
            }
            SCIP_CALL( SCIPfreeSol(scip, &copysol) );
         }
      }
      else
      {
         SCIP_CALL( SCIPfreeTransform(scip) );
         assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
      }

      ++iteration;

      /* update batchsize */
      SCIP_CALL( updateBatchsize(scip, initbatchsize, maxbatchsize, iteration, FALSE, batchingfactor, batchingoffset, batchupdateinterval, &batchsize) );

      assert( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM );
   }

   SCIPiisSetSubscipInfeasible(iis, !feasible);
   if( !silent )
      SCIPiisfinderInfoMessage(iis, FALSE);

   /* Release any cons not in the IS */
   for( i = 0; i < nconss; ++i )
   {
      if( !inIS[order[i]] )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &conss[order[i]]) );
      }
   }

   SCIPfreeBlockMemoryArray(scip, &order, nconss);
   SCIPfreeBlockMemoryArray(scip, &inIS, nconss);
   SCIPfreeBlockMemoryArray(scip, &conss, nconss);

   return SCIP_OKAY;
}

/** perform a greedy addition or deletion algorithm to obtain an infeasible subsystem (IS).
 *
 *  This is the generation method for the greedy IIS finder rule.
 *  Depending on the parameter choices, constraints are either greedily added from an empty problem,
 *  or deleted from a complete problem. In the case of constraints being added, this is done until the problem
 *  becomes infeasible, after which one can then begin deleting constraints. In the case of deleting constraints,
 *  this is done until no more constraints (or batches of constraints) can be deleted without making
 *  the problem feasible.
 */
static
SCIP_RETCODE execIISfinderGreedy(
   SCIP_IIS*             iis,                /**< IIS data structure */
   SCIP_IISFINDERDATA*   iisfinderdata,      /**< IIS finder data */
   SCIP_RESULT*          result              /**< pointer to store the result of the IIS finder run. SCIP_DIDNOTFIND if the algorithm failed, otherwise SCIP_SUCCESS. */
   )
{
   SCIP* scip = SCIPiisGetSubscip(iis);
   SCIP_Real timelim;
   SCIP_Longint nodelim;
   SCIP_Bool removebounds;
   SCIP_Bool silent;
   SCIP_Bool alldeletionssolved = TRUE;

   assert( scip != NULL );
   assert( iisfinderdata != NULL );
   assert( result != NULL );

   SCIP_CALL( SCIPgetRealParam(scip, "iis/time", &timelim) );
   SCIP_CALL( SCIPgetLongintParam(scip, "iis/nodes", &nodelim) );
   SCIP_CALL( SCIPgetBoolParam(scip, "iis/removebounds", &removebounds) );
   SCIP_CALL( SCIPgetBoolParam(scip, "iis/silent", &silent) );

   *result = SCIP_SUCCESS;

   if( iisfinderdata->additive )
   {
      if( !silent )
      {
         SCIPdebugMsg(scip, "----- STARTING GREEDY ADDITION ALGORITHM -----\n");
      }
      SCIP_CALL( additionFilterBatch(iis, timelim, nodelim, silent, iisfinderdata->timelimperiter,
            iisfinderdata->nodelimperiter, iisfinderdata->dynamicreordering, iisfinderdata->initbatchsize,
            iisfinderdata->initrelbatchsize, iisfinderdata->maxbatchsize, iisfinderdata->maxrelbatchsize,
            iisfinderdata->batchingfactor, iisfinderdata->batchingoffset, iisfinderdata->batchupdateinterval) );
      SCIPiisSetSubscipIrreducible(iis, FALSE);
      if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) )
         return SCIP_OKAY;
   }
   else
   {
      if( !silent )
      {
         SCIPdebugMsg(scip, "----- STARTING GREEDY DELETION ALGORITHM -----\n");
      }
      SCIP_CALL( deletionFilterBatch(iis, timelim, nodelim, removebounds, silent, iisfinderdata->timelimperiter,
            iisfinderdata->nodelimperiter, iisfinderdata->conservative, iisfinderdata->initbatchsize,
            iisfinderdata->initrelbatchsize, iisfinderdata->maxbatchsize, iisfinderdata->maxrelbatchsize,
            iisfinderdata->batchingfactor, iisfinderdata->batchingoffset, iisfinderdata->batchupdateinterval,
            iisfinderdata->detectcomponents, iisfinderdata->componentminsize, &alldeletionssolved) );
      if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) )
         return SCIP_OKAY;
   }

   if( iisfinderdata->delafteradd && iisfinderdata->additive )
   {
      if( !silent )
      {
         SCIPdebugMsg(scip, "----- STARTING GREEDY DELETION ALGORITHM FOLLOWING COMPLETED ADDITION ALGORITHM -----\n");
      }
      SCIP_CALL( deletionFilterBatch(iis, timelim, nodelim, removebounds, silent, iisfinderdata->timelimperiter,
            iisfinderdata->nodelimperiter, iisfinderdata->conservative, iisfinderdata->initbatchsize,
            iisfinderdata->initrelbatchsize, iisfinderdata->maxbatchsize, iisfinderdata->maxrelbatchsize,
            iisfinderdata->batchingfactor, iisfinderdata->batchingoffset, iisfinderdata->batchupdateinterval,
            iisfinderdata->detectcomponents, iisfinderdata->componentminsize, &alldeletionssolved) );
      if( SCIPiisGetTime(iis) >= timelim || ( nodelim != -1 && SCIPiisGetNNodes(iis) >= nodelim ) )
         return SCIP_OKAY;
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of IIS finder
 */


/** copy method for IIS finder plugin (called when SCIP copies plugins) */
static
SCIP_DECL_IISFINDERCOPY(iisfinderCopyGreedy)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(iisfinder != NULL);

   SCIP_STRINGEQ( SCIPiisfinderGetName(iisfinder), IISFINDER_NAME, SCIP_INVALIDCALL );

   /* call inclusion method of IIS finder */
   SCIP_CALL( SCIPincludeIISfinderGreedy(scip) );

   return SCIP_OKAY;
}

/** destructor of IIS finder to free user data (called when SCIP is exiting) */
/**! [SnippetIISfinderFreeGreedy] */
static
SCIP_DECL_IISFINDERFREE(iisfinderFreeGreedy)
{  /*lint --e{715}*/
   SCIP_IISFINDERDATA* iisfinderdata;

   iisfinderdata = SCIPiisfinderGetData(iisfinder);

   SCIPfreeBlockMemory(scip, &iisfinderdata);

   SCIPiisfinderSetData(iisfinder, NULL);

   return SCIP_OKAY;
}
/**! [SnippetIISfinderFreeGreedy] */

/** IIS finder generation method of IIS */
static
SCIP_DECL_IISFINDEREXEC(iisfinderExecGreedy)
{  /*lint --e{715}*/
   SCIP_IISFINDERDATA* iisfinderdata;

   assert(iisfinder != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTFIND;

   iisfinderdata = SCIPiisfinderGetData(iisfinder);
   assert(iisfinderdata != NULL);

   SCIP_CALL( execIISfinderGreedy(iis, iisfinderdata, result) );

   return SCIP_OKAY;
}


/*
 * IIS finder specific interface methods
 */

/** creates the greedy IIS finder and includes it in SCIP */
SCIP_RETCODE SCIPincludeIISfinderGreedy(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_IISFINDERDATA* iisfinderdata;
   SCIP_IISFINDER* iisfinder;

   /* create greedy IIS finder data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &iisfinderdata) );
   BMSclearMemory(iisfinderdata);

   SCIP_CALL( SCIPincludeIISfinderBasic(scip, &iisfinder, IISFINDER_NAME, IISFINDER_DESC, IISFINDER_PRIORITY,
         IISFINDER_ENABLE, iisfinderExecGreedy, iisfinderdata) );

   assert(iisfinder != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetIISfinderCopy(scip, iisfinder, iisfinderCopyGreedy) );
   SCIP_CALL( SCIPsetIISfinderFree(scip, iisfinder, iisfinderFreeGreedy) );

   /* add greedy IIS finder parameters */
   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/timelimperiter",
         "time limit of optimization process for each individual subproblem",
         &iisfinderdata->timelimperiter, FALSE, DEFAULT_TIMELIMPERITER, 0.0, SCIP_INVALID/10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddLongintParam(scip,
         "iis/" IISFINDER_NAME "/nodelimperiter",
         "node limit of optimization process for each individual subproblem",
         &iisfinderdata->nodelimperiter, FALSE, DEFAULT_NODELIMPERITER, -1L, SCIP_LONGINT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/additive",
         "should an additive constraint approach be used instead of deletion",
         &iisfinderdata->additive, FALSE, DEFAULT_ADDITIVE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/conservative",
         "should an unsolved problem (by e.g. user interrupt, node limit, time limit) be considered feasible when deleting constraints",
         &iisfinderdata->conservative, TRUE, DEFAULT_CONSERVATIVE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/delafteradd",
         "should the deletion routine be performed after the addition routine (in the case of additive)",
         &iisfinderdata->delafteradd, TRUE, DEFAULT_DELAFTERADD, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/dynamicreordering",
         "should satisfied constraints outside the batch of an intermediate solve be added during the additive method",
         &iisfinderdata->dynamicreordering, TRUE, DEFAULT_DYNAMICREORDERING, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "iis/" IISFINDER_NAME "/initbatchsize",
         "the initial batchsize for the first iteration, ignored if initrelbatchsize is positive",
         &iisfinderdata->initbatchsize, FALSE, DEFAULT_INITBATCHSIZE, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/initrelbatchsize",
         "the initial batchsize relative to the original problem for the first iteration (0.0: use initbatchsize)",
         &iisfinderdata->initrelbatchsize, FALSE, DEFAULT_INITRELBATCHSIZE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "iis/" IISFINDER_NAME "/maxbatchsize",
         "the maximum batchsize per iteration",
         &iisfinderdata->maxbatchsize, TRUE, DEFAULT_MAXBATCHSIZE, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/maxrelbatchsize",
         "the maximum batchsize relative to the original problem per iteration",
         &iisfinderdata->maxrelbatchsize, TRUE, DEFAULT_MAXRELBATCHSIZE, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/batchingfactor",
         "the factor with which the batchsize is multiplied in every update",
         &iisfinderdata->batchingfactor, TRUE, DEFAULT_BATCHINGFACTOR, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "iis/" IISFINDER_NAME "/batchingoffset",
         "the offset which is added to the multiplied batchsize in every update",
         &iisfinderdata->batchingoffset, TRUE, DEFAULT_BATCHINGOFFSET, SCIP_REAL_MIN, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "iis/" IISFINDER_NAME "/batchupdateinterval",
         "the number of iterations to run with a constant batchsize before updating (1: always update)",
         &iisfinderdata->batchupdateinterval, TRUE, DEFAULT_BATCHUPDATEINTERVAL, 1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "iis/" IISFINDER_NAME "/detectcomponents",
         "should the deletion filter detect and delete disconnected components",
         &iisfinderdata->detectcomponents, FALSE, DEFAULT_DETECTCOMPONENTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "iis/" IISFINDER_NAME "/componentminsize",
         "number of constraints a component must have at least to be detected",
         &iisfinderdata->componentminsize, FALSE, DEFAULT_COMPONENTMINSIZE, 1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** perform the greedy deletion algorithm with singleton batches to obtain an irreducible infeasible subsystem (IIS) */
SCIP_RETCODE SCIPiisGreedyMakeIrreducible(
   SCIP_IIS*             iis                 /**< IIS data structure */
   )
{
   SCIP* scip = SCIPiisGetSubscip(iis);
   SCIP_CONS** conss;
   SCIP_CONS* imagecons;
   SCIP_HASHMAP* invconssmap = NULL;
   SCIP_Real timelim;
   SCIP_Longint nodelim;
   SCIP_Bool isstandalone;
   SCIP_Bool removebounds;
   SCIP_Bool silent;
   SCIP_Bool alldeletionssolved = TRUE;
   int nconss;
   int c;

   assert( scip != NULL );

   if( !SCIPiisIsSubscipInfeasible(iis) )
   {
      SCIPerrorMessage("infeasible problem required\n");
      return SCIP_INVALIDDATA;
   }

   nconss = SCIPgetNOrigConss(scip);

   /* if this function is called by a user outside of iisfinder.c::SCIPiisGenerate(), build inverse constraints hashmap */
   isstandalone = !SCIPhashmapIsEmpty(iis->conssmap);
   if( isstandalone )
   {
      conss = SCIPgetOrigConss(scip);
      SCIP_CALL( SCIPhashmapCreate(&invconssmap, SCIPblkmem(scip), nconss) );
      for( c = 0; c < nconss; ++c )
      {
         imagecons = SCIPhashmapGetImage(iis->conssmap, conss[c]);
         assert(imagecons != NULL);
         SCIP_CALL( SCIPhashmapInsert(invconssmap, imagecons, conss[c]) );
      }
      SCIP_CALL( SCIPhashmapRemoveAll(iis->conssmap) );
   }

   /* get relevant parameters */
   SCIP_CALL( SCIPgetRealParam(scip, "iis/time", &timelim) );
   SCIP_CALL( SCIPgetLongintParam(scip, "iis/nodes", &nodelim) );
   SCIP_CALL( SCIPgetBoolParam(scip, "iis/removebounds", &removebounds) );
   SCIP_CALL( SCIPgetBoolParam(scip, "iis/silent", &silent) );

   /* make irreducible by running the deletion filter with singleton batches */
   SCIP_CALL( deletionFilterBatch(iis, timelim, nodelim, removebounds, silent,
         DEFAULT_TIMELIMPERITER, DEFAULT_NODELIMPERITER, TRUE, 1, 0.0, DEFAULT_MAXBATCHSIZE, DEFAULT_MAXRELBATCHSIZE,
         DEFAULT_BATCHINGFACTOR, DEFAULT_BATCHINGOFFSET, DEFAULT_BATCHUPDATEINTERVAL, DEFAULT_DETECTCOMPONENTS,
         DEFAULT_COMPONENTMINSIZE, &alldeletionssolved) );
   if( alldeletionssolved && SCIPiisGetTime(iis) < timelim && ( nodelim == -1 || SCIPiisGetNNodes(iis) < nodelim ) )
      SCIPiisSetSubscipIrreducible(iis, TRUE);

   /* recreate main constraints hashmap */
   if( isstandalone )
   {
      assert(invconssmap != NULL);
      nconss = SCIPgetNOrigConss(scip);
      conss = SCIPgetOrigConss(scip);
      for( c = 0; c < nconss; ++c )
      {
         imagecons = SCIPhashmapGetImage(invconssmap, conss[c]);
         assert(imagecons != NULL);
         SCIP_CALL( SCIPhashmapInsert(iis->conssmap, imagecons, conss[c]) );
      }
      SCIPhashmapFree(&invconssmap);
   }

   return SCIP_OKAY;
}
