/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   compr_weakcompr.c
 * @brief  weakcompr tree compression
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/compr_weakcompr.h"
#include "scip/compress.h"
#include "scip/reopt.h"
#include "scip/type_reopt.h"

#define COMPR_NAME             "weakcompr"
#define COMPR_DESC             "reduce the search frontier to k+1 or max{2, |C|+1} nodes."
#define COMPR_DISPCHAR         'W'
#define COMPR_PRIORITY         1000
#define COMPR_MINDEPTH         0
#define COMPR_MINNNODES        50
#define COMPR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

#define DEFAUL_MEM_REPR        10
/*
 * Data structures
 */

/** tree compression data */
struct SCIP_ComprData
{
   /* representative data */
   SCIP_REOPTNODE**      representatives;         /**< list of representatives */
   int                   nrepresentatives;        /**< number of representatives */
   int                   allocmemrepr;            /**< allocated memory for representatives */

   /* statistics */
   SCIP_Real             rate;                    /**< compression rate */
   SCIP_Real             loi;                     /**< loss of information */
   int                   nnodes;                  /**< number of nodes after compression */

   /* parameter */
   SCIP_Bool             convertconss;            /**< convert added logic-or constraints of size k into k nodes */
   int                   min_leaves;              /**< minimal number of nodes to compress */
   int                   numstartnodes;           /**< size of the compression, -1: log(#variables) */
};


/*
 * Local methods
 */

/**
 * subroutine for quicksort.
 */
static
int partition(
   SCIP*                 scip,                    /**< SCIP data structure */
   int*                  childids,                /**< array of child ids */
   int                   left,                    /**< first position */
   int                   right                    /**< last position */
   )
{
   int pivot;
   int i;
   int j;
   int t;

   pivot = left;
   i = left;
   j = right+1;

   while( TRUE )
   {
      do
      {
         ++i;
      } while( i <= right && SCIPisGT(scip, SCIPgetReoptNodeLb(scip, childids[i]), SCIPgetReoptNodeLb(scip, childids[pivot])) );

      do
      {
         --j;
      } while( j > 0 && SCIPisLT(scip, SCIPgetReoptNodeLb(scip, childids[j]), SCIPgetReoptNodeLb(scip, childids[pivot])) );

      if( i >= j )
         break;

      t = childids[i];
      childids[i] = childids[j];
      childids[j] = t;
   }

   t = childids[left];
   childids[left] = childids[j];
   childids[j] = t;

   return j;
}

/**
 * quicksprt implementation.
 * sort the ids of child nodes by there dual bound of the last iteration
 */
static
void sortIDs(
   SCIP*                 scip,                    /**< SCIP data structure */
   int*                  childids,                /**< array of child ids */
   int                   left,                    /**< first position */
   int                   right                    /**< last position */
   )
{
   int j;

   if( left < right )
   {
      j = partition(scip, childids, left, right);
      sortIDs(scip, childids, left, j-1);
      sortIDs(scip, childids, j+1, right);
   }
}

/**
 * check of enough memory is allocated and reallocate of necessary.
 */
static
SCIP_RETCODE checkMemSize(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_COMPRDATA*       comprdata,               /**< compression data */
   int                   nrepresentatives         /**< number of representatives */
   )
{
   assert(scip != NULL);
   assert(comprdata != NULL);

   if( comprdata->allocmemrepr < nrepresentatives )
   {
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &comprdata->representatives, comprdata->allocmemrepr, nrepresentatives) );
      comprdata->allocmemrepr = nrepresentatives;
   }

   return SCIP_OKAY;
}

/**
 * try to find a good representation
 */
static
SCIP_RETCODE constructCompression(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_COMPR*           compr,                   /**< compression method */
   SCIP_COMPRDATA*       comprdata,               /**< compression data */
   SCIP_RESULT*          result                   /**< result pointer */
   )
{
   SCIP_NODE* currentnode;
   SCIP_VAR**** conss_var;
   SCIP_VAR*** vars;
   SCIP_Real*** conss_val;
   SCIP_Real** vals;
   SCIP_BOUNDTYPE** bounds;
   int** conss_nvars;
   int* leaveids;
   int* nconss;
   int* nvars;
   int nids;
   int nleaveids;
   int depth;
   int k;
   int size;

   assert(scip != NULL);
   assert(comprdata != NULL);

   *result = SCIP_DIDNOTRUN;

   depth = 0;

   /* calculate the size of the representation */
   size = comprdata->numstartnodes == -1 ? log10(SCIPgetNBinVars(scip))/log10(2.0) : comprdata->numstartnodes;
   size = 1; /* @ todo: fix this */

   currentnode = SCIPgetStage(scip) <= SCIP_STAGE_PRESOLVED ? NULL : SCIPgetCurrentNode(scip);

   if( SCIPgetStage(scip) <= SCIP_STAGE_PRESOLVED )
      nleaveids = SCIPgetReoptNLeaves(scip, currentnode);
   else
   {
      assert(currentnode != NULL);
      nleaveids = SCIPgetReoptNLeaves(scip, currentnode);
      depth = SCIPnodeGetDepth(currentnode);
   }

   if( size > 1 ) /* && !comprdata->convertconss ) */
      return SCIP_OKAY;

   SCIPdebugMessage(">> start <%s> at node %llu (nleaves: %d, depth: %d)\n", COMPR_NAME,
         SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED ? 0 : SCIPnodeGetNumber(SCIPgetCurrentNode(scip)),
         nleaveids, depth);

   if( SCIPcomprGetMinnnodes(compr) > nleaveids )
   {
      SCIPdebugMessage("-> skip compression (min. leaves = %d)\n", SCIPcomprGetMinnnodes(compr));
      return SCIP_OKAY;
   }

   if( SCIPcomprGetMindepth(compr) > depth )
   {
      SCIPdebugMessage("-> skip compression (min. depth = %d)\n", SCIPcomprGetMindepth(compr));
      return SCIP_OKAY;
   }

   if( size > nleaveids )
   {
      SCIPdebugMessage("-> skip compression (k = %d, nleaves = %d)\n", size, nleaveids);
      return SCIP_OKAY;
   }

   SCIPdebugMessage("-> try compression with %d node(s)\n", comprdata->numstartnodes);

   *result = SCIP_DIDNOTFIND;

   if( comprdata->convertconss && size > 1)
   {
      printf("WARNING: converting constraints into nodes is currently not implemented. run with k = 1");
      size = 1;
   }

   /* collect the nodes to compress */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &leaveids, nleaveids) );

   SCIP_CALL( SCIPgetReoptLeaveIDs(scip, currentnode, leaveids, nleaveids, &nids) );
   assert(nids == nleaveids);

   /* sort the ids */
   sortIDs(scip, leaveids, 0, nleaveids-1);

   /* allocate memory to store the old tree */
   SCIP_CALL( SCIPallocBufferArray(scip, &vars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &vals, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &bounds, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conss_var, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conss_val, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &conss_nvars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nvars, size) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nconss, size) );

   /* get data of nodes */
   for(k = 0; k < size; k++)
   {
      int mem_vars;
      int mem_conss;
      int nvars2;
      int nafterdualvars;

      mem_vars = SCIPgetNBinVars(scip);

      /* allocate memory */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars[k], mem_vars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &vals[k], mem_vars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &bounds[k], mem_vars) );

      /* get the branching path */
      SCIPgetReoptnodePath(scip, leaveids[k], vars[k], vals[k], bounds[k], mem_vars, &nvars2, &nafterdualvars);

      /* reallocate memory */
      if( mem_vars < nvars2 + nafterdualvars )
      {
         mem_vars = nvars2 + nafterdualvars;
         SCIP_CALL( SCIPreallocBufferArray(scip, &vars[k], mem_vars) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &vals[k], mem_vars) );
         SCIP_CALL( SCIPreallocBufferArray(scip, &bounds[k], mem_vars) );

         /* get the branching path */
         SCIPgetReoptnodePath(scip, leaveids[k], vars[k], vals[k], bounds[k], mem_vars, &nvars2, &nafterdualvars);
      }

      nvars[k] = nvars2 + nafterdualvars;

      /* get the constraints */
      mem_conss = SCIPgetReoptnodeNConss(scip, leaveids[k]);

      SCIP_CALL( SCIPallocBufferArray(scip, &conss_var[k], mem_conss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &conss_val[k], mem_conss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &conss_nvars[k], mem_conss) );

      SCIPgetReoptnodeConss(scip, leaveids[k], conss_var[k], conss_val[k], mem_conss, &nconss[k], conss_nvars[k]);
      assert(mem_conss == nconss[k]);

#ifdef SCIP_DEBUG
      for(c = 0; c < mem_conss; c++)
         assert(conss_nvars[k][c] <= SCIPgetNBinVars(scip));
#endif

      SCIPdebugMessage("-> use node at id %d, %d vars, %d conss, lowerbound = %.g\n", leaveids[k], nvars[k],
            SCIPgetReoptnodeNConss(scip, leaveids[k]), SCIPgetReoptNodeLb(scip, leaveids[k]));
   }

   /* @ todo: convert constraints into nodes */

   /* perform the compression */
   if( size == 1 )
   {
      int pos_repr_fix;
      int r;

      assert(comprdata->nrepresentatives == 0);

      pos_repr_fix = 1;

      /* calculate the number of representatives */
      comprdata->nrepresentatives = (nvars[0] > 0 ? 2 : 1);
      comprdata->nrepresentatives += nconss[0];

      /* check memory size */
      SCIP_CALL( checkMemSize(scip, comprdata, comprdata->nrepresentatives) );
      assert(comprdata->nrepresentatives <= comprdata->allocmemrepr);

      /* initialize the representatives */
      for(r = 0; r < comprdata->nrepresentatives; r++)
      {
         SCIP_CALL( SCIPallocMemory(scip, &comprdata->representatives[r]) );
      }

      SCIPinitilizeRepresentation(comprdata->representatives, comprdata->nrepresentatives);

      /* create 2 candidates for the fixed variables */
      if( nvars[0] >= 1 )
      {
         int v;

         assert(pos_repr_fix < comprdata->nrepresentatives);

         /* create a representative at position 1 with fixed branching path */
         assert(SCIPreoptnodeGetNVars(comprdata->representatives[pos_repr_fix]) == 0);
         for(r = pos_repr_fix; r < comprdata->nrepresentatives; r++)
         {
            /* copy the branching path to all representatives */
            assert(comprdata->representatives[r] != NULL);

            for(v = 0; v < nvars[0]; v++)
            {
               SCIP_CALL( SCIPaddReoptnodeVar(scip, comprdata->representatives[r], vars[0][v],
                     vals[0][v], SCIPisFeasEQ(scip, vals[0][v], 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER) );
            }
         }

         /* create a representative at position 0 with an added constraint corresponding
          * to the branching path of the node*/
         assert(comprdata->representatives[pos_repr_fix-1] != NULL);
         SCIP_CALL( SCIPaddReoptnodeCons(scip, comprdata->representatives[pos_repr_fix-1], vars[0], vals[0], nvars[0], REOPT_CONSTYPE_STRBRANCHED) );

      }

      assert(0 <= pos_repr_fix && pos_repr_fix < comprdata->nrepresentatives);

      /* create nconss[0] nodes for the added constraints */
      for(k = 0; k < nconss[0]; k++)
      {
         int v;

         assert(pos_repr_fix < comprdata->nrepresentatives);

         /* create a node with fixed bounds corresponding to constraint at position k */

         /* fix the branching path */
         for(v = 0; v < conss_nvars[0][k]; v++)
         {
            SCIP_CALL( SCIPaddReoptnodeVar(scip, comprdata->representatives[pos_repr_fix], conss_var[0][k][v], conss_val[0][k][v],
                  SCIPisFeasEQ(scip, conss_val[0][k][v], 1) ? SCIP_BOUNDTYPE_LOWER : SCIP_BOUNDTYPE_UPPER) );
         }

         /* add this constraint to all further representatives */
         for(r = pos_repr_fix+1; r < comprdata->nrepresentatives; r++)
         {
            SCIP_CALL( SCIPaddReoptnodeCons(scip, comprdata->representatives[r], conss_var[0][k], conss_val[0][k],
                  conss_nvars[0][k], REOPT_CONSTYPE_STRBRANCHED) );
         }

         pos_repr_fix++;
      }
   }
   else if( size > 1 )
   {
      /* @ todo*/
   }

   *result = SCIP_SUCCESS;

   SCIPdebugMessage("-> found representation of size %d.\n", comprdata->nrepresentatives);

   /* free memory */
   for(k = size-1; k >= 0; k--)
   {
      SCIPfreeBufferArray(scip, &conss_nvars[k]);
      SCIPfreeBufferArray(scip, &conss_val[k]);
      SCIPfreeBufferArray(scip, &conss_var[k]);
      SCIPfreeBufferArray(scip, &bounds[k]);
      SCIPfreeBufferArray(scip, &vals[k]);
      SCIPfreeBufferArray(scip, &vars[k]);
   }

   SCIPfreeBufferArray(scip, &nconss);
   SCIPfreeBufferArray(scip, &nvars);
   SCIPfreeBufferArray(scip, &conss_nvars);
   SCIPfreeBufferArray(scip, &conss_val);
   SCIPfreeBufferArray(scip, &conss_var);
   SCIPfreeBufferArray(scip, &bounds);
   SCIPfreeBufferArray(scip, &vals);
   SCIPfreeBufferArray(scip, &vars);

   SCIPfreeBlockMemoryArray(scip, &leaveids, nleaveids);

   return SCIP_OKAY;
}

/**
 * apply the stored representation to the reopttree
 */
static
SCIP_RETCODE applyCompression(
   SCIP*                 scip,                    /**< SCIP data structure */
   SCIP_COMPR*           compr,                   /**< compression method */
   SCIP_COMPRDATA*       comprdata,               /**< compression data */
   SCIP_RESULT*          result                   /**< result pointer */
   )
{
   SCIP_Bool success;
   int old_nnodes;
   int r;

   assert(scip != NULL);
   assert(compr != NULL);
   assert(comprdata != NULL);

   *result = SCIP_DIDNOTRUN;

   if( comprdata->nrepresentatives == 0 )
      return SCIP_OKAY;

   /* set references to the root node */
   for(r = 0; r < comprdata->nrepresentatives; r++)
      SCIPreoptnodeSetParentID(comprdata->representatives[r], 0);

   old_nnodes = SCIPgetNReoptNodeIDs(scip, NULL)+1;

   success = FALSE;
   SCIP_CALL( SCIPsetReoptCompression(scip, comprdata->representatives, comprdata->nrepresentatives, &success) );

   if( success )
   {
      comprdata->nnodes = SCIPgetNReoptNodeIDs(scip, NULL)+1;
      comprdata->rate = ((SCIP_Real) comprdata->nnodes)/old_nnodes;

      SCIPcomprUpdateRate(compr, comprdata->rate);
      SCIPcomprUpdateNNodes(compr, comprdata->nnodes);
      *result = SCIP_SUCCESS;
   }
   return SCIP_OKAY;
}

/*
 * Callback methods of tree compression
 */

/** copy method for tree compression plugins (called when SCIP copies plugins) */
#define comprCopyWeakcompr NULL

/** destructor of tree compression to free user data (called when SCIP is exiting) */
static
SCIP_DECL_COMPRFREE(comprFreeWeakcompr)
{
   SCIP_COMPRDATA* comprdata;

   assert(scip != NULL);
   assert(compr != NULL);

   comprdata = SCIPcomprGetData(compr);
   assert(comprdata != NULL);

   SCIPfreeMemory(scip, &comprdata);
   SCIPcomprSetData(compr, NULL);

   return SCIP_OKAY;
}


/** initialization method of tree compression (called after problem was transformed) */
static
SCIP_DECL_COMPRINIT(comprInitWeakcompr)
{
   SCIP_COMPRDATA* comprdata;

   assert(scip != NULL);
   assert(compr != NULL);

   comprdata = SCIPcomprGetData(compr);
   assert(comprdata != NULL);

   SCIPdebugMessage(">> initializing <%s>\n", COMPR_NAME);

   comprdata->allocmemrepr = DEFAUL_MEM_REPR;
   comprdata->nrepresentatives = 0;
   comprdata->rate = 0.0;
   comprdata->loi = 0.0;
   comprdata->nnodes = 0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &comprdata->representatives, comprdata->allocmemrepr) );

   return SCIP_OKAY;
}


/** deinitialization method of tree compression (called before transformed problem is freed) */
#define comprExitWeakcompr NULL

/** solving process initialization method of tree compression (called when branch and bound process is about to begin) */
#define comprInitsolWeakcompr NULL

/** solving process deinitialization method of tree compression (called before branch and bound process data is freed) */
#define comprExitsolWeakcompr NULL


/** execution method of tree compression */
static
SCIP_DECL_COMPREXEC(comprExecWeakcompr)
{
   SCIP_COMPRDATA* comprdata;

   assert(SCIPcomprIsInitialized(compr));

   comprdata = SCIPcomprGetData(compr);
   assert(comprdata != NULL);

   /* try to find a representation */
   SCIP_CALL( constructCompression(scip, compr, comprdata, result) );

   assert(*result == SCIP_DIDNOTRUN || *result == SCIP_DIDNOTFIND || *result == SCIP_SUCCESS);

   /* apply the representation, if some was found */
   if( *result == SCIP_SUCCESS )
   {
      SCIP_CALL( applyCompression(scip, compr, comprdata, result) );
      assert(*result == SCIP_DIDNOTRUN || *result == SCIP_SUCCESS);

      SCIPdebugMessage("->%s apply compression.\n", *result == SCIP_DIDNOTRUN ? " did not" : "");
   }

   return SCIP_OKAY;
}


/*
 * tree compression specific interface methods
 */

/** creates the weakcompr tree compression and includes it in SCIP */
SCIP_RETCODE SCIPincludeComprWeakcompr(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_COMPRDATA* comprdata;
   SCIP_COMPR* compr;

   /* create weakcompr tree compression data */
   SCIP_CALL( SCIPallocMemory(scip, &comprdata) );
   assert(comprdata != NULL);

   /* include tree compression */
   SCIP_CALL( SCIPincludeComprBasic(scip, &compr,COMPR_NAME, COMPR_DESC,
         COMPR_DISPCHAR, COMPR_PRIORITY, COMPR_MINDEPTH, COMPR_MINNNODES,
         COMPR_USESSUBSCIP, comprExecWeakcompr, comprdata) );

   assert(compr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetComprInit(scip, compr, comprInitWeakcompr) );
   SCIP_CALL( SCIPsetComprFree(scip, compr, comprFreeWeakcompr) );

   /* add weakcompr tree compression parameters */
   SCIP_CALL( SCIPaddBoolParam(scip, "compression/"COMPR_NAME"/convertconss", "convert constraints into nodes", &comprdata->convertconss, FALSE, FALSE, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "compression/"COMPR_NAME"/numstartnodes", "size of the compression (-1: k = log(#variables))", &comprdata->numstartnodes, FALSE, 5, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
