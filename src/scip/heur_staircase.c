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
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_staircase.c
 * @brief  start heuristic for decompositions with staircase structure
 * @author Christine Tawfik
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/heur_staircase.h"
#include "scip/decomp.h"

#define HEUR_NAME             "staircase"
#define HEUR_DESC             "start heuristic for decompositions with staircase structure"
#define HEUR_DISPCHAR         'D'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      TRUE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopyStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopyStaircase NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_HEURFREE(heurFreeStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreeStaircase NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitStaircase NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitStaircase NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolStaircase NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolStaircase)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolStaircase NULL
#endif

/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecStaircase)
{  /*lint --e{715}*/

   SCIP_DECOMPSTORE* decompstore;
   SCIP_DECOMP** decomps;
   int ndecomps;
   int d;

   /* query if there is a user decomposition (see methods in decomp.h) */
   decompstore = SCIPgetDecompstore(scip);
   assert(decompstore != NULL);

   decomps = SCIPdecompstoreGetDecomps(decompstore);
   ndecomps = SCIPdecompstoreGetNDecomps(decompstore);

   /* staircase heuristic requires a decomposition to work with */
   if( ndecomps == 0 )
      return SCIP_OKAY;

   SCIPdebugMsg(scip, "Staircase heuristic works on %d decomposition%s\n", ndecomps, ndecomps == 1 ? "" : "s");



   /* author bzfhende
    *
    * TODO check if the decomposition has the appropriate staircase structure
    *    - no linking constraints
    *    - all linking variables link exactly 2 named blocks of constraints
    *    - the block connectivity graph is a path
    *
    *    if all the above assumptions are true, we should have an order of the blocks for the rolling horizon (if graph
    *    is a path, there are exactly two orders depending on which degree 1 node we start from)
    */

   /* todo comment */
   for( d = 0; d < ndecomps; ++d )
   {
      SCIP_DECOMP* decomp;
      SCIP_CONS** conss;
      SCIP_VAR** vars;
      SCIP_VAR** consvars;
      SCIP_VAR** consvars2;
      SCIP_DIGRAPH* blocklinkingvargraph;
      SCIP_DIGRAPH* blockgraph;
      SCIP_QUEUE* queue;
      SCIP_Bool success, success2;
      SCIP_Bool varblocksflag = TRUE;
      SCIP_Bool cycleflag;
      int* conslabels;
      int* varlabels;
      int* linkvaridx;
      int* succnodes;
      int* visited;
      int* parent;
      int linkingvarlabel1 = -1;
      int linkingvarlabel2 = -1;
      int nlinkingconss = 0;
      int nlinkingvars = 0;
      int ncontlinkingvars = 0;
      int nconss;
      int nvars;
      int nconsvars, nconsvars2;
      int ndiscvars = 0;
      int ndiscblks = 0;
      int i;
      int v;
      int j;
      int l;
      int n;
      int root;
      int succ1;
      int succ2;
      int nsucc;
      int nnodes;

      nconss = SCIPgetNConss(scip);
      conss = SCIPgetConss(scip);
      nvars = SCIPgetNVars(scip);
      vars = SCIPgetVars(scip);

      decomp = decomps[d];

      SCIP_CALL( SCIPallocBufferArray(scip, &conslabels, nconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &varlabels, nvars) );
      SCIP_CALL( SCIPallocBufferArray(scip, &linkvaridx, nvars) );

      for( i = 0; i < nvars; ++i )
    	  linkvaridx[i] = -1;

      SCIPdecompGetConsLabels(decomp, conss, conslabels, nconss);

      SCIPdecompGetVarsLabels(decomp, vars, varlabels, nvars);

      SCIPdebugMsg(scip, "This decomposition has %d blocks.\n", SCIPdecompGetNBlocks(decomp));

      /* count number of linking constraints */
      for( i = 0; i < nconss; ++i )
      {
         if( conslabels[i] == SCIP_DECOMP_LINKCONS )
         {
            ++nlinkingconss;
            SCIPdebugMsg(scip, "Constraint %s is a linking %sconstraint %s\n", SCIPconsGetName(conss[i]), SCIPconsIsConflict(conss[i]) == TRUE ? "(conflict) " : "", SCIPconsIsChecked(conss[i]) == TRUE ? " and it is CHECKED." : "");

            SCIP_CALL( SCIPprintCons(scip, conss[i], NULL) );
         }
      }

      SCIPdebugMsg(scip, "Decomposition %d has %d linking constraints\n", d, nlinkingconss);

      /* count number of linking variables and
       *
       * create a unique mapping of all linking variables to 0,..,nlinkingvars -1 and store it in array linkvaridx */

      for( v = 0; v < nvars; ++v )
      {
         if( varlabels[v] == SCIP_DECOMP_LINKVAR )
         {
            linkvaridx[v] = nlinkingvars;
            assert(linkvaridx[v] == nlinkingvars);
            assert(SCIPvarGetProbindex(vars[v]) == v);
            if( SCIPvarGetType(vars[v]) == SCIP_VARTYPE_CONTINUOUS )
               ++ncontlinkingvars;
            ++nlinkingvars;
         }
      }
      SCIPdebugMsg(scip, "Decomposition %d has %d linking variables (%d continuous)\n", d, nlinkingvars, ncontlinkingvars);

      /* author bzftawfi
       *
       * TODO construct a block DS and a graph where:
       * Nodes corresponds to blocks
       * Edges correspond to linking variables
       *
       * */

      /* author bzftawfi
       *
       * TODO check the following:
       * The graph should contain no cycle of more than 2 edges;
       * 2-edge cycles are acceptable since the nodes could be connected by several edges
       * no cycles detection entails that only consecutive blocks could be linked
       * The graph should be connected
       *
       * Starting from the origin/destination node, query the next reachable node and save its information -> constructing the components of the sub-MIP
       *
       * Required information:
       * identifying the linking variables associated with a certain constraint and whether its coef. is nnz
       * identifying that a linking variables appears in exactly two blocks
       * if the previous step is successful, get the number of its up and down locks
       *
       * */

      /* future TODO:
       *
       * act upon the existence of linking constraints; special treatment for the constraints that do not violate the staircase structure
       * act upon the negative check of the negative check of the connectedness of linking variables to exactly 2 blocks
       * Currently: aggregation is not allowed while running SCIP
       * Activate a check for all blocks - within the rolling horizon - of the avtive/inactive variables; if only inactive variables exist, this block could be skipped as a stage
       *
       * */


      /* Check 1: using a bipartite graph, check the connectedness of each linking variable to 2 blocks
       * Checking procedure in O(no. of nonzeros in constraint matrix)  */

      /* create a bipartite graph consisting of 2 sets of nodes corresponding to the decomposition blocks and the linking variables, respectively */
      SCIP_CALL( SCIPcreateDigraph(scip, &blocklinkingvargraph, SCIPdecompGetNBlocks(decomp) + nlinkingvars) );

      for( i = 0; i < nconss; ++i )
      {
    	  int requiredsize;
    	  /* linking constraints are skipped in this checking step */
    	  if( conslabels[i] == SCIP_DECOMP_LINKCONS )
    		  continue;
    	  SCIP_CALL(SCIPgetConsNVars(scip, conss[i], &nconsvars, &success) );
    	  if (!success)
    		  continue;
    	  SCIP_CALL(SCIPallocBufferArray(scip, &consvars, nconsvars) );
    	  SCIP_CALL(SCIPgetConsVars(scip, conss[i], consvars, nconsvars, &success));
    	  if(!success)
    		  continue;

    	  /* re-transform given variables to active variables*/
    	  SCIP_CALL( SCIPgetActiveVars(scip, consvars, &nconsvars, nconsvars, &requiredsize) );
    	  assert(requiredsize <= nconsvars);
    	  SCIPdecompGetVarsLabels(decomp, consvars, varlabels, nconsvars);


		  for( j = 0; j < nconsvars; ++j )
		  {
			  assert(consvars[j] != NULL);

			  if( varlabels[j] == SCIP_DECOMP_LINKVAR )
			  {
				  int linkingvarnodeidx = linkvaridx[SCIPvarGetProbindex(consvars[j])];

				  SCIPdigraphAddArcSafe(blocklinkingvargraph, SCIPdecompGetNBlocks(decomp) + linkingvarnodeidx, conslabels[i] - 2, NULL);
				  SCIPdigraphAddArcSafe(blocklinkingvargraph, conslabels[i] - 2, SCIPdecompGetNBlocks(decomp) + linkingvarnodeidx, NULL);
			  }
		  }
      }
      SCIPdebugMsg(scip, "Total no. of edges are equal %d and it should be equal %d \n", SCIPdigraphGetNArcs(blocklinkingvargraph)/2, nlinkingvars*2);


      success = TRUE;
	  for( n = SCIPdecompGetNBlocks(decomp); n < SCIPdigraphGetNNodes(blocklinkingvargraph); ++n )
	  {
		  nsucc = (int)SCIPdigraphGetNSuccessors(blocklinkingvargraph, n);
		  if( nsucc != 2 && nsucc != 0)
		  {
			  SCIPdebugMsg(scip, "Node %d has %d successors \n", n, nsucc);
			  success = FALSE;
		  }
		  if (nsucc == 0)
			  ++ndiscvars;
	  }
	  for( n = 0; n < SCIPdecompGetNBlocks(decomp); ++n)
	  {
		  nsucc = (int) SCIPdigraphGetNSuccessors(blocklinkingvargraph, n);
		  if (nsucc == 0)
			  ++ndiscblks;
	  }
      SCIPdebugMsg(scip, "Check #1: Connectedness of linking variables to 2 blocks: %s\n", success == FALSE ? " FAIL " : " SUCCESS ");
      SCIPdebugMsg(scip, "There are %d disconnected linking variables and %d disconnected blocks. \n", ndiscvars, ndiscblks);

      /* Check 2: Using the information from the previous bi-partite graph, build a block-decomposition graph
       * and check for CONNECTIVITY and CYCLES
       */

      /* create a block graph consisting of nodes corresponding to the decomposition blocks and arcs corresponding to linking variables */
      SCIP_CALL( SCIPcreateDigraph(scip, &blockgraph, SCIPdecompGetNBlocks(decomp)) );

      SCIPdebugMsg(scip, "Building the block decomposition graph... \n");

      for( n = SCIPdecompGetNBlocks(decomp); n < SCIPdigraphGetNNodes(blocklinkingvargraph); ++n )
      {
    	  nsucc = (int) SCIPdigraphGetNSuccessors(blocklinkingvargraph, n);
    	  //SCIP_CALL( SCIPallocBufferArray(scip, &succnodes, nsucc) );
    	  succnodes = (int*) SCIPdigraphGetSuccessors(blocklinkingvargraph, n);
    	  for ( succ1 = 0; succ1 < nsucc; ++succ1 )
    	  {
    		  for( succ2 = 0; succ2 < nsucc; ++succ2 )
    		  {
    			  if( succnodes[succ1] != succnodes[succ2] )
    				  SCIPdigraphAddArcSafe(blockgraph, succnodes[succ1], succnodes[succ2], NULL);
    		  }
    	  }
    	  //SCIPfreeBufferArray(scip, &succnodes);
      }

      SCIPdebugMsg(scip, "The block decomposition graph contains %d edges. \n", SCIPdigraphGetNArcs(blockgraph)/2);

      /*
       * Checking for the graph connectedness and the existence of cycles using a BFS implementation
       */

      //SCIPdigraphAddArcSafe(blockgraph, 8, 6, NULL);

      /* Mark all vertices as unvisited (0 - unvisited; 1 - visited) */
      SCIP_CALL( SCIPallocBufferArray(scip, &visited, SCIPdigraphGetNNodes(blockgraph)) );
      for( n = 0; n < SCIPdigraphGetNNodes(blockgraph); ++n )
      {
         visited[n] = 0;
      }

      /* Create an Int array to keep track of the nodes' parents in the search tree */
      SCIP_CALL( SCIPallocBufferArray(scip, &parent, SCIPdigraphGetNNodes(blockgraph)) );

      /* Create a queue of the visited vertices in a BFS order (initially having the root) */
      /* todo: make the root choice random (?) */
      root = SCIPdigraphGetNNodes(blockgraph) / 2;
      SCIP_CALL( SCIPqueueCreate(&queue, SCIPdigraphGetNNodes(blockgraph), 2.0) );
      SCIP_CALL( SCIPqueueInsertUInt(queue, root) );
      visited[root] = 1;
      parent[root] = -1; /* the root has no parent */
      nnodes = 0;
      cycleflag = FALSE;

      while( !SCIPqueueIsEmpty(queue) )
      {
    	  ++nnodes;
    	  n = SCIPqueueRemoveUInt(queue);
    	  nsucc = (int) SCIPdigraphGetNSuccessors(blockgraph, n);
    	  //SCIP_CALL( SCIPallocBufferArray(scip, &succnodes, nsucc) );
    	  succnodes = (int*) SCIPdigraphGetSuccessors(blockgraph, n);
    	  for ( succ1 = 0; succ1 < nsucc; ++succ1 )
    	  {
    		  if(visited[succnodes[succ1]] == 0)
    		  {
    			  visited[succnodes[succ1]] = 1;
    			  parent[succnodes[succ1]] = n;
    			  SCIP_CALL( SCIPqueueInsertUInt(queue, succnodes[succ1]) );
    		  }
    		  else
    			  if(visited[succnodes[succ1]] == 1 && parent[n] != succnodes[succ1])
    				  cycleflag = TRUE;
    	  }
    	  //SCIPfreeBufferArray(scip, &succnodes);

      }

      SCIPdebugMsg(scip, "Check #2: The block decomposition graph is%s and contains %s cycle. \n",  nnodes == SCIPdigraphGetNNodes(blockgraph) ? " CONNECTED" : " NOT CONNECTED", cycleflag == TRUE ? "A" : "NO");

      /* get the number of connected components */
      SCIP_CALL( SCIPdigraphComputeUndirectedComponents(blockgraph, -1, NULL, NULL) );
      int ncomponents = SCIPdigraphGetNComponents(blockgraph);
      SCIPdebugMsg(scip, "The number of connected components = %d. \n", ncomponents);


      SCIPdigraphFree(&blocklinkingvargraph);
      SCIPdigraphFree(&blockgraph);
      SCIPqueueFree(&queue);
      SCIPfreeBufferArray(scip, &linkvaridx);
      SCIPfreeBufferArray(scip, &varlabels);
      SCIPfreeBufferArray(scip, &conslabels);
      SCIPfreeBufferArray(scip, &visited);
      SCIPfreeBufferArray(scip, &parent);
   }

   SCIPerrorMessage("method of staircase primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   /* author bzfhende
    *
    * TODO create or query a start solution
    *
    * - incumbent solution, if available
    * - zero solution (always available, perhaps infeasible)
    * - rounded LP relaxation solution (rounding based on locks or fractionality)
    * - partial solution (possible user input, see reader_mst.h)
    * - ...
    */

   /* author bzfhende
    *
    * TODO heuristic parameters
    *
    * - width of the horizon (previous blocks, subsequent blocks (with integer restrictions intact), subsequent relaxed blocks)
    *
    */

   /* author bzfhende
    *
    * TODO create a working solution by copying the start solution. This working solution is iteratively updated in the
    * main loop
    */


   /* author bzfhende
    *
    * TODO main loop: apply the rolling horizon algorithm.
    *
    * - in each iteration, formulate the corresponding sub-MIP in an auxiliary SCIP instance (have a look at the methods in scip_copy.c)
    * - solve this sub-MIP (have a look at heur_rins.c)
    * - update the working solution with the result from the sub-MIP
    */

   /* author bzfhende
    *
    * TODO if working solution is feasible, add it to SCIP
    */

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the staircase primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurStaircase(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create staircase primal heuristic data */
   heurdata = NULL;

   heur = NULL;

   /* include primal heuristic */
#if 0
   /* use SCIPincludeHeur() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyStaircase, heurFreeStaircase, heurInitStaircase, heurExitStaircase, heurInitsolStaircase, heurExitsolStaircase, heurExecStaircase,
         heurdata) );
#else
   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecStaircase, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyStaircase) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeStaircase) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitStaircase) );
   SCIP_CALL( SCIPsetHeurExit(scip, heur, heurExitStaircase) );
   SCIP_CALL( SCIPsetHeurInitsol(scip, heur, heurInitsolStaircase) );
   SCIP_CALL( SCIPsetHeurExitsol(scip, heur, heurExitsolStaircase) );
#endif

   /* add staircase primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
