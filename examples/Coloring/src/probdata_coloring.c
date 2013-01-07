/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_coloring.c
 * @brief  problem data for vertex coloring algorithm
 * @author Gerald Gamrath
 *
 * This file implements the problem data for the coloring algorithm.
 *
 * The problem data contains the original graph, preprocessing information, the preprocessed graph,
 * the array with the node-constraints, and an array with all stable sets and corresponding
 * variables.
 *
 * The preprocessing deletes nodes that have a lower degree than the size of a maximum clique.
 * Additionally, it also deletes nodes that have a dominated neighborhood. For further information,
 * look at the documentation for the method preprocessGraph().
 *
 * The deleted nodes and the relation between the nodes of the original graph and the nodes of the
 * preprocessed graph are stored in order to convert a solution of the preprocessed problem to a
 * solution for the original graph and vice versa.
 *
 * Each variable has a pointer of type SCIP_VARDATA* that is used in this case to store an integer
 * representing the number of the stable set. With the aid of this, the corresponding stable set can
 * be found in the array returned by COLORprobGetStableSets().  This array contains all stable sets
 * and is also used to check whether a stable set found by the pricer is really new. This can be
 * done by calling COLORprobStableSetIsNew(). All sets are sorted decreasingly with respect to the
 * indices of the nodes. New candidates should also be sorted that way.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "probdata_coloring.h"

#define EVENTHDLR_NAME         "probdatavardeleted"
#define EVENTHDLR_DESC         "event handler for variable deleted event"

struct SCIP_ProbData
{
   TCLIQUE_GRAPH*   graph;              /* the preprocessed graph */
   SCIP_CONS**      constraints;        /* array of added constraints */

   /* stable set / variable - information*/
   int**            stablesets;         /* array of stable sets */
   int*             stablesetlengths;   /* length of the array in stablesets */
   int              maxstablesets;      /* length of array stablesets */
   int              nstablesets;        /* number of stable sets saved in array stablesets */
   SCIP_VAR**       stablesetvars;      /* variables belonging to stable sets */

   /* preprocessing information */
   TCLIQUE_GRAPH*   oldgraph;           /* the original graph */
   int*             deletednodes;       /* array of nodes which were deleted during preprocessing */
   int*             new2oldnode;        /* new2oldnode[i] = j iff node i in the (preprocessed) graph corresponds to node j in the old graph*/
};


/*
 * Local methods
 */

/** 
 *  Preprocessing of the graph, using 2 methods in order to find redundant nodes
 *  that can be deleted and easily colored later.
 *
 *  Foundation of these methods is the computation of a maximum clique C with M nodes. 
 *  After this computation, the following two steps are repeated until no node was deleted
 *  in the last iteration:
 *  
 *  1: Low-Degree:
 *  Iterativly delete all nodes v in the graph G with degree d(v) < M ( don't delete nodes of C )
 *
 *  2: Dominated Neighbourhood:
 *  If the neighbourhood of one node v is part of the neighbourhood of another node w, v can 
 *  be deleted, since it can later get the same color as w.
 */
static 
SCIP_RETCODE preprocessGraph( 
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   SCIP_PROBDATA*   probdata;           /* the problemdata */
   SCIP_Bool        changed;            /* was the graph changed in the last round of preprocessing? */
   SCIP_Bool        dominates;          /* is the neighbourhood of one node dominated by the neigbourhood of another one?*/
   int*             maxcliquenodes;     /* pointer to store nodes of the maximum weight clique */
   int              nmaxcliquenodes;    /* number of nodes in the maximum weight clique */
   TCLIQUE_WEIGHT   maxcliqueweight;    /* weight of the maximum weight clique */
   TCLIQUE_STATUS   status;             /* status of clique-computation */
   TCLIQUE_GRAPH*   currgraph;          /* the current, not preprocessed graph in each step */
   int currnnodes;                      /* the current number of nodes in each step */ 
   int actnewnode;                      /* the number of nodes yet marked for beeing in the graph in the next round */
   int* newnodes;                       /* the nodes that will be in the graph in the next round */
   int* degrees;                        /* the degrees of the nodes */
   int myround;                         /* the number of the current round */
   int ndeletednodes;                   /* the total number of deleted nodes */
   int nnodesdeleteddegreethisround;    /* the number of nodes deleted due to low degree in the current round */
   int nnodesdeletedneighbourthisround; /* the number of nodes deleted due to neighbourhood in the current round */
   int*  firstedge;                     /* pointer for the edges in the graph */
   int*  lastedge;                      /* pointer for the edges in the graph */
   int i;
   int j;
   char opt;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   printf("\npreprocessing...\n");

   /* copy the old graph */
   tcliqueCreate(&currgraph);
   assert(currgraph != NULL);
   tcliqueAddNode(currgraph, tcliqueGetNNodes(probdata->oldgraph)-1, 0);
   for ( i = 0; i < tcliqueGetNNodes(probdata->oldgraph); i++ )
   {
      /* get adjacent nodes for node i */
      firstedge = tcliqueGetFirstAdjedge(probdata->oldgraph, i);
      lastedge = tcliqueGetLastAdjedge(probdata->oldgraph, i);
      while ( firstedge <= lastedge )
      {
         if ( *firstedge > i )
         {
            tcliqueAddEdge(currgraph, i, *firstedge);
         }
         firstedge++;
      }
   }
   tcliqueFlush(currgraph);
   currnnodes = tcliqueGetNNodes(probdata->oldgraph);

   /* get memory for array of deletednodes */  
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->deletednodes), COLORprobGetOriginalNNodes(scip)) );

   /* get memory for array new2oldnodes */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->new2oldnode), COLORprobGetOriginalNNodes(scip)) );  

   SCIP_CALL( SCIPallocBufferArray(scip, &newnodes, COLORprobGetOriginalNNodes(scip)) );
   SCIP_CALL( SCIPallocBufferArray(scip, &maxcliquenodes, COLORprobGetOriginalNNodes(scip)) );
   
   for ( i = 0; i < currnnodes; i++ )
   {
      /* set weights of the nodes to 1 */
      tcliqueChangeWeight(currgraph, i, 1);
      /* every node in the graph represents the same node in the original graph */
      probdata->new2oldnode[i] = i;
   } 

   /* compute maximum clique */
   tcliqueMaxClique(NULL, NULL, NULL, NULL, currgraph, NULL, NULL, maxcliquenodes,
      &nmaxcliquenodes, &maxcliqueweight, 0, 0, 50000, 0, INT_MAX, -1, NULL, &status);
   opt = ( status == TCLIQUE_OPTIMAL ? ' ' : '*' );
   printf("size of the maximum clique: %d%c \n", nmaxcliquenodes, opt);

   ndeletednodes = 0;
   nnodesdeleteddegreethisround = 1;
   nnodesdeletedneighbourthisround = 1;
   myround = 0;

   /* main loop */
   while ( (nnodesdeleteddegreethisround > 0) || (nnodesdeletedneighbourthisround > 0) ) 
   {
      myround++;
      nnodesdeleteddegreethisround = 0;
      nnodesdeletedneighbourthisround = 0;
      changed = TRUE;

      /* node degree deletion loop */
      while ( changed == TRUE )
      {
         changed = FALSE;
         actnewnode = 0;
         degrees = tcliqueGetDegrees(currgraph);
         for (i = 0; i < currnnodes; i++)
         {
            /* degree is low, node can be deleted */
            if ( (degrees[i] < nmaxcliquenodes ) 
               && (!COLORprobIsNodeInArray(probdata->new2oldnode[i], maxcliquenodes, nmaxcliquenodes)) )
            {
               probdata->deletednodes[ndeletednodes] = probdata->new2oldnode[i];
               changed = TRUE;
               nnodesdeleteddegreethisround++;
               ndeletednodes++;
            }
            /* node will be in the new graph, because degree is not low enought for deletion or it is in the maxClique */
            else
            {
               newnodes[actnewnode] = probdata->new2oldnode[i];
               actnewnode++;
            }
         }

         /* at least one node was deletet, create new graph ( tclique doesn't support deletion of nodes) */
         if ( changed )
         {
            assert(actnewnode+ndeletednodes == COLORprobGetOriginalNNodes(scip));
            /* create new current graph */
            tcliqueFree(&currgraph);
            tcliqueCreate(&currgraph);
            assert(currgraph != NULL);
            tcliqueAddNode(currgraph, actnewnode-1, 0);
            for ( i = 0; i < actnewnode; i++ )
            {
               /* get adjacent nodes for node newnodes[i] */
               firstedge = tcliqueGetFirstAdjedge(probdata->oldgraph, newnodes[i]);
               lastedge = tcliqueGetLastAdjedge(probdata->oldgraph, newnodes[i]);
               while ( firstedge <= lastedge )
               {
                  /* try to find a node in newnodes which corresponds 
                     to the node in the old graph, that is the end-node of the edge */
                  for ( j = i+1; j < actnewnode; j++ )
                  {
                     if ( *firstedge == newnodes[j] )
                     {
                        tcliqueAddEdge(currgraph, i, j);
                        break;
                     }
                  }
                  firstedge++;
               }
            }
            tcliqueFlush(currgraph);
            /* update currnnodes */
            currnnodes = tcliqueGetNNodes(currgraph);
            /* update new2oldnodes */
            for ( i = 0; i < actnewnode; i++ )
            {
               probdata->new2oldnode[i] = newnodes[i];
            }
            /* set the corresponding old node to -1 for all nodes not in the current graph (only for bug-finding) */
            for ( i = actnewnode; i < COLORprobGetOriginalNNodes(scip); i++ )
            {
               probdata->new2oldnode[i] = -1;
            }
         }
      } /* end node degree deletion loop */

      /* set changed to TRUE for getting into the while-loop */
      changed = TRUE;
      /* loop for finding dominated neighbourhoods */
      while ( changed == TRUE )
      {
         changed = FALSE;
         actnewnode = 0;
         /* i is the node which is checked for beeing dominated */
         for ( i = 0; i < currnnodes; i++ )
         {
            assert(!COLORprobIsNodeInArray(probdata->new2oldnode[i], probdata->deletednodes, ndeletednodes));

            /* i must be not in the clique and not yet deleted */
            if ( (!COLORprobIsNodeInArray(probdata->new2oldnode[i], maxcliquenodes, nmaxcliquenodes)) )
            {
               /* j is the node for which is checked whether it dominates i */
               for ( j = 0; j < currnnodes; j++ )
               {
                  /* i must be distinct from j, there must be no edge between i and j, 
                     j may not have been deleted due to degree in the last round */
                  if ( (j != i) && !tcliqueIsEdge(currgraph, i, j) 
                     && (!COLORprobIsNodeInArray(probdata->new2oldnode[j], probdata->deletednodes, ndeletednodes)) )
                     /** @todo only check nodes deleted in the last round */
                  {
                     /* check whether nodes adjacent to i are also adjacent to j <-> j dominates i */
                     dominates = TRUE;
                     /* get adjacent nodes for node i in currgraph */
                     firstedge = tcliqueGetFirstAdjedge(currgraph, i);
                     lastedge = tcliqueGetLastAdjedge(currgraph, i);
                     while ( firstedge <= lastedge )
                     {
                        /* check whether (j,firstedge) is in currgraph, if not, j doesn't dominate i */
                        if ( !tcliqueIsEdge(currgraph, j, *firstedge) )
                        {
                           dominates = FALSE;
                           break;
                        }
                        firstedge++;
                     }
                     if ( dominates )
                     {
                        probdata->deletednodes[ndeletednodes] = probdata->new2oldnode[i];
                        changed = TRUE;
                        ndeletednodes++;
                        nnodesdeletedneighbourthisround++;
                        break; /* for j, because we already now that i is dominated and deleted i */
                     }
                  }
               } /* end for j */
               
               /* if i is dominated by no other node and thus not deleted, 
                  put it into newnodes, so that it is in the next graph */
               if ( ndeletednodes == 0 || probdata->deletednodes[ndeletednodes-1] != probdata->new2oldnode[i])
               {
                  newnodes[actnewnode] = probdata->new2oldnode[i];
                  actnewnode++;
               }
            }
            /* if i is in the maxClique und was thus not deleted,
               put it into newnodes, so that it is in the next graph */
            else
            {
               newnodes[actnewnode] = probdata->new2oldnode[i];
               actnewnode++;
            }
         } /*end for i */

         /* at least one node was deletet, create new graph ( tclique doesn't support deletion of nodes) */
         if ( changed )
         {
            assert(actnewnode+ndeletednodes == COLORprobGetOriginalNNodes(scip));
            /* create new current graph */
            tcliqueFree(&currgraph);
            tcliqueCreate(&currgraph);
            assert(currgraph != NULL);
            tcliqueAddNode(currgraph, actnewnode-1, 0);
            for ( i = 0; i < actnewnode; i++ )
            {
               /* get adjacent nodes for node newnodes[i] */
               firstedge = tcliqueGetFirstAdjedge(probdata->oldgraph, newnodes[i]);
               lastedge = tcliqueGetLastAdjedge(probdata->oldgraph, newnodes[i]);
               while ( firstedge <= lastedge )
               {
                  /* try to find a node in newnodes which corresponds 
                     to the node in the old graph, that is the end-node of the edge */
                  for ( j = i+1; j < actnewnode; j++ )
                  {
                     if ( *firstedge == newnodes[j] )
                     {
                        tcliqueAddEdge(currgraph, i, j);
                        break;
                     }
                  }
                  firstedge++;
               }
            }
            tcliqueFlush(currgraph);

            /* update currnnodes */
            currnnodes = tcliqueGetNNodes(currgraph);

            /* update new2oldnodes */
            for ( i = 0; i < actnewnode; i++ )
            {
               probdata->new2oldnode[i] = newnodes[i];
            }

            /* set the corresponding old node to -1 for all nodes not in the current graph (only for bug-finding) */
            for ( i = actnewnode; i < COLORprobGetOriginalNNodes(scip); i++ )
            {
               probdata->new2oldnode[i] = -1;
            }
         }
      } /* end of loop for finding dominated neighbourhoods */
      
      printf("Round %d of preprocessing:\n", myround);
      printf("   deleted low degree vertices: %d\n", nnodesdeleteddegreethisround);
      printf("   deleted almost cliques:      %d\n", nnodesdeletedneighbourthisround);

   }

   for ( i = ndeletednodes; i < COLORprobGetOriginalNNodes(scip); i++ )
   {
      probdata->deletednodes[i] = -1;
   }

   printf("preprocessing overall deleted vertices: %d\n\n", ndeletednodes);
   /* copy preprocessed graph into problem data */
   probdata->graph = currgraph;

   SCIPfreeBufferArray(scip, &newnodes);
   SCIPfreeBufferArray(scip, &maxcliquenodes);

   return SCIP_OKAY;
}




/*
 * Callback methods of probdata
 */

/** transforms the problem */
static
SCIP_DECL_PROBTRANS(probtransColoring)
{
   int i;
   int j;
   int* firstedge;
   int* lastedge;

   assert(scip != NULL);
   assert(sourcedata != NULL);
   assert(targetdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, targetdata) );
   tcliqueCreate(&((*targetdata)->graph));              /* create the transformed graph */

   (*targetdata)->maxstablesets = sourcedata->maxstablesets;        /* copy length of array sets */
   (*targetdata)->nstablesets = sourcedata->nstablesets;            /* copy number of sets saved in array sets */
   (*targetdata)->oldgraph = sourcedata->oldgraph;      /* copy link to original graph */

   /* allocate memory for sets and lenghts of the sets */
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->stablesets), sourcedata->maxstablesets) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->stablesetlengths), sourcedata->maxstablesets) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->stablesetvars), sourcedata->maxstablesets) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->deletednodes), tcliqueGetNNodes(sourcedata->oldgraph)) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &((*targetdata)->new2oldnode), tcliqueGetNNodes(sourcedata->oldgraph)) );  

   for ( i = 0; i < tcliqueGetNNodes(sourcedata->oldgraph); i++ )
   {
      (*targetdata)->deletednodes[i] = sourcedata->deletednodes[i];
      (*targetdata)->new2oldnode[i] = sourcedata->new2oldnode[i];
   }

   /* copy stablesetlengths and stablesets */
   for ( i = 0; i < sourcedata->nstablesets; i++ )
   {
      assert(sourcedata->stablesetvars[i] != NULL);
      (*targetdata)->stablesetlengths[i] = sourcedata->stablesetlengths[i];
      SCIP_CALL( SCIPtransformVar(scip, sourcedata->stablesetvars[i], &((*targetdata)->stablesetvars[i])) );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*targetdata)->stablesets[i]),  sourcedata->stablesetlengths[i]) );
      for ( j = 0; j <sourcedata->stablesetlengths[i]; j++ )
      {
         (*targetdata)->stablesets[i][j] =  sourcedata->stablesets[i][j];
      }
   }

   /* create array for constraints */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(*targetdata)->constraints, tcliqueGetNNodes(sourcedata->graph)) );
   /* tranform constraints */
   SCIP_CALL( SCIPtransformConss(scip, tcliqueGetNNodes(sourcedata->graph), sourcedata->constraints,
         (*targetdata)->constraints) );
   /* copy the graph */
   tcliqueAddNode((*targetdata)->graph, tcliqueGetNNodes(sourcedata->graph)-1, 0);
   for ( i = 0; i < tcliqueGetNNodes(sourcedata->graph); i++ )
   {
      /* get adjacent nodes for node i */
      firstedge = tcliqueGetFirstAdjedge(sourcedata->graph, i);
      lastedge = tcliqueGetLastAdjedge(sourcedata->graph, i);
      while ( firstedge <= lastedge )
      {
         if ( *firstedge > i )
         {
            tcliqueAddEdge((*targetdata)->graph, i, *firstedge);
         }
         firstedge++;
      }
   }
   tcliqueFlush((*targetdata)->graph);

   return SCIP_OKAY;
}


/** deletes the transformed problem */
static
SCIP_DECL_PROBDELTRANS(probdeltransColoring)
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* relese constraints and free array for constraints */
   for ( i = 0; i < tcliqueGetNNodes((*probdata)->graph); i++)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->constraints[i])) );
   }
   SCIPfreeMemoryArray(scip, &((*probdata)->constraints));

   /* free the arrays for the stable sets and relese the related variables */
   for ( i = (*probdata)->nstablesets-1; i >= 0; i-- )
   {
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->stablesets[i]), (*probdata)->stablesetlengths[i]);
      SCIP_CALL( SCIPreleaseVar(scip, &((*probdata)->stablesetvars[i])) );
   }

   SCIPfreeMemoryArray(scip, &((*probdata)->new2oldnode));
   SCIPfreeMemoryArray(scip, &((*probdata)->deletednodes));
   SCIPfreeMemoryArray(scip, &((*probdata)->stablesetvars));
   SCIPfreeMemoryArray(scip, &((*probdata)->stablesetlengths));
   SCIPfreeMemoryArray(scip, &((*probdata)->stablesets));

   tcliqueFree(&(*probdata)->graph);
   SCIPfreeMemory(scip, probdata);
   return SCIP_OKAY;
}

static
SCIP_DECL_PROBDELORIG(probdelorigColoring)
{
   int i;

   assert(probdata != NULL);
   assert(*probdata != NULL);
  
   SCIPfreeMemoryArray(scip, &((*probdata)->new2oldnode));
   SCIPfreeMemoryArray(scip, &((*probdata)->deletednodes));  

   for ( i = (*probdata)->nstablesets-1; i >= 0; i-- )
   {
      SCIPfreeBlockMemoryArray(scip, &((*probdata)->stablesets[i]), (*probdata)->stablesetlengths[i]);
      SCIP_CALL( SCIPreleaseVar(scip, &((*probdata)->stablesetvars[i])) );
   }
   SCIPfreeMemoryArray(scip, &((*probdata)->stablesetvars));
   SCIPfreeMemoryArray(scip, &((*probdata)->stablesetlengths));
   SCIPfreeMemoryArray(scip, &((*probdata)->stablesets));

   /* release Constraints */
   for ( i = 0; i < tcliqueGetNNodes((*probdata)->graph); i++ )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &((*probdata)->constraints[i])) );
   }
   SCIPfreeMemoryArray(scip, &((*probdata)->constraints));   

   /* free memory used for graph */
   tcliqueFree(&((*probdata)->graph));
   tcliqueFree(&((*probdata)->oldgraph));

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}


/*
 * Callback methods of event handler
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecProbdatavardeleted)
{
   SCIP_VAR* var;
   SCIP_PROBDATA* probdata;
   int idx;

   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARDELETED);
   var = SCIPeventGetVar(event);
   probdata = (SCIP_PROBDATA*) eventdata;

   assert(probdata != NULL);
   assert(var != NULL);

   /* get index of variable in stablesets array */
   idx = (int)(size_t) SCIPvarGetData(var);

   SCIPdebugMessage("remove variable %s [%d] from list of stable sets\n", SCIPvarGetName(var), idx);

   assert(probdata->stablesetvars[idx] == var);

   /* remove variable from stablesets array and release it */
   SCIPfreeBlockMemoryArray(scip, &(probdata->stablesets[idx]), probdata->stablesetlengths[idx]);
   SCIP_CALL( SCIPreleaseVar(scip, &(probdata->stablesetvars[idx])) );

   /* move all subsequent variables to the front */
   for( ; idx < probdata->nstablesets - 1; idx++)
   {
      probdata->stablesets[idx] = probdata->stablesets[idx + 1];
      probdata->stablesetlengths[idx] = probdata->stablesetlengths[idx + 1];
      probdata->stablesetvars[idx] = probdata->stablesetvars[idx + 1];
      SCIPvarSetData(probdata->stablesetvars[idx], (SCIP_VARDATA*) (size_t) idx);
   }

   probdata->nstablesets--;

   return SCIP_OKAY;
}



/*
 * probdata specific interface methods
 */

/** sets up the problem data */
SCIP_RETCODE SCIPcreateProbColoring(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name,               /**< problem name */           
   int                   nnodes,             /**< number of nodes */
   int                   nedges,             /**< number of edges */
   int**                 edges               /**< array with start- and endpoints of the edges */
   )
{
   int i;
   SCIP_PROBDATA* probdata = NULL;
#ifndef NDEBUG
   TCLIQUE_Bool graphcreated;
#endif

   assert(nnodes > 0);  /* at least one node */
   assert(nedges >= 0); /* no negative number of edges */

   printf("Creating problem: %s \n", name);
   
   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, &probdata) );

   /* create graph */
#ifndef NDEBUG
   graphcreated = tcliqueCreate(&((probdata)->oldgraph));
   assert(graphcreated);
#else
   tcliqueCreate(&((probdata)->oldgraph));
#endif
   
   /* add all nodes from 0 to nnodes-1 */
   tcliqueAddNode((probdata)->oldgraph, nnodes-1, 0);

   /* add all edges, first into cache, then flush to add all of them to the graph */
   for ( i = 0; i < nedges; i++ )
   {
      assert((edges[i][0] > 0) && (edges[i][0] <= nnodes));
      assert((edges[i][1] > 0) && (edges[i][1] <= nnodes));

      tcliqueAddEdge((probdata)->oldgraph, edges[i][0]-1, edges[i][1]-1);
   }
   tcliqueFlush((probdata)->oldgraph);

   /* create constraints */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->constraints), nnodes) );

   /* at the beginning memory for 2 sets */
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->stablesets), 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->stablesetlengths), 2) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(probdata->stablesetvars), 2) );

   probdata->maxstablesets = 2;
   probdata->nstablesets = 0;

   /* include variable deleted event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC,
         eventExecProbdatavardeleted, NULL) );

   /* create problem in SCIP */
   SCIP_CALL( SCIPcreateProb(scip, name, probdelorigColoring, probtransColoring, probdeltransColoring, 
         NULL, NULL, NULL, probdata) );

   SCIP_CALL( preprocessGraph(scip) );

   return SCIP_OKAY;
}


/* ----------------------------------- external methods -------------------------- */

/** returns the number of stable sets / variables */
int COLORprobGetNStableSets(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->nstablesets;
}


/** prints all stable sets to standart output */
void COLORprobPrintStableSets(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;
   int i;
   int j;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   for ( i = 0; i < probdata->nstablesets; i++ )
   {
      printf( "Set %d: ", i);
      for ( j = 0; j < probdata->stablesetlengths[i]; j++ )
      {
         printf("%d, ", probdata->stablesets[i][j]+1);
      }
      printf("ub = %f", SCIPvarGetUbLocal(probdata->stablesetvars[i]));
      printf(", inLP = %d", SCIPvarIsInLP(probdata->stablesetvars[i]));
      //assert( (SCIPvarGetUbLocal(probdata->stablesetvars[i]) == 0) ||  SCIPvarIsInLP(probdata->stablesetvars[i]) );
      printf("\n");
   }
}


/** prints the requested stable set to standart output */
void COLORprobPrintStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setnumber           /**< the number of the requested set */
   )
{
   SCIP_PROBDATA* probdata;
   int i;
   int j;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   i = setnumber;
   printf( "Set %d: ", i);
   for ( j = 0; j < probdata->stablesetlengths[i]; j++ )
   {
      printf("%d, ", probdata->stablesets[i][j]+1);
   }
   if ( probdata->stablesetvars[i] != NULL )
      printf("ub = %f", SCIPvarGetUbLocal(probdata->stablesetvars[i]));
   printf("\n");
}




/** adds a variable that belongs to a given stable set */
SCIP_RETCODE COLORprobAddVarForStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,           /**< index of the stable set */
   SCIP_VAR*             var                 /**< pointer to the variable */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert((setindex >= 0) && (setindex < probdata->nstablesets));

   /* catch variable deleted event on the variable to update the stablesetvars array in the problem data */
   SCIP_CALL( SCIPcatchVarEvent(scip, var, SCIP_EVENTTYPE_VARDELETED, SCIPfindEventhdlr(scip, EVENTHDLR_NAME),
         (SCIP_EVENTDATA*) probdata, NULL) );

   probdata->stablesetvars[setindex] = var;

   return SCIP_OKAY;
}


/** gets the variable belonging to a given stable set */
SCIP_VAR* COLORprobGetVarForStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex            /**< index of the stable set */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert ( (setindex >= 0) && (setindex < probdata->nstablesets));

   return probdata->stablesetvars[setindex];
}


/** checks whether a node is in a given stable set, returns true iff it is */
SCIP_Bool COLORprobIsNodeInStableSet( 
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,           /**< index of the stable set */
   int                   node                /**< number of the node */
   )
{
   SCIP_PROBDATA* probdata;
   int l;
   int u;
   int m;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   l = 0;
   u = probdata->stablesetlengths[setindex]-1;
   while ( l <= u )
   {
      m = (l+u)/2;
      if ( probdata->stablesets[setindex][m] == node )
      {
         return TRUE;
      }
      if ( probdata->stablesets[setindex][m] > node )
      {
         l = m+1;
      }
      if ( probdata->stablesets[setindex][m] < node )
      {
         u = m-1;
      }
   }
   return FALSE;
}


/** checks whether the first set is equal to the second set, both sets have to be sorted in a decreasing way */
SCIP_Bool COLORprobStableSetsAreEqual(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  set1,               /**< array of nodes in the first set */ 
   int                   nset1nodes,         /**< number of nodes in the first set */
   int*                  set2,               /**< array of nodes in the second set */ 
   int                   nset2nodes          /**< number of nodes in the second set */
   )
{
   
   int i;

   assert(scip != NULL);
   assert(set1 != NULL && set2 != NULL);
   assert(nset1nodes > 0 && nset2nodes > 0);

   if ( nset1nodes != nset2nodes )
   {
      return FALSE;
   }
   for ( i = 0; i < nset1nodes; i++ )
   {
      if ( set1[i] != set2[i] )
      {
         return FALSE;
      }
   }
   return TRUE;

}


/** checks whether the given stable set is new
    returns TRUE if the stable is new, 
            FALSE if it is equal to an already existing stable set */
SCIP_Bool COLORprobStableSetIsNew(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  stablesetnodes,     /**< array of nodes in the stable set */
   int                   nstablesetnodes     /**< number of nodes in the stable set */
   )
{
   SCIP_PROBDATA* probdata; 
   int i;
   
   assert(stablesetnodes != NULL);
   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* sort the set */
   SCIPsortDownInt(stablesetnodes, nstablesetnodes);

   for ( i = 0; i < COLORprobGetNStableSets(scip); i++ )
   {
      if ( COLORprobStableSetsAreEqual(scip, stablesetnodes, nstablesetnodes, 
               probdata->stablesets[i], 
               probdata->stablesetlengths[i]) )
      {
         return FALSE;
      }
   }

   return TRUE;
}

/** adds a new stable set, the set must be sorted descendingly, 
 *  attention: you need to check whether it is new before adding it
 */
SCIP_RETCODE COLORprobAddNewStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int*                  stablesetnodes,     /**< array of nodes in the stable set */
   int                   nstablesetnodes,    /**< number of nodes in the stable set */
   int*                  setindex            /**< return value: index of the stable set */
   )
{
   SCIP_PROBDATA* probdata; 
   int newsize;
   int i;
   
   assert(stablesetnodes != NULL);
   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   /* the set should be sorted descendingly */
#ifndef NDEBUG
   for ( i = 0; i < nstablesetnodes-2; i++ )
   {
      assert(stablesetnodes[i]>stablesetnodes[i+1]);
   }
#endif
 
   /* ensure that array is big enough */
   if ( (probdata->nstablesets + 1) > probdata->maxstablesets)
   {
      newsize = 2* probdata->maxstablesets;
      assert(newsize >  probdata->nstablesets + 1);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->stablesets), newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->stablesetlengths), newsize) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &(probdata->stablesetvars), newsize) );
      probdata->maxstablesets = newsize;
      SCIPdebugMessage("Set-array resized: %d --> %d\n", newsize/2, newsize);
   }

   /* alloc memory for the new stable set */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->stablesets[probdata->nstablesets]), nstablesetnodes) );
   probdata->stablesetlengths[probdata->nstablesets] = nstablesetnodes;
   probdata->stablesetvars[probdata->nstablesets] = NULL;
   for ( i = 0; i < nstablesetnodes; i++ )
   {
      assert(stablesetnodes[i] >= 0);
      probdata->stablesets[probdata->nstablesets][i] = stablesetnodes[i];
   }
   *setindex = probdata->nstablesets;

   probdata->nstablesets++;

   return SCIP_OKAY;
}


/** returns the stable set with the given index */
void COLORprobGetStableSet(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   setindex,            /**< index of the stable set */
   int**                 stableset,          /**< return value: pointer to the stable set */
   int*                  nelements           /**< return value: number of elements in the stable set */
   )
{
   SCIP_PROBDATA* probdata; 

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   *stableset = probdata->stablesets[setindex];
   *nelements = probdata->stablesetlengths[setindex];
}


/** returns all stable sets */
void COLORprobGetStableSets(
   SCIP*                 scip,               /**< SCIP data structure */
   int***                stablesets,         /**< return value: pointer to the stable sets */
   int**                 nelements,          /**< return value: number of elements in the stable sets */
   int*                  nstablesets         /**< return value: number of sets */
   )
{
   SCIP_PROBDATA* probdata; 

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   *stablesets = probdata->stablesets;
   *nelements = probdata->stablesetlengths;
   *nstablesets = probdata->nstablesets;
}


/** returns the number of nodes in the (preprocessed) graph */
int COLORprobGetNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return tcliqueGetNNodes(probdata->graph);
}


/** returns the number of nodes in the original graph */
int COLORprobGetOriginalNNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return tcliqueGetNNodes(probdata->oldgraph);
}


/** returns the (preprocessed) graph */
TCLIQUE_GRAPH* COLORprobGetGraph(
   SCIP*                 scip                /**< SCIP data structure */
   )
{

   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->graph;
}


/** returns the original graph */
TCLIQUE_GRAPH* COLORprobGetOriginalGraph(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->oldgraph;
}


/** returns the array of nodes deleted during preprocessing, length = COLORprobGetOriginalNNodes(), filled with -1 at the end */
int* COLORprobGetDeletedNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->deletednodes;
}


/** returns the array in which for every node in the preprocessed graph, the related node in the original graph is saved */
int* COLORprobGetOriginalNodesForNewNodes(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->new2oldnode;
}



/** returns the node in the preprocessed graph, that belongs to the given node, returns -1 if node was deleted */
int COLORprobGetNewNodeForOriginalNode(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node                /**< a node in the original graph */
   )
{
   SCIP_PROBDATA* probdata;
   int i;
   
   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);

   assert(probdata != NULL);
   assert(node >= 0 && node < COLORprobGetOriginalNNodes(scip));
   
   for ( i = 0; i < COLORprobGetOriginalNNodes(scip); i++ )
   {
      if ( probdata->new2oldnode[i] == node )
         return i;
      if ( probdata->new2oldnode[i] == -1 )
         return -1;
   }
   return -1;
   
}


/** returns all node-constraints */
SCIP_CONS** COLORprobGetConstraints(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   return probdata->constraints;
}


/** returns the node-constraint belonging to a given node */
SCIP_CONS* COLORprobGetConstraint(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   node                /**< number of the node, for which this constraint assures coloring */
   )
{
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);
   assert(node >= 0 && node < tcliqueGetNNodes(probdata->graph));

   return probdata->constraints[node];
}


/** computes the complementary graph for a given graph and stores it in the given pointer */
void COLORprobGetComplementaryGraph(
   SCIP*                 scip,                /**< SCIP data structure */
   TCLIQUE_GRAPH*        graph,               /**< the given graph */
   TCLIQUE_GRAPH*        cgraph               /**< the complementary graph is saved in here */
   )
{   
   int nnodes;
   int i;
   int j;
   int*  firstedge;
   int*  lastedge;

   assert(scip != NULL);
   assert(graph != NULL);
   assert(cgraph != NULL);

   /* get number of nodes */
   nnodes = tcliqueGetNNodes(graph);
   assert(nnodes > 0);

   /* add all nodes from 0 to nnodes-1 */
   tcliqueAddNode(cgraph, nnodes-1, 0);
   /* add edge between i and j iff there is no edge between i and j in old graph */
   /* assumption: all edges are undirected, (i,j) exists --> (j,i) exists */
   for ( i = 0; i < nnodes; i++ )
   {
      firstedge = tcliqueGetFirstAdjedge(graph, i);
      lastedge = tcliqueGetLastAdjedge(graph, i);
      for ( j = 0; j < *firstedge && j < i; j++ )
      {
         tcliqueAddEdge(cgraph, i, j);
      }
      while ( firstedge < lastedge )
      {
         for ( j = *firstedge+1; j < *(firstedge+1) && j < i; j++ )
         {
            tcliqueAddEdge(cgraph, i, j);
         }
         firstedge++;
      }
      for ( j = (*lastedge)+1; j < COLORprobGetNNodes(scip) && j < i; j++ )
      {
         tcliqueAddEdge(cgraph, i, j);
      }
   }
   tcliqueFlush(cgraph);
   for ( i = 0; i < COLORprobGetNNodes(scip); i++ )
   {
      for ( j = i+1; j < COLORprobGetNNodes(scip); j++ )
      {
         assert((tcliqueIsEdge(graph, i, j) && !tcliqueIsEdge(cgraph, i, j))
            || (!tcliqueIsEdge(graph, i, j) && tcliqueIsEdge(cgraph, i, j)));
      }
   }
}


/** creates all node-constraints and saves them in an array */
SCIP_RETCODE COLORprobSetUpArrayOfCons(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONS** constraints;
   int nnodes;
   int i; 

   assert(scip != NULL);

   constraints = COLORprobGetConstraints(scip);
   assert(constraints != NULL);
   nnodes = COLORprobGetNNodes(scip);
   for ( i = 0; i < nnodes; i++ )
   {
      char consname[SCIP_MAXSTRLEN];
     
      /* create the constraint */
      sprintf(consname, "Node-Constraint%d", i+1);

      SCIP_CALL( SCIPcreateConsSetcover(scip, &constraints[i], consname, 0, NULL, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE) );
      SCIP_CALL( SCIPaddCons(scip, constraints[i]) );
   }

   return SCIP_OKAY;
}


/** checks whether the given node is in the given array */
SCIP_Bool COLORprobIsNodeInArray(
   int                   node,               /**< the number of the node */
   int*                  arraynodes,         /**< the nodes of the maximum stableset */
   int                   narraynodes         /**< number of nodes in the maximum stableset */
   )
{
   int i;

   assert(arraynodes != NULL);

   for ( i = 0; i < narraynodes; i++ )
   {
      if ( arraynodes[i] == node )
      {
         return TRUE;
      }
   }
   return FALSE;
}

/** checks whether the given two given arrays are equal */
SCIP_Bool COLORprobEqualSortedArrays(
   int*                  array1nodes,         /**< the nodes of the first set */
   int                   narray1nodes,        /**< number of nodes in the first set */
   int*                  array2nodes,         /**< the nodes of the second set */
   int                   narray2nodes         /**< number of nodes in the second set */
   )
{
   int i;

   assert(array1nodes != NULL);
   assert(narray1nodes > 0);
   assert(array2nodes != NULL);
   assert(narray2nodes > 0);

   if ( narray1nodes != narray2nodes )
   {
      return FALSE;
   }
   for ( i = 0; i < narray1nodes; i++ )
   {
      if ( array1nodes[i] != array2nodes[i] )
      {
         return FALSE;
      }
   }
   return TRUE;
}
