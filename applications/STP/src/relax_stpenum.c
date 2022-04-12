/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_stpenum.c
 * @brief  Steiner tree enumeration relaxator
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>
#include "relax_stpenum.h"
#include "probdata_stp.h"
#include "prop_stp.h"
#include "mst.h"
#include "solstp.h"
#include "relax_stpdp.h"

#define RELAX_NAME             "stpenum"
#define RELAX_DESC             "enumeration relaxator for STP"
#define RELAX_PRIORITY         10
#define RELAX_FREQ             1

#define ENUM_MAXNSTEINVERTS    8

/*
 * Data structures
 */


/*
 * Local methods
 */


/** collects and returns all non-terminals */
static
STP_Vectype(int) getSteinerVertices(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph               /**< graph data structure */
)
{
   STP_Vectype(int) steinverts = NULL;
   const int nnodes = graph_get_nNodes(graph);

   for( int i = 0; i < nnodes; i++ )
   {
      if( !Is_term(graph->term[i]) && graph->grad[i] > 0 )
         StpVecPushBack(scip, steinverts, i);
   }

   assert(StpVecGetSize(steinverts) <= ENUM_MAXNSTEINVERTS);

   return steinverts;
}


/** are all terminals visited by MST? */
static
SCIP_Bool allTermsAreVisited(
   const GRAPH*          graph,              /**< graph data structure */
   const MST*            mst                 /**< MST */
)
{
   const int nnodes = graph_get_nNodes(graph);
   const int* const terms = graph->term;
   const int* const nodes_prededge = mst->nodes_predEdge;

   for( int i = 0; i < nnodes; i++ )
   {
      if( !Is_term(terms[i]) )
         continue;

      if( nodes_prededge[i] == UNKNOWN && i != graph->source )
      {
         SCIPdebugMessage("terminal vertex %d not visited, solution is invalid! \n", i);
         return FALSE;
      }
   }

   return TRUE;
}


/** do enumeration */
static
SCIP_RETCODE enumExec(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int* RESTRICT         edges_solstat,      /**< solution edges */
   SCIP_Real*            obj                 /**< optimal objective value */
)
{
   MST* mst;
   STP_Bool* nodes_isActive;
   STP_Vectype(int) steinverts = getSteinerVertices(scip, graph);
   const int nnodes = graph_get_nNodes(graph);
   const uint32_t nsteinverts = (uint32_t) StpVecGetSize(steinverts);
   const uint32_t powsize = (uint32_t) pow(2.0, nsteinverts);
   *obj = FARAWAY;

   SCIP_CALL( mst_init(scip, graph, &mst) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &nodes_isActive, graph->knots) );

#ifndef WITH_UG
   printf("solving problem with enumeration (%u Steiner vertices) \n", nsteinverts);
#endif

   for( int i = 0; i < nnodes; i++ )
      nodes_isActive[i] = (graph->grad[i] > 0);
   nodes_isActive[graph->source] = TRUE;

   /* loop over all subsets of Steiner vertices */
   for( uint32_t bitset = 0; bitset < powsize; bitset++ )
   {
      for( uint32_t i = 0; i < nsteinverts; i++ )
      {
         assert(nodes_isActive[steinverts[i]]);
         if( (bitset & (1 << i)) )
            nodes_isActive[steinverts[i]] = FALSE;
      }

      mst_computeOnMarked(graph, nodes_isActive, graph->source, mst);

      if( allTermsAreVisited(graph, mst) )
      {
         const SCIP_Real obj_curr = mst_getObj(graph, mst);

         if( LT(obj_curr, *obj) )
         {
            *obj = obj_curr;
            SCIPdebugMessage("updating MST objective to %f \n", obj_curr);
            mst_getSoledges(graph, mst, edges_solstat);
            assert(EQ(obj_curr, solstp_getObj(graph, edges_solstat, 0.0)));
            assert(solstp_isValid(scip, graph, edges_solstat));
         }
      }

      /* reset */
      for( uint32_t i = 0; i < nsteinverts; i++ )
         nodes_isActive[steinverts[i]] = TRUE;
   }

   StpVecFree(scip, steinverts);
   mst_free(scip, &mst);
   SCIPfreeMemoryArray(scip, &nodes_isActive);

   assert(LT(*obj, FARAWAY));
   assert(solstp_isValid(scip, graph, edges_solstat));

   return SCIP_OKAY;
}


/** is enumeration promising? */
static
SCIP_Bool enumIsPromising(
   const GRAPH*          graph               /**< graph data structure */
)
{
   int nsteinverts = 0;
   const int nnodes = graph_get_nNodes(graph);
   const int* const terms = graph->term;

   for( int i = 0; i < nnodes; i++ )
   {
      if( !Is_term(terms[i]) )
      {
         nsteinverts++;
         if( nsteinverts > ENUM_MAXNSTEINVERTS )
         {
            SCIPdebugMessage("enumeration of Steiner vertices is not promising \n");
            return FALSE;
         }
      }
   }

   return TRUE;
}



/*
 * Callback methods of relaxator
 */


/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeStpenum)
{  /*lint --e{715}*/
   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecStpenum)
{  /*lint --e{715}*/
   GRAPH* graph;
   int* edges_solstat;
   SCIP_Bool success;
   SCIP_Longint graphnodenumber;
   SCIP_Bool probisinfeas;
   SCIP_Real offset = 0.0;

   *lowerbound = -SCIPinfinity(scip);
   *result = SCIP_DIDNOTRUN;

   if( !graph_typeIsSpgLike(SCIPprobdataGetGraph2(scip)) )
      return SCIP_OKAY;

   if( SCIPStpDpRelaxIsActive(scip) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPStpPropGetGraph(scip, &graph, &graphnodenumber, &probisinfeas, &offset) );

   if( probisinfeas )
   {
      SCIPdebugMessage("STP enumeration relaxator found infeasibility \n");
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   assert(graphnodenumber == SCIPnodeGetNumber(SCIPgetCurrentNode(scip)));

   if( !enumIsPromising(graph) )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocMemoryArray(scip, &edges_solstat, graph->edges) );
   SCIP_CALL( enumExec(scip, graph, edges_solstat, lowerbound) );

   SCIP_CALL( solstp_addSolToProb(scip, graph, edges_solstat, NULL, &success) );
   SCIPfreeMemoryArray(scip, &edges_solstat);

   *lowerbound += offset;
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods
 */


/** is using the relaxator promising? */
SCIP_Bool SCIPStpEnumRelaxIsPromising(
   const GRAPH*          graph               /**< graph */
   )
{
   assert(graph);

   return enumIsPromising(graph);
}


/** Solve instance by enumeration. Only call when promising. */
SCIP_RETCODE SCIPStpEnumRelaxComputeSol(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph,              /**< graph data structure */
   int* RESTRICT         edges_solstat       /**< solution edges */
)
{
   SCIP_Real obj;

   if( !SCIPStpEnumRelaxIsPromising(graph) )
   {
      SCIPerrorMessage("cannot call enumeration algorithm with too many non-terminals! \n");
      return SCIP_ERROR;
   }

   SCIP_CALL( enumExec(scip, graph, edges_solstat, &obj) );
   SCIPdebugMessage("enumeration found solution with objective %f \n", obj);

   return SCIP_OKAY;
}


/** creates the stp relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxStpenum(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata = NULL;
   SCIP_RELAX* relax = NULL;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecStpenum, relaxdata) );
   assert(relax != NULL);

   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeStpenum) );

   return SCIP_OKAY;
}
