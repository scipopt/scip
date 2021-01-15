/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_stp.c
 * @brief  Steiner tree relaxator
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "relax_stp.h"
#include "prop_stp.h"
#include "dualascent.h"

#define RELAX_NAME             "stp"
#define RELAX_DESC             "relaxator for STP"
#define RELAX_PRIORITY         0
#define RELAX_FREQ             1


/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP_Longint          lastnodenumber;     /**< node number of last call */
   SCIP_Bool             isActive;
};


/*
 * Local methods
 */




/** computes lower bound */
static
SCIP_RETCODE runDualAscent(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real*            lowerbound
   )
{
   const SCIP_Bool mw = (graph->stp_type == STP_MWCSP);
   const SCIP_Bool pc = (graph->stp_type == STP_PCSPG);
   const SCIP_Bool doAscendPrune = !FALSE; // todo?

   if( doAscendPrune )
   {
      SCIP_CALL( graph_path_init(scip, graph) );
   }

   // int dod

   if( pc || mw )
   {
      SCIP_CALL( dualascent_execPcMw(scip, graph, NULL, lowerbound, FALSE, doAscendPrune, 1) );
   }
   else
   {
      DAPARAMS daparams = { .addcuts = FALSE, .ascendandprune = doAscendPrune, .root = graph->source,
                   .is_pseudoroot = FALSE, .damaxdeviation = 0.1 };

      if( graph->stp_type == STP_DCSTP )
      {
         SCIP_CALL( dualascent_execDegCons(scip, graph, NULL, &daparams, NULL, lowerbound) );
      }
      else
      {
         SCIP_CALL( dualascent_exec(scip, graph, NULL, &daparams, NULL, lowerbound) );
      }
   }

   if( doAscendPrune )
      graph_path_exit(scip, graph);

   return SCIP_OKAY;
}

/*
 * Callback methods of relaxator
 */


/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeStp)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata);

   SCIPfreeMemory(scip, &relaxdata);

   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolStp)
{  /*lint --e{715}*/

   SCIP_Bool usedarelax;
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata);

   SCIP_CALL( SCIPgetBoolParam(scip, "stp/usedarelax", &usedarelax) );

   relaxdata->lastnodenumber = -1;
   relaxdata->isActive = usedarelax;

   if( usedarelax )
   {
      SCIP_CALL( SCIPsetIntParam(scip, "lp/solvefreq", -1) );
   }

   return SCIP_OKAY;
}



/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolStp)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}



/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecStp)
{  /*lint --e{715}*/
   GRAPH* graph;
   SCIP_RELAXDATA* const relaxdata = SCIPrelaxGetData(relax);
   SCIP_Longint nodenumber;
   SCIP_Longint graphnodenumber;
   SCIP_Bool probisinfeas;

   *lowerbound = -SCIPinfinity(scip);
   *result = SCIP_DIDNOTRUN;

   if( !relaxdata->isActive )
      return SCIP_OKAY;

   nodenumber = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
   if( nodenumber == relaxdata->lastnodenumber )
   {
      relaxdata->lastnodenumber = nodenumber;
      printf("same node number (%lld), don't execute relaxator \n", nodenumber);

      return SCIP_OKAY;
   }

   relaxdata->lastnodenumber = nodenumber;

   SCIP_CALL( SCIPStpPropGetGraph(scip, &graph, &graphnodenumber, &probisinfeas) );

   if( probisinfeas )
   {
      printf("STP relaxator found infeasibility \n");
      *result = SCIP_CUTOFF;
   }
   else
   {
      assert(graph);
      assert(graphnodenumber == nodenumber);
      // todo there might be a problem for pc/mw due to the offset from reomving nodes!
      assert(!graph_pc_isPcMw(graph) && "todo: adapt offset");

      SCIP_CALL( runDualAscent(scip, graph, lowerbound) );

      *lowerbound += SCIPprobdataGetOffset(scip);

      printf("Stp lower bound = %f \n", *lowerbound);
      *result = SCIP_SUCCESS;
   }

   // todo possibly compute solution by ascend-prune and add it

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods
 */

/** is the relaxator active? */
SCIP_Bool SCIPStpRelaxIsActive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax = SCIPfindRelax(scip, "stp");
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);

   assert(relaxdata);

   return relaxdata->isActive;
}


/** creates the stp relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecStp, relaxdata) );
   assert(relax != NULL);

   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitsolStp) );
   SCIP_CALL( SCIPsetRelaxExitsol(scip, relax, relaxExitsolStp) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeStp) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "stp/usedarelax",
         "use dual ascent bounds only?",
         NULL, FALSE, FALSE, NULL, NULL) );

   relaxdata->lastnodenumber = -1;
   relaxdata->isActive = FALSE;

   return SCIP_OKAY;
}
