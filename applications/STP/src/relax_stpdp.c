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

/**@file   relax_stpdp.c
 * @brief  Steiner tree relaxator
 * @author Daniel Rehfeldt
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "relax_stpdp.h"
#include "probdata_stp.h"
#include "solstp.h"
#include "dpterms.h"
#include "dpborder.h"
#ifndef NDEBUG
#include "substpsolver.h"
#endif

#define RELAX_NAME             "stpdp"
#define RELAX_DESC             "DP relaxator for STP"
#define RELAX_PRIORITY         100000000
#define RELAX_FREQ             1

enum DP_TYPE {dp_border, dp_terms};


/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   DPBORDER*             dpborder;           /**< DP border */
   enum DP_TYPE          mode;               /**< DP algo to be used */
   SCIP_Bool             isActive;           /**< is the relaxator being used? */
};


/*
 * Local methods
 */


/** solves problem with border-FPT DP */
static
SCIP_RETCODE solveWithDpBorder(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   DPBORDER*             dpborder,           /**< DP border algorithm data structure */
   SCIP_Real*            obj,                /**< pointer to the objective value (OUT) */
   SCIP_Bool*            wasSolved           /**< pointer to mark whether problem was solved (OUT) */
   )
{
   int* soledges;
   SCIP_Bool success;

   SCIP_CALL( SCIPallocMemoryArray(scip, &soledges, graph->edges) );

   printf("solving problem with DPB... ");
   graph_printInfo(graph);

   SCIP_CALL( dpborder_solve(scip, graph, dpborder, soledges, wasSolved) );

   if( !(*wasSolved) )
   {
      printf("...aborted \n");

      SCIPfreeMemoryArray(scip, &soledges);

      return SCIP_OKAY;
   }

   SCIP_CALL( solstp_addSolToProb(scip, graph, soledges, NULL, &success) );

   printf("...finished \n");

   *obj = solstp_getObj(graph, soledges, 0.0);

   SCIPfreeMemoryArray(scip, &soledges);

   return SCIP_OKAY;
}


/** solves problem with terminals-FPT DP */
static
SCIP_RETCODE solveWithDpTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real*            obj,                /**< pointer to the objective value (OUT) */
   SCIP_Bool*            wasSolved           /**< pointer to mark whether problem was solved (OUT) */
   )
{
   int* soledges;
   SCIP_Bool success;

   SCIP_CALL( SCIPallocMemoryArray(scip, &soledges, graph->edges) );

   printf("solving problem with DP... \n");

   SCIP_CALL( dpterms_solve(scip, graph, soledges, wasSolved) );

   if( !(*wasSolved) )
   {
      printf("...aborted \n");

      SCIPfreeMemoryArray(scip, &soledges);

      return SCIP_OKAY;
   }

   SCIP_CALL( solstp_addSolToProb(scip, graph, soledges, NULL, &success) );

   printf("...finished \n");

   *obj = solstp_getObj(graph, soledges, 0.0);

   SCIPfreeMemoryArray(scip, &soledges);

#ifndef NDEBUG
   {
      SCIP_Real obj_bc;
      SCIP_CALL( substpsolver_getObjFromGraph(scip, graph, &obj_bc) );
      assert(EQ(obj_bc, *obj));
   }
#endif

   return SCIP_OKAY;
}


/*
 * Callback methods of relaxator
 */


/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreeStpdp)
{  /*lint --e{715}*/
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata);

   if( relaxdata->dpborder )
      dpborder_free(scip, &(relaxdata->dpborder));

   SCIPfreeMemory(scip, &relaxdata);
   SCIPrelaxSetData(relax, NULL);

   return SCIP_OKAY;
}


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolStpdp)
{  /*lint --e{715}*/

   //SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);
   //assert(relaxdata);

   return SCIP_OKAY;
}



/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolStpdp)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecStpdp)
{  /*lint --e{715}*/
   GRAPH* graph;
   SCIP_RELAXDATA* const relaxdata = SCIPrelaxGetData(relax);
   SCIP_Bool success;

   *lowerbound = -SCIPinfinity(scip);
   *result = SCIP_DIDNOTRUN;

   if( !relaxdata->isActive )
      return SCIP_OKAY;

   graph = SCIPprobdataGetGraph2(scip);
   assert(graph);

   if( relaxdata->mode == dp_terms )
   {
      SCIP_CALL( solveWithDpTerms(scip, graph, lowerbound, &success) );
   }
   else
   {
      assert(relaxdata->mode == dp_border);
      assert(relaxdata->dpborder);

      SCIP_CALL( solveWithDpBorder(scip, graph, relaxdata->dpborder, lowerbound, &success) );
   }

   if( !success )
   {
      assert(SCIPisStopped(scip));
      return SCIP_OKAY;
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** is using the relaxator promising? */
SCIP_Bool SCIPStpDpRelaxIsPromising(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph               /**< graph */
   )
{
   SCIP_Bool dpbHasPotential;
   SCIP_RELAX* relax = SCIPfindRelax(scip, "stpdp");
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);

   assert(graph);
   assert(relaxdata);
   assert(!relaxdata->dpborder);

   if( dpterms_isPromisingPartly(graph) )
   {
      relaxdata->mode = dp_terms;
      return TRUE;
   }

   SCIP_CALL( dpborder_init(scip, graph, &(relaxdata->dpborder)) );
   SCIP_CALL( dpborder_probePotential(scip, graph, relaxdata->dpborder, &dpbHasPotential) );

   if( dpbHasPotential )
   {
      relaxdata->mode = dp_border;
      return TRUE;
   }

   return FALSE;
}


/** activates */
SCIP_RETCODE SCIPStpDpRelaxActivate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax = SCIPfindRelax(scip, "stpdp");
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);

   assert(relaxdata);

   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/TM/initruns", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "presolving/maxrounds", 0) );
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/trivial/freq", -1) );
   SCIP_CALL( SCIPsetIntParam(scip, "propagating/maxroundsroot", 0) );

   relaxdata->isActive = TRUE;

   return SCIP_OKAY;
}


/** is active? */
SCIP_Bool SCIPStpDpRelaxIsActive(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax = SCIPfindRelax(scip, "stpdp");
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);

   assert(relaxdata);

   return relaxdata->isActive;
}


/** creates the relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxStpdp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecStpdp, relaxdata) );
   assert(relax != NULL);

   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitsolStpdp) );
   SCIP_CALL( SCIPsetRelaxExitsol(scip, relax, relaxExitsolStpdp) );
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreeStpdp) );

   relaxdata->isActive = FALSE;
   relaxdata->mode = dp_terms;
   relaxdata->dpborder = NULL;

   return SCIP_OKAY;
}
