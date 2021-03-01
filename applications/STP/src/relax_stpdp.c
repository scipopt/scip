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

#define RELAX_NAME             "stpdp"
#define RELAX_DESC             "DP relaxator for STP"
#define RELAX_PRIORITY         100000000
#define RELAX_FREQ             1

/*
 * Data structures
 */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP_Bool             isActive;           /**< is the relaxator being used? */
};


/*
 * Local methods
 */


/** solves problem with terminals-FPT DP */
static
SCIP_RETCODE solveWithDpTerms(
   SCIP*                 scip,               /**< SCIP data structure */
   GRAPH*                graph,              /**< the graph */
   SCIP_Real*            obj
   )
{
   int* soledges;
   SCIP_Bool success;

   SCIP_CALL( SCIPallocMemoryArray(scip, &soledges, graph->edges) );

   printf("solving problem with DP... \n");

   SCIP_CALL( dpterms_solve(scip, graph, soledges) );
   SCIP_CALL( solstp_addSolToProb(scip, graph, soledges, NULL, &success) );
   assert(success);

   printf("...finished \n");

   *obj = solstp_getObj(graph, soledges, 0.0);

   SCIPfreeMemoryArray(scip, &soledges);

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

   *lowerbound = -SCIPinfinity(scip);
   *result = SCIP_DIDNOTRUN;

   if( !relaxdata->isActive )
      return SCIP_OKAY;

   graph = SCIPprobdataGetGraph2(scip);
   assert(graph);

   SCIP_CALL( solveWithDpTerms(scip, graph, lowerbound) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/*
 * Interface methods
 */


/** is using the relaxator promising? */
SCIP_Bool SCIPStpDpRelaxIsPromising(
   SCIP*                 scip,               /**< SCIP data structure */
   const GRAPH*          graph               /**< graph */
   )
{
   /*
   SCIP_RELAX* relax = SCIPfindRelax(scip, "stpdp");
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);

   assert(relaxdata);
*/
   assert(graph);
   assert(!dpterms_isPromising(graph)); // todo!


   return dpterms_isPromising(graph);
}


/** activates */
SCIP_RETCODE SCIPStpDpRelaxActivate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAX* relax = SCIPfindRelax(scip, "stpdp");
   SCIP_RELAXDATA* relaxdata = SCIPrelaxGetData(relax);

   assert(relaxdata);

   // todo want to check which DP algo to use!

   // todo maybe also turn SCIP presolving off
   SCIP_CALL( SCIPsetIntParam(scip, "heuristics/TM/initruns", 0) );

   relaxdata->isActive = TRUE;

   return SCIP_OKAY;
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

   return SCIP_OKAY;
}
