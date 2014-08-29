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

/**@file   relax_reoptsolvelp.c
 * @brief  reoptsolvelp relaxator
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/relax_reoptsolvelp.h"
#include "scip/scip.h"
#include "scip/branch_nodereopt.h"
#include "scip/struct_tree.h"


#define RELAX_NAME             "reoptsolvelp"
#define RELAX_DESC             "relaxator to check if the LP should be solved"
#define RELAX_PRIORITY         1
#define RELAX_FREQ             1


/*
 * Data structures
 */

/* TODO: fill in the necessary relaxator data */

/** relaxator data */
struct SCIP_RelaxData
{
   SCIP_Bool             reopt;
   SCIP_Longint          lastnode;
};

/*
 * Local methods
 */

/*
 * Callback methods of relaxator
 */

/** copy method for relaxator plugins (called when SCIP copies plugins) */
#define relaxCopyreoptsolvelp NULL

/** destructor of relaxator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_RELAXFREE(relaxFreereoptsolvelp)
{
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL );

   SCIPfreeMemory(scip, &relaxdata);

   return SCIP_OKAY;
}

/** initialization method of relaxator (called after problem was transformed) */
#define relaxInitreoptsolvelp NULL

/** deinitialization method of relaxator (called before transformed problem is freed) */
#define relaxExitreoptsolvelp NULL

/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolreoptsolvelp)
{
   SCIP_RELAXDATA* relaxdata;

   assert(scip != NULL );
   assert(relax != NULL );

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL );

   relaxdata->lastnode = -1;

   SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/enable", &relaxdata->reopt) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
#define relaxExitsolreoptsolvelp NULL

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecreoptsolvelp)
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_NODE* node;
   SCIP_Bool solvelp;

   assert(scip != NULL);
   assert(relax != NULL);

   relaxdata = SCIPrelaxGetData(relax);
   assert(relaxdata != NULL);

   node = SCIPgetCurrentNode(scip);

   if( !relaxdata->reopt )
      return SCIP_OKAY;

   /** check is the node come from reoptimiziation */
   if( SCIPgetRootNode(scip) == node || (SCIPnodeGetReoptID(node) > -1 && relaxdata->lastnode != SCIPnodeGetNumber(node)) )
   {
      relaxdata->lastnode = SCIPnodeGetNumber(node);
      SCIP_CALL( SCIPbranchruleNodereoptSolveLP(scip, node, &solvelp) );
      SCIPsetFocusnodeLP(scip, solvelp);
   }

   return SCIP_OKAY;
}

/*
 * relaxator specific interface methods
 */

/** creates the reoptsolvelp relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxReoptsolvelp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   /* create reoptsolvelp relaxator data */
   SCIP_CALL( SCIPallocMemory(scip, &relaxdata) );

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecreoptsolvelp, relaxdata) );

   assert(relax != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetRelaxFree(scip, relax, relaxFreereoptsolvelp) );
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitsolreoptsolvelp) );

   return SCIP_OKAY;
}
