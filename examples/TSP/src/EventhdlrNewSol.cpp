/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: EventhdlrNewSol.cpp,v 1.1 2005/03/03 16:43:34 bzfberth Exp $"

/**@file   EventhdlrNewSol.cpp
 * @brief  event handler for new solutions in TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include "EventhdlrNewSol.h"
#include "TSPProbData.h"
#include "gminucut.h"


using namespace tsp;
using namespace std;

     
/** destructor of event handler to free user data (called when SCIP is exiting) */
RETCODE EventhdlrNewSol::scip_free(
   SCIP*         scip,               /**< SCIP data structure */
   EVENTHDLR*    eventhdlr           /**< the event handler itself */
   )
{
   return SCIP_OKAY;
}
   
/** initialization method of event handler (called after problem was transformed) */
RETCODE EventhdlrNewSol::scip_init(
   SCIP*         scip,               /**< SCIP data structure */
   EVENTHDLR*    eventhdlr           /**< the event handler itself */
   )
{
   return SCIP_OKAY;
}
   
/** deinitialization method of event handler (called before transformed problem is freed) */
RETCODE EventhdlrNewSol::scip_exit(
   SCIP*         scip,               /**< SCIP data structure */
   EVENTHDLR*    eventhdlr           /**< the event handler itself */
   )
{
   return SCIP_OKAY;
}
   
/** solving process initialization method of event handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The event handler may use this call to initialize its branch and bound specific data.
 *
 */
RETCODE EventhdlrNewSol::scip_initsol(
   SCIP*         scip,               /**< SCIP data structure */
   EVENTHDLR*    eventhdlr           /**< the event handler itself */
   )
{
   CHECK_OKAY( SCIPcatchEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );
   return SCIP_OKAY;
}
   
/** solving process deinitialization method of event handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The event handler should use this call to clean up its branch and bound data.
 */
RETCODE EventhdlrNewSol::scip_exitsol(
   SCIP*         scip,               /**< SCIP data structure */
   EVENTHDLR*    eventhdlr           /**< the event handler itself */
   )
{
   CHECK_OKAY( SCIPdropEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1) );
   return SCIP_OKAY;
}
   
/** frees specific constraint data */
RETCODE EventhdlrNewSol::scip_delete(
   SCIP*         scip,               /**< SCIP data structure */
   EVENTHDLR*    eventhdlr,          /**< the event handler itself */
   EVENTDATA**   eventdata           /**< pointer to the event data to free */
   )
{
   return SCIP_OKAY;
}

/** execution method of event handler
 *
 *  Processes the event. The method is called every time an event occurs, for which the event handler
 *  is responsible. Event handlers may declare themselves resposible for events by calling the
 *  corresponding SCIPcatch...() method. This method creates an event filter object to point to the
 *  given event handler and event data.
 */
RETCODE EventhdlrNewSol::scip_exec(
   SCIP*         scip,               /**< SCIP data structure */
   EVENTHDLR*    eventhdlr,          /**< the event handler itself */
   EVENT*        event,              /**< event to process */
   EVENTDATA*    eventdata           /**< user data for the event */
   ) 
{
   SOL* sol;  
   sol = SCIPgetBestSol(scip);
   TSPProbData* probdata = dynamic_cast<TSPProbData*>(SCIPgetObjProbData(scip));
   GRAPH* graph = probdata->getGraph();
   GRAPHNODE* node = &graph->nodes[0];
   GRAPHEDGE* lastedge = NULL;
   int lockwaits = 0;

   /* wait for lock file to disappear */
   while( SCIPfileExists("temp.tour.lock") && lockwaits < 10 )
   {
      /* wait one second and try again */
      sleep(1);
      lockwaits++;
   }
   if( SCIPfileExists("temp.tour.lock") )
   {
      warningMessage("cannot output tour in file, because lockfile <temp.tour.lock> is still existing\n");
      return SCIP_OKAY;
   }

   /* create lock file */
   ofstream lockfile("temp.tour.lock");
   lockfile << "lock" << endl;
   lockfile.close();

   // create output file which can be read by TSPViewer 
   ofstream filedata("temp.tour");
   filedata << graph->nnodes << endl;
   do
   {
      // output the number of nodes
      filedata << node->id << " " << node->x << " " << node->y << endl;
      GRAPHEDGE* edge = node->first_edge;
     
      // output the id and the coordinates of the nodes in the order they appear in the tour
      if( edge != NULL )
      {
         do
         {
            if( edge->back != lastedge && SCIPgetSolVal(scip, sol, edge->var) > 0.5 )
            {
               node = edge->adjac;
               lastedge = edge;
               break;
            }
            edge = edge->next;
         }
         while ( edge != node->first_edge );
         assert( edge != node->first_edge );
      }
   }
   while ( node != &graph->nodes[0] );

   filedata.close();

   /* delete lock file */
   unlink("temp.tour.lock");

   return SCIP_OKAY;
}

