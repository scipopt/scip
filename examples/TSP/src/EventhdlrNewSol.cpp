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
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   EventhdlrNewSol.cpp
 * @brief  event handler for new solutions in TSP
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>
#include <fstream>
#include <iostream>
#ifndef _MSC_VER
#include <unistd.h>
#else
#include <windows.h>
#define sleep Sleep
#endif

#include "objscip/objscip.h"
#include "EventhdlrNewSol.h"
#include "ProbDataTSP.h"
#include "GomoryHuTree.h"

using namespace tsp;
using namespace std;


/** destructor of event handler to free user data (called when SCIP is exiting) */
SCIP_DECL_EVENTFREE(EventhdlrNewSol::scip_free)
{
   return SCIP_OKAY;
}


/** initialization method of event handler (called after problem was transformed) */
SCIP_DECL_EVENTINIT(EventhdlrNewSol::scip_init)
{

   int lockwaits = 0;
   while( SCIPfileExists("temp.tour.lock") && lockwaits < 10 )
   {
      /* wait one second and try again */
      sleep(1);
      lockwaits++;
   }
   if( SCIPfileExists("temp.tour.lock") )
   {
      SCIPwarningMessage(scip, "cannot reset, because lockfile <temp.tour.lock> is still existing\n");
      return SCIP_OKAY;
   }

   /* create lock file */
   ofstream lockfile("temp.tour.lock");
   lockfile << "lock" << endl;
   lockfile.close();

   // create output file which can be read by TSPViewer
   ofstream filedata("temp.tour");
   filedata << "RESET" << endl;
   filedata.close();

   /* delete lock file */
   unlink("temp.tour.lock");
   sleep(1); /* wait for the Java TSPViewer */

   return SCIP_OKAY;
}


/** deinitialization method of event handler (called before transformed problem is freed) */
SCIP_DECL_EVENTEXIT(EventhdlrNewSol::scip_exit)
{
   return SCIP_OKAY;
}


/** solving process initialization method of event handler (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The event handler may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_DECL_EVENTINITSOL(EventhdlrNewSol::scip_initsol)
{
   SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );
   return SCIP_OKAY;
}


/** solving process deinitialization method of event handler (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The event handler should use this call to clean up its branch and bound data.
 */
SCIP_DECL_EVENTEXITSOL(EventhdlrNewSol::scip_exitsol)
{
   SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1) );
   return SCIP_OKAY;
}


/** frees specific constraint data */
SCIP_DECL_EVENTDELETE(EventhdlrNewSol::scip_delete)
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
SCIP_DECL_EVENTEXEC(EventhdlrNewSol::scip_exec)
{
   SCIP_SOL* sol;
   sol = SCIPgetBestSol(scip);
   ProbDataTSP* probdata = dynamic_cast<ProbDataTSP*>(SCIPgetObjProbData(scip));
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
      SCIPwarningMessage(scip, "cannot output tour in file, because lockfile <temp.tour.lock> is still existing\n");
      return SCIP_OKAY;
   }

   /* create lock file */
   ofstream lockfile("temp.tour.lock");
   lockfile << "lock" << endl;
   lockfile.close();

   // create output file which can be read by TSPViewer
   ofstream filedata("temp.tour");
   filedata << graph->nnodes << endl;

   SCIP_HEUR* heur = SCIPgetSolHeur(scip, sol);
   if ( heur == NULL)
      filedata << "relaxation" << endl;
   else
      filedata << SCIPheurGetName(heur) << endl;

   filedata << SCIPgetSolOrigObj(scip,sol) << endl;
   do
   {
      // output the number of nodes
      filedata << node->id << " " << node->x << " " << node->y << endl;
      GRAPHEDGE* edge = node->first_edge;

      // output the id and the coordinates of the nodes in the order they appear in the tour

      while( edge != NULL)
      {
         if( edge->back != lastedge && SCIPgetSolVal(scip, sol, edge->var) > 0.5 )
         {
            node = edge->adjac;
            lastedge = edge;
            break;
         }
         edge = edge->next;
      }
   }
   while ( node != &graph->nodes[0] );

   filedata.close();

   /* delete lock file */
   unlink("temp.tour.lock");
   sleep(1); /* wait for the Java TSPViewer */

   return SCIP_OKAY;
}
