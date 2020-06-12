/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   examples/TSP/src/cppmain.cpp
 * @brief  main file for C++ TSP example using SCIP as a callable library
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 * This is an example of using SCIP to solve the TSP problem on undirected graphs. See the doxygen documentation for an
 * explanation.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <iostream>

/* include SCIP components */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* include TSP specific components */
#include "ReaderTSP.h"
#include "ConshdlrSubtour.h"
#include "HeurFarthestInsert.h"
#include "Heur2opt.h"
#include "HeurFrats.h"
#include "EventhdlrNewSol.h"

using namespace scip;
using namespace tsp;
using namespace std;

/** creates and runs a SCIP instance with default and TSP plugins */
static
SCIP_RETCODE runSCIP(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP* scip = NULL;


   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* we explicitly enable the use of a debug solution for this main SCIP instance */
   SCIPenableDebugSol(scip);

   /* include TSP specific plugins */
   SCIP_CALL( SCIPincludeObjReader(scip, new ReaderTSP(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjConshdlr(scip, new ConshdlrSubtour(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjEventhdlr(scip, new EventhdlrNewSol(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjHeur(scip, new HeurFarthestInsert(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjHeur(scip, new Heur2opt(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjHeur(scip, new HeurFrats(scip), TRUE) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, "sciptsp.set") );


   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method starting TSP code */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
