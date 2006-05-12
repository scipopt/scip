/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cppmain.cpp,v 1.8 2006/05/12 08:56:38 bzfpfend Exp $"

/**@file   cppmain.cpp
 * @brief  main file for C++ TSP example using SCIP as a callable library
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <iostream>

/* include SCIP components */
#include "objscip/objscip.h"


extern "C"
{
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"
}

/* include TSP specific components */
#include "ReaderTSP.h"
#include "ConshdlrSubtour.h"
#include "HeurFarthestInsert.h"
#include "Heur2opt.h"
#include "EventhdlrNewSol.h"

using namespace scip;
using namespace tsp;
using namespace std;

static
SCIP_RETCODE runSCIP(
   int                   argc,
   char**                argv
   )
{
   SCIP* scip = NULL;


   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include TSP specific plugins */
   SCIP_CALL( SCIPincludeObjReader(scip, new ReaderTSP(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjConshdlr(scip, new ConshdlrSubtour(), TRUE) ); 
   SCIP_CALL( SCIPincludeObjEventhdlr(scip, new EventhdlrNewSol(), TRUE) );
   SCIP_CALL( SCIPincludeObjHeur(scip, new HeurFarthestInsert(), TRUE) );
   SCIP_CALL( SCIPincludeObjHeur(scip, new Heur2opt(), TRUE) );
   
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

int
main(
   int                   argc,
   char**                argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode, stderr);
      return -1;
   }

   return 0;
}
