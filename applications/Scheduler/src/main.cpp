/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   Scheduler/src/main.cpp
 * @brief  Main file for C++ compilation
 * @author Stefan Heinz
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"

#include "cons_optcumulative.h"
#include "heur_optcumulative.h"
#include "heur_listscheduling.h"
#include "reader_cmin.h"
#include "reader_sch.h"
#include "reader_sm.h"
#include "reader_rcp.h"

/** runs the shell */
static
SCIP_RETCODE runShell(
   int                        argc,               /**< number of shell parameters */
   char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
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

   /* include default plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include problem reader */
   SCIP_CALL( SCIPincludeReaderCmin(scip) );
   SCIP_CALL( SCIPincludeReaderSch(scip) );
   SCIP_CALL( SCIPincludeReaderSm(scip) );
   SCIP_CALL( SCIPincludeReaderRcp(scip) );

   /* include problem specific heuristic */
   SCIP_CALL( SCIPincludeHeurListScheduling(scip) );
   SCIP_CALL( SCIPincludeHeurOptcumulative(scip) );

   /* include cumulative constraint handler with optional activities */
   SCIP_CALL( SCIPincludeConshdlrOptcumulative(scip) );

#ifdef WITH_CPOPTIMIZER
   SCIP_CALL( SCIPsetSolveCumulative(scip, cpoptimizer) );
#endif

   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );

   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   /* check block memory */
   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method */
int main(
   int                   argc,          /**< number of arguments */
   char**                argv           /**< string array with arguments */
   )
{
  SCIP_RETCODE retcode;

  retcode = runShell(argc, argv, "scip.set");

  if( retcode != SCIP_OKAY )
  {
     SCIPprintError(retcode);
     return -1;
  }

  return 0;
}
