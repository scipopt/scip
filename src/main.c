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

/**@file   scip/src/main.c
 * @ingroup OTHER_CFILES
 * @brief  main file for C compilation
 * @author Tobias Achterberg
 */

/* global todos: */

/**@todo pricing for pseudo solutions */
/**@todo unboundness detection in presolving -> convert problem into feasibility problem to decide
 *       unboundness/infeasibility */
/**@todo variable event PSSOLCHANGED, update pseudo activities in constraints to speed up checking of pseudo solutions */
/**@todo branching rule acting as a filter by temporary changing the branching priority of variables and returning
 *       SCIP_DIDNOTFIND to let the next branching rule select the branching variable */
/**@todo try to not use the first but the shortest constraint as reason for a deduction */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <signal.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/interrupt.h"

/** callback function for handling signals */ /*lint -e715*/
static
void handleSigterm(
   int                   signum              /**< signal code */
   )
{ /*lint --e{715}*/
   /* Calling the following function is not directly async-signal-safe, since it counts how often the function is
    * called. For achieving the goal of terminating this seems unproblematic. */
   SCIPtryTerminate(); /*lint !e2761*/
}

/** main method starting SCIP */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode;

   (void)signal(SIGTERM, handleSigterm);
   /* run interactive shell */
   retcode = SCIPrunShell(argc, argv, "scip.set");

   /* evaluate retrun code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
