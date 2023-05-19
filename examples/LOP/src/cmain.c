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

/**@file   examples/LOP/src/cmain.c
 * @brief  main file for linear ordering example
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include "cons_lop.h"
#include "reader_lop.h"

/** define macro to print error message and exit */
#define SCIP_CALL_ERROR(x)   do                                                                               \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) != SCIP_OKAY )                                                 \
                          {                                                                                   \
                             SCIPprintError(_restat_);                                                        \
                             return -1;                                                                       \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )


/** main function, which starts the solution of the linear ordering problem */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP* scip = NULL;

   /* initialize SCIP */
   SCIP_CALL_ERROR( SCIPcreate(&scip) );

   /* output version information */
   SCIPinfoMessage(scip, NULL, "Solving the linear ordering problem using SCIP.\n");
   SCIPinfoMessage(scip, NULL, "\n");

   /* include default SCIP plugins */
   SCIP_CALL_ERROR( SCIPincludeDefaultPlugins(scip) );

   /* include linear ordering constraint handler */
   SCIP_CALL_ERROR( SCIPincludeConshdlrLOP(scip) );

   /* include linear ordering file reader */
   SCIP_CALL_ERROR( SCIPincludeReaderLOP(scip) );

   /* Process command line arguments */
   SCIP_CALL_ERROR( SCIPprocessShellArguments(scip, argc, argv, "scip.set") );

   /* free SCIP */
   SCIP_CALL_ERROR( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return 0;
}
