/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   LOP/src/cmain.c
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
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
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
