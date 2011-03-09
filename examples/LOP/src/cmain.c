/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  main file for linear ordering example
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <scip/scip.h>
#include <scip/scipdefplugins.h>

#include "probdata_lop.h"
#include "cons_linearordering.h"


/** read parameters from a file */
static
SCIP_RETCODE readParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< parameter file name, or NULL */
   )
{
   if (filename != NULL)
   {
      if( SCIPfileExists(filename) )
      {
         printf("reading parameter file <%s>\n", filename);
         SCIP_CALL( SCIPreadParams(scip, filename) );
      }
      else
      {
         printf("parameter file <%s> not found.\n", filename);
	 return SCIP_NOFILE;
      }
   }
   else
   {
      printf("parameter file empty.\n");
      return SCIP_NOFILE;
   }

   return SCIP_OKAY;
}



/** main function, which starts the solution of the linear ordering problem */
int main(
   int    argc,
   char** argv
   )
{
   SCIP* scip = NULL;

   /* check paramters */
   if (argc < 2 || argc > 3)
   {
      printf("usage: %s <LOP file> [<parameter file>]\n", argv[0]);
      return 1;
   }

   /* output version information */
   printf("Solving the linear ordering problem using SCIP.\n\n");

   SCIPprintVersion(NULL);
   printf("\n");

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   /* include linear ordering constraint handler */
   SCIP_CALL( SCIPincludeConshdlrLinearOrdering(scip) );

   /* read parameters if requested */
   if( argc == 3 )
      SCIP_CALL( readParams(scip, argv[2]) );

   /* read problem data */
   SCIP_CALL( LOPcreateProb(scip, argv[1]) );

   /* generate linear ordering model */
   SCIP_CALL( LOPgenerateModel(scip) );

   /* print model */
   if ( LOPgetNElements(scip) <= 10 )
   {
      SCIP_CALL( SCIPprintOrigProblem(scip, NULL, NULL, FALSE) );
      SCIPinfoMessage(scip, NULL, "\n");
   }

   /* solve the model */
   SCIP_CALL( SCIPsolve(scip) );

   /* print statistics */
   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   /* SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );*/
   SCIP_CALL( LOPevalSolution(scip) );

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return 0;
}
