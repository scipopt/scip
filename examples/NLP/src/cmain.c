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
 * @brief  main file for C++ TSP example using SCIP as a callable library
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 * Description...
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* include SCIP components */
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/* include NLP specific components */
#include "cons_polynomial.h"

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

   /* include NLP specific plugins */
   SCIP_CALL( SCIPincludeConshdlrPolynomial(scip) ); 

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, "scipnlp.set") );


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
