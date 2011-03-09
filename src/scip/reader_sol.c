/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_sol.c
 * @ingroup FILEREADERS 
 * @brief  file reader for primal solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/reader_sol.h"


#define READER_NAME             "solreader"
#define READER_DESC             "file reader for primal solutions"
#define READER_EXTENSION        "sol"

/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeSol NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadSol)
{  /*lint --e{715}*/
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   assert(result != NULL);

   if( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
   {
      SCIPwarningMessage("reading of solution file is only possible after a problem was created\n");
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if( SCIPgetStage(scip) == SCIP_STAGE_SOLVED )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL,
         "primal solution from solution file <%s> was ignored - problem is already solved to optimality\n",
         filename);
      *result = SCIP_SUCCESS;
      return SCIP_OKAY;
   }

   /* transform the problem such that adding primal solutions is possible */
   SCIP_CALL( SCIPtransformProb(scip) );

   /* read the solution and add it to the solution pool */
   SCIP_CALL( SCIPreadSol(scip, filename) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


/** problem writing method of reader */
#define readerWriteSol NULL


/*
 * sol file reader specific interface methods
 */

/** includes the sol file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderSol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create sol reader data */
   readerdata = NULL;

   /* include sol reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeSol, readerReadSol, readerWriteSol, readerdata) );

   return SCIP_OKAY;
}

