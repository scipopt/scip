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
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_stp.c
 * @brief  Steiner Tree problem reader file reader
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_stp.h"
#include "reader_stp.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "stpreader"
#define READER_DESC             "file reader for steiner tree data format"
#define READER_EXTENSION        "stp"


#define   DEFAULT_COMPCENTRAL  1
#define   DEFAULT_EMITGRAPH    FALSE
#define   DEFAULT_REDUCTION    4

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadStp)
{
   SCIP_RETCODE          retcode;
   SCIP_PROBDATA*        probdata;
   char                  mode;

   *result = SCIP_DIDNOTRUN;

   /* get solving mode parameter */
   SCIP_CALL( SCIPgetCharParam(scip, "stp/mode", &mode) );

   retcode = SCIPprobdataCreate(scip, filename);

   if( retcode == SCIP_READERROR )
      return SCIP_READERROR;

   SCIP_CALL( retcode );

   probdata = SCIPgetProbData(scip);
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT ||  probdata == NULL ) /*|| SCIPprobdataGetGraph(SCIPgetProbData(scip)) == NULL )*/
      return SCIP_READERROR;
   else if(SCIPprobdataGetGraph(probdata) != NULL && mode == 'p'){
      printf("activate pricer \n");
      SCIP_CALL( SCIPsetBoolParam(scip, "lp/disablecutoff", TRUE) );
      SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "stp")) );
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the stp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* valid values for user parameter 'stp/mode' */
   char vals [3] = {'c', 'f', 'p'};

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadStp) );

   /* include user parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "stp/compcentral",
         "Comp. Central Term: 0 disable, 1 max. degree, 2 min. dist. sum to all terminals, 3 min. max. dist., 4 min. dist to all nodes",
         NULL, FALSE, DEFAULT_COMPCENTRAL, 0, 4, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "stp/reduction",
         "Reduction: 0 disable, 5 maximum",
         NULL, FALSE, DEFAULT_REDUCTION, 0, 5, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "stp/emitgraph",
         "Emit graph",
         NULL, FALSE, DEFAULT_EMITGRAPH, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "stp/bigt",
         "use 'T' model", NULL, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "stp/printGraph",
         "print the graph before and after the presolving", NULL, FALSE, FALSE, NULL, NULL) );

   SCIP_CALL( SCIPaddCharParam( scip,
	 "stp/mode",
         "Solving mode: 'c'ut, 'f'low ,'p'rice",
         NULL, FALSE, 'c', vals, NULL, NULL) );


   return SCIP_OKAY;
}

/**@} */
