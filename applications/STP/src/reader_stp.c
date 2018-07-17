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

/**@file   reader_stp.c
 * @brief  Steiner tree problem file reader
 * @author Gerald Gamrath
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 * @author Michael Winkler
 *
 * This file implements the reader used to read and write Steiner tree problems.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_stp.h"
#include "reader_stp.h"
#include "grph.h"


/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "stpreader"
#define READER_DESC             "file reader for steiner tree data format"
#define READER_EXTENSION        "stp"

#define   DEFAULT_COMPCENTRAL  1             /**< selection type for the root (for undirected STPs) */
#define   DEFAULT_EMITGRAPH    FALSE         /**< emit graph? */
#define   DEFAULT_COUNTPRESOLTIME  TRUE      /**< count presolving time as part of overall solution time? */
#define   DEFAULT_REDUCTION    2             /**< reduction mode to apply */
#define   DEFAULT_SYMCONS      2             /**< symmetry constraints */
#define   DEFAULT_CYCLECONS    2             /**< cycle constraints */
#define   DEFAULT_MINELIMS     3             /**< minimal number of eliminations to be achieved for reiteration of reduction methods */
#define   DEFAULT_PRETIMELIMIT -1.0          /**< presolving time limit */

#define STP_MODES "cfp" /**< valid values for user parameter 'stp/mode' */

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyStp)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderStp(scip) );

   return SCIP_OKAY;
}

/** problem reading method of the reader */
static
SCIP_DECL_READERREAD(readerReadStp)
{  /*lint --e{715}*/
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
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || probdata == NULL )
      return SCIP_READERROR;
   else if(SCIPprobdataGetGraph(probdata) != NULL && mode == 'p')
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "activate pricer\n");
      SCIP_CALL( SCIPactivatePricer(scip, SCIPfindPricer(scip, "stp")) );
   }

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** problem writing method of the reader */
static
SCIP_DECL_READERWRITE(readerWriteStp)
{  /*lint --e{715}*/
   const GRAPH* graph;
   SCIP_Real offset;

   /* get the graph of the problem */
   graph = SCIPprobdataGetGraph(probdata);

   /* get the offset of the problem */
   offset = SCIPprobdataGetOffset(scip);

   /* save the graph in a .stp file */
   SCIPwriteStp(scip, graph, file, offset);

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** include user parameters */
SCIP_RETCODE SCIPStpReaderIncludeParams(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   /* include user parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "stp/compcentral",
         "Comp. Central Term: 0 disable, 1 max. degree, 2 min. dist. sum to all terminals, 3 min. max. dist., 4 min. dist to all nodes",
         NULL, FALSE, DEFAULT_COMPCENTRAL, 0, 4, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "stp/reduction",
         "Reduction: 0 disable, 1 diminish, 2 default",
         NULL, FALSE, DEFAULT_REDUCTION, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "stp/usesymcons",
         "Use symmetry constraints (PC, MW): 0 never, 1 always, 2 problem specific",
         NULL, FALSE, DEFAULT_SYMCONS, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "stp/usecyclecons",
         "Use 2-cycle constraints (PC): 0 never, 1 always, 2 problem specific",
         NULL, FALSE, DEFAULT_CYCLECONS, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
         "stp/minelims",
         "minimal number of eliminations per reduction method",
         NULL, FALSE, DEFAULT_MINELIMS, 0, 10000, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip,
         "stp/pretimelimit",
         "presolving time limit",
         NULL, FALSE, DEFAULT_PRETIMELIMIT, -1.0, SCIPinfinity(scip), NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "stp/countpresoltime",
         "count presolving time to solving time?",
         NULL, FALSE, DEFAULT_COUNTPRESOLTIME, NULL, NULL) );

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

   SCIP_CALL( SCIPaddCharParam(scip,
    "stp/mode",
         "Solving mode: 'c'ut, 'f'low ,'p'rice",
         NULL, FALSE, 'c', STP_MODES, NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
         "stp/logfile",
         "log file in DIMACS challenge format; use_probname for using problem name",
         NULL, FALSE, "",
         NULL, NULL) );

   SCIP_CALL( SCIPaddStringParam(scip,
         "stp/intlogfile",
         "log file for intermediate solutions; use_probname for using problem name",
         NULL, FALSE, "",
         NULL, NULL) );

   return SCIP_OKAY;
}


/** includes the stp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderStp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create reader data */
   readerdata = NULL;

   /* include reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyStp) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadStp) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteStp) );

   SCIP_CALL( SCIPStpReaderIncludeParams(scip) );

   return SCIP_OKAY;
}

/**@} */
