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

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadStp)
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;

   *result = SCIP_DIDNOTRUN;

   retcode = SCIPprobdataCreate(scip, filename);

   if( retcode == SCIP_READERROR )
      return SCIP_READERROR;

   SCIP_CALL( retcode );

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetProbData(scip) == NULL ) /*|| SCIPprobdataGetGraph(SCIPgetProbData(scip)) == NULL )*/
      return SCIP_READERROR;

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

   /* create binpacking reader data */
   readerdata = NULL;

   /* include binpacking reader */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );
   assert(reader != NULL);

   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadStp) );

   return SCIP_OKAY;
}

/**@} */
