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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_gr.c
 * @brief  Steiner tree problem file reader
 * @author Daniel Rehfeldt
 *
 * This file implements the reader used to read and write Steiner tree problems.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "probdata_stp.h"
#include "reader_gr.h"
#include "grph.h"


/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "grreader"
#define READER_DESC             "file reader for Steiner tree data format (gr)"
#define READER_EXTENSION        "gr"


/**@} */


/**@name Callback methods
 *
 * @{
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyGr)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderGr(scip) );

   return SCIP_OKAY;
}

/** problem reading method of the reader */
static
SCIP_DECL_READERREAD(readerReadGr)
{  /*lint --e{715}*/
   SCIP_RETCODE          retcode;
   SCIP_PROBDATA*        probdata;

   *result = SCIP_DIDNOTRUN;

   retcode = SCIPprobdataCreate(scip, filename);

   if( retcode == SCIP_READERROR )
      return SCIP_READERROR;

   SCIP_CALL( retcode );

   probdata = SCIPgetProbData(scip);
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || probdata == NULL )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/**@} */


/**@name Interface methods
 *
 * @{
 */


/** includes the gr file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderGr(
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

   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyGr) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadGr) );

   return SCIP_OKAY;
}

/**@} */
