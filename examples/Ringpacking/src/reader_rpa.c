/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_rpa.c
 * @brief  Ringpacking problem reader
 * @author Benjamin Mueller
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "reader_rpa.h"
#include "probdata_rpa.h"

/**@name Reader properties
 *
 * @{
 */

#define READER_NAME             "rpareader"
#define READER_DESC             "file reader for binpacking data format"
#define READER_EXTENSION        "rpa"

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadRpa)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   SCIP_Real* rints;
   SCIP_Real* rexts;
   int* demands;
   SCIP_Bool error;
   char name[SCIP_MAXSTRLEN];
   char buffer[SCIP_MAXSTRLEN];
   SCIP_Real width;
   SCIP_Real height;
   SCIP_Real r_int;
   SCIP_Real r_ext;
   int demand;
   int ntypes;
   int nread;
   int lineno;
   int i;

   *result = SCIP_DIDNOTRUN;
   width = -1.0;
   height = -1.0;

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   lineno = 0;
   sprintf(name, "++ uninitialized ++");
   ntypes = 0;

   /* read problem dimension */
   if( !SCIPfeof(file) )
   {
      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;
      lineno++;

      /* parse instance name line */
      nread = sscanf(buffer, "%s", name);
      if( nread != 1 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         return SCIP_READERROR;
      lineno++;

      /* parse dimension line */
      nread = sscanf(buffer, "%d %" SCIP_REAL_FORMAT " %" SCIP_REAL_FORMAT "\n", &ntypes, &width, &height);
      if( nread < 3 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }
   }

   SCIPdebugMessage("instance name = %s\n", name);
   SCIPdebugMessage("width = %e height = %e\n", MAX(width,height), MIN(width,height));

   /* allocate buffer memory for storing the demands, rints, rexts */
   SCIP_CALL( SCIPallocBufferArray(scip, &demands, ntypes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rints, ntypes) );
   SCIP_CALL( SCIPallocBufferArray(scip, &rexts, ntypes) );

   /* ring types */
   r_int = 0.0;
   r_ext = 0.0;
   demand = 0;
   i = 0;
   error = FALSE;

   while( !SCIPfeof(file) && !error )
   {
      /* get next line */
      if( SCIPfgets(buffer, (int)sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* parse the line */
      nread = sscanf(buffer, "%d %" SCIP_REAL_FORMAT " %" SCIP_REAL_FORMAT "\n", &demand, &r_int, &r_ext);
      if( nread == 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      if( r_int > r_ext )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: internal radius is greater than the external one\n");
         error = TRUE;
         break;
      }

      if( demand <= 0 )
      {
         SCIPwarningMessage(scip, "invalid input line %d in file <%s>: demand has to be positive\n");
         error = TRUE;
         break;
      }

      demands[i] = demand;
      rints[i] = r_int;
      rexts[i] = r_ext;
      ++i;

      if( i == ntypes )
         break;
   }

   if( i < ntypes )
   {
      SCIPwarningMessage(scip, "found %d different types of rings, needed %d\n", i, ntypes);
      error = TRUE;
   }

   if( !SCIPisPositive(scip, width) || !SCIPisPositive(scip, height) )
   {
      SCIPwarningMessage(scip, "non-positive width and height = (%f, %f)!\n", width, height);
      error = TRUE;
   }

   if( !error )
   {
      /* sort rings by their external radii */
      SCIPsortDownRealRealInt(rexts, rints, demands, ntypes);

      /* create and set problem data */
      SCIP_CALL( SCIPprobdataCreate(scip, filename, demands, rints, rexts, ntypes, MAX(width,height), MIN(width,height)) );
      SCIP_CALL( SCIPprobdataSetupProblem(scip) );
   }

   (void)SCIPfclose(file);
   SCIPfreeBufferArray(scip, &rints);
   SCIPfreeBufferArray(scip, &rexts);
   SCIPfreeBufferArray(scip, &demands);

   if( error )
      return SCIP_READERROR;

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** includes the rpa file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderRpa(
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

   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadRpa) );

   return SCIP_OKAY;
}

/**@} */
