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
#pragma ident "@(#) $Id: reader_bpa.c,v 1.1 2010/09/09 12:21:14 bzfheinz Exp $"

/**@file   reader_bpa.h
 * @brief  binpacking problem reader file reader
 * @author Timo Berthold
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/cons_setppc.h"

#include "probdata_binpacking.h"
#include "reader_bpa.h"

#define READER_NAME             "bpareader"
#define READER_DESC             "file reader for binpacking data format"
#define READER_EXTENSION        "bpa"


/*
 * Data structures
 */

/** data for bpa reader */
struct SCIP_ReaderData
{
};

/*
 * Callback methods of reader
 */

/** copy method for reader plugins (called when SCIP copies plugins) */
static
SCIP_DECL_READERCOPY(readerCopyBpa)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(reader != NULL);
   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);

   /* call inclusion method of reader */
   SCIP_CALL( SCIPincludeReaderBpa(scip) );
 
   return SCIP_OKAY;
}

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeBpa NULL

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadBpa)
{  /*lint --e{715}*/
   SCIP_FILE* file;
   SCIP_Longint* weights;
   int* ids;
   SCIP_Bool error;

   char name[SCIP_MAXSTRLEN];
   char format[16];
   char buffer[SCIP_MAXSTRLEN];
   int capacity;
   int nitems;
   int bestsolvalue;
   int nread;
   int weight;
   int nweights;
   int lineno;

   *result = SCIP_DIDNOTRUN;

   /* open file */
   file = SCIPfopen(filename, "r");
   if( file == NULL )
   {
      SCIPerrorMessage("cannot open file <%s> for reading\n", filename);
      SCIPprintSysError(filename);
      return SCIP_NOFILE;
   }

   lineno = 0;

   /* read problem name */
   if( !SCIPfeof(file) )
   {
      /* get next line */
      if( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
         return SCIP_READERROR;
      lineno++;

      /* parse dimension line */
      sprintf(format, "%%%ds\n", SCIP_MAXSTRLEN);
      nread = sscanf(buffer, format, name);
      if( nread == 0 )
      {
         SCIPwarningMessage("invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      SCIPdebugMessage("problem name <%s>\n", name);
   }
   
   /* read problem dimension */
   if( !SCIPfeof(file) )
   {
      /* get next line */
      if( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
         return SCIP_READERROR;
      lineno++;

      /* parse dimension line */
      nread = sscanf(buffer, "%d %d %d\n", &capacity, &nitems, &bestsolvalue);
      if( nread < 2 )
      {
         SCIPwarningMessage("invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         return SCIP_READERROR;
      }

      SCIPdebugMessage("capacity = <%d>, number of items = <%d>, best known solution = <%d>\n", capacity, nitems, bestsolvalue);
   }
   

   /* allocate buffer memory for storing the weights and ids temporary */
   SCIP_CALL( SCIPallocBufferArray(scip, &weights, nitems) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ids, nitems) );
   
   /* pasre weights */
   nweights = 0;
   error = FALSE;
   
   while( !SCIPfeof(file) && !error )
   {
      /* get next line */
      if( SCIPfgets(buffer, sizeof(buffer), file) == NULL )
         break;
      lineno++;

      /* parse the line */
      nread = sscanf(buffer, "%d\n", &weight);
      if( nread == 0 )
      {
         SCIPwarningMessage("invalid input line %d in file <%s>: <%s>\n", lineno, filename, buffer);
         error = TRUE;
         break;
      }

      SCIPdebugMessage("found weight %d <%d>\n", nweights, weight);
      weights[nweights] = weight;
      ids[nweights] = nweights;
      nweights++;

      if( nweights == nitems )
         break;
   }

   if( nweights < nitems )
   {
      SCIPwarningMessage("set nitems from <%d> to <%d> since the file <%s> only contains <%d> weights\n", nitems, weights, filename, weights);
      nitems = nweights;
   }
   
   if( !error )
   {
      /* create a new problem in SCIP */  
      SCIP_CALL( SCIPprobdataCreate(scip, name, ids, weights, nitems, (SCIP_Longint)capacity) );
   }

   (void)SCIPfclose(file);
   SCIPfreeBufferArray(scip, &ids);
   SCIPfreeBufferArray(scip, &weights);

   if( error )
      return SCIP_READERROR;
   
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

#define readerWriteBpa NULL


/*
 * reader specific interface methods
 */

/** includes the bpa file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderBpa(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create binpacking reader data */
   readerdata = NULL;
   
   /* include binpacking reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyBpa, readerFreeBpa, readerReadBpa, readerWriteBpa, readerdata) );

   /* add bpa reader parameters */
   /* TODO: (optional) add reader specific parameters with SCIPaddTypeParam() here */
   
   return SCIP_OKAY;
}
