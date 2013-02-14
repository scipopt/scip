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

/**@file   reader_xyz.c
 * @brief  XYZ file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/reader_xyz.h"


#define READER_NAME             "xyzreader"
#define READER_DESC             "xyz file reader"
#define READER_EXTENSION        "xyz"


/*
 * Data structures
 */

/* TODO: (optional) fill in the necessary reader data */

/** data for xyz reader */
struct SCIP_ReaderData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of reader
 */


/** copy method for reader plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_READERCOPY(readerCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define readerCopyXyz NULL
#endif

/** destructor of reader to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_READERFREE(readerFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerFreeXyz NULL
#endif


/** problem reading method of reader */
#if 0
static
SCIP_DECL_READERREAD(readerReadXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerReadXyz NULL
#endif


#if 0
/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerWriteXyz NULL
#endif


/*
 * reader specific interface methods
 */

/** includes the xyz file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   SCIP_READER* reader;

   /* create xyz reader data */
   readerdata = NULL;
   /* TODO: (optional) create reader specific data here */
   
   reader = NULL;

   /* include reader */
#if 0
   /* use SCIPincludeReader() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyXyz, readerFreeXyz, readerReadXyz, readerWriteXyz, readerdata) );
#else
   /* use SCIPincludeReaderBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeReaderBasic(scip, &reader, READER_NAME, READER_DESC, READER_EXTENSION, readerdata) );

   assert(reader != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetReaderCopy(scip, reader, readerCopyXyz) );
   SCIP_CALL( SCIPsetReaderFree(scip, reader, readerFreeXyz) );
   SCIP_CALL( SCIPsetReaderRead(scip, reader, readerReadXyz) );
   SCIP_CALL( SCIPsetReaderWrite(scip, reader, readerWriteXyz) );
#endif

   /* add xyz reader parameters */
   /* TODO: (optional) add reader specific parameters with SCIPaddTypeParam() here */
   
   return SCIP_OKAY;
}
