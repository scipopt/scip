/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   reader_xxx.c
 * @brief  XXX file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "reader_xxx.h"


#define READER_NAME             "xxxreader"
#define READER_DESC             "xxx file reader"
#define READER_EXTENSION        "xxx"




/*
 * Data structures
 */

/* TODO: (optional) fill in the necessary reader data */

/** data for xxx reader */
struct ReaderData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */





/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#if 0
static
DECL_READERFREE(readerFreeXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx reader not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerFreeXxx NULL
#endif


/** problem reading method of reader */
static
DECL_READERREAD(readerReadXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx reader not implemented yet");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}




/*
 * reader specific interface methods
 */

/** includes the xxx file reader in SCIP */
RETCODE SCIPincludeReaderXxx(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   READERDATA* readerdata;

   /* create xxx reader data */
   readerdata = NULL;
   /* TODO: (optional) create reader specific data here */
   
   /* include xxx reader */
   CHECK_OKAY( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
                  readerFreeXxx, readerReadXxx, readerdata) );

   /* add xxx reader parameters */
   /* TODO: (optional) add reader specific parameters with SCIPaddTypeParam() here */
   
   return SCIP_OKAY;
}
