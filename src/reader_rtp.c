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

/**@file   reader_rtp.c
 * @brief  RTP file reader
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "reader_rtp.h"
#include "cons_bitvar.h"
#include "cons_bitarith.h"


#define READER_NAME             "rtpreader"
#define READER_DESC             "rtp file reader"
#define READER_EXTENSION        "rtp"




/*
 * Data structures
 */

/* TODO: (optional) fill in the necessary reader data */

/** data for rtp reader */
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
DECL_READERFREE(readerFreeRtp)
{
   errorMessage("method of rtp reader not implemented yet");
   abort();

   return SCIP_OKAY;
}
#else
#define readerFreeRtp NULL
#endif


/** problem reading method of reader */
static
DECL_READERREAD(readerReadRtp)
{
   errorMessage("method of rtp reader not implemented yet");
   abort();

   return SCIP_OKAY;
}




/*
 * reader specific interface methods
 */

/** includes the rtp file reader in SCIP */
RETCODE SCIPincludeReaderRtp(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   READERDATA* readerdata;

   /* create rtp reader data */
   readerdata = NULL;
   /* TODO: (optional) create reader specific data here */
   
   /* include rtp reader */
   CHECK_OKAY( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
                  readerFreeRtp, readerReadRtp, readerdata) );

   /* add rtp reader parameters */
   /* TODO: (optional) add reader specific parameters with SCIPaddTypeParam() here */
   
   return SCIP_OKAY;
}
