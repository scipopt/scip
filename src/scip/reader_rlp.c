/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_rlp.c,v 1.2 2007/10/31 10:32:20 bzfheinz Exp $"

/**@file   reader_rlp.h
 * @brief  RLP file reader (LP format with generic variables and row names)
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/reader_lp.h"
#include "scip/reader_rlp.h"

#define READER_NAME             "rlpreader"
#define READER_DESC             "file reader for MIPs in ILOG's RLP file format"
#define READER_EXTENSION        "rlp"


/*
 * Callback methods of reader
 */

/** destructor of reader to free user data (called when SCIP is exiting) */
#define readerFreeRlp NULL


/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadRlp)
{  /*lint --e{715}*/

   SCIP_CALL( SCIPreadLp(scip, reader, filename, result) );

   return SCIP_OKAY;
}


/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteRlp)
{  /*lint --e{715}*/
   
   /* print statistics as comment to file */
   SCIPmessageFPrintInfo(file, "\\ STATISTICS\n");
   SCIPmessageFPrintInfo(file, "\\   Problem name     : %s\n", name);
   SCIPmessageFPrintInfo(file, "\\   Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      nvars, nbinvars, nintvars, nimplvars, ncontvars);
   
   SCIPmessageFPrintInfo(file, "\\   Constraints      : %d initial, %d maximal\n", startnconss, maxnconss);
   
   SCIP_CALL( SCIPwriteLp(scip, file, TRUE, objsense, vars, nvars, conss, nconss, result) );
   
   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** includes the lp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderRlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create lp reader data */
   readerdata = NULL;
   
   /* include lp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerFreeRlp, readerReadRlp, readerWriteRlp, readerdata) );
   
   return SCIP_OKAY;
}
