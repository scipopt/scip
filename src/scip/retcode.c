/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   retcode.c
 * @brief  methods for return codes for SCIP methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>

#include "scip/retcode.h"
#include "scip/message.h"


/** prints error message for return code */
void SCIPretcodePrint(
   FILE*                 file,               /**< file stream to write error message */
   SCIP_RETCODE          retcode             /**< SCIP return code causing the error */
   )
{
   switch( retcode )
   {
   case SCIP_OKAY:
      SCIPmessageFPrintInfo(file, "normal termination");
      break;
   case SCIP_ERROR:
      SCIPmessageFPrintInfo(file, "unspecified error");
      break;
   case SCIP_NOMEMORY:
      SCIPmessageFPrintInfo(file, "insufficient memory error");
      break;
   case SCIP_READERROR:
      SCIPmessageFPrintInfo(file, "read error");
      break;
   case SCIP_WRITEERROR:
      SCIPmessageFPrintInfo(file, "write error");
      break;
   case SCIP_NOFILE:
      SCIPmessageFPrintInfo(file, "file not found error");
      break;
   case SCIP_FILECREATEERROR:
      SCIPmessageFPrintInfo(file, "cannot create file");
      break;
   case SCIP_LPERROR:
      SCIPmessageFPrintInfo(file, "error in LP solver");
      break;
   case SCIP_NOPROBLEM:
      SCIPmessageFPrintInfo(file, "no problem exists");
      break;
   case SCIP_INVALIDCALL:
      SCIPmessageFPrintInfo(file, "method cannot be called at this time in solution process");
      break;
   case SCIP_INVALIDDATA:
      SCIPmessageFPrintInfo(file, "method cannot be called with this type of data");
      break;
   case SCIP_INVALIDRESULT:
      SCIPmessageFPrintInfo(file, "method returned an invalid result code");
      break;
   case SCIP_PLUGINNOTFOUND:
      SCIPmessageFPrintInfo(file, "a required plugin was not found");
      break;
   case SCIP_PARAMETERUNKNOWN:
      SCIPmessageFPrintInfo(file, "the parameter with the given name was not found");
      break;
   case SCIP_PARAMETERWRONGTYPE:
      SCIPmessageFPrintInfo(file, "the parameter is not of the expected type");
      break;
   case SCIP_PARAMETERWRONGVAL:
      SCIPmessageFPrintInfo(file, "the value is invalid for the given parameter");
      break;
   case SCIP_KEYALREADYEXISTING:
      SCIPmessageFPrintInfo(file, "the given key is already existing in table");
      break;
   case SCIP_MAXDEPTHLEVEL:
      SCIPmessageFPrintInfo(file, "maximal branching depth level exceeded");
      break;
   default:
      SCIPmessageFPrintInfo(file, "unknown error code");
      break;
   }
}
