/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: retcode.c,v 1.12 2003/11/21 10:35:39 bzfpfend Exp $"

/**@file   retcode.c
 * @brief  return codes for SCIP methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "retcode.h"


/** prints error message for return code */
void SCIPretcodePrint(
   FILE*            file,               /**< file stream to write error message */
   RETCODE          retcode             /**< SCIP return code causing the error */
   )
{
   assert(file != NULL);

   switch( retcode )
   {
   case SCIP_OKAY:
      fprintf(file, "normal termination");
      break;
   case SCIP_ERROR:
      fprintf(file, "unspecified error");
      break;
   case SCIP_NOMEMORY:
      fprintf(file, "insufficient memory error");
      break;
   case SCIP_READERROR:
      fprintf(file, "file read error");
      break;
   case SCIP_WRITEERROR:
      fprintf(file, "file write error");
      break;
   case SCIP_NOFILE:
      fprintf(file, "file not found error");
      break;
   case SCIP_FILECREATEERROR:
      fprintf(file, "cannot create file");
      break;
   case SCIP_LPERROR:
      fprintf(file, "error in LP solver");
      break;
   case SCIP_NOPROBLEM:
      fprintf(file, "no problem exists");
      break;
   case SCIP_INVALIDCALL:
      fprintf(file, "method cannot be called at this time in solution process");
      break;
   case SCIP_INVALIDDATA:
      fprintf(file, "method cannot be called with this type of data");
      break;
   case SCIP_INVALIDRESULT:
      fprintf(file, "method returned an invalid result code");
      break;
   case SCIP_PLUGINNOTFOUND:
      fprintf(file, "a required plugin was not found");
      break;
   case SCIP_PARAMETERUNKNOWN:
      fprintf(file, "the parameter with the given name was not found");
      break;
   case SCIP_PARAMETERWRONGTYPE:
      fprintf(file, "the parameter is not of the expected type");
      break;
   case SCIP_PARAMETERWRONGVAL:
      fprintf(file, "the value is invalid for the given parameter");
      break;
   case SCIP_KEYALREADYEXISTING:
      fprintf(file, "the given key is already existing in table");
      break;
   case SCIP_PARSEERROR:
      fprintf(file, "invalid input given to the parser");
      break;
   case SCIP_MAXDEPTHLEVEL:
      fprintf(file, "maximal branching depth level exceeded");
      break;
   default:
      fprintf(file, "unknown error code");
      break;
   }
}
