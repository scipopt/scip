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

/**@file   retcode.c
 * @brief  return codes for SCIP methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "retcode.h"


/** prints error message for return code */
void SCIPretcodePrint(
   FILE*            errout,             /**< file stream to write error message */
   RETCODE          retcode             /**< SCIP return code causing the error */
   )
{
   assert(errout != NULL);

   switch( retcode )
   {
   case SCIP_OKAY:
      fprintf(errout, "normal termination");
      break;
   case SCIP_ERROR:
      fprintf(errout, "unspecified error");
      break;
   case SCIP_NOMEMORY:
      fprintf(errout, "insufficient memory error");
      break;
   case SCIP_READERR:
      fprintf(errout, "file read error");
      break;
   case SCIP_NOFILE:
      fprintf(errout, "file not found error");
      break;
   case SCIP_LPERROR:
      fprintf(errout, "error in LP solver");
      break;
   case SCIP_NOPROBLEM:
      fprintf(errout, "no problem exists");
      break;
   case SCIP_INVALIDCALL:
      fprintf(errout, "method cannot be called at this time in solution process");
      break;
   case SCIP_INVALIDDATA:
      fprintf(errout, "method cannot be called with this type of data");
      break;
   case SCIP_INVALIDRESULT:
      fprintf(errout, "method returned an invalid result code");
      break;
   case SCIP_PLUGINNOTFOUND:
      fprintf(errout, "a required plugin was not found");
      break;
   default:
      fprintf(errout, "unknown error code");
      break;
   }
}
