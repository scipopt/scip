/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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


void SCIPretcodePrint(                  /**< prints error message for return code */
   FILE*            errout,             /**< file stream to write error message */
   RETCODE          retcode             /**< SCIP return code causing the error */
   )
{
   assert(errout != NULL);

   switch( retcode )
   {
   case SCIP_DIDNOTRUN:
      fprintf(errout, "the method was not executed");
      break;
   case SCIP_FAILURE:
      fprintf(errout, "the method was executed, but did not have success to find anything");
      break;
   case SCIP_BRANCHED:
      fprintf(errout, "the method created a branching");
      break;
   case SCIP_REDUCEDDOM:
      fprintf(errout, "the method reduced the domain of a variable");
      break;
   case SCIP_SEPARATED:
      fprintf(errout, "the method added a cutting plane");
      break;
   case SCIP_INFEASIBLE:
      fprintf(errout, "an infeasibility was detected");
      break;
   case SCIP_FEASIBLE:
      fprintf(errout, "no infeasibility could be found");
      break;
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
   default:
      fprintf(errout, "unknown error code");
      break;
   }
}
