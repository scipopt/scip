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

/**@file   retcode.h
 * @brief  return codes for SCIP methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __RETCODE_H__
#define __RETCODE_H__

#include <stdio.h>

/** return codes for SCIP methods: negative return codes are errors */
enum Retcode
{
   SCIP_DIDNOTRUN   =   8,              /**< the method was not executed */
   SCIP_FAILURE     =   7,              /**< the method was executed, but did not have success to find anything */
   SCIP_BRANCHED    =   6,              /**< the method created a branching */
   SCIP_REDUCEDDOM  =   5,              /**< the method reduced the domain of a variable */
   SCIP_SEPARATED   =   4,              /**< the method added a cutting plane */
   SCIP_UNBOUNDED   =   3,              /**< an unboundness was detected */
   SCIP_INFEASIBLE  =   2,              /**< an infeasibility was detected */
   SCIP_FEASIBLE    =   1,              /**< no infeasibility could be found */
   SCIP_OKAY        =   0,              /**< normal termination */
   SCIP_ERROR       =  -1,              /**< ERROR: unspecified error */
   SCIP_NOMEMORY    =  -2,              /**< ERROR: insufficient memory error */
   SCIP_READERR     =  -3,              /**< ERROR: file read error */
   SCIP_NOFILE      =  -4,              /**< ERROR: file not found error */
   SCIP_LPERROR     =  -5,              /**< ERROR: error in LP solver */
   SCIP_NOPROBLEM   =  -6,              /**< ERROR: no problem exists */
   SCIP_INVALIDCALL =  -7,              /**< ERROR: method cannot be called at this time in solution process */
   SCIP_INVALIDDATA =  -8               /**< ERROR: error in input data */
};
typedef enum Retcode RETCODE;           /**< return code for SCIP method */



extern
void SCIPretcodePrint(                  /**< prints error message for return code */
   FILE*            errout,             /**< file stream to write error message */
   RETCODE          retcode             /**< SCIP return code causing the error */
   );

#endif
