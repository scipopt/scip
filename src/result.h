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

/**@file   result.h
 * @brief  result codes for SCIP callback methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __RESULT_H__
#define __RESULT_H__

/** result codes for SCIP callback methods */
enum Result
{
   SCIP_DIDNOTRUN   =   0,            /**< the method was not executed */
   SCIP_DIDNOTFIND  =   1,            /**< the method was executed, but did not have success to find anything */
   SCIP_FEASIBLE    =   2,            /**< no infeasibility could be found */
   SCIP_INFEASIBLE  =   3,            /**< an infeasibility was detected */
   SCIP_UNBOUNDED   =   4,            /**< an unboundness was detected */
   SCIP_SEPARATED   =   5,            /**< the method added a cutting plane */
   SCIP_REDUCEDDOM  =   6,            /**< the method reduced the domain of a variable */
   SCIP_BRANCHED    =   7,            /**< the method created a branching */
   SCIP_SUCCESS     =   8,            /**< the method was successfully executed */  
   SCIP_FAILURE     =   9             /**< the method's execution failed */
};
typedef enum Result RESULT;           /**< result codes for SCIP callback methods */



#endif
