/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_result.h,v 1.7 2005/01/21 09:17:12 bzfpfend Exp $"

/**@file   type_result.h
 * @brief  result codes for SCIP callback methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_RESULT_H__
#define __TYPE_RESULT_H__

/** result codes for SCIP callback methods */
enum Result
{
   SCIP_DIDNOTRUN   =   1,            /**< the method was not executed */
   SCIP_DELAYED     =   2,            /**< the method was not executed, but should be called again later */
   SCIP_DIDNOTFIND  =   3,            /**< the method was executed, but failed finding anything */
   SCIP_FEASIBLE    =   4,            /**< no infeasibility could be found */
   SCIP_INFEASIBLE  =   5,            /**< an infeasibility was detected */
   SCIP_UNBOUNDED   =   6,            /**< an unboundness was detected */
   SCIP_CUTOFF      =   7,            /**< the current node is infeasible and can be cut off */
   SCIP_SEPARATED   =   8,            /**< the method added a cutting plane */
   SCIP_REDUCEDDOM  =   9,            /**< the method reduced the domain of a variable */
   SCIP_CONSADDED   =  10,            /**< the method added a constraint */
   SCIP_CONSCHANGED =  11,            /**< the method changed a constraint */
   SCIP_BRANCHED    =  12,            /**< the method created a branching */
   SCIP_SOLVELP     =  13,            /**< the current node's LP must be solved */
   SCIP_FOUNDSOL    =  14,            /**< the method found a feasible primal solution */
   SCIP_SUCCESS     =  15,            /**< the method was successfully executed */  
   SCIP_SUSPENDED   =  16             /**< the method interrupted its execution, but can continue if needed */
};
typedef enum Result RESULT;           /**< result codes for SCIP callback methods */



#endif
