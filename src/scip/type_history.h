/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_history.h,v 1.4 2004/04/15 10:41:27 bzfpfend Exp $"

/**@file   type_history.h
 * @brief  type definitions for branching and inference history
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_HISTORY_H__
#define __TYPE_HISTORY_H__


/** branching direction for branching on variables */
enum BranchDir
{
   SCIP_BRANCHDIR_DOWNWARDS = 0,        /**< downwards branching: decreasing upper bound */
   SCIP_BRANCHDIR_UPWARDS   = 1         /**< upwards branching: increasing lower bound */
};
typedef enum BranchDir BRANCHDIR;       /**< branching direction for branching on variables */

typedef struct History HISTORY;         /**< branching and inference history information for single variable */


#endif
