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
#pragma ident "@(#) $Id: struct_primal.h,v 1.10 2005/02/14 13:35:52 bzfpfend Exp $"

/**@file   struct_primal.h
 * @brief  datastructures for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_PRIMAL_H__
#define __STRUCT_PRIMAL_H__


#include "scip/def.h"
#include "scip/type_sol.h"
#include "scip/type_primal.h"



/** primal data and solution storage */
struct Primal
{
   Longint          nsolsfound;         /**< number of primal CIP solutions found up to now */
   Longint          nbestsolsfound;     /**< number of new best primal CIP solutions found up to now */
   Real             upperbound;         /**< upper (primal) bound of CIP: objective value of best solution or user bound */
   Real             cutoffbound;        /**< upper bound for better primal solutions (if objective value is always
                                         *   integral, cutoffbound is equal to ceil(upperbound) - 1.0 (+eps) */
   SOL**            sols;               /**< primal CIP solutions */
   SOL**            existingsols;       /**< all existing primal solutions (feasible and infeasible) */
   SOL*             currentsol;         /**< internal solution for temporarily storing the current solution */
   int              solssize;           /**< size of sols array */
   int              nsols;              /**< number of primal CIP solutions stored in sols array */
   int              existingsolssize;   /**< size of existingsols array */
   int              nexistingsols;      /**< number of primal CIP solutions stored in existingsols array */
};


#endif
