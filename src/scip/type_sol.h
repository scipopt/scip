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
#pragma ident "@(#) $Id: type_sol.h,v 1.3 2004/02/05 14:12:45 bzfpfend Exp $"

/**@file   type_sol.h
 * @brief  type definitions for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_SOL_H__
#define __TYPE_SOL_H__


/** origin of solution: where to retrieve uncached elements */
enum SolOrigin
{
   SCIP_SOLORIGIN_ZERO      = 0,        /**< all non-cached elements in solution are equal to zero */
   SCIP_SOLORIGIN_LPSOL     = 1,        /**< all non-cached elements in solution are equal to current LP solution */
   SCIP_SOLORIGIN_PSEUDOSOL = 2         /**< all non-cached elements in solution are equal to current pseudo solution */
};
typedef enum SolOrigin SOLORIGIN;

typedef struct Sol SOL;                 /**< primal CIP solution */


#endif
