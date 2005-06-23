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
#pragma ident "@(#) $Id: heur_octane.h,v 1.1 2005/06/23 16:30:04 bzfberth Exp $"

/**@file   heur_octane.h
 * @brief  octane primal heuristic based on Balas, Ceria, Dawande, Margot, and Pataki
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __HEUR_OCTANE_H__
#define __HEUR_OCTANE_H__


#include "scip/scip.h"


/** creates the octane primal heuristic and includes it in SCIP */
extern
RETCODE SCIPincludeHeurOctane(
   SCIP*            scip                /**< SCIP data structure */
   );

#endif
