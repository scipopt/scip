/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_octane.h,v 1.6.2.1 2009/06/19 07:53:43 bzfwolte Exp $"

/**@file   heur_octane.h
 * @brief  octane primal heuristic based on Balas, Ceria, Dawande, Margot, and Pataki
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_OCTANE_H__
#define __SCIP_HEUR_OCTANE_H__


#include "scip/scip.h"


/** creates the octane primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurOctane(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
