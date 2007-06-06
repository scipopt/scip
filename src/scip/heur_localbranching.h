/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_localbranching.h,v 1.3 2007/06/06 11:25:17 bzfpfend Exp $"

/**@file   heur_localbranching.h
 * @brief  localbranching primal heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_LOCALBRANCHING_H__
#define __SCIP_HEUR_LOCALBRANCHING_H__


#include "scip/scip.h"


/** creates local branching primal heuristic and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeHeurLocalbranching(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
