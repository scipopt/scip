/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   branch_ryanfoster.h
 * @brief  Ryan/Foster branching rule
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file implements the Ryan/Foster branching rule. For more details see \ref BRANCHING page
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_RYANFOSTER_H__
#define __SCIP_BRANCH_RYANFOSTER_H__


#include "scip/scip.h"

/** creates the ryanfoster branching rule and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeBranchruleRyanFoster(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
