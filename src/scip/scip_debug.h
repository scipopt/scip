/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_debug.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for debugging
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Leona Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_DEBUG_H__
#define __SCIP_SCIP_DEBUG_H__


#include "scip/def.h"
#include "scip/type_scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup DebugSolutionMethods
 *
 * @{
 */

/** enable debug solution mechanism
 *
 *  the debug solution mechanism allows to trace back the invalidation of
 *  a debug solution during the solution process of SCIP. It must be explicitly
 *  enabled for the SCIP data structure.
 *
 *  @see debug.h for more information on debug solution mechanism
 */
SCIP_EXPORT
void SCIPenableDebugSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** disable solution debugging mechanism
 *
 *  @see debug.h for more information on debug solution mechanism
 */
SCIP_EXPORT
void SCIPdisableDebugSol(
   SCIP*                 scip                /**< SCIP data structure */
   );

/**@} */

#ifdef __cplusplus
}
#endif

#endif
