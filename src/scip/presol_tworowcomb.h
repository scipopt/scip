/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   presol_tworowcomb.h
 * @ingroup PRESOLVERS
 * @brief  derive variable bounds from convex combination of two rows
 * @author Dieter Weninger
 * @author Patrick Gemander
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PRESOL_TWOROWCOMB_H__
#define __SCIP_PRESOL_TWOROWCOMB_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the tworowcombine presolver and includes it in SCIP
 *
 * @ingroup PresolverIncludes
 */
extern
SCIP_RETCODE SCIPincludePresolTworowcomb(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

/**@addtogroup PRESOLVERS
 *
 * @{
 */

/* TODO place other public methods in this group to facilitate navigation through the documentation */

/* @} */

#endif
