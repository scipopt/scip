/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heuristics.h
 * @brief  methods commonly used by primal heuristics
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEURISTICS_H__
#define __SCIP_HEURISTICS_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** get a sub-SCIP copy of the transformed problem */
extern
SCIP_RETCODE SCIPheuristicsInstanceCopy(
   SCIP*                 sourcescip,         /**< source SCIP data structure */
   SCIP*                 subscip,            /**< sub-SCIP used by the heuristic */
   SCIP_HASHMAP*         varmap,             /**< a hashmap to store the mapping of source variables to the corresponding
                                              *   target variables */
   const char*           suffix,             /**< suffix for the problem name */
   SCIP_VAR**            fixedvars,          /**< source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Real*            fixedvals,          /**< array of fixing values for target SCIP variables, or NULL */
   int                   nfixedvars,         /**< number of source variables whose copies should be fixed in the target SCIP environment, or NULL */
   SCIP_Bool             uselprows,          /**< should the linear relaxation of the problem defined by LP rows be copied? */
   SCIP_Bool             copycuts,           /**< should cuts be copied (only if uselprows == FALSE) */
   SCIP_Bool*            success             /**< was the copying successful? */
   );

#ifdef __cplusplus
}
#endif

#endif
