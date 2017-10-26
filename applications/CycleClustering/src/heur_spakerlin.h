/*
 * heur_spaswitch.h
 *
 *  Created on: Feb 23, 2016
 *      Author: bzfeifle
 */

#ifndef APPLICATIONS_SPARSEAPPROX_SRC_HEUR_SPAKERLIN_H_
#define APPLICATIONS_SPARSEAPPROX_SRC_HEUR_SPAKERLIN_H_

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_spaswitch.h
 * @ingroup PRIMALHEURISTICS
 * @brief  Improvement heuristic that trades bin-variables between clusters
 * @author Leon Eifler
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the oneopt primal heuristic and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeHeurSpakerlin(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif




#endif /* APPLICATIONS_SPARSEAPPROX_SRC_HEUR_SPAKERLIN_H_ */
