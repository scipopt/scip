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

/**@file   bendersdefcuts.c
 * @brief  default cuts for Benders' decomposition
 * @author Stephen J. Maher
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BENDERSDEFCUTS_H__
#define __SCIP_BENDERSDEFCUTS_H__

/* include header files here, such that the user only has to include bendersdefcuts.h */
#include "scip/benderscut_feas.h"
#include "scip/benderscut_feasalt.h"
#include "scip/benderscut_int.h"
#include "scip/benderscut_nogood.h"
#include "scip/benderscut_opt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes default Benders' decomposition cuts plugins into SCIP and the associated Benders' decomposition */
SCIP_EXPORT
SCIP_RETCODE SCIPincludeBendersDefaultCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition struture */
   );

#ifdef __cplusplus
}
#endif

#endif
