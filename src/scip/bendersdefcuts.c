/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   bendersdefcuts.c
 * @brief  default cuts for Benders' decomposition
 * @author Stephen J. Maher
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/bendersdefcuts.h"

/** includes default SCIP plugins into SCIP */
SCIP_RETCODE SCIPincludeBendersDefaultCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_BENDERS*         benders             /**< Benders' decomposition struture */
   )
{
   SCIP_CALL( SCIPincludeBenderscutOpt(scip, benders) );
   SCIP_CALL( SCIPincludeBenderscutFeas(scip, benders) );

   return SCIP_OKAY;
}
