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

/**@file   spaplugins.c
 * @brief  SCIP plugins for Sparse Approximation of Transition Networks
 * @author Leon Eifler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "spaplugins.h"

#include "scip/debug.h"
#include "scip/scipdefplugins.h"
#include "event_newsol.h"
#include "sepa_sparseapprox.h"

/** includes default plugins for coloring into SCIP */
SCIP_RETCODE SCIPincludeSpaPlugins(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPincludeReaderSpa(scip) );
   SCIP_CALL( SCIPincludeHeurSpaGreedy(scip) );
   SCIP_CALL( SCIPincludeHeurSpaswitch(scip) );
   SCIP_CALL( SCIPincludeEventHdlrNewsol(scip) );
   SCIP_CALL( SCIPincludeSepaSparseApprox(scip) );

   SCIP_CALL( SCIPaddRealParam(scip,"coherence_bound","lower bound to within-cluster coherence", NULL, FALSE, 0.1, 0.0, 1.0, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "edge_representation", "true, if the edge represantation should be used. Otherwise the bin represantation is used", NULL, FALSE, TRUE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "irrev_relaxation", "should the abs-value be inside the sum (relaxing the constraint by triangle-inequality)", NULL, FALSE, FALSE, NULL, NULL) );

   return SCIP_OKAY;
}
