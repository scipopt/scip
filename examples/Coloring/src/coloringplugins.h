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

/**@file   coloringplugins.h
 * @brief  SCIP plugins for coloring
 * @author Gerald Gamrath
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIPCOLORINGPLUGINS_H__
#define __SCIP_SCIPCOLORINGPLUGINS_H__


#include "scip/scip.h"

/* include header files here, such that the user only has to include
 * coloringplugins.h
 */
#include "scip/branch_allfullstrong.h"
#include "scip/branch_fullstrong.h"
#include "scip/branch_inference.h"
#include "scip/branch_mostinf.h"
#include "scip/branch_leastinf.h"
#include "scip/branch_pscost.h"
#include "scip/branch_random.h"
#include "scip/branch_relpscost.h"

#include "scip/cons_integral.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"

#include "scip/dialog_default.h"
#include "scip/disp_default.h"

#include "scip/heur_actconsdiving.h"
#include "scip/heur_coefdiving.h"
#include "scip/heur_crossover.h"
#include "scip/heur_dins.h"
#include "scip/heur_feaspump.h"
#include "scip/heur_fixandinfer.h"
#include "scip/heur_fracdiving.h"
#include "scip/heur_guideddiving.h"
#include "scip/heur_intdiving.h"
#include "scip/heur_intshifting.h"
#include "scip/heur_linesearchdiving.h"
#include "scip/heur_localbranching.h"
#include "scip/heur_mutation.h"
#include "scip/heur_objpscostdiving.h"
#include "scip/heur_octane.h"
#include "scip/heur_oneopt.h"
#include "scip/heur_pscostdiving.h"
#include "scip/heur_rens.h"
#include "scip/heur_rins.h"
#include "scip/heur_rootsoldiving.h"
#include "scip/heur_rounding.h"
#include "scip/heur_shifting.h"
#include "scip/heur_simplerounding.h"
#include "scip/heur_trivial.h"
#include "scip/heur_trysol.h"
#include "scip/heur_twoopt.h"
#include "scip/heur_undercover.h"
#include "scip/heur_veclendiving.h"
#include "scip/heur_zirounding.h"

#include "scip/nodesel_bfs.h"
#include "scip/nodesel_dfs.h"
#include "scip/nodesel_estimate.h"
#include "scip/nodesel_hybridestim.h"
#include "scip/nodesel_restartdfs.h"

#include "scip/prop_pseudoobj.h"
#include "scip/prop_rootredcost.h"

#include "scip/scipshell.h"

/** project plugins */
#include "cons_storeGraph.h"
#include "branch_coloring.h"
#include "branch_strongcoloring.h"
#include "reader_col.h"
#include "reader_csol.h"
#include "heur_init.h"
#include "pricer_coloring.h"

#ifdef __cplusplus
extern "C" {
#endif

/** includes default SCIP plugins into SCIP */
extern
SCIP_RETCODE SCIPincludeColoringPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
