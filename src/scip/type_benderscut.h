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

/**@file   type_benderscut.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for Benders' decomposition cuts
 * @author Stephen J. Maher
 *
 *  This file defines the interface for Benders' decomposition cuts implemented in C.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_BENDERSCUT_H__
#define __SCIP_TYPE_BENDERSCUT_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "scip/type_timing.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Benderscut SCIP_BENDERSCUT;              /**< Benders' decomposition cuts */
typedef struct SCIP_BenderscutData SCIP_BENDERSCUTDATA;      /**< locally defined Benders' decomposition cuts data */


/** copy method for compression plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benderscut           : the compression technique itself
 */
#define SCIP_DECL_BENDERSCUTCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERSCUT* benderscut)

/** destructor of Benders' decomposition cuts to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benderscut           : the compression technique itself
 */
#define SCIP_DECL_BENDERSCUTFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERSCUT* benderscut)

/** initialization method of Benders' decomposition cuts (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benderscut           : the compression technique itself
 */
#define SCIP_DECL_BENDERSCUTINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERSCUT* benderscut)

/** deinitialization method of Benders' decomposition cuts (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benderscut           : the compression technique itself
 */
#define SCIP_DECL_BENDERSCUTEXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERSCUT* benderscut)

/** solving process initialization method of Benders' decomposition cutsc (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The tree benderscutession may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benderscut           : the compression technique itself
 */
#define SCIP_DECL_BENDERSCUTINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERSCUT* benderscut)

/** solving process deinitialization method of Benders' decomposition cuts (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The Benders' decomposition cuts should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benderscut           : the compression technique itself
 */
#define SCIP_DECL_BENDERSCUTEXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERSCUT* benderscut)

/** execution method of Benders' decomposition cuts technique
 *
 *  Try to compress the current search tree. The method is called in the node processing loop.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - benders         : the Benders' decomposition structure
 *  - benderscut      : the Benders' cut structure
 *  - probnumber      : the number of the subproblem from which the cut is generated
 *  - result          : pointer to store the result of the cut algorithm
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define SCIP_DECL_BENDERSCUTEXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_BENDERS* benders, SCIP_BENDERSCUT* benderscut,\
   SCIP_SOL* sol, int probnumber, SCIP_RESULT* result)

#ifdef __cplusplus
}
#endif

#endif
