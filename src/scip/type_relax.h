/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_relax.h,v 1.8 2005/07/15 17:20:24 bzfpfend Exp $"

/**@file   type_relax.h
 * @brief  type definitions for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_RELAX_H__
#define __SCIP_TYPE_RELAX_H__


typedef struct Relax RELAX;             /**< relaxator */
typedef struct RelaxData RELAXDATA;     /**< locally defined relaxator data */



/** destructor of relaxator to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define DECL_RELAXFREE(x) RETCODE x (SCIP* scip, RELAX* relax)

/** initialization method of relaxator (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define DECL_RELAXINIT(x) RETCODE x (SCIP* scip, RELAX* relax)

/** deinitialization method of relaxator (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define DECL_RELAXEXIT(x) RETCODE x (SCIP* scip, RELAX* relax)

/** solving process initialization method of relaxator (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The relaxator may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define DECL_RELAXINITSOL(x) RETCODE x (SCIP* scip, RELAX* relax)

/** solving process deinitialization method of relaxator (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The relaxator should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 */
#define DECL_RELAXEXITSOL(x) RETCODE x (SCIP* scip, RELAX* relax)

/** execution method of relaxator
 *
 *  The method is called in the node processing loop. It solves the current subproblem's relaxation.
 *  Like the LP relaxation, the relaxator should only operate on COLUMN variables.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - relax           : the relaxator itself
 *  - result          : pointer to store the result of the relaxation call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_SUCCESS    : the relaxator solved the relaxation and should not be called again on the same relaxation
 *  - SCIP_SUSPENDED  : the relaxator interrupted its solving process to wait for additional input (e.g. cutting
 *                      planes); however, it is able to continue the solving in order to improve the dual bound
 *  - SCIP_DIDNOTRUN  : the relaxator was skipped
 */
#define DECL_RELAXEXEC(x) RETCODE x (SCIP* scip, RELAX* relax, RESULT* result)


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_scip.h"


#endif
