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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_sepa.h,v 1.9 2005/02/04 14:27:24 bzfpfend Exp $"

/**@file   type_sepa.h
 * @brief  type definitions for separators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_SEPA_H__
#define __TYPE_SEPA_H__


typedef struct Sepa SEPA;               /**< separator */
typedef struct SepaData SEPADATA;       /**< locally defined separator data */



/** destructor of separator to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define DECL_SEPAFREE(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** initialization method of separator (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define DECL_SEPAINIT(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** deinitialization method of separator (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define DECL_SEPAEXIT(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** solving process initialization method of separator (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The separator may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define DECL_SEPAINITSOL(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** solving process deinitialization method of separator (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The separator should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define DECL_SEPAEXITSOL(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** execution method of separator
 *
 *  Searches for cutting planes. The method is called in the LP solving loop.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 *  - result          : pointer to store the result of the separation call
 *
 *  possible return values for *result (if more than one applies, the first in the list should be used):
 *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
 *  - SCIP_CONSADDED  : an additional constraint was generated
 *  - SCIP_REDUCEDDOM : a variable's domain was reduced
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 */
#define DECL_SEPAEXEC(x) RETCODE x (SCIP* scip, SEPA* sepa, RESULT* result)


#include "def.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_scip.h"


#endif
