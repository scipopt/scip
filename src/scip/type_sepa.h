/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_sepa.h,v 1.2 2003/12/08 13:24:54 bzfpfend Exp $"

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

/** initialization method of separator (called when problem solving starts)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define DECL_SEPAINIT(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** deinitialization method of separator (called when problem solving exits)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 */
#define DECL_SEPAEXIT(x) RETCODE x (SCIP* scip, SEPA* sepa)

/** execution method of separator
 *
 *  Searches for cutting planes. The method is called in the LP solving loop.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sepa            : the separator itself
 *  - result          : pointer to store the result of the separation call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : at least one unmodifiable row is infeasible in the variable's bounds -> node is infeasible
 *  - SCIP_SEPARATED  : a cutting plane was generated
 *  - SCIP_REDUCEDDOM : no cutting plane was generated, but at least one domain was reduced
 *  - SCIP_CONSADDED  : no cutting plane or domain reductions, but at least one additional constraint was generated
 *  - SCIP_DIDNOTFIND : the separator searched, but did not find a feasible cutting plane
 *  - SCIP_DIDNOTRUN  : the separator was skipped
 */
#define DECL_SEPAEXEC(x) RETCODE x (SCIP* scip, SEPA* sepa, RESULT* result)


#include "def.h"
#include "type_retcode.h"
#include "type_result.h"
#include "type_scip.h"


#endif
