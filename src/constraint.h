/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   constraint.h
 * @brief  datastructures and methods for managing constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONSTRAINT_H__
#define __CONSTRAINT_H__

#include "scip.h"
#include "retcode.h"
#include "lp.h"

typedef struct ConsHdlr CONSHDLR;       /**< constraint handler for a specific constraint type */
typedef struct Constraint CONSTRAINT;   /**< constraint data structure */
typedef struct ConsList CONSLIST;       /**< list of constraints */
typedef void CONSDATA;                  /**< constraint type specific data; default is void */

/** initialization method of constraint handler
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSINIT(x) RETCODE x (SCIP* scip, CONSHDLR* self)

/** deinitialization method of constraint handler
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSEXIT(x) RETCODE x (SCIP* scip, CONSHDLR* self)

/** feasibility check method of constraint handler
 *  possible return values:
 *    SCIP_SUCCESS: constraint is feasible
 *    SCIP_FAILURE: constraint is infeasible
 *    neg. values : error codes
 */
#define DECL_CONSCHCK(x) RETCODE x (SCIP* scip, CONSHDLR* self, CONSTRAINT* constraint, double* psol)

/** domain propagation method of constraint handler
 *  possible return values:
 *    SCIP_SUCCESS: propagator searched and found domain reductions
 *    SCIP_FAILURE: propagator searched and did not find domain reductions
 *    SCIP_OKAY   : propagator was skipped
 *    neg. values : error codes
 */
#define DECL_CONSPROP(x) RETCODE x (SCIP* scip, CONSHDLR* self, CONSTRAINT* constraint, VARDOM* vardom)



#endif
