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


typedef struct ConsHdlr CONSHDLR;       /**< constraint handler for a specific constraint type */
typedef struct Cons CONS;               /**< constraint data structure */
typedef struct ConsList CONSLIST;       /**< list of constraints */
typedef void CONSHDLRDATA;              /**< constraint handler data */
typedef void CONSDATA;                  /**< constraint type specific data; default is void */


#include "scip.h"
#include "retcode.h"
#include "mem.h"
#include "lp.h"

/** initialization method of constraint handler
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSINIT(x) RETCODE x (CONSHDLR* self, SCIP* scip, MEMHDR* memhdr)

/** deinitialization method of constraint handler
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSEXIT(x) RETCODE x (CONSHDLR* self, SCIP* scip, MEMHDR* memhdr)

/** frees specific constraint data
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSFREE(x) RETCODE x (CONSHDLR* self, MEMHDR* memhdr, CONSDATA* consdata)

/** feasibility check method of constraint handler
 *  possible return values:
 *    SCIP_SUCCESS: constraint is feasible
 *    SCIP_FAILURE: constraint is infeasible
 *    neg. values : error codes
 */
#define DECL_CONSCHCK(x) RETCODE x (CONSHDLR* self, SCIP* scip, MEMHDR* memhdr, CONS* cons, Real* psol)

/** domain propagation method of constraint handler
 *  possible return values:
 *    SCIP_SUCCESS: propagator searched and found domain reductions
 *    SCIP_FAILURE: propagator searched and did not find domain reductions
 *    SCIP_OKAY   : propagator was skipped
 *    neg. values : error codes
 */
#define DECL_CONSPROP(x) RETCODE x (CONSHDLR* self, SCIP* scip, MEMHDR* memhdr, CONS* cons)



/** constraint handler */
struct ConsHdlr
{
   const char*      name;               /**< name of constraint handler */
   DECL_CONSINIT((*consinit));          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit));          /**< deinitialise constraint handler */
   DECL_CONSFREE((*consfree));          /**< frees specific constraint data */
   DECL_CONSCHCK((*conschck));          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop));          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata;       /**< constraint handler data */
};



extern
RETCODE SCIPconsCreate(                 /**< creates a constraint */
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   Bool             original,           /**< belongs constraint to the original problem formulation? */
   Bool             model,              /**< is constraint necessary for feasibility? */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata            /**< data for this specific constraint */
   );

extern
void SCIPconsFree(                      /**< frees a constraint */
   CONS**           cons,               /**< constraint to free */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
void SCIPconsCapture(                   /**< increases usage counter of constraint */
   CONS*            cons                /**< constraint */
   );

extern
void SCIPconsRelease(                   /**< decreases usage counter of constraint, and frees memory if necessary */
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
RETCODE SCIPconslistAdd(                /**< adds constraint to a list of constraints and captures it */
   CONSLIST**       conslist,           /**< constraint list to extend */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   );

extern
void SCIPconslistFreePart(              /**< partially unlinks and releases the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory */
   CONSLIST*        firstkeep           /**< first constraint list entry to keep */
   );

extern
void SCIPconslistFree(                  /**< unlinks and releases all the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr              /**< block memory */
   );

#endif
