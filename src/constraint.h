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
#define DECL_CONSPROP(x) RETCODE x (CONSHDLR* self, SCIP* scip, MEMHDR* memhdr, CONS* cons, DOM* dom)



extern
CONS* SCIPconsCreate(                   /**< creates a constraint */
   MEM*             mem,                /**< block memory buffers */
   Bool             model,              /**< is constraint necessary for feasibility? */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata            /**< data for this specific constraint */
   );

extern
void SCIPconsFree(                      /**< frees a constraint */
   CONS**           cons,               /**< constraint to free */
   MEM*             mem                 /**< block memory buffers */
   );

extern
RETCODE SCIPconslistAdd(                /**< adds constraint to a list of constraints */
   CONSLIST**       conslist,           /**< constraint list to extend */
   MEM*             mem,                /**< block memory buffers */
   CONS*            cons                /**< constraint to add */
   );

extern
void SCIPconslistFreePart(              /**< partially unlinks and deletes the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEM*             mem,                /**< block memory buffers */
   CONSLIST*        firstkeep           /**< first constraint list entry to keep */
      );

extern
void SCIPconslistFree(                  /**< unlinks and deletes all the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEM*             mem                 /**< block memory buffers */
   );

#endif
