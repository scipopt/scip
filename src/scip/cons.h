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

/**@file   cons.h
 * @brief  datastructures and methods for managing constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_H__
#define __CONS_H__


typedef struct ConsHdlr CONSHDLR;       /**< constraint handler for a specific constraint type */
typedef struct Cons CONS;               /**< constraint data structure */
typedef struct ConsList CONSLIST;       /**< list of constraints */
typedef void CONSHDLRDATA;              /**< constraint handler data */
typedef struct ConsData CONSDATA;       /**< locally defined constraint type specific data */

/** initialization method of constraint handler
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSINIT(x) RETCODE x (CONSHDLR* self, SCIP* scip)

/** deinitialization method of constraint handler
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSEXIT(x) RETCODE x (CONSHDLR* self, SCIP* scip)

/** frees specific constraint data
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSFREE(x) RETCODE x (CONSHDLR* self, SCIP* scip, CONSDATA** consdata)

/** transforms constraint data into data belonging to the transformed problem
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_CONSTRAN(x) RETCODE x (CONSHDLR* self, SCIP* scip, CONSDATA* sourcedata, CONSDATA** targetdata)

/** feasibility check method of constraint handler
 *  possible return values:
 *    SCIP_SUCCESS: constraint is feasible
 *    SCIP_FAILURE: constraint is infeasible
 *    neg. values : error codes
 */
#define DECL_CONSCHCK(x) RETCODE x (CONSHDLR* self, SCIP* scip, CONS* cons, Real* psol)

/** domain propagation method of constraint handler
 *  possible return values:
 *    SCIP_SUCCESS: propagator searched and found domain reductions
 *    SCIP_FAILURE: propagator searched and did not find domain reductions
 *    SCIP_OKAY   : propagator was skipped
 *    neg. values : error codes
 */
#define DECL_CONSPROP(x) RETCODE x (CONSHDLR* self, SCIP* scip, CONS* cons)


#include "scip.h"
#include "retcode.h"
#include "mem.h"
#include "lp.h"


/** linked list of constraints */
struct ConsList
{
   CONS*            cons;               /**< pointer to constraint data structure */
   CONSLIST*        next;               /**< next list entry */
};



extern
RETCODE SCIPconshdlrCreate(             /**< creates a constraint handler */
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSFREE((*consfree)),          /**< frees specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transforms constraint data into data belonging to the transformed problem */
   DECL_CONSCHCK((*conschck)),          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop)),          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   );

extern
RETCODE SCIPconshdlrFree(               /**< frees memory of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer to constraint handler data structure */
   );

extern
RETCODE SCIPconshdlrInit(               /**< initializes constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPconshdlrExit(               /**< calls exit method of constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
const char* SCIPconshdlrGetName(        /**< gets name of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handlert */
   );

extern
Bool SCIPconshdlrIsInitialized(         /**< is constraint handler initialized? */
   CONSHDLR*        conshdlr            /**< constraint handlert */
   );

extern
RETCODE SCIPconsCreate(                 /**< creates a constraint */
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   Bool             ismodel,            /**< is constraint necessary for feasibility? */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata            /**< data for this specific constraint */
   );

extern
RETCODE SCIPconsFree(                   /**< frees a constraint */
   CONS**           cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPconslistAdd(                /**< adds constraint to a list of constraints and captures it */
   CONSLIST**       conslist,           /**< constraint list to extend */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   );

extern
RETCODE SCIPconslistFreePart(           /**< partially unlinks and frees the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONSLIST*        firstkeep           /**< first constraint list entry to keep */
   );

extern
RETCODE SCIPconslistFree(               /**< unlinks and frees all the constraints in the list */
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPconsTransform(              /**< copies original constraint into transformed constraint */
   CONS*            origcons,           /**< original constraint */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONS**           transcons           /**< pointer to transformed constraint */
   );

#endif
