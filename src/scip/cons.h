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
typedef struct ConsHdlrData CONSHDLRDATA;/**< constraint handler data */
typedef struct ConsData CONSDATA;       /**< locally defined constraint type specific data */


/** destructor of constraint handler to free user data (called when SCIP is exiting)
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 */
#define DECL_CONSFREE(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip)

/** initialization method of constraint handler (called at problem creation)
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 */
#define DECL_CONSINIT(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip)

/** deinitialization method of constraint handler (called at problem destruction)
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 */
#define DECL_CONSEXIT(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip)

/** frees specific constraint data
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    consdata        : pointer to the constraint data to free
 */
#define DECL_CONSDELE(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip, CONSDATA** consdata)

/** transforms constraint data into data belonging to the transformed problem
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    sourcedata      : constraint data to transform
 *    targetdata      : pointer to constraint data where to store transformed data
 */
#define DECL_CONSTRAN(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip, CONSDATA* sourcedata, CONSDATA** targetdata)

/** separation method of constraint handler
 *
 *  Separates all constraints of the constraint handler. The method is called in the LP solution loop,
 *  which means that a valid LP solution exists.
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    result          : pointer to store the result of the separation call
 *
 *  possible return values for *result:
 *    SCIP_SEPARATED  : at least one cutting plane was generated
 *    SCIP_DIDNOTFIND : the separator searched, but didn't found a cutting plane
 *    SCIP_DIDNOTRUN  : the separator was skipped
 */
#define DECL_CONSSEPA(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip, CONS** conss, int nconss, RESULT* result)

/** constraint enforcing method of constraint handler for LP solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was solved.
 *  The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching, reducing a variable's domain to exclude the solution or separating the solution with a valid
 *  cutting plane.
 *
 *  The enforcing methods of the active constraint handlers are called in decreasing order of their enforcing
 *  priorities until the first constraint handler returned with the value SCIP_BRANCHED, SCIP_REDUCEDDOM, or
 *  SCIP_SEPARATED. 
 *  The integrality constraint handler has an enforcing priority of zero. A constraint handler which can
 *  (or wants) to enforce its constraints only for integral solutions should have a negative enforcing priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to incorporate its own branching strategy even on non-integral
 *  solutions must have an enforcing priority greater than zero (e.g. the SOS-constraint incorporates
 *  SOS-branching on non-integral solutions).
 *  If the solution is integral and one of the constraints of the constraint handler is violated, the
 *  constraint handler has to branch, to reduce a variable's domain, or to create a cutting plane -- otherwise,
 *  the infeasibility cannot be resolved.
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result:
 *    SCIP_BRANCHED   : at least one constraint is infeasible, and branching was applied
 *    SCIP_REDUCEDDOM : at least one constraint is infeasible, and a domain was reduced
 *    SCIP_SEPARATED  : at least one constraint is infeasible, and a cutting plane was generated
 *    SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *    SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define DECL_CONSENLP(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip, CONS** conss, int nconss, RESULT* result)

/** constraint enforcing method of constraint handler for pseudo solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was not solved.
 *  The pseudo solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching or reducing a variable's domain to exclude the solution. Separation is not possible, since the
 *  LP is not processed at the current node. All LP informations like LP solution, slack values, or reduced costs
 *  are invalid and must not be accessed.
 *
 *  Like in the enforcing method for LP solutions, the enforcing methods of the active constraint handlers are
 *  called in decreasing order of their enforcing priorities until the first constraint handler returned with
 *  the value SCIP_BRANCHED or SCIP_REDUCEDDOM.
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result:
 *    SCIP_BRANCHED   : at least one constraint is infeasible, and branching was applied
 *    SCIP_REDUCEDDOM : at least one constraint is infeasible, and a domain was reduced
 *    SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *    SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define DECL_CONSENPS(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip, CONS** conss, int nconss, RESULT* result)

/** feasibility check method of constraint handler for integral solutions
 *
 *  The given solution has to be checked for feasibility.
 *  
 *  The check methods of the active constraint handlers are called in decreasing order of their check
 *  priorities until the first constraint handler returned with the result SCIP_INFEASIBLE.
 *  The integrality constraint handler has a check priority of zero. A constraint handler which can
 *  (or wants) to check its constraints only for integral solutions should have a negative check priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to check feasibility even on non-integral solutions must have a
 *  check priority greater than zero (e.g. if the check is much faster than testing all variables for
 *  integrality).
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    sol             : the solution to check feasibility for
 *    result          : pointer to store the result of the feasibility checking call
 *
 *  possible return values for *result:
 *    SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
 *    SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define DECL_CONSCHCK(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip, CONS** conss, int nconss, SOL* sol, RESULT* result)

/** domain propagation method of constraint handler
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *    SCIP_REDUCEDDOM : at least one domain reduction was found
 *    SCIP_DIDNOTFIND : the propagator searched and did not find any domain reductions
 *    SCIP_DIDNOTRUN  : the propagator was skipped
 */
#define DECL_CONSPROP(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip, CONS** conss, int nconss, RESULT* result)



#include "scip.h"
#include "retcode.h"
#include "result.h"
#include "mem.h"
#include "lp.h"
#include "sol.h"



/** linked list of constraints */
struct ConsList
{
   CONS*            cons;               /**< pointer to constraint data structure */
   CONSLIST*        next;               /**< next list entry */
};



/*
 * Constraint handler methods
 */

DECL_SORTPTRCOMP(SCIPconshdlrCompSepa); /**< compares two constraint handlers w. r. to their separation priority */
DECL_SORTPTRCOMP(SCIPconshdlrCompEnfo); /**< compares two constraint handlers w. r. to their enforcing priority */
DECL_SORTPTRCOMP(SCIPconshdlrCompChck); /**< compares two constraint handlers w. r. to their feasibility check priority */

extern
RETCODE SCIPconshdlrCreate(             /**< creates a constraint handler */
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority,       /**< priority of the constraint handler for checking infeasibility */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE((*consfree)),          /**< destructor of constraint handler */
   DECL_CONSINIT((*consinit)),          /**< initialise constraint handler */
   DECL_CONSEXIT((*consexit)),          /**< deinitialise constraint handler */
   DECL_CONSDELE((*consdele)),          /**< free specific constraint data */
   DECL_CONSTRAN((*constran)),          /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA((*conssepa)),          /**< separate cutting planes */
   DECL_CONSENLP((*consenlp)),          /**< enforcing constraints for LP solutions */
   DECL_CONSENPS((*consenps)),          /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHCK((*conschck)),          /**< check feasibility of primal solution */
   DECL_CONSPROP((*consprop)),          /**< propagate variable domains */
   CONSHDLRDATA*    conshdlrdata        /**< constraint handler data */
   );

extern
RETCODE SCIPconshdlrFree(               /**< calls destructor and frees memory of constraint handler */
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
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
RETCODE SCIPconshdlrSeparate(           /**< calls separator method of constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
RETCODE SCIPconshdlrEnforceLPSol(       /**< calls enforcing method of constraint handler for LP solutions */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
RETCODE SCIPconshdlrEnforcePseudoSol(   /**< calls enforcing method of constraint handler for pseudo solutions */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
RETCODE SCIPconshdlrCheck(              /**< calls feasibility check method of constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   SOL*             sol,                /**< primal CIP solution */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
RETCODE SCIPconshdlrPropagate(          /**< calls propagation method of constraint handler */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

extern
const char* SCIPconshdlrGetName(        /**< gets name of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

extern
CONSHDLRDATA* SCIPconshdlrGetData(      /**< gets user data of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

extern
void SCIPconshdlrSetData(               /**< sets user data of constraint handler; user has to free old data in advance! */
   CONSHDLR*        conshdlr,           /**< constraint handler */
   CONSHDLRDATA*    conshdlrdata        /**< new constraint handler user data */
   );

extern
CONS** SCIPconshdlrGetConss(            /**< gets constraints array of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

extern
int SCIPconshdlrGetNConss(              /**< gets number of constraints in constraints array of constraint handler */
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

extern
Bool SCIPconshdlrIsInitialized(         /**< is constraint handler initialized? */
   CONSHDLR*        conshdlr            /**< constraint handler */
   );



/*
 * Constraint methods
 */

extern
RETCODE SCIPconsCreate(                 /**< creates and captures a constraint */
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             original,           /**< is constraint belonging to the original problem? */
   Bool             model               /**< is constraint necessary for feasibility? */
   );

extern
RETCODE SCIPconsFree(                   /**< frees a constraint */
   CONS**           cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   );

extern
void SCIPconsCapture(                   /**< increases usage counter of constraint */
   CONS*            cons                /**< constraint */
   );

extern
RETCODE SCIPconsRelease(                /**< decreases usage counter of constraint, and frees memory if necessary */
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPconsActivate(               /**< activates constraint */
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPconsDeactivate(             /**< deactivates constraint */
   CONS*            cons                /**< constraint */
   );

extern
RETCODE SCIPconsTransform(              /**< copies original constraint into transformed constraint, that is captured */
   CONS**           transcons,          /**< pointer to store the transformed constraint */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONS*            origcons            /**< original constraint */
   );

extern
const char* SCIPconsGetName(            /**< returns the name of the constraint */
   CONS*            cons                /**< constraint */
   );

extern
CONSHDLR* SCIPconsGetConsHdlr(          /**< returns the constraint handler of the constraint */
   CONS*            cons                /**< constraint */
   );

extern
CONSDATA* SCIPconsGetConsdata(          /**< returns the constraint data field of the constraint */
   CONS*            cons                /**< constraint */
   );

extern
Bool SCIPconsIsOriginal(                /**< returns TRUE iff constraint is belonging to original problem */
   CONS*            cons                /**< constraint */
   );

extern
Bool SCIPconsIsModel(                   /**< returns TRUE iff constraint is necessary for feasibility */
   CONS*            cons                /**< constraint */
   );


/*
 * Hash functions
 */

extern
DECL_HASHGETKEY(SCIPhashGetKeyCons);    /**< gets the key (i.e. the name) of the given constraint */



/*
 * Constraint list methods
 */

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


#endif
