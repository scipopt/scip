/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
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
typedef struct ConsHdlrData CONSHDLRDATA; /**< constraint handler data */
typedef struct ConsData CONSDATA;       /**< locally defined constraint type specific data */


/** destructor of constraint handler to free user data (called when SCIP is exiting)
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 */
#define DECL_CONSFREE(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip)

/** initialization method of constraint handler (called when problem solving starts)
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 */
#define DECL_CONSINIT(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip)

/** deinitialization method of constraint handler (called when problem solving exits)
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
 *    SCIP_CUTOFF     : at least one constraint is infeasible, and it cannot be resolved -> node is infeasible
 *    SCIP_BRANCHED   : at least one constraint is infeasible, and branching was applied to resolve infeasibility
 *    SCIP_REDUCEDDOM : at least one constraint is infeasible, and a domain was reduced to resolve infeasibility
 *    SCIP_SEPARATED  : at least one constraint is infeasible, and a cutting plane was generated to resolve infeasibility
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
 *    SCIP_CUTOFF     : at least one constraint is infeasible, and it cannot be resolved -> node is infeasible
 *    SCIP_BRANCHED   : at least one constraint is infeasible, and branching was applied to resolve infeasibility
 *    SCIP_REDUCEDDOM : at least one constraint is infeasible, and a domain was reduced to resolve infeasibility
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
 *  In some cases, integrality conditions or rows in actual LP don't have to be checked, because their
 *  feasibility is already checked or implicitly given. In these cases, 'chckintegrality' or
 *  'chcklprows' is FALSE.
 *
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    sol             : the solution to check feasibility for
 *    chckintegrality : has integrality to be checked?
 *    chcklprows      : have current LP rows to be checked?
 *    result          : pointer to store the result of the feasibility checking call
 *
 *  possible return values for *result:
 *    SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
 *    SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define DECL_CONSCHCK(x) RETCODE x (CONSHDLR* conshdlr, SCIP* scip, CONS** conss, int nconss, SOL* sol, \
                                    Bool chckintegrality, Bool chcklprows, RESULT* result)

/** domain propagation method of constraint handler
 *  input:
 *    conshdlr        : the constraint handler itself
 *    scip            : SCIP main data structure
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : at least one constraint is infeasible for the actual domains -> node is infeasible
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

/** compares two constraint handlers w. r. to their separation priority */
extern
DECL_SORTPTRCOMP(SCIPconshdlrCompSepa);

/** compares two constraint handlers w. r. to their enforcing priority */
extern
DECL_SORTPTRCOMP(SCIPconshdlrCompEnfo);

/** compares two constraint handlers w. r. to their feasibility check priority */
extern
DECL_SORTPTRCOMP(SCIPconshdlrCompChck);

/** creates a constraint handler */
extern
RETCODE SCIPconshdlrCreate(
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              chckpriority,       /**< priority of the constraint handler for checking infeasibility */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
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

/** calls destructor and frees memory of constraint handler */
extern
RETCODE SCIPconshdlrFree(
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes constraint handler */
extern
RETCODE SCIPconshdlrInit(
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls exit method of constraint handler */
extern
RETCODE SCIPconshdlrExit(
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls separator method of constraint handler */
extern
RETCODE SCIPconshdlrSeparate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls enforcing method of constraint handler for LP solutions */
extern
RETCODE SCIPconshdlrEnforceLPSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls enforcing method of constraint handler for pseudo solutions */
extern
RETCODE SCIPconshdlrEnforcePseudoSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls feasibility check method of constraint handler */
extern
RETCODE SCIPconshdlrCheck(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   SOL*             sol,                /**< primal CIP solution */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls propagation method of constraint handler */
extern
RETCODE SCIPconshdlrPropagate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   int              actdepth,           /**< depth of active node; -1 if preprocessing domain propagation */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** adds constraint to constraint handler's problem constraint array */
extern
RETCODE SCIPconshdlrAddProbCons(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< model constraint of initial problem to add */
   );

/** gets name of constraint handler */
extern
const char* SCIPconshdlrGetName(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets user data of constraint handler */
extern
CONSHDLRDATA* SCIPconshdlrGetData(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** sets user data of constraint handler; user has to free old data in advance! */
extern
void SCIPconshdlrSetData(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   CONSHDLRDATA*    conshdlrdata        /**< new constraint handler user data */
   );

/** gets constraints array of constraint handler */
extern
CONS** SCIPconshdlrGetConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints in constraints array of constraint handler */
extern
int SCIPconshdlrGetNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets checking priority of constraint handler */
extern
int SCIPconshdlrGetChckPriority(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets propagation frequency of constraint handler */
extern
int SCIPconshdlrGetPropFreq(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** is constraint handler initialized? */
extern
Bool SCIPconshdlrIsInitialized(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );



/*
 * Constraint methods
 */

/** creates and captures a constraint */
extern
RETCODE SCIPconsCreate(
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             original,           /**< is constraint belonging to the original problem? */
   Bool             model               /**< is constraint necessary for feasibility? */
   );

/** frees a constraint */
extern
RETCODE SCIPconsFree(
   CONS**           cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   );

/** increases usage counter of constraint */
extern
void SCIPconsCapture(
   CONS*            cons                /**< constraint */
   );

/** decreases usage counter of constraint, and frees memory if necessary */
extern
RETCODE SCIPconsRelease(
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** activates constraint */
extern
RETCODE SCIPconsActivate(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   );

/** deactivates constraint */
extern
RETCODE SCIPconsDeactivate(
   CONS*            cons                /**< constraint */
   );

/** copies original constraint into transformed constraint, that is captured */
extern
RETCODE SCIPconsTransform(
   CONS**           transcons,          /**< pointer to store the transformed constraint */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONS*            origcons            /**< original constraint */
   );

/** returns the name of the constraint */
extern
const char* SCIPconsGetName(
   CONS*            cons                /**< constraint */
   );

/** returns the constraint handler of the constraint */
extern
CONSHDLR* SCIPconsGetConsHdlr(
   CONS*            cons                /**< constraint */
   );

/** returns the constraint data field of the constraint */
extern
CONSDATA* SCIPconsGetConsData(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is belonging to original problem */
extern
Bool SCIPconsIsOriginal(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is necessary for feasibility */
extern
Bool SCIPconsIsModel(
   CONS*            cons                /**< constraint */
   );


/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given constraint */
extern
DECL_HASHGETKEY(SCIPhashGetKeyCons);



/*
 * Constraint list methods
 */

/** adds constraint to a list of constraints and captures it */
extern
RETCODE SCIPconslistAdd(
   CONSLIST**       conslist,           /**< constraint list to extend */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   );

/** partially unlinks and frees the constraints in the list */
extern
RETCODE SCIPconslistFreePart(
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONSLIST*        firstkeep           /**< first constraint list entry to keep */
   );

/** unlinks and frees all the constraints in the list */
extern
RETCODE SCIPconslistFree(
   CONSLIST**       conslist,           /**< constraint list to delete from */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
   );


#endif
