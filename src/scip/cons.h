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
typedef struct ConsHdlrData CONSHDLRDATA; /**< constraint handler data */
typedef struct ConsData CONSDATA;       /**< locally defined constraint type specific data */
typedef struct ConsSetChg CONSSETCHG;   /**< tracks additions and removals of the set of active constraints */
typedef struct ConsSetChgDyn CONSSETCHGDYN; /**< dynamic size attachment for constraint set change data */



/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 */
#define DECL_CONSFREE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr)

/** initialization method of constraint handler (called when problem solving starts)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 */
#define DECL_CONSINIT(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr)

/** deinitialization method of constraint handler (called when problem solving exits)
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 */
#define DECL_CONSEXIT(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr)

/** frees specific constraint data
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    consdata        : pointer to the constraint data to free
 */
#define DECL_CONSDELETE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONSDATA** consdata)

/** transforms constraint data into data belonging to the transformed problem
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    sourcecons      : source constraint to transform
 *    targetcons      : pointer to store created target constraint
 */
#define DECL_CONSTRANS(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* sourcecons, CONS** targetcons)

/** separation method of constraint handler
 *
 *  Separates all constraints of the constraint handler. The method is called in the LP solution loop,
 *  which means that a valid LP solution exists.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    nusefulconss    : number of useful (non-obsolete) constraints to process
 *    result          : pointer to store the result of the separation call
 *
 *  possible return values for *result:
 *    SCIP_SEPARATED  : at least one cutting plane was generated
 *    SCIP_CONSADDED  : at least one additional constraint was generated
 *    SCIP_DIDNOTFIND : the separator searched, but didn't found a cutting plane
 *    SCIP_DIDNOTRUN  : the separator was skipped
 */
#define DECL_CONSSEPA(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss, int nusefulconss, \
                                    RESULT* result)

/** constraint enforcing method of constraint handler for LP solutions
 *
 *  The method is called at the end of the node processing loop for a node where the LP was solved.
 *  The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
 *  branching, reducing a variable's domain to exclude the solution or separating the solution with a valid
 *  cutting plane.
 *
 *  The enforcing methods of the active constraint handlers are called in decreasing order of their enforcing
 *  priorities until the first constraint handler returned with the value SCIP_BRANCHED, SCIP_REDUCEDDOM,
 *  SCIP_SEPARATED, or SCIP_CONSADDED.
 *  The integrality constraint handler has an enforcing priority of zero. A constraint handler which can
 *  (or wants) to enforce its constraints only for integral solutions should have a negative enforcing priority
 *  (e.g. the alldiff-constraint can only operate on integral solutions).
 *  A constraint handler which wants to incorporate its own branching strategy even on non-integral
 *  solutions must have an enforcing priority greater than zero (e.g. the SOS-constraint incorporates
 *  SOS-branching on non-integral solutions).
 *  If the solution is integral and one of the constraints of the constraint handler is violated, the
 *  constraint handler has to branch, to reduce a variable's domain, to create a cutting plane, or to add an
 *  additional constraint that cuts off the solution -- otherwise, the infeasibility cannot be resolved.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
 *  method should process the useful constraints first, and only if no violation was found, the remaining
 *  nconss - nusefulconss constraints must be enforced.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    nusefulconss    : number of useful (non-obsolete) constraints to process
 *    result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : at least one constraint is infeasible, and it cannot be resolved -> node is infeasible
 *    SCIP_BRANCHED   : at least one constraint is infeasible, and branching was applied to resolve infeasibility
 *    SCIP_REDUCEDDOM : at least one constraint is infeasible, and a domain was reduced to resolve infeasibility
 *    SCIP_SEPARATED  : at least one constraint is infeasible, and a cutting plane was generated to resolve infeasibility
 *    SCIP_CONSADDED  : at least one constraint is infeasible, and a constraint was generated to resolve infeasibility
 *    SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *    SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define DECL_CONSENFOLP(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss, int nusefulconss, \
                                      RESULT* result)

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
 *  the value SCIP_BRANCHED, SCIP_REDUCEDDOM, or SCIP_CONSADDED.
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
 *  method should process the useful constraints first, and only if no violation was found, the remaining
 *  nconss - nusefulconss constraints must be enforced.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    nusefulconss    : number of useful (non-obsolete) constraints to process
 *    result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : at least one constraint is infeasible, and it cannot be resolved -> node is infeasible
 *    SCIP_BRANCHED   : at least one constraint is infeasible, and branching was applied to resolve infeasibility
 *    SCIP_REDUCEDDOM : at least one constraint is infeasible, and a domain was reduced to resolve infeasibility
 *    SCIP_CONSADDED  : at least one constraint is infeasible, and a constraint was generated to resolve infeasibility
 *    SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *    SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define DECL_CONSENFOPS(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss, int nusefulconss, \
                                      RESULT* result)

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
 *  feasibility is already checked or implicitly given. In these cases, 'checkintegrality' or
 *  'checklprows' is FALSE.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    sol             : the solution to check feasibility for
 *    checkintegrality: has integrality to be checked?
 *    checklprows     : have current LP rows to be checked?
 *    result          : pointer to store the result of the feasibility checking call
 *
 *  possible return values for *result:
 *    SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
 *    SCIP_FEASIBLE   : all constraints of the handler are feasible
 */
#define DECL_CONSCHECK(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss, SOL* sol, \
                                     Bool checkintegrality, Bool checklprows, RESULT* result)

/** domain propagation method of constraint handler
 *
 *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The propagation
 *  method should process only the useful constraints in most runs, and only occasionally the remaining
 *  nconss - nusefulconss constraints.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    conss           : array of constraints to process
 *    nconss          : number of constraints to process
 *    nusefulconss    : number of useful (non-obsolete) constraints to process
 *    result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *    SCIP_CUTOFF     : at least one constraint is infeasible for the actual domains -> node is infeasible
 *    SCIP_REDUCEDDOM : at least one domain reduction was found
 *    SCIP_DIDNOTFIND : the propagator searched and did not find any domain reductions
 *    SCIP_DIDNOTRUN  : the propagator was skipped
 */
#define DECL_CONSPROP(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss, int nusefulconss, \
                                    RESULT* result)

/** constraint enabling notification method of constraint handler
 *
 *  This method is always called after a constraint of the constraint handler was enabled. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    cons            : the constraint that has been activated
 */
#define DECL_CONSENABLE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* cons)

/** constraint disabling notification method of constraint handler
 *
 *  This method is always called before a constraint of the constraint handler is disabled. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *    scip            : SCIP main data structure
 *    conshdlr        : the constraint handler itself
 *    cons            : the constraint that will be deactivated
 */
#define DECL_CONSDISABLE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* cons)



#include "scip.h"
#include "retcode.h"
#include "result.h"
#include "mem.h"
#include "lp.h"
#include "sol.h"
#include "tree.h"



/** constraint data structure */
struct Cons
{
   char*            name;               /**< name of the constraint */
   CONSHDLR*        conshdlr;           /**< constraint handler for this constraint */
   CONSDATA*        consdata;           /**< data for this specific constraint */
   NODE*            node;               /**< node where this constraint was created, or NULL if it's a global constraint */
   int              nuses;              /**< number of times, this constraint is referenced */
   int              age;                /**< age of constraint: number of successive times, the constraint was irrelevant */
   int              sepaconsspos;       /**< position of constraint in the handler's sepaconss array */
   int              enfoconsspos;       /**< position of constraint in the handler's enfoconss array */
   int              checkconsspos;      /**< position of constraint in the handler's checkconss array */
   int              propconsspos;       /**< position of constraint in the handler's propconss array */
   int              arraypos;           /**< position of constraint in the node's/problem's addedconss/conss array */
   unsigned int     separate:1;         /**< TRUE iff constraint should be separated during LP processing */
   unsigned int     enforce:1;          /**< TRUE iff constraint should be enforced during node processing */
   unsigned int     check:1;            /**< TRUE iff constraint should be checked for feasibility */
   unsigned int     propagate:1;        /**< TRUE iff constraint should be propagated during node processing */
   unsigned int     original:1;         /**< TRUE iff constraint belongs to original problem */
   unsigned int     active:1;           /**< TRUE iff constraint is active in the active node */
   unsigned int     enabled:1;          /**< TRUE iff constraint is enforced, separated, and propagated in active node */
   unsigned int     obsolete:1;         /**< TRUE iff constraint is too seldomly used and therefore obsolete */
   unsigned int     update:1;           /**< TRUE iff constraint has to be updated in update phase */
   unsigned int     updateactivate:1;   /**< TRUE iff constraint has to be activated in update phase */
   unsigned int     updatedeactivate:1; /**< TRUE iff constraint has to be deactivated in update phase */
   unsigned int     updateenable:1;     /**< TRUE iff constraint has to be enabled in update phase */
   unsigned int     updatedisable:1;    /**< TRUE iff constraint has to be disabled in update phase */
   unsigned int     updateobsolete:1;   /**< TRUE iff obsolete status of constraint has to be updated in update phase */
};

/** tracks additions and removals of the set of active constraints */
struct ConsSetChg
{
   CONS**           addedconss;         /**< constraints added to the set of active constraints */
   CONS**           disabledconss;      /**< constraints disabled in the set of active constraints */
   int              naddedconss;        /**< number of added constraints */
   int              ndisabledconss;     /**< number of disabled constraints */
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
DECL_SORTPTRCOMP(SCIPconshdlrCompCheck);

/** creates a constraint handler */
extern
RETCODE SCIPconshdlrCreate(
   CONSHDLR**       conshdlr,           /**< pointer to constraint handler data structure */
   const char*      name,               /**< name of constraint handler */
   const char*      desc,               /**< description of constraint handler */
   int              sepapriority,       /**< priority of the constraint handler for separation */
   int              enfopriority,       /**< priority of the constraint handler for constraint enforcing */
   int              checkpriority,      /**< priority of the constraint handler for checking infeasibility */
   int              sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
   int              propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
   Bool             needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
   DECL_CONSFREE    ((*consfree)),      /**< destructor of constraint handler */
   DECL_CONSINIT    ((*consinit)),      /**< initialise constraint handler */
   DECL_CONSEXIT    ((*consexit)),      /**< deinitialise constraint handler */
   DECL_CONSDELETE  ((*consdelete)),    /**< free specific constraint data */
   DECL_CONSTRANS   ((*constrans)),     /**< transform constraint data into data belonging to the transformed problem */
   DECL_CONSSEPA    ((*conssepa)),      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   DECL_CONSENABLE  ((*consenable)),    /**< enabling notification method */
   DECL_CONSDISABLE ((*consdisable)),   /**< disabling notification method */
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

/** calls separator method of constraint handler to separate all constraints added after last conshdlrReset() call */
extern
RETCODE SCIPconshdlrSeparate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   int              actdepth,           /**< depth of active node */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls enforcing method of constraint handler for LP solution for all constraints added after last
 *  conshdlrReset() call
 */
extern
RETCODE SCIPconshdlrEnforceLPSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls enforcing method of constraint handler for pseudo solution for all constraints added after last
 *  conshdlrReset() call
 */
extern
RETCODE SCIPconshdlrEnforcePseudoSol(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls feasibility check method of constraint handler */
extern
RETCODE SCIPconshdlrCheck(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** calls propagation method of constraint handler */
extern
RETCODE SCIPconshdlrPropagate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   int              actdepth,           /**< depth of active node; -1 if preprocessing domain propagation */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** resets separation to start with first constraint in the next call */
extern
void SCIPconshdlrResetSepa(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** resets enforcement to start with first constraint in the next call */
extern
void SCIPconshdlrResetEnfo(
   CONSHDLR*        conshdlr            /**< constraint handler */
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

/** gets number of active constraints of constraint handler */
extern
int SCIPconshdlrGetNActiveConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of enabled constraints of constraint handler */
extern
int SCIPconshdlrGetNEnabledConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets checking priority of constraint handler */
extern
int SCIPconshdlrGetCheckPriority(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets separation frequency of constraint handler */
extern
int SCIPconshdlrGetSepaFreq(
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

/** creates and captures a constraint
 *  Warning! If a constraint is marked to be checked for feasibility but not to be enforced, a LP or pseudo solution
 *  may be declared feasible even if it violates this particular constraint.
 *  This constellation should only be used, if no LP or pseudo solution can violate the constraint -- e.g. if a
 *  local constraint is redundant due to the variable's local bounds.
 */
extern
RETCODE SCIPconsCreate(
   CONS**           cons,               /**< pointer to constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const char*      name,               /**< name of constraint */
   CONSHDLR*        conshdlr,           /**< constraint handler for this constraint */
   CONSDATA*        consdata,           /**< data for this specific constraint */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             original            /**< is constraint belonging to the original problem? */
   );

/** frees constraint data of a constraint, leaving the constraint itself as a zombie constraint */
extern
RETCODE SCIPconsFreeData(
   CONS*            cons,               /**< constraint to free */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set                 /**< global SCIP settings */
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

/** globally removes constraint from all subproblems; removes constraint from the addedconss array of the node, where it
 *  was created, or from the problem, if it was a problem constraint;
 *  the method must not be called for local check-constraint (i.e. constraints, that locally ensure feasibility);
 *  the constraint data is freed, and if the constraint is no longer used, it is freed completely
 */
extern
RETCODE SCIPconsDelete(
   CONS*            cons,               /**< constraint to delete */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   );

/** copies original constraint into transformed constraint, that is captured */
extern
RETCODE SCIPconsTransform(
   CONS**           transcons,          /**< pointer to store the transformed constraint */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   CONS*            origcons            /**< original constraint */
   );

/** activates constraint or marks constraint to be activated in next update */
extern
RETCODE SCIPconsActivate(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   );

/** deactivates constraint or marks constraint to be deactivated in next update */
extern
RETCODE SCIPconsDeactivate(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   );

/** enables constraint's separation, enforcing, and propagation capabilities or marks them to be enabled in next update */
extern
RETCODE SCIPconsEnable(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   );

/** disables constraint's separation, enforcing, and propagation capabilities or marks them to be disabled in next update */
extern
RETCODE SCIPconsDisable(
   CONS*            cons,               /**< constraint */
   const SET*       set                 /**< global SCIP settings */
   );

/** increases age of constraint; should be called in constraint separation, if no cut was found for this constraint,
 *  in constraint enforcing, if constraint was feasible, and in constraint propagation, if no domain reduction was
 *  deduced;
 *  if it's age exceeds the constraint age limit, makes constraint obsolete or marks constraint to be made obsolete
 *  in next update,
 */
extern
RETCODE SCIPconsIncAge(
   CONS*            cons,               /**< constraint */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   );

/** resets age of constraint to zero; should be called in constraint separation, if a cut was found for this constraint,
 *  in constraint enforcing, if the constraint was violated, and in constraint propagation, if a domain reduction was
 *  deduced;
 *  if it was obsolete, makes constraint useful again or marks constraint to be made useful again in next update
 */
extern
RETCODE SCIPconsResetAge(
   CONS*            cons,               /**< constraint */
   SET*             set                 /**< global SCIP settings */
   );

/** returns the name of the constraint */
extern
const char* SCIPconsGetName(
   CONS*            cons                /**< constraint */
   );

/** returns the constraint handler of the constraint */
extern
CONSHDLR* SCIPconsGetHdlr(
   CONS*            cons                /**< constraint */
   );

/** returns the constraint data field of the constraint */
extern
CONSDATA* SCIPconsGetData(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be separated during LP processing */
extern
Bool SCIPconsIsSeparated(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be enforced during node processing */
extern
Bool SCIPconsIsEnforced(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be checked for feasibility */
extern
Bool SCIPconsIsChecked(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be propagated during node processing */
extern
Bool SCIPconsIsPropagated(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is belonging to original problem */
extern
Bool SCIPconsIsOriginal(
   CONS*            cons                /**< constraint */
   );




/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given constraint */
extern
DECL_HASHGETKEY(SCIPhashGetKeyCons);




/*
 * Constraint set change methods
 */

/** frees fixed size constraint set change data and releases all included constraints */
extern
RETCODE SCIPconssetchgFree(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** deactivates, deletes, and releases constraint from the addedconss array of the constraint set change data */
extern
RETCODE SCIPconssetchgDelAddedCons(
   CONSSETCHG*      conssetchg,         /**< constraint set change to delete constraint from */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to delete from addedconss array */
   );

/** applies constraint set change */
extern
RETCODE SCIPconssetchgApply(
   CONSSETCHG*      conssetchg,         /**< constraint set change to apply */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** undoes constraint set change */
extern
RETCODE SCIPconssetchgUndo(
   CONSSETCHG*      conssetchg,         /**< constraint set change to undo */
   const SET*       set                 /**< global SCIP settings */
   );




/*
 * dynamic size attachment methods for constraint set changes
 */

/** creates a dynamic size attachment for a constraint set change data structure */
extern
RETCODE SCIPconssetchgdynCreate(
   CONSSETCHGDYN**  conssetchgdyn,      /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   );

/** frees a dynamic size attachment for a constraint set change data structure */
extern
void SCIPconssetchgdynFree(
   CONSSETCHGDYN**  conssetchgdyn,      /**< pointer to dynamic size attachment */
   MEMHDR*          memhdr              /**< block memory */
   );

/** attaches dynamic size information to constraint set change data */
extern
void SCIPconssetchgdynAttach(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamic size information */
   CONSSETCHG**     conssetchg          /**< pointer to static constraint set change */
   );

/** detaches dynamic size information and shrinks constraint set change data */
extern
RETCODE SCIPconssetchgdynDetach(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamic size information */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** frees attached constraint set change data and detaches dynamic size attachment */
extern
RETCODE SCIPconssetchgdynDiscard(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamically sized constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** adds constraint addition to constraint set changes, and captures constraint */
extern
RETCODE SCIPconssetchgdynAddAddedCons(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamically sized constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   NODE*            node,               /**< node that the constraint set change belongs to */
   CONS*            cons                /**< added constraint */
   );

/** adds constraint disabling to constraint set changes, and captures constraint */
extern
RETCODE SCIPconssetchgdynAddDisabledCons(
   CONSSETCHGDYN*   conssetchgdyn,      /**< dynamically sized constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< disabled constraint */
   );

/** gets pointer to constraint set change data the dynamic size information references */
extern
CONSSETCHG** SCIPconssetchgdynGetConssetchgPtr(
   CONSSETCHGDYN*   conssetchgdyn       /**< dynamically sized constraint set change data structure */
   );


#endif
