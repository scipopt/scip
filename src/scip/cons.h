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



/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 */
#define DECL_CONSFREE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr)

/** initialization method of constraint handler (called when problem solving starts)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 */
#define DECL_CONSINIT(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr)

/** deinitialization method of constraint handler (called when problem solving exits)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 */
#define DECL_CONSEXIT(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr)

/** frees specific constraint data
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
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
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - sourcecons      : source constraint to transform
 *  - targetcons      : pointer to store created target constraint
 */ 
#define DECL_CONSTRANS(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* sourcecons, CONS** targetcons)

/** LP initialization method of constraint handler
 *
 *  Puts the LP relaxations of all "initial" constraints into the LP. The method should scan the constraints
 *  array for constraints that are marked initial via calls to SCIPconsIsInitial() and put the LP relaxation
 *  of all initial constraints to the LP with calls to SCIPaddCut().
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 */
#define DECL_CONSINITLP(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss)

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
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - result          : pointer to store the result of the separation call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> node is infeasible
 *  - SCIP_SEPARATED  : at least one cutting plane was generated
 *  - SCIP_REDUCEDDOM : no cutting plane was generated, but at least one domain was reduced
 *  - SCIP_CONSADDED  : no cutting plane or domain reductions, but at least one additional constraint was generated
 *  - SCIP_DIDNOTFIND : the separator searched, but didn't found a cutting plane
 *  - SCIP_DIDNOTRUN  : the separator was skipped
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
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : at least one constraint is infeasible, and it cannot be resolved -> node is infeasible
 *  - SCIP_SEPARATED  : a cutting plane was generated to resolve an infeasibility
 *  - SCIP_REDUCEDDOM : no cutting plane was generated, but at least one domain was reduced to resolve an infeasibility
 *  - SCIP_CONSADDED  : no cutting plane or domain reductions, but a constraint was generated to resolve an infeasibility
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
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
 *  If the pseudo solution's objective value is lower than the lower bound of the node, it cannot be feasible
 *  and the enforcing method may skip it's check and set *result to SCIP_DIDNOTRUN. However, it can also process
 *  it's constraints and return any other possible result code.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - objinfeasible   : is the solution infeasible anyway due to violating lower objective bound?
 *  - result          : pointer to store the result of the enforcing call
 *
 *  possible return values for *result:
 *  - SCIP_DIDNOTRUN  : the enforcement was skipped (only possible, if objinfeasible is true)
 *  - SCIP_CUTOFF     : at least one constraint is infeasible, and it cannot be resolved -> node is infeasible
 *  - SCIP_REDUCEDDOM : at least one domain was reduced to resolve an infeasibility
 *  - SCIP_CONSADDED  : no domain reductions, but a constraint was generated to resolve an infeasibility
 *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
 *  - SCIP_SOLVELP    : at least one constraint is infeasible, and this can only be resolved by solving the LP
 *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
 */ 
#define DECL_CONSENFOPS(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss, int nusefulconss, \
                                      Bool objinfeasible, RESULT* result)

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
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - sol             : the solution to check feasibility for
 *  - checkintegrality: has integrality to be checked?
 *  - checklprows     : have current LP rows to be checked?
 *  - result          : pointer to store the result of the feasibility checking call
 *
 *  possible return values for *result:
 *  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
 *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
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
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nusefulconss    : number of useful (non-obsolete) constraints to process
 *  - result          : pointer to store the result of the propagation call
 *
 *  possible return values for *result:
 *  - SCIP_CUTOFF     : at least one constraint is infeasible for the actual domains -> node is infeasible
 *  - SCIP_REDUCEDDOM : at least one domain reduction was found
 *  - SCIP_DIDNOTFIND : the propagator searched and did not find any domain reductions
 *  - SCIP_DIDNOTRUN  : the propagator was skipped
 */
#define DECL_CONSPROP(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss, int nusefulconss, \
                                    RESULT* result)

/** presolving method of constraint handler
 *
 *  The presolver should go through the variables and constraints and tighten the domains or
 *  constraints. Each tightening should increase the given total number of changes.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - conss           : array of constraints to process
 *  - nconss          : number of constraints to process
 *  - nrounds         : number of presolving rounds already done
 *  - nnewfixedvars   : number of variables fixed since the last call to the presolving method
 *  - nnewaggrvars    : number of variables aggregated since the last call to the presolving method
 *  - nnewchgvartypes : number of variable type changes since the last call to the presolving method
 *  - nnewchgbds      : number of variable bounds tightend since the last call to the presolving method
 *  - nnewholes       : number of domain holes added since the last call to the presolving method
 *  - nnewdelconss    : number of deleted constraints since the last call to the presolving method
 *  - nnewupgdconss   : number of upgraded constraints since the last call to the presolving method
 *  - nnewchgcoefs    : number of changed coefficients since the last call to the presolving method
 *  - nnewchgsides    : number of changed left or right hand sides since the last call to the presolving method
 *
 *  input/output:
 *  - nfixedvars      : pointer to total number of variables fixed of all presolvers
 *  - naggrvars       : pointer to total number of variables aggregated of all presolvers
 *  - nchgvartypes    : pointer to total number of variable type changes of all presolvers
 *  - nchgbds         : pointer to total number of variable bounds tightend of all presolvers
 *  - naddholes       : pointer to total number of domain holes added of all presolvers
 *  - ndelconss       : pointer to total number of deleted constraints of all presolvers
 *  - nupgdconss      : pointer to total number of upgraded constraints of all presolvers
 *  - nchgcoefs       : pointer to total number of changed coefficients of all presolvers
 *  - nchgsides       : pointer to total number of changed left/right hand sides of all presolvers
 *
 *  output:
 *  - result          : pointer to store the result of the presolving call
 *
 *  possible return values for *result:
 *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
 *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
 *  - SCIP_SUCCESS    : the presolver found a reduction
 *  - SCIP_DIDNOTFIND : the presolver searched, but didn't found a presolving change
 *  - SCIP_DIDNOTRUN  : the presolver was skipped
 */
#define DECL_CONSPRESOL(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS** conss, int nconss, int nrounds, \
   int nnewfixedvars, int nnewaggrvars, int nnewchgvartypes, int nnewchgbds, int nnewholes, \
   int nnewdelconss, int nnewupgdconss, int nnewchgcoefs, int nnewchgsides,                 \
   int* nfixedvars, int* naggrvars, int* nchgvartypes, int* nchgbds, int* naddholes,        \
   int* ndelconss, int* nupgdconss, int* nchgcoefs, int* nchgsides, RESULT* result)


/** conflict variable resolving method of constraint handler
 *
 *  This method is called during conflict analysis. If the conflict handler wants to support conflict analysis,
 *  it should call SCIPinferBinVar() in domain propagation in order to fix binary variables to deduced values.
 *  In this call, the handler provides the constraint, that deduced the variable's assignment. The conflict
 *  variable resolving method must then be implemented, to provide the "reasons" for the variable assignment,
 *  i.e. the fixed binary variables, that forced the constraint to set the conflict variable to its current
 *  value. The variables that form the reason of the assignment must be provided by calls to SCIPaddConflictVar().
 *
 *  For example, the logicor constraint c = "x or y or z" fixes variable z to TRUE, if both, x and y, are assigned
 *  to FALSE. It uses SCIPinferBinVar(scip, z, 1.0, c) to apply this assignment. In the conflict analysis, the
 *  constraint handler may be asked to resolve variable z. With a call to SCIPvarGetLbLocal(z), the handler can find
 *  out, that variable z is currently assigned to TRUE, and should call SCIPaddConflictVar(scip, x) and
 *  SCIPaddConflictVar(scip, y) to tell SCIP, that the assignments to x and y are the reason for the deduction of z.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that deduced the assignment of the conflict variable
 *  - infervar        : the binary conflict variable that has to be resolved
 */
#define DECL_CONSRESCVAR(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* cons, VAR* infervar)

/** constraint activation notification method of constraint handler
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 *
 *  This method is always called after a constraint of the constraint handler was activated. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that has been activated
 */
#define DECL_CONSACTIVE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* cons)

/** constraint deactivation notification method of constraint handler
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 *
 *  This method is always called before a constraint of the constraint handler is deactivated. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that will be deactivated
 */
#define DECL_CONSDEACTIVE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* cons)

/** constraint enabling notification method of constraint handler
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 *
 *  This method is always called after a constraint of the constraint handler was enabled. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that has been enabled
 */
#define DECL_CONSENABLE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* cons)

/** constraint disabling notification method of constraint handler
 *
 *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
 *  the corresponding bound change event was not yet processed.
 *
 *  This method is always called before a constraint of the constraint handler is disabled. The constraint
 *  handler may use this call to update his own (statistical) data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - conshdlr        : the constraint handler itself
 *  - cons            : the constraint that will be disabled
 */
#define DECL_CONSDISABLE(x) RETCODE x (SCIP* scip, CONSHDLR* conshdlr, CONS* cons)



#include "scip.h"
#include "retcode.h"
#include "result.h"
#include "mem.h"
#include "lp.h"
#include "sol.h"
#include "tree.h"
#include "sepastore.h"



/** constraint data structure */
struct Cons
{
   char*            name;               /**< name of the constraint */
   CONSHDLR*        conshdlr;           /**< constraint handler for this constraint */
   CONSDATA*        consdata;           /**< data for this specific constraint */
   CONSSETCHG*      addconssetchg;      /**< constraint change that added constraint to current subproblem, or NULL if
                                         *   constraint is from global problem */
   int              addarraypos;        /**< position of constraint in the conssetchg's/prob's addedconss/conss array */
   int              consspos;           /**< position of constraint in the handler's conss array */
   int              sepaconsspos;       /**< position of constraint in the handler's sepaconss array */
   int              enfoconsspos;       /**< position of constraint in the handler's enfoconss array */
   int              checkconsspos;      /**< position of constraint in the handler's checkconss array */
   int              propconsspos;       /**< position of constraint in the handler's propconss array */
   int              nuses;              /**< number of times, this constraint is referenced */
   int              age;                /**< age of constraint: number of successive times, the constraint was irrelevant */
   unsigned int     initial:1;          /**< TRUE iff LP relaxation of constraint should be in initial LP, if possible */
   unsigned int     separate:1;         /**< TRUE iff constraint should be separated during LP processing */
   unsigned int     enforce:1;          /**< TRUE iff constraint should be enforced during node processing */
   unsigned int     check:1;            /**< TRUE iff constraint should be checked for feasibility */
   unsigned int     propagate:1;        /**< TRUE iff constraint should be propagated during node processing */
   unsigned int     local:1;            /**< TRUE iff constraint is only valid locally */
   unsigned int     modifiable:1;       /**< TRUE iff constraint is modifiable (subject to column generation) */
   unsigned int     removeable:1;       /**< TRUE iff constraint should be removed from the LP due to aging or cleanup */
   unsigned int     original:1;         /**< TRUE iff constraint belongs to original problem */
   unsigned int     active:1;           /**< TRUE iff constraint is active in the active node */
   unsigned int     enabled:1;          /**< TRUE iff constraint is enforced, separated, and propagated in active node */
   unsigned int     obsolete:1;         /**< TRUE iff constraint is too seldomly used and therefore obsolete */
   unsigned int     deleted:1;          /**< TRUE iff constraint was globally deleted */
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
   int              addedconsssize;     /**< size of added constraints array */
   int              naddedconss;        /**< number of added constraints */
   int              disabledconsssize;  /**< size of disabled constraints array */
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
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
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
   DECL_CONSINITLP  ((*consinitlp)),    /**< initialize LP with relaxations of "initial" constraints */
   DECL_CONSSEPA    ((*conssepa)),      /**< separate cutting planes */
   DECL_CONSENFOLP  ((*consenfolp)),    /**< enforcing constraints for LP solutions */
   DECL_CONSENFOPS  ((*consenfops)),    /**< enforcing constraints for pseudo solutions */
   DECL_CONSCHECK   ((*conscheck)),     /**< check feasibility of primal solution */
   DECL_CONSPROP    ((*consprop)),      /**< propagate variable domains */
   DECL_CONSPRESOL  ((*conspresol)),    /**< presolving method */
   DECL_CONSRESCVAR ((*consrescvar)),   /**< conflict variable resolving method */
   DECL_CONSACTIVE  ((*consactive)),    /**< activation notification method */
   DECL_CONSDEACTIVE((*consdeactive)),  /**< deactivation notification method */
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

/** calls LP initialization method of constraint handler to separate all initial constraints */
extern
RETCODE SCIPconshdlrInitLP(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   SEPASTORE*       sepastore           /**< separation storage */
   );

/** calls separator method of constraint handler to separate all constraints added after last conshdlrReset() call */
extern
RETCODE SCIPconshdlrSeparate(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   SEPASTORE*       sepastore,          /**< separation storage */
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
   SEPASTORE*       sepastore,          /**< separation storage */
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
   Bool             objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
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

/** calls presolving method of constraint handler */
extern
RETCODE SCIPconshdlrPresolve(
   CONSHDLR*        conshdlr,           /**< constraint handler */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   int              nrounds,            /**< number of presolving rounds already done */
   int*             nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*             naggrvars,          /**< pointer to total number of variables aggregated of all presolvers */
   int*             nchgvartypes,       /**< pointer to total number of variable type changes of all presolvers */
   int*             nchgbds,            /**< pointer to total number of variable bounds tightend of all presolvers */
   int*             naddholes,          /**< pointer to total number of domain holes added of all presolvers */
   int*             ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   int*             nupgdconss,         /**< pointer to total number of upgraded constraints of all presolvers */
   int*             nchgcoefs,          /**< pointer to total number of changed coefficients of all presolvers */
   int*             nchgsides,          /**< pointer to total number of changed left/right hand sides of all presolvers */
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
int SCIPconshdlrGetNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of enabled constraints of constraint handler */
extern
int SCIPconshdlrGetNEnabledConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for presolving in this constraint handler */
extern
Real SCIPconshdlrGetPresolTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for separation in this constraint handler */
extern
Real SCIPconshdlrGetSepaTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for LP enforcement in this constraint handler */
extern
Real SCIPconshdlrGetEnfoLPTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for pseudo enforcement in this constraint handler */
extern
Real SCIPconshdlrGetEnfoPSTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets time in seconds used for propagation in this constraint handler */
extern
Real SCIPconshdlrGetPropTime(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's separation method */
extern
Longint SCIPconshdlrGetNSepaCalls(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's LP enforcing method */
extern
Longint SCIPconshdlrGetNEnfoLPCalls(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's pseudo enforcing method */
extern
Longint SCIPconshdlrGetNEnfoPSCalls(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of calls to the constraint handler's propagation method */
extern
Longint SCIPconshdlrGetNPropCalls(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets total number of cuts found by this constraint handler */
extern
Longint SCIPconshdlrGetNCutsFound(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of branchings performed by this constraint handler */
extern
Longint SCIPconshdlrGetNBranchings(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets maximum number of active constraints of constraint handler existing at the same time */
extern
int SCIPconshdlrGetMaxNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** resets maximum number of active constraints to current number of active constraints */
extern
void SCIPconshdlrResetNMaxNConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variables fixed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNFixedVars(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variables aggregated in presolving method of constraint handler */
extern
int SCIPconshdlrGetNAggrVars(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of variable types changed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNVarTypes(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of bounds changed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNChgBds(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of holes added to domains of variables in presolving method of constraint handler */
extern
int SCIPconshdlrGetNAddHoles(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints deleted in presolving method of constraint handler */
extern
int SCIPconshdlrGetNDelConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraints upgraded in presolving method of constraint handler */
extern
int SCIPconshdlrGetNUpgdConss(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of coefficients changed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNChgCoefs(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** gets number of constraint sides changed in presolving method of constraint handler */
extern
int SCIPconshdlrGetNChgSides(
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

/** needs constraint handler a constraint to be called? */
extern
Bool SCIPconshdlrNeedsCons(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** does the constraint handler perform presolving? */
extern
Bool SCIPconshdlrDoesPresolve(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** is constraint handler initialized? */
extern
Bool SCIPconshdlrIsInitialized(
   CONSHDLR*        conshdlr            /**< constraint handler */
   );



/*
 * Constraint set change methods
 */

/** frees constraint set change data and releases all included constraints */
extern
RETCODE SCIPconssetchgFree(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** adds constraint addition to constraint set changes, and captures constraint */
extern
RETCODE SCIPconssetchgAddAddedCons(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   NODE*            node,               /**< node that the constraint set change belongs to */
   CONS*            cons                /**< added constraint */
   );

/** adds constraint disabling to constraint set changes, and captures constraint */
extern
RETCODE SCIPconssetchgAddDisabledCons(
   CONSSETCHG**     conssetchg,         /**< pointer to constraint set change data structure */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< disabled constraint */
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
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
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
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable,         /**< should the constraint be removed from the LP due to aging or cleanup? */
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

/** globally removes constraint from all subproblems; removes constraint from the constraint set change data of the
 *  node, where it was created, or from the problem, if it was a problem constraint
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

/** resolves the given conflict var, that was deduced by the given constraint, by putting all "reason" variables
 *  leading to the deduction into the conflict queue with calls to SCIPaddConflictVar()
 */
extern
RETCODE SCIPconsResolveConflictVar(
   CONS*            cons,               /**< constraint that deduced the assignment */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< conflict variable, that was deduced by the constraint */
   );

/** checks single constraint for feasibility of the given solution */
extern
RETCODE SCIPconsCheck(
   CONS*            cons,               /**< constraint to check */
   const SET*       set,                /**< global SCIP settings */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   );

#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

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

/** returns TRUE iff constraint is active in the current node */
extern
Bool SCIPconsIsActive(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is enabled in the current node */
extern
Bool SCIPconsIsEnabled(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff the LP relaxation of constraint should be in the initial LP */
extern
Bool SCIPconsIsInitial(
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

/** returns TRUE iff constraint is only locally valid */
extern
Bool SCIPconsIsLocal(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is modifiable (subject to column generation) */
extern
Bool SCIPconsIsModifiable(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint should be removed from the LP due to aging or cleanup */
extern
Bool SCIPconsIsRemoveable(
   CONS*            cons                /**< constraint */
   );

/** returns TRUE iff constraint is belonging to original problem */
extern
Bool SCIPconsIsOriginal(
   CONS*            cons                /**< constraint */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPconsGetName(cons)           (cons)->name
#define SCIPconsGetHdlr(cons)           (cons)->conshdlr
#define SCIPconsGetData(cons)           (cons)->consdata
#define SCIPconsIsActive(cons)          ((cons)->updateactivate || ((cons)->active && !(cons)->updatedeactivate))
#define SCIPconsIsEnabled(cons)         ((cons)->updateenable || ((cons)->enabled && !(cons)->updatedisable))
#define SCIPconsIsInitial(cons)         (cons)->initial
#define SCIPconsIsSeparated(cons)       (cons)->separate
#define SCIPconsIsEnforced(cons)        (cons)->enforce
#define SCIPconsIsChecked(cons)         (cons)->check
#define SCIPconsIsPropagated(cons)      (cons)->propagate
#define SCIPconsIsLocal(cons)           (cons)->local
#define SCIPconsIsModifiable(cons)      (cons)->modifiable
#define SCIPconsIsRemoveable(cons)      (cons)->removeable
#define SCIPconsIsOriginal(cons)        (cons)->original

#endif




/*
 * Hash functions
 */

/** gets the key (i.e. the name) of the given constraint */
extern
DECL_HASHGETKEY(SCIPhashGetKeyCons);


#endif
