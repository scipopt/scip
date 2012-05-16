/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   objconshdlr.h
 * @brief  C++ wrapper for constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_OBJCONSHDLR_H__
#define __SCIP_OBJCONSHDLR_H__


#include <cassert>
#include <cstring>

#include "scip/scip.h"
#include "objscip/objprobcloneable.h"

namespace scip
{

/**
 *  @brief C++ wrapper for constraint handlers
 *
 *  This class defines the interface for constraint handlers implemented in C++. Note that there are pure virtual
 *  functions (these have to be implemented). These functions are: scip_trans(), scip_enfolp(), scip_enfops(),
 *  scip_check(), and scip_lock().
 *
 *  - \ref CONS "Instructions for implementing a constraint handler"
 *  - \ref CONSHDLRS "List of available constraint handlers"
 *  - \ref type_cons.h "Corresponding C interface"
 */
class ObjConshdlr : public ObjProbCloneable
{
public:
   /*lint --e{1540}*/

   /** SCIP data structure */
   SCIP* scip_;

   /** name of the constraint handler */
   char* scip_name_;
   
   /** description of the constraint handler */
   char* scip_desc_;
   
   /** default separation priority of the constraint handler */
   const int scip_sepapriority_;

   /** default enforcing priority of the constraint handler */
   const int scip_enfopriority_;

   /** default checking priority of the constraint handler */
   const int scip_checkpriority_;

   /** default separation frequency of the constraint handler */
   const int scip_sepafreq_;

   /** default propagation frequency of the constraint handler */
   const int scip_propfreq_;

   /** default frequency of the constraint handler for eager evaluations in separation, propagation and enforcement */
   const int scip_eagerfreq_;

   /** maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
   const int scip_maxprerounds_;

   /** should separation method be delayed, if other separators found cuts? */
   const SCIP_Bool scip_delaysepa_;

   /** should propagation method be delayed, if other propagators found reductions? */
   const SCIP_Bool scip_delayprop_;

   /** should presolving method be delayed, if other presolvers found reductions? */
   const SCIP_Bool scip_delaypresol_;

   /** should the constraint handler be skipped, if no constraints are available? */
   const SCIP_Bool scip_needscons_;

   /** positions in the node solving loop where propagation method of constraint handler should be executed */
   const unsigned int scip_timingmask_;

   /** default constructor */
   ObjConshdlr(
      SCIP*              scip,               /**< SCIP data structure */
      const char*        name,               /**< name of constraint handler */
      const char*        desc,               /**< description of constraint handler */
      int                sepapriority,       /**< priority of the constraint handler for separation */
      int                enfopriority,       /**< priority of the constraint handler for constraint enforcing */
      int                checkpriority,      /**< priority of the constraint handler for checking infeasibility (and propagation) */
      int                sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
      int                propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
      int                eagerfreq,          /**< frequency for using all instead of only the useful constraints in separation,
                                              *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
      int                maxprerounds,       /**< maximal number of presolving rounds the constraint handler participates in (-1: no limit) */
      SCIP_Bool          delaysepa,          /**< should separation method be delayed, if other separators found cuts? */
      SCIP_Bool          delayprop,          /**< should propagation method be delayed, if other propagators found reductions? */
      SCIP_Bool          delaypresol,        /**< should presolving method be delayed, if other presolvers found reductions? */
      SCIP_Bool          needscons,          /**< should the constraint handler be skipped, if no constraints are available? */
      unsigned int       timingmask          /**< positions in the node solving loop where propagation method of constraint handler should be executed */
      )
      : scip_(scip),
        scip_name_(0),
        scip_desc_(0),
        scip_sepapriority_(sepapriority),
        scip_enfopriority_(enfopriority),
        scip_checkpriority_(checkpriority),
        scip_sepafreq_(sepafreq),
        scip_propfreq_(propfreq),
        scip_eagerfreq_(eagerfreq),
        scip_maxprerounds_(maxprerounds),
        scip_delaysepa_(delaysepa),
        scip_delayprop_(delayprop),
        scip_delaypresol_(delaypresol),
        scip_needscons_(needscons),
        scip_timingmask_(timingmask)
   {
      /* the macro SCIPduplicateMemoryArray does not need the first argument: */
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_name_, name, std::strlen(name)+1) );
      SCIP_CALL_ABORT( SCIPduplicateMemoryArray(scip_, &scip_desc_, desc, std::strlen(desc)+1) );
   }

   /** destructor */
   virtual ~ObjConshdlr()
   {
      /* the macro SCIPfreeMemoryArray does not need the first argument: */
      /*lint --e{64}*/
      SCIPfreeMemoryArray(scip_, &scip_name_);
      SCIPfreeMemoryArray(scip_, &scip_desc_);
   }

   /** destructor of constraint handler to free user data (called when SCIP is exiting) */
   virtual SCIP_RETCODE scip_free(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr            /**< the constraint handler itself */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** initialization method of constraint handler (called after problem has been transformed) */
   virtual SCIP_RETCODE scip_init(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints in transformed problem */
      int                nconss              /**< number of constraints in transformed problem */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** deinitialization method of constraint handler (called before transformed problem is freed) */
   virtual SCIP_RETCODE scip_exit(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints in transformed problem */
      int                nconss              /**< number of constraints in transformed problem */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** presolving initialization method of constraint handler (called when presolving is about to begin)
    *
    *  This method is called when the presolving process is about to begin, even if presolving is turned off.
    *  The constraint handler may use this call to initialize its presolving data, or to modify its constraints
    *  before the presolving process begins.
    *  Necessary constraint modifications that have to be performed even if presolving is turned off should be done here
    *  or in the presolving deinitialization call.
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_FEASIBLE   : no infeasibility nor unboundness could be found
    */
   virtual SCIP_RETCODE scip_initpre(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints in transformed problem */
      int                nconss,             /**< number of constraints in transformed problem */
      SCIP_Bool          isunbounded,        /**< was unboundedness already detected */
      SCIP_Bool          isinfeasible,       /**< was infeasibility already detected */
      SCIP_RESULT*       result              /**< pointer to store the result of the callback method */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);

      *result = SCIP_FEASIBLE;

      return SCIP_OKAY;
   }
   
   /** presolving deinitialization method of constraint handler (called after presolving has been finished)
    *
    *  This method is called after the presolving has been finished, even if presolving is turned off.
    *  The constraint handler may use this call e.g. to clean up its presolving data, or to finally modify its constraints
    *  before the branch and bound process begins.
    *  Necessary constraint modifications that have to be performed even if presolving is turned off should be done here.
    *  Besides necessary modifications and clean up, no time consuming operations should be done.
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_FEASIBLE   : no infeasibility nor unboundness could be found
    */
   virtual SCIP_RETCODE scip_exitpre(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< final array of constraints in transformed problem */
      int                nconss,             /**< final number of constraints in transformed problem */
      SCIP_Bool          isunbounded,        /**< was unboundedness already detected */
      SCIP_Bool          isinfeasible,       /**< was infeasibility already detected */
      SCIP_RESULT*       result              /**< pointer to store the result of the callback method */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);

      *result = SCIP_FEASIBLE;

      return SCIP_OKAY;
   }
   
   /** solving process initialization method of constraint handler (called when branch and bound process is about to begin)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  The constraint handler may use this call to initialize its branch and bound specific data.
    */
   virtual SCIP_RETCODE scip_initsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints of the constraint handler */
      int                nconss              /**< number of constraints of the constraint handler */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** solving process deinitialization method of constraint handler (called before branch and bound process data is freed)
    *
    *  This method is called before the branch and bound process is freed.
    *  The constraint handler should use this call to clean up its branch and bound data, in particular to release
    *  all LP rows that it has created or captured.
    */
   virtual SCIP_RETCODE scip_exitsol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints of the constraint handler */
      int                nconss              /**< number of constraints of the constraint handler */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }
   
   /** frees specific constraint data
    *
    *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
    *  the corresponding bound change event was not yet processed.
    */
   virtual SCIP_RETCODE scip_delete(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< the constraint belonging to the constraint data */
      SCIP_CONSDATA**    consdata            /**< pointer to the constraint data to free */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** transforms constraint data into data belonging to the transformed problem */
   virtual SCIP_RETCODE scip_trans(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         sourcecons,         /**< source constraint to transform */
      SCIP_CONS**        targetcons          /**< pointer to store created target constraint */
      ) = 0;

   /** LP initialization method of constraint handler (called before the initial LP relaxation at a node is solved)
    *
    *  Puts the LP relaxations of all "initial" constraints into the LP. The method should add a canonic LP relaxation
    *  of all given constraints to the LP with calls to SCIPaddCut().
    */
   virtual SCIP_RETCODE scip_initlp(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss              /**< number of constraints to process */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** separation method of constraint handler for LP solution
    *
    *  Separates all constraints of the constraint handler. The method is called in the LP solution loop,
    *  which means that a valid LP solution exists.
    *
    *  The first nusefulconss constraints are the ones that are identified to likely be violated. The separation
    *  method should process only the useful constraints in most runs, and only occasionally the remaining
    *  nconss - nusefulconss constraints.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
    *  - SCIP_DIDNOTRUN  : the separator was skipped
    *  - SCIP_DELAYED    : the separator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_sepalp(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** separation method of constraint handler for arbitrary primal solution
    *
    *  Separates all constraints of the constraint handler. The method is called outside the LP solution loop (e.g., by
    *  a relaxator or a primal heuristic), which means that there is no valid LP solution.
    *  Instead, the method should produce cuts that separate the given solution.
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
    *  method should process only the useful constraints in most runs, and only occasionally the remaining
    *  nconss - nusefulconss constraints.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_DIDNOTFIND : the separator searched, but did not find domain reductions, cutting planes, or cut constraints
    *  - SCIP_DIDNOTRUN  : the separator was skipped
    *  - SCIP_DELAYED    : the separator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_sepasol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_SOL*          sol,                /**< primal solution that should be separated */
      SCIP_RESULT*       result              /**< pointer to store the result of the separation call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** constraint enforcing method of constraint handler for LP solutions
    *
    *  The method is called at the end of the node processing loop for a node where the LP was solved.
    *  The LP solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
    *  branching, reducing a variable's domain to exclude the solution or separating the solution with a valid
    *  cutting plane.
    *
    *  The enforcing methods of the active constraint handlers are called in decreasing order of their enforcing
    *  priorities until the first constraint handler returned with the value SCIP_CUTOFF, SCIP_SEPARATED,
    *  SCIP_REDUCEDDOM, SCIP_CONSADDED, or SCIP_BRANCHED.
    *  The integrality constraint handler has an enforcing priority of zero. A constraint handler which can
    *  (or wants) to enforce its constraints only for integral solutions should have a negative enforcing priority
    *  (e.g. the alldiff-constraint can only operate on integral solutions).
    *  A constraint handler which wants to incorporate its own branching strategy even on non-integral
    *  solutions must have an enforcing priority greater than zero (e.g. the SOS-constraint incorporates
    *  SOS-branching on non-integral solutions).
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
    *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
    *  be enforced, if no violation was found in the useful constraints.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_SEPARATED  : a cutting plane was generated
    *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
    *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
    *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
    */
   virtual SCIP_RETCODE scip_enfolp(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
      SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
      ) = 0;

   /** constraint enforcing method of constraint handler for pseudo solutions
    *
    *  The method is called at the end of the node processing loop for a node where the LP was not solved.
    *  The pseudo solution has to be checked for feasibility. If possible, an infeasibility should be resolved by
    *  branching, reducing a variable's domain to exclude the solution or adding an additional constraint.
    *  Separation is not possible, since the LP is not processed at the current node. All LP informations like
    *  LP solution, slack values, or reduced costs are invalid and must not be accessed.
    *
    *  Like in the enforcing method for LP solutions, the enforcing methods of the active constraint handlers are
    *  called in decreasing order of their enforcing priorities until the first constraint handler returned with
    *  the value SCIP_CUTOFF, SCIP_REDUCEDDOM, SCIP_CONSADDED, SCIP_BRANCHED, or SCIP_SOLVELP.
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The enforcing
    *  method should process the useful constraints first. The other nconss - nusefulconss constraints should only
    *  be enforced, if no violation was found in the useful constraints.
    *
    *  If the pseudo solution's objective value is lower than the lower bound of the node, it cannot be feasible
    *  and the enforcing method may skip it's check and set *result to SCIP_DIDNOTRUN. However, it can also process
    *  its constraints and return any other possible result code.
    *
    *  possible return values for *result (if more than one applies, the first in the list should be used):
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_CONSADDED  : an additional constraint was generated
    *  - SCIP_REDUCEDDOM : a variable's domain was reduced
    *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
    *  - SCIP_SOLVELP    : at least one constraint is infeasible, and this can only be resolved by solving the SCIP_LP
    *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
    *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
    *  - SCIP_DIDNOTRUN  : the enforcement was skipped (only possible, if objinfeasible is true)
    */
   virtual SCIP_RETCODE scip_enfops(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_Bool          solinfeasible,      /**< was the solution already declared infeasible by a constraint handler? */
      SCIP_Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
      SCIP_RESULT*       result              /**< pointer to store the result of the enforcing call */
      ) = 0;

   /** feasibility check method of constraint handler for primal solutions
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
    *  In some cases, integrality conditions or rows of the current LP do not have to be checked, because their
    *  feasibility is already checked or implicitly given. In these cases, 'checkintegrality' or
    *  'checklprows' is FALSE.
    *
    *  possible return values for *result:
    *  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
    *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
    */
   virtual SCIP_RETCODE scip_check(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      SCIP_SOL*          sol,                /**< the solution to check feasibility for */
      SCIP_Bool          checkintegrality,   /**< has integrality to be checked? */
      SCIP_Bool          checklprows,        /**< have current LP rows to be checked? */
      SCIP_Bool          printreason,        /**< should the reason for the violation be printed? */
      SCIP_RESULT*       result              /**< pointer to store the result of the feasibility checking call */
      ) = 0;

   /** domain propagation method of constraint handler
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The propagation
    *  method should process only the useful constraints in most runs, and only occasionally the remaining
    *  nconss - nusefulconss constraints.
    *
    *  possible return values for *result:
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_REDUCEDDOM : at least one domain reduction was found
    *  - SCIP_DIDNOTFIND : the propagator searched, but did not find any domain reductions
    *  - SCIP_DIDNOTRUN  : the propagator was skipped
    *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_prop(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< number of constraints to process */
      int                nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      SCIP_PROPTIMING    proptiming,         /**< current point in the node solving process */
      SCIP_RESULT*       result              /**< pointer to store the result of the propagation call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** presolving method of constraint handler
    *
    *  The presolver should go through the variables and constraints and tighten the domains or
    *  constraints. Each tightening should increase the given total number of changes.
    *
    *  @note the counters state the changes since the last call including the changes of this presolving method during
    *        its call
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_SUCCESS    : the presolving method found a reduction
    *  - SCIP_DIDNOTFIND : the presolving method searched, but did not find a presolving change
    *  - SCIP_DIDNOTRUN  : the presolving method was skipped
    *  - SCIP_DELAYED    : the presolving method was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_presol(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints to process */
      int                nconss,             /**< no. of constraints to process */
      int                nrounds,            /**< no. of presolving rounds already done */
      int                nnewfixedvars,      /**< no. of variables fixed since last call to presolving method */
      int                nnewaggrvars,       /**< no. of variables aggregated since last call to presolving method */
      int                nnewchgvartypes,    /**< no. of variable type changes since last call to presolving method */
      int                nnewchgbds,         /**< no. of variable bounds tightend since last call to presolving method */
      int                nnewholes,          /**< no. of domain holes added since last call to presolving method */
      int                nnewdelconss,       /**< no. of deleted constraints since last call to presolving method */
      int                nnewaddconss,       /**< no. of added constraints since last call to presolving method */
      int                nnewupgdconss,      /**< no. of upgraded constraints since last call to presolving method */
      int                nnewchgcoefs,       /**< no. of changed coefficients since last call to presolving method */
      int                nnewchgsides,       /**< no. of changed left or right hand sides since last call to presolving method */
      int*               nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
      int*               naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
      int*               nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
      int*               nchgbds,            /**< pointer to count total number of variable bounds tightend of all presolvers */
      int*               naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
      int*               ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
      int*               naddconss,          /**< pointer to count total number of added constraints of all presolvers */
      int*               nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
      int*               nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
      int*               nchgsides,          /**< pointer to count total number of changed sides of all presolvers */
      SCIP_RESULT*       result              /**< pointer to store the result of the presolving call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** propagation conflict resolving method of constraint handler
    *
    *  This method is called during conflict analysis. If the constraint handler wants to support conflict analysis,
    *  it should call SCIPinferVarLbCons() or SCIPinferVarUbCons() in domain propagation instead of SCIPchgVarLb() or
    *  SCIPchgVarUb() in order to deduce bound changes on variables.
    *  In the SCIPinferVarLbCons() and SCIPinferVarUbCons() calls, the handler provides the constraint, that deduced the
    *  variable's bound change, and an integer value "inferinfo" that can be arbitrarily chosen.
    *  The propagation conflict resolving method can then be implemented, to provide a "reasons" for the bound
    *  changes, i.e. the bounds of variables at the time of the propagation, that forced the constraint to set the
    *  conflict variable's bound to its current value. It can use the "inferinfo" tag to identify its own propagation
    *  rule and thus identify the "reason" bounds. The bounds that form the reason of the assignment must then be provided
    *  by calls to SCIPaddConflictLb(), SCIPaddConflictUb(), SCIPaddConflictBd(), SCIPaddConflictRelaxedLb(),
    *  SCIPaddConflictRelaxedUb(), SCIPaddConflictRelaxedBd(), and/or SCIPaddConflictBinvar() in the propagation conflict
    *  resolving method.
    *
    *  For example, the logicor constraint c = "x or y or z" fixes variable z to TRUE (i.e. changes the lower bound of z
    *  to 1.0), if both, x and y, are assigned to FALSE (i.e. if the upper bounds of these variables are 0.0). It uses
    *  SCIPinferVarLbCons(scip, z, 1.0, c, 0) to apply this assignment (an inference information tag is not needed by the
    *  constraint handler and is set to 0).
    *  In the conflict analysis, the constraint handler may be asked to resolve the lower bound change on z with
    *  constraint c, that was applied at a time given by a bound change index "bdchgidx".
    *  With a call to SCIPvarGetLbAtIndex(z, bdchgidx, TRUE), the handler can find out, that the lower bound of
    *  variable z was set to 1.0 at the given point of time, and should call SCIPaddConflictUb(scip, x, bdchgidx) and
    *  SCIPaddConflictUb(scip, y, bdchgidx) to tell SCIP, that the upper bounds of x and y at this point of time were
    *  the reason for the deduction of the lower bound of z.
    *
    *  possible return values for *result:
    *  - SCIP_SUCCESS    : the conflicting bound change has been successfully resolved by adding all reason bounds
    *  - SCIP_DIDNOTFIND : the conflicting bound change could not be resolved and has to be put into the conflict set
    *
    *  @note it is sufficient to explain/resolve the relaxed bound
    */
   virtual SCIP_RETCODE scip_resprop(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< the constraint that deduced the bound change of the conflict variable */
      SCIP_VAR*          infervar,           /**< the conflict variable whose bound change has to be resolved */
      int                inferinfo,          /**< the user information passed to the corresponding SCIPinferVarLbCons()
                                              *   or SCIPinferVarUbCons() call */
      SCIP_BOUNDTYPE     boundtype,          /**< the type of the changed bound (lower or upper bound) */
      SCIP_BDCHGIDX*     bdchgidx,           /**< the index of the bound change, representing the point of time where the
                                              *   change took place */
      SCIP_Real          relaxedbd,          /**< the relaxed bound which is sufficient to be explained */
      SCIP_RESULT*       result              /**< pointer to store the result of the propagation conflict resolving resolving
                                              *   call */
      )
   {  /*lint --e{715}*/
      assert(result != NULL);
      *result = SCIP_DIDNOTFIND;
      return SCIP_OKAY;
   }

   /** variable rounding lock method of constraint handler
    *
    *  This method is called, after a constraint is added or removed from the transformed problem.
    *  It should update the rounding locks of all associated variables with calls to SCIPaddVarLocks(),
    *  depending on the way, the variable is involved in the constraint:
    *  - If the constraint may get violated by decreasing the value of a variable, it should call
    *    SCIPaddVarLocks(scip, var, nlockspos, nlocksneg), saying that rounding down is potentially rendering the
    *    (positive) constraint infeasible and rounding up is potentially rendering the negation of the constraint
    *    infeasible.
    *  - If the constraint may get violated by increasing the value of a variable, it should call
    *    SCIPaddVarLocks(scip, var, nlocksneg, nlockspos), saying that rounding up is potentially rendering the
    *    constraint's negation infeasible and rounding up is potentially rendering the constraint itself
    *    infeasible.
    *  - If the constraint may get violated by changing the variable in any direction, it should call
    *    SCIPaddVarLocks(scip, var, nlockspos + nlocksneg, nlockspos + nlocksneg).
    *
    *  Consider the linear constraint "3x -5y +2z <= 7" as an example. The variable rounding lock method of the
    *  linear constraint handler should call SCIPaddVarLocks(scip, x, nlocksneg, nlockspos), 
    *  SCIPaddVarLocks(scip, y, nlockspos, nlocksneg) and SCIPaddVarLocks(scip, z, nlocksneg, nlockspos) to tell SCIP,
    *  that rounding up of x and z and rounding down of y can destroy the feasibility of the constraint, while rounding
    *  down of x and z and rounding up of y can destroy the feasibility of the constraint's negation "3x -5y +2z > 7".
    *  A linear constraint "2 <= 3x -5y +2z <= 7" should call
    *  SCIPaddVarLocks(scip, ..., nlockspos + nlocksneg, nlockspos + nlocksneg) on all variables, since rounding in both
    *  directions of each variable can destroy both the feasibility of the constraint and it's negation
    *  "3x -5y +2z < 2  or  3x -5y +2z > 7".
    *
    *  If the constraint itself contains other constraints as sub constraints (e.g. the "or" constraint concatenation
    *  "c(x) or d(x)"), the rounding lock methods of these constraints should be called in a proper way.
    *  - If the constraint may get violated by the violation of the sub constraint c, it should call
    *    SCIPaddConsLocks(scip, c, nlockspos, nlocksneg), saying that infeasibility of c may lead to infeasibility of
    *    the (positive) constraint, and infeasibility of c's negation (i.e. feasibility of c) may lead to infeasibility
    *    of the constraint's negation (i.e. feasibility of the constraint).
    *  - If the constraint may get violated by the feasibility of the sub constraint c, it should call
    *    SCIPaddConsLocks(scip, c, nlocksneg, nlockspos), saying that infeasibility of c may lead to infeasibility of
    *    the constraint's negation (i.e. feasibility of the constraint), and infeasibility of c's negation (i.e. feasibility
    *    of c) may lead to infeasibility of the (positive) constraint.
    *  - If the constraint may get violated by any change in the feasibility of the sub constraint c, it should call
    *    SCIPaddConsLocks(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg).
    *
    *  Consider the or concatenation "c(x) or d(x)". The variable rounding lock method of the or constraint handler
    *  should call SCIPaddConsLocks(scip, c, nlockspos, nlocksneg) and SCIPaddConsLocks(scip, d, nlockspos, nlocksneg)
    *  to tell SCIP, that infeasibility of c and d can lead to infeasibility of "c(x) or d(x)".
    *
    *  As a second example, consider the equivalence constraint "y <-> c(x)" with variable y and constraint c. The
    *  constraint demands, that y == 1 if and only if c(x) is satisfied. The variable lock method of the corresponding
    *  constraint handler should call SCIPaddVarLocks(scip, y, nlockspos + nlocksneg, nlockspos + nlocksneg) and
    *  SCIPaddConsLocks(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg), because any modification to the
    *  value of y or to the feasibility of c can alter the feasibility of the equivalence constraint.
    */
   virtual SCIP_RETCODE scip_lock(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< the constraint that should lock rounding of its variables, or NULL if the
                                              *   constraint handler does not need constraints */
      int                nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
      int                nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
      ) = 0;

   /** constraint activation notification method of constraint handler
    *
    *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
    *  the corresponding bound change event was not yet processed.
    *
    *  This method is always called after a constraint of the constraint handler was activated. The constraint
    *  handler may use this call to update his own (statistical) data.
    */
   virtual SCIP_RETCODE scip_active(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons                /**< the constraint that has been activated */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint deactivation notification method of constraint handler
    *
    *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
    *  the corresponding bound change event was not yet processed.
    *
    *  This method is always called before a constraint of the constraint handler is deactivated. The constraint
    *  handler may use this call to update his own (statistical) data.
    */
   virtual SCIP_RETCODE scip_deactive(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons                /**< the constraint that has been deactivated */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint enabling notification method of constraint handler
    *
    *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
    *  the corresponding bound change event was not yet processed.
    *
    *  This method is always called after a constraint of the constraint handler was enabled. The constraint
    *  handler may use this call to update his own (statistical) data.
    */
   virtual SCIP_RETCODE scip_enable(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons                /**< the constraint that has been enabled */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint disabling notification method of constraint handler
    *
    *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
    *  the corresponding bound change event was not yet processed.
    *
    *  This method is always called before a constraint of the constraint handler is disabled. The constraint
    *  handler may use this call to update his own (statistical) data.
    */
   virtual SCIP_RETCODE scip_disable(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons                /**< the constraint that has been disabled */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }


   /** variable deletion method of constraint handler
    *
    *  This method goes through all constraints of the constraint handler and deletes all variables
    *  that were marked for deletion by SCIPdelVar().
    *
    *  possible return values for *result:
    *  - SCIP_CUTOFF     : the node is infeasible in the variable's bounds and can be cut off
    *  - SCIP_REDUCEDDOM : at least one domain reduction was found
    *  - SCIP_DIDNOTFIND : the propagator searched, but did not find any domain reductions
    *  - SCIP_DIDNOTRUN  : the propagator was skipped
    *  - SCIP_DELAYED    : the propagator was skipped, but should be called again
    */
   virtual SCIP_RETCODE scip_delvars(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        conss,              /**< array of constraints in transformed problem */
      int                nconss              /**< number of constraints in transformed problem */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint display method of constraint handler
    *
    *  The constraint handler should output a representation of the constraint into the given text file.
    */
   virtual SCIP_RETCODE scip_print(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< the constraint that should be displayed */
      FILE*              file                /**< the text file to store the information into */
      )
   {  /*lint --e{715}*/
      if ( file == NULL )
	 fprintf(stdout, "constraint handler <%s> does not support printing constraints\n", SCIPconshdlrGetName(conshdlr));
      else
	 fprintf(file, "constraint handler <%s> does not support printing constraints\n", SCIPconshdlrGetName(conshdlr));
      return SCIP_OKAY;
   }

   /** constraint copying method of constraint handler
    *
    *  The constraint handler can provide a copy method which copies a constraint from one SCIP data structure into a other
    *  SCIP data structure. If a copy of a constraint is created the constraint has to be captured (The capture is usually
    *  already done due to the creation of the constraint).
    * 
    *  If the copy process was a one to one the valid pointer can set to TRUE. Otherwise, you have to set this pointer to
    *  FALSE. In case all problem defining objects (constraint handlers and variable pricers) return a valid TRUE for all
    *  their copying calls, SCIP assumes that it is a overall one to one copy of the original instance. In this case any
    *  reductions made in the copied SCIP instance can be transfer to the original SCIP instance. If the valid pointer is
    *  set to TRUE and it was not one to one copy, it might happen that optimal solutions are cut off.
    *
    *  To get copy of variable in the target SCIP you should use the function SCIPgetVarCopy(). 
    */
   virtual SCIP_RETCODE scip_copy(
      SCIP*              scip,               /**< target SCIP data structure */
      SCIP_CONS**        cons,               /**< pointer to store the created target constraint */
      const char*        name,               /**< name of constraint, or NULL if the name of the source constraint should be used */
      SCIP*              sourcescip,         /**< source SCIP data structure */
      SCIP_CONSHDLR*     sourceconshdlr,     /**< source constraint handler of the source SCIP */
      SCIP_CONS*         sourcecons,         /**< source constraint of the source SCIP */
      SCIP_HASHMAP*      varmap,             /**< a SCIP_HASHMAP mapping variables of the source SCIP to corresponding
                                              *   variables of the target SCIP */
      SCIP_HASHMAP*      consmap,            /**< a hashmap to store the mapping of source constraints to the corresponding
                                              *   target constraints, or NULL */
      SCIP_Bool          initial,            /**< should the LP relaxation of constraint be in the initial LP? */
      SCIP_Bool          separate,           /**< should the constraint be separated during LP processing? */
      SCIP_Bool          enforce,            /**< should the constraint be enforced during node processing? */
      SCIP_Bool          check,              /**< should the constraint be checked for feasibility? */
      SCIP_Bool          propagate,          /**< should the constraint be propagated during node processing? */
      SCIP_Bool          local,              /**< is constraint only valid locally? */
      SCIP_Bool          modifiable,         /**< is constraint modifiable (subject to column generation)? */
      SCIP_Bool          dynamic,            /**< is constraint subject to aging? */
      SCIP_Bool          removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
      SCIP_Bool          stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
      SCIP_Bool          global,             /**< create a global or a local copy? */
      SCIP_Bool*         valid               /**< pointer to store whether the copying was valid or not */
      )
   {  /*lint --e{715}*/
      *valid = FALSE;
      return SCIP_OKAY;
   }
   
   /** constraint parsing method of constraint handler
    *
    *  The constraint handler should be able to parse the output created by the display method (SCIP_DECL_CONSDISABLE) 
    *  and to create constraint out of it.
    */
   virtual SCIP_RETCODE scip_parse(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS**        cons,               /**< pointer to store the created constraint */
      const char*        name,               /**< name of the constraint */
      const char*        str,                /**< string to parse */
      SCIP_Bool          initial,            /**< should the LP relaxation of constraint be in the initial LP? */
      SCIP_Bool          separate,           /**< should the constraint be separated during LP processing? */
      SCIP_Bool          enforce,            /**< should the constraint be enforced during node processing? */
      SCIP_Bool          check,              /**< should the constraint be checked for feasibility? */
      SCIP_Bool          propagate,          /**< should the constraint be propagated during node processing? */
      SCIP_Bool          local,              /**< is constraint only valid locally? */
      SCIP_Bool          modifiable,         /**< is constraint modifiable (subject to column generation)? */
      SCIP_Bool          dynamic,            /**< is constraint subject to aging? */
      SCIP_Bool          removable,          /**< should the relaxation be removed from the LP due to aging or cleanup? */
      SCIP_Bool          stickingatnode,     /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node? */
      SCIP_Bool*         success             /**< pointer to store whether the parsing was successful or not */
      )
   {  /*lint --e{715}*/
      return SCIP_OKAY;
   }

   /** constraint method of constraint handler which returns the variables (if possible)
    *
    *  The constraint handler can (this callback is optional) provide this callback to return the variables which are
    *  involved in that particular constraint. If this not possible, the variables should be copyied into the variables
    *  array and the success pointers has to be set to TRUE. Otherwise the success has to be set FALSE or the callback
    *  should not be implemented.
    */
   virtual SCIP_RETCODE scip_getvars(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< constraint for which the variables are wanted */
      SCIP_VAR**         vars,               /**< array to store/copy the involved variable of the constraint */
      int                varssize,           /**< available slots in vars array which is needed to check if the array is large enough */
      SCIP_Bool*         success             /**< pointer to store whether the variables are successfully copied */
      )
   {  /*lint --e{715}*/

      (*success) = FALSE;

      return SCIP_OKAY;
   }

   /** constraint method of constraint handler which returns the number of variables (if possible)
    *
    *  The constraint handler can (this callback is optional) provide this callback to return the number variable which
    *  are involved in that particular constraint. If this not possible, the success pointers has to be set to FALSE or
    *  the callback should not be implemented.
    */
   virtual SCIP_RETCODE scip_getnvars(
      SCIP*              scip,               /**< SCIP data structure */
      SCIP_CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      SCIP_CONS*         cons,               /**< constraint for which the number of variables is wanted */
      int*               nvars,              /**< pointer to store the number of variables */
      SCIP_Bool*         success             /**< pointer to store whether the constraint successfully returned the number of variables */
      )
   {  /*lint --e{715}*/

      (*nvars) = 0;
      (*success) = FALSE;

      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the constraint handler for the given constraint handler object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       MyConshdlr* myconshdlr = new MyConshdlr(...);
 *       SCIP_CALL( SCIPincludeObjConshdlr(scip, &myconshdlr, FALSE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );
 *       delete myconshdlr;    // delete conshdlr AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       SCIP_CALL( SCIPcreate(&scip) );
 *       ...
 *       SCIP_CALL( SCIPincludeObjConshdlr(scip, new MyConshdlr(...), TRUE) );
 *       ...
 *       SCIP_CALL( SCIPfree(&scip) );  // destructor of MyConshdlr is called here
 */
extern
SCIP_RETCODE SCIPincludeObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   scip::ObjConshdlr*    objconshdlr,        /**< constraint handler object */
   SCIP_Bool             deleteobject        /**< should the constraint handler object be deleted when conshdlr is freed? */
   );

/** returns the conshdlr object of the given name, or 0 if not existing */
extern
scip::ObjConshdlr* SCIPfindObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of constraint handler */
   );

/** returns the conshdlr object for the given constraint handler */
extern
scip::ObjConshdlr* SCIPgetObjConshdlr(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler */
   );

#endif
