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
#pragma ident "@(#) $Id: objconshdlr.h,v 1.5 2003/12/08 13:32:05 bzfpfend Exp $"

/**@file   objconshdlr.h
 * @brief  C++ wrapper for constraint handlers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __OBJCONSHDLR_H__
#define __OBJCONSHDLR_H__


#include <cassert>

extern "C" 
{
#include "scip.h"
}


namespace scip
{

/** C++ wrapper object for constraint handlers */
class ObjConshdlr
{
public:
   /** name of the constraint handler */
   const char* const scip_name_;
   
   /** description of the constraint handler */
   const char* const scip_desc_;
   
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

   /** should the constraint handler be skipped, if no constraints are available? */
   const Bool scip_needscons_;

   /** default constructor */
   ObjConshdlr(
      const char*   name,               /**< name of constraint handler */
      const char*   desc,               /**< description of constraint handler */
      int           sepapriority,       /**< priority of the constraint handler for separation */
      int           enfopriority,       /**< priority of the constraint handler for constraint enforcing */
      int           checkpriority,      /**< priority of the constraint handler for checking infeasibility */
      int           sepafreq,           /**< frequency for separating cuts; zero means to separate only in the root node */
      int           propfreq,           /**< frequency for propagating domains; zero means only preprocessing propagation */
      Bool          needscons           /**< should the constraint handler be skipped, if no constraints are available? */
      )
      : scip_name_(name),
        scip_desc_(desc),
        scip_sepapriority_(sepapriority),
        scip_enfopriority_(enfopriority),
        scip_checkpriority_(checkpriority),
        scip_sepafreq_(sepafreq),
        scip_propfreq_(propfreq),
        scip_needscons_(needscons)
   {
   }

   /** destructor of constraint handler to free user data (called when SCIP is exiting) */
   virtual RETCODE scip_free(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr            /**< the constraint handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** initialization method of constraint handler (called when problem solving starts) */
   virtual RETCODE scip_init(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr            /**< the constraint handler itself */
      )
   {
      return SCIP_OKAY;
   }
   
   /** deinitialization method of constraint handler (called when problem solving exits) */
   virtual RETCODE scip_exit(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr            /**< the constraint handler itself */
      )
   {
      return SCIP_OKAY;
   }

   /** solving start notification method of constraint handler (called when presolving was finished)
    *
    *  This method is called when the presolving was finished and the branch and bound process is about to begin.
    *  It is called even when presolving is turned off.
    *  The constraint handler may use this call e.g. to clean up its presolving data, or to finally modify its constraints
    *  before the branch and bound process begins.
    *  Necessary constraint modifications that have to be performed even if presolving is turned off should be done here.
    */
   virtual RETCODE scip_solstart(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< final array of constraints in transformed problem */
      int           nconss              /**< final number of constraints in transformed problem */
      )
   {
      return SCIP_OKAY;
   }
   
   /** frees specific constraint data
    *
    *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
    *  the corresponding bound change event was not yet processed.
    */
   virtual RETCODE scip_delete(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONSDATA**    consdata            /**< pointer to the constraint data to free */
      )
   {
      return SCIP_OKAY;
   }

   /** transforms constraint data into data belonging to the transformed problem */
   virtual RETCODE scip_trans(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         sourcecons,         /**< source constraint to transform */
      CONS**        targetcons          /**< pointer to store created target constraint */
      ) = 0;

   /** LP initialization method of constraint handler
    *
    *  Puts the LP relaxations of all "initial" constraints into the LP. The method should scan the constraints
    *  array for constraints that are marked initial via calls to SCIPconsIsInitial() and put the LP relaxation
    *  of all initial constraints to the LP with calls to SCIPaddCut().
    */
   virtual RETCODE scip_initlp(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss              /**< number of constraints to process */
      )
   {
      return SCIP_OKAY;
   }

   /** separation method of constraint handler
    *
    *  Separates all constraints of the constraint handler. The method is called in the LP solution loop,
    *  which means that a valid LP solution exists.
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The separation
    *  method should process only the useful constraints in most runs, and only occasionally the remaining
    *  nconss - nusefulconss constraints.
    *
    *  possible return values for *result:
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> node is infeasible
    *  - SCIP_SEPARATED  : at least one cutting plane was generated
    *  - SCIP_REDUCEDDOM : no cutting plane was generated, but at least one domain was reduced
    *  - SCIP_CONSADDED  : no cutting plane or domain reductions, but at least one additional constraint was generated
    *  - SCIP_DIDNOTFIND : the separator searched, but did not find a cutting plane
    *  - SCIP_DIDNOTRUN  : the separator was skipped
    */
   virtual RETCODE scip_sepa(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      RESULT*       result              /**< pointer to store the result of the separation call */
      )
   {
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
    *  possible return values for *result:
    *  - SCIP_CUTOFF     : at least one constraint is infeasible, and it cannot be resolved -> node is infeasible
    *  - SCIP_SEPARATED  : a cutting plane was generated to resolve an infeasibility
    *  - SCIP_REDUCEDDOM : no cutting plane was generated, but at least one domain was reduced to resolve an infeasibility
    *  - SCIP_CONSADDED  : no cutting plane or domain reductions, but a constraint was generated to resolve an infeasibility
    *  - SCIP_BRANCHED   : no changes were made to the problem, but a branching was applied to resolve an infeasibility
    *  - SCIP_INFEASIBLE : at least one constraint is infeasible, but it was not resolved
    *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
    */
   virtual RETCODE scip_enfolp(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      RESULT*       result              /**< pointer to store the result of the enforcing call */
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
   virtual RETCODE scip_enfops(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      Bool          objinfeasible,      /**< is the solution infeasible anyway due to violating lower objective bound? */
      RESULT*       result              /**< pointer to store the result of the enforcing call */
      ) = 0;

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
    *  In some cases, integrality conditions or rows of the actual LP don't have to be checked, because their
    *  feasibility is already checked or implicitly given. In these cases, 'checkintegrality' or
    *  'checklprows' is FALSE.
    *
    *  possible return values for *result:
    *  - SCIP_INFEASIBLE : at least one constraint of the handler is infeasible
    *  - SCIP_FEASIBLE   : all constraints of the handler are feasible
    */
   virtual RETCODE scip_check(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      SOL*          sol,                /**< the solution to check feasibility for */
      Bool          checkintegrality,   /**< has integrality to be checked? */
      Bool          checklprows,        /**< have current LP rows to be checked? */
      RESULT*       result              /**< pointer to store the result of the feasibility checking call */
      ) = 0;

   /** domain propagation method of constraint handler
    *
    *  The first nusefulconss constraints are the ones, that are identified to likely be violated. The propagation
    *  method should process only the useful constraints in most runs, and only occasionally the remaining
    *  nconss - nusefulconss constraints.
    *
    *  possible return values for *result:
    *  - SCIP_CUTOFF     : at least one constraint is infeasible for the actual domains -> node is infeasible
    *  - SCIP_REDUCEDDOM : at least one domain reduction was found
    *  - SCIP_DIDNOTFIND : the propagator searched, but did not find any domain reductions
    *  - SCIP_DIDNOTRUN  : the propagator was skipped
    */
   virtual RETCODE scip_prop(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< number of constraints to process */
      int           nusefulconss,       /**< number of useful (non-obsolete) constraints to process */
      RESULT*       result              /**< pointer to store the result of the propagation call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** presolving method of constraint handler
    *
    *  The presolver should go through the variables and constraints and tighten the domains or
    *  constraints. Each tightening should increase the given total number of changes.
    *
    *  possible return values for *result:
    *  - SCIP_UNBOUNDED  : at least one variable is not bounded by any constraint in obj. direction -> problem is unbounded
    *  - SCIP_CUTOFF     : at least one constraint is infeasible in the variable's bounds -> problem is infeasible
    *  - SCIP_SUCCESS    : the presolver found a reduction
    *  - SCIP_DIDNOTFIND : the presolver searched, but did not find a presolving change
    *  - SCIP_DIDNOTRUN  : the presolver was skipped
    */
   virtual RETCODE scip_presol(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS**        conss,              /**< array of constraints to process */
      int           nconss,             /**< no. of constraints to process */
      int           nrounds,            /**< no. of presolving rounds already done */
      int           nnewfixedvars,      /**< no. of variables fixed since last call to presolving method */
      int           nnewaggrvars,       /**< no. of variables aggregated since last call to presolving method */
      int           nnewchgvartypes,    /**< no. of variable type changes since last call to presolving method */
      int           nnewchgbds,         /**< no. of variable bounds tightend since last call to presolving method */
      int           nnewholes,          /**< no. of domain holes added since last call to presolving method */
      int           nnewdelconss,       /**< no. of deleted constraints since last call to presolving method */
      int           nnewupgdconss,      /**< no. of upgraded constraints since last call to presolving method */
      int           nnewchgcoefs,       /**< no. of changed coefficients since last call to presolving method */
      int           nnewchgsides,       /**< no. of changed left or right hand sides since last call to presolving method */
      int*          nfixedvars,         /**< pointer to count total number of variables fixed of all presolvers */
      int*          naggrvars,          /**< pointer to count total number of variables aggregated of all presolvers */
      int*          nchgvartypes,       /**< pointer to count total number of variable type changes of all presolvers */
      int*          nchgbds,            /**< pointer to count total number of variable bounds tightend of all presolvers */
      int*          naddholes,          /**< pointer to count total number of domain holes added of all presolvers */
      int*          ndelconss,          /**< pointer to count total number of deleted constraints of all presolvers */
      int*          nupgdconss,         /**< pointer to count total number of upgraded constraints of all presolvers */
      int*          nchgcoefs,          /**< pointer to count total number of changed coefficients of all presolvers */
      int*          nchgsides,          /**< pointer to count total number of changed sides of all presolvers */
      RESULT*       result              /**< pointer to store the result of the presolving call */
      )
   {
      assert(result != NULL);
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

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
    *  to FALSE. It uses SCIPinferBinVar(scip, z, TRUE, c) to apply this assignment. In the conflict analysis, the
    *  constraint handler may be asked to resolve variable z with constraint c. With a call to SCIPvarGetLbLocal(z), 
    *  the handler can find out, that variable z is currently assigned to TRUE, and should call 
    *  SCIPaddConflictVar(scip, x) and SCIPaddConflictVar(scip, y) to tell SCIP, that the assignments to x and y were
    *  the reason for the deduction of z.
    */
   virtual RETCODE scip_rescvar(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons,               /**< the constraint that deduced the assignment of the conflict variable */
      VAR*          infervar            /**< the binary conflict variable that has to be resolved */
      )
   {
      errorMessage("conflict variable resolving method of constraint handler <%s> not implemented\n",
         SCIPconshdlrGetName(conshdlr));
      return SCIP_INVALIDCALL;
   }

   /** variable rounding lock method of constraint handler
    *
    *  This method is called, after a constraint is added to the transformed problem. It should lock the rounding
    *  of all associated variables with calls to SCIPvarLock(), depending on the way, the variable is involved
    *  in the constraint:
    *  - If the constraint may get violated by decreasing the value of a variable, it should call
    *    SCIPvarLock(var, nlockspos, nlocksneg), saying that rounding down is potentially rendering the
    *    (positive) constraint infeasible and rounding up is potentially rendering the negation of the constraint
    *    infeasible.
    *  - If the constraint may get violated by increasing the value of a variable, it should call
    *    SCIPvarLock(var, nlocksneg, nlockspos), saying that rounding down is potentially rendering the
    *    constraint's negation infeasible and rounding up is potentially rendering the constraint itself
    *    infeasible.
    *  - If the constraint may get violated by changing the variable in any direction, it should call
    *    SCIPvarLock(var, nlockspos + nlocksneg, nlockspos + nlocksneg).
    *
    *  Consider the linear constraint "3x -5y +2z <= 7" as an example. The variable rounding lock method of the
    *  linear constraint handler should call SCIPvarLock(x, nlocksneg, nlockspos), 
    *  SCIPvarLock(y, nlockspos, nlocksneg) and SCIPvarLock(z, nlocksneg, nlockspos) to tell SCIP, that
    *  rounding up of x and z and rounding down of y can destroy the feasibility of the constraint, while rounding
    *  down of x and z and rounding up of y can destroy the feasibility of the constraint's negation "3x -5y +2z > 7".
    *  A linear constraint "2 <= 3x -5y +2z <= 7" should call
    *  SCIPvarLock(..., nlockspos + nlocksneg, nlockspos + nlocksneg) on all variables, since rounding in both
    *  direction of each variable can destroy both the feasibility of the constraint and it's negation
    *  "3x -5y +2z < 2  or  3x -5y +2z > 7".
    *
    *  If the constraint itself contains other constraints as sub constraints (e.g. the "or" constraint concatenation
    *  "c(x) or d(x)"), the rounding lock methods of these constraints should be called in a proper way.
    *  - If the constraint may get violated by the violation of the sub constraint c, it should call
    *    SCIPlockConsVars(scip, c, nlockspos, nlocksneg), saying that infeasibility of c may lead to infeasibility of
    *    the (positive) constraint, and infeasibility of c's negation (i.e. feasibility of c) may lead to infeasibility
    *    of the constraint's negation (i.e. feasibility of the constraint).
    *  - If the constraint may get violated by the feasibility of the sub constraint c, it should call
    *    SCIPlockConsVars(scip, c, nlocksneg, nlockspos), saying that infeasibility of c may lead to infeasibility of
    *    the constraint's negation (i.e. feasibility of the constraint), and infeasibility of c's negation (i.e. feasibility
    *    of c) may lead to infeasibility of the (positive) constraint.
    *  - If the constraint may get violated by any change in the feasibility of the sub constraint c, it should call
    *    SCIPlockConsVars(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg).
    *
    *  Consider the or concatenation "c(x) or d(x)". The variable rounding lock method of the or constraint handler
    *  should call SCIPlockConsVars(scip, c, nlockspos, nlocksneg) and SCIPlockConsVars(scip, d, nlockspos, nlocksneg)
    *  to tell SCIP, that infeasibility of c and d can lead to infeasibility of "c(x) or d(x)".
    *
    *  As a second example, consider the equivalence constraint "y <-> c(x)" with variable y and constraint c. The
    *  constraint demands, that y == 1 if and only if c(x) is satisfied. The variable lock method of the corresponding
    *  constraint handler should call SCIPvarLock(y, nlockspos + nlocksneg, nlockspos + nlocksneg) and
    *  SCIPlockConsVars(scip, c, nlockspos + nlocksneg, nlockspos + nlocksneg), because any modification to the
    *  value of y or to the feasibility of c can alter the feasibility of the equivalence constraint.
    */
   virtual RETCODE scip_lock(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons,               /**< the constraint that should lock rounding of its variables */
      int           nlockspos,          /**< no. of times, the roundings should be locked for the constraint */
      int           nlocksneg           /**< no. of times, the roundings should be locked for the constraint's negation */
      ) = 0;

   /** variable rounding unlock method of constraint handler
    *
    *  This method is called, before a constraint is deleted from the transformed problem. It should unlock the rounding
    *  of all associated variables with calls to SCIPvarUnlock(), depending on the way, the variable is involved
    *  in the constraint:
    *  - If the constraint may get violated by decreasing the value of a variable, it should call
    *    SCIPvarUnlock(var, nunlockpos, nunlockneg).
    *  - If the constraint may get violated by increasing the value of a variable, it should call
    *    SCIPvarUnlock(var, nunlockneg, nunlockpos).
    *  - If the constraint may get violated by changing the variable in any direction, it should call
    *    SCIPvarUnlock(var, nunlockpos + nunlockneg, nunlockpos + nunlockneg).
    *
    *  If the constraint itself contains other constraints as sub constraints, the rounding lock methods of these
    *  constraints should be called in a proper way.
    *  - If the constraint may get violated by the violation of the sub constraint c, it should call
    *    SCIPunlockConsVars(scip, c, nunlockspos, nunlocksneg).
    *  - If the constraint may get violated by the feasibility of the sub constraint c, it should call
    *    SCIPunlockConsVars(scip, c, nunlocksneg, nunlockspos).
    *  - If the constraint may get violated by any change in the feasibility of the sub constraint c, it should call
    *    SCIPunlockConsVars(scip, c, nunlockspos + nunlocksneg, nunlockspos + nunlocksneg).
    *
    *  The unlocking method should exactly undo all lockings performed in the locking method of the constraint handler.
    */
   virtual RETCODE scip_unlock(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons,               /**< the constraint that should unlock rounding of its variables */
      int           nlockspos,          /**< no. of times, the roundings should be unlocked for the constraint */
      int           nlocksneg           /**< no. of times, the roundings should be unlocked for the constraint's negation */
      ) = 0;

   /** constraint activation notification method of constraint handler
    *
    *  WARNING! There may exist unprocessed events. For example, a variable's bound may have been already changed, but
    *  the corresponding bound change event was not yet processed.
    *
    *  This method is always called after a constraint of the constraint handler was activated. The constraint
    *  handler may use this call to update his own (statistical) data.
    */
   virtual RETCODE scip_active(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons                /**< the constraint that has been activated */
      )
   {
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
   virtual RETCODE scip_deactive(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons                /**< the constraint that has been deactivated */
      )
   {
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
   virtual RETCODE scip_enable(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons                /**< the constraint that has been enabled */
      )
   {
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
   virtual RETCODE scip_disable(
      SCIP*         scip,               /**< SCIP data structure */
      CONSHDLR*     conshdlr,           /**< the constraint handler itself */
      CONS*         cons                /**< the constraint that has been disabled */
      )
   {
      return SCIP_OKAY;
   }
};

} /* namespace scip */


   
/** creates the constraint handler for the given constraint handler object and includes it in SCIP
 *
 *  The method should be called in one of the following ways:
 *
 *   1. The user is resposible of deleting the object:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       MyConshdlr* myconshdlr = new MyConshdlr(...);
 *       CHECK_OKAY( SCIPincludeObjConshdlr(scip, &myconshdlr, FALSE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );
 *       delete myconshdlr;    // delete conshdlr AFTER SCIPfree() !
 *
 *   2. The object pointer is passed to SCIP and deleted by SCIP in the SCIPfree() call:
 *       CHECK_OKAY( SCIPcreate(&scip) );
 *       ...
 *       CHECK_OKAY( SCIPincludeObjConshdlr(scip, new MyConshdlr(...), TRUE) );
 *       ...
 *       CHECK_OKAY( SCIPfree(&scip) );  // destructor of MyConshdlr is called here
 */
extern
RETCODE SCIPincludeObjConshdlr(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjConshdlr* objconshdlr,      /**< constraint handler object */
   Bool             deleteobject        /**< should the constraint handler object be deleted when conshdlr is freed? */
   );

#endif
