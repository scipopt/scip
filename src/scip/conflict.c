/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: conflict.c,v 1.77 2004/11/17 12:53:48 bzfpfend Exp $"

/**@file   conflict.c
 * @brief  methods and datastructures for conflict analysis
 * @author Tobias Achterberg
 *
 * This file implements a conflict analysis method like the one used in modern
 * SAT solvers like zchaff. The algorithm works as follows:
 * 
 * Given is a set of bound changes that are not allowed being applied simultaneously,
 * because they render the current node infeasible (e.g. because a single constraint
 * is infeasible in the these bounds, or because the LP relaxation is infeasible).
 * The goal is to deduce a clause on binary variables -- a conflict clause -- representing
 * the "reason" for this conflict, i.e. the branching decisions or the deductions
 * (applied e.g. in domain propagation) that lead to the conflict. This clause
 * can then be added to the constraint set to help cutting off similar parts
 * of the branch and bound tree, that would lead to the same conflict.
 * A conflict clause can also be generated, if the conflict was detected by a
 * locally valid constraint. In this case, the resulting conflict clause is
 * also locally valid in the same depth as the conflict detecting constraint.
 *
 *  1. Put all the given bound changes to a priority queue, which is ordered,
 *     such that the bound change that was applied last due to branching or deduction
 *     is at the top of the queue. The variables in the queue are always active
 *     problem variables. Because non-binary variables must not exist in the final
 *     conflict clause, they are resolved first and put on the priority queue prior
 *     to the binary variables.
 *     Create an empty conflict set.
 *  2. Remove the top bound change b from the priority queue.
 *  3. (a) If the remaining queue is non-empty, and bound change b' (the one that is now
 *         on the top of the queue) was applied at the same depth level as b, and if
 *         b was a deduction with known inference reason, and if the inference constraint's
 *         valid depth is smaller or equal to the conflict detecting constraint's valid
 *         depth:
 *          - Resolve bound change b by asking the constraint that infered the
 *            bound change to put all the bound changes on the priority queue, that
 *            lead to the deduction of b.
 *            Note that these bound changes have at most the same inference depth
 *            level as b, and were deduced earlier than b.
 *     (b) Otherwise, the bound change b was a branching decision or a deduction with
 *         missing inference reason, or the inference constraint's validity is more local
 *         than the one of the conflict detecing constraint.
 *          - If b was a bound change on a non-binary variable, abort -- unresolved
 *            bound changes on non-binary variables cannot be handled, because the
 *            final conflict set must consist of only binary variables.
 *          - Otherwise, put the binary variable that was changed, or the negation of it
 *            into the conflict set, depending on which of them is currently fixed to
 *            FALSE (i.e., the conflict set consists of literals that cannot be FALSE
 *            altogether at the same time).
 *            Note that if the bound change was a branching, all deduced bound changes
 *            remaining in the priority queue have smaller inference depth level than b,
 *            since deductions are always applied after the branching decisions. However,
 *            there is the possibility, that b was a deduction, where the inference
 *            reason was not given or the inference constraint was too local.
 *            With this lack of information, we must treat the deduced bound change like
 *            a branching, and there may exist other deduced bound changes of the same
 *            inference depth level in the priority queue.
 *  4. If priority queue is non-empty, goto step 2.
 *  5. The conflict set represents the conflict clause saying that at least one
 *     of the binary conflict variables must be set to TRUE.
 *     The conflict set is then passed to the conflict handlers, that may create 
 *     a corresponding constraint (e.g. a logicor constraint) out of these conflict
 *     variables and add it to the problem.
 *
 * If all deduced bound changes come with (global) inference information, depending on
 * the conflict analyzing strategy, the resulting conflict set has the following property:
 *  - 1-FirstUIP: In the depth level where the conflict was found, at most one variable
 *    assigned at that level is member of the conflict set. This conflict variable is the
 *    first unique implication point of its depth level (FUIP).
 *  - All-FirstUIP: For each depth level, at most one variable assigned at that level is
 *    member of the conflict set. This conflict variable is the first unique implication
 *    point of its depth level (FUIP).
 *
 * The user has to do the following to get the conflict analysis running in its
 * current implementation:
 *  - A constraint handler or propagator supporting the conflict analysis must implement
 *    the CONSRESPROP/PROPRESPROP call, that processes a bound change inference b and puts all
 *    the reason bounds leading to the application of b with calls to
 *    SCIPaddConflictBound() on the conflict queue (algorithm step 3.(a)).
 *  - If the current bounds lead to a deduction of a bound change (e.g. in domain
 *    propagation), a constraint handler should call SCIPinferVarLbCons() or
 *    SCIPinferVarUbCons(), thus providing the constraint that infered the bound change.
 *    A propagator should call SCIPinferVarLbProp() or SCIPinferVarUbProp() instead,
 *    thus providing a pointer to itself.
 *  - If (in the current bounds) an infeasibility is detected, the constraint handler or
 *    propagator should
 *     1. call SCIPinitConflictAnalysis() to initialize the conflict queue,
 *     2. call SCIPaddConflictBound() for each bound that lead to the conflict,
 *     3. call SCIPanalyzeConflictCons() or SCIPanalyzeConflict() to analyze the conflict
 *        and add an appropriate conflict constraint.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "stat.h"
#include "clock.h"
#include "vbc.h"
#include "lpi.h"
#include "misc.h"
#include "paramset.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "tree.h"
#include "scip.h"
#include "conflict.h"
#include "cons.h"
#include "prop.h"

#include "struct_conflict.h"



/*
 * Conflict Handler
 */

/** compares two conflict handlers w. r. to their priority */
DECL_SORTPTRCOMP(SCIPconflicthdlrComp)
{  /*lint --e{715}*/
   return ((CONFLICTHDLR*)elem2)->priority - ((CONFLICTHDLR*)elem1)->priority;
}

/** method to call, when the priority of a conflict handler was changed */
static
DECL_PARAMCHGD(paramChgdConflicthdlrPriority)
{  /*lint --e{715}*/
   PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetConflicthdlrPriority() to mark the conflicthdlrs unsorted */
   CHECK_OKAY( SCIPsetConflicthdlrPriority(scip, (CONFLICTHDLR*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a conflict handler */
RETCODE SCIPconflicthdlrCreate(
   CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of conflict handler */
   const char*      desc,               /**< description of conflict handler */
   int              priority,           /**< priority of the conflict handler */
   DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   )
{
   char paramname[MAXSTRLEN];
   char paramdesc[MAXSTRLEN];

   assert(conflicthdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   ALLOC_OKAY( allocMemory(conflicthdlr) );
   ALLOC_OKAY( duplicateMemoryArray(&(*conflicthdlr)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*conflicthdlr)->desc, desc, strlen(desc)+1) );
   (*conflicthdlr)->priority = priority;
   (*conflicthdlr)->conflictfree = conflictfree;
   (*conflicthdlr)->conflictinit = conflictinit;
   (*conflicthdlr)->conflictexit = conflictexit;
   (*conflicthdlr)->conflictexec = conflictexec;
   (*conflicthdlr)->conflicthdlrdata = conflicthdlrdata;
   (*conflicthdlr)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "conflict/%s/priority", name);
   sprintf(paramdesc, "priority of conflict handler <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
         &(*conflicthdlr)->priority, priority, INT_MIN, INT_MAX, 
         paramChgdConflicthdlrPriority, (PARAMDATA*)(*conflicthdlr)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of conflict handler */
RETCODE SCIPconflicthdlrFree(
   CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conflicthdlr != NULL);
   assert(*conflicthdlr != NULL);
   assert(!(*conflicthdlr)->initialized);
   assert(SCIPgetStage(scip) == SCIP_STAGE_INIT);

   /* call destructor of conflict handler */
   if( (*conflicthdlr)->conflictfree != NULL )
   {
      CHECK_OKAY( (*conflicthdlr)->conflictfree(scip, *conflicthdlr) );
   }

   freeMemoryArray(&(*conflicthdlr)->name);
   freeMemoryArray(&(*conflicthdlr)->desc);
   freeMemory(conflicthdlr);

   return SCIP_OKAY;
}

/** calls init method of conflict handler */
RETCODE SCIPconflicthdlrInit(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conflicthdlr != NULL);

   if( conflicthdlr->initialized )
   {
      errorMessage("conflict handler <%s> already initialized\n", conflicthdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call initialization method of conflict handler */
   if( conflicthdlr->conflictinit != NULL )
   {
      CHECK_OKAY( conflicthdlr->conflictinit(scip, conflicthdlr) );
   }
   conflicthdlr->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of conflict handler */
RETCODE SCIPconflicthdlrExit(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(conflicthdlr != NULL);

   if( !conflicthdlr->initialized )
   {
      errorMessage("conflict handler <%s> not initialized\n", conflicthdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call deinitialization method of conflict handler */
   if( conflicthdlr->conflictexit != NULL )
   {
      CHECK_OKAY( conflicthdlr->conflictexit(scip, conflicthdlr) );
   }
   conflicthdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** calls execution method of conflict handler */
RETCODE SCIPconflicthdlrExec(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP*            scip,               /**< SCIP data structure */   
   NODE*            node,               /**< node to add conflict constraint to */
   VAR**            conflictvars,       /**< variables of the conflict set */
   int              nconflictvars,      /**< number of variables in the conflict set */
   Bool             local,              /**< is the conflict set only valid locally? */
   Bool             resolved,           /**< is the conflict set already used to create a constraint? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conflicthdlr != NULL);
   assert(conflictvars != NULL || nconflictvars == 0);
   assert(result != NULL);

   /* call solution start method of conflict handler */
   *result = SCIP_DIDNOTRUN;
   if( conflicthdlr->conflictexec != NULL )
   {
      CHECK_OKAY( conflicthdlr->conflictexec(scip, conflicthdlr, node, conflictvars, nconflictvars, local, resolved, 
            result) );

      if( *result != SCIP_CONSADDED
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         errorMessage("execution method of conflict handler <%s> returned invalid result <%d>\n", 
            conflicthdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }
   
   return SCIP_OKAY;
}

/** gets user data of conflict handler */
CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->conflicthdlrdata;
}

/** sets user data of conflict handler; user has to free old data in advance! */
void SCIPconflicthdlrSetData(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflicthdlrdata = conflicthdlrdata;
}

/** gets name of conflict handler */
const char* SCIPconflicthdlrGetName(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->name;
}

/** gets description of conflict handler */
const char* SCIPconflicthdlrGetDesc(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->desc;
}

/** gets priority of conflict handler */
int SCIPconflicthdlrGetPriority(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->priority;
}

/** sets priority of conflict handler */
void SCIPconflicthdlrSetPriority(
   CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the conflict handler */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);
   
   conflicthdlr->priority = priority;
   set->conflicthdlrssorted = FALSE;
}

/** is conflict handler initialized? */
Bool SCIPconflicthdlrIsInitialized(
   CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->initialized;
}




/*
 * Propagation Conflict Analysis
 */

/** resizes conflictvars array to be able to store at least num constraints */
static
RETCODE conflictEnsureConflictvarsMem(
   CONFLICT*        conflict,           /**< conflict analysis data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(conflict != NULL);

   if( num > conflict->conflictvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&conflict->conflictvars, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&conflict->conflictvardepths, newsize) );
      conflict->conflictvarssize = newsize;
   }
   assert(num <= conflict->conflictvarssize);

   return SCIP_OKAY;
}

/** puts given locally fixed binary variable or its negation into conflict set, depending which one is fixed to zero */
static
RETCODE conflictAddConflictVar(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR*             var,                /**< variable that should enter the conflict set */
   Bool             literalvalue,       /**< current value of the literal, that should enter the conflict set */
   Bool             temporary           /**< should the variable be added only temporary to the conflict set? */
   )
{
   int fixdepth;
   int i;

   assert(conflict != NULL);
   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(SCIPsetIsEQ(set, SCIPvarGetLbLP(var), SCIPvarGetUbLP(var)));
   assert(SCIPvarGetLbLP(var) > 0.5 || SCIPvarGetUbLP(var) < 0.5);

   /* get the depth level, at which the variable was fixed */
   fixdepth = SCIPvarGetLastBdchgDepth(var);
   assert(fixdepth != -1); /* variables fixed in preprocessing shouldn't occur here */
   if( fixdepth == -2 )
   {
      /* the variable is not fixed at the focus node: it was fixed due to diving, probing, or strong branching */
      fixdepth = INT_MAX;
   }

   /* choose variable or its negation, such that the literal is currently fixed to the desired value */
   if( (SCIPvarGetLbLP(var) > 0.5) != literalvalue )
   {
      /* variable is fixed to the opposite value -> use the negated variable as conflict literal */
      CHECK_OKAY( SCIPvarNegate(var, memhdr, set, stat, &var) );
   }
   assert(SCIPsetIsEQ(set, SCIPvarGetLbLP(var), SCIPvarGetUbLP(var)));
   assert((SCIPvarGetLbLP(var) > 0.5) == literalvalue);

   debugMessage("putting variable <%s> fixed to %d at depth %d to conflict set (temporary: %d)\n",
      SCIPvarGetName(var), (SCIPvarGetLbLP(var) > 0.5), fixdepth, temporary);

   if( var->conflictsetcount == conflict->count )
   {
      /* variable is already member of the conflict set */
#ifndef NDEBUG
      for( i = 0; i < conflict->nconflictvars + conflict->ntmpconflictvars && conflict->conflictvars[i] != var; ++i )
      {}
      assert(i < conflict->nconflictvars + conflict->ntmpconflictvars);
      assert(conflict->conflictvars[i] == var);
#endif
      return SCIP_OKAY;
   }

   /* put variable to conflict set, sorted by non-increasing fixing depths */
   CHECK_OKAY( conflictEnsureConflictvarsMem(conflict, set, conflict->nconflictvars + conflict->ntmpconflictvars + 1) );
   if( temporary )
   {
      /* sort the variable into the temporary part of the conflict vars array */
      for( i = conflict->nconflictvars + conflict->ntmpconflictvars;
           i > conflict->nconflictvars && conflict->conflictvardepths[i-1] < fixdepth; --i )
      {
         conflict->conflictvars[i] = conflict->conflictvars[i-1];
         conflict->conflictvardepths[i] = conflict->conflictvardepths[i-1];
      }
      conflict->conflictvars[i] = var;
      conflict->conflictvardepths[i] = fixdepth;
      conflict->ntmpconflictvars++;
   }
   else
   {
      assert(conflict->ntmpconflictvars == 0);

      /* sort the variable into the conflict vars array */
      for( i = conflict->nconflictvars; i > 0 && conflict->conflictvardepths[i-1] < fixdepth; --i )
      {
         conflict->conflictvars[i] = conflict->conflictvars[i-1];
         conflict->conflictvardepths[i] = conflict->conflictvardepths[i-1];
      }
      conflict->conflictvars[i] = var;
      conflict->conflictvardepths[i] = fixdepth;
      conflict->nconflictvars++;
   }

   /* mark the active problem variable to belong to the current conflict */
   var->conflictsetcount = conflict->count;

   return SCIP_OKAY;
}

/** removes all temporary variables from the conflict set */
static
void conflictRemoveTemporaryConflictVars(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   int i;

   assert(conflict != NULL);
   assert(conflict->nconflictvars >= 0);
   assert(conflict->nconflictvars == 0 || conflict->conflictvardepths != NULL);

   for( i = conflict->nconflictvars; i < conflict->nconflictvars + conflict->ntmpconflictvars; ++i )
      conflict->conflictvars[i]->conflictsetcount = 0;
   conflict->ntmpconflictvars = 0;
}

/** after the conflict set is computated, gets the depth at which now a propagation could be performed with this
 *  conflict set; this is the depth of the next to last fixed variable
 */
static
int conflictGetPropagateDepth(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   int lastdepth;
   int lastbutonedepth;

   assert(conflict != NULL);

#ifndef NDEBUG
   {
      int i;
      for( i = 0; i < conflict->nconflictvars; ++i )
      {
         assert(conflict->conflictvardepths[i] >= 0);
         assert(conflict->conflictvardepths[i] == SCIPvarGetLastBdchgDepth(conflict->conflictvars[i])
            || (conflict->conflictvardepths[i] == INT_MAX && SCIPvarGetLastBdchgDepth(conflict->conflictvars[i]) == -2));
         assert(i == conflict->nconflictvars-1 || conflict->conflictvardepths[i] >= conflict->conflictvardepths[i+1]);
      }
      for( i = conflict->nconflictvars; i < conflict->nconflictvars + conflict->ntmpconflictvars; ++i )
      {
         assert(conflict->conflictvardepths[i] >= 0);
         assert(conflict->conflictvardepths[i] == SCIPvarGetLastBdchgDepth(conflict->conflictvars[i])
            || (conflict->conflictvardepths[i] == INT_MAX && SCIPvarGetLastBdchgDepth(conflict->conflictvars[i]) == -2));
         assert(i == conflict->nconflictvars + conflict->ntmpconflictvars - 1
            || conflict->conflictvardepths[i] >= conflict->conflictvardepths[i+1]);
      }
   }
#endif
   
   /* the last fixed variable could already be fixed to the opposite value in the depth of the previously
    * fixed conflict variable by using this conflict clause
    */
   lastdepth = 0;
   lastbutonedepth = 0;
   if( conflict->nconflictvars >= 1 )
      lastdepth = conflict->conflictvardepths[0];
   if( conflict->nconflictvars >= 2 )
      lastbutonedepth = conflict->conflictvardepths[1];
   if( conflict->ntmpconflictvars >= 1 )
   {
      if( conflict->conflictvardepths[conflict->nconflictvars] > lastdepth )
      {
         lastbutonedepth = lastdepth;
         lastdepth = conflict->conflictvardepths[conflict->nconflictvars];
      }
      else if( conflict->conflictvardepths[conflict->nconflictvars] > lastbutonedepth )
         lastbutonedepth = conflict->conflictvardepths[conflict->nconflictvars];
   }
   if( conflict->ntmpconflictvars >= 2 )
   {
      if( conflict->conflictvardepths[conflict->nconflictvars+1] > lastdepth )
      {
         lastbutonedepth = lastdepth;
         lastdepth = conflict->conflictvardepths[conflict->nconflictvars+1];
      }
      else if( conflict->conflictvardepths[conflict->nconflictvars+1] > lastbutonedepth )
         lastbutonedepth = conflict->conflictvardepths[conflict->nconflictvars+1];
   }
   assert(lastdepth >= lastbutonedepth);

   return lastbutonedepth;
}

/** compares two conflict set entries, such that bound changes infered later are
 *  ordered prior to ones that were infered earlier
 */
static
DECL_SORTPTRCOMP(conflictBdchginfoComp)
{  /*lint --e{715}*/
   BDCHGINFO* bdchginfo1;
   BDCHGINFO* bdchginfo2;
   
   bdchginfo1 = (BDCHGINFO*)elem1;
   bdchginfo2 = (BDCHGINFO*)elem2;
   assert(bdchginfo1 != NULL);
   assert(bdchginfo2 != NULL);

   assert((SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo1)) == SCIP_VARTYPE_BINARY)
      == (SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo2)) == SCIP_VARTYPE_BINARY));

   if( !SCIPbdchgidxIsEarlierNonNull(SCIPbdchginfoGetIdx(bdchginfo1), SCIPbdchginfoGetIdx(bdchginfo2)) )
      return -1;
   else
      return +1;
}

/** creates conflict analysis data for propagation conflicts */
RETCODE SCIPconflictCreate(
   CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflict != NULL);

   ALLOC_OKAY( allocMemory(conflict) );

   CHECK_OKAY( SCIPclockCreate(&(*conflict)->propanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*conflict)->lpanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*conflict)->sbanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPclockCreate(&(*conflict)->pseudoanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
#if 0
   CHECK_OKAY( SCIPlpiCreate(&(*conflict)->lpi, "LPconflict") );
#endif
   CHECK_OKAY( SCIPpqueueCreate(&(*conflict)->binbdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp) );
   CHECK_OKAY( SCIPpqueueCreate(&(*conflict)->nonbinbdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac, 
         conflictBdchginfoComp) );
   (*conflict)->conflictvars = NULL;
   (*conflict)->conflictvardepths = NULL;
   (*conflict)->conflictvarssize = 0;
   (*conflict)->nconflictvars = 0;
   (*conflict)->ntmpconflictvars = 0;
   (*conflict)->count = 0;
   (*conflict)->npropcalls = 0;
   (*conflict)->npropconfclauses = 0;
   (*conflict)->npropconfliterals = 0;
   (*conflict)->npropreconvclauses = 0;
   (*conflict)->npropreconvliterals = 0;
   (*conflict)->nlpcalls = 0;
   (*conflict)->nlpconfclauses = 0;
   (*conflict)->nlpconfliterals = 0;
   (*conflict)->nlpreconvclauses = 0;
   (*conflict)->nlpreconvliterals = 0;
   (*conflict)->nlpiterations = 0;
   (*conflict)->nsbcalls = 0;
   (*conflict)->nsbconfclauses = 0;
   (*conflict)->nsbconfliterals = 0;
   (*conflict)->nsbreconvclauses = 0;
   (*conflict)->nsbreconvliterals = 0;
   (*conflict)->nsbiterations = 0;
   (*conflict)->npseudocalls = 0;
   (*conflict)->npseudoconfclauses = 0;
   (*conflict)->npseudoconfliterals = 0;
   (*conflict)->npseudoreconvclauses = 0;
   (*conflict)->npseudoreconvliterals = 0;

   return SCIP_OKAY;
}

/** frees conflict analysis data for propagation conflicts */
RETCODE SCIPconflictFree(
   CONFLICT**       conflict            /**< pointer to conflict analysis data */
   )
{
   assert(conflict != NULL);
   assert(*conflict != NULL);


   SCIPclockFree(&(*conflict)->propanalyzetime);
   SCIPclockFree(&(*conflict)->lpanalyzetime);
   SCIPclockFree(&(*conflict)->sbanalyzetime);
   SCIPclockFree(&(*conflict)->pseudoanalyzetime);
#if 0
   CHECK_OKAY( SCIPlpiFree(&(*conflict)->lpi) );
#endif
   SCIPpqueueFree(&(*conflict)->binbdchgqueue);
   SCIPpqueueFree(&(*conflict)->nonbinbdchgqueue);
   freeMemoryArrayNull(&(*conflict)->conflictvars);
   freeMemoryArrayNull(&(*conflict)->conflictvardepths);
   freeMemory(conflict);


   return SCIP_OKAY;
}

/** initializes the propagation conflict analysis by clearing the conflict candidate queue */
RETCODE SCIPconflictInit(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   debugMessage("initializing conflict analysis\n");

   /* clear the conflict candidate queues and the conflict set */
   SCIPpqueueClear(conflict->binbdchgqueue);
   SCIPpqueueClear(conflict->nonbinbdchgqueue);
   conflict->nconflictvars = 0;
   conflict->ntmpconflictvars = 0;

   /* increase the conflict set counter, such that variables of new conflict set are labeled with this new counter */
   conflict->count++;

   return SCIP_OKAY;
}

/** adds given bound change information to the conflict candidate queue */
static
RETCODE conflictAddBdchginfo(
   CONFLICT*        conflict,           /**< conflict analysis data */
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(conflict != NULL);
   assert(bdchginfo != NULL);

   /* put candidate in the appropriate priority queue */
   if( SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) == SCIP_VARTYPE_BINARY )
   {
      CHECK_OKAY( SCIPpqueueInsert(conflict->binbdchgqueue, (void*)bdchginfo) );
   }
   else
   {
      CHECK_OKAY( SCIPpqueueInsert(conflict->nonbinbdchgqueue, (void*)bdchginfo) );
   }

   return SCIP_OKAY;
}

/** adds variable's bound to conflict candidate queue */
RETCODE SCIPconflictAddBound(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR*             var,                /**< problem variable */
   BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   BDCHGINFO* bdchginfo;
   Real scalar;
   Real constant;

   assert(conflict != NULL);
   assert(var != NULL);

   /* get active problem variable */
   scalar = 1.0;
   constant = 0.0;
   CHECK_OKAY( SCIPvarGetProbvarSum(&var, &scalar, &constant) );

   /* we can ignore fixed variables */
   if( var == NULL )
      return SCIP_OKAY;

   assert(!SCIPsetIsZero(set, scalar));

   /* if the scalar of the aggregation is negative, we have to switch the bound type */
   if( scalar < 0.0 )
      boundtype = (boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BOUNDTYPE_UPPER : SCIP_BOUNDTYPE_LOWER);

   /* get bound change information */
   bdchginfo = SCIPvarGetBdchgInfo(var, boundtype, bdchgidx, FALSE);

   /* if bound of variable was not changed, we can ignore the conflicting bound */
   if( bdchginfo == NULL )
      return SCIP_OKAY;

   debugMessage(" -> adding bound <%s> %s %g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] to candidates\n",
      SCIPvarGetName(var),
      boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
      boundtype == SCIP_BOUNDTYPE_LOWER ?
      SCIPvarGetLbAtIndex(var, bdchgidx, FALSE) : SCIPvarGetUbAtIndex(var, bdchgidx, FALSE),
      SCIPvarGetStatus(var), SCIPvarGetType(var),
      SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo), 
      SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
      : (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
         ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
         : (SCIPbdchginfoGetInferProp(bdchginfo) != NULL ? SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo))
            : "none")),
      SCIPbdchginfoGetChgtype(bdchginfo) != SCIP_BOUNDCHGTYPE_BRANCHING ? SCIPbdchginfoGetInferInfo(bdchginfo) : -1);

   /* the local bound change may be resolved and has to be put on the candidate queue;
    * we even put bound changes without inference information on the queue in order to automatically
    * eliminate multiple insertions of the same bound change
    */
   assert(SCIPbdchginfoGetVar(bdchginfo) == var);
   assert(SCIPbdchginfoGetBoundtype(bdchginfo) == boundtype);
   assert(SCIPbdchginfoGetDepth(bdchginfo) >= 0);
   assert(SCIPbdchginfoGetPos(bdchginfo) >= 0);
   assert(SCIPbdchgidxIsEarlier(SCIPbdchginfoGetIdx(bdchginfo), bdchgidx));

   /* put bound change information into priority queue */
   CHECK_OKAY( conflictAddBdchginfo(conflict, bdchginfo) );

   return SCIP_OKAY;
}

/** removes and returns next conflict analysis candidate from the candidate queues */
static
BDCHGINFO* conflictRemoveCand(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   BDCHGINFO* bdchginfo;

   assert(conflict != NULL);

   bdchginfo = (BDCHGINFO*)(SCIPpqueueRemove(conflict->nonbinbdchgqueue));
   if( bdchginfo == NULL )
      bdchginfo = (BDCHGINFO*)(SCIPpqueueRemove(conflict->binbdchgqueue));

   return bdchginfo;
}

/** returns next conflict analysis candidate from the candidate queues without removing it */
static
BDCHGINFO* conflictFirstCand(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   BDCHGINFO* bdchginfo;

   assert(conflict != NULL);

   bdchginfo = (BDCHGINFO*)(SCIPpqueueFirst(conflict->nonbinbdchgqueue));
   if( bdchginfo == NULL )
      bdchginfo = (BDCHGINFO*)(SCIPpqueueFirst(conflict->binbdchgqueue));

   return bdchginfo;
}

/** calls the conflict handlers in order to create a conflict clause */
static
RETCODE conflictAddClause(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   int              validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   int              maxsize,            /**< maximal size of conflict set */
   Bool*            success,            /**< pointer to store whether the conflict set is valid */
   int*             nliterals           /**< pointer to store the number of literals in the generated clause */
   )
{
   BDCHGINFO** bdchginfos;
   int nbdchginfos;
   int currentdepth;
   int d;
   int i;

   assert(conflict != NULL);
   assert(conflict->nconflictvars >= 0);
   assert(conflict->nconflictvars == 0 || conflict->conflictvars != NULL);
   assert(conflict->nconflictvars == 0 || conflict->conflictvardepths != NULL);
   assert(conflict->ntmpconflictvars == 0);
   assert(SCIPpqueueNElems(conflict->nonbinbdchgqueue) == 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(success != NULL);
   assert(nliterals != NULL);

   *success = FALSE;
   *nliterals = 0;

   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(0 <= validdepth && validdepth <= currentdepth);

   /* temporarily move remaining bound changes from the queue into the conflict set */
   bdchginfos = (BDCHGINFO**)SCIPpqueueElems(conflict->binbdchgqueue);
   nbdchginfos = SCIPpqueueNElems(conflict->binbdchgqueue);
   debugMessage("adding %d variables from the queue as temporary conflict variables\n", nbdchginfos);
   for( i = 0; i < nbdchginfos; ++i )
   {
      CHECK_OKAY( conflictAddConflictVar(conflict, memhdr, set, stat, SCIPbdchginfoGetVar(bdchginfos[i]), FALSE, TRUE) );
   }

#ifndef NDEBUG
   {
      int v;
         
      for( v = 0; v < conflict->nconflictvars + conflict->ntmpconflictvars; ++v )
      {
         assert(conflict->conflictvars[v] != NULL);
         assert(conflict->conflictvars[v]->conflictsetcount == conflict->count);
      }
   }
#endif

#if 0 /*???????????????????? is this really true? */
   /* even if all branching variables up to a certain depth d are member of the conflict, we don't want to delete
    * those variables from the conflict clause and attach the conflict clause locally to the node in depth d,
    * because those branching decisions may be changed into inferences due to repropagation of nodes higher in the
    * tree
    */
   d = validdepth+1;
#else
   /* identify the depth, at which the conflict clause should be added:
    * - if the branching rule operates on variables only, and if all branching variables up to a certain
    *   depth level are member of the conflict, the conflict clause can only be violated in the subtree
    *   of the node at that depth, because in all other nodes, at least one of these branching variables
    *   takes a different value, such that the conflict clause is feasible
    * - if there is at least one branching variable in a node, we assume, that this branching was performed
    *   on variables, and that the siblings of this node are disjunct w.r.t. the branching variables' fixings
    * - the variables in the conflict set are labeled with the current conflict set counter
    * - we have to add the conflict clause at least in the valid depth of the initial conflict set,
    *   so we start searching at the first branching after this depth level, i.e. validdepth+1
    */
   for( d = validdepth+1; d <= currentdepth; ++d )
   {
      NODE* node;
      VAR* var;
      BOUNDCHG* boundchgs;
      BOUNDCHGTYPE boundchgtype;
      int nboundchgs;
      int nbranchingvars;
      int b;

      node = tree->path[d];
      assert(node != NULL);

      /* if node has no domain changes, the branching was not performed on variables */
      if( node->domchg == NULL )
         break;

      /* scan the bound changes of the current depth level: branchings are always first entries */
      boundchgs = node->domchg->domchgbound.boundchgs;
      nboundchgs = node->domchg->domchgbound.nboundchgs;
      nbranchingvars = 0;
      for( b = 0; b < nboundchgs; ++b )
      {
         var = boundchgs[b].var;
         boundchgtype = (BOUNDCHGTYPE)boundchgs[b].boundchgtype;

         debugMessage(" -> depth %d, bound change %d: <%s> %s %g (branching: %d, conflict set: %d)\n", 
            d, b, SCIPvarGetName(var), boundchgs[b].boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            boundchgs[b].newbound, (boundchgtype == SCIP_BOUNDCHGTYPE_BRANCHING),
            (var->conflictsetcount == conflict->count));

         if( boundchgtype != SCIP_BOUNDCHGTYPE_BRANCHING )
            break;
         if( var->conflictsetcount != conflict->count )
         {
            /* the branching variable is not member of the conflict set: abort at this depth level */
            nbranchingvars = 0;
            break;
         }
         nbranchingvars++;
      }

      /* break, if a branching variable not in the conflict set was found, or no branching variables exist,
       * which probably means, that this branching was not operating on variables
       */
      if( nbranchingvars == 0 )
         break;
   }
#endif

   /* now, d is the depth level of the first node, that has non-variable branching or a branching variable
    * not in the conflict set; this means, the siblings of the node may benefit from the conflict clause,
    * and the clause should be added to the node's parent, i.e. at depth level d-1
    */
   d--;
   assert(validdepth <= d && d <= currentdepth);
   debugMessage(" -> conflict with %d literals found at depth %d is active in depth %d and valid in depth %d\n", 
      conflict->nconflictvars + conflict->ntmpconflictvars, currentdepth, d, validdepth);
   
   /* if all branching variables are in the conflict set, the conflict clause is of no use */
   if( d < currentdepth )
   {
      VAR** conflictset;

      /* if the conflict clause is only valid locally, use only those variables from the conflict set, that are unfixed
       * at the local node
       */
      if( validdepth > 0 )
      {
         CHECK_OKAY( SCIPsetAllocBufferArray(set, &conflictset, conflict->nconflictvars + conflict->ntmpconflictvars) );
         *nliterals = 0;
         for( i = 0; i < conflict->nconflictvars + conflict->ntmpconflictvars; ++i )
         {
            if( conflict->conflictvardepths[i] > d )
            {
               conflictset[*nliterals] = conflict->conflictvars[i];
               (*nliterals)++;
            }
         }
      }
      else
      {
         conflictset = conflict->conflictvars;
         *nliterals = conflict->nconflictvars + conflict->ntmpconflictvars;
      }
      debugMessage(" -> final conflict clause has %d literals\n", *nliterals);

      /* if no conflict variables exist, the node and its sub tree in the conflict clause's depth can be 
       * cut off completely
       */
      if( *nliterals == 0 )
      {
         debugMessage("empty conflict clause in depth %d cuts off sub tree at depth %d\n", currentdepth, validdepth);
      
         SCIPnodeCutoff(tree->path[d], set, stat, tree);
         *success = TRUE;
      }
      else if( *nliterals <= maxsize )
      {
         RESULT result;
         int propdepth;
         int h;

         /* sort conflict handlers by priority */
         SCIPsetSortConflicthdlrs(set);
      
         /* call conflict handlers to create a conflict constraint */
         for( h = 0; h < set->nconflicthdlrs; ++h )
         {
            CHECK_OKAY( SCIPconflicthdlrExec(set->conflicthdlrs[h], set->scip, tree->path[d], 
                  conflictset, *nliterals, (validdepth > 0), *success, &result) );
            *success = *success || (result == SCIP_CONSADDED);
            debugMessage(" -> calling conflict handler <%s> (prio=%d) to create conflict clause with %d literals returned result %d\n",
               SCIPconflicthdlrGetName(set->conflicthdlrs[h]), SCIPconflicthdlrGetPriority(set->conflicthdlrs[h]),
               *nliterals, result);
         }

         if( *success )
         {
#ifdef DEBUG
            debugMessage(" -> conflict clause (valid:%d, active:%d):", validdepth, d);
            for( i = 0; i < *nliterals; ++i )
               printf(" <%s>", SCIPvarGetName(conflictset[i]));
            printf("\n");
#endif

            SCIPvbcFoundConflict(stat->vbc, stat, tree->path[currentdepth]);

            /* reactivate propagation on the node at depth level of the last but one fixed conflict variable,
             * such that the last fixed conflict variable can be deduced to the opposite value
             */
            if( set->conf_repropagate )
            {
               propdepth = conflictGetPropagateDepth(conflict);
               if( propdepth < INT_MAX )
               {
                  assert(0 <= propdepth && propdepth < tree->pathlen);
                  assert(tree->path[propdepth]->depth == propdepth);
                  SCIPnodePropagateAgain(tree->path[propdepth], set, stat, tree);
                  
                  debugMessage("marked node %p in depth %d to be repropagated due to conflict found in depth %d with %d literals\n", 
                     tree->path[propdepth], propdepth, currentdepth, *nliterals);
               }
            }
         }
      }

      /* free temporary buffer for the conflict set */
      if( validdepth > 0 )
         SCIPsetFreeBufferArray(set, &conflictset);
   }

   /* remove temporary conflict variables */
   conflictRemoveTemporaryConflictVars(conflict);

   return SCIP_OKAY;
}

/** tries to resolve given bound change
 *   - bound changes on non-binary variables are resolved in any case, but if this is done with a
 *     locally valid constraint, the valid depth level of the conflict set is updated
 *   - binary resolutions on local constraints are only be applied, if the constraint is valid at the 
 *     current minimal valid depth level, because this depth level is the topmost level to add the conflict
 *     clause to anyways
 */
static
RETCODE conflictResolveBound(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   BDCHGINFO*       bdchginfo,          /**< bound change to resolve */
   int*             validdepth,         /**< pointer to update the minimal depth level at which the conflict is valid */
   Bool*            resolved            /**< pointer to store whether the bound change was resolved */
   )
{
   VAR* actvar;
   CONS* infercons;
   PROP* inferprop;
   RESULT result;

   assert(conflict != NULL);
   assert(validdepth != NULL);
   assert(resolved != NULL);

   *resolved = FALSE;

   actvar = SCIPbdchginfoGetVar(bdchginfo);
   assert(actvar != NULL);
   assert(SCIPvarIsActive(actvar));

#ifdef DEBUG
   {
      int v;
         
      debugMessage("processing next conflicting bound (depth: %d, valid depth: %d, bdchgtype: %s, vartype: %d): [<%s> %s %g]\n",
         SCIPbdchginfoGetDepth(bdchginfo), *validdepth,
         SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
         : SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER ? "cons" : "prop",
         SCIPvarGetType(actvar), SCIPvarGetName(actvar), 
         SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(bdchginfo));
      debugMessage(" - conflict set             :");
      for( v = 0; v < conflict->nconflictvars; ++v )
         printf(" <%s>[%d]", SCIPvarGetName(conflict->conflictvars[v]), conflict->conflictvardepths[v]);
      printf("\n");
      debugMessage(" - candidate queue (non-bin):");
      for( v = 0; v < SCIPpqueueNElems(conflict->nonbinbdchgqueue); ++v )
      {
         BDCHGINFO* info = (BDCHGINFO*)(SCIPpqueueElems(conflict->nonbinbdchgqueue)[v]);
         printf(" [%d:<%s> %s %g]", SCIPbdchginfoGetDepth(info), SCIPvarGetName(SCIPbdchginfoGetVar(info)),
            SCIPbdchginfoGetBoundtype(info) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(info));
      }
      printf("\n");
      debugMessage(" - candidate queue (bin)    :");
      for( v = 0; v < SCIPpqueueNElems(conflict->binbdchgqueue); ++v )
      {
         BDCHGINFO* info = (BDCHGINFO*)(SCIPpqueueElems(conflict->binbdchgqueue)[v]);
         printf(" [%d:<%s> %s %g]", SCIPbdchginfoGetDepth(info), SCIPvarGetName(SCIPbdchginfoGetVar(info)),
            SCIPbdchginfoGetBoundtype(info) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(info));
      }
      printf("\n");
   }
#endif

   /* check, if the bound change can and should be resolved:
    *  - bound changes on non-binary variables have to be resolved in any case, but if this is done with a
    *    locally valid constraint, the valid depth level of the conflict set has to be updated
    *  - binary resolutions on local constraints should only be applied, if the constraint is valid at the 
    *    current minimal valid depth level (which is initialized with the valid depth level of the initial 
    *    conflict set), because this depth level is the topmost level to add the conflict clause to anyways
    */
   switch( SCIPbdchginfoGetChgtype(bdchginfo) )
   {
   case SCIP_BOUNDCHGTYPE_CONSINFER:
      infercons = SCIPbdchginfoGetInferCons(bdchginfo);
      assert(infercons != NULL);

      if( SCIPvarGetType(actvar) != SCIP_VARTYPE_BINARY
         || SCIPconsIsGlobal(infercons)
         || (SCIPconsIsActive(infercons) && SCIPconsGetActiveDepth(infercons) <= *validdepth) )
      {
         VAR* infervar;
         int inferinfo;
         BOUNDTYPE inferboundtype;
         BDCHGIDX* bdchgidx;

         /* resolve bound change by asking the constraint that infered the bound to put all bounds that were
          * the reasons for the conflicting bound change on the priority queue
          */
         infervar = SCIPbdchginfoGetInferVar(bdchginfo);
         inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);
         inferboundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         assert(infervar != NULL);
         assert(SCIPvarGetType(actvar) == SCIP_VARTYPE_BINARY
            || SCIPvarGetType(infervar) != SCIP_VARTYPE_BINARY);

         debugMessage("resolving bound <%s> %s %g [status:%d, type:%d, depth:%d, pos:%d]: <%s> %s %g [cons:<%s>(%s), info:%d]\n",
            SCIPvarGetName(actvar), 
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo),
            SCIPvarGetStatus(actvar), SCIPvarGetType(actvar),
            SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar), 
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPvarGetBdAtIndex(infervar, inferboundtype, bdchgidx, TRUE),
            SCIPconsGetName(infercons),
            SCIPconsIsGlobal(infercons) ? "global" : "local",
            inferinfo);
                  
         CHECK_OKAY( SCIPconsResolvePropagation(infercons, set, infervar, inferinfo, inferboundtype,
               bdchgidx, &result) );
         *resolved = (result == SCIP_SUCCESS);
         
         if( *resolved && SCIPconsIsLocal(infercons) )
            *validdepth = MAX(*validdepth, SCIPconsGetActiveDepth(infercons));
      }
      break;

   case SCIP_BOUNDCHGTYPE_PROPINFER:
      inferprop = SCIPbdchginfoGetInferProp(bdchginfo);
      if( inferprop != NULL )
      {
         VAR* infervar;
         int inferinfo;
         BOUNDTYPE inferboundtype;
         BDCHGIDX* bdchgidx;

         /* resolve bound change by asking the propagator that infered the bound to put all bounds that were
          * the reasons for the conflicting bound change on the priority queue
          */
         infervar = SCIPbdchginfoGetInferVar(bdchginfo);
         inferinfo = SCIPbdchginfoGetInferInfo(bdchginfo);
         inferboundtype = SCIPbdchginfoGetInferBoundtype(bdchginfo);
         bdchgidx = SCIPbdchginfoGetIdx(bdchginfo);
         assert(infervar != NULL);
         assert(SCIPvarGetType(actvar) == SCIP_VARTYPE_BINARY
            || SCIPvarGetType(infervar) != SCIP_VARTYPE_BINARY);

         debugMessage("resolving bound <%s> %s %g [status:%d, depth:%d, pos:%d]: <%s> %s %g [prop:<%s>, info:%d]\n",
            SCIPvarGetName(actvar), 
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo),
            SCIPvarGetStatus(actvar), SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar), 
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPvarGetBdAtIndex(infervar, inferboundtype, bdchgidx, TRUE),
            SCIPpropGetName(inferprop), inferinfo);
                  
         CHECK_OKAY( SCIPpropResolvePropagation(inferprop, set, infervar, inferinfo, inferboundtype,
               bdchgidx, &result) );
         *resolved = (result == SCIP_SUCCESS);
      }
      break;

   case SCIP_BOUNDCHGTYPE_BRANCHING:
      assert(!(*resolved));
      break;

   default:
      errorMessage("invalid bound change type <%d>\n", SCIPbdchginfoGetChgtype(bdchginfo));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** if only one conflicting bound change of the last depth level was used, and if this can be resolved,
 *  creates GRASP-like reconvergence clauses in the conflict graph up to the branching variable of the depth level
 */
static
RETCODE conflictCreateReconvergenceClauses(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   int              validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   int              maxsize,            /**< maximal size of conflict set */
   BDCHGINFO*       firstuip,           /**< first UIP of conflict graph */
   int*             nreconvclauses,     /**< pointer to store the number of generated reconvergence clauses */
   int*             nreconvliterals     /**< pointer to store the number of literals generated reconvergence clauses */
   )
{
   BDCHGINFO* uip;
   int firstuipdepth;

   assert(conflict != NULL);
   assert(firstuip != NULL);
   assert(nreconvclauses != NULL);
   assert(nreconvliterals != NULL);

   firstuipdepth = SCIPbdchginfoGetDepth(firstuip);

   /* for each succeeding UIP pair of the last depth level, create one reconvergence clause */
   uip = firstuip;
   while( uip != NULL && SCIPbdchginfoGetDepth(uip) == SCIPbdchginfoGetDepth(firstuip) )
   {
      BDCHGINFO* bdchginfo;
      BDCHGINFO* nextuip;
      VAR* var;
      int nresolutions;

      debugMessage("creating reconvergence clause for UIP <%s> in depth %d\n", 
         SCIPvarGetName(SCIPbdchginfoGetVar(uip)), SCIPbdchginfoGetDepth(uip));

      /* initialize conflict data */
      CHECK_OKAY( SCIPconflictInit(conflict) );

      /* put the variable of first UIP into the conflict set, using the literal that is currently fixed to TRUE */
      var = SCIPbdchginfoGetVar(uip);
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      CHECK_OKAY( SCIPvarNegate(var, memhdr, set, stat, &var) );
      CHECK_OKAY( conflictAddConflictVar(conflict, memhdr, set, stat, var, TRUE, FALSE) );
      
      /* put current UIP into priority queue */
      CHECK_OKAY( conflictAddBdchginfo(conflict, uip) );

      /* resolve the queue until the next UIP is reached */
      bdchginfo = conflictFirstCand(conflict);
      nextuip = NULL;
      nresolutions = 0;
      while( bdchginfo != NULL && conflict->nconflictvars < maxsize && validdepth < firstuipdepth )
      {
         BDCHGINFO* nextbdchginfo;
         int bdchgdepth;

         /* remove currently processed candidate and get next conflicting bound from the conflict candidate queue */
         assert(bdchginfo == conflictFirstCand(conflict));
         bdchginfo = conflictRemoveCand(conflict);
         nextbdchginfo = conflictFirstCand(conflict);
         bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
         assert(bdchginfo != NULL);
         assert(nextbdchginfo == NULL
            || SCIPbdchginfoGetDepth(bdchginfo) >= SCIPbdchginfoGetDepth(nextbdchginfo)
            || (SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) != SCIP_VARTYPE_BINARY
               && SCIPvarGetType(SCIPbdchginfoGetVar(nextbdchginfo)) == SCIP_VARTYPE_BINARY));
         assert(SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) != SCIP_VARTYPE_BINARY || bdchgdepth <= firstuipdepth);
         assert(SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) != SCIP_VARTYPE_BINARY
            || SCIPpqueueFirst(conflict->nonbinbdchgqueue) == NULL);
         assert(SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) != SCIP_VARTYPE_BINARY
            || SCIPpqueueNElems(conflict->nonbinbdchgqueue) == 0);

         /* bound changes that are higher in the tree than the valid depth of the conflict can be ignored;
          * multiple insertions of the same bound change can be ignored
          */
         if( bdchgdepth > validdepth && bdchginfo != nextbdchginfo )
         {
            VAR* actvar;
            Bool resolved;
            
            actvar = SCIPbdchginfoGetVar(bdchginfo);
            assert(actvar != NULL);
            assert(SCIPvarIsActive(actvar));
            
            /* check if we have to resolve the bound change in this depth level
             *  - a bound change on non-binary variables has to be resolved in any case,
             *  - the starting uip has to be resolved
             *  - a bound change on binary variables should be resolved, if it is in the fuip's depth level and not the
             *    next uip (i.e., if it is not the last bound change in the fuip's depth level)
             */
            resolved = FALSE;
            if( SCIPvarGetType(actvar) != SCIP_VARTYPE_BINARY
               || bdchginfo == uip
               || (bdchgdepth == firstuipdepth
                  && nextbdchginfo != NULL
                  && SCIPbdchginfoGetDepth(nextbdchginfo) == bdchgdepth) )
            {
               CHECK_OKAY( conflictResolveBound(conflict, memhdr, set, stat, tree, bdchginfo, &validdepth, &resolved) );
            }

            if( resolved )
               nresolutions++;
            else if( SCIPvarGetType(actvar) != SCIP_VARTYPE_BINARY )
            {
               /* non-binary variables cannot enter the conflict clause: we have to make the reconvergence clause
                * local, s.t. the unresolved non-binary bound change is active in the whole sub tree of the
                * reconvergence clause
                */
               assert(bdchgdepth >= validdepth);
               validdepth = bdchgdepth;
            }
            else if( bdchginfo != uip )
            {
               assert(SCIPsetIsEQ(set, SCIPbdchginfoGetNewbound(bdchginfo), SCIPvarGetLbLP(actvar)));
               assert(SCIPsetIsEQ(set, SCIPbdchginfoGetNewbound(bdchginfo), SCIPvarGetUbLP(actvar)));
               assert(SCIPvarGetLbLP(actvar) > 0.5 || SCIPvarGetUbLP(actvar) < 0.5);
               assert(conflict->nconflictvars >= 1); /* the starting UIP is already member of the conflict set */

               /* if this is the first variable of the conflict set besides the current starting UIP, it is the next
                * UIP (or the first unresolvable bound change)
                */
               if( bdchgdepth == firstuipdepth && conflict->nconflictvars == 1 )
               {
                  assert(nextuip == NULL);
                  nextuip = bdchginfo;
               }

               /* put variable into the conflict set, using the literal that is currently fixed to FALSE */
               CHECK_OKAY( conflictAddConflictVar(conflict, memhdr, set, stat, actvar, FALSE, FALSE) );
               assert(conflict->nconflictvars >= 2);
            }
            else
               assert(conflictFirstCand(conflict) == NULL); /* the starting UIP was not resolved */
         }

         /* get next conflicting bound from the conflict candidate queue (this need not to be nextbdchginfo, because
          * due to resolving the bound changes, a non-binary variable could be added to the queue which must be
          * resolved before nextbdchginfo
          */
         bdchginfo = conflictFirstCand(conflict);
      }
      assert(nextuip != uip);

      /* if only one propagation was resolved, the reconvergence clause is already member of the constraint set
       * (it is exactly the clause that produced the propagation)
       */
      if( nextuip != NULL && nresolutions >= 2
         && bdchginfo == NULL && conflict->nconflictvars < maxsize && validdepth < firstuipdepth )
      {
         int nlits;
         Bool success;

         assert(SCIPbdchginfoGetDepth(nextuip) == SCIPbdchginfoGetDepth(uip));
         assert(SCIPvarGetType(SCIPbdchginfoGetVar(nextuip)) == SCIP_VARTYPE_BINARY);

         debugMessage("creating reconvergence clause from UIP <%s> to UIP <%s> in depth %d with %d literals after %d resolutions\n", 
            SCIPvarGetName(SCIPbdchginfoGetVar(uip)), SCIPvarGetName(SCIPbdchginfoGetVar(nextuip)),
            SCIPbdchginfoGetDepth(uip), conflict->nconflictvars, nresolutions);

         /* call the conflict handlers to create a conflict clause */
         CHECK_OKAY( conflictAddClause(conflict, memhdr, set, stat, tree, validdepth, maxsize, &success, &nlits) );
         if( success )
         {
            (*nreconvclauses)++;
            (*nreconvliterals) += nlits;
         }
      }

      uip = nextuip;
   }

   return SCIP_OKAY;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set
 */
static
RETCODE conflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   int              validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   int              maxsize,            /**< maximal size of conflict set */
   Bool             mustresolve,        /**< should the conflict set only be used, if a resolution was applied? */
   int*             nclauses,           /**< pointer to store the number of generated conflict clauses */
   int*             nliterals,          /**< pointer to store the number of literals in generated conflict clauses */
   int*             nreconvclauses,     /**< pointer to store the number of generated reconvergence clauses */
   int*             nreconvliterals     /**< pointer to store the number of literals generated reconvergence clauses */
   )
{
   BDCHGINFO* bdchginfo;
   BDCHGINFO* firstuip;
   int currentdepth;
   int resolvedepth;
   int nresolutions;
   int lastclausenresolutions;
   int lastclauseresoldepth;

   assert(conflict != NULL);
   assert(conflict->nconflictvars >= 0);
   assert(set != NULL);
   assert(0 <= validdepth && validdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(nclauses != NULL);
   assert(nliterals != NULL);
   assert(nreconvclauses != NULL);
   assert(nreconvliterals != NULL);

   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);

   resolvedepth = ((set->conf_fuiplevels >= 0 && set->conf_fuiplevels < currentdepth)
      ? currentdepth - set->conf_fuiplevels : 0);
   assert(0 <= resolvedepth && resolvedepth <= currentdepth);

   debugMessage("analyzing conflict with %d+%d conflict candidates and starting conflict set of size %d in depth %d (maxsize=%d, resolvedepth=%d)\n",
      SCIPpqueueNElems(conflict->binbdchgqueue), SCIPpqueueNElems(conflict->nonbinbdchgqueue),
      conflict->nconflictvars, currentdepth, maxsize, resolvedepth);

   *nclauses = 0;
   *nliterals = 0;
   *nreconvclauses = 0;
   *nreconvliterals = 0;

   /* if the conflict set is only valid at the current node and not higher in the tree, no useful conflict clause
    * can be found
    */
   if( validdepth == currentdepth )
      return SCIP_OKAY;

   /* process all bound changes in the conflict candidate queue */
   nresolutions = 0;
   lastclausenresolutions = (mustresolve ? 0 : -1);
   lastclauseresoldepth = (mustresolve ? currentdepth : INT_MAX);
   bdchginfo = conflictFirstCand(conflict);
   firstuip = NULL;
   while( bdchginfo != NULL && conflict->nconflictvars < maxsize && validdepth < currentdepth )
   {
      BDCHGINFO* nextbdchginfo;
      int bdchgdepth;

      /* resolve next bound change in queue */
      bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
      assert(0 <= bdchgdepth && bdchgdepth <= currentdepth);
      assert(SCIPvarIsActive(SCIPbdchginfoGetVar(bdchginfo)));
      assert(bdchgdepth < tree->pathlen);
      assert(tree->path[bdchgdepth] != NULL);
      assert(tree->path[bdchgdepth]->domchg != NULL);
      assert(SCIPbdchginfoGetPos(bdchginfo) < tree->path[bdchgdepth]->domchg->domchgbound.nboundchgs);
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].var
         == SCIPbdchginfoGetVar(bdchginfo));
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].newbound
         == SCIPbdchginfoGetNewbound(bdchginfo));
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].boundtype
         == SCIPbdchginfoGetBoundtype(bdchginfo));

      /* create intermediate conflict clause */
      if( nresolutions > lastclausenresolutions
         && (set->conf_interclauses == -1 || *nclauses < set->conf_interclauses)
         && validdepth < currentdepth
         && SCIPpqueueNElems(conflict->nonbinbdchgqueue) == 0
         && bdchgdepth < lastclauseresoldepth )
      {
         int nlits;
         Bool success;
         
         /* call the conflict handlers to create a conflict clause */
         debugMessage("creating intermediate clause after %d resolutions up to depth %d (valid at depth %d): %d conflict vars, %d vars in queue\n",
            nresolutions, bdchgdepth, validdepth, conflict->nconflictvars, SCIPpqueueNElems(conflict->binbdchgqueue));
         CHECK_OKAY( conflictAddClause(conflict, memhdr, set, stat, tree, validdepth, maxsize, &success, &nlits) );
         lastclausenresolutions = nresolutions;
         lastclauseresoldepth = bdchgdepth;
         if( success )
         {
            (*nclauses)++;
            (*nliterals) += nlits;
         }
      }

      /* remove currently processed candidate and get next conflicting bound from the conflict candidate queue */
      assert(bdchginfo == conflictFirstCand(conflict));
      bdchginfo = conflictRemoveCand(conflict);
      nextbdchginfo = conflictFirstCand(conflict);
      assert(bdchginfo != NULL);
      assert(nextbdchginfo == NULL
         || SCIPbdchginfoGetDepth(bdchginfo) >= SCIPbdchginfoGetDepth(nextbdchginfo)
         || (SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) != SCIP_VARTYPE_BINARY
            && SCIPvarGetType(SCIPbdchginfoGetVar(nextbdchginfo)) == SCIP_VARTYPE_BINARY));

      /* we don't need to resolve bound changes that are already active in the valid depth of the current conflict set,
       * because the conflict clause can only be added locally at the valid depth, and all bound changes applied in this
       * depth or earlier can be removed from the conflict clause, since they are already applied in the constraint's 
       * subtree;
       * if the next bound change on the remaining queue is equal to the current bound change,
       * this is a multiple insertion in the conflict candidate queue and we can ignore the current
       * bound change
       */
      if( bdchgdepth > validdepth && bdchginfo != nextbdchginfo )
      {
         VAR* actvar;
         Bool resolved;

         actvar = SCIPbdchginfoGetVar(bdchginfo);
         assert(actvar != NULL);
         assert(SCIPvarIsActive(actvar));

         /* check if we want to resolve the bound change in this depth level
          *  - bound changes on non-binary variables have to be resolved in any case,
          *  - bound changes on binary variables should be resolved, if
          *     (i)  we must apply at least one resolution and didn't resolve a bound change yet, or
          *     (ii) their depth level is at least equal to the minimal resolving depth, and
          *          they are not the last remaining conflicting bound change in their depth level
          */
         resolved = FALSE;
         if( SCIPvarGetType(actvar) != SCIP_VARTYPE_BINARY
            || (mustresolve && nresolutions == 0)
            || (bdchgdepth >= resolvedepth
               && nextbdchginfo != NULL
               && SCIPbdchginfoGetDepth(nextbdchginfo) == bdchgdepth) )
         {
            CHECK_OKAY( conflictResolveBound(conflict, memhdr, set, stat, tree, bdchginfo, &validdepth, &resolved) );
         }

         if( resolved )
            nresolutions++;
         else if( SCIPvarGetType(actvar) != SCIP_VARTYPE_BINARY )
         {
            /* non-binary variables cannot enter the conflict clause: we have to make the conflict clause local, s.t.
             * the unresolved non-binary bound change is active in the whole sub tree of the conflict clause
             */
            assert(bdchgdepth >= validdepth);
            validdepth = bdchgdepth;

            debugMessage("couldn't resolve bound change on <%s> -> new valid depth: %d\n",
               SCIPvarGetName(actvar), validdepth);
         }
         else
         {
            assert(SCIPsetIsEQ(set, SCIPbdchginfoGetNewbound(bdchginfo), SCIPvarGetLbLP(actvar)));
            assert(SCIPsetIsEQ(set, SCIPbdchginfoGetNewbound(bdchginfo), SCIPvarGetUbLP(actvar)));
            assert(SCIPvarGetLbLP(actvar) > 0.5 || SCIPvarGetUbLP(actvar) < 0.5);

            /* if this is the first variable of the conflict set, it is the first UIP (or the first unresolvable bound
             * change)
             */
            if( conflict->nconflictvars == 0 )
               firstuip = bdchginfo;

            /* put variable into the conflict set, using the literal that is currently fixed to FALSE */
            CHECK_OKAY( conflictAddConflictVar(conflict, memhdr, set, stat, actvar, FALSE, FALSE) );
         }
      }

      /* get next conflicting bound from the conflict candidate queue (this need not to be nextbdchginfo, because
       * due to resolving the bound changes, a non-binary variable could be added to the queue which must be
       * resolved before nextbdchginfo
       */
      bdchginfo = conflictFirstCand(conflict);
   }

   /* check, if a valid conflict set was found */
   if( bdchginfo == NULL
      && nresolutions > lastclausenresolutions
      && validdepth < currentdepth
      && (!mustresolve || nresolutions > 0 || conflict->nconflictvars == 0) )
   {
      int nlits;
      Bool success;

      /* call the conflict handlers to create a conflict clause */
      CHECK_OKAY( conflictAddClause(conflict, memhdr, set, stat, tree, validdepth, maxsize, &success, &nlits) );
      if( success )
      {
         (*nclauses)++;
         (*nliterals) += nlits;
      }
   }

   /* produce reconvergence clauses defined by succeeding UIP's of the last depth level */
   if( set->conf_reconvclauses && firstuip != NULL && SCIPbdchginfoHasInferenceReason(firstuip) )
   {
      CHECK_OKAY( conflictCreateReconvergenceClauses(conflict, memhdr, set, stat, tree, validdepth, maxsize, firstuip,
            nreconvclauses, nreconvliterals) );
   }

   return SCIP_OKAY;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  updates statistics for propagation conflict analysis
 */
RETCODE SCIPconflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   int              validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   int maxsize;
   int nclauses;
   int nliterals;
   int nreconvclauses;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if propagation conflict analysis is enabled */
   if( !set->conf_useprop )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   debugMessage("analyzing conflict after infeasible propagation\n");

   /* calculate the maximal size of the conflict set */
   maxsize = (int)(set->conf_maxvarsfac * prob->nbinvars);
   maxsize = MAX(maxsize, set->conf_minmaxvars);
   if( maxsize < 2 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->propanalyzetime, set);

   conflict->npropcalls++;

   /* analyze the conflict set, and create a conflict constraint on success */
   CHECK_OKAY( conflictAnalyze(conflict, memhdr, set, stat, tree, validdepth, maxsize, TRUE,
         &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
   conflict->npropconfclauses += nclauses;
   conflict->npropconfliterals += nliterals;
   conflict->npropreconvclauses += nreconvclauses;
   conflict->npropreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nclauses > 0);

   /* stop timing */
   SCIPclockStop(conflict->propanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing propagation conflicts */
Real SCIPconflictGetPropTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->propanalyzetime);
}

/** gets number of calls to propagation conflict analysis */
Longint SCIPconflictGetNPropCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropcalls;
}

/** gets number of conflict clauses detected in propagation conflict analysis */
Longint SCIPconflictGetNPropConflictClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfclauses;
}

/** gets total number of literals in conflict clauses created in propagation conflict analysis */
Longint SCIPconflictGetNPropConflictLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfliterals;
}

/** gets number of reconvergence clauses detected in propagation conflict analysis */
Longint SCIPconflictGetNPropReconvergenceClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvclauses;
}

/** gets total number of literals in reconvergence clauses created in propagation conflict analysis */
Longint SCIPconflictGetNPropReconvergenceLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvliterals;
}




/*
 * Infeasible LP Conflict Analysis
 */

#if 0
/** Generate the alternative lp from the current lp */
/**@todo Correct hard coded 1000.0 */
static
RETCODE generateAltLP(
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   LPI*             alt_lpi,            /**< alternative lp on output */
   int*             alt_nrowvars,       /**< number of variables in the alt. LP on output belonging to original rows */
   int*             alt_ncolvars,       /**< number of variables in the alt. LP on output belonging to original columns */
   ROW**            alt_rowvars,        /**< array to store the original rows belonging to the variables in the alt. LP */
   COL**            alt_colvars,        /**< array to store the original cols belonging to the variables in the alt. LP */
   Bool*            alt_islower         /**< array to store whether the alt. LP variable belongs to the lower bound */
   )
{
   ROW** rows;
   COL** cols;
   ROW* row;
   COL* col;
   Real* vals;
   Real* rhs;
   int* inds;
   Real obj;
   Real lb;
   Real lpi_infinity;
   Real sum_rhs;    /* sum of the absolute nonzeros in the last row */
   int beg;
   int cnt;
   int nrows;
   int ncols;
   int i;
   int j;

   assert(lp != NULL);
   assert(alt_rowvars != NULL);
   assert(alt_colvars != NULL);
   assert(alt_nrowvars != NULL);
   assert(alt_ncolvars != NULL);
   assert(alt_islower != NULL);

   lpi_infinity = SCIPlpiInfinity(alt_lpi);
   sum_rhs = 0.0;
   *alt_nrowvars = 0;
   *alt_ncolvars = 0;

   /* get original LP data */
   ncols = SCIPlpGetNCols(lp);
   nrows = SCIPlpGetNRows(lp);
   cols = SCIPlpGetCols(lp);
   rows = SCIPlpGetRows(lp);

   /* allocate memory buffer for storing a single column of the alternative LP */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &inds, ncols + 1) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &vals, ncols + 1) );

   /* add empty rows to the alternative LP */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &rhs, ncols + 1) );
   clearMemoryArray(rhs, ncols + 1); /* clear rhs vector; the rhs of the last row will be modified later */
   clearMemoryArray(inds, ncols + 1); /* use inds array for storing "beg"; no matrix entries are used -> all zeros */
   CHECK_OKAY( SCIPlpiAddRows(alt_lpi, ncols+1, rhs, rhs, NULL, 0, inds, inds, vals) );
   SCIPsetFreeBufferArray(set, &rhs);

   /* set up main part */
   /**@todo Check for equality constraints and add only one column corr. to unbounded variable */
   for( i = 0; i < nrows; ++i )
   {
      COL** rowcols;
      Real* rowvals;
      Real rowlhs;
      Real rowrhs;
      int rowlen;

      row = rows[i];
      assert(row != NULL);

      rowlen = SCIProwGetNLPNonz(row);
      rowcols = SCIProwGetCols(row);
      rowvals = SCIProwGetVals(row);
      rowlhs = SCIProwGetLhs(row);
      rowrhs = SCIProwGetRhs(row);

      /* left hand side */
      if( !SCIPsetIsInfinity(set, -rowlhs) )
      {
         cnt = 0;
         for( j = 0; j < rowlen; ++j )
         {
            assert(SCIPcolIsInLP(rowcols[j]));
            inds[cnt] = SCIPcolGetLPPos(rowcols[j]);
            vals[cnt] = -rowvals[j];
            ++cnt;
         }
         if( !SCIPsetIsZero(set, rowlhs) )
         {
            inds[cnt] = ncols;
            vals[cnt] = -(rowlhs - SCIProwGetConstant(row));
            sum_rhs += REALABS(rowlhs);
            ++cnt;
         }
         if( row->local )
            obj = 1000.0;
         else
            obj = 0.0;
         lb = 0.0;
         beg = 0;
         CHECK_OKAY( SCIPlpiAddCols(alt_lpi, 1, &obj, &lb, &lpi_infinity, NULL, cnt, &beg, inds, vals) );
         alt_rowvars[*alt_nrowvars] = row;
         ++(*alt_nrowvars);
      }

      /* right hand side */
      if( !SCIPsetIsInfinity(set, rowrhs) )
      {
         cnt = 0;
         for( j = 0; j < rowlen; ++j )
         {
            assert(SCIPcolIsInLP(rowcols[j]));
            inds[cnt] = SCIPcolGetLPPos(rowcols[j]);
            vals[cnt] = rowvals[j];
            ++cnt;
         }
         if( !SCIPsetIsZero(set, rowrhs) )
         {
            inds[cnt] = ncols;
            vals[cnt] = rowrhs - SCIProwGetConstant(row);
            sum_rhs += REALABS(rowrhs);
            ++cnt;
         }
         if( row->local )
            obj = 1000.0;
         else
            obj = 0.0;
         lb = 0.0;
         beg = 0;
         CHECK_OKAY( SCIPlpiAddCols(alt_lpi, 1, &obj, &lb, &lpi_infinity, NULL, cnt, &beg, inds, vals) );
         alt_rowvars[*alt_nrowvars] = row;
         ++(*alt_nrowvars);
      }
   }
  
   /* add column for objective function */
   cnt = 0;
   for( j = 0; j < ncols; ++j )
   {
      col = cols[j];
      assert( col != NULL);
  
      if( !SCIPsetIsZero(set, col->obj) )
      {
         inds[cnt] = j;
         vals[cnt] = col->obj;
         ++cnt;
      }
   }
   if( !SCIPsetIsZero(set, lp->cutoffbound) )
   {
      inds[cnt] = ncols;
      vals[cnt] = lp->cutoffbound;
      ++cnt;
   }
   obj = 0.0;
   lb = 0.0;
   beg = 0;
   CHECK_OKAY( SCIPlpiAddCols(alt_lpi, 1, &obj, &lb, &lpi_infinity, NULL, cnt, &beg, inds, vals) );
   alt_rowvars[*alt_nrowvars] = NULL;
   ++(*alt_nrowvars);

   /* set up part for bounds */
   for( j = 0; j < ncols; ++j )
   {
      VAR* colvar;
      Real collb;
      Real colub;

      col = cols[j];
      assert( col != NULL);

      colvar = SCIPcolGetVar(col);
      collb = SCIPcolGetLb(col);
      colub = SCIPcolGetUb(col);

      /* lower bounds */
      if( !SCIPsetIsInfinity(set, -collb) )
      {
         cnt = 0;
         inds[cnt] = j;
         vals[cnt] = -1.0;
         ++cnt;

         if( !SCIPsetIsZero(set, collb) )
         {
            inds[cnt] = ncols;
            vals[cnt] = -collb;
            sum_rhs += REALABS(collb);
            ++cnt;
         }
         if( !SCIPsetIsEQ(set, SCIPvarGetLbGlobal(colvar), collb) )
         {
            if( SCIPvarGetType(colvar) == SCIP_VARTYPE_BINARY )
               obj = 1.0 + SCIPvarGetFixDepth(colvar);
            else
               obj = 1000.0;
         }
         else
            obj = 0.0;
         lb = 0.0;
         beg = 0;
         CHECK_OKAY( SCIPlpiAddCols(alt_lpi, 1, &obj, &lb, &lpi_infinity, NULL, cnt, &beg, inds, vals) );
         alt_colvars[*alt_ncolvars] = col;
         alt_islower[*alt_ncolvars] = TRUE;
         ++(*alt_ncolvars);
      }

      /* upper bounds */
      if( !SCIPsetIsInfinity(set, colub) )
      {
         cnt = 0;
         inds[cnt] = j;
         vals[cnt] = 1.0;
         ++cnt;

         if( !SCIPsetIsZero(set, colub) )
         {
            inds[cnt] = ncols;
            vals[cnt] = colub;
            sum_rhs += REALABS(colub);
            ++cnt;
         }
         if( !SCIPsetIsEQ(set, SCIPvarGetUbGlobal(colvar), colub) )
         {
            if( SCIPvarGetType(colvar) == SCIP_VARTYPE_BINARY )
               obj = 1.0 + SCIPvarGetFixDepth(colvar);
            else
               obj = 1000.0;
         }
         else
            obj = 0.0;
         lb = 0.0;
         beg = 0;
         CHECK_OKAY( SCIPlpiAddCols(alt_lpi, 1, &obj, &lb, &lpi_infinity, NULL, cnt, &beg, inds, vals) );
         alt_colvars[*alt_ncolvars] = col;
         alt_islower[*alt_ncolvars] = FALSE;
         ++(*alt_ncolvars);
      }
   }

   /* adjust rhs of last row */
   sum_rhs = - sum_rhs;
   CHECK_OKAY( SCIPlpiChgSides(alt_lpi, 1, &ncols, &sum_rhs, &sum_rhs) );

   /* free space */
   SCIPsetFreeBufferArray(set, &vals);
   SCIPsetFreeBufferArray(set, &inds);

   return SCIP_OKAY;
}

/** analyzes an infeasible LP trying to create a conflict set in the LP conflict analysis data structure */
static
RETCODE conflictAnalyzeLPAltpoly(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   int              maxsize,            /**< maximal size of conflict set */
   int*             iterations,         /**< pointer to store the total number of used LP iterations */
   Bool*            success             /**< pointer to store whether the conflict set is valid */
   )
{
   COL** alt_colvars;
   ROW** alt_rowvars;
   Bool* alt_islower;  
   int alt_ncolvars;
   int alt_nrowvars;
   int ncols;
   int nrows;
   int iter;

   assert(conflict != NULL);
   assert(success != NULL);
  
   *iterations = 0;
   *success = FALSE;

   warningMessage("this is not tested for infeasible strong branchings and infeasible dives\n");

   ncols = SCIPlpGetNCols(lp);
   nrows = SCIPlpGetNRows(lp);

   /* get temporary memory for alternative LP information */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &alt_rowvars, 2*nrows + 1) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &alt_colvars, 2*ncols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &alt_islower, 2*ncols) );

   /* create alternative LP */
   CHECK_OKAY( SCIPlpiClear(conflict->lpi) );
   CHECK_OKAY( generateAltLP(set, lp, conflict->lpi, 
         &alt_nrowvars, &alt_ncolvars, alt_rowvars, alt_colvars, alt_islower) );

   /* solve alternative LP with primal simplex, because the LP is primal feasible but not necessary dual feasible */
   CHECK_OKAY( SCIPlpiSolvePrimal(conflict->lpi) );

   /* update the LP iteration count */
   CHECK_OKAY( SCIPlpiGetIntpar(conflict->lpi, SCIP_LPPAR_LPITER, &iter) );
   (*iterations) += iter;

   /* check, if the LP is stable and we have a primal feasible solution */
   if( SCIPlpiIsStable(conflict->lpi) && SCIPlpiIsPrimalFeasible(conflict->lpi) )
   {
      COL* col;
      Real* psol;
      Bool valid;
      int i;
      int j;

      debugMessage("found IIS\n");
     
      /* get primal solution of alternative LP: the non-zero variables define an IIS in the original LP */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &psol, alt_nrowvars + alt_ncolvars) );
      CHECK_OKAY( SCIPlpiGetSol(conflict->lpi, NULL, psol, NULL, NULL, NULL) ); 

      valid = TRUE;

      /* scan the rows of the IIS for local cuts (the global upper bound is denoted by alt_rowvars[i] == NULL) */
      for( i = 0; i < alt_nrowvars; ++i )
      {
         if( alt_rowvars[i] != NULL && alt_rowvars[i]->local && SCIPsetIsPositive(set, psol[i]) )
         {
            debugMessage(" -> IIS includes local cut <%s>\n", SCIProwGetName(alt_rowvars[i]));
            valid = FALSE;
            break;
         }
      }

      if( valid )
      {
         /* initialize conflict data */
         CHECK_OKAY( SCIPconflictInit(conflict) );

         /* scan the column's bounds of the IIS:
          *  - if bound is an unchanged bound, we can ignore it, because it is valid in the root node
          *  - if bound is a changed bound of a binary variable, we have to include the variable in the conflict set
          *  - if bound is a changed bound of a non-binary variable, we must abort (conflict set can only include binaries)
          */
         for( j = 0; j < alt_ncolvars; ++j )
         {
            /* check, if the specific bound belongs to the IIS */
            if( SCIPsetIsPositive(set, psol[j + alt_nrowvars]) )
            {
               col = alt_colvars[j];
               /* check, if the specific bound differs from the root LP */
               if( (alt_islower[j] && !SCIPsetIsEQ(set, SCIPvarGetLbGlobal(col->var), col->lb))
                  || (!alt_islower[j] && !SCIPsetIsEQ(set, SCIPvarGetUbGlobal(col->var), col->ub)) )
               {
                  /* check, if the bound belongs to a binary variable */
                  if( SCIPvarGetType(col->var) == SCIP_VARTYPE_BINARY )
                  {
                     /* put binary variable into the conflict set */
                     CHECK_OKAY( SCIPconflictAddVar(conflict, col->var) );
                     debugMessage(" -> adding <%s> to IIS\n", SCIPvarGetName(col->var));
                  }
                  else
                  {
                     /* we cannot put a non-binary variable into the conflict set: we have to abort */
                     debugMessage(" -> IIS includes changed bound of non-binary variable <%s>\n", 
                        SCIPvarGetName(col->var));
                     valid = FALSE;
                     break;
                  }
               }      
            }
         }

         if( valid )
         {
            /* analyze the conflict set, and create a conflict constraint on success */
            CHECK_OKAY( conflictAnalyze(conflict, memhdr, set, stat, tree, 0, maxsize, FALSE, success) );
         }
      }
      
      /* free memory buffer for primal solution of alternative LP */
      SCIPsetFreeBufferArray(set, &psol);      
   }
   else
   {
      infoMessage(set->disp_verblevel, SCIP_VERBLEVEL_FULL, 
         "(node %lld) alternative LP is unstable or primal infeasible\n", stat->nnodes);
   }

   /* free memory buffers for alternative LP information */
   SCIPsetFreeBufferArray(set, &alt_islower);
   SCIPsetFreeBufferArray(set, &alt_colvars);
   SCIPsetFreeBufferArray(set, &alt_rowvars);
  
   return SCIP_OKAY;
}
#endif

/** returns, whether bound change information can be used in conflict analysis:
 *  - if variable is binary, or
 *  - if bound change information has valid inference data
 */
static
Bool bdchginfoIsUseable(
   BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return (SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) == SCIP_VARTYPE_BINARY
      || SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
      || (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_PROPINFER
         && SCIPbdchginfoGetInferProp(bdchginfo) != NULL));
}

/** ensures, that side change arrays can store at least num entries */
static
RETCODE ensureSidechgsSize(
   SET*             set,                /**< global SCIP settings */
   int**            sidechginds,        /**< pointer to side change index array */
   Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*             sidechgssize,       /**< pointer to size of side change arrays */
   int              num                 /**< minimal number of entries to be able to store in side change arrays */
   )
{
   assert(sidechginds != NULL);
   assert(sidechgoldlhss != NULL);
   assert(sidechgoldrhss != NULL);
   assert(sidechgnewlhss != NULL);
   assert(sidechgnewrhss != NULL);
   assert(sidechgssize != NULL);

   if( num > *sidechgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechginds, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechgoldlhss, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechgoldrhss, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechgnewlhss, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechgnewrhss, newsize) );
      *sidechgssize = newsize;
   }
   assert(num <= *sidechgssize);

   return SCIP_OKAY;
}

/** adds removal of row's side to side change arrays; finite sides are only replaced by near infinite sides, such
 *  that the row's sense in the LP solver is not changed
 */
static
RETCODE addSideRemoval(
   SET*             set,                /**< global SCIP settings */
   ROW*             row,                /**< LP row to change the sides for */
   Real             lpiinfinity,        /**< value treated as infinity in LP solver */
   int**            sidechginds,        /**< pointer to side change index array */
   Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*             sidechgssize,       /**< pointer to size of side change arrays */
   int*             nsidechgs           /**< pointer to number of used slots in side change arrays */
   )
{
   Real lhs;
   Real rhs;
   Real constant;

   assert(sidechginds != NULL);
   assert(sidechgoldlhss != NULL);
   assert(sidechgoldrhss != NULL);
   assert(sidechgnewlhss != NULL);
   assert(sidechgnewrhss != NULL);
   assert(sidechgssize != NULL);
   assert(nsidechgs != NULL);

   lhs = SCIProwGetLhs(row);
   rhs = SCIProwGetRhs(row);
   constant = SCIProwGetConstant(row);
   assert(!SCIPsetIsInfinity(set, -lhs) || !SCIPsetIsInfinity(set, rhs));

   /* get memory to store additional side change */
   CHECK_OKAY( ensureSidechgsSize(set, sidechginds, sidechgoldlhss, sidechgoldrhss, sidechgnewlhss, sidechgnewrhss,
         sidechgssize, (*nsidechgs)+1) );
   assert(*nsidechgs < *sidechgssize);
   assert(*sidechginds != NULL);
   assert(*sidechgoldlhss != NULL);
   assert(*sidechgoldrhss != NULL);
   assert(*sidechgnewlhss != NULL);
   assert(*sidechgnewrhss != NULL);
   
   /* store side change */
   (*sidechginds)[*nsidechgs] = SCIProwGetLPPos(row);
   if( SCIPsetIsInfinity(set, -lhs) )
   {
      (*sidechgoldlhss)[*nsidechgs] = -lpiinfinity;
      (*sidechgnewlhss)[*nsidechgs] = -lpiinfinity;
   }
   else
   {
      (*sidechgoldlhss)[*nsidechgs] = lhs - constant;
      (*sidechgnewlhss)[*nsidechgs] = -lpiinfinity/2;
   }
   if( SCIPsetIsInfinity(set, rhs) )
   {
      (*sidechgoldrhss)[*nsidechgs] = lpiinfinity;
      (*sidechgnewrhss)[*nsidechgs] = lpiinfinity;
   }
   else
   {
      (*sidechgoldrhss)[*nsidechgs] = rhs - constant;
      (*sidechgnewrhss)[*nsidechgs] = lpiinfinity/2;
   }
   (*nsidechgs)++;

   return SCIP_OKAY;
}

/** ensures, that bound change arrays can store at least num entries */
static
RETCODE ensureBdchgsSize(
   SET*             set,                /**< global SCIP settings */
   int**            sidechginds,        /**< pointer to side change index array */
   Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*             sidechgssize,       /**< pointer to size of side change arrays */
   int              num                 /**< minimal number of entries to be able to store in side change arrays */
   )
{
   assert(sidechginds != NULL);
   assert(sidechgoldlhss != NULL);
   assert(sidechgoldrhss != NULL);
   assert(sidechgnewlhss != NULL);
   assert(sidechgnewrhss != NULL);
   assert(sidechgssize != NULL);

   if( num > *sidechgssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechginds, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechgoldlhss, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechgoldrhss, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechgnewlhss, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, sidechgnewrhss, newsize) );
      *sidechgssize = newsize;
   }
   assert(num <= *sidechgssize);

   return SCIP_OKAY;
}

/** inserts variable's new bounds into bound change arrays */
static
RETCODE addBdchg(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to change the LP bounds for */
   Real             newlb,              /**< new lower bound */
   Real             newub,              /**< new upper bound */
   int**            bdchginds,          /**< pointer to bound change index array */
   Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array */
   Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array */
   Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array */
   Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array */
   int*             bdchgssize,         /**< pointer to size of bound change arrays */
   int*             nbdchgs             /**< pointer to number of used slots in bound change arrays */
   )
{
   assert(newlb <= newub);
   assert(bdchginds != NULL);
   assert(bdchgoldlbs != NULL);
   assert(bdchgoldubs != NULL);
   assert(bdchgnewlbs != NULL);
   assert(bdchgnewubs != NULL);
   assert(bdchgssize != NULL);
   assert(nbdchgs != NULL);
   assert(*nbdchgs <= *bdchgssize);

   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
   {
      COL* col;
      int c;
      
      col = SCIPvarGetCol(var);
      c = SCIPcolGetLPPos(col);
      if( c >= 0 )
      {
         CHECK_OKAY( ensureBdchgsSize(set, bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, 
               bdchgssize, (*nbdchgs)+1) );
         assert(*bdchginds != NULL);
         assert(*bdchgoldlbs != NULL);
         assert(*bdchgoldubs != NULL);
         assert(*bdchgnewlbs != NULL);
         assert(*bdchgnewubs != NULL);

         (*bdchginds)[*nbdchgs] = c;
         (*bdchgoldlbs)[*nbdchgs] = SCIPvarGetLbLP(var);
         (*bdchgoldubs)[*nbdchgs] = SCIPvarGetUbLP(var);
         (*bdchgnewlbs)[*nbdchgs] = newlb;
         (*bdchgnewubs)[*nbdchgs] = newub;
         (*nbdchgs)++;
      }
   }

   return SCIP_OKAY;
}

/** ensures, that candidate array can store at least num entries */
static
RETCODE ensureCandsSize(
   SET*             set,                /**< global SCIP settings */
   VAR***           cands,              /**< pointer to candidate array */
   Real**           candscores,         /**< pointer to candidate score array */
   Real**           newbounds,          /**< pointer to candidate new bounds array */
   Real**           proofactdeltas,     /**< pointer to candidate proof delta array */
   int*             candssize,          /**< pointer to size of array */
   int              num                 /**< minimal number of candidates to store in array */
   )
{
   assert(cands != NULL);
   assert(candssize != NULL);
   
   if( num > *candssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      CHECK_OKAY( SCIPsetReallocBufferArray(set, cands, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, candscores, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, newbounds, newsize) );
      CHECK_OKAY( SCIPsetReallocBufferArray(set, proofactdeltas, newsize) );
      *candssize = newsize;
   }
   assert(num <= *candssize);

   return SCIP_OKAY;
}

/** adds variable to candidate list, if the current best bound corresponding to the proof coefficient is local;
 *  returns the array position in the candidate list, where the new candidate was inserted, or -1 if the
 *  variable can relaxed to global bounds immediately without increasing the proof's activity;
 *  the candidates are sorted with respect to the following two criteria:
 *  - prefer bound changes that have been applied deeper in the tree, to get a more global conflict
 *  - prefer variables with small farkas coefficient to get rid of as many bound changes as possible
 */
static
RETCODE addCand(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to candidate array */
   int              lbchginfopos,       /**< positions of currently active lower bound change infos in variable's array */
   int              ubchginfopos,       /**< positions of currently active upper bound change infos in variable's array */
   Real             proofcoef,          /**< coefficient of variable in infeasibility/bound proof */
   Real             prooflhs,           /**< left hand side of infeasibility/bound proof */
   Real             proofact,           /**< activity of infeasibility/bound proof row */
   VAR***           cands,              /**< pointer to candidate array for undoing bound changes */
   Real**           candscores,         /**< pointer to candidate score array for undoing bound changes */
   Real**           newbounds,          /**< pointer to candidate new bounds array for undoing bound changes */
   Real**           proofactdeltas,     /**< pointer to proof activity increase array for undoing bound changes */
   int*             candssize,          /**< pointer to size of cands arrays */
   int*             ncands,             /**< pointer to count number of candidates in bound change list */
   int              firstcand,          /**< position of first unprocessed bound change candidate */
   int*             insertpos           /**< pointer to store insertion position, or -1 if not inserted */
   )
{
   Real oldbound;
   Real newbound;
   Real proofactdelta;
   Real score;
   int depth;
   int i;
   Bool useable;

   assert(set != NULL);
   assert(var != NULL);
   assert(-1 <= lbchginfopos && lbchginfopos <= var->nlbchginfos);
   assert(-1 <= ubchginfopos && ubchginfopos <= var->nubchginfos);
   assert(SCIPsetIsGT(set, prooflhs, proofact));
   assert(cands != NULL);
   assert(candscores != NULL);
   assert(newbounds != NULL);
   assert(proofactdeltas != NULL);
   assert(candssize != NULL);
   assert(ncands != NULL);
   assert(*ncands <= *candssize);
   assert(0 <= firstcand && firstcand <= *ncands);
   assert(insertpos != NULL);

   *insertpos = -1;

   /* in the infeasibility or dual bound proof, the variable's bound is chosen to maximize the proof's activity */
   if( SCIPsetIsPositive(set, proofcoef) )
   {
      /* if bound is global, nothing has to be done */
      if( ubchginfopos == -1 )
         return SCIP_OKAY;

      /* calculate the difference of current bound to the previous bound the variable was set to */
      if( ubchginfopos == var->nubchginfos )
      {
         /* current bound is the strong branching or diving bound: this can only happen on binary variables,
          * because the non-binary bound changes due to strong branching and diving were already removed in
          * the main loop of the infeasible LP analysis
          */
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         assert(SCIPvarGetUbLP(var) < 0.5);
         assert(SCIPvarGetUbLocal(var) > 0.5);
         oldbound = 0.0;
         newbound = 1.0;
         depth = 1000;
         useable = TRUE;
      }
      else
      {
         /* current bound is the result of a local bound change */
         oldbound = var->ubchginfos[ubchginfopos].newbound;
         newbound = var->ubchginfos[ubchginfopos].oldbound;
         depth = var->ubchginfos[ubchginfopos].bdchgidx.depth;
         useable = bdchginfoIsUseable(&var->ubchginfos[ubchginfopos]);
      }
   }
   else if( SCIPsetIsNegative(set, proofcoef) )
   {
      /* if bound is global, nothing has to be done */
      if( lbchginfopos == -1 )
         return SCIP_OKAY;

      /* calculate the difference of current bound to the previous bound the variable was set to */
      if( lbchginfopos == var->nlbchginfos )
      {
         /* current bound is the strong branching or diving bound: this can only happen on binary variables,
          * because the non-binary bound changes due to strong branching and diving were already removed in
          * the main loop of the infeasible LP analysis
          */
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         assert(SCIPvarGetLbLP(var) > 0.5);
         assert(SCIPvarGetLbLocal(var) < 0.5);
         oldbound = 1.0;
         newbound = 0.0;
         depth = 1000;
         useable = TRUE;
      }
      else
      {
         /* current bound is the result of a local bound change */
         oldbound = var->lbchginfos[lbchginfopos].newbound;
         newbound = var->lbchginfos[lbchginfopos].oldbound;
         depth = var->lbchginfos[lbchginfopos].bdchgidx.depth;
         useable = bdchginfoIsUseable(&var->lbchginfos[lbchginfopos]);
      }
   }
   else
      return SCIP_OKAY;

   /* calculate the increase in the proof's activity */
   proofactdelta = (newbound - oldbound)*proofcoef;
   assert(proofactdelta > 0.0);

   /* if the bound change is not useable in conflict analysis, we have to undo it */
   if( !useable )
      score = SCIPsetInfinity(set);
   else
   {
      /* calculate score for undoing the bound change */
      score = 1.0 - proofactdelta/(prooflhs - proofact);
      score = MAX(score, 0.0) + 1e-6;
      score *= depth;
      score = MIN(score, SCIPsetInfinity(set)/2);
   }
   
   /* get enough memory to store new candidate */
   CHECK_OKAY( ensureCandsSize(set, cands, candscores, newbounds, proofactdeltas, candssize, (*ncands)+1) );
   assert(*cands != NULL);
   assert(*candscores != NULL);
   assert(*newbounds != NULL);
   assert(*proofactdeltas != NULL);
    
   /* insert variable in candidate list without touching the already processed candidates */
   for( i = *ncands; i > firstcand && score > (*candscores)[i-1]; --i )
   {
      (*cands)[i] = (*cands)[i-1];
      (*candscores)[i] = (*candscores)[i-1];
      (*newbounds)[i] = (*newbounds)[i-1];
      (*proofactdeltas)[i] = (*proofactdeltas)[i-1];
   }
   (*cands)[i] = var;
   (*candscores)[i] = score;
   (*newbounds)[i] = newbound;
   (*proofactdeltas)[i] = proofactdelta;
   (*ncands)++;

   *insertpos = i;

   return SCIP_OKAY;
}

/** undos bound changes on variables, still leaving the given infeasibility proof valid */
static
RETCODE undoBdchgsProof(
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   Real*            proofcoefs,         /**< coefficients in infeasibility proof */
   Real             prooflhs,           /**< left hand side of proof */
   Real             proofact,           /**< current activity of proof */
   Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*             lbchginfoposs,      /**< positions of currently active lower bound change infos in variables' arrays */
   int*             ubchginfoposs,      /**< positions of currently active upper bound change infos in variables' arrays */
   int**            bdchginds,          /**< pointer to bound change index array, or NULL */
   Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*             bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*             nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again, or NULL */
   )
{
   VAR** vars;
   VAR** cands;
   Real* candscores;
   Real* newbounds;
   Real* proofactdeltas;
   int nvars;
   int ncands;
   int candssize;
   int pos;
   int v;
   int i;

   assert(prob != NULL);
   assert(proofcoefs != NULL);
   assert(SCIPsetIsFeasGT(set, prooflhs, proofact));
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);

   if( resolve != NULL )
      *resolve = FALSE;

   vars = prob->vars;
   nvars = prob->nvars;
   assert(nvars == 0 || vars != NULL);

   /* calculate the order in which the bound changes are tried to be undone, and relax all bounds if this doesn't
    * increase the proof's activity
    */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &cands, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &candscores, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &newbounds, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &proofactdeltas, nvars) );
   ncands = 0;
   candssize = nvars;
   for( v = 0; v < nvars; ++v )
   {
      /* ignore variables already relaxed to global bounds */
      if( lbchginfoposs[v] == -1 && ubchginfoposs[v] == -1 )
         continue;

      /* add variable to candidate list */
      CHECK_OKAY( addCand(set, vars[v], lbchginfoposs[v], ubchginfoposs[v], proofcoefs[v], prooflhs, proofact,
            &cands, &candscores, &newbounds, &proofactdeltas, &candssize, &ncands, 0, &pos) );
      assert(-1 <= pos && pos < ncands);
      if( pos == -1 )
      {
         assert(SCIPsetIsZero(set, proofcoefs[v])
            || (SCIPsetIsPositive(set, proofcoefs[v]) && ubchginfoposs[v] == -1)
            || (SCIPsetIsNegative(set, proofcoefs[v]) && lbchginfoposs[v] == -1));

         /* variable can be relaxed to global bounds */
         debugMessage(" -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g\n",
            SCIPvarGetName(vars[v]), curvarlbs[v], curvarubs[v], SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]),
            proofcoefs[v], prooflhs, proofact);
         curvarlbs[v] = SCIPvarGetLbGlobal(vars[v]);
         curvarubs[v] = SCIPvarGetUbGlobal(vars[v]);
         lbchginfoposs[v] = -1;
         ubchginfoposs[v] = -1;
         if( nbdchgs != NULL )
         {
            CHECK_OKAY( addBdchg(set, vars[v], curvarlbs[v], curvarubs[v], 
                  bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs) );
         }
      }
   }

   /* try to undo remaining local bound changes while still keeping the proof row violated:
    * bound changes can be undone, if prooflhs > proofact + proofactdelta;
    * afterwards, the current proof activity has to be updated
    */
   for( i = 0; i < ncands; ++i )
   {
      assert(proofactdeltas[i] > 0.0);

      if( SCIPsetIsGT(set, prooflhs, proofact + proofactdeltas[i]) )
      {
         v = SCIPvarGetProbindex(cands[i]);
         assert(0 <= v && v < nvars);
         assert(lbchginfoposs[v] >= 0 || ubchginfoposs[v] >= 0);

         debugMessage(" -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g + %g\n",
            SCIPvarGetName(cands[i]), curvarlbs[v], curvarubs[v],
            proofcoefs[v] > 0.0 ? curvarlbs[v] : newbounds[i],
            proofcoefs[v] > 0.0 ? newbounds[i] : curvarubs[v],
            proofcoefs[v], prooflhs, proofact, proofactdeltas[i]);

         assert((SCIPsetIsPositive(set, proofcoefs[v]) && SCIPsetIsGT(set, newbounds[i], curvarubs[v]))
            || (SCIPsetIsNegative(set, proofcoefs[v]) && SCIPsetIsLT(set, newbounds[i], curvarlbs[v])));
         assert((SCIPsetIsPositive(set, proofcoefs[v])
               && SCIPsetIsEQ(set, proofactdeltas[i], (newbounds[i] - curvarubs[v])*proofcoefs[v]))
            || (SCIPsetIsNegative(set, proofcoefs[v])
               && SCIPsetIsEQ(set, proofactdeltas[i], (newbounds[i] - curvarlbs[v])*proofcoefs[v])));
         assert(!SCIPsetIsZero(set, proofcoefs[v]));

         if( proofcoefs[v] > 0.0 )
         {
            assert(ubchginfoposs[v] >= 0);
            curvarubs[v] = newbounds[i];
            ubchginfoposs[v]--;
         }
         else
         {
            assert(lbchginfoposs[v] >= 0);
            curvarlbs[v] = newbounds[i];
            lbchginfoposs[v]--;
         }
         if( nbdchgs != NULL )
         {
            CHECK_OKAY( addBdchg(set, cands[i], curvarlbs[v], curvarubs[v], 
                  bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs) );
         }
         proofact += proofactdeltas[i];
         if( resolve != NULL )
            *resolve = TRUE;

         /* insert the new local bound of the variable into the candidate list */
         CHECK_OKAY( addCand(set, cands[i], lbchginfoposs[v], ubchginfoposs[v], proofcoefs[v], prooflhs, proofact,
               &cands, &candscores, &newbounds, &proofactdeltas, &candssize, &ncands, i+1, &pos) );
         assert(pos == -1 || (i < pos && pos < ncands));
      }
   }

   /* free the buffer for the sorted bound change candidates */
   SCIPsetFreeBufferArray(set, &proofactdeltas);
   SCIPsetFreeBufferArray(set, &newbounds);
   SCIPsetFreeBufferArray(set, &candscores);
   SCIPsetFreeBufferArray(set, &cands);

   return SCIP_OKAY;
}

/** analyzes an infeasible LP and undoes additional bound changes while staying infeasible */
static
RETCODE undoBdchgsDualfarkas(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< LP data */
   Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*             lbchginfoposs,      /**< positions of currently active lower bound change infos in variables' arrays */
   int*             ubchginfoposs,      /**< positions of currently active upper bound change infos in variables' arrays */
   int**            bdchginds,          /**< pointer to bound change index array, or NULL */
   Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*             bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*             nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again */
   )
{
   LPI* lpi;
   ROW** rows;
   VAR** vars;
   ROW* row;
   VAR* var;
   Real* dualfarkas;
   Real* farkascoefs;
   Real farkaslhs;
   Real farkasact;
   int nrows;
   int nvars;
   int r;
   int v;
   int i;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(valid != NULL);
   assert(resolve != NULL);

   debugMessage("undoing bound changes in infeasible LP: cutoff=%g, depth=%d\n", 
      lp->cutoffbound, SCIPgetDepth(set->scip));

   *valid = FALSE;
   *resolve = FALSE;

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);

   /* get LP rows and problem variables */
   rows = SCIPlpGetRows(lp);
   nrows = SCIPlpGetNRows(lp);
   vars = prob->vars;
   nvars = prob->nvars;
   assert(nrows == 0 || rows != NULL);
   assert(nrows == lp->nlpirows);

   /* allocate temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &dualfarkas, nrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &farkascoefs, nvars) );

   /* get dual farkas values of rows */
   CHECK_OKAY( SCIPlpiGetDualfarkas(lpi, dualfarkas) );

   /* calculate the farkas row */
   clearMemoryArray(farkascoefs, nvars);
   farkaslhs = 0.0;
   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore local rows and rows with farkas value 0.0 */
      if( !row->local && !SCIPsetIsZero(set, dualfarkas[r]) )
      {
#ifndef NDEBUG
         {
            Real lpilhs;
            Real lpirhs;

            CHECK_OKAY( SCIPlpiGetSides(lpi, r, r, &lpilhs, &lpirhs) );
            assert((SCIPsetIsInfinity(set, -lpilhs) && SCIPsetIsInfinity(set, -row->lhs))
               || SCIPsetIsEQ(set, lpilhs, row->lhs - row->constant));
            assert((SCIPsetIsInfinity(set, lpirhs) && SCIPsetIsInfinity(set, row->rhs))
               || SCIPsetIsEQ(set, lpirhs, row->rhs - row->constant));
         }
#endif

         /* add row coefficients to farkas row */
         for( i = 0; i < row->len; ++i )
         {
            v = SCIPvarGetProbindex(SCIPcolGetVar(row->cols[i]));
            assert(0 <= v && v < nvars);
            farkascoefs[v] += dualfarkas[r] * row->vals[i];
         }

         /* add row side to farkas row lhs: dualfarkas > 0 -> lhs, dualfarkas < 0 -> rhs */
         if( dualfarkas[r] > 0.0 )
         {
            assert(!SCIPsetIsInfinity(set, -row->lhs));
            farkaslhs += dualfarkas[r] * (row->lhs - row->constant);
         }
         else
         {
            assert(!SCIPsetIsInfinity(set, row->rhs));
            farkaslhs += dualfarkas[r] * (row->rhs - row->constant);
         }
      }
      else
      {
         /*debugMessage(" -> ignoring %s row <%s> with dual farkas value %.10f (lhs=%g, rhs=%g)\n",
           row->local ? "local" : "global", SCIProwGetName(row), dualfarkas[r], 
           row->lhs - row->constant, row->rhs - row->constant);*/
      }
   }

   /* calculate the current farkas activity, always using the best bound w.r.t. the farkas coefficient */
   farkasact = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(SCIPvarGetProbindex(var) == v);

      /* ignore coefs close to 0.0 */
      if( SCIPsetIsZero(set, farkascoefs[v]) )
      {
         /*debugMessage(" -> ignoring zero farkas coefficient <%s> [%g,%g]: %.10f\n",
           SCIPvarGetName(var), curvarlbs[v], curvarubs[v], farkascoefs[v]);*/
         farkascoefs[v] = 0.0;
      }
      else if( farkascoefs[v] > 0.0 )
      {
         assert((SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && SCIPcolGetLPPos(SCIPvarGetCol(var)) >= 0)
            || !SCIPsetIsPositive(set, curvarubs[v]));
         if( SCIPsetIsInfinity(set, curvarubs[v]) )
            goto TERMINATE;
         farkasact += farkascoefs[v] * curvarubs[v];
      }
      else
      {
         assert((SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && SCIPcolGetLPPos(SCIPvarGetCol(var)) >= 0)
            || !SCIPsetIsNegative(set, curvarlbs[v]));
         if( SCIPsetIsInfinity(set, -curvarlbs[v]) )
            goto TERMINATE;
         farkasact += farkascoefs[v] * curvarlbs[v];
      }
   }
   debugMessage(" -> farkaslhs=%g, farkasact=%g\n", farkaslhs, farkasact);

   /* check, if the farkas row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, farkaslhs, farkasact) )
   {
      /* undo bound changes while keeping the infeasibility proof valid */
      CHECK_OKAY( undoBdchgsProof(set, prob, farkascoefs, farkaslhs, farkasact, 
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
            bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs, resolve) );

      *valid = TRUE;
   }

 TERMINATE:

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &farkascoefs);
   SCIPsetFreeBufferArray(set, &dualfarkas);

   return SCIP_OKAY;
}

/** analyzes an LP exceeding the objective limit and undoes additional bound changes while staying beyond the
 *  objective limit
 */
static
RETCODE undoBdchgsDualsol(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< LP data */
   Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*             lbchginfoposs,      /**< positions of currently active lower bound change infos in variables' arrays */
   int*             ubchginfoposs,      /**< positions of currently active upper bound change infos in variables' arrays */
   int**            bdchginds,          /**< pointer to bound change index array, or NULL */
   Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*             bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*             nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again */
   )
{
   LPI* lpi;
   ROW** rows;
   VAR** vars;
   ROW* row;
   VAR* var;
   Real* primsols;
   Real* dualsols;
   Real* redcosts;
   Real* dualcoefs;
   Real* varredcosts;
   Real duallhs;
   Real dualact;
   Real duallhsdelta;
   Real dualactdelta;
   int nrows;
   int ncols;
   int nvars;
   int r;
   int v;
   int i;

   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(curvarlbs != NULL);
   assert(curvarubs != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(valid != NULL);
   assert(resolve != NULL);

   *valid = FALSE;
   *resolve = FALSE;

   debugMessage("undoing bound changes in LP exceeding cutoff: cutoff=%g, depth=%d\n", 
      lp->cutoffbound, SCIPgetDepth(set->scip));

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);

   /* get LP rows and problem variables */
   rows = SCIPlpGetRows(lp);
   nrows = SCIPlpGetNRows(lp);
   ncols = SCIPlpGetNCols(lp);
   vars = prob->vars;
   nvars = prob->nvars;
   assert(nrows == 0 || rows != NULL);
   assert(nrows == lp->nlpirows);

   /* get temporary memory */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &primsols, ncols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &dualsols, nrows) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &redcosts, ncols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &dualcoefs, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &varredcosts, nvars) );

   /* get solution from LPI */
   CHECK_OKAY( SCIPlpiGetSol(lpi, NULL, primsols, dualsols, NULL, redcosts) );

   /* Let y be the dual solution and r be the reduced cost vector. Let z be defined as
    *    z_i := y_i if i is a global row,
    *    z_i := 0   if i is a local row.
    * Define the set X := {x | lhs <= Ax <= rhs, lb <= x <= ub, c^Tx <= c*}, with c* being the current primal bound.
    * Then the following inequalities are valid for all x \in X:
    *                                 - c* <= -c^Tx
    *   <=>                     z^TAx - c* <= (z^TA - c^T) x
    *   <=>                     z^TAx - c* <= (y^TA - c^T - (y-z)^TA) x
    *   <=>                     z^TAx - c* <= (-r^T - (y-z)^TA) x         (dual feasibility of (y,r): y^TA + r^T == c^T)
    * Because lhs <= Ax <= rhs and lb <= x <= ub, the inequality can be relaxed to give
    *     min{z^Tq | lhs <= q <= rhs} - c* <= max{(-r^T - (y-z)^TA) x | lb <= x <= ub}, or X = {}.
    *
    * The resulting dual row is:  z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub},
    * where lhs, rhs, lb, and ub are selected in order to maximize the feasibility of the row.
    */

   clearMemoryArray(dualcoefs, nvars);

   /* use a slightly tighter cutoff bound, because solutions with equal objective value should also be declared
    * infeasible
    */
   duallhs = -(lp->cutoffbound - SCIPsetSumepsilon(set));
   dualact = 0.0;

   /* dual row: z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub}
    * process rows: add z^T{lhs,rhs} to the dual row's left hand side, and -(y-z)^TA to the dual row's coefficients
    */
   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore dual solution values of 0.0 (in this case: y_i == z_i == 0) */
      if( SCIPsetIsFeasZero(set, dualsols[r]) )
         continue;

      /* check dual feasibility */
      if( (SCIPsetIsInfinity(set, -row->lhs) && SCIPsetIsFeasPositive(set, dualsols[r]))
         || (SCIPsetIsInfinity(set, row->rhs) && SCIPsetIsFeasNegative(set, dualsols[r])) )
      {
         debugMessage(" -> infeasible dual solution %g in row <%s>: lhs=%g, rhs=%g\n",
            dualsols[r], SCIProwGetName(row), row->lhs, row->rhs);
         goto TERMINATE;
      }

      /* local rows add up to the dual row's coefficients (because z_i == 0 => -(y_i - z_i) == -y_i),
       * global rows add up to the dual row's left hand side (because z_i == y_i != 0)
       */
      if( row->local )
      {
         /* add -y_i A_i to coefficients of dual row */
         for( i = 0; i < row->len; ++i )
         {
            v = SCIPvarGetProbindex(SCIPcolGetVar(row->cols[i]));
            assert(0 <= v && v < nvars);
            dualcoefs[v] -= dualsols[r] * row->vals[i];
         }
         /*debugMessage(" -> local row <%s>: dual=%g\n", SCIProwGetName(row), dualsols[r]);*/
      }
      else
      {
         /* add minimal value to dual row's left hand side: z_i == y_i > 0 -> lhs, z_i == y_i < 0 -> rhs */
         if( dualsols[r] > 0.0 )
         {
            assert(!SCIPsetIsInfinity(set, -row->lhs));
            duallhs += dualsols[r] * (row->lhs - row->constant);
         }
         else
         {
            assert(!SCIPsetIsInfinity(set, row->rhs));
            duallhs += dualsols[r] * (row->rhs - row->constant);
         }
         /*debugMessage(" -> global row <%s>[%g,%g]: dual=%g -> duallhs=%g\n", 
           SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant, dualsols[r], duallhs);*/
      }
   }

   /* dual row: z^T{lhs,rhs} - c* <= (-r^T - (y-z)^TA){lb,ub}
    * process variables: subtract reduced costs from dual row's coefficients, and calculate current maximal dual
    *                    activity by multiplying the resultant coefficient with lb or ub
    */
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(SCIPvarGetProbindex(var) == v);

      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         /* reduced costs for loose variables are equal to the objective value */
         varredcosts[v] = SCIPvarGetObj(var);
      }
      else
      {
         COL* col;
         int c;

         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         col = SCIPvarGetCol(var);
         c = SCIPcolGetLPPos(col);
         assert(c == -1 || col == lp->cols[c]);
         assert(c == -1 || col == lp->lpicols[c]);

         /* get reduced costs from LPI, or calculate it manually if the column is not in current LP */
         varredcosts[v] = (c >= 0 ? redcosts[c] : SCIPcolCalcRedcost(col, dualsols));

         /* check dual feasibility */
         if( (SCIPsetIsGT(set, primsols[c], curvarlbs[v]) && SCIPsetIsFeasPositive(set, varredcosts[v]))
            || (SCIPsetIsLT(set, primsols[c], curvarubs[v]) && SCIPsetIsFeasNegative(set, varredcosts[v])) )
         {
            debugMessage(" -> infeasible reduced costs %g in var <%s>: lb=%g, ub=%g\n",
               varredcosts[v], SCIPvarGetName(var), curvarlbs[v], curvarubs[v]);
            goto TERMINATE;
         }
      }

      /* subtract reduced costs from dual row's coefficients */
      dualcoefs[v] -= varredcosts[v];

      /* add maximal value to dual row's activity: dualcoef > 0 -> ub, dualcoef < 0 -> lb */
      if( dualcoefs[v] > 0.0 )
      {
         if( SCIPsetIsInfinity(set, curvarubs[v]) )
            goto TERMINATE;
         dualact += dualcoefs[v] * curvarubs[v];
      }
      else
      {
         if( SCIPsetIsInfinity(set, -curvarlbs[v]) )
            goto TERMINATE;
         dualact += dualcoefs[v] * curvarlbs[v];
      }
      /*debugMessage(" -> col <%s> [%g,%g]: redcost=%g, coef=%g -> activity=%g\n", 
        SCIPvarGetName(var), curvarlbs[v], curvarubs[v], varredcosts[v], dualcoefs[v], dualact);*/
   }
   debugMessage(" -> final dual values: lhs=%g, act=%g\n", duallhs, dualact);

   /* check, if the dual row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, duallhs, dualact) )
   {
      /* undo bound changes while keeping the infeasibility proof valid */
      CHECK_OKAY( undoBdchgsProof(set, prob, dualcoefs, duallhs, dualact, 
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
            bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs, resolve) );

      *valid = TRUE;
   }

 TERMINATE:

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &varredcosts);
   SCIPsetFreeBufferArray(set, &dualcoefs);
   SCIPsetFreeBufferArray(set, &redcosts);
   SCIPsetFreeBufferArray(set, &dualsols);
   SCIPsetFreeBufferArray(set, &primsols);

   return SCIP_OKAY;
}

/** applies conflict analysis starting with given bound changes, that could not be undone during previous 
 *  infeasibility analysis
 */
static
RETCODE conflictAnalyzeRemainingBdchgs(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   int              maxsize,            /**< maximal size of conflict set */
   int*             lbchginfoposs,      /**< positions of currently active lower bound change infos in variables' arrays */
   int*             ubchginfoposs,      /**< positions of currently active upper bound change infos in variables' arrays */
   int*             nclauses,           /**< pointer to store the number of generated conflict clauses */
   int*             nliterals,          /**< pointer to store the number of literals in generated conflict clauses */
   int*             nreconvclauses,     /**< pointer to store the number of generated reconvergence clauses */
   int*             nreconvliterals     /**< pointer to store the number of literals generated reconvergence clauses */
   )
{
   VAR** vars;
   VAR* var;
   int nvars;
   int v;
   int nbdchgs;

   assert(prob != NULL);
   assert(lbchginfoposs != NULL);
   assert(ubchginfoposs != NULL);
   assert(nclauses != NULL);
   assert(nliterals != NULL);
   assert(nreconvclauses != NULL);
   assert(nreconvliterals != NULL);

   *nclauses = 0;
   *nliterals = 0;
   *nreconvclauses = 0;
   *nreconvliterals = 0;

   vars = prob->vars;
   nvars = prob->nvars;
   assert(nvars == 0 || vars != NULL);

   /* initialize conflict data */
   CHECK_OKAY( SCIPconflictInit(conflict) );
   
   /* add remaining bound changes to conflict queue; don't try to analyse the conflict, if the initial number of
    * remaining bound changes is too large (10 times larger than maximal conflict size)
    */
   debugMessage("initial conflict set after undoing bound changes:\n");
   nbdchgs = 0;
   for( v = 0; v < nvars && nbdchgs < 10*maxsize; ++v )
   {
      var = vars[v];
      assert(var != NULL);
      assert(-1 <= lbchginfoposs[v] && lbchginfoposs[v] <= var->nlbchginfos);
      assert(-1 <= ubchginfoposs[v] && ubchginfoposs[v] <= var->nubchginfos);

      if( lbchginfoposs[v] == var->nlbchginfos || ubchginfoposs[v] == var->nubchginfos )
      {
         /* the strong branching or diving bound stored in the column is responsible for the conflict:
          * it cannot be resolved and therefore has to be directly put into the conflict set;
          * the variable must be binary, because non-binary strong branching or diving bound changes
          * are already undone
          */
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         assert(SCIPvarGetLbLocal(var) < 0.5);
         assert(SCIPvarGetUbLocal(var) > 0.5);
         assert((lbchginfoposs[v] == var->nlbchginfos) != (ubchginfoposs[v] == var->nubchginfos));
         assert((lbchginfoposs[v] == var->nlbchginfos) == (SCIPvarGetLbLP(var) > 0.5));
         assert((ubchginfoposs[v] == var->nubchginfos) == (SCIPvarGetUbLP(var) < 0.5));

         /* put variable into the conflict set, using the literal that is currently fixed to FALSE */
         debugMessage("   fixed: <%s> == %g [status: %d, type: %d, dive/strong]\n",
            SCIPvarGetName(var), SCIPvarGetLbLP(var), SCIPvarGetStatus(var), SCIPvarGetType(var));
         CHECK_OKAY( conflictAddConflictVar(conflict, memhdr, set, stat, var, FALSE, FALSE) );
         nbdchgs++;
      }
      else
      {
         /* put remaining bound changes into conflict candidate queue */
         if( lbchginfoposs[v] >= 0 )
         {
            debugMessage("   queue: <%s> >= %g [status: %d, type: %d, depth: %d, pos:%d, chgtype: %d]\n", 
               SCIPvarGetName(var), SCIPbdchginfoGetNewbound(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPvarGetStatus(var), SCIPvarGetType(var),
               SCIPbdchginfoGetDepth(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPbdchginfoGetPos(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPbdchginfoGetChgtype(&var->lbchginfos[lbchginfoposs[v]]));
            CHECK_OKAY( conflictAddBdchginfo(conflict, &var->lbchginfos[lbchginfoposs[v]]) );
            nbdchgs++;
         }
         if( ubchginfoposs[v] >= 0 )
         {
            debugMessage("   queue: <%s> <= %g [status: %d, type: %d, depth: %d, pos:%d, chgtype: %d]\n", 
               SCIPvarGetName(var), SCIPbdchginfoGetNewbound(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPvarGetStatus(var), SCIPvarGetType(var),
               SCIPbdchginfoGetDepth(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPbdchginfoGetPos(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPbdchginfoGetChgtype(&var->ubchginfos[ubchginfoposs[v]]));
            CHECK_OKAY( conflictAddBdchginfo(conflict, &var->ubchginfos[ubchginfoposs[v]]) );
            nbdchgs++;
         }
      }
   }

   if( v == nvars )
   {
      /* analyze the conflict set, and create conflict constraints on success */
      CHECK_OKAY( conflictAnalyze(conflict, memhdr, set, stat, tree, 0, maxsize, FALSE, 
            nclauses, nliterals, nreconvclauses, nreconvliterals) );
   }

   return SCIP_OKAY;
}

/** actually performs analyzation of infeasible LP */
static
RETCODE conflictAnalyzeLP(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   Bool             diving,             /**< are we in strong branching or diving mode? */
   int*             iterations,         /**< pointer to store the total number of LP iterations used */
   int*             nclauses,           /**< pointer to store the number of generated conflict clauses */
   int*             nliterals,          /**< pointer to store the number of literals in generated conflict clauses */
   int*             nreconvclauses,     /**< pointer to store the number of generated reconvergence clauses */
   int*             nreconvliterals     /**< pointer to store the number of literals generated reconvergence clauses */
   )
{
   LPI* lpi;
   VAR** vars;
   Real* curvarlbs;
   Real* curvarubs;
   int* lbchginfoposs;
   int* ubchginfoposs;
   Real objval;
   int nvars;
   int maxsize;
   int v;
   int currentdepth;
   Bool valid;

   assert(set != NULL);
   assert(set->nactivepricers == 0); /* conflict analysis is not possible with unknown variables */
   assert(stat != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(iterations != NULL);
   assert(nclauses != NULL);
   assert(nliterals != NULL);
   assert(nreconvclauses != NULL);
   assert(nreconvliterals != NULL);

   *iterations = 0;
   *nclauses = 0;
   *nliterals = 0;
   *nreconvclauses = 0;
   *nreconvliterals = 0;

   /* conflict analysis only makes sense, if binary variables exist */
   if( prob->nbinvars == 0 )
      return SCIP_OKAY;

   /* calculate the maximal size of the conflict set */
   maxsize = (int)(set->conf_maxvarsfac * prob->nbinvars);
   maxsize = MAX(maxsize, set->conf_minmaxvars);
   if( maxsize < 2 )
      return SCIP_OKAY;

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);
   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsOptimal(lpi));
   if( SCIPlpiIsOptimal(lpi) )
   {
      CHECK_OKAY( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiuobjlim )
         return SCIP_OKAY;
   }

   currentdepth = SCIPtreeGetCurrentDepth(tree);

   debugMessage("analyzing conflict on infeasible LP (infeasible: %d, objlimexc: %d, optimal:%d) in depth %d\n",
      SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsOptimal(lpi), currentdepth);

   /* get active problem variables */
   vars = prob->vars;
   nvars = prob->nvars;

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information 
    * positions in variable's bound change information arrays
    */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &lbchginfoposs, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &ubchginfoposs, nvars) );

   /* the following algorithm is used to find a subset of changed bounds leading to an infeasible LP:
    * 1. undo all strong branching or diving bound changes on non-binary variables
    *    -> update lb/ubchginfoposs arrays
    *    -> store changes in bdchg and curvarlbs/ubs arrays
    *    -> apply changes to the LPI
    * 2. resolve LP
    * 3. call undoBdchgsDualfarkas() or undoBdchgsDualsol()
    *    -> update lb/ubchginfoposs arrays
    *    -> store additional changes in bdchg and curvarlbs/ubs arrays
    *    -> apply additional changes to the LPI
    * 4. if additional bound changes were undone, goto 2
    * 5. redo all bound changes in the LPI to restore the LPI to its original state
    * 6. analyze conflict
    *    -> put remaining changed bounds (see lb/ubchginfoposs arrays) into starting conflict set
    */
   
   /* get current bounds and current positions in lb/ubchginfos arrays of binary variables */
   valid = FALSE;
   for( v = 0; v < prob->nbinvars; ++v )
   {
      VAR* var;
      Real lb;
      Real ub;

      var = vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      
      curvarlbs[v] = SCIPvarGetLbLP(var);
      curvarubs[v] = SCIPvarGetUbLP(var);
      lbchginfoposs[v] = var->nlbchginfos-1;
      ubchginfoposs[v] = var->nubchginfos-1;

      /* check, if last bound changes were due to strong branching or diving */
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      if( !SCIPsetIsEQ(set, curvarlbs[v], lb) )
         lbchginfoposs[v] = var->nlbchginfos;
      if( !SCIPsetIsEQ(set, curvarubs[v], ub) )
         ubchginfoposs[v] = var->nubchginfos;
      assert(diving || lbchginfoposs[v] < var->nlbchginfos);
      assert(diving || ubchginfoposs[v] < var->nubchginfos);

      valid = valid || (v < prob->nbinvars && (curvarlbs[v] > 0.5 || curvarubs[v] < 0.5));
      /*debugMessage(" -> binary variable <%s> [%g,%g], bdchg=(%d,%d)\n", SCIPvarGetName(var), curvarlbs[v], curvarubs[v],
        lbchginfoposs[v], ubchginfoposs[v]);*/
   }

   /* conflict analysis only makes sense, if fixed binary variables exist */
   if( valid )
   {
      ROW** rows;
      COL** cols;
      int* sidechginds;
      Real* sidechgoldlhss;
      Real* sidechgoldrhss;
      Real* sidechgnewlhss;
      Real* sidechgnewrhss;
      int* bdchginds;
      Real* bdchgoldlbs;
      Real* bdchgoldubs;
      Real* bdchgnewlbs;
      Real* bdchgnewubs;
      Real lpiinfinity;
      Bool resolve;
      int sidechgssize;
      int nsidechgs;
      int bdchgssize;
      int nbdchgs;
      int lastnbdchgs;
      int nrows;
      int ncols;
      int iter;
      int nloops;
      int r;

      /* get infinity value of LP solver */
      lpiinfinity = SCIPlpiInfinity(lpi);

      /* get LP data */
      rows = SCIPlpGetRows(lp);
      nrows = SCIPlpGetNRows(lp);
      cols = SCIPlpGetCols(lp);
      ncols = SCIPlpGetNCols(lp);
      assert(nrows == 0 || rows != NULL);
      assert(ncols == 0 || cols != NULL);

      /* temporarily remove objective limit in LP solver */
      CHECK_OKAY( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, lpiinfinity) );

      /* get temporary memory for remembering side and bound changes on LPI rows and columns */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &sidechginds, nrows) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &sidechgoldlhss, nrows) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &sidechgoldrhss, nrows) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &sidechgnewlhss, nrows) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &sidechgnewrhss, nrows) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchginds, ncols) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchgoldlbs, ncols) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchgoldubs, ncols) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchgnewlbs, ncols) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchgnewubs, ncols) );
      sidechgssize = nrows;
      bdchgssize = ncols;
      nsidechgs = 0;
      nbdchgs = 0;
      lastnbdchgs = 0;

      /* remove all local rows by setting their sides to infinity;
       * finite sides are only changed to near infinity, such that the row's sense in the LP solver
       * is not affected (e.g. CPLEX cannot handle free rows)
       */
      for( r = 0 ; r < nrows; ++r )
      {
         assert(SCIProwGetLPPos(rows[r]) == r);

         if( SCIProwIsLocal(rows[r]) )
         {
            debugMessage(" -> removing local row <%s> [%g,%g]\n", 
               SCIProwGetName(rows[r]), SCIProwGetLhs(rows[r]), SCIProwGetRhs(rows[r]));
            CHECK_OKAY( addSideRemoval(set, rows[r], lpiinfinity, &sidechginds, &sidechgoldlhss, &sidechgoldrhss,
                  &sidechgnewlhss, &sidechgnewrhss, &sidechgssize, &nsidechgs) );
         }
      }

      /* apply changes of local rows to the LP solver */
      if( nsidechgs > 0 )
      {
         CHECK_OKAY( SCIPlpiChgSides(lpi, nsidechgs, sidechginds, sidechgnewlhss, sidechgnewrhss) );
      }

      /* get current bounds and current positions in lb/ubchginfos arrays of non-binary variables */
      for( v = prob->nbinvars; v < nvars; ++v )
      {
         VAR* var;

         var = vars[v];
         assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);

         curvarlbs[v] = SCIPvarGetLbLP(var);
         curvarubs[v] = SCIPvarGetUbLP(var);
         lbchginfoposs[v] = var->nlbchginfos-1;
         ubchginfoposs[v] = var->nubchginfos-1;
         
         /*debugMessage(" -> non-binary variable <%s> [%g,%g], bdchg=(%d,%d)\n", SCIPvarGetName(var), curvarlbs[v], curvarubs[v],
           lbchginfoposs[v], ubchginfoposs[v]);*/

         /* bound change has to be undone, if it is a strong branching/diving bound change */
         assert(diving || SCIPsetIsEQ(set, curvarlbs[v], SCIPvarGetLbLocal(var)));
         assert(diving || SCIPsetIsEQ(set, curvarubs[v], SCIPvarGetUbLocal(var)));
         if( diving )
         {
            Real lb;
            Real ub;
            Bool changed;

            lb = SCIPvarGetLbLocal(var);
            ub = SCIPvarGetUbLocal(var);
            changed = FALSE;
            if( !SCIPsetIsEQ(set, curvarlbs[v], lb) )
            {
               if( lbchginfoposs[v] >= 0 )
                  curvarlbs[v] = var->lbchginfos[lbchginfoposs[v]].newbound;
               else
                  curvarlbs[v] = SCIPvarGetLbGlobal(var);
               changed = TRUE;
            }
            if( !SCIPsetIsEQ(set, curvarubs[v], ub) )
            {
               if( ubchginfoposs[v] >= 0 )
                  curvarubs[v] = var->ubchginfos[ubchginfoposs[v]].newbound;
               else
                  curvarubs[v] = SCIPvarGetUbGlobal(var);
               changed = TRUE;
            }

            /* insert bound change in LPI bound change arrays */
            if( changed )
            {
               CHECK_OKAY( addBdchg(set, var, curvarlbs[v], curvarubs[v], 
                     &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs) );
               debugMessage(" -> undoing strong branching/diving bound change on non-binary variable <%s>: [%g,%g], bdchg=(%d,%d)\n", 
                  SCIPvarGetName(var), curvarlbs[v], curvarubs[v], lbchginfoposs[v], ubchginfoposs[v]);
            }
         }
      }

#if 0 /*?????????????????????????*/
      /* undo current depth's bound changes on non-binary variables without inference information */
      for( v = prob->nbinvars; v < nvars; ++v )
      {
         int lbpos;
         int ubpos;
         Bool changed;

         var = vars[v];
         assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
         assert(-1 <= lbchginfoposs[v] && lbchginfoposs[v] <= var->nlbchginfos);
         assert(-1 <= ubchginfoposs[v] && ubchginfoposs[v] <= var->nubchginfos);

         /* bound change has to be undone, if it is a strong branching/diving bound change, or
          * if no inference information is available in the corresponding bound change information data
          */

         /* calculate position of new active lower bound change */
         lbpos = lbchginfoposs[v];
         while( lbpos == var->nlbchginfos
            || (lbpos >= 0
               && SCIPbdchginfoGetDepth(&var->lbchginfos[lbpos]) == currentdepth
               && !bdchginfoIsUseable(&var->lbchginfos[lbpos])) )
         {
            lbpos--;
         }
         assert(-1 <= lbpos && lbpos < var->nlbchginfos);
            
         /* calculate position of new active upper bound change */
         ubpos = ubchginfoposs[v];
         while( ubpos == var->nubchginfos
            || (ubpos >= 0
               && SCIPbdchginfoGetDepth(&var->ubchginfos[ubpos]) == currentdepth
               && !bdchginfoIsUseable(&var->ubchginfos[ubpos])) )
         {
            ubpos--;
         }
         assert(-1 <= ubpos && ubpos < var->nubchginfos);

         /* change bounds and bound change positions if necessary */
         changed = FALSE;
         if( lbpos < lbchginfoposs[v] )
         {
            lbchginfoposs[v] = lbpos;
            if( lbpos >= 0 )
               curvarlbs[v] = var->lbchginfos[lbpos].newbound;
            else
               curvarlbs[v] = SCIPvarGetLbGlobal(var);
            changed = TRUE;
         }
         if( ubpos < ubchginfoposs[v] )
         {
            ubchginfoposs[v] = ubpos;
            if( ubpos >= 0 )
               curvarubs[v] = var->ubchginfos[ubpos].newbound;
            else
               curvarubs[v] = SCIPvarGetUbGlobal(var);
            changed = TRUE;
         }

         /* insert bound change in LPI bound change arrays */
         if( changed )
         {
            CHECK_OKAY( addBdchg(set, var, curvarlbs[v], curvarubs[v], 
                  &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs) );
            debugMessage(" -> undoing bound change on non-binary variable <%s> without inference info: [%g,%g], bdchg=(%d,%d)\n", 
               SCIPvarGetName(var), curvarlbs[v], curvarubs[v], lbchginfoposs[v], ubchginfoposs[v]);
         }
      }
#endif

      /* undo as many additional bound changes as possible */
      nloops = 0;
      do
      {
         nloops++;
         resolve = FALSE;

         debugMessage("infeasible LP conflict analysis loop %d (changed col bounds: %d)\n", nloops, nbdchgs);

#if 1 /*?????????????????????*/
         /* undo bound changes on non-binary variables in current depth without inference information */
         if( !diving )
         {
            for( v = prob->nbinvars; v < nvars; ++v )
            {
               VAR* var;
               int lbpos;
               int ubpos;
               Bool changed;

               var = vars[v];
               assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
               assert(-1 <= lbchginfoposs[v] && lbchginfoposs[v] < var->nlbchginfos);
               assert(-1 <= ubchginfoposs[v] && ubchginfoposs[v] < var->nubchginfos);

               /* bound change has to be undone, if it is a strong branching/diving bound change, or
                * if no inference information is available in the corresponding bound change information data
                */

               /* calculate position of new active lower bound change */
               lbpos = lbchginfoposs[v];
               while( lbpos >= 0
                  && SCIPbdchginfoGetDepth(&var->lbchginfos[lbpos]) == currentdepth
                  && !bdchginfoIsUseable(&var->lbchginfos[lbpos]) )
               {
                  lbpos--;
               }
               assert(-1 <= lbpos && lbpos < var->nlbchginfos);
            
               /* calculate position of new active upper bound change */
               ubpos = ubchginfoposs[v];
               while( ubpos >= 0
                  && SCIPbdchginfoGetDepth(&var->ubchginfos[ubpos]) == currentdepth
                  && !bdchginfoIsUseable(&var->ubchginfos[ubpos]) )
               {
                  ubpos--;
               }
               assert(-1 <= ubpos && ubpos < var->nubchginfos);

               /* change bounds and bound change positions if necessary */
               changed = FALSE;
               if( lbpos < lbchginfoposs[v] )
               {
                  lbchginfoposs[v] = lbpos;
                  if( lbpos >= 0 )
                     curvarlbs[v] = var->lbchginfos[lbpos].newbound;
                  else
                     curvarlbs[v] = SCIPvarGetLbGlobal(var);
                  changed = TRUE;
               }
               if( ubpos < ubchginfoposs[v] )
               {
                  ubchginfoposs[v] = ubpos;
                  if( ubpos >= 0 )
                     curvarubs[v] = var->ubchginfos[ubpos].newbound;
                  else
                     curvarubs[v] = SCIPvarGetUbGlobal(var);
                  changed = TRUE;
               }

               /* insert bound change in LPI bound change arrays */
               if( changed )
               {
                  CHECK_OKAY( addBdchg(set, var, curvarlbs[v], curvarubs[v], 
                        &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs) );
                  debugMessage(" -> undoing bound change on non-binary variable <%s> without inference info: [%g,%g], bdchg=(%d,%d)\n", 
                     SCIPvarGetName(var), curvarlbs[v], curvarubs[v], lbchginfoposs[v], ubchginfoposs[v]);
               }
            }
         }
#endif

         /* apply bound changes to the LP solver */
         assert(nbdchgs >= lastnbdchgs);
         if( nbdchgs > lastnbdchgs )
         {
            debugMessage(" -> applying %d bound changes to the LP solver (total: %d)\n", nbdchgs - lastnbdchgs, nbdchgs);
            CHECK_OKAY( SCIPlpiChgBounds(lpi, nbdchgs - lastnbdchgs,
                  &bdchginds[lastnbdchgs], &bdchgnewlbs[lastnbdchgs], &bdchgnewubs[lastnbdchgs]) );
            lastnbdchgs = nbdchgs;
         }
         
         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve LP */
         CHECK_OKAY( SCIPlpiSolveDual(lpi) );

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* evaluate result */
         if( SCIPlpiIsOptimal(lpi) )
         {
            CHECK_OKAY( SCIPlpiGetObjval(lpi, &objval) );
            valid = (objval >= lp->lpiuobjlim);
         }
         else
            valid = SCIPlpiIsPrimalInfeasible(lpi);

         /* count number of LP iterations */
         CHECK_OKAY( SCIPlpiGetIterations(lpi, &iter) );
         (*iterations) += iter;
         stat->nconflictlps++;
         stat->nconflictlpiterations += iter;
         debugMessage(" -> resolved LP in %d iterations (total: %lld) (valid: %d, infeasible:%d)\n",
            iter, stat->nconflictlpiterations, valid, SCIPlpiIsPrimalInfeasible(lpi));

         if( valid )
         {
            /* undo additional bound changes */
            if( SCIPlpiIsPrimalInfeasible(lpi) )
            {
               CHECK_OKAY( undoBdchgsDualfarkas(set, stat, prob, lp, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
                     &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
                     &valid, &resolve) );
            }
            else
            {
               assert(SCIPlpiIsOptimal(lpi));
               CHECK_OKAY( undoBdchgsDualsol(set, stat, prob, lp, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
                     &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
                     &valid, &resolve) );
            }
         }
         assert(!resolve || valid);
         assert(!resolve || nbdchgs > lastnbdchgs);
         debugMessage(" -> finished infeasible LP conflict analysis loop %d (iter: %d, nbdchgs: %d)\n",
            nloops, iter, nbdchgs - lastnbdchgs);
      }
      while( resolve && nloops < set->conf_maxlploops );
      debugMessage("finished undoing bound changes after %d loops (valid=%d, nbdchgs: %d)\n",
         nloops, valid, nbdchgs);

      /* reset variables to local bounds */
      if( nbdchgs > 0 )
      {
         CHECK_OKAY( SCIPlpiChgBounds(lpi, nbdchgs, bdchginds, bdchgoldlbs, bdchgoldubs) );
      }

      /* reset changes of local rows */
      if( nsidechgs > 0 )
      {
         CHECK_OKAY( SCIPlpiChgSides(lpi, nsidechgs, sidechginds, sidechgoldlhss, sidechgoldrhss) );
      }

      /* reset objective limit in LP solver */
      CHECK_OKAY( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, lp->lpiuobjlim) );

      /* analyze the conflict starting with remaining bound changes */
      if( valid )
      {
         CHECK_OKAY( conflictAnalyzeRemainingBdchgs(conflict, memhdr, set, stat, prob, tree,
               maxsize, lbchginfoposs, ubchginfoposs, nclauses, nliterals, nreconvclauses, nreconvliterals) );
      }

      /* free temporary memory */
      SCIPsetFreeBufferArray(set, &sidechgnewrhss);
      SCIPsetFreeBufferArray(set, &sidechgnewlhss);
      SCIPsetFreeBufferArray(set, &sidechgoldrhss);
      SCIPsetFreeBufferArray(set, &sidechgoldlhss);
      SCIPsetFreeBufferArray(set, &sidechginds);
      SCIPsetFreeBufferArray(set, &bdchgnewubs);
      SCIPsetFreeBufferArray(set, &bdchgnewlbs);
      SCIPsetFreeBufferArray(set, &bdchgoldubs);
      SCIPsetFreeBufferArray(set, &bdchgoldlbs);
      SCIPsetFreeBufferArray(set, &bdchginds);
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &ubchginfoposs);
   SCIPsetFreeBufferArray(set, &lbchginfoposs);
   SCIPsetFreeBufferArray(set, &curvarubs);
   SCIPsetFreeBufferArray(set, &curvarlbs);

   return SCIP_OKAY;
}

/** analyzes an infeasible LP to find out the bound changes on binary variables that were responsible for the
 *  infeasibility;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible LP conflict analysis
 */
RETCODE SCIPconflictAnalyzeLP(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   int iterations;
   int nclauses;
   int nliterals;
   int nreconvclauses;
   int nreconvliterals;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->conf_uselp )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;
   
   /* LP conflict analysis is only possible, if all variables are known */
   if( set->nactivepricers > 0 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->lpanalyzetime, set);
   conflict->nlpcalls++;

   /* perform conflict analysis */
   CHECK_OKAY( conflictAnalyzeLP(conflict, memhdr, set, stat, prob, tree, lp, SCIPlpDiving(lp),
         &iterations, &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
   conflict->nlpiterations += iterations;
   conflict->nlpconfclauses += nclauses;
   conflict->nlpconfliterals += nliterals;
   conflict->nlpreconvclauses += nreconvclauses;
   conflict->nlpreconvliterals += nreconvliterals;
   if( success != NULL )
      *success = (nclauses > 0);

   /* stop timing */
   SCIPclockStop(conflict->lpanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing infeasible LP conflicts */
Real SCIPconflictGetLPTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->lpanalyzetime);
}

/** gets number of calls to infeasible LP conflict analysis */
Longint SCIPconflictGetNLPCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpcalls;
}

/** gets number of conflict clauses detected in infeasible LP conflict analysis */
Longint SCIPconflictGetNLPConflictClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpconfclauses;
}

/** gets total number of literals in conflict clauses created in infeasible LP conflict analysis */
Longint SCIPconflictGetNLPConflictLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpconfliterals;
}

/** gets number of reconvergence clauses detected in infeasible LP conflict analysis */
Longint SCIPconflictGetNLPReconvergenceClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpreconvclauses;
}

/** gets total number of literals in reconvergence clauses created in infeasible LP conflict analysis */
Longint SCIPconflictGetNLPReconvergenceLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpreconvliterals;
}

/** gets number of LP iterations in infeasible LP conflict analysis */
Longint SCIPconflictGetNLPIterations(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpiterations;
}




/*
 * infeasible strong branching conflict analysis
 */

/** analyses infeasible strong branching sub problems for conflicts */
RETCODE SCIPconflictAnalyzeStrongbranch(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   COL*             col,                /**< LP column with at least one infeasible strong branching subproblem */
   Bool*            downconflict,       /**< pointer to store whether a conflict clause was created for an infeasible
                                         *   downwards branch, or NULL */
   Bool*            upconflict          /**< pointer to store whether a conflict clause was created for an infeasible
                                         *   upwards branch, or NULL */
   )
{
   VAR* var;
   int* cstat;
   int* rstat;
   Real oldlb;
   Real oldub;
   Real newlb;
   Real newub;
   int iter;
   int nclauses;
   int nliterals;
   int nreconvclauses;
   int nreconvliterals;

   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(col != NULL);
   assert((SCIPsetIsGE(set, col->strongbranchdown, lp->cutoffbound) && SCIPsetCeil(set, col->primsol-1.0) >= col->lb - 0.5)
      || (SCIPsetIsGE(set, col->strongbranchup, lp->cutoffbound) && SCIPsetFloor(set, col->primsol+1.0) <= col->ub + 0.5));

   if( downconflict != NULL )
      *downconflict = FALSE;
   if( upconflict != NULL )
      *upconflict = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->conf_usesb )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* LP conflict analysis is only possible, if all variables are known */
   if( set->nactivepricers > 0 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->sbanalyzetime, set);
   conflict->nsbcalls++;

   /* get temporary memory for storing current LP basis */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &cstat, lp->nlpicols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &rstat, lp->nlpirows) );

   /* get current LP basis */
   CHECK_OKAY( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

   /* remember old bounds */
   oldlb = col->lb;
   oldub = col->ub;

   var = SCIPcolGetVar(col);

   /* is down branch infeasible? */
   if( col->strongbranchdown >= lp->cutoffbound )
   {
      newub = SCIPsetCeil(set, col->primsol-1.0);
      if( newub >= col->lb - 0.5 )
      {
         debugMessage("analyzing conflict on infeasible downwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPtreeGetCurrentDepth(tree));

         /* change the upper bound */
         col->ub = newub;
         CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         CHECK_OKAY( SCIPlpiSolveDual(lp->lpi) );
            
         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* count number of LP iterations */
         CHECK_OKAY( SCIPlpiGetIterations(lp->lpi, &iter) );
         stat->nconflictlps++;
         stat->nconflictlpiterations += iter;
         conflict->nsbiterations += iter;
         debugMessage(" -> resolved downwards strong branching LP in %d iterations\n", iter);

         /* perform conflict analysis on infeasible LP */
         CHECK_OKAY( conflictAnalyzeLP(conflict, memhdr, set, stat, prob, tree, lp, TRUE, 
               &iter, &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
         conflict->nsbiterations += iter;
         conflict->nsbconfclauses += nclauses;
         conflict->nsbconfliterals += nliterals;
         conflict->nsbreconvclauses += nreconvclauses;
         conflict->nsbreconvliterals += nreconvliterals;
         if( downconflict != NULL )
            *downconflict = (nclauses > 0);

         /* reset the upper bound */
         col->ub = oldub;
         CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );
            
         /* reset LP basis */
         CHECK_OKAY( SCIPlpiSetBase(lp->lpi, cstat, rstat) );
      }
   }

   /* is up branch infeasible? */
   if( col->strongbranchup >= lp->cutoffbound )
   {
      newlb = SCIPsetFloor(set, col->primsol+1.0);
      if( newlb <= col->ub + 0.5 )
      {
         debugMessage("analyzing conflict on infeasible upwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPtreeGetCurrentDepth(tree));

         /* change the lower bound */
         col->lb = newlb;
         CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );
            
         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         CHECK_OKAY( SCIPlpiSolveDual(lp->lpi) );
            
         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* count number of LP iterations */
         CHECK_OKAY( SCIPlpiGetIterations(lp->lpi, &iter) );
         stat->nconflictlps++;
         stat->nconflictlpiterations += iter;
         conflict->nsbiterations += iter;
         debugMessage(" -> resolved upwards strong branching LP in %d iterations\n", iter);

         /* perform conflict analysis on infeasible LP */
         CHECK_OKAY( conflictAnalyzeLP(conflict, memhdr, set, stat, prob, tree, lp, TRUE, 
               &iter, &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
         conflict->nsbiterations += iter;
         conflict->nsbconfclauses += nclauses;
         conflict->nsbconfliterals += nliterals;
         conflict->nsbreconvclauses += nreconvclauses;
         conflict->nsbreconvliterals += nreconvliterals;
         if( upconflict != NULL )
            *upconflict = (nclauses > 0);

         /* reset the lower bound */
         col->lb = oldlb;
         CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );
            
         /* reset LP basis */
         CHECK_OKAY( SCIPlpiSetBase(lp->lpi, cstat, rstat) );
      }
   }

   /* free temporary memory for storing current LP basis */
   SCIPsetFreeBufferArray(set, &rstat);
   SCIPsetFreeBufferArray(set, &cstat);

   /* stop timing */
   SCIPclockStop(conflict->sbanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing infeasible strong branching conflicts */
Real SCIPconflictGetStrongbranchTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->sbanalyzetime);
}

/** gets number of calls to infeasible strong branching conflict analysis */
Longint SCIPconflictGetNStrongbranchCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbcalls;
}

/** gets number of conflict clauses detected in infeasible strong branching conflict analysis */
Longint SCIPconflictGetNStrongbranchConflictClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfclauses;
}

/** gets total number of literals in conflict clauses created in infeasible strong branching conflict analysis */
Longint SCIPconflictGetNStrongbranchConflictLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfliterals;
}

/** gets number of reconvergence clauses detected in infeasible strong branching conflict analysis */
Longint SCIPconflictGetNStrongbranchReconvergenceClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvclauses;
}

/** gets total number of literals in reconvergence clauses created in infeasible strong branching conflict analysis */
Longint SCIPconflictGetNStrongbranchReconvergenceLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvliterals;
}

/** gets number of LP iterations in infeasible strong branching conflict analysis */
Longint SCIPconflictGetNStrongbranchIterations(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbiterations;
}




/*
 * pseudo solution conflict analysis
 */

/** analyzes a pseudo solution with objective value exceeding the current cutoff to find out the bound changes on 
 *  variables that were responsible for the objective value degradation;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for pseudo solution conflict analysis
 */
RETCODE SCIPconflictAnalyzePseudo(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< LP data */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   VAR** vars;
   VAR* var;
   Real* curvarlbs;
   Real* curvarubs;
   int* lbchginfoposs;
   int* ubchginfoposs;
   Real* pseudocoefs;
   Real pseudolhs;
   Real pseudoact;
   int maxsize;
   int nvars;
   int v;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(lp != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if pseudo solution conflict analysis is enabled */
   if( !set->conf_usepseudo )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* calculate the maximal size of the conflict set */
   maxsize = (int)(set->conf_maxvarsfac * prob->nbinvars);
   maxsize = MAX(maxsize, set->conf_minmaxvars);
   if( maxsize < 2 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->pseudoanalyzetime, set);

   conflict->npseudocalls++;

   debugMessage("analyzing pseudo solution (obj: %g) that exceeds objective limit (%g)\n",
      SCIPlpGetPseudoObjval(lp, set), lp->cutoffbound);

   vars = prob->vars;
   nvars = prob->nvars;
   assert(nvars == 0 || vars != NULL);

   /* The current primal bound c* gives an upper bound for the current pseudo objective value:
    *   min{c^T x | lb <= x <= ub} <= c*.
    * We have to transform this row into a >= inequality in order to use methods above:
    *                          -c* <= max{-c^T x | lb <= x <= ub}.
    * In the local subproblem, this row is violated. We want to undo bound changes while still keeping the
    * row violated.
    */

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information 
    * positions in variable's bound change information arrays
    */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &lbchginfoposs, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &ubchginfoposs, nvars) );

   /* get temporary memory for infeasibility proof coefficients */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &pseudocoefs, nvars) );

   /* use a slightly tighter cutoff bound, because solutions with equal objective value should also be declared
    * infeasible
    */
   pseudolhs = -(lp->cutoffbound - SCIPsetSumepsilon(set));

   /* store the objective values as infeasibility proof coefficients, and recalculate the pseudo activity */
   pseudoact = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      pseudocoefs[v] = -SCIPvarGetObj(var);
      curvarlbs[v] = SCIPvarGetLbLocal(var);
      curvarubs[v] = SCIPvarGetUbLocal(var);
      if( pseudocoefs[v] > 0.0 )
         pseudoact += pseudocoefs[v] * curvarubs[v];
      else
         pseudoact += pseudocoefs[v] * curvarlbs[v];
      lbchginfoposs[v] = var->nlbchginfos-1;
      ubchginfoposs[v] = var->nubchginfos-1;
   }
   assert(SCIPsetIsFeasEQ(set, pseudoact, -SCIPlpGetPseudoObjval(lp, set)));
   debugMessage("  -> recalculated pseudo infeasibility proof:  %g <= %g\n", pseudolhs, pseudoact);

   /* check, if the pseudo row is still violated (after recalculation of pseudo activity) */
   if( SCIPsetIsFeasGT(set, pseudolhs, pseudoact) )
   {
      int nclauses;
      int nliterals;
      int nreconvclauses;
      int nreconvliterals;

      /* undo bound changes without destroying the infeasibility proof */
      CHECK_OKAY( undoBdchgsProof(set, prob, pseudocoefs, pseudolhs, pseudoact, 
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* analyze conflict on remaining bound changes */
      CHECK_OKAY( conflictAnalyzeRemainingBdchgs(conflict, memhdr, set, stat, prob, tree,
            maxsize, lbchginfoposs, ubchginfoposs, &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
      conflict->npseudoconfclauses += nclauses;
      conflict->npseudoconfliterals += nliterals;
      conflict->npseudoreconvclauses += nreconvclauses;
      conflict->npseudoreconvliterals += nreconvliterals;
      if( success != NULL )
         *success = (nclauses > 0);
   }

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &pseudocoefs);
   SCIPsetFreeBufferArray(set, &curvarubs);
   SCIPsetFreeBufferArray(set, &curvarlbs);
   SCIPsetFreeBufferArray(set, &ubchginfoposs);
   SCIPsetFreeBufferArray(set, &lbchginfoposs);

   /* stop timing */
   SCIPclockStop(conflict->pseudoanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing pseudo solution conflicts */
Real SCIPconflictGetPseudoTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->pseudoanalyzetime);
}

/** gets number of calls to pseudo solution conflict analysis */
Longint SCIPconflictGetNPseudoCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudocalls;
}

/** gets number of conflict clauses detected in pseudo solution conflict analysis */
Longint SCIPconflictGetNPseudoConflictClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfclauses;
}

/** gets total number of literals in conflict clauses created in pseudo solution conflict analysis */
Longint SCIPconflictGetNPseudoConflictLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfliterals;
}

/** gets number of reconvergence clauses detected in pseudo solution conflict analysis */
Longint SCIPconflictGetNPseudoReconvergenceClauses(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvclauses;
}

/** gets total number of literals in reconvergence clauses created in pseudo solution conflict analysis */
Longint SCIPconflictGetNPseudoReconvergenceLiterals(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvliterals;
}

