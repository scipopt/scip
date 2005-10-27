/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: conflict.c,v 1.107 2005/10/27 16:36:47 bzfpfend Exp $"

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

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/clock.h"
#include "scip/vbc.h"
#include "scip/lpi.h"
#include "scip/misc.h"
#include "scip/paramset.h"
#include "scip/lp.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/tree.h"
#include "scip/scip.h"
#include "scip/conflict.h"
#include "scip/cons.h"
#include "scip/prop.h"
#include "scip/debug.h"

#include "scip/struct_conflict.h"



#define CLAUSESCORE(clause) (-(clause)->nvars - 100*(clause)->repropdepth - 1000*(clause)->validdepth)



/*
 * Conflict Handler
 */

/** compares two conflict handlers w. r. to their priority */
SCIP_DECL_SORTPTRCOMP(SCIPconflicthdlrComp)
{  /*lint --e{715}*/
   return ((SCIP_CONFLICTHDLR*)elem2)->priority - ((SCIP_CONFLICTHDLR*)elem1)->priority;
}

/** method to call, when the priority of a conflict handler was changed */
static
SCIP_DECL_PARAMCHGD(paramChgdConflicthdlrPriority)
{  /*lint --e{715}*/
   SCIP_PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetConflicthdlrPriority() to mark the conflicthdlrs unsorted */
   SCIP_CALL( SCIPsetConflicthdlrPriority(scip, (SCIP_CONFLICTHDLR*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a conflict handler */
SCIP_RETCODE SCIPconflicthdlrCreate(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of conflict handler */
   const char*           desc,               /**< description of conflict handler */
   int                   priority,           /**< priority of the conflict handler */
   SCIP_DECL_CONFLICTFREE((*conflictfree)),  /**< destructor of conflict handler */
   SCIP_DECL_CONFLICTINIT((*conflictinit)),  /**< initialize conflict handler */
   SCIP_DECL_CONFLICTEXIT((*conflictexit)),  /**< deinitialize conflict handler */
   SCIP_DECL_CONFLICTINITSOL((*conflictinitsol)),/**< solving process initialization method of conflict handler */
   SCIP_DECL_CONFLICTEXITSOL((*conflictexitsol)),/**< solving process deinitialization method of conflict handler */
   SCIP_DECL_CONFLICTEXEC((*conflictexec)),  /**< conflict processing method of conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< conflict handler data */
   )
{
   char paramname[SCIP_MAXSTRLEN];
   char paramdesc[SCIP_MAXSTRLEN];

   assert(conflicthdlr != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   SCIP_ALLOC( BMSallocMemory(conflicthdlr) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*conflicthdlr)->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&(*conflicthdlr)->desc, desc, strlen(desc)+1) );
   (*conflicthdlr)->priority = priority;
   (*conflicthdlr)->conflictfree = conflictfree;
   (*conflicthdlr)->conflictinit = conflictinit;
   (*conflicthdlr)->conflictexit = conflictexit;
   (*conflicthdlr)->conflictinitsol = conflictinitsol;
   (*conflicthdlr)->conflictexitsol = conflictexitsol;
   (*conflicthdlr)->conflictexec = conflictexec;
   (*conflicthdlr)->conflicthdlrdata = conflicthdlrdata;
   (*conflicthdlr)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "conflict/%s/priority", name);
   sprintf(paramdesc, "priority of conflict handler <%s>", name);
   SCIP_CALL( SCIPsetAddIntParam(set, blkmem, paramname, paramdesc,
         &(*conflicthdlr)->priority, priority, INT_MIN, INT_MAX,
         paramChgdConflicthdlrPriority, (SCIP_PARAMDATA*)(*conflicthdlr)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of conflict handler */
SCIP_RETCODE SCIPconflicthdlrFree(
   SCIP_CONFLICTHDLR**   conflicthdlr,       /**< pointer to conflict handler data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(*conflicthdlr != NULL);
   assert(!(*conflicthdlr)->initialized);
   assert(set != NULL);

   /* call destructor of conflict handler */
   if( (*conflicthdlr)->conflictfree != NULL )
   {
      SCIP_CALL( (*conflicthdlr)->conflictfree(set->scip, *conflicthdlr) );
   }

   BMSfreeMemoryArray(&(*conflicthdlr)->name);
   BMSfreeMemoryArray(&(*conflicthdlr)->desc);
   BMSfreeMemory(conflicthdlr);

   return SCIP_OKAY;
}

/** calls init method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrInit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   if( conflicthdlr->initialized )
   {
      SCIPerrorMessage("conflict handler <%s> already initialized\n", conflicthdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call initialization method of conflict handler */
   if( conflicthdlr->conflictinit != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictinit(set->scip, conflicthdlr) );
   }
   conflicthdlr->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExit(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   if( !conflicthdlr->initialized )
   {
      SCIPerrorMessage("conflict handler <%s> not initialized\n", conflicthdlr->name);
      return SCIP_INVALIDCALL;
   }

   /* call deinitialization method of conflict handler */
   if( conflicthdlr->conflictexit != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictexit(set->scip, conflicthdlr) );
   }
   conflicthdlr->initialized = FALSE;

   return SCIP_OKAY;
}

/** informs conflict handler that the branch and bound process is being started */
SCIP_RETCODE SCIPconflicthdlrInitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   /* call solving process initialization method of conflict handler */
   if( conflicthdlr->conflictinitsol != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictinitsol(set->scip, conflicthdlr) );
   }

   return SCIP_OKAY;
}

/** informs conflict handler that the branch and bound process data is being freed */
SCIP_RETCODE SCIPconflicthdlrExitsol(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   /* call solving process deinitialization method of conflict handler */
   if( conflicthdlr->conflictexitsol != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictexitsol(set->scip, conflicthdlr) );
   }

   return SCIP_OKAY;
}

/** calls execution method of conflict handler */
SCIP_RETCODE SCIPconflicthdlrExec(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_NODE*            node,               /**< node to add conflict constraint to */
   SCIP_NODE*            validnode,          /**< node at which the constraint is valid */
   SCIP_VAR**            conflictvars,       /**< variables of the conflict set */
   int                   nconflictvars,      /**< number of variables in the conflict set */
   SCIP_Bool             resolved,           /**< is the conflict set already used to create a constraint? */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);
   assert(conflictvars != NULL || nconflictvars == 0);
   assert(result != NULL);

   /* call solution start method of conflict handler */
   *result = SCIP_DIDNOTRUN;
   if( conflicthdlr->conflictexec != NULL )
   {
      SCIP_CALL( conflicthdlr->conflictexec(set->scip, conflicthdlr, node, validnode, conflictvars, nconflictvars,
            (SCIPnodeGetDepth(validnode) > 0), set->conf_dynamic, set->conf_removeable, resolved, result) );

      if( *result != SCIP_CONSADDED
         && *result != SCIP_DIDNOTFIND
         && *result != SCIP_DIDNOTRUN )
      {
         SCIPerrorMessage("execution method of conflict handler <%s> returned invalid result <%d>\n",
            conflicthdlr->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** gets user data of conflict handler */
SCIP_CONFLICTHDLRDATA* SCIPconflicthdlrGetData(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->conflicthdlrdata;
}

/** sets user data of conflict handler; user has to free old data in advance! */
void SCIPconflicthdlrSetData(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_CONFLICTHDLRDATA* conflicthdlrdata   /**< new conflict handler user data */
   )
{
   assert(conflicthdlr != NULL);

   conflicthdlr->conflicthdlrdata = conflicthdlrdata;
}

/** gets name of conflict handler */
const char* SCIPconflicthdlrGetName(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->name;
}

/** gets description of conflict handler */
const char* SCIPconflicthdlrGetDesc(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->desc;
}

/** gets priority of conflict handler */
int SCIPconflicthdlrGetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->priority;
}

/** sets priority of conflict handler */
void SCIPconflicthdlrSetPriority(
   SCIP_CONFLICTHDLR*    conflicthdlr,       /**< conflict handler */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the conflict handler */
   )
{
   assert(conflicthdlr != NULL);
   assert(set != NULL);

   conflicthdlr->priority = priority;
   set->conflicthdlrssorted = FALSE;
}

/** is conflict handler initialized? */
SCIP_Bool SCIPconflicthdlrIsInitialized(
   SCIP_CONFLICTHDLR*    conflicthdlr        /**< conflict handler */
   )
{
   assert(conflicthdlr != NULL);

   return conflicthdlr->initialized;
}




/*
 * Conflict Clause Storage
 */

/** creates a clause */
static
SCIP_RETCODE clauseCreate(
   SCIP_CLAUSE**         clause,             /**< pointer to store the clause */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_VAR**            vars,               /**< literals of the clause */
   int                   nvars,              /**< number of literals in the clause */
   int                   validdepth,         /**< depth in the tree where the clause is valid */
   int                   insertdepth         /**< minimal depth where the clause should be added */
   )
{
   SCIP_VAR** cvars;
   int i;
   int maxdepth[2];

   assert(clause != NULL);
   assert(nvars == 0 || vars != NULL);
   assert(validdepth <= insertdepth);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, clause) );
   SCIP_ALLOC( BMSallocBlockMemoryArray(blkmem, &(*clause)->vars, nvars) );
   (*clause)->nvars = nvars;
   (*clause)->validdepth = validdepth;
   (*clause)->insertdepth = insertdepth;

   /* sort the variables, and identify the depth of the last and last but one fixing */
   maxdepth[0] = validdepth;
   maxdepth[1] = validdepth;
   cvars = (*clause)->vars;
   for( i = 0; i < nvars; ++i )
   {
      int idx;
      int j;
      int depth;

      idx = SCIPvarGetIndex(vars[i]);
      for( j = i; j > 0 && idx < SCIPvarGetIndex(cvars[j-1]); --j )
         cvars[j] = cvars[j-1];
      assert(j == 0 || idx > SCIPvarGetIndex(cvars[j-1]));
      cvars[j] = vars[i];

      depth = SCIPvarGetLastBdchgDepth(vars[i]);
      assert(depth == -2 || depth >= 0);

      if( depth == -2 )
         depth = INT_MAX; /* variable is currently not fixed -> strong branching, diving, or probing bound change */
      if( depth > maxdepth[0] )
      {
         maxdepth[1] = maxdepth[0];
         maxdepth[0] = depth;
      }
      else if( depth > maxdepth[1] )
         maxdepth[1] = depth;
   }
   assert(maxdepth[0] >= maxdepth[1]);

   (*clause)->conflictdepth = maxdepth[0];
   (*clause)->repropdepth = maxdepth[1];

   (*clause)->score = CLAUSESCORE(*clause);

   return SCIP_OKAY;
}

/** frees a clause */
static
SCIP_RETCODE clauseFree(
   SCIP_CLAUSE**         clause,             /**< pointer to the clause */
   BMS_BLKMEM*           blkmem              /**< block memory of transformed problem */
   )
{
   assert(clause != NULL);
   assert(*clause != NULL);

   BMSfreeBlockMemoryArray(blkmem, &(*clause)->vars, (*clause)->nvars);
   BMSfreeBlockMemory(blkmem, clause);

   return SCIP_OKAY;
}

/** checks whether the first clause is redundant to the second one */
static
SCIP_Bool clauseIsRedundant(
   SCIP_CLAUSE*          clause1,            /**< first conflict clause */
   SCIP_CLAUSE*          clause2             /**< second conflict clause */
   )
{
   int i1;
   int i2;

   assert(clause1 != NULL);
   assert(clause2 != NULL);

   /* if clause1 has smaller validdepth, it is definitely not redundant to clause2 */
   if( clause1->validdepth < clause2->validdepth )
      return FALSE;

   /* check, if all literals in clause2 are also present in clause1;
    * we can stop immediately, if more literals are remaining in clause2 than in clause1
    */
   for( i1 = 0, i2 = 0; i2 < clause2->nvars && clause1->nvars - i1 < clause2->nvars - i2; ++i1, ++i2 )
   {
      int idx;

      assert(i2 == 0 || SCIPvarGetIndex(clause2->vars[i2-1]) < SCIPvarGetIndex(clause2->vars[i2]));

      idx = SCIPvarGetIndex(clause2->vars[i2]);
      for( ; i1 < clause1->nvars && SCIPvarGetIndex(clause1->vars[i1]) < idx; ++i1 )
      {
         /* while scanning clause1, check consistency */
         assert(i1 == 0 || SCIPvarGetIndex(clause1->vars[i1-1]) < SCIPvarGetIndex(clause1->vars[i1]));
      }
      if( i1 >= clause1->nvars || SCIPvarGetIndex(clause1->vars[i1]) > idx )
         return FALSE;
   }

   return (i2 == clause2->nvars);
}

#ifdef SCIP_DEBUG
/** prints a conflict clause to the screen */
static
void clausePrint(
   SCIP_CLAUSE*          clause              /**< conflict clause */
   )
{
   int i;

   assert(clause != NULL);
   for( i = 0; i < clause->nvars; ++i )
      SCIPdebugPrintf(" <%s>[%d]", SCIPvarGetName(clause->vars[i]), SCIPvarGetLastBdchgDepth(clause->vars[i]));
   SCIPdebugPrintf("\n");
}
#endif

/** resizes clauses array to be able to store at least num entries */
static
SCIP_RETCODE conflictEnsureClausesMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflict != NULL);
   assert(set != NULL);

   if( num > conflict->clausessize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->clauses, newsize) );
      conflict->clausessize = newsize;
   }
   assert(num <= conflict->clausessize);

   return SCIP_OKAY;
}

/** inserts conflict clause into sorted clauses array */
static
SCIP_RETCODE conflictInsertClause(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR**            vars,               /**< literals of the clause */
   int                   nvars,              /**< number of literals in the clause */
   int                   validdepth,         /**< depth in the tree where the clause is valid */
   int                   insertdepth         /**< minimal depth where the clause should be added */
   )
{
   SCIP_CLAUSE* clause;
   int pos;
   int i;
   int j;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(validdepth <= insertdepth);
   assert(set->conf_allowlocal || validdepth == 0);

   /* create a out of the given conflict set */
   SCIP_CALL( clauseCreate(&clause, blkmem, vars, nvars, validdepth, insertdepth) );
   assert(clause != NULL);
   assert(0 <= clause->validdepth && clause->validdepth <= clause->insertdepth);

   /* if we apply repropagations, the clause should be inserted at most at its repropdepth */
   if( set->conf_repropagate )
      clause->insertdepth = MIN(clause->insertdepth, clause->repropdepth);
   else
      clause->repropdepth = INT_MAX;
   assert(clause->insertdepth <= clause->repropdepth);

   SCIPdebugMessage("inserting conflict clause (valid: %d, insert: %d, conf: %d, reprop: %d):",
      clause->validdepth, clause->insertdepth, clause->conflictdepth, clause->repropdepth);
   SCIPdebug(clausePrint(clause));

   /* check, if clause is redundant to a better clause */
   for( pos = 0; pos < conflict->nclauses && clause->score < conflict->clauses[pos]->score; ++pos )
   {
      /* check if clause is redundant with respect to clauses[pos] */
      if( clauseIsRedundant(clause, conflict->clauses[pos]) )
      {
         SCIP_CALL( clauseFree(&clause, blkmem) );
         return SCIP_OKAY;
      }
      /**@todo like in sepastore.c: calculate overlap between clauses -> large overlap reduces score */
      /*????????????????????????*/
   }

   /* insert clause into the sorted clauses array*/
   SCIP_CALL( conflictEnsureClausesMem(conflict, set, conflict->nclauses + 1) );
   for( i = conflict->nclauses; i > pos; --i )
   {
      assert(clause->score >= conflict->clauses[i-1]->score);
      conflict->clauses[i] = conflict->clauses[i-1];
   }
   conflict->clauses[pos] = clause;
   conflict->nclauses++;

   /* remove worse clauses that are redundant to the new clause */
   for( i = pos+1, j = pos+1; i < conflict->nclauses; ++i )
   {
      if( clauseIsRedundant(conflict->clauses[i], clause) )
      {
         SCIP_CALL( clauseFree(&conflict->clauses[i], blkmem) );
      }
      else
      {
         assert(j <= i);
         conflict->clauses[j] = conflict->clauses[i];
         j++;
      }
   }
   assert(j <= conflict->nclauses);
   conflict->nclauses = j;

   return SCIP_OKAY;
}

/** calculates the maximal size of conflict clauses to be used */
static
int conflictCalcMaxsize(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob                /**< problem data */
   )
{
   int maxsize;

   assert(set != NULL);
   assert(prob != NULL);

   maxsize = (int)(set->conf_maxvarsfac * prob->nbinvars);
   maxsize = MAX(maxsize, set->conf_minmaxvars);

   return maxsize;
}

#if 0 /*??????????????????????????*/
/** adds the collected conflict clauses to the corresponding nodes; the best set->conf_maxclauses clauses are added
 *  to the node of their validdepth; the remaining clauses are added to the node of their repropdepth in order to just
 *  trigger the deduction but not get introduced to a more global subtree
 */
SCIP_RETCODE SCIPconflictFlushClauses(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   int i;
   int nclausesused;
   int currentdepth;
   int cutoffdepth;
   int repropdepth;
   int maxclauses;
   int maxsize;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);

   /* is there anything to do? */
   if( conflict->nclauses == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of conflict clauses to accept, and the maximal size of each accepted conflict clause */
   maxclauses = (set->conf_maxclauses == -1 ? INT_MAX : set->conf_maxclauses);
   maxsize = conflictCalcMaxsize(set, prob);

   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);

   SCIPdebugMessage("flushing %d conflict clauses at depth %d (maxclauses: %d, maxsize: %d)\n",
      conflict->nclauses, currentdepth, maxclauses, maxsize);

   /* mark the current node to have produced conflict clauses in the VBC tool output */
   SCIPvbcFoundConflict(stat->vbc, stat, tree->path[currentdepth]);

   /* insert the conflict clauses at the corresponding nodes */
   nclausesused = 0;
   cutoffdepth = currentdepth+1;
   repropdepth = currentdepth+1;
   for( i = 0; i < conflict->nclauses && nclausesused < maxclauses; ++i )
   {
      SCIP_CLAUSE* clause;
      int insertdepth;

      clause = conflict->clauses[i];
      assert(clause != NULL);
      assert(0 <= clause->validdepth && clause->validdepth <= clause->insertdepth && clause->insertdepth <= currentdepth);
      assert(clause->insertdepth <= clause->repropdepth);
      assert(clause->repropdepth <= currentdepth || clause->repropdepth == INT_MAX);
      assert(clause->conflictdepth <= currentdepth || clause->conflictdepth == INT_MAX); /* INT_MAX for dive/strong */

      /* ignore clauses that are only valid at a node that was already cut off */
      if( clause->validdepth >= cutoffdepth )
      {
         SCIPdebugMessage("ignoring clause with validdepth %d >= cutoffdepth %d\n", clause->validdepth, cutoffdepth);
         continue;
      }

      /* if no conflict variables exist, the node and its sub tree in the conflict clause's valid depth can be
       * cut off completely
       */
      if( clause->nvars == 0 )
      {
         SCIPdebugMessage("empty conflict clause in depth %d cuts off sub tree at depth %d\n",
            currentdepth, clause->validdepth);

         SCIPnodeCutoff(tree->path[clause->validdepth], set, stat, tree);
         cutoffdepth = clause->validdepth;
         continue;
      }

      /* if the clause is too long, use the clause only as triggering clause
       * (i.e. insert it at the node of its repropdepth)
       */
      if( clause->nvars > maxsize )
         insertdepth = clause->repropdepth;
      else
         insertdepth = clause->insertdepth;

      if( insertdepth < cutoffdepth && (insertdepth < currentdepth || clause->conflictdepth == INT_MAX) )
      {
         SCIP_Bool success;
         int h;

         assert(insertdepth <= currentdepth);

         /* sort conflict handlers by priority */
         SCIPsetSortConflicthdlrs(set);

         /* call conflict handlers to create a conflict constraint */
         success = FALSE;
         for( h = 0; h < set->nconflicthdlrs; ++h )
         {
            SCIP_RESULT result;

            SCIP_CALL( SCIPconflicthdlrExec(set->conflicthdlrs[h], set, tree->path[insertdepth],
                  tree->path[clause->validdepth], clause->vars, clause->nvars, success, &result) );
            if( result == SCIP_CONSADDED )
            {
               success = TRUE;
               if( insertdepth > 0 )
               {
                  conflict->nappliedlocclauses++;
                  conflict->nappliedlocliterals += clause->nvars;
               }
               else
               {
                  conflict->nappliedglbclauses++;
                  conflict->nappliedglbliterals += clause->nvars;
               }
            }
            SCIPdebugMessage(" -> calling conflict handler <%s> (prio=%d) to create conflict clause with %d literals returned result %d\n",
               SCIPconflicthdlrGetName(set->conflicthdlrs[h]), SCIPconflicthdlrGetPriority(set->conflicthdlrs[h]),
               clause->nvars, result);
         }

         if( success )
         {
#ifdef SCIP_DEBUG
            SCIPdebugMessage(" -> conflict clause added at depth %d (valid:%d, conf:%d, reprop:%d, len:%d):",
               insertdepth, clause->validdepth, clause->conflictdepth, clause->repropdepth, clause->nvars);
            SCIPdebug(clausePrint(clause));
#endif

            repropdepth = MIN(repropdepth, clause->repropdepth);
            nclausesused++;
         }
      }
   }

   /* reactivate propagation on the first node where one of the new conflict clauses trigger a deduction */
   if( set->conf_repropagate && repropdepth < cutoffdepth && repropdepth < currentdepth )
   {
      assert(0 <= repropdepth && repropdepth < tree->pathlen);
      assert(tree->path[repropdepth]->depth == repropdepth);
      SCIPnodePropagateAgain(tree->path[repropdepth], set, stat, tree);

      SCIPdebugMessage("marked node %p in depth %d to be repropagated due to conflicts found in depth %d\n",
         tree->path[repropdepth], repropdepth, currentdepth);
   }

   /* free the conflict storage */
   for( i = 0; i < conflict->nclauses; ++i )
   {
      SCIP_CALL( clauseFree(&conflict->clauses[i], blkmem) );
   }
   conflict->nclauses = 0;

   return SCIP_OKAY;
}
#else
/** adds the given clause as conflict constraint to the problem */
static
SCIP_RETCODE conflictAddClauseCons(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_CLAUSE*          clause,             /**< clause to add to the tree */
   int                   insertdepth,        /**< depth level at which the clause should be added */
   SCIP_Bool*            success             /**< pointer to store whether the addition was successful */
   )
{
   int h;

   assert(conflict != NULL);
   assert(tree != NULL);
   assert(tree->path != NULL);
   assert(clause != NULL);
   assert(clause->validdepth <= insertdepth);
   assert(success != NULL);

   /* sort conflict handlers by priority */
   SCIPsetSortConflicthdlrs(set);

   /* call conflict handlers to create a conflict constraint */
   *success = FALSE;
   for( h = 0; h < set->nconflicthdlrs; ++h )
   {
      SCIP_RESULT result;

      SCIP_CALL( SCIPconflicthdlrExec(set->conflicthdlrs[h], set, tree->path[insertdepth],
            tree->path[clause->validdepth], clause->vars, clause->nvars, *success, &result) );
      if( result == SCIP_CONSADDED )
      {
         *success = TRUE;
         if( insertdepth > 0 )
         {
            conflict->nappliedlocclauses++;
            conflict->nappliedlocliterals += clause->nvars;
         }
         else
         {
            conflict->nappliedglbclauses++;
            conflict->nappliedglbliterals += clause->nvars;
         }
      }
      SCIPdebugMessage(" -> calling conflict handler <%s> (prio=%d) to create conflict clause with %d literals returned result %d\n",
         SCIPconflicthdlrGetName(set->conflicthdlrs[h]), SCIPconflicthdlrGetPriority(set->conflicthdlrs[h]),
         clause->nvars, result);
   }

   return SCIP_OKAY;
}

/** adds the collected conflict clauses to the corresponding nodes; the best set->conf_maxclauses clauses are added
 *  to the node of their validdepth; additionally (if not yet added, and if repropagation is activated), the clause that
 *  triggers the earliest repropagation is added to the node of its validdepth
 */
SCIP_RETCODE SCIPconflictFlushClauses(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree                /**< branch and bound tree */
   )
{
   SCIP_CLAUSE* repropclause;
   int nclausesused;
   int currentdepth;
   int cutoffdepth;
   int repropdepth;
   int maxclauses;
   int maxsize;
   int i;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(tree != NULL);

   /* is there anything to do? */
   if( conflict->nclauses == 0 )
      return SCIP_OKAY;

   /* calculate the maximal number of conflict clauses to accept, and the maximal size of each accepted conflict clause */
   maxclauses = (set->conf_maxclauses == -1 ? INT_MAX : set->conf_maxclauses);
   maxsize = conflictCalcMaxsize(set, prob);

   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);

   SCIPdebugMessage("flushing %d conflict clauses at depth %d (maxclauses: %d, maxsize: %d)\n",
      conflict->nclauses, currentdepth, maxclauses, maxsize);

   /* mark the current node to have produced conflict clauses in the VBC tool output */
   SCIPvbcFoundConflict(stat->vbc, stat, tree->path[currentdepth]);

   /* insert the conflict clauses at the corresponding nodes */
   nclausesused = 0;
   cutoffdepth = currentdepth+1;
   repropdepth = currentdepth+1;
   repropclause = NULL;
   for( i = 0; i < conflict->nclauses && nclausesused < maxclauses; ++i )
   {
      SCIP_CLAUSE* clause;

      clause = conflict->clauses[i];
      assert(clause != NULL);
      assert(0 <= clause->validdepth && clause->validdepth <= clause->insertdepth && clause->insertdepth <= currentdepth);
      assert(clause->insertdepth <= currentdepth);
      assert(clause->insertdepth <= clause->repropdepth);
      assert(clause->repropdepth <= currentdepth || clause->repropdepth == INT_MAX);
      assert(clause->conflictdepth <= currentdepth || clause->conflictdepth == INT_MAX); /* INT_MAX for dive/strong */

      /* ignore clauses that are only valid at a node that was already cut off */
      if( clause->insertdepth >= cutoffdepth )
      {
         SCIPdebugMessage(" -> ignoring clause with insertdepth %d >= cutoffdepth %d\n", clause->validdepth, cutoffdepth);
         continue;
      }

      /* ignore clauses that are only valid in the current node (which is cut off); if the conflict analysis was applied
       * in strong branching, probing, or diving, the current node is not cut off -> the conflict clause is useful
       */
      if( clause->insertdepth >= currentdepth && clause->conflictdepth != INT_MAX )
      {
         SCIPdebugMessage(" -> ignoring clause with insertdepth %d >= currentdepth %d\n", clause->validdepth, currentdepth);
         continue;
      }

      /* if no conflict variables exist, the node and its sub tree in the conflict clause's valid depth can be
       * cut off completely
       */
      if( clause->nvars == 0 )
      {
         SCIPdebugMessage(" -> empty conflict clause in depth %d cuts off sub tree at depth %d\n",
            currentdepth, clause->validdepth);

         SCIPnodeCutoff(tree->path[clause->validdepth], set, stat, tree);
         cutoffdepth = clause->validdepth;
         continue;
      }

      /* if the clause is too long, use the clause only if it decreases the repropagation depth */
      if( clause->nvars > maxsize )
      {
         SCIPdebugMessage(" -> clause is too long: %d > %d literals\n", clause->nvars, maxsize);
         if( clause->repropdepth < repropdepth )
         {
            repropdepth = clause->repropdepth;
            repropclause = clause;
         }
      }
      else
      {
         SCIP_Bool success;

         /* call conflict handlers to create a conflict constraint */
         SCIP_CALL( conflictAddClauseCons(conflict, set, stat, prob, tree, clause, clause->insertdepth, &success) );

         if( success )
         {
            SCIPdebugMessage(" -> conflict clause %d/%d added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf:%d, reprop:%d, len:%d):\n",
               nclausesused+1, maxclauses, SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
               clause->insertdepth, clause->validdepth, clause->conflictdepth, clause->repropdepth, clause->nvars);
            SCIPdebug(clausePrint(clause));

            if( clause->repropdepth <= repropdepth )
            {
               repropdepth = clause->repropdepth;
               repropclause = NULL;
            }
            nclausesused++;
         }
      }
   }

   /* reactivate propagation on the first node where one of the new conflict clauses trigger a deduction */
   if( set->conf_repropagate && repropdepth < cutoffdepth && repropdepth < currentdepth )
   {
      assert(0 <= repropdepth && repropdepth < tree->pathlen);
      assert(tree->path[repropdepth]->depth == repropdepth);

      /* if the conflict constraint of smallest repropagation depth was not yet added, insert it now */
      if( repropclause != NULL )
      {
         SCIP_Bool success;

         SCIP_CALL( conflictAddClauseCons(conflict, set, stat, prob, tree, repropclause, repropdepth, &success) );
#ifdef SCIP_DEBUG
         if( success )
         {
            SCIPdebugMessage(" -> additional reprop conflict clause added (cdpt:%d, fdpt:%d, insert:%d, valid:%d, conf:%d, reprop:%d, len:%d):\n",
               SCIPtreeGetCurrentDepth(tree), SCIPtreeGetFocusDepth(tree),
               repropclause->insertdepth, repropclause->validdepth, repropclause->conflictdepth,
               repropclause->repropdepth, repropclause->nvars);
            SCIPdebug(clausePrint(repropclause));
         }
#endif
      }

      /* mark the node in the repropdepth to be propagated again */
      SCIPnodePropagateAgain(tree->path[repropdepth], set, stat, tree);

      SCIPdebugMessage("marked node %p in depth %d to be repropagated due to conflicts found in depth %d\n",
         tree->path[repropdepth], repropdepth, currentdepth);
   }

   /* free the conflict storage */
   for( i = 0; i < conflict->nclauses; ++i )
   {
      SCIP_CALL( clauseFree(&conflict->clauses[i], blkmem) );
   }
   conflict->nclauses = 0;

   return SCIP_OKAY;
}
#endif

/** returns the current number of conflict clauses in the conflict clause storage */
int SCIPconflictGetNClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nclauses;
}

/** returns the total number of conflict clauses that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbclauses + conflict->nappliedlocclauses;
}

/** returns the total number of literals in conflict clauses that were added to the problem */
SCIP_Longint SCIPconflictGetNAppliedLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbliterals + conflict->nappliedlocliterals;
}

/** returns the total number of conflict clauses that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbclauses;
}

/** returns the total number of literals in conflict clauses that were added globally to the problem */
SCIP_Longint SCIPconflictGetNAppliedGlobalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedglbliterals;
}

/** returns the total number of conflict clauses that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedlocclauses;
}

/** returns the total number of literals in conflict clauses that were added locally to the problem */
SCIP_Longint SCIPconflictGetNAppliedLocalLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nappliedlocliterals;
}




/*
 * Propagation Conflict Analysis
 */

/** resizes conflictvars array to be able to store at least num constraints */
static
SCIP_RETCODE conflictEnsureConflictvarsMem(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimal number of slots in array */
   )
{
   assert(conflict != NULL);

   if( num > conflict->conflictvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&conflict->conflictvars, newsize) );
      conflict->conflictvarssize = newsize;
   }
   assert(num <= conflict->conflictvarssize);

   return SCIP_OKAY;
}

/** puts given locally fixed binary variable or its negation into conflict set, depending which one is fixed to zero */
static
SCIP_RETCODE conflictAddConflictVar(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_VAR*             var,                /**< variable that should enter the conflict set */
   SCIP_Bool             literalvalue,       /**< current value of the literal, that should enter the conflict set */
   SCIP_Bool             temporary           /**< should the variable be added only temporary to the conflict set? */
   )
{
   assert(conflict != NULL);
   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(SCIPsetIsEQ(set, SCIPvarGetLbLP(var), SCIPvarGetUbLP(var)));
   assert(SCIPvarGetLbLP(var) > 0.5 || SCIPvarGetUbLP(var) < 0.5);

   /* choose variable or its negation, such that the literal is currently fixed to the desired value */
   if( (SCIPvarGetLbLP(var) > 0.5) != literalvalue )
   {
      /* variable is fixed to the opposite value -> use the negated variable as conflict literal */
      SCIP_CALL( SCIPvarNegate(var, blkmem, set, stat, &var) );
   }
   assert(SCIPsetIsEQ(set, SCIPvarGetLbLP(var), SCIPvarGetUbLP(var)));
   assert((SCIPvarGetLbLP(var) > 0.5) == literalvalue);

   SCIPdebugMessage("putting variable <%s> fixed to %d at depth %d to conflict set (temporary: %d)\n",
      SCIPvarGetName(var), (SCIPvarGetLbLP(var) > 0.5), SCIPvarGetLastBdchgDepth(var), temporary);

   if( var->conflictsetcount == conflict->count )
   {
      /* variable is already member of the conflict set */
#ifndef NDEBUG
      {
         int i;
         for( i = 0; i < conflict->nconflictvars + conflict->ntmpconflictvars && conflict->conflictvars[i] != var; ++i )
         {}
         assert(i < conflict->nconflictvars + conflict->ntmpconflictvars);
         assert(conflict->conflictvars[i] == var);
      }
#endif
      return SCIP_OKAY;
   }

   /* put variable to conflict set, sorted by non-increasing fixing depths */
   SCIP_CALL( conflictEnsureConflictvarsMem(conflict, set, conflict->nconflictvars + conflict->ntmpconflictvars + 1) );
   if( temporary )
   {
      /* put the variable into the temporary part of the conflict vars array */
      conflict->conflictvars[conflict->nconflictvars + conflict->ntmpconflictvars] = var;
      conflict->ntmpconflictvars++;
   }
   else
   {
      assert(conflict->ntmpconflictvars == 0);

      /* put the variable into the conflict vars array */
      conflict->conflictvars[conflict->nconflictvars] = var;
      conflict->nconflictvars++;
   }

   /* mark the active problem variable to belong to the current conflict */
   var->conflictsetcount = conflict->count;

   return SCIP_OKAY;
}

/** removes all temporary variables from the conflict set */
static
void conflictRemoveTemporaryConflictVars(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   int i;

   assert(conflict != NULL);
   assert(conflict->nconflictvars >= 0);

   for( i = conflict->nconflictvars; i < conflict->nconflictvars + conflict->ntmpconflictvars; ++i )
      conflict->conflictvars[i]->conflictsetcount = 0;
   conflict->ntmpconflictvars = 0;
}

/** compares two conflict set entries, such that bound changes infered later are
 *  ordered prior to ones that were infered earlier
 */
static
SCIP_DECL_SORTPTRCOMP(conflictBdchginfoComp)
{  /*lint --e{715}*/
   SCIP_BDCHGINFO* bdchginfo1;
   SCIP_BDCHGINFO* bdchginfo2;

   bdchginfo1 = (SCIP_BDCHGINFO*)elem1;
   bdchginfo2 = (SCIP_BDCHGINFO*)elem2;
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
SCIP_RETCODE SCIPconflictCreate(
   SCIP_CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   SCIP_SET*             set                 /**< global SCIP settings */
   )
{
   assert(conflict != NULL);

   SCIP_ALLOC( BMSallocMemory(conflict) );

   SCIP_CALL( SCIPclockCreate(&(*conflict)->propanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->lpanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->sbanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPclockCreate(&(*conflict)->pseudoanalyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   SCIP_CALL( SCIPpqueueCreate(&(*conflict)->binbdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp) );
   SCIP_CALL( SCIPpqueueCreate(&(*conflict)->nonbinbdchgqueue, set->mem_arraygrowinit, set->mem_arraygrowfac,
         conflictBdchginfoComp) );
   (*conflict)->clauses = NULL;
   (*conflict)->conflictvars = NULL;
   (*conflict)->conflictvarssize = 0;
   (*conflict)->nconflictvars = 0;
   (*conflict)->ntmpconflictvars = 0;
   (*conflict)->count = 0;
   (*conflict)->clausessize = 0;
   (*conflict)->nclauses = 0;
   (*conflict)->nappliedglbclauses = 0;
   (*conflict)->nappliedglbliterals = 0;
   (*conflict)->nappliedlocclauses = 0;
   (*conflict)->nappliedlocliterals = 0;
   (*conflict)->npropcalls = 0;
   (*conflict)->npropsuccess = 0;
   (*conflict)->npropconfclauses = 0;
   (*conflict)->npropconfliterals = 0;
   (*conflict)->npropreconvclauses = 0;
   (*conflict)->npropreconvliterals = 0;
   (*conflict)->nlpcalls = 0;
   (*conflict)->nlpsuccess = 0;
   (*conflict)->nlpconfclauses = 0;
   (*conflict)->nlpconfliterals = 0;
   (*conflict)->nlpreconvclauses = 0;
   (*conflict)->nlpreconvliterals = 0;
   (*conflict)->nlpiterations = 0;
   (*conflict)->nsbcalls = 0;
   (*conflict)->nsbsuccess = 0;
   (*conflict)->nsbconfclauses = 0;
   (*conflict)->nsbconfliterals = 0;
   (*conflict)->nsbreconvclauses = 0;
   (*conflict)->nsbreconvliterals = 0;
   (*conflict)->nsbiterations = 0;
   (*conflict)->npseudocalls = 0;
   (*conflict)->npseudosuccess = 0;
   (*conflict)->npseudoconfclauses = 0;
   (*conflict)->npseudoconfliterals = 0;
   (*conflict)->npseudoreconvclauses = 0;
   (*conflict)->npseudoreconvliterals = 0;

   return SCIP_OKAY;
}

/** frees conflict analysis data for propagation conflicts */
SCIP_RETCODE SCIPconflictFree(
   SCIP_CONFLICT**       conflict            /**< pointer to conflict analysis data */
   )
{
   assert(conflict != NULL);
   assert(*conflict != NULL);
   assert((*conflict)->nclauses == 0);

   SCIPclockFree(&(*conflict)->propanalyzetime);
   SCIPclockFree(&(*conflict)->lpanalyzetime);
   SCIPclockFree(&(*conflict)->sbanalyzetime);
   SCIPclockFree(&(*conflict)->pseudoanalyzetime);
   SCIPpqueueFree(&(*conflict)->binbdchgqueue);
   SCIPpqueueFree(&(*conflict)->nonbinbdchgqueue);
   BMSfreeMemoryArrayNull(&(*conflict)->clauses);
   BMSfreeMemoryArrayNull(&(*conflict)->conflictvars);
   BMSfreeMemory(conflict);

   return SCIP_OKAY;
}

/** initializes the propagation conflict analysis by clearing the conflict candidate queue */
SCIP_RETCODE SCIPconflictInit(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   SCIPdebugMessage("initializing conflict analysis\n");

   /* clear the conflict candidate queues and the conflict set */
   SCIPpqueueClear(conflict->binbdchgqueue);
   SCIPpqueueClear(conflict->nonbinbdchgqueue);
   conflict->nconflictvars = 0;
   conflict->ntmpconflictvars = 0;

   /* increase the conflict set counter, such that variables of new conflict set are labeled with this new counter */
   conflict->count++;
   if( conflict->count == 0 ) /* make sure, 0 is not a valid conflict counter (may happen due to integer overflow) */
      conflict->count = 1;

   return SCIP_OKAY;
}

/** adds given bound change information to the conflict candidate queue */
static
SCIP_RETCODE conflictAddBdchginfo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(conflict != NULL);
   assert(bdchginfo != NULL);

   /* put candidate in the appropriate priority queue */
   if( SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) == SCIP_VARTYPE_BINARY )
   {
      SCIP_CALL( SCIPpqueueInsert(conflict->binbdchgqueue, (void*)bdchginfo) );
   }
   else
   {
      SCIP_CALL( SCIPpqueueInsert(conflict->nonbinbdchgqueue, (void*)bdchginfo) );
   }

   return SCIP_OKAY;
}

/** adds variable's bound to conflict candidate queue */
SCIP_RETCODE SCIPconflictAddBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_VAR*             var,                /**< problem variable */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound that was changed: lower or upper bound */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_Real scalar;
   SCIP_Real constant;

   assert(conflict != NULL);
   assert(stat != NULL);
   assert(var != NULL);

   /* get active problem variable */
   scalar = 1.0;
   constant = 0.0;
   SCIP_CALL( SCIPvarGetProbvarSum(&var, &scalar, &constant) );

   /* we can ignore fixed variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED )
      return SCIP_OKAY;

   assert(!SCIPsetIsZero(set, scalar));

   /* if the scalar of the aggregation is negative, we have to switch the bound type */
   if( scalar < 0.0 )
      boundtype = SCIPboundtypeOpposite(boundtype);

   /* if the variable is multi-aggregated, add the bounds of all aggregation variables */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR )
   {
      SCIP_VAR** vars;
      SCIP_Real* scalars;
      int nvars;
      int i;

      vars = SCIPvarGetMultaggrVars(var);
      scalars = SCIPvarGetMultaggrScalars(var);
      nvars = SCIPvarGetMultaggrNVars(var);
      for( i = 0; i < nvars; ++i )
      {
         SCIP_CALL( SCIPconflictAddBound(conflict, set, stat, vars[i],
               (scalars[i] < 0.0 ? SCIPboundtypeOpposite(boundtype) : boundtype), bdchgidx) );
      }

      return SCIP_OKAY;
   }
   assert(SCIPvarIsActive(var));

   /* get bound change information */
   bdchginfo = SCIPvarGetBdchgInfo(var, boundtype, bdchgidx, FALSE);

   /* if bound of variable was not changed, we can ignore the conflicting bound */
   if( bdchginfo == NULL )
      return SCIP_OKAY;

   SCIPdebugMessage(" -> adding bound <%s> %s %g [status:%d, type:%d, depth:%d, pos:%d, reason:<%s>, info:%d] to candidates\n",
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
   SCIP_CALL( conflictAddBdchginfo(conflict, bdchginfo) );

   /* update the conflict score of the variable */
   SCIPvarIncConflictScore(var, boundtype == SCIP_BOUNDTYPE_LOWER ? SCIP_BRANCHDIR_UPWARDS : SCIP_BRANCHDIR_DOWNWARDS,
      stat->conflictscoreweight);

   return SCIP_OKAY;
}

/** removes and returns next conflict analysis candidate from the candidate queues */
static
SCIP_BDCHGINFO* conflictRemoveCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;

   assert(conflict != NULL);

   bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->nonbinbdchgqueue));
   if( bdchginfo == NULL )
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueRemove(conflict->binbdchgqueue));

   return bdchginfo;
}

/** returns next conflict analysis candidate from the candidate queues without removing it */
static
SCIP_BDCHGINFO* conflictFirstCand(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   SCIP_BDCHGINFO* bdchginfo;

   assert(conflict != NULL);

   bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->nonbinbdchgqueue));
   if( bdchginfo == NULL )
      bdchginfo = (SCIP_BDCHGINFO*)(SCIPpqueueFirst(conflict->binbdchgqueue));

   return bdchginfo;
}

/** identifies the depth, at which the conflict clause should be added:
 *  - at first, the insert depth is increased until the number of remaining unfixed literals in the clause
 *    does not exceed the parameter conflict/minmaxvars anymore
 *  - if the branching rule operates on variables only, and if all branching variables up to a certain
 *    depth level are member of the conflict, the conflict clause can only be violated in the subtree
 *    of the node at that depth, because in all other nodes, at least one of these branching variables
 *    takes a different value, such that the conflict clause is feasible
 *  - if there is at least one branching variable in a node, we assume, that this branching was performed
 *    on variables, and that the siblings of this node are disjunct w.r.t. the branching variables' fixings
 *  - the variables in the conflict set are labeled with the current conflict set counter
 *  - we have to add the conflict clause at least in the valid depth of the initial conflict set,
 *    so we start searching at the first branching after this depth level, i.e. validdepth+1
 */
static
int conflictCalcInsertDepth(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth          /**< minimal depth level at which the initial conflict set is valid */
   )
{
   int currentdepth;
   int insertdepth;

   assert(conflict != NULL);
   assert(tree != NULL);

   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);

   /* the clause must not be inserted prior to its valid depth */
   insertdepth = validdepth;

   /* increase insert depth until the number of remaining unfixed literals in the clause does not exceed the
    * parameter conflict/minmaxvars anymore
    */
   if( set->conf_maxunfixed >= 0 )
   {
      int nunfixedlits;

      nunfixedlits = conflict->nconflictvars + conflict->ntmpconflictvars;
      if( nunfixedlits > set->conf_maxunfixed )
      {
         int* nlits;
         int i;

         SCIP_CALL( SCIPsetAllocBufferArray(set, &nlits, currentdepth+2) );
         BMSclearMemoryArray(nlits, currentdepth+2);

         for( i = 0; i < conflict->nconflictvars + conflict->ntmpconflictvars; ++i )
         {
            int fixdepth;

            /* get the depth level, at which the variable was fixed */
            fixdepth = SCIPvarGetLastBdchgDepth(conflict->conflictvars[i]);
            assert(fixdepth != -1); /* variables fixed in preprocessing shouldn't occur here */
            if( fixdepth == -2 )
            {
               /* the variable is not fixed at the focus node: it was fixed due to diving, probing, or strong branching */
               fixdepth = currentdepth+1;
            }
            nlits[fixdepth]++;
         }
         for( i = 0; i <= currentdepth && nunfixedlits > set->conf_maxunfixed; ++i )
            nunfixedlits -= nlits[i];
         insertdepth = MAX(insertdepth, i-1);

         SCIPsetFreeBufferArray(set, &nlits);
      }
   }

   /* skip additional depth levels where branching on the conflict variables was applied */
   for( insertdepth++; insertdepth <= currentdepth; ++insertdepth )
   {
      SCIP_NODE* node;
      SCIP_VAR* var;
      SCIP_BOUNDCHG* boundchgs;
      SCIP_BOUNDCHGTYPE boundchgtype;
      int nboundchgs;
      int nbranchingvars;
      int b;

      node = tree->path[insertdepth];
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
         boundchgtype = (SCIP_BOUNDCHGTYPE)boundchgs[b].boundchgtype;

         SCIPdebugMessage(" -> depth %d, bound change %d: <%s> %s %g (branching: %d, conflict set: %d)\n",
            insertdepth, b, SCIPvarGetName(var), (SCIP_BOUNDTYPE)boundchgs[b].boundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
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

   /* now, insertdepth is the depth level of the first node, that has non-variable branching or a branching variable
    * not in the conflict set; this means, the siblings of the node may benefit from the conflict clause,
    * and the clause should be added to the node's parent, i.e. at depth level insertdepth-1
    */
   insertdepth--;
   assert(validdepth <= insertdepth && insertdepth <= currentdepth);

   return insertdepth;
}

/** calls the conflict handlers in order to create a conflict clause */
static
SCIP_RETCODE conflictAddClause(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success,            /**< pointer to store whether the conflict set is valid */
   int*                  nliterals           /**< pointer to store the number of literals in the generated clause */
   )
{
   SCIP_BDCHGINFO** bdchginfos;
   int nbdchginfos;
   int currentdepth;
   int focusdepth;
   int insertdepth;
   int i;

   assert(conflict != NULL);
   assert(conflict->nconflictvars >= 0);
   assert(conflict->nconflictvars == 0 || conflict->conflictvars != NULL);
   assert(conflict->ntmpconflictvars == 0);
   assert(SCIPpqueueNElems(conflict->nonbinbdchgqueue) == 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(tree != NULL);
   assert(success != NULL);
   assert(nliterals != NULL);

   *success = FALSE;
   *nliterals = 0;

   /* check, whether local conflicts are allowed */
   if( !set->conf_allowlocal && validdepth > 0 )
      return SCIP_OKAY;

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);
   assert(0 <= validdepth && validdepth <= currentdepth);

   /* temporarily move remaining bound changes from the queue into the conflict set */
   bdchginfos = (SCIP_BDCHGINFO**)SCIPpqueueElems(conflict->binbdchgqueue);
   nbdchginfos = SCIPpqueueNElems(conflict->binbdchgqueue);
   SCIPdebugMessage("adding %d variables from the queue as temporary conflict variables\n", nbdchginfos);
   for( i = 0; i < nbdchginfos; ++i )
   {
      SCIP_CALL( conflictAddConflictVar(conflict, blkmem, set, stat, SCIPbdchginfoGetVar(bdchginfos[i]), FALSE, TRUE) );
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

   /* calculate the depth, at which the clause should be inserted */
   insertdepth = conflictCalcInsertDepth(conflict, set, tree, validdepth);
   assert(validdepth <= insertdepth && insertdepth <= currentdepth);
   SCIPdebugMessage(" -> conflict with %d literals found at depth %d is active in depth %d and valid in depth %d\n",
      conflict->nconflictvars + conflict->ntmpconflictvars, currentdepth, insertdepth, validdepth);

   /* if all branching variables are in the conflict set, the conflict clause is of no use;
    * don't use clauses that are only valid in the probing path but not in the problem tree
    */
   if( insertdepth < currentdepth && insertdepth <= focusdepth )
   {
      *nliterals = conflict->nconflictvars + conflict->ntmpconflictvars;
      SCIPdebugMessage(" -> final conflict clause has %d literals\n", *nliterals);

      /* check conflict clause on debugging solution */
      SCIP_CALL( SCIPdebugCheckConflict(blkmem, set, tree->path[validdepth],
            conflict->conflictvars, *nliterals) ); /*lint !e506 !e774*/

      /* insert clause to the clause storage */
      SCIP_CALL( conflictInsertClause(conflict, blkmem, set, conflict->conflictvars, *nliterals,
            validdepth, insertdepth) );
      *success = TRUE;
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
SCIP_RETCODE conflictResolveBound(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change to resolve */
   int*                  validdepth,         /**< pointer to update the minimal depth level at which the conflict is valid */
   SCIP_Bool*            resolved            /**< pointer to store whether the bound change was resolved */
   )
{
   SCIP_VAR* actvar;
   SCIP_CONS* infercons;
   SCIP_PROP* inferprop;
   SCIP_RESULT result;

   assert(conflict != NULL);
   assert(validdepth != NULL);
   assert(resolved != NULL);

   *resolved = FALSE;

   actvar = SCIPbdchginfoGetVar(bdchginfo);
   assert(actvar != NULL);
   assert(SCIPvarIsActive(actvar));

#ifdef SCIP_DEBUG
   {
      int v;

      SCIPdebugMessage("processing next conflicting bound (depth: %d, valid depth: %d, bdchgtype: %s [%s], vartype: %d): [<%s> %s %g]\n",
         SCIPbdchginfoGetDepth(bdchginfo), *validdepth,
         SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "branch"
         : SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER ? "cons" : "prop",
         SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_BRANCHING ? "-"
         : SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
         ? SCIPconsGetName(SCIPbdchginfoGetInferCons(bdchginfo))
         : SCIPbdchginfoGetInferProp(bdchginfo) == NULL ? "-"
         : SCIPpropGetName(SCIPbdchginfoGetInferProp(bdchginfo)),
         SCIPvarGetType(actvar), SCIPvarGetName(actvar),
         SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
         SCIPbdchginfoGetNewbound(bdchginfo));
      SCIPdebugMessage(" - conflict set             :");
      for( v = 0; v < conflict->nconflictvars; ++v )
         SCIPdebugPrintf(" <%s>[%d]", SCIPvarGetName(conflict->conflictvars[v]),
            SCIPvarGetLastBdchgDepth(conflict->conflictvars[v]));
      SCIPdebugPrintf("\n");
      SCIPdebugMessage(" - candidate queue (non-bin):");
      for( v = 0; v < SCIPpqueueNElems(conflict->nonbinbdchgqueue); ++v )
      {
         SCIP_BDCHGINFO* info = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->nonbinbdchgqueue)[v]);
         SCIPdebugPrintf(" [%d:<%s> %s %g]", SCIPbdchginfoGetDepth(info), SCIPvarGetName(SCIPbdchginfoGetVar(info)),
            SCIPbdchginfoGetBoundtype(info) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(info));
      }
      SCIPdebugPrintf("\n");
      SCIPdebugMessage(" - candidate queue (bin)    :");
      for( v = 0; v < SCIPpqueueNElems(conflict->binbdchgqueue); ++v )
      {
         SCIP_BDCHGINFO* info = (SCIP_BDCHGINFO*)(SCIPpqueueElems(conflict->binbdchgqueue)[v]);
         SCIPdebugPrintf(" [%d:<%s> %s %g]", SCIPbdchginfoGetDepth(info), SCIPvarGetName(SCIPbdchginfoGetVar(info)),
            SCIPbdchginfoGetBoundtype(info) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(info));
      }
      SCIPdebugPrintf("\n");
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
         || SCIPconsGetValidDepth(infercons) <= *validdepth )
      {
         SCIP_VAR* infervar;
         int inferinfo;
         SCIP_BOUNDTYPE inferboundtype;
         SCIP_BDCHGIDX* bdchgidx;

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

         SCIPdebugMessage("resolving bound <%s> %s %g [status:%d, type:%d, depth:%d, pos:%d]: <%s> %s %g [cons:<%s>(%s), info:%d]\n",
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

         SCIP_CALL( SCIPconsResolvePropagation(infercons, set, infervar, inferinfo, inferboundtype,
               bdchgidx, &result) );
         *resolved = (result == SCIP_SUCCESS);

         if( *resolved && SCIPconsIsLocal(infercons) )
         {
            int vdepth;

            vdepth = SCIPconsGetValidDepth(infercons);
            *validdepth = MAX(*validdepth, vdepth);
         }
      }
      break;

   case SCIP_BOUNDCHGTYPE_PROPINFER:
      inferprop = SCIPbdchginfoGetInferProp(bdchginfo);
      if( inferprop != NULL )
      {
         SCIP_VAR* infervar;
         int inferinfo;
         SCIP_BOUNDTYPE inferboundtype;
         SCIP_BDCHGIDX* bdchgidx;

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

         SCIPdebugMessage("resolving bound <%s> %s %g [status:%d, depth:%d, pos:%d]: <%s> %s %g [prop:<%s>, info:%d]\n",
            SCIPvarGetName(actvar),
            SCIPbdchginfoGetBoundtype(bdchginfo) == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPbdchginfoGetNewbound(bdchginfo),
            SCIPvarGetStatus(actvar), SCIPbdchginfoGetDepth(bdchginfo), SCIPbdchginfoGetPos(bdchginfo),
            SCIPvarGetName(infervar),
            inferboundtype == SCIP_BOUNDTYPE_LOWER ? ">=" : "<=",
            SCIPvarGetBdAtIndex(infervar, inferboundtype, bdchgidx, TRUE),
            SCIPpropGetName(inferprop), inferinfo);

         SCIP_CALL( SCIPpropResolvePropagation(inferprop, set, infervar, inferinfo, inferboundtype,
               bdchgidx, &result) );
         *resolved = (result == SCIP_SUCCESS);
      }
      break;

   case SCIP_BOUNDCHGTYPE_BRANCHING:
      assert(!(*resolved));
      break;

   default:
      SCIPerrorMessage("invalid bound change type <%d>\n", SCIPbdchginfoGetChgtype(bdchginfo));
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** if only one conflicting bound change of the last depth level was used, and if this can be resolved,
 *  creates GRASP-like reconvergence clauses in the conflict graph up to the branching variable of the depth level
 */
static
SCIP_RETCODE conflictCreateReconvergenceClauses(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_BDCHGINFO*       firstuip,           /**< first UIP of conflict graph */
   int*                  nreconvclauses,     /**< pointer to store the number of generated reconvergence clauses */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence clauses */
   )
{
   SCIP_BDCHGINFO* uip;
   int firstuipdepth;
   int maxvaliddepth;

   assert(conflict != NULL);
   assert(firstuip != NULL);
   assert(nreconvclauses != NULL);
   assert(nreconvliterals != NULL);

   firstuipdepth = SCIPbdchginfoGetDepth(firstuip);

   /* check, whether local conflicts are allowed */
   maxvaliddepth = (set->conf_allowlocal ? firstuipdepth-1 : 0);

   /* for each succeeding UIP pair of the last depth level, create one reconvergence clause */
   uip = firstuip;
   while( uip != NULL && SCIPbdchginfoGetDepth(uip) == SCIPbdchginfoGetDepth(firstuip) )
   {
      SCIP_BDCHGINFO* bdchginfo;
      SCIP_BDCHGINFO* nextuip;
      SCIP_VAR* var;
      int nresolutions;

      SCIPdebugMessage("creating reconvergence clause for UIP <%s> in depth %d\n",
         SCIPvarGetName(SCIPbdchginfoGetVar(uip)), SCIPbdchginfoGetDepth(uip));

      /* initialize conflict data */
      SCIP_CALL( SCIPconflictInit(conflict) );

      /* put the variable of first UIP into the conflict set, using the literal that is currently fixed to TRUE */
      var = SCIPbdchginfoGetVar(uip);
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      SCIP_CALL( SCIPvarNegate(var, blkmem, set, stat, &var) );
      SCIP_CALL( conflictAddConflictVar(conflict, blkmem, set, stat, var, TRUE, FALSE) );

      /* put current UIP into priority queue */
      SCIP_CALL( conflictAddBdchginfo(conflict, uip) );

      /* resolve the queue until the next UIP is reached */
      bdchginfo = conflictFirstCand(conflict);
      nextuip = NULL;
      nresolutions = 0;
      while( bdchginfo != NULL && validdepth <= maxvaliddepth )
      {
         SCIP_BDCHGINFO* nextbdchginfo;
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
            SCIP_VAR* actvar;
            SCIP_Bool resolved;

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
               SCIP_CALL( conflictResolveBound(conflict, set, bdchginfo, &validdepth, &resolved) );
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
               SCIP_CALL( conflictAddConflictVar(conflict, blkmem, set, stat, actvar, FALSE, FALSE) );
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
      if( nextuip != NULL && nresolutions >= 2 && bdchginfo == NULL && validdepth <= maxvaliddepth )
      {
         int nlits;
         SCIP_Bool success;

         assert(SCIPbdchginfoGetDepth(nextuip) == SCIPbdchginfoGetDepth(uip));
         assert(SCIPvarGetType(SCIPbdchginfoGetVar(nextuip)) == SCIP_VARTYPE_BINARY);

         SCIPdebugMessage("creating reconvergence clause from UIP <%s> to UIP <%s> in depth %d with %d literals after %d resolutions\n",
            SCIPvarGetName(SCIPbdchginfoGetVar(uip)), SCIPvarGetName(SCIPbdchginfoGetVar(nextuip)),
            SCIPbdchginfoGetDepth(uip), conflict->nconflictvars, nresolutions);

         /* call the conflict handlers to create a conflict clause */
         SCIP_CALL( conflictAddClause(conflict, blkmem, set, stat, tree, validdepth, &success, &nlits) );
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
SCIP_RETCODE conflictAnalyze(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool             mustresolve,        /**< should the conflict set only be used, if a resolution was applied? */
   int*                  nclauses,           /**< pointer to store the number of generated conflict clauses */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict clauses */
   int*                  nreconvclauses,     /**< pointer to store the number of generated reconvergence clauses */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence clauses */
   )
{
   SCIP_BDCHGINFO* bdchginfo;
   SCIP_BDCHGINFO* firstuip;
   int focusdepth;
   int currentdepth;
   int maxvaliddepth;
   int resolvedepth;
   int nresolutions;
   int lastclausenresolutions;
   int lastclauseresoldepth;

   assert(conflict != NULL);
   assert(conflict->nconflictvars >= 0);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(0 <= validdepth && validdepth <= SCIPtreeGetCurrentDepth(tree));
   assert(nclauses != NULL);
   assert(nliterals != NULL);
   assert(nreconvclauses != NULL);
   assert(nreconvliterals != NULL);

   focusdepth = SCIPtreeGetFocusDepth(tree);
   currentdepth = SCIPtreeGetCurrentDepth(tree);
   assert(currentdepth == tree->pathlen-1);
   assert(focusdepth <= currentdepth);

   resolvedepth = ((set->conf_fuiplevels >= 0 && set->conf_fuiplevels <= currentdepth)
      ? currentdepth - set->conf_fuiplevels + 1 : 0);
   assert(0 <= resolvedepth && resolvedepth <= currentdepth + 1);

   /* if we must resolve at least one bound change, find the first UIP at least in the last depth level */
   if( mustresolve )
      resolvedepth = MIN(resolvedepth, currentdepth);

   SCIPdebugMessage("analyzing conflict with %d+%d conflict candidates and starting conflict set of size %d in depth %d (resolvedepth=%d)\n",
      SCIPpqueueNElems(conflict->binbdchgqueue), SCIPpqueueNElems(conflict->nonbinbdchgqueue),
      conflict->nconflictvars, currentdepth, resolvedepth);

   *nclauses = 0;
   *nliterals = 0;
   *nreconvclauses = 0;
   *nreconvliterals = 0;

   /* check, whether local conflicts are allowed; however, don't generate clauses that are only valid in the
    * probing path and not in the problem tree (i.e. that exceed the focusdepth)
    */
   maxvaliddepth = (set->conf_allowlocal ? MIN(currentdepth-1, focusdepth) : 0);
   if( validdepth > maxvaliddepth )
      return SCIP_OKAY;

   /* process all bound changes in the conflict candidate queue */
   nresolutions = 0;
   lastclausenresolutions = (mustresolve ? 0 : -1);
   lastclauseresoldepth = (mustresolve ? currentdepth : INT_MAX);
   bdchginfo = conflictFirstCand(conflict);
   firstuip = NULL;
   while( bdchginfo != NULL && validdepth <= maxvaliddepth )
   {
      SCIP_BDCHGINFO* nextbdchginfo;
      int bdchgdepth;

      /* resolve next bound change in queue */
      bdchgdepth = SCIPbdchginfoGetDepth(bdchginfo);
      assert(0 <= bdchgdepth && bdchgdepth <= currentdepth);
      assert(SCIPvarIsActive(SCIPbdchginfoGetVar(bdchginfo)));
      assert(bdchgdepth < tree->pathlen);
      assert(tree->path[bdchgdepth] != NULL);
      assert(tree->path[bdchgdepth]->domchg != NULL);
      assert(SCIPbdchginfoGetPos(bdchginfo) < (int)tree->path[bdchgdepth]->domchg->domchgbound.nboundchgs);
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].var
         == SCIPbdchginfoGetVar(bdchginfo));
      assert(tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].newbound
         == SCIPbdchginfoGetNewbound(bdchginfo)); /*lint !e777*/
      assert((SCIP_BOUNDTYPE)tree->path[bdchgdepth]->domchg->domchgbound.boundchgs[SCIPbdchginfoGetPos(bdchginfo)].boundtype
         == SCIPbdchginfoGetBoundtype(bdchginfo));

      /* create intermediate conflict clause */
      if( nresolutions > lastclausenresolutions
         && (set->conf_interclauses == -1 || *nclauses < set->conf_interclauses)
         && validdepth <= maxvaliddepth
         && SCIPpqueueNElems(conflict->nonbinbdchgqueue) == 0
         && bdchgdepth < lastclauseresoldepth )
      {
         int nlits;
         SCIP_Bool success;

         /* call the conflict handlers to create a conflict clause */
         SCIPdebugMessage("creating intermediate clause after %d resolutions up to depth %d (valid at depth %d): %d conflict vars, %d vars in queue\n",
            nresolutions, bdchgdepth, validdepth, conflict->nconflictvars, SCIPpqueueNElems(conflict->binbdchgqueue));
         SCIP_CALL( conflictAddClause(conflict, blkmem, set, stat, tree, validdepth, &success, &nlits) );
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
         SCIP_VAR* actvar;
         SCIP_Bool resolved;

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
            SCIP_CALL( conflictResolveBound(conflict, set, bdchginfo, &validdepth, &resolved) );
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

            SCIPdebugMessage("couldn't resolve bound change on <%s> -> new valid depth: %d\n",
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
            SCIP_CALL( conflictAddConflictVar(conflict, blkmem, set, stat, actvar, FALSE, FALSE) );
         }
      }

      /* get next conflicting bound from the conflict candidate queue (this needs not to be nextbdchginfo, because
       * due to resolving the bound changes, a non-binary variable could be added to the queue which must be
       * resolved before nextbdchginfo
       */
      bdchginfo = conflictFirstCand(conflict);
   }

   /* check, if a valid conflict set was found */
   if( bdchginfo == NULL
      && nresolutions > lastclausenresolutions
      && validdepth <= maxvaliddepth
      && (!mustresolve || nresolutions > 0 || conflict->nconflictvars == 0) )
   {
      int nlits;
      SCIP_Bool success;

      /* call the conflict handlers to create a conflict clause */
      SCIP_CALL( conflictAddClause(conflict, blkmem, set, stat, tree, validdepth, &success, &nlits) );
      if( success )
      {
         (*nclauses)++;
         (*nliterals) += nlits;
      }
   }

   /* produce reconvergence clauses defined by succeeding UIP's of the last depth level */
   if( set->conf_reconvclauses && firstuip != NULL && SCIPbdchginfoHasInferenceReason(firstuip) )
   {
      SCIP_CALL( conflictCreateReconvergenceClauses(conflict, blkmem, set, stat, tree, validdepth, firstuip,
            nreconvclauses, nreconvliterals) );
   }

   /* increase the conflict score weight for history updates of future conflict reasons */
   if( stat->nnodes > stat->lastconflictnode )
   {
      assert(0.0 < set->conf_scorefac && set->conf_scorefac <= 1.0);
      stat->conflictscoreweight /= set->conf_scorefac;
      assert(stat->conflictscoreweight > 0.0);

      /* if the conflict score for the next conflict exceeds 1000.0, rescale all history conflict scores */
      if( stat->conflictscoreweight >= 1000.0 )
      {
         int v;

         for( v = 0; v < prob->nvars; ++v )
            SCIPvarScaleConflictScores(prob->vars[v], 1.0/stat->conflictscoreweight);
         stat->conflictscoreweight = 1.0;
      }
      stat->lastconflictnode = stat->nnodes;
   }

   return SCIP_OKAY;
}

/** analyzes conflicting bound changes that were added with calls to SCIPconflictAddBound(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set;
 *  updates statistics for propagation conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyze(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int                   validdepth,         /**< minimal depth level at which the initial conflict set is valid */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
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

   /* check, if the conflict clause will get too large with high probability */
   if( conflict->nconflictvars + SCIPpqueueNElems(conflict->binbdchgqueue)
      + SCIPpqueueNElems(conflict->nonbinbdchgqueue) >= 2*conflictCalcMaxsize(set, prob) )
      return SCIP_OKAY;

   SCIPdebugMessage("analyzing conflict after infeasible propagation in depth %d\n", SCIPtreeGetCurrentDepth(tree));

   /* start timing */
   SCIPclockStart(conflict->propanalyzetime, set);

   conflict->npropcalls++;

   /* analyze the conflict set, and create a conflict constraint on success */
   SCIP_CALL( conflictAnalyze(conflict, blkmem, set, stat, prob, tree, validdepth, TRUE,
         &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
   conflict->npropsuccess += (nclauses > 0);
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
SCIP_Real SCIPconflictGetPropTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->propanalyzetime);
}

/** gets number of calls to propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropcalls;
}

/** gets number of calls to propagation conflict analysis that yield at least one conflict clause */
SCIP_Longint SCIPconflictGetNPropSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropsuccess;
}

/** gets number of conflict clauses detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfclauses;
}

/** gets total number of literals in conflict clauses created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconfliterals;
}

/** gets number of reconvergence clauses detected in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvclauses;
}

/** gets total number of literals in reconvergence clauses created in propagation conflict analysis */
SCIP_Longint SCIPconflictGetNPropReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropreconvliterals;
}




/*
 * Infeasible LP Conflict Analysis
 */

/** returns whether bound change has a valid reason that can be resolved in conflict analysis */
static
SCIP_Bool bdchginfoIsResolvable(
   SCIP_BDCHGINFO*       bdchginfo           /**< bound change information */
   )
{
   assert(bdchginfo != NULL);

   return (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_CONSINFER
      || (SCIPbdchginfoGetChgtype(bdchginfo) == SCIP_BOUNDCHGTYPE_PROPINFER
         && SCIPbdchginfoGetInferProp(bdchginfo) != NULL));
}

/** returns whether bound change information can be used in conflict analysis:
 *  - if variable is binary, or
 *  - if bound change information has valid inference data
 */
static
void bdchginfoGetUsability(
   SCIP_BDCHGINFO*       bdchginfo,          /**< bound change information */
   SCIP_Bool*            resolvable,         /**< has bound change a valid reason that can be resolved in conflict analysis? */
   SCIP_Bool*            usable              /**< is bound change binary or resolveable? */
   )
{
   assert(bdchginfo != NULL);
   assert(resolvable != NULL);
   assert(usable != NULL);

   *resolvable = bdchginfoIsResolvable(bdchginfo);
   *usable = *resolvable || (SCIPvarGetType(SCIPbdchginfoGetVar(bdchginfo)) == SCIP_VARTYPE_BINARY);
}

/** ensures, that side change arrays can store at least num entries */
static
SCIP_RETCODE ensureSidechgsSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int**                 sidechginds,        /**< pointer to side change index array */
   SCIP_Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   SCIP_Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   SCIP_Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   SCIP_Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*                  sidechgssize,       /**< pointer to size of side change arrays */
   int                   num                 /**< minimal number of entries to be able to store in side change arrays */
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
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechginds, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldrhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewrhss, newsize) );
      *sidechgssize = newsize;
   }
   assert(num <= *sidechgssize);

   return SCIP_OKAY;
}

/** adds removal of row's side to side change arrays; finite sides are only replaced by near infinite sides, such
 *  that the row's sense in the LP solver is not changed
 */
static
SCIP_RETCODE addSideRemoval(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_ROW*             row,                /**< LP row to change the sides for */
   SCIP_Real             lpiinfinity,        /**< value treated as infinity in LP solver */
   int**                 sidechginds,        /**< pointer to side change index array */
   SCIP_Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   SCIP_Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   SCIP_Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   SCIP_Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*                  sidechgssize,       /**< pointer to size of side change arrays */
   int*                  nsidechgs           /**< pointer to number of used slots in side change arrays */
   )
{
   SCIP_Real lhs;
   SCIP_Real rhs;
   SCIP_Real constant;

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
   SCIP_CALL( ensureSidechgsSize(set, sidechginds, sidechgoldlhss, sidechgoldrhss, sidechgnewlhss, sidechgnewrhss,
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
SCIP_RETCODE ensureBdchgsSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   int**                 sidechginds,        /**< pointer to side change index array */
   SCIP_Real**           sidechgoldlhss,     /**< pointer to side change old left hand sides array */
   SCIP_Real**           sidechgoldrhss,     /**< pointer to side change old right hand sides array */
   SCIP_Real**           sidechgnewlhss,     /**< pointer to side change new left hand sides array */
   SCIP_Real**           sidechgnewrhss,     /**< pointer to side change new right hand sides array */
   int*                  sidechgssize,       /**< pointer to size of side change arrays */
   int                   num                 /**< minimal number of entries to be able to store in side change arrays */
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
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechginds, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgoldrhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewlhss, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, sidechgnewrhss, newsize) );
      *sidechgssize = newsize;
   }
   assert(num <= *sidechgssize);

   return SCIP_OKAY;
}

/** inserts variable's new bounds into bound change arrays */
static
SCIP_RETCODE addBdchg(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR*             var,                /**< variable to change the LP bounds for */
   SCIP_Real             newlb,              /**< new lower bound */
   SCIP_Real             newub,              /**< new upper bound */
   int**                 bdchginds,          /**< pointer to bound change index array */
   SCIP_Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array */
   SCIP_Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array */
   SCIP_Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array */
   SCIP_Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array */
   int*                  bdchgssize,         /**< pointer to size of bound change arrays */
   int*                  nbdchgs             /**< pointer to number of used slots in bound change arrays */
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
      SCIP_COL* col;
      int c;

      col = SCIPvarGetCol(var);
      c = SCIPcolGetLPPos(col);
      if( c >= 0 )
      {
         SCIP_CALL( ensureBdchgsSize(set, bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs,
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
SCIP_RETCODE ensureCandsSize(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_VAR***           cands,              /**< pointer to candidate array */
   SCIP_Real**           candscores,         /**< pointer to candidate score array */
   SCIP_Real**           newbounds,          /**< pointer to candidate new bounds array */
   SCIP_Real**           proofactdeltas,     /**< pointer to candidate proof delta array */
   int*                  candssize,          /**< pointer to size of array */
   int                   num                 /**< minimal number of candidates to store in array */
   )
{
   assert(cands != NULL);
   assert(candssize != NULL);

   if( num > *candssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_CALL( SCIPsetReallocBufferArray(set, cands, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, candscores, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, newbounds, newsize) );
      SCIP_CALL( SCIPsetReallocBufferArray(set, proofactdeltas, newsize) );
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
SCIP_RETCODE addCand(
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_VAR*             var,                /**< variable to add to candidate array */
   int                   lbchginfopos,       /**< positions of currently active lower bound change infos in variable's array */
   int                   ubchginfopos,       /**< positions of currently active upper bound change infos in variable's array */
   SCIP_Real             proofcoef,          /**< coefficient of variable in infeasibility/bound proof */
   SCIP_Real             prooflhs,           /**< left hand side of infeasibility/bound proof */
   SCIP_Real             proofact,           /**< activity of infeasibility/bound proof row */
   SCIP_VAR***           cands,              /**< pointer to candidate array for undoing bound changes */
   SCIP_Real**           candscores,         /**< pointer to candidate score array for undoing bound changes */
   SCIP_Real**           newbounds,          /**< pointer to candidate new bounds array for undoing bound changes */
   SCIP_Real**           proofactdeltas,     /**< pointer to proof activity increase array for undoing bound changes */
   int*                  candssize,          /**< pointer to size of cands arrays */
   int*                  ncands,             /**< pointer to count number of candidates in bound change list */
   int                   firstcand,          /**< position of first unprocessed bound change candidate */
   int*                  insertpos           /**< pointer to store insertion position, or -1 if not inserted */
   )
{
   SCIP_Real oldbound;
   SCIP_Real newbound;
   SCIP_Real proofactdelta;
   SCIP_Real score;
   int depth;
   int i;
   SCIP_Bool resolvable;
   SCIP_Bool usable;

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
         depth = currentdepth+1;
         resolvable = FALSE;
         usable = TRUE;
      }
      else
      {
         /* current bound is the result of a local bound change */
         bdchginfoGetUsability(&var->ubchginfos[ubchginfopos], &resolvable, &usable);
         depth = var->ubchginfos[ubchginfopos].bdchgidx.depth;
         oldbound = var->ubchginfos[ubchginfopos].newbound;
         newbound = var->ubchginfos[ubchginfopos].oldbound;
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
         depth = currentdepth+1;
         resolvable = FALSE;
         usable = TRUE;
      }
      else
      {
         /* current bound is the result of a local bound change */
         bdchginfoGetUsability(&var->lbchginfos[lbchginfopos], &resolvable, &usable);
         depth = var->lbchginfos[lbchginfopos].bdchgidx.depth;
         oldbound = var->lbchginfos[lbchginfopos].newbound;
         newbound = var->lbchginfos[lbchginfopos].oldbound;
      }
   }
   else
      return SCIP_OKAY;

   /* calculate the increase in the proof's activity */
   proofactdelta = (newbound - oldbound)*proofcoef;
   assert(proofactdelta > 0.0);

   /* if the bound change is not usable in conflict analysis, we have to undo it */
   if( !usable )
      score = (depth+1) * SCIPsetInfinity(set);
   else
   {
      /* calculate score for undoing the bound change */
      score = 1.0 - proofactdelta/(prooflhs - proofact);
      score = MAX(score, 0.0);
      score += (SCIP_Real)(depth+1)/(SCIP_Real)(currentdepth+1);
      if( !resolvable )
         score += 10.0;
   }

   /* get enough memory to store new candidate */
   SCIP_CALL( ensureCandsSize(set, cands, candscores, newbounds, proofactdeltas, candssize, (*ncands)+1) );
   assert(*cands != NULL);
   assert(*candscores != NULL);
   assert(*newbounds != NULL);
   assert(*proofactdeltas != NULL);

   SCIPdebugMessage(" -> local <%s> %s %g, relax <%s> %s %g, dpt=%d, resolve=%d, use=%d, delta=%g, score=%g\n",
      SCIPvarGetName(var), proofcoef > 0.0 ? "<=" : ">=", oldbound,
      SCIPvarGetName(var), proofcoef > 0.0 ? "<=" : ">=", newbound,
      depth, resolvable, usable, proofactdelta, score);

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
SCIP_RETCODE undoBdchgsProof(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            proofcoefs,         /**< coefficients in infeasibility proof */
   SCIP_Real             prooflhs,           /**< left hand side of proof */
   SCIP_Real             proofact,           /**< current activity of proof */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change infos in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change infos in variables' arrays */
   int**                 bdchginds,          /**< pointer to bound change index array, or NULL */
   SCIP_Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   SCIP_Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   SCIP_Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   SCIP_Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*                  bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*                  nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   SCIP_Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR** cands;
   SCIP_Real* candscores;
   SCIP_Real* newbounds;
   SCIP_Real* proofactdeltas;
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
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cands, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &candscores, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &newbounds, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &proofactdeltas, nvars) );
   ncands = 0;
   candssize = nvars;
   for( v = 0; v < nvars; ++v )
   {
      /* ignore variables already relaxed to global bounds */
      if( lbchginfoposs[v] == -1 && ubchginfoposs[v] == -1 )
         continue;

      /* add variable to candidate list */
      SCIP_CALL( addCand(set, currentdepth, vars[v], lbchginfoposs[v], ubchginfoposs[v], proofcoefs[v],
            prooflhs, proofact, &cands, &candscores, &newbounds, &proofactdeltas, &candssize, &ncands, 0, &pos) );
      assert(-1 <= pos && pos < ncands);
      if( pos == -1 )
      {
         assert(SCIPsetIsZero(set, proofcoefs[v])
            || (SCIPsetIsPositive(set, proofcoefs[v]) && ubchginfoposs[v] == -1)
            || (SCIPsetIsNegative(set, proofcoefs[v]) && lbchginfoposs[v] == -1));

         /* variable can be relaxed to global bounds */
         SCIPdebugMessage(" -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g\n",
            SCIPvarGetName(vars[v]), curvarlbs[v], curvarubs[v], SCIPvarGetLbGlobal(vars[v]), SCIPvarGetUbGlobal(vars[v]),
            proofcoefs[v], prooflhs, proofact);
         curvarlbs[v] = SCIPvarGetLbGlobal(vars[v]);
         curvarubs[v] = SCIPvarGetUbGlobal(vars[v]);
         lbchginfoposs[v] = -1;
         ubchginfoposs[v] = -1;
         if( nbdchgs != NULL )
         {
            SCIP_CALL( addBdchg(set, vars[v], curvarlbs[v], curvarubs[v],
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

         SCIPdebugMessage(" -> relaxing variable <%s>[%g,%g] to [%g,%g]: proofcoef=%g, %g <= %g + %g\n",
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
            SCIP_CALL( addBdchg(set, cands[i], curvarlbs[v], curvarubs[v],
                  bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs) );
         }
         proofact += proofactdeltas[i];
         if( resolve != NULL )
            *resolve = TRUE;

         /* insert the new local bound of the variable into the candidate list */
         SCIP_CALL( addCand(set, currentdepth, cands[i], lbchginfoposs[v], ubchginfoposs[v], proofcoefs[v],
               prooflhs, proofact, &cands, &candscores, &newbounds, &proofactdeltas, &candssize, &ncands, i+1, &pos) );
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
SCIP_RETCODE undoBdchgsDualfarkas(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change infos in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change infos in variables' arrays */
   int**                 bdchginds,          /**< pointer to bound change index array, or NULL */
   SCIP_Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   SCIP_Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   SCIP_Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   SCIP_Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*                  bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*                  nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   SCIP_Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   SCIP_Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again */
   )
{
   SCIP_LPI* lpi;
   SCIP_ROW** rows;
   SCIP_VAR** vars;
   SCIP_ROW* row;
   SCIP_VAR* var;
   SCIP_Real* dualfarkas;
   SCIP_Real* farkascoefs;
   SCIP_Real farkaslhs;
   SCIP_Real farkasact;
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

   SCIPdebugMessage("undoing bound changes in infeasible LP: cutoff=%g, depth=%d\n",
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
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualfarkas, nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &farkascoefs, nvars) );

   /* get dual farkas values of rows */
   SCIP_CALL( SCIPlpiGetDualfarkas(lpi, dualfarkas) );

   /* calculate the farkas row */
   BMSclearMemoryArray(farkascoefs, nvars);
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
            SCIP_Real lpilhs;
            SCIP_Real lpirhs;

            SCIP_CALL( SCIPlpiGetSides(lpi, r, r, &lpilhs, &lpirhs) );
            assert((SCIPsetIsInfinity(set, -lpilhs) && SCIPsetIsInfinity(set, -row->lhs))
               || SCIPsetIsEQ(set, lpilhs, row->lhs - row->constant));
            assert((SCIPsetIsInfinity(set, lpirhs) && SCIPsetIsInfinity(set, row->rhs))
               || SCIPsetIsEQ(set, lpirhs, row->rhs - row->constant));
         }
#endif

         /* add row side to farkas row lhs: dualfarkas > 0 -> lhs, dualfarkas < 0 -> rhs */
         if( dualfarkas[r] > 0.0 )
         {
            /* check if sign of dual farkas value is valid */
            if( SCIPsetIsInfinity(set, -row->lhs) )
               continue;
            farkaslhs += dualfarkas[r] * (row->lhs - row->constant);
         }
         else
         {
            /* check if sign of dual farkas value is valid */
            if( SCIPsetIsInfinity(set, row->rhs) )
               continue;
            farkaslhs += dualfarkas[r] * (row->rhs - row->constant);
         }
         SCIPdebugMessage(" -> farkaslhs: %g<%s>[%g,%g] -> %g\n", dualfarkas[r], SCIProwGetName(row),
            row->lhs - row->constant, row->rhs - row->constant, farkaslhs);

         /* add row coefficients to farkas row */
         for( i = 0; i < row->len; ++i )
         {
            v = SCIPvarGetProbindex(SCIPcolGetVar(row->cols[i]));
            assert(0 <= v && v < nvars);
            farkascoefs[v] += dualfarkas[r] * row->vals[i];
         }
      }
#ifdef SCIP_DEBUG
      else if( !SCIPsetIsZero(set, dualfarkas[r]) )
      {
         SCIPdebugMessage(" -> ignoring %s row <%s> with dual farkas value %.10f (lhs=%g, rhs=%g)\n",
            row->local ? "local" : "global", SCIProwGetName(row), dualfarkas[r],
            row->lhs - row->constant, row->rhs - row->constant);
      }
#endif
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
            || !SCIPsetIsPositive(set, SCIPvarGetUbLP(var)));
         if( SCIPsetIsInfinity(set, curvarubs[v]) )
            goto TERMINATE;
         farkasact += farkascoefs[v] * curvarubs[v];
         SCIPdebugMessage(" -> farkasact: %g<%s>[%g,%g] -> %g\n", farkascoefs[v], SCIPvarGetName(var),
            curvarlbs[v], curvarubs[v], farkasact);
      }
      else
      {
         assert((SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN && SCIPcolGetLPPos(SCIPvarGetCol(var)) >= 0)
            || !SCIPsetIsNegative(set, SCIPvarGetLbLP(var)));
         if( SCIPsetIsInfinity(set, -curvarlbs[v]) )
            goto TERMINATE;
         farkasact += farkascoefs[v] * curvarlbs[v];
         SCIPdebugMessage(" -> farkasact: %g<%s>[%g,%g] -> %g\n", farkascoefs[v], SCIPvarGetName(var),
            curvarlbs[v], curvarubs[v], farkasact);
      }
   }
   SCIPdebugMessage(" -> farkaslhs=%g, farkasact=%g\n", farkaslhs, farkasact);

   /* check, if the farkas row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, farkaslhs, farkasact) )
   {
      /* undo bound changes while keeping the infeasibility proof valid */
      SCIP_CALL( undoBdchgsProof(set, prob, currentdepth, farkascoefs, farkaslhs, farkasact,
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
            bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, bdchgssize, nbdchgs, resolve) );

      *valid = TRUE;

      /* resolving does not make sense: the old dual ray is still valid -> resolving will not change the solution */
      *resolve = FALSE;
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
SCIP_RETCODE undoBdchgsDualsol(
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_LP*              lp,                 /**< LP data */
   int                   currentdepth,       /**< current depth in the tree */
   SCIP_Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   SCIP_Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change infos in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change infos in variables' arrays */
   int**                 bdchginds,          /**< pointer to bound change index array, or NULL */
   SCIP_Real**           bdchgoldlbs,        /**< pointer to bound change old lower bounds array, or NULL */
   SCIP_Real**           bdchgoldubs,        /**< pointer to bound change old upper bounds array, or NULL */
   SCIP_Real**           bdchgnewlbs,        /**< pointer to bound change new lower bounds array, or NULL */
   SCIP_Real**           bdchgnewubs,        /**< pointer to bound change new upper bounds array, or NULL */
   int*                  bdchgssize,         /**< pointer to size of bound change arrays, or NULL */
   int*                  nbdchgs,            /**< pointer to number of used slots in bound change arrays, or NULL */
   SCIP_Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   SCIP_Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again */
   )
{
   SCIP_LPI* lpi;
   SCIP_ROW** rows;
   SCIP_VAR** vars;
   SCIP_ROW* row;
   SCIP_VAR* var;
   SCIP_Real* primsols;
   SCIP_Real* dualsols;
   SCIP_Real* redcosts;
   SCIP_Real* dualcoefs;
   SCIP_Real* varredcosts;
   SCIP_Real duallhs;
   SCIP_Real dualact;
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

   SCIPdebugMessage("undoing bound changes in LP exceeding cutoff: cutoff=%g, depth=%d\n",
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
   SCIP_CALL( SCIPsetAllocBufferArray(set, &primsols, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualsols, nrows) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &redcosts, ncols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &dualcoefs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &varredcosts, nvars) );

   /* get solution from LPI */
   SCIP_CALL( SCIPlpiGetSol(lpi, NULL, primsols, dualsols, NULL, redcosts) );
#ifdef SCIP_DEBUG
   {
      SCIP_Real objval;
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      SCIPdebugMessage(" -> LP objval: %g\n", objval);
   }
#endif

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

   BMSclearMemoryArray(dualcoefs, nvars);

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
      if( SCIPsetIsZero(set, dualsols[r]) )
         continue;

      /* check dual feasibility */
      if( (SCIPsetIsInfinity(set, -row->lhs) && dualsols[r] > 0.0) || (SCIPsetIsInfinity(set, row->rhs) && dualsols[r] < 0.0) )
      {
         SCIPdebugMessage(" -> infeasible dual solution %g in row <%s>: lhs=%g, rhs=%g\n",
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
         SCIPdebugMessage(" -> local row <%s>: dual=%g\n", SCIProwGetName(row), dualsols[r]);
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
         SCIPdebugMessage(" -> global row <%s>[%g,%g]: dual=%g -> duallhs=%g\n",
            SCIProwGetName(row), row->lhs - row->constant, row->rhs - row->constant, dualsols[r], duallhs);
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
         SCIP_COL* col;
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
            SCIPdebugMessage(" -> infeasible reduced costs %g in var <%s>: lb=%g, ub=%g\n",
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
   SCIPdebugMessage(" -> final dual values: lhs=%g, act=%g\n", duallhs, dualact);

   /* check, if the dual row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, duallhs, dualact) )
   {
      /* undo bound changes while keeping the infeasibility proof valid */
      SCIP_CALL( undoBdchgsProof(set, prob, currentdepth, dualcoefs, duallhs, dualact,
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
SCIP_RETCODE conflictAnalyzeRemainingBdchgs(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   int*                  lbchginfoposs,      /**< positions of currently active lower bound change infos in variables' arrays */
   int*                  ubchginfoposs,      /**< positions of currently active upper bound change infos in variables' arrays */
   int*                  nclauses,           /**< pointer to store the number of generated conflict clauses */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict clauses */
   int*                  nreconvclauses,     /**< pointer to store the number of generated reconvergence clauses */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence clauses */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   int nvars;
   int v;
   int nbdchgs;
   int maxsize;

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

   maxsize = 2*conflictCalcMaxsize(set, prob);

   /* initialize conflict data */
   SCIP_CALL( SCIPconflictInit(conflict) );

   /* add remaining bound changes to conflict queue */
   SCIPdebugMessage("initial conflict set after undoing bound changes:\n");
   nbdchgs = 0;
   for( v = 0; v < nvars && nbdchgs < maxsize; ++v )
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
         SCIPdebugMessage("   fixed: <%s> == %g [status: %d, type: %d, dive/strong]\n",
            SCIPvarGetName(var), SCIPvarGetLbLP(var), SCIPvarGetStatus(var), SCIPvarGetType(var));
         SCIP_CALL( conflictAddConflictVar(conflict, blkmem, set, stat, var, FALSE, FALSE) );
         nbdchgs++;
      }
      else
      {
         /* put remaining bound changes into conflict candidate queue */
         if( lbchginfoposs[v] >= 0 )
         {
            SCIPdebugMessage("   queue: <%s> >= %g [status: %d, type: %d, depth: %d, pos:%d, chgtype: %d]\n",
               SCIPvarGetName(var), SCIPbdchginfoGetNewbound(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPvarGetStatus(var), SCIPvarGetType(var),
               SCIPbdchginfoGetDepth(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPbdchginfoGetPos(&var->lbchginfos[lbchginfoposs[v]]),
               SCIPbdchginfoGetChgtype(&var->lbchginfos[lbchginfoposs[v]]));
            SCIP_CALL( conflictAddBdchginfo(conflict, &var->lbchginfos[lbchginfoposs[v]]) );
            nbdchgs++;
         }
         if( ubchginfoposs[v] >= 0 )
         {
            SCIPdebugMessage("   queue: <%s> <= %g [status: %d, type: %d, depth: %d, pos:%d, chgtype: %d]\n",
               SCIPvarGetName(var), SCIPbdchginfoGetNewbound(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPvarGetStatus(var), SCIPvarGetType(var),
               SCIPbdchginfoGetDepth(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPbdchginfoGetPos(&var->ubchginfos[ubchginfoposs[v]]),
               SCIPbdchginfoGetChgtype(&var->ubchginfos[ubchginfoposs[v]]));
            SCIP_CALL( conflictAddBdchginfo(conflict, &var->ubchginfos[ubchginfoposs[v]]) );
            nbdchgs++;
         }
      }
   }

   if( v == nvars )
   {
      /* analyze the conflict set, and create conflict constraints on success */
      SCIP_CALL( conflictAnalyze(conflict, blkmem, set, stat, prob, tree, 0, FALSE,
            nclauses, nliterals, nreconvclauses, nreconvliterals) );
   }

   return SCIP_OKAY;
}

/** actually performs analyzation of infeasible LP */
static
SCIP_RETCODE conflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool             diving,             /**< are we in strong branching or diving mode? */
   int*                  iterations,         /**< pointer to store the total number of LP iterations used */
   int*                  nclauses,           /**< pointer to store the number of generated conflict clauses */
   int*                  nliterals,          /**< pointer to store the number of literals in generated conflict clauses */
   int*                  nreconvclauses,     /**< pointer to store the number of generated reconvergence clauses */
   int*                  nreconvliterals     /**< pointer to store the number of literals generated reconvergence clauses */
   )
{
   SCIP_LPI* lpi;
   SCIP_VAR** vars;
   SCIP_Real* curvarlbs;
   SCIP_Real* curvarubs;
   int* lbchginfoposs;
   int* ubchginfoposs;
   SCIP_Real objval;
   int nvars;
   int v;
   int currentdepth;
   SCIP_Bool valid;

   assert(conflict != NULL);
   assert(conflict->nclauses == 0);
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

   /* get LP solver interface */
   lpi = SCIPlpGetLPI(lp);
   assert(SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsOptimal(lpi));
   assert(SCIPlpiIsPrimalInfeasible(lpi) || !SCIPlpDivingObjChanged(lp));

   if( SCIPlpiIsObjlimExc(lpi) )
   {
      assert(!SCIPlpDivingObjChanged(lp));

      /* make sure, a dual feasible solution exists, that exceeds the objective limit;
       * With FASTMIP setting, CPLEX does not apply the final pivot to reach the dual solution exceeding the objective
       * limit. Therefore, we have turn off FASTMIP and resolve the problem.
       */
      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiuobjlim )
      {
         SCIP_RETCODE retcode;
         int oldfastmip;
         int iter;

         /* get current value for FASTMIP */
         retcode = SCIPlpiGetIntpar(lpi, SCIP_LPPAR_FASTMIP, &oldfastmip);
	 if( retcode != SCIP_PARAMETERUNKNOWN )
	 {
	    SCIP_CALL( retcode );
	    if( !oldfastmip )
	       return SCIP_OKAY;
	 }
         else
            return SCIP_OKAY;

         /* temporarily disable FASTMIP in LP solver */
         SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FASTMIP, FALSE) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve LP */
         retcode = SCIPlpiSolveDual(lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         valid = (retcode != SCIP_LPERROR);
         if( valid )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lpi, &iter) );
            (*iterations) += iter;
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            SCIPdebugMessage(" -> resolved objlim exceeding LP in %d iterations (total: %lld) (infeasible:%d, objlim: %d, optimal:%d)\n",
               iter, stat->nconflictlpiterations, SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi),
               SCIPlpiIsOptimal(lpi));
            valid = (SCIPlpiIsObjlimExc(lpi) || SCIPlpiIsPrimalInfeasible(lpi) || SCIPlpiIsOptimal(lpi));
         }

         /* reinstall FASTMIP in LP solver */
         SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FASTMIP, oldfastmip) );

         /* abort, if the LP produced an error */
         if( !valid )
            return SCIP_OKAY;
      }
   }

   if( SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi) )
   {
      assert(!SCIPlpDivingObjChanged(lp));

      SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
      if( objval < lp->lpiuobjlim )
      {
         SCIPdebugMessage(" -> LP does not exceed the cutoff bound: obj=%g, cutoff=%g\n", objval, lp->lpiuobjlim);
         return SCIP_OKAY;
      }
      else
      {
         SCIPdebugMessage(" -> LP exceeds the cutoff bound: obj=%g, cutoff=%g\n", objval, lp->lpiuobjlim);
      }
   }

   currentdepth = SCIPtreeGetCurrentDepth(tree);

   SCIPdebugMessage("analyzing conflict on infeasible LP (infeasible: %d, objlimexc: %d, optimal:%d) in depth %d\n",
      SCIPlpiIsPrimalInfeasible(lpi), SCIPlpiIsObjlimExc(lpi), SCIPlpiIsOptimal(lpi), currentdepth);
#ifdef SCIP_DEBUG
   {
      SCIP_Real uobjlim;

      SCIP_CALL( SCIPlpiGetRealpar(lpi, SCIP_LPPAR_UOBJLIM, &uobjlim) );
      SCIPdebugMessage(" -> objective limit in LP solver: %g (in LP: %g)\n", uobjlim, lp->lpiuobjlim);
   }
#endif

   /* get active problem variables */
   vars = prob->vars;
   nvars = prob->nvars;

   /* get temporary memory for remembering variables' current bounds and corresponding bound change information
    * positions in variable's bound change information arrays
    */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbchginfoposs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubchginfoposs, nvars) );

   /* the following algorithm is used to find a subset of changed bounds leading to an infeasible LP:
    * 1. undo all strong branching or diving bound changes on non-binary variables
    *    -> update lb/ubchginfoposs arrays
    *    -> store changes in bdchg and curvarlbs/ubs arrays
    *    -> apply changes to the SCIP_LPI
    * 2. undo bound changes on non-binary variables in current depth without inference information
    *    -> update lb/ubchginfoposs arrays
    *    -> store additional changes in bdchg and curvarlbs/ubs arrays
    *    -> apply additional changes to the SCIP_LPI
    * 3. resolve SCIP_LP
    * 4. call undoBdchgsDualfarkas() or undoBdchgsDualsol()
    *    -> update lb/ubchginfoposs arrays
    *    -> store additional changes in bdchg and curvarlbs/ubs arrays
    *    -> apply additional changes to the SCIP_LPI
    * 5. if additional bound changes were undone, goto 2
    * 6. redo all bound changes in the LPI to restore the LPI to its original state
    * 7. analyze conflict
    *    -> put remaining changed bounds (see lb/ubchginfoposs arrays) into starting conflict set
    */

   /* get current bounds and current positions in lb/ubchginfos arrays of binary variables */
   valid = FALSE;
   for( v = 0; v < prob->nbinvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real lb;
      SCIP_Real ub;

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
      int* bdchginds;
      SCIP_Real* bdchgoldlbs;
      SCIP_Real* bdchgoldubs;
      SCIP_Real* bdchgnewlbs;
      SCIP_Real* bdchgnewubs;
      SCIP_Bool resolve;
      SCIP_Bool solvelp;
      int bdchgssize;
      int nbdchgs;
      int lastnbdchgs;
      int ncols;

      /* get number of columns in the LP */
      ncols = SCIPlpGetNCols(lp);

      /* get temporary memory for remembering bound changes on LPI columns */
      SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchginds, ncols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchgoldlbs, ncols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchgoldubs, ncols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchgnewlbs, ncols) );
      SCIP_CALL( SCIPsetAllocBufferArray(set, &bdchgnewubs, ncols) );
      bdchgssize = ncols;
      nbdchgs = 0;
      lastnbdchgs = 0;

      /* get current bounds and current positions in lb/ubchginfos arrays of non-binary variables */
      for( v = prob->nbinvars; v < nvars; ++v )
      {
         SCIP_VAR* var;

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
            SCIP_Real lb;
            SCIP_Real ub;
            SCIP_Bool changed;

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
               SCIP_CALL( addBdchg(set, var, curvarlbs[v], curvarubs[v],
                     &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs) );
               SCIPdebugMessage(" -> undoing strong branching/diving bound change on non-binary variable <%s>: [%g,%g], bdchg=(%d,%d)\n",
                  SCIPvarGetName(var), curvarlbs[v], curvarubs[v], lbchginfoposs[v], ubchginfoposs[v]);
            }
         }
      }

      /* undo as many bound changes as possible with the current LP solution */
      if( SCIPlpiIsPrimalInfeasible(lpi) )
      {
         SCIP_CALL( undoBdchgsDualfarkas(set, prob, lp, currentdepth, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
               &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
               &valid, &resolve) );
      }
      else
      {
         assert(SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi));
         SCIP_CALL( undoBdchgsDualsol(set, prob, lp, currentdepth, curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs,
               &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
               &valid, &resolve) );
      }

      /* check if we want to solve the LP and if all columns are present in the LP */
      solvelp = (set->conf_maxlploops > 0 && SCIPprobAllColsInLP(prob, set, lp));

      if( valid && resolve && solvelp )
      {
	 SCIP_RETCODE retcode;
         SCIP_ROW** rows;
         int* sidechginds;
         SCIP_Real* sidechgoldlhss;
         SCIP_Real* sidechgoldrhss;
         SCIP_Real* sidechgnewlhss;
         SCIP_Real* sidechgnewrhss;
         SCIP_Real lpiinfinity;
         int sidechgssize;
         int nsidechgs;
         int nrows;
         int nloops;
         int oldfastmip;
         int r;

         /* get infinity value of LP solver */
         lpiinfinity = SCIPlpiInfinity(lpi);

         /* disable FASTMIP setting */
         retcode = SCIPlpiGetIntpar(lpi, SCIP_LPPAR_FASTMIP, &oldfastmip);
	 if( retcode != SCIP_PARAMETERUNKNOWN )
	 {
	    SCIP_CALL( retcode );
	    SCIP_CALL( SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FASTMIP, FALSE) );
	 }

         /* temporarily remove objective limit in LP solver */
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, lpiinfinity) );

         /* get LP rows */
         rows = SCIPlpGetRows(lp);
         nrows = SCIPlpGetNRows(lp);
         assert(nrows == 0 || rows != NULL);

         /* get temporary memory for remembering side changes on LPI rows */
         SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechginds, nrows) );
         SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgoldlhss, nrows) );
         SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgoldrhss, nrows) );
         SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgnewlhss, nrows) );
         SCIP_CALL( SCIPsetAllocBufferArray(set, &sidechgnewrhss, nrows) );
         sidechgssize = nrows;
         nsidechgs = 0;

         /* remove all local rows by setting their sides to infinity;
          * finite sides are only changed to near infinity, such that the row's sense in the LP solver
          * is not affected (e.g. CPLEX cannot handle free rows)
          */
         for( r = 0 ; r < nrows; ++r )
         {
            assert(SCIProwGetLPPos(rows[r]) == r);

            if( SCIProwIsLocal(rows[r]) )
            {
               SCIPdebugMessage(" -> removing local row <%s> [%g,%g]\n",
                  SCIProwGetName(rows[r]), SCIProwGetLhs(rows[r]), SCIProwGetRhs(rows[r]));
               SCIP_CALL( addSideRemoval(set, rows[r], lpiinfinity, &sidechginds, &sidechgoldlhss, &sidechgoldrhss,
                     &sidechgnewlhss, &sidechgnewrhss, &sidechgssize, &nsidechgs) );
            }
         }

         /* apply changes of local rows to the LP solver */
         if( nsidechgs > 0 )
         {
            SCIP_CALL( SCIPlpiChgSides(lpi, nsidechgs, sidechginds, sidechgnewlhss, sidechgnewrhss) );
         }

         /* undo as many additional bound changes as possible by resolving the LP */
         assert(valid);
         assert(resolve);
         nloops = 0;
         while( valid && resolve && nloops < set->conf_maxlploops )
         {
            int iter;

            nloops++;
            resolve = FALSE;

            SCIPdebugMessage("infeasible LP conflict analysis loop %d (changed col bounds: %d)\n", nloops, nbdchgs);

            /* apply bound changes to the LP solver */
            assert(nbdchgs >= lastnbdchgs);
            if( nbdchgs > lastnbdchgs )
            {
               SCIPdebugMessage(" -> applying %d bound changes to the LP solver (total: %d)\n", nbdchgs - lastnbdchgs, nbdchgs);
               SCIP_CALL( SCIPlpiChgBounds(lpi, nbdchgs - lastnbdchgs,
                     &bdchginds[lastnbdchgs], &bdchgnewlbs[lastnbdchgs], &bdchgnewubs[lastnbdchgs]) );
               lastnbdchgs = nbdchgs;
            }

            /* start LP timer */
            SCIPclockStart(stat->conflictlptime, set);

            /* resolve LP */
            retcode = SCIPlpiSolveDual(lpi);

            /* stop LP timer */
            SCIPclockStop(stat->conflictlptime, set);

            /* check return code of LP solving call */
            if( retcode == SCIP_LPERROR )
            {
               valid = FALSE;
               break;
            }
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lpi, &iter) );
            (*iterations) += iter;
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            SCIPdebugMessage(" -> resolved LP in %d iterations (total: %lld) (infeasible:%d)\n",
               iter, stat->nconflictlpiterations, SCIPlpiIsPrimalInfeasible(lpi));

            /* evaluate result */
            if( SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi) )
            {
               SCIP_CALL( SCIPlpiGetObjval(lpi, &objval) );
               valid = (objval >= lp->lpiuobjlim && !SCIPlpDivingObjChanged(lp));
            }
            else
               valid = SCIPlpiIsPrimalInfeasible(lpi);

            if( valid )
            {
               /* undo additional bound changes */
               if( SCIPlpiIsPrimalInfeasible(lpi) )
               {
                  SCIP_CALL( undoBdchgsDualfarkas(set, prob, lp, currentdepth, curvarlbs, curvarubs,
                        lbchginfoposs, ubchginfoposs,
                        &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
                        &valid, &resolve) );
               }
               else
               {
                  assert(SCIPlpiIsOptimal(lpi) || SCIPlpiIsObjlimExc(lpi));
                  SCIP_CALL( undoBdchgsDualsol(set, prob, lp, currentdepth, curvarlbs, curvarubs,
                        lbchginfoposs, ubchginfoposs,
                        &bdchginds, &bdchgoldlbs, &bdchgoldubs, &bdchgnewlbs, &bdchgnewubs, &bdchgssize, &nbdchgs,
                        &valid, &resolve) );
               }
            }
            assert(!resolve || valid);
            assert(!resolve || nbdchgs > lastnbdchgs);
            SCIPdebugMessage(" -> finished infeasible LP conflict analysis loop %d (iter: %d, nbdchgs: %d)\n",
               nloops, iter, nbdchgs - lastnbdchgs);
         }

         SCIPdebugMessage("finished undoing bound changes after %d loops (valid=%d, nbdchgs: %d)\n",
            nloops, valid, nbdchgs);

         /* reset variables to local bounds */
         if( nbdchgs > 0 )
         {
            SCIP_CALL( SCIPlpiChgBounds(lpi, nbdchgs, bdchginds, bdchgoldlbs, bdchgoldubs) );
         }

         /* reset changes of local rows */
         if( nsidechgs > 0 )
         {
            SCIP_CALL( SCIPlpiChgSides(lpi, nsidechgs, sidechginds, sidechgoldlhss, sidechgoldrhss) );
         }

         /* reset objective limit in LP solver */
         SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, lp->lpiuobjlim) );

         /* reinstall FASTMIP in LP solver */
         retcode = SCIPlpiSetIntpar(lpi, SCIP_LPPAR_FASTMIP, oldfastmip);
	 if( retcode != SCIP_PARAMETERUNKNOWN )
	 {
	    SCIP_CALL( retcode );
	 }

         /* free temporary memory */
         SCIPsetFreeBufferArray(set, &sidechgnewrhss);
         SCIPsetFreeBufferArray(set, &sidechgnewlhss);
         SCIPsetFreeBufferArray(set, &sidechgoldrhss);
         SCIPsetFreeBufferArray(set, &sidechgoldlhss);
         SCIPsetFreeBufferArray(set, &sidechginds);
      }

      /* analyze the conflict starting with remaining bound changes */
      if( valid )
      {
         SCIP_CALL( conflictAnalyzeRemainingBdchgs(conflict, blkmem, set, stat, prob, tree,
               lbchginfoposs, ubchginfoposs, nclauses, nliterals, nreconvclauses, nreconvliterals) );
      }

      /* free temporary memory */
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

   /* flush conflict clause storage */
   SCIP_CALL( SCIPconflictFlushClauses(conflict, blkmem, set, stat, prob, tree) );

   return SCIP_OKAY;
}

/** analyzes an infeasible LP to find out the bound changes on binary variables that were responsible for the
 *  infeasibility;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible LP conflict analysis
 */
SCIP_RETCODE SCIPconflictAnalyzeLP(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
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

   SCIPdebugMessage("analyzing conflict on infeasible LP in depth %d (solstat: %d, objchanged: %d)\n",
      SCIPtreeGetCurrentDepth(tree), SCIPlpGetSolstat(lp), SCIPlpDivingObjChanged(lp));

   /* start timing */
   SCIPclockStart(conflict->lpanalyzetime, set);
   conflict->nlpcalls++;

   /* perform conflict analysis */
   SCIP_CALL( conflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, SCIPlpDiving(lp),
         &iterations, &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
   conflict->nlpsuccess += (nclauses > 0);
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
SCIP_Real SCIPconflictGetLPTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->lpanalyzetime);
}

/** gets number of calls to infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNLPCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpcalls;
}

/** gets number of calls to infeasible LP conflict analysis that yield at least one conflict clause */
SCIP_Longint SCIPconflictGetNLPSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpsuccess;
}

/** gets number of conflict clauses detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNLPConflictClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpconfclauses;
}

/** gets total number of literals in conflict clauses created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNLPConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpconfliterals;
}

/** gets number of reconvergence clauses detected in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNLPReconvergenceClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpreconvclauses;
}

/** gets total number of literals in reconvergence clauses created in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNLPReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpreconvliterals;
}

/** gets number of LP iterations in infeasible LP conflict analysis */
SCIP_Longint SCIPconflictGetNLPIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpiterations;
}




/*
 * infeasible strong branching conflict analysis
 */

/** analyses infeasible strong branching sub problems for conflicts */
SCIP_RETCODE SCIPconflictAnalyzeStrongbranch(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory buffers */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< dynamic problem statistics */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_COL*             col,                /**< LP column with at least one infeasible strong branching subproblem */
   SCIP_Bool*            downconflict,       /**< pointer to store whether a conflict clause was created for an infeasible
                                              *   downwards branch, or NULL */
   SCIP_Bool*            upconflict          /**< pointer to store whether a conflict clause was created for an infeasible
                                              *   upwards branch, or NULL */
   )
{
   SCIP_VAR* var;
   int* cstat;
   int* rstat;
   SCIP_Real oldlb;
   SCIP_Real oldub;
   SCIP_Real newlb;
   SCIP_Real newub;
   int iter;
   int nclauses;
   int nliterals;
   int nreconvclauses;
   int nreconvliterals;

   assert(stat != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(SCIPprobAllColsInLP(prob, set, lp));
   assert(col != NULL);
   assert((col->sbdownvalid && SCIPsetIsGE(set, col->sbdown, lp->cutoffbound)
         && SCIPsetFeasCeil(set, col->primsol-1.0) >= col->lb - 0.5)
      || (col->sbupvalid && SCIPsetIsGE(set, col->sbup, lp->cutoffbound)
         && SCIPsetFeasFloor(set, col->primsol+1.0) <= col->ub + 0.5));

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

   /* get temporary memory for storing current LP basis */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &cstat, lp->nlpicols) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &rstat, lp->nlpirows) );

   /* get current LP basis */
   SCIP_CALL( SCIPlpiGetBase(lp->lpi, cstat, rstat) );

   /* remember old bounds */
   oldlb = col->lb;
   oldub = col->ub;

   var = SCIPcolGetVar(col);

   /* is down branch infeasible? */
   if( col->sbdownvalid && SCIPsetIsGE(set, col->sbdown, lp->cutoffbound) )
   {
      newub = SCIPsetFeasCeil(set, col->primsol-1.0);
      if( newub >= col->lb - 0.5 )
      {
         SCIP_RETCODE retcode;

         SCIPdebugMessage("analyzing conflict on infeasible downwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPtreeGetCurrentDepth(tree));

         conflict->nsbcalls++;

         /* change the upper bound */
         col->ub = newub;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         retcode = SCIPlpiSolveDual(lp->lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode != SCIP_LPERROR )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lp->lpi, &iter) );
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            conflict->nsbiterations += iter;
            SCIPdebugMessage(" -> resolved downwards strong branching LP in %d iterations\n", iter);

            /* perform conflict analysis on infeasible LP */
            SCIP_CALL( conflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, TRUE,
                  &iter, &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
            conflict->nsbsuccess += (nclauses > 0);
            conflict->nsbiterations += iter;
            conflict->nsbconfclauses += nclauses;
            conflict->nsbconfliterals += nliterals;
            conflict->nsbreconvclauses += nreconvclauses;
            conflict->nsbreconvliterals += nreconvliterals;
            if( downconflict != NULL )
               *downconflict = (nclauses > 0);
         }

         /* reset the upper bound */
         col->ub = oldub;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* reset LP basis */
         SCIP_CALL( SCIPlpiSetBase(lp->lpi, cstat, rstat) );
      }
   }

   /* is up branch infeasible? */
   if( col->sbupvalid && SCIPsetIsGE(set, col->sbup, lp->cutoffbound) )
   {
      newlb = SCIPsetFeasFloor(set, col->primsol+1.0);
      if( newlb <= col->ub + 0.5 )
      {
         SCIP_RETCODE retcode;

         SCIPdebugMessage("analyzing conflict on infeasible upwards strongbranch for variable <%s>[%g,%g] in depth %d\n",
            SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPtreeGetCurrentDepth(tree));

         conflict->nsbcalls++;

         /* change the lower bound */
         col->lb = newlb;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* start LP timer */
         SCIPclockStart(stat->conflictlptime, set);

         /* resolve the LP */
         retcode = SCIPlpiSolveDual(lp->lpi);

         /* stop LP timer */
         SCIPclockStop(stat->conflictlptime, set);

         /* check return code of LP solving call */
         if( retcode != SCIP_LPERROR )
         {
            SCIP_CALL( retcode );

            /* count number of LP iterations */
            SCIP_CALL( SCIPlpiGetIterations(lp->lpi, &iter) );
            stat->nconflictlps++;
            stat->nconflictlpiterations += iter;
            conflict->nsbiterations += iter;
            SCIPdebugMessage(" -> resolved upwards strong branching LP in %d iterations\n", iter);

            /* perform conflict analysis on infeasible LP */
            SCIP_CALL( conflictAnalyzeLP(conflict, blkmem, set, stat, prob, tree, lp, TRUE,
                  &iter, &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
            conflict->nsbsuccess += (nclauses > 0);
            conflict->nsbiterations += iter;
            conflict->nsbconfclauses += nclauses;
            conflict->nsbconfliterals += nliterals;
            conflict->nsbreconvclauses += nreconvclauses;
            conflict->nsbreconvliterals += nreconvliterals;
            if( upconflict != NULL )
               *upconflict = (nclauses > 0);
         }

         /* reset the lower bound */
         col->lb = oldlb;
         SCIP_CALL( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );

         /* reset LP basis */
         SCIP_CALL( SCIPlpiSetBase(lp->lpi, cstat, rstat) );
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
SCIP_Real SCIPconflictGetStrongbranchTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->sbanalyzetime);
}

/** gets number of calls to infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbcalls;
}

/** gets number of calls to infeasible strong branching conflict analysis that yield at least one conflict clause */
SCIP_Longint SCIPconflictGetNStrongbranchSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbsuccess;
}

/** gets number of conflict clauses detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfclauses;
}

/** gets total number of literals in conflict clauses created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconfliterals;
}

/** gets number of reconvergence clauses detected in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvclauses;
}

/** gets total number of literals in reconvergence clauses created in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbreconvliterals;
}

/** gets number of LP iterations in infeasible strong branching conflict analysis */
SCIP_Longint SCIPconflictGetNStrongbranchIterations(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
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
SCIP_RETCODE SCIPconflictAnalyzePseudo(
   SCIP_CONFLICT*        conflict,           /**< conflict analysis data */
   BMS_BLKMEM*           blkmem,             /**< block memory of transformed problem */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics */
   SCIP_PROB*            prob,               /**< problem data */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_LP*              lp,                 /**< LP data */
   SCIP_Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real* curvarlbs;
   SCIP_Real* curvarubs;
   int* lbchginfoposs;
   int* ubchginfoposs;
   SCIP_Real* pseudocoefs;
   SCIP_Real pseudolhs;
   SCIP_Real pseudoact;
   int nvars;
   int v;

   assert(conflict != NULL);
   assert(conflict->nclauses == 0);
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

   SCIPdebugMessage("analyzing pseudo solution (obj: %g) that exceeds objective limit (%g)\n",
      SCIPlpGetPseudoObjval(lp, set), lp->cutoffbound);

   /* start timing */
   SCIPclockStart(conflict->pseudoanalyzetime, set);
   conflict->npseudocalls++;

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
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &lbchginfoposs, nvars) );
   SCIP_CALL( SCIPsetAllocBufferArray(set, &ubchginfoposs, nvars) );

   /* get temporary memory for infeasibility proof coefficients */
   SCIP_CALL( SCIPsetAllocBufferArray(set, &pseudocoefs, nvars) );

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
   SCIPdebugMessage("  -> recalculated pseudo infeasibility proof:  %g <= %g\n", pseudolhs, pseudoact);

   /* check, if the pseudo row is still violated (after recalculation of pseudo activity) */
   if( SCIPsetIsFeasGT(set, pseudolhs, pseudoact) )
   {
      int nclauses;
      int nliterals;
      int nreconvclauses;
      int nreconvliterals;

      /* undo bound changes without destroying the infeasibility proof */
      SCIP_CALL( undoBdchgsProof(set, prob, SCIPtreeGetCurrentDepth(tree), pseudocoefs, pseudolhs, pseudoact,
            curvarlbs, curvarubs, lbchginfoposs, ubchginfoposs, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL) );

      /* analyze conflict on remaining bound changes */
      SCIP_CALL( conflictAnalyzeRemainingBdchgs(conflict, blkmem, set, stat, prob, tree,
            lbchginfoposs, ubchginfoposs, &nclauses, &nliterals, &nreconvclauses, &nreconvliterals) );
      conflict->npseudosuccess += (nclauses > 0);
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

   /* flush conflict clause storage */
   SCIP_CALL( SCIPconflictFlushClauses(conflict, blkmem, set, stat, prob, tree) );

   /* stop timing */
   SCIPclockStop(conflict->pseudoanalyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing pseudo solution conflicts */
SCIP_Real SCIPconflictGetPseudoTime(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->pseudoanalyzetime);
}

/** gets number of calls to pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoCalls(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudocalls;
}

/** gets number of calls to pseudo solution conflict analysis that yield at least one conflict clause */
SCIP_Longint SCIPconflictGetNPseudoSuccess(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudosuccess;
}

/** gets number of conflict clauses detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfclauses;
}

/** gets total number of literals in conflict clauses created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoConflictLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconfliterals;
}

/** gets number of reconvergence clauses detected in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceClauses(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvclauses;
}

/** gets total number of literals in reconvergence clauses created in pseudo solution conflict analysis */
SCIP_Longint SCIPconflictGetNPseudoReconvergenceLiterals(
   SCIP_CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoreconvliterals;
}
