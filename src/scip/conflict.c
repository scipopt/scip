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
#pragma ident "@(#) $Id: conflict.c,v 1.47 2004/08/06 08:18:01 bzfpfend Exp $"

/**@file   conflict.c
 * @brief  methods and datastructures for conflict analysis
 * @author Tobias Achterberg
 *
 * This file implements a conflict analysis method like the one used in modern
 * SAT solvers like zchaff. The algorithm works as follows:
 * 
 * Given is a set of binary variables that are not allowed being set to their
 * current values simultaneously (e.g. a single clause, i.e. a logicor constraint),
 * rendering the current node infeasible.
 * The goal is to deduce a different clause -- a conflict clause -- representing
 * the "reason" for this conflict, i.e. the branching decisions or the deductions
 * (applied e.g. in domain propagation) that lead to the conflict. This clause
 * can then be added to the constraint set to help cutting off similar parts
 * of the branch and bound tree, that would lead to the same conflict.
 *
 *  1. Put all the given variables to a priority queue, which is ordered,
 *     such that the variable that was fixed last due to branching or deduction
 *     is at the top of the queue. The variables in the queue are always active
 *     problem variables.
 *     Create an empty conflict set.
 *  2. Remove the top variable v from the priority queue.
 *  3. (a) If the remaining queue is non-empty, and variable w (the one that is now
 *         on the top of the queue) was fixed at the same depth level as v, and if
 *         the assignment to v was a deduction with known inference reason, and if
 *         the inference constraint is globally valid:
 *          - Resolve variable v by asking the constraint that infered the
 *            assignment of v to put all the variables that when assigned to their
 *            current values lead to the deduction of v's current value on the
 *            priority queue.
 *            Note that these variables have at most the same inference depth
 *            level as variable v, and were deduced earlier than v.
 *     (b) Otherwise, the assignment to variable v was a branching decision or
 *         a deduction with missing inference reason.
 *          - Put v or the negation of v into the conflict set, depending on which
 *            of the ones is currently fixed to FALSE.
 *            Note that if v was a branching variable, all deduced variables remaining
 *            in the priority queue have smaller inference depth level than v, since
 *            deductions are always applied after the branching decisions. However,
 *            there is the possibility, that v was a deduced variable, where the
 *            inference reason was not given. With this lack of information, we must
 *            treat the deduced variable like a branching variable, and there may
 *            exist other deduced variables of the same inference depth level in
 *            the priority queue.
 *  4. If priority queue is non-empty, goto step 2.
 *  5. The conflict set represents the conflict clause saying that at least one
 *     of the conflict variables must be set to TRUE.
 *     The conflict set is then passed to the conflict handlers, that may create 
 *     a corresponding constraint (e.g. a logicor constraint) out of these conflict
 *     variables and add it to the problem.
 *
 * If all deduced variables come with inference information, the resulting
 * conflict set has the property, that for each depth level at most one variable
 * assigned at that level is member of the conflict set. This conflict variable
 * is the first unique implication point of its depth level (FUIP).
 *
 * The user has to do the following to get the conflict analysis running in its
 * current implementation:
 *  - A constraint handler supporting the conflict analysis must implement
 *    the CONSRESCVAR call, that processes a candidate variable v and puts all
 *    the reason variables leading to the assignment of v with a call to
 *    SCIPaddConflictVar() on the conflict queue (algorithm step 3.(a)).
 *  - If the current assignment of the binary variables lead to a deduction of
 *    a different binary variable, the constraint handler should call
 *    SCIPinferBinVar(), thus providing the constraint that infered the
 *    assignment to SCIP.
 *  - If an (with the current assignment) infeasible globally valid constraint
 *    is detected, the constraint handler should
 *     1. call SCIPinitConflictAnalysis() to initialize the conflict queue,
 *     2. call SCIPaddConflictVar() for each variable in the infeasible
 *        constraint,
 *     3. call SCIPanalyzeConflict() to analyze the conflict and add an
 *        appropriate conflict constraint.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "stat.h"
#include "clock.h"
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

#include "struct_conflict.h"



/*
 * dynamic memory arrays
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
      conflict->conflictvarssize = newsize;
   }
   assert(num <= conflict->conflictvarssize);

   return SCIP_OKAY;
}




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
   assert(SCIPstage(scip) == SCIP_STAGE_INIT);

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
      errorMessage("Conflict handler <%s> already initialized\n", conflicthdlr->name);
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
      errorMessage("Conflict handler <%s> not initialized\n", conflicthdlr->name);
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
      CHECK_OKAY( conflicthdlr->conflictexec(scip, conflicthdlr, node, conflictvars, nconflictvars, resolved, result) );

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

/** compares two binary variables w.r.t. their inference information, such that variables infered later are
 *  ordered prior to variables infered earlier
 */
static
DECL_SORTPTRCOMP(conflictVarComp)
{  /*lint --e{715}*/
   VAR* var1;
   VAR* var2;
   
   var1 = (VAR*)elem1;
   var2 = (VAR*)elem2;
   assert(var1 != NULL);
   assert(SCIPvarGetType(var1) == SCIP_VARTYPE_BINARY);
   assert(var2 != NULL);
   assert(SCIPvarGetType(var2) == SCIP_VARTYPE_BINARY);

   if( SCIPvarGetFixDepth(var1) > SCIPvarGetFixDepth(var2) )
      return -1;
   else if( SCIPvarGetFixDepth(var1) < SCIPvarGetFixDepth(var2) )
      return +1;
   else if( SCIPvarGetFixIndex(var1) > SCIPvarGetFixIndex(var2) )
      return -1;
   else if( SCIPvarGetFixIndex(var1) < SCIPvarGetFixIndex(var2) )
      return +1;
   else if( SCIPvarGetIndex(var1) > SCIPvarGetIndex(var2) )
      return +1;
   else if( SCIPvarGetIndex(var1) < SCIPvarGetIndex(var2) )
      return -1;
   else
      return 0;
}

/* forward reference */
static
Real getUbLocal(
   VAR*             var                 /**< active or negated problem variable */
   );

/** gets current lower bound of variable; if the variable is a COLUMN, the column's bound is returned, because
 *  it can be temporarily different from the variable's bound due to diving or strong branching
 */
static
Real getLbLocal(
   VAR*             var                 /**< active or negated problem variable */
   )
{
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIPcolGetLb(SCIPvarGetCol(var));
   else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      return 1.0 - getUbLocal(SCIPvarGetNegatedVar(var));
   }
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);
      return SCIPvarGetLbLocal(var);
   }
}

/** gets current upper bound of variable; if the variable is a COLUMN, the column's bound is returned, because
 *  it can be temporarily different from the variable's bound due to diving or strong branching
 */
static
Real getUbLocal(
   VAR*             var                 /**< active or negated problem variable */
   )
{
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      return SCIPcolGetUb(SCIPvarGetCol(var));
   else if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED )
   {
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      return 1.0 - getLbLocal(SCIPvarGetNegatedVar(var));
   }
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED);
      return SCIPvarGetUbLocal(var);
   }
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
   CHECK_OKAY( SCIPpqueueCreate(&(*conflict)->varqueue, set->memgrowinit, set->memgrowfac, conflictVarComp) );
   (*conflict)->conflictvars = NULL;
   (*conflict)->conflictvarssize = 0;
   (*conflict)->nconflictvars = 0;
   (*conflict)->count = 0;
   (*conflict)->npropcalls = 0;
   (*conflict)->npropconflicts = 0;
   (*conflict)->nlpcalls = 0;
   (*conflict)->nlpconflicts = 0;
   (*conflict)->nlpiterations = 0;
   (*conflict)->nsbcalls = 0;
   (*conflict)->nsbconflicts = 0;
   (*conflict)->nsbiterations = 0;
   (*conflict)->npseudocalls = 0;
   (*conflict)->npseudoconflicts = 0;

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
   SCIPpqueueFree(&(*conflict)->varqueue);
   freeMemoryArrayNull(&(*conflict)->conflictvars);
   freeMemory(conflict);


   return SCIP_OKAY;
}

/** initializes the propagation conflict analysis by clearing the conflict variable candidate queue */
RETCODE SCIPconflictInit(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   debugMessage("initializing conflict analysis\n");

   SCIPpqueueClear(conflict->varqueue);
   conflict->nconflictvars = 0;

   return SCIP_OKAY;
}

/** adds variable to conflict variable candidates */
RETCODE SCIPconflictAddVar(
   CONFLICT*        conflict,           /**< conflict analysis data */
   VAR*             var                 /**< problem variable */
   )
{
   VAR* actvar;

   assert(conflict != NULL);
   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(getLbLocal(var) > 0.5 || getUbLocal(var) < 0.5);

   debugMessage(" -> adding variable <%s>[%g,%g] [status:%d, depth:%d, num:%d, cons:<%s>, info:%d] to conflict candidates\n",
      SCIPvarGetName(var), getLbLocal(var), getUbLocal(var), SCIPvarGetStatus(var),
      SCIPvarGetFixDepth(var), SCIPvarGetFixIndex(var), 
      SCIPvarGetInferCons(var) == NULL ? "null" : SCIPconsGetName(SCIPvarGetInferCons(var)), SCIPvarGetInferInfo(var));

   /* get active problem variable */
   actvar = SCIPvarGetProbvar(var);

   /* we can ignore fixed variables */
   if( actvar != NULL )
   {
      assert(SCIPvarGetType(actvar) == SCIP_VARTYPE_BINARY);
      assert(getLbLocal(actvar) > 0.5 || getUbLocal(actvar) < 0.5);
      
#ifdef DEBUG
      if( actvar != var )
      {
         debugMessage("   -> active conflict candidate is <%s>[%g,%g] [status:%d]\n",
            SCIPvarGetName(actvar), getLbLocal(actvar), getUbLocal(actvar), SCIPvarGetStatus(actvar));
      }
#endif

      /* put candidate in priority queue */
      CHECK_OKAY( SCIPpqueueInsert(conflict->varqueue, (void*)actvar) );
   }

   return SCIP_OKAY;
}

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set
 */
static
RETCODE conflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   int              maxsize,            /**< maximal size of conflict set */
   Bool             mustresolve,        /**< should the conflict set only be used, if a resolution was applied? */
   Bool*            success             /**< pointer to store whether the conflict set is valid */
   )
{
   VAR* var;
   VAR* nextvar;
   VAR* infervar;
   CONS* infercons;
   RESULT result;
   int nresolutions;

   assert(conflict != NULL);
   assert(tree != NULL);
   assert(success != NULL);

   debugMessage("analyzing conflict with %d conflict candidates (maxsize=%d)\n",
      SCIPpqueueNElems(conflict->varqueue), maxsize);

   *success = FALSE;

   /* check, if there is something to analyze */
   if( SCIPpqueueNElems(conflict->varqueue) == 0 )
   {
      errorMessage("no conflict variables to analyze\n");
      return SCIP_INVALIDDATA;
   }

   /* increase the conflict set counter, such that variables of new conflict set are labeled with this new counter */
   conflict->count++;

   /* clear the conflict set */
   conflict->nconflictvars = 0;
   nresolutions = 0;

   /* process all variables in the conflict candidate queue */
   var = (VAR*)(SCIPpqueueRemove(conflict->varqueue));
   while( var != NULL && conflict->nconflictvars < maxsize )
   {
#ifdef DEBUG
      {
         CONS* infercons = SCIPvarGetInferCons(var);
         int v;

         debugMessage("processing variable <%s>[%g,%g] [status:%d, depth:%d, num:%d, cons:<%s>(%s), info:%d]\n",
            SCIPvarGetName(var), getLbLocal(var), getUbLocal(var), SCIPvarGetStatus(var),
            SCIPvarGetFixDepth(var), SCIPvarGetFixIndex(var), 
            infercons == NULL ? "null" : SCIPconsGetName(infercons),
            infercons == NULL ? "--" : (SCIPconsIsGlobal(infercons) ? "global" : "local"),
            SCIPvarGetInferInfo(var));
         debugMessage(" - conflict set   :");
         for( v = 0; v < conflict->nconflictvars; ++v )
            printf(" <%s>", SCIPvarGetName(conflict->conflictvars[v]));
         printf("\n");
         debugMessage(" - candidate queue:");
         for( v = 0; v < SCIPpqueueNElems(conflict->varqueue); ++v )
            printf(" <%s>", SCIPvarGetName((VAR*)(SCIPpqueueElems(conflict->varqueue)[v])));
         printf("\n");
      }
#endif

      assert(SCIPsetIsEQ(set, getLbLocal(var), getUbLocal(var)));

      /* if the first variable on the remaining queue is equal to the current variable,
       * this is a multiple insertion in the conflict candidate queue and we can ignore the current
       * variable
       */
      nextvar = (VAR*)(SCIPpqueueFirst(conflict->varqueue));
      assert(nextvar == NULL || SCIPvarGetFixDepth(var) >= SCIPvarGetFixDepth(nextvar));

      if( var != nextvar )
      {
         Bool resolved = FALSE;

         /* check, if the variable can and should be resolved */
         infercons = SCIPvarGetInferCons(var);
         if( nextvar != NULL
            && infercons != NULL
            && SCIPconsIsGlobal(infercons)
            && SCIPvarGetFixDepth(var) == SCIPvarGetFixDepth(nextvar) )
         {
            /* resolve variable v by asking the constraint that infered the
             * assignment of v to put all the variables that when assigned to their
             * current values lead to the deduction of v on the priority queue.
             */
            infervar = SCIPvarGetInferVar(var);
            assert(infervar != NULL);
            assert(SCIPsetIsEQ(set, getLbLocal(infervar), getUbLocal(infervar)));
            debugMessage("resolving variable <%s>[%g,%g] [status:%d, depth:%d, num:%d, cons:<%s>(%s), info:%d]\n",
               SCIPvarGetName(var), getLbLocal(var), getUbLocal(var), SCIPvarGetStatus(var),
               SCIPvarGetFixDepth(var), SCIPvarGetFixIndex(var), 
               SCIPconsGetName(infercons), SCIPconsIsGlobal(infercons) ? "global" : "local",
               SCIPvarGetInferInfo(var));
            CHECK_OKAY( SCIPconsResolveConflictVar(infercons, set, SCIPvarGetInferVar(var), &result) );
            resolved = (result == SCIP_SUCCESS);
         }

         if( resolved )
            nresolutions++;
         else
         {
            assert(SCIPsetIsFeasEQ(set, getLbLocal(var), getUbLocal(var)));

            /* mark the active problem variable to belong to the current conflict */
            var->conflictsetcount = conflict->count;

            /* choose between variable or its negation, such that the literal is fixed to FALSE */
            if( SCIPsetIsEQ(set, getLbLocal(var), 1.0) )
            {
               /* variable is fixed to TRUE -> use the negation */
               CHECK_OKAY( SCIPvarNegate(var, memhdr, set, stat, &var) );
            }
            assert(SCIPsetIsFeasEQ(set, getLbLocal(var), getUbLocal(var)));
            assert(SCIPsetIsFeasZero(set, getLbLocal(var)));

            debugMessage("putting variable <%s> to conflict set\n", SCIPvarGetName(var));

            /* put variable in the conflict set */
            CHECK_OKAY( conflictEnsureConflictvarsMem(conflict, set, conflict->nconflictvars+1) );
            conflict->conflictvars[conflict->nconflictvars] = var;
            conflict->nconflictvars++;
         }
      }

      /* get next variable from the conflict candidate queue */
      var = (VAR*)(SCIPpqueueRemove(conflict->varqueue));
   }

   /* if a valid conflict set was found (where at least one resolution was applied), call the conflict handlers */
   if( var == NULL && (!mustresolve || nresolutions >= 1) )
   {
      NODE* node;
      BOUNDCHG* boundchgs;
      RESULT result;
      int nboundchgs;
      int nbranchingvars;
      int d;
      int b;
      int h;

      assert(conflict->conflictvars != NULL);
      assert(conflict->nconflictvars > 0);

      debugMessage(" -> found conflict set of size %d/%d after %d resolutions\n", 
         conflict->nconflictvars, maxsize, nresolutions);

      /* identify the depth, at which the conflict clause should be added:
       * - if the branching rule operates on variables only, and if all branching variables up to a certain
       *   depth level are member of the conflict, the conflict clause can only be violated in the subtree
       *   of the node at that depth, because in all other nodes, at least one of these branching variables
       *   takes a different value, such that the conflict clause is feasible
       * - if there is at least one branching variable in a node, we assume, that this branching was performed
       *   on variables, and that the siblings of this node are disjunct w.r.t. the branching variables' fixings
       * - the variables in the conflict set are labeled with the current conflict set counter
       * - there are no branching variables in the root node, so we can start with depth level 1
       */
      for( d = 1; d < tree->pathlen; ++d )
      {
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
            if( SCIPvarGetBoundchgType(var) != SCIP_BOUNDCHGTYPE_BRANCHING )
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

      /* now, d is the depth level of the first node, that has non-variable branching or a branching variable
       * not in the conflict set; this means, the siblings of the node may benefit from the conflict clause,
       * and the clause should be added to the node's parent, i.e. at depth level d-1
       */
      assert(1 <= d && d <= tree->pathlen);
      debugMessage(" -> conflict found at depth %d is active in depth %d\n", tree->pathlen-1, d-1);

      /* if all branching variables are in the conflict set, the conflict clause is of no use */
      if( d == tree->pathlen )
         return SCIP_OKAY;

      /* sort conflict handlers by priority */
      SCIPsetSortConflicthdlrs(set);

      /* call conflict handlers to create a conflict constraint */
      for( h = 0; h < set->nconflicthdlrs; ++h )
      {
         CHECK_OKAY( SCIPconflicthdlrExec(set->conflicthdlrs[h], set->scip, tree->path[d-1], 
                        conflict->conflictvars, conflict->nconflictvars, *success, &result) );
         *success = *success || (result == SCIP_CONSADDED);
      }
   }

   return SCIP_OKAY;
}

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
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
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   Bool valid;
   int maxsize;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if propagation conflict analysis is enabled */
   if( !set->usepropconflict )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* calculate the maximal size of the conflict set */
   maxsize = (int)(set->maxconfvarsfac * prob->nbinvars);
   maxsize = MAX(maxsize, set->minmaxconfvars);
   if( maxsize < 2 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->propanalyzetime, set);

   conflict->npropcalls++;

   /* analyze the conflict set, and create a conflict constraint on success */
   CHECK_OKAY( conflictAnalyze(conflict, memhdr, set, stat, tree, maxsize, TRUE, &valid) );

   /* check, if a conflict constraint was created */
   if( valid )
   {
      conflict->npropconflicts++;
      if( success != NULL )
         *success = TRUE;
   }

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

/** gets number of valid conflicts detected in propagation conflict analysis */
Longint SCIPconflictGetNPropConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npropconflicts;
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
            sum_rhs += ABS(rowlhs);
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
            sum_rhs += ABS(rowrhs);
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
            sum_rhs += ABS(collb);
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
            sum_rhs += ABS(colub);
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
   if( SCIPlpiIsStable(conflict->lpi)
      && (SCIPlpiIsOptimal(conflict->lpi) || SCIPlpiIsPrimalUnbounded(conflict->lpi)) )
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
            CHECK_OKAY( conflictAnalyze(conflict, memhdr, set, stat, tree, maxsize, FALSE, success) );
         }
      }
      
      /* free memory buffer for primal solution of alternative LP */
      SCIPsetFreeBufferArray(set, &psol);      
   }
   else
   {
      infoMessage(set->verblevel, SCIP_VERBLEVEL_FULL, 
         "(node %lld) alternative LP is unstable or primal infeasible\n", stat->nnodes);
   }

   /* free memory buffers for alternative LP information */
   SCIPsetFreeBufferArray(set, &alt_islower);
   SCIPsetFreeBufferArray(set, &alt_colvars);
   SCIPsetFreeBufferArray(set, &alt_rowvars);
  
   return SCIP_OKAY;
}
#endif

/** analyzes an infeasible LP and unfixes additional binary variables while staying infeasible */
static
RETCODE unfixBinariesDualfarkas(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< LP data */
   Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*             nbdchgs,            /**< pointer to count total number of LP column bound changes */
   int*             bdchginds,          /**< array to store indices of bound changed columns */
   Real*            bdchgoldlbs,        /**< array to store old lower bounds of bound changed columns */
   Real*            bdchgoldubs,        /**< array to store old upper bounds of bound changed columns */
   Real*            bdchgnewlbs,        /**< array to store new lower bounds of bound changed columns */
   Real*            bdchgnewubs,        /**< array to store new upper bounds of bound changed columns */
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
   assert(valid != NULL);
   assert(resolve != NULL);

   debugMessage("unfixing binary variables in infeasible LP: cutoff=%g, depth=%d\n", 
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
      if( !row->local && !SCIPsetIsFeasZero(set, dualfarkas[r]) )
      {
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
   }

   /* calculate the current farkas activity, always using the best bound w.r.t. the farkas coefficient */
   farkasact = 0.0;
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      assert(SCIPvarGetProbindex(var) == v);

      /* ignore coefs close to 0.0 */
      if( SCIPsetIsZero(set, farkascoefs[v]) )
         farkascoefs[v] = 0.0;
      else if( farkascoefs[v] > 0.0 )
      {
         if( SCIPsetIsInfinity(set, curvarubs[v]) )
            goto TERMINATE;
         farkasact += farkascoefs[v] * curvarubs[v];
      }
      else
      {
         if( SCIPsetIsInfinity(set, -curvarlbs[v]) )
            goto TERMINATE;
         farkasact += farkascoefs[v] * curvarlbs[v];
      }
   }
   debugMessage("farkaslhs=%g, farkasact=%g\n", farkaslhs, farkasact);

   /* check, if the farkas row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, farkaslhs, farkasact) )
   {
      VAR** binvars;
      Real* binvarscores;
      Real score;
      Real farkascoef;
      int nfixings;
      Bool unfix;

      /* calculate the order in which the binary variables are tried to be unfixed:
       *  - prefer variables that have been fixed deeper in the tree, to get a more global conflict
       *  - prefer variables with small farkas coefficient to get rid of as many variables as possible
       */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &binvars, prob->nbinvars) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &binvarscores, prob->nbinvars) );
      nfixings = 0;
      for( v = 0; v < prob->nbinvars; ++v )
      {
         var = vars[v];

         /* ignore already unfixed variables */
         if( curvarlbs[v] < 0.5 && curvarubs[v] > 0.5 )
            continue;

         /* calculate score of variable fixing */
         farkascoef = ABS(farkascoefs[v]);
         score = 1.0 - farkascoef/(farkaslhs - farkasact);
         score = MAX(score, 0.0) + 1e-6;
         score *= SCIPvarGetFixDepth(var);

         /* insert fixed variable in candidate list */
         for( i = nfixings; i > 0 && score > binvarscores[i-1]; --i )
         {
            binvars[i] = binvars[i-1];
            binvarscores[i] = binvarscores[i-1];
         }
         binvars[i] = var;
         binvarscores[i] = score;
         nfixings++;
      }

      /* try to unfix binary variables and still keep the farkas row violated:
       * binaries fixed to 0.0 can be unfixed
       *  - if farkas value is positive and farkasact + farkascoef < farkaslhs (and farkasact has to be updated)
       *  - if farkas value is non-positive (because 0.0 is best farkas bound anyways)
       * binaries fixed to 1.0 can be unfixed
       *  - if farkas value is negative and farkasact - farkascoef < farkaslhs (and farkasact has to be updated)
       *  - if farkas value is non-negative (because 1.0 is best farkas bound anyways)
       */
      for( i = 0; i < nfixings; ++i )
      {
         var = binvars[i];
         v = SCIPvarGetProbindex(var);
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         assert(0 <= v && v < prob->nbinvars);

         debugMessage(" <%s> [%g,%g]: farkascoef=%g, %g <= %g\n",
            SCIPvarGetName(var), curvarlbs[v], curvarubs[v], farkascoefs[v], farkaslhs, farkasact);

         unfix = FALSE;
         if( curvarubs[v] < 0.5 )
         {
            /* variable is fixed to 0.0 */
            if( SCIPsetIsPositive(set, farkascoefs[v]) )
            {
               if( SCIPsetIsFeasGT(set, farkaslhs, farkasact + farkascoefs[v]) )
               {
                  /* we can unfix the variable but have to increase the farkas activity */
                  farkasact += farkascoefs[v];
                  unfix = TRUE;
                  *resolve = TRUE;
               }
            }
            else
               unfix = TRUE;
         }
         else
         {
            /* variable is fixed to 1.0 */
            assert(curvarlbs[v] > 0.5);
            if( SCIPsetIsNegative(set, farkascoefs[v]) )
            {
               if( SCIPsetIsFeasGT(set, farkaslhs, farkasact - farkascoefs[v]) )
               {
                  /* we can unfix the variable but have to increase the farkas activity */
                  farkasact -= farkascoefs[v];
                  unfix = TRUE;
                  *resolve = TRUE;
               }
            }
            else
               unfix = TRUE;
         }

         if( unfix )
         {
            if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
            {
               COL* col;
               int c;
               
               col = SCIPvarGetCol(var);
               c = SCIPcolGetLPPos(col);
               if( c >= 0 )
               {
                  assert(c < lp->ncols);
                  bdchginds[*nbdchgs] = c;
                  bdchgoldlbs[*nbdchgs] = curvarlbs[v];
                  bdchgoldubs[*nbdchgs] = curvarubs[v];
                  bdchgnewlbs[*nbdchgs] = 0.0;
                  bdchgnewubs[*nbdchgs] = 1.0;
                  (*nbdchgs)++;
               }
            }
            curvarlbs[v] = 0.0;
            curvarubs[v] = 1.0;
            debugMessage("  -> unfixed variable -> %g <= %g\n", farkaslhs, farkasact);
         }
      }

      /* free the buffer for the sorted binary variables */
      SCIPsetFreeBufferArray(set, &binvarscores);
      SCIPsetFreeBufferArray(set, &binvars);

      *valid = TRUE;
   }

 TERMINATE:

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &farkascoefs);
   SCIPsetFreeBufferArray(set, &dualfarkas);

   return SCIP_OKAY;
}

/** analyzes an LP exceeding the objective limit and unfixes additional binary variables while staying beyond the
 *  objective limit
 */
static
RETCODE unfixBinariesDualsol(
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< LP data */
   Real*            curvarlbs,          /**< current lower bounds of active problem variables */
   Real*            curvarubs,          /**< current upper bounds of active problem variables */
   int*             nbdchgs,            /**< pointer to count total number of LP column bound changes */
   int*             bdchginds,          /**< array to store indices of bound changed columns */
   Real*            bdchgoldlbs,        /**< array to store old lower bounds of bound changed columns */
   Real*            bdchgoldubs,        /**< array to store old upper bounds of bound changed columns */
   Real*            bdchgnewlbs,        /**< array to store new lower bounds of bound changed columns */
   Real*            bdchgnewubs,        /**< array to store new upper bounds of bound changed columns */
   Bool*            valid,              /**< pointer to store whether the unfixings are valid */
   Bool*            resolve             /**< pointer to store whether the changed LP should be resolved again */
   )
{
   LPI* lpi;
   ROW** rows;
   VAR** vars;
   ROW* row;
   COL* col;
   VAR* var;
   Real* primsols;
   Real* dualsols;
   Real* redcosts;
   Real* dualcoef;
   Real* varredcosts;
   Real duallhs;
   Real dualact;
   Real duallhsdelta;
   Real dualactdelta;
   int nrows;
   int ncols;
   int nvars;
   int r;
   int c;
   int v;
   int i;
   Bool havelocal;

   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(valid != NULL);
   assert(resolve != NULL);

   *valid = FALSE;
   *resolve = FALSE;

   debugMessage("unfixing binary variables in LP exceeding cutoff: cutoff=%g, depth=%d\n", 
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
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &dualcoef, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &varredcosts, nvars) );

   /* get solution from LPI */
   CHECK_OKAY( SCIPlpiGetSol(lpi, NULL, primsols, dualsols, NULL, redcosts) );

   /* calculate the dual row: y'^T{lhs,rhs} + r^T{lb,ub} - c* <= (y'^T A + r^T - c^T) x
    *  - y'i := 0 for local rows, y'i := yi for global rows
    *  - if yi > 0, use left hand side, if yi < 0, use right hand side
    *  - if ri > 0, use lower bound,    if ri < 0, use upper bound
    *  - c* is the primal bound (the LP's cutoffbound)
    *
    * The row coefficients transform into:
    *  y'^T A + r^T - c^T
    *   = y^T A + r^T - c^T - z^T A     (with zi := 0 for global, zi := yi for local rows)
    *   = -z^T A                        (because y^T A + r^T == c^T due to dual feasibility)
    *
    * The resulting dual row is: y'^T{lhs,rhs} + r^T{lb,ub} - c* <= -z^T A x
    */
   clearMemoryArray(dualcoef, nvars);
   duallhs = -lp->cutoffbound;

   /* calculate y'^T{lhs,rhs} and -z^T A x */
   havelocal = FALSE;
   for( r = 0; r < nrows; ++r )
   {
      row = rows[r];
      assert(row != NULL);
      assert(row->len == 0 || row->cols != NULL);
      assert(row->len == 0 || row->vals != NULL);
      assert(row == lp->lpirows[r]);

      /* ignore dual solution values of 0.0 */
      if( SCIPsetIsFeasZero(set, dualsols[r]) )
         continue;

      /* check dual feasibility */
      if( (SCIPsetIsInfinity(set, -row->lhs) && SCIPsetIsFeasPositive(set, dualsols[r]))
         || (SCIPsetIsInfinity(set, row->rhs) && SCIPsetIsFeasNegative(set, dualsols[r])) )
      {
         debugMessage("infeasible dual solution %g in row <%s>: lhs=%g, rhs=%g\n",
            dualsols[r], SCIProwGetName(row), row->lhs, row->rhs);
         goto TERMINATE;
      }

      /* local rows add up (negatively) to the dual row, global rows add up to the left hand side */
      if( row->local )
      {
         /* add negated local row coefficients to dual row */
         for( i = 0; i < row->len; ++i )
         {
            v = SCIPvarGetProbindex(SCIPcolGetVar(row->cols[i]));
            assert(0 <= v && v < nvars);
            dualcoef[v] -= dualsols[r] * row->vals[i];
         }
         havelocal = TRUE;
         debugMessage(" local row <%s>: dual=%g\n", SCIProwGetName(row), dualsols[r]);
      }
      else
      {
         /* add row side to dual row lhs: dualsol > 0 -> lhs, dualsol < 0 -> rhs */
         if( dualsols[r] > 0.0 )
         {
            assert(!SCIPsetIsInfinity(set, -row->lhs));
            duallhs += dualsols[r] * (row->lhs - row->constant);
            debugMessage(" global row <%s>: lhs=%g, dual=%g -> %g\n", 
               SCIProwGetName(row), row->lhs - row->constant, dualsols[r], duallhs);
         }
         else
         {
            assert(!SCIPsetIsInfinity(set, row->rhs));
            duallhs += dualsols[r] * (row->rhs - row->constant);
            debugMessage(" global row <%s>: rhs=%g, dual=%g -> %g\n", 
               SCIProwGetName(row), row->rhs - row->constant, dualsols[r], duallhs);
         }
      }
   }

   /* calculate r^T{lb,ub}
    * loose variables have reduced costs r_i = c_i, s.t. the LP's loose objective value is equal to
    * r^T{lb,ub} for the loose variables
    */
   duallhs += lp->looseobjval;
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];

      if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE )
      {
         /* reduced costs for loose variables are equal to the objective value */
         varredcosts[v] = SCIPvarGetObj(var);
         continue;
      }

      col = SCIPvarGetCol(var);
      c = SCIPcolGetLPPos(col);
      assert(c == -1 || col == lp->cols[c]);
      assert(c == -1 || col == lp->lpicols[c]);

      v = SCIPvarGetProbindex(var);
      assert(0 <= v && v < nvars);
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, curvarlbs[v], SCIPvarGetLbGlobal(var)));
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, curvarubs[v], SCIPvarGetUbGlobal(var)));

      /* get reduced costs from LPI, or calculate it manually if the column is not in current LP */
      if( c == -1 )
         varredcosts[v] = SCIPcolCalcRedcost(col, dualsols);
      else
      {
#ifndef NDEBUG
         varredcosts[v] = SCIPcolCalcRedcost(col, dualsols);
         assert(SCIPsetIsSumEQ(set, varredcosts[v], redcosts[c]));
#endif
         varredcosts[v] = redcosts[c];
      }

      /* ignore reduced costs close to 0.0 */
      if( SCIPsetIsFeasZero(set, varredcosts[v]) )
      {
         varredcosts[v] = 0.0;
         continue;
      }

      /* check dual feasibility */
      if( (SCIPsetIsGT(set, primsols[c], curvarlbs[v]) && SCIPsetIsFeasPositive(set, varredcosts[v]))
         || (SCIPsetIsLT(set, primsols[c], curvarubs[v]) && SCIPsetIsFeasNegative(set, varredcosts[v])) )
      {
         debugMessage("infeasible reduced costs %g in var <%s>: lb=%g, ub=%g\n",
            varredcosts[v], SCIPvarGetName(var), curvarlbs[v], curvarubs[v]);
         goto TERMINATE;
      }

      /* add column's bound to dual row lhs: redcost > 0 -> lb, redcost < 0 -> ub */
      if( varredcosts[v] > 0.0 )
      {
         if( SCIPsetIsInfinity(set, -curvarlbs[v]) )
            goto TERMINATE;
         duallhs += varredcosts[v] * curvarlbs[v];
      }
      else
      {
         if( SCIPsetIsInfinity(set, curvarubs[v]) )
            goto TERMINATE;
         duallhs += varredcosts[v] * curvarubs[v];
      }
      debugMessage(" col <%s> [%g,%g]: redcost=%g -> %g\n", 
         SCIPvarGetName(var), curvarlbs[v], curvarubs[v], varredcosts[v], duallhs);
   }
   debugMessage("duallhs=%g, localrows=%d\n", duallhs, havelocal);

   /* calculate the current dual activity, always using the best bound w.r.t. the dual coefficient;
    * if no local rows are included, the dual row is empty
    */
   dualact = 0.0;
   if( havelocal )
   {
      for( v = 0; v < nvars; ++v )
      {
         var = vars[v];
         assert(SCIPvarGetProbindex(var) == v);
         
         /* ignore coefs close to 0.0 */
         if( SCIPsetIsZero(set, dualcoef[v]) )
            dualcoef[v] = 0.0;
         else if( dualcoef[v] > 0.0 )
            dualact += dualcoef[v] * curvarubs[v];
         else
            dualact += dualcoef[v] * curvarlbs[v];
      }
   }
   debugMessage("duallhs=%g, dualact=%g, cutoffbound=%g\n", duallhs, dualact, lp->cutoffbound);

   /* check, if the dual row is still violated (using current bounds and ignoring local rows) */
   if( SCIPsetIsFeasGT(set, duallhs, dualact) )
   {
      VAR** binvars;
      Real* binvarscores;
      Real score;
      Real redcost;
      int nfixings;
      int sign;

      /* calculate the order in which the binary variables are tried to be unfixed:
       *  - prefer variables that have been fixed deeper in the tree, to get a more global conflict
       *  - prefer variables with small reduced costs to get rid of as many variables as possible
       */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &binvars, prob->nbinvars) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &binvarscores, prob->nbinvars) );
      nfixings = 0;
      for( v = 0; v < prob->nbinvars; ++v )
      {
         var = prob->vars[v];

         /* ignore already unfixed variables */
         if( curvarlbs[v] < 0.5 && curvarubs[v] > 0.5 )
            continue;

         /* calculate score of variable fixing */
         redcost = ABS(varredcosts[v]);
         score = 1.0 - redcost/(duallhs - dualact);
         score = MAX(score, 0.0) + 1e-6;
         score *= SCIPvarGetFixDepth(var);

         /* insert fixed variable in candidate list */
         for( i = nfixings; i > 0 && score > binvarscores[i-1]; --i )
         {
            binvars[i] = binvars[i-1];
            binvarscores[i] = binvarscores[i-1];
         }
         binvars[i] = var;
         binvarscores[i] = score;
         nfixings++;
      }

      /* try to unfix binary variables and still keep the dual row  y'^T{lhs,rhs} + r^T{lb,ub} - c* <= -z^T A x  violated:
       * binary variables i fixed to x_i == 0.0:
       *  - set duallhsdelta := r_i if r_i < 0.0,             and duallhsdelta := 0.0, if r_i >= 0.0
       *  - set dualactdelta := dualcoef, if dualcoef > 0.0,  and dualactdelta := 0.0, if dualcoef <= 0.0
       * binary variables i fixed to x_i == 1.0:
       *  - set duallhsdelta := -r_i if r_i > 0.0,            and duallhsdelta := 0.0, if r_i <= 0.0
       *  - set dualactdelta := -dualcoef, if dualcoef < 0.0, and dualactdelta := 0.0, if dualcoef >= 0.0
       * if duallhs + duallhsdelta > dualact + dualactdelta:
       *  - variable can be unfixed
       *  - update duallhs and dualact with duallhsdelta and dualactdelta
       * if duallhs + duallhsdelta <= dualact + dualactdelta:
       *  - variable cannot be unfixed and has to be put in the conflict set
       */
      for( i = 0; i < nfixings; ++i )
      {
         var = binvars[i];
         v = SCIPvarGetProbindex(var);
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         assert(0 <= v && v < prob->nbinvars);

         debugMessage(" <%s> [%g,%g]: score=%g, infdepth=%d, redcost=%g, dualcoef=%g, %g <= %g\n",
            SCIPvarGetName(var), curvarlbs[v], curvarubs[v], binvarscores[i], SCIPvarGetFixDepth(var),
            varredcosts[v], dualcoef[v], duallhs, dualact);

         /* variable fixed to lower bound: sign = +1, upper bound: sign = -1 */
         sign = curvarubs[v] < 0.5 ? +1 : -1;
         
         /* calculate duallhsdelta and dualactdelta for unfixing the variable */
         duallhsdelta = sign * varredcosts[v];
         duallhsdelta = MIN(duallhsdelta, 0.0);
         dualactdelta = sign * dualcoef[v];
         dualactdelta = MAX(dualactdelta, 0.0);
         assert(SCIPsetIsLE(set, duallhsdelta, dualactdelta));

         /* check, if variable can be unfixed */
         if( SCIPsetIsFeasGT(set, duallhs + duallhsdelta, dualact + dualactdelta) )
         {
            /* we can unfix the variable and have to update the dual row's left hand side and activity */
            duallhs += duallhsdelta;
            dualact += dualactdelta;
            if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
            {
               col = SCIPvarGetCol(var);
               c = SCIPcolGetLPPos(col);
               if( c >= 0 )
               {
                  assert(c < ncols);
                  bdchginds[*nbdchgs] = c;
                  bdchgoldlbs[*nbdchgs] = curvarlbs[v];
                  bdchgoldubs[*nbdchgs] = curvarubs[v];
                  bdchgnewlbs[*nbdchgs] = 0.0;
                  bdchgnewubs[*nbdchgs] = 1.0;
                  (*nbdchgs)++;
               }
            }
            curvarlbs[v] = 0.0;
            curvarubs[v] = 1.0;
            *resolve = *resolve || SCIPsetIsLT(set, duallhsdelta, dualactdelta);
            debugMessage("  -> unfixed variable <%s> -> %g <= %g\n", SCIPvarGetName(var), duallhs, dualact);
         }
      }

      /* free the buffer for the sorted binary variables */
      SCIPsetFreeBufferArray(set, &binvarscores);
      SCIPsetFreeBufferArray(set, &binvars);

      *valid = TRUE;
   }

 TERMINATE:

   /* free temporary memory */
   SCIPsetFreeBufferArray(set, &varredcosts);
   SCIPsetFreeBufferArray(set, &dualcoef);
   SCIPsetFreeBufferArray(set, &redcosts);
   SCIPsetFreeBufferArray(set, &dualsols);
   SCIPsetFreeBufferArray(set, &primsols);

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
   int*             iterations,         /**< pointer to store the total number of LP iterations used */
   Bool*            success             /**< pointer to store whether a conflict constraint was created */
   )
{
   LPI* lpi;
   COL** cols;
   VAR** vars;
   VAR* var;
   int* bdchginds;
   Real* bdchgoldlbs;
   Real* bdchgoldubs;
   Real* bdchgnewlbs;
   Real* bdchgnewubs;
   Real* curvarlbs;
   Real* curvarubs;
   Real loclb;
   Real locub;
   Real objval;
   Bool valid;
   Bool resolve;
   int nbdchgs;
   int lastnbdchgs;
   int ncols;
   int nvars;
   int maxsize;
   int iter;
   int nloops;
   int v;

   assert(set != NULL);
   assert(prob != NULL);
   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(iterations != NULL);
   assert(success != NULL);

   debugMessage("analyzing conflict on infeasible LP\n");

   *iterations = 0;
   *success = FALSE;

   /* conflict analysis only makes sense, if binary variables exist */
   if( prob->nbinvars == 0 )
      return SCIP_OKAY;

   /* calculate the maximal size of the conflict set */
   maxsize = (int)(set->maxconfvarsfac * prob->nbinvars);
   maxsize = MAX(maxsize, set->minmaxconfvars);
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

   /* get LP columns */
   cols = SCIPlpGetCols(lp);
   ncols = SCIPlpGetNCols(lp);

   /* get active problem variables */
   vars = prob->vars;
   nvars = prob->nvars;

   /* get temporary memory for remembering variables' current bounds */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &curvarlbs, nvars) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &curvarubs, nvars) );

   /* the following algorithm is used to find a subset of fixed binary variables leading to an infeasible LP:
    * 1. unfix all non-binary variables
    *    -> store changes in bdchg and curvarlbs/ubs arrays
    *    -> apply changes to the LPI
    * 2. resolve LP
    * 3. call unfixBinariesDualfarkas() or unfixBinariesDualsol()
    *    -> store additional changes in bdchg and curvarlbs/ubs arrays
    *    -> apply additional changes to the LPI
    * 4. if additional unfixings were found, goto 2
    * 5. undo all bound changes in the LPI
    * 6. analyze conflict
    *    -> put remaining fixed binaries (see curvarlbs/ubs arrays) into starting conflict set
    */

   /* get current bounds of binary variables */
   valid = FALSE;
   for( v = 0; v < prob->nbinvars; ++v )
   {
      var = vars[v];
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);

      curvarlbs[v] = getLbLocal(var);
      curvarubs[v] = getUbLocal(var);
      valid = valid || (curvarlbs[v] > 0.5 || curvarubs[v] < 0.5);
      debugMessage("binary variable <%s> [%g,%g]\n", SCIPvarGetName(var), curvarlbs[v], curvarubs[v]);
   }

   /* conflict analysis only makes sense, if fixed binary variables exist */
   if( valid )
   {
      /* temporarily remove objective limit in LP solver */
      CHECK_OKAY( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, set->infinity) );

      /* get temporary memory for remembering bound changes on LPI columns */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchginds, ncols) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchgoldlbs, ncols) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchgoldubs, ncols) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchgnewlbs, ncols) );
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &bdchgnewubs, ncols) );

      /* get current bounds of non-binary variables, and unfix them:
       * compare local and global bounds of non-binary variables and remember differences
       */
      nbdchgs = 0;
      for( v = prob->nbinvars; v < nvars; ++v )
      {
         var = vars[v];
         assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);

         loclb = getLbLocal(var);
         locub = getUbLocal(var);
         curvarlbs[v] = SCIPvarGetLbGlobal(var);
         curvarubs[v] = SCIPvarGetUbGlobal(var);

         debugMessage("non-binary <%s> [%g,%g] -> [%g,%g]\n", SCIPvarGetName(var), 
            loclb, locub, curvarlbs[v], curvarubs[v]);

         if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
         {
            COL* col;
            int c;
         
            col = SCIPvarGetCol(var);
            c = SCIPcolGetLPPos(col);
            if( c >= 0 && (!SCIPsetIsEQ(set, curvarlbs[v], loclb) || !SCIPsetIsEQ(set, curvarubs[v], locub)) )
            {
               assert(c < ncols);
               assert(cols[c] == col);
               assert(lp->lpicols[c] == col);
               bdchginds[nbdchgs] = c;
               bdchgoldlbs[nbdchgs] = loclb;
               bdchgoldubs[nbdchgs] = locub;
               bdchgnewlbs[nbdchgs] = curvarlbs[v];
               bdchgnewubs[nbdchgs] = curvarubs[v];
               nbdchgs++;
            }
         }
      }

      /* apply changes of non-binary variables to the LP solver */
      if( nbdchgs > 0 )
      {
         CHECK_OKAY( SCIPlpiChgBounds(lpi, nbdchgs, bdchginds, bdchgnewlbs, bdchgnewubs) );
      }

      /* unfix binary variables until no more unfixings could be found */
      nloops = 0;
      do
      {
         nloops++;
         resolve = FALSE;
         lastnbdchgs = nbdchgs;

         /* resolve LP */
         CHECK_OKAY( SCIPlpiSolveDual(lpi) );
      
         /* evaluate result */
         if( SCIPlpiIsOptimal(lpi) )
         {
            CHECK_OKAY( SCIPlpiGetObjval(lpi, &objval) );
            valid = (objval >= lp->lpiuobjlim);
         }
         else
            valid = SCIPlpiIsPrimalInfeasible(lpi);

         /* count number of LP iterations */
         CHECK_OKAY( SCIPlpiGetIntpar(lpi, SCIP_LPPAR_LPITER, &iter) );
         (*iterations) += iter;

         if( valid )
         {
            /* unfix additional binary variables */
            if( SCIPlpiIsPrimalInfeasible(lpi) )
            {
               CHECK_OKAY( unfixBinariesDualfarkas(set, stat, prob, lp, curvarlbs, curvarubs,
                     &nbdchgs, bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, &valid, &resolve) );
            }
            else
            {
               assert(SCIPlpiIsOptimal(lpi));
               CHECK_OKAY( unfixBinariesDualsol(set, stat, prob, lp, curvarlbs, curvarubs,
                     &nbdchgs, bdchginds, bdchgoldlbs, bdchgoldubs, bdchgnewlbs, bdchgnewubs, &valid, &resolve) );
            }

            /* apply additional changes of binary variables to the LP solver */
            if( valid && nbdchgs > lastnbdchgs )
            {
               CHECK_OKAY( SCIPlpiChgBounds(lpi, nbdchgs - lastnbdchgs, 
                     &bdchginds[lastnbdchgs], &bdchgnewlbs[lastnbdchgs], &bdchgnewubs[lastnbdchgs]) );
            }
         }
         assert(nbdchgs >= lastnbdchgs);
         assert(!resolve || (nbdchgs > lastnbdchgs && valid));
      }
      while( resolve && nloops < 100 );

      /* reset variables to local bounds */
      CHECK_OKAY( SCIPlpiChgBounds(lpi, nbdchgs, bdchginds, bdchgoldlbs, bdchgoldubs) );
   
      /* reset objective limit in LP solver */
      CHECK_OKAY( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_UOBJLIM, lp->lpiuobjlim) );

      if( valid )
      {
         /* initialize conflict data */
         CHECK_OKAY( SCIPconflictInit(conflict) );

         /* add remaining fixed binary variables to conflit set */
         debugMessage("initial conflict set after LP analysis:\n");
         for( v = 0; v < nvars; ++v )
         {
            var = vars[v];
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, curvarlbs[v], SCIPvarGetLbGlobal(var)));
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPsetIsEQ(set, curvarubs[v], SCIPvarGetUbGlobal(var)));

            if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY )
            {
               if( curvarlbs[v] > 0.5 || curvarubs[v] < 0.5 )
               {
                  debugMessage("   <%s> [%g,%g]\n", SCIPvarGetName(var), curvarlbs[v], curvarubs[v]);
                  CHECK_OKAY( SCIPconflictAddVar(conflict, var) );
               }
            }
         }

         /* analyze the conflict set, and create a conflict constraint on success */
         CHECK_OKAY( conflictAnalyze(conflict, memhdr, set, stat, tree, maxsize, FALSE, success) );
      }

      /* free temporary memory */
      SCIPsetFreeBufferArray(set, &bdchgnewubs);
      SCIPsetFreeBufferArray(set, &bdchgnewlbs);
      SCIPsetFreeBufferArray(set, &bdchgoldubs);
      SCIPsetFreeBufferArray(set, &bdchgoldlbs);
      SCIPsetFreeBufferArray(set, &bdchginds);
   }

   /* free temporary memory */
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
   Bool found;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(lp != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->uselpconflict )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->lpanalyzetime, set);
   conflict->nlpcalls++;

   /* perform conflict analysis */
   CHECK_OKAY( conflictAnalyzeLP(conflict, memhdr, set, stat, prob, tree, lp, &iterations, &found) );
   conflict->nlpiterations += iterations;
   if( found )
   {
      conflict->nlpconflicts++;
      if( success != NULL )
         *success = TRUE;
   }

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

/** gets number of valid conflicts detected in infeasible LP conflict analysis */
Longint SCIPconflictGetNLPConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nlpconflicts;
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
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   int* cstat;
   int* rstat;
   Real oldlb;
   Real oldub;
   Real newlb;
   Real newub;
   Bool found;
   int iterations;

   assert(lp != NULL);
   assert(lp->flushed);
   assert(lp->solved);
   assert(col != NULL);
   assert((col->strongbranchdown >= lp->cutoffbound && SCIPsetCeil(set, col->primsol-1.0) >= col->lb - 0.5)
      || (col->strongbranchup >= lp->cutoffbound && SCIPsetFloor(set, col->primsol+1.0) <= col->ub + 0.5));
   
   if( success != NULL )
      *success = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->usesbconflict )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
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

   /* is down branch infeasible? */
   if( col->strongbranchdown >= lp->cutoffbound )
   {
      newub = SCIPsetCeil(set, col->primsol-1.0);
      if( newub >= col->lb - 0.5 )
      {
         /* change the upper bound */
         col->ub = newub;
         CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );
            
         /* resolve the LP */
         CHECK_OKAY( SCIPlpiSolveDual(lp->lpi) );
            
         /* perform conflict analysis on infeasible LP */
         CHECK_OKAY( conflictAnalyzeLP(conflict, memhdr, set, stat, prob, tree, lp, &iterations, &found) );
         conflict->nsbiterations += iterations;
         if( found )
         {
            conflict->nsbconflicts++;
            if( success != NULL )
               *success = TRUE;
         }

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
         /* change the lower bound */
         col->lb = newlb;
         CHECK_OKAY( SCIPlpiChgBounds(lp->lpi, 1, &col->lpipos, &col->lb, &col->ub) );
            
         /* resolve the LP */
         CHECK_OKAY( SCIPlpiSolveDual(lp->lpi) );
            
         /* perform conflict analysis on infeasible LP */
         CHECK_OKAY( conflictAnalyzeLP(conflict, memhdr, set, stat, prob, tree, lp, &iterations, &found) );
         conflict->nsbiterations += iterations;
         if( found )
         {
            conflict->nsbconflicts++;
            if( success != NULL )
               *success = TRUE;
         }

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

/** gets number of valid conflicts detected in infeasible strong branching conflict analysis */
Longint SCIPconflictGetNStrongbranchConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nsbconflicts;
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

/** analyzes a pseudo solution with objective value exceeding the current cutoff to find out the bound changes on binary
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
   VAR* var;
   Real pseudoobjval;
   Real obj;
   Bool valid;
   int maxsize;
   int v;

   assert(conflict != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(prob != NULL);
   assert(lp != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if pseudo solution conflict analysis is enabled */
   if( !set->usepseudoconflict )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* calculate the maximal size of the conflict set */
   maxsize = (int)(set->maxconfvarsfac * prob->nbinvars);
   maxsize = MAX(maxsize, set->minmaxconfvars);
   if( maxsize < 2 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->pseudoanalyzetime, set);

   conflict->npseudocalls++;
   valid = FALSE;

   /* if any non-binary variable is not on its best bound in current pseudo solution, we have to abort */
   for( v = prob->nbinvars; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      assert(SCIPvarGetType(var) != SCIP_VARTYPE_BINARY);
      obj = SCIPvarGetObj(var);
      if( (obj > 0.0 && SCIPvarGetLbLocal(var) != SCIPvarGetLbGlobal(var))
         || (obj < 0.0 && SCIPvarGetUbLocal(var) != SCIPvarGetUbGlobal(var)) )
         break;
   }

   if( v == prob->nvars )
   {
      pseudoobjval = SCIPlpGetPseudoObjval(lp, set);
      assert(pseudoobjval >= lp->cutoffbound);
      
      debugMessage("analyzing pseudo solution (obj: %g) that exceeds objective limit (%g)\n",
         pseudoobjval, lp->cutoffbound);

      /* initialize conflict data */
      CHECK_OKAY( SCIPconflictInit(conflict) );
      
      /* add all binary variables, that are not on their best bound in current pseudo solution;
       * however, the variable need not to be added, if after unfixing the pseudo objective value is still too large
       */
      for( v = 0; v < prob->nbinvars; ++v )
      {
         var = prob->vars[v];
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
         obj = SCIPvarGetObj(var);
         if( obj > 0.0 && SCIPvarGetLbLocal(var) > 0.5 )
         {
            if( pseudoobjval - obj >= lp->cutoffbound )
               pseudoobjval -= obj;
            else
            {
               debugMessage("adding <%s>[%g,%g] (obj=%g) to conflict set\n", 
                  SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), obj);
               CHECK_OKAY( SCIPconflictAddVar(conflict, var) );
            }
         }
         else if( obj < 0.0 && SCIPvarGetUbLocal(var) < 0.5 )
         {
            if( pseudoobjval + obj >= lp->cutoffbound )
               pseudoobjval += obj;
            else
            {
               debugMessage("adding <%s>[%g,%g] (obj=%g) to conflict set\n", 
                  SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), obj);
               CHECK_OKAY( SCIPconflictAddVar(conflict, var) );
            }
         }
      }
      
      /* analyze the conflict set, and create a conflict constraint on success */
      CHECK_OKAY( conflictAnalyze(conflict, memhdr, set, stat, tree, maxsize, TRUE, &valid) );
      
      /* check, if a valid conflict was found */
      if( valid )
      {
         conflict->npseudoconflicts++;
         if( success != NULL )
            *success = TRUE;
      }
   }

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

/** gets number of valid conflicts detected in pseudo solution conflict analysis */
Longint SCIPconflictGetNPseudoConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->npseudoconflicts;
}

