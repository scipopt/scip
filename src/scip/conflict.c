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
#pragma ident "@(#) $Id: conflict.c,v 1.23 2004/01/22 14:42:26 bzfpfend Exp $"

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
 * of the branch-and-bound tree, that would lead to the same conflict.
 *
 *  1. Put all the given variables to a priority queue, which is ordered,
 *     such that the variable that was fixed last due to branching or deduction
 *     is at the top of the queue. The variables in the queue are represented
 *     by themselves or their negation in order to have all variable's representants
 *     currently fixed to FALSE.
 *     Create an empty conflict set.
 *  2. Remove the top variable v from the priority queue.
 *  3. (a) If the remaining queue is non-empty, and variable w (the one that is now
 *         on the top of the queue) was fixed at the same depth level as v, and if
 *         the assignment to v was a deduction with known inference reason, and if
 *         the inference constraint is globally valid:
 *          - resolve variable v by asking the constraint that infered the
 *            assignment of v to put all the variables that when assigned to their
 *            current values lead to the deduction of v on the priority queue.
 *            Note that these variables have at most the same inference depth
 *            level as variable v, and were deduced earlier than v.
 *            If the variable's current assignment is TRUE, it is automatically
 *            replaced by its negation in the priority queue.
 *     (b) Otherwise, the assignment to variable v was a branching decision or
 *         a deduction with missing inference reason.
 *          - Put v into the conflict set.
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
   const SET*       set,                /**< global SCIP settings */
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
      CHECK_OKAY( conflicthdlr->conflictexec(scip, conflicthdlr, conflictvars, nconflictvars, resolved, result) );

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
DECL_SORTPTRCOMP(conflictVarCmp)
{  /*lint --e{715}*/
   VAR* var1;
   VAR* var2;
   
   var1 = (VAR*)elem1;
   var2 = (VAR*)elem2;
   assert(var1 != NULL);
   assert(SCIPvarGetType(var1) == SCIP_VARTYPE_BINARY);
   assert(var2 != NULL);
   assert(SCIPvarGetType(var2) == SCIP_VARTYPE_BINARY);

   if( SCIPvarGetInferDepth(var1) > SCIPvarGetInferDepth(var2) )
      return -1;
   else if( SCIPvarGetInferDepth(var1) < SCIPvarGetInferDepth(var2) )
      return +1;
   else if( SCIPvarGetInferIndex(var1) > SCIPvarGetInferIndex(var2) )
      return -1;
   else if( SCIPvarGetInferIndex(var1) < SCIPvarGetInferIndex(var2) )
      return +1;
   else if( SCIPvarGetIndex(var1) > SCIPvarGetIndex(var2) )
      return +1;
   else if( SCIPvarGetIndex(var1) < SCIPvarGetIndex(var2) )
      return -1;
   else
      return 0;
}

/** creates conflict analysis data for propagation conflicts */
RETCODE SCIPconflictCreate(
   CONFLICT**       conflict,           /**< pointer to conflict analysis data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(conflict != NULL);

   ALLOC_OKAY( allocMemory(conflict) );

   CHECK_OKAY( SCIPclockCreate(&(*conflict)->analyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPpqueueCreate(&(*conflict)->varqueue, set->memgrowinit, set->memgrowfac, conflictVarCmp) );
   (*conflict)->conflictvars = NULL;
   (*conflict)->conflictvarssize = 0;
   (*conflict)->nconflictvars = 0;
   (*conflict)->ncalls = 0;
   (*conflict)->nconflicts = 0;
   
   return SCIP_OKAY;
}

/** frees conflict analysis data for propagation conflicts */
RETCODE SCIPconflictFree(
   CONFLICT**       conflict            /**< pointer to conflict analysis data */
   )
{
   assert(conflict != NULL);
   assert(*conflict != NULL);

   SCIPclockFree(&(*conflict)->analyzetime);
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

   SCIPpqueueClear(conflict->varqueue);
   conflict->nconflictvars = 0;

   return SCIP_OKAY;
}

/** adds variable to conflict variable candidates */
RETCODE SCIPconflictAddVar(
   CONFLICT*        conflict,           /**< conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   VAR*             var                 /**< problem variable */
   )
{
   assert(conflict != NULL);
   assert(var != NULL);
   assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
   assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   debugMessage("adding variable <%s>[%g,%g] [status:%d, depth:%d, num:%d, cons:%p] to conflict candidates\n",
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetStatus(var),
      SCIPvarGetInferDepth(var), SCIPvarGetInferIndex(var), SCIPvarGetInferCons(var));

   /* get active problem variable */
   var = SCIPvarGetProbvar(var);

   /* we can ignore fixed variables */
   if( var != NULL )
   {
      assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
      assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
      
      /* choose between variable or its negation, such that the literal is fixed to FALSE */
      if( SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), 1.0) )
      {
         /* variable is fixed to TRUE -> use the negation */
         CHECK_OKAY( SCIPvarNegate(var, memhdr, set, stat, &var) );
      }
      debugMessage(" -> active conflict candidate is <%s>[%g,%g] [status:%d]\n",
         SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetStatus(var));
      
      /* put candidate in priority queue */
      CHECK_OKAY( SCIPpqueueInsert(conflict->varqueue, (void*)var) );
   }

   return SCIP_OKAY;
}

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set
 */
static
RETCODE conflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   SET*             set,                /**< global SCIP settings */
   int              maxsize,            /**< maximal size of conflict set */
   Bool*            success             /**< pointer to store whether the conflict set is valid */
   )
{
   VAR* var;
   VAR* nextvar;
   VAR* infervar;
   CONS* infercons;

   assert(conflict != NULL);
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

   /* clear the conflict set */
   conflict->nconflictvars = 0;

   /* process all variables in the conflict candidate queue */
   var = (VAR*)(SCIPpqueueRemove(conflict->varqueue));
   while( var != NULL && conflict->nconflictvars < maxsize )
   {
#ifdef DEBUG
      {
         int v;
         debugMessage("processing conflict var <%s>\n", SCIPvarGetName(var));
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

      assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

      /* if the first variable on the remaining queue is equal to the actual variable,
       * this is a multiple insertion in the conflict candidate queue and we can ignore the current
       * variable
       */
      nextvar = (VAR*)(SCIPpqueueFirst(conflict->varqueue));
      if( var != nextvar )
      {
         /* check, if the variable can and should be resolved */
         infercons = SCIPvarGetInferCons(var);
         if( nextvar != NULL
            && infercons != NULL
            && SCIPconsIsGlobal(infercons)
            && SCIPvarGetInferDepth(var) == SCIPvarGetInferDepth(nextvar) )
         {
            /* resolve variable v by asking the constraint that infered the
             * assignment of v to put all the variables that when assigned to their
             * current values lead to the deduction of v on the priority queue.
             */
            infervar = SCIPvarGetInferVar(var);
            assert(infervar != NULL);
            assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(infervar), SCIPvarGetUbLocal(infervar)));
            debugMessage("resolving conflict var <%s>: constraint <%s> infered <%s> == %g at depth %d, num %d\n",
               SCIPvarGetName(var), SCIPconsGetName(infercons), SCIPvarGetName(infervar), SCIPvarGetLbLocal(infervar),
               SCIPvarGetInferDepth(var), SCIPvarGetInferIndex(var));
            CHECK_OKAY( SCIPconsResolveConflictVar(infercons, set, SCIPvarGetInferVar(var)) );
         }
         else
         {
            assert(nextvar == NULL || infercons == NULL || SCIPvarGetInferDepth(var) > SCIPvarGetInferDepth(nextvar));
            assert(nextvar == NULL || SCIPvarGetInferDepth(var) >= SCIPvarGetInferDepth(nextvar));
            assert(SCIPsetIsFeasEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
            assert(SCIPsetIsFeasZero(set, SCIPvarGetLbLocal(var)));

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

   /* if a valid conflict set was found, call the conflict handlers */
   if( var == NULL )
   {
      RESULT result;
      int h;

      assert(conflict->conflictvars != NULL);
      assert(conflict->nconflictvars > 0);

      debugMessage(" -> found conflict set of size %d/%d\n", conflict->nconflictvars, maxsize);

      /* sort conflict handlers by priority */
      SCIPsetSortConflicthdlrs(set);

      /* call conflict handlers to create a conflict constraint */
      for( h = 0; h < set->nconflicthdlrs; ++h )
      {
         CHECK_OKAY( SCIPconflicthdlrExec(set->conflicthdlrs[h], set->scip, 
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
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
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
   maxsize = (int)(set->maxconfvarsfac * prob->nbin);
   maxsize = MAX(maxsize, set->minmaxconfvars);
   if( maxsize < 2 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(conflict->analyzetime, set);

   conflict->ncalls++;

   /* analyze the conflict set, and create a conflict constraint on success */
   CHECK_OKAY( conflictAnalyze(conflict, set, maxsize, &valid) );

   /* check, if a conflict constraint was created */
   if( valid )
   {
      conflict->nconflicts++;
      if( success != NULL )
         *success = TRUE;
   }

   /* stop timing */
   SCIPclockStop(conflict->analyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing propagation conflicts */
Real SCIPconflictGetTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->analyzetime);
}

/** gets number of calls to propagation conflict analysis */
Longint SCIPconflictGetNCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ncalls;
}

/** gets number of valid conflicts detected in propagation conflict analysis */
Longint SCIPconflictGetNConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nconflicts;
}




/*
 * Infeasible LP Conflict Analysis
 */

/** Generate the alternative lp from the current lp */
/**@todo Correct hard coded 1000.0 */
static
RETCODE lpGenerateAltLP(
   LP*              lp,                 /**< LP data */
   const SET*       set,                /**< global SCIP settings */
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

   lpi_infinity = SCIPlpiInfinity(lp->lpi);
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
      row = rows[i];
      assert( row != NULL );

      /* left hand side */
      if( !SCIPsetIsInfinity(set, -row->lhs) )
      {
         cnt = 0;
         for( j = 0; j < row->len; ++j )
         {
            inds[cnt] = row->cols[j]->lppos;
            vals[cnt] = -row->vals[j];
            ++cnt;
         }
         if( !SCIPsetIsZero(set, row->lhs) )
         {
            inds[cnt] = ncols;
            vals[cnt] = -(row->lhs - row->constant);
            sum_rhs += ABS(row->lhs);
            ++cnt;
         }
         if( row->local )
            obj = 2000.0;
         else
            obj = 0.0;
         lb = 0.0;
         beg = 0;
         CHECK_OKAY( SCIPlpiAddCols(alt_lpi, 1, &obj, &lb, &lpi_infinity, NULL, cnt, &beg, inds, vals) );
         alt_rowvars[*alt_nrowvars] = row;
         ++(*alt_nrowvars);
      }

      /* right hand side */
      if( !SCIPsetIsInfinity(set, row->rhs) )
      {
         cnt = 0;
         for( j = 0; j < row->len; ++j )
         {
            inds[cnt] = row->cols[j]->lppos;
            vals[cnt] = row->vals[j];
            ++cnt;
         }
         if( !SCIPsetIsZero(set, row->rhs) )
         {
            inds[cnt] = ncols;
            vals[cnt] = row->rhs - row->constant;
            sum_rhs += ABS(row->rhs);
            ++cnt;
         }
         if( row->local )
            obj = 3000.0;
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
      col = cols[j];
      assert( col != NULL);
      /* lower bounds */
      if( !SCIPsetIsInfinity(set, -col->lb) )
      {
         cnt = 0;
         inds[cnt] = j;
         vals[cnt] = -1.0;
         ++cnt;

         if( !SCIPsetIsZero(set, col->lb) )
         {
            inds[cnt] = ncols;
            vals[cnt] = -col->lb;
            sum_rhs += ABS(col->lb);
            ++cnt;
         }
         if( !SCIPsetIsEQ(set, SCIPvarGetLbGlobal(col->var), col->lb) )
         {
            if( SCIPvarGetType(col->var) == SCIP_VARTYPE_BINARY )
               obj = 1.0;
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
      if( !SCIPsetIsInfinity(set, col->ub) )
      {
         cnt = 0;
         inds[cnt] = j;
         vals[cnt] = 1.0;
         ++cnt;

         if( !SCIPsetIsZero(set, col->ub) )
         {
            inds[cnt] = ncols;
            vals[cnt] = col->ub;
            sum_rhs += ABS(col->ub);
            ++cnt;
         }
         if( !SCIPsetIsEQ(set, SCIPvarGetUbGlobal(col->var), col->ub) )
         {
            if( SCIPvarGetType(col->var) == SCIP_VARTYPE_BINARY )
               obj = 1.0;
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



/** creates conflict analysis data for infeasible LP conflicts */
RETCODE SCIPlpconflictCreate(
   LPCONFLICT**     lpconflict          /**< pointer to LP conflict analysis data */
   )
{
   assert(lpconflict != NULL);

   ALLOC_OKAY( allocMemory(lpconflict) );

   CHECK_OKAY( SCIPclockCreate(&(*lpconflict)->analyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPlpiCreate(&(*lpconflict)->lpi, "LPconflict") );
   (*lpconflict)->ncalls = 0;
   (*lpconflict)->nconflicts = 0;
   (*lpconflict)->nlpiterations = 0;
   
   return SCIP_OKAY;
}

/** frees conflict analysis data for infeasible LP conflicts */
RETCODE SCIPlpconflictFree(
   LPCONFLICT**     lpconflict          /**< pointer to LP conflict analysis data */
   )
{
   assert(lpconflict != NULL);
   assert(*lpconflict != NULL);

   SCIPclockFree(&(*lpconflict)->analyzetime);
   CHECK_OKAY( SCIPlpiFree(&(*lpconflict)->lpi) );
   freeMemory(lpconflict);

   return SCIP_OKAY;
}

/** analyzes an infeasible LP trying to create a conflict set in the LP conflict analysis data structure */
static
RETCODE lpconflictAnalyze(
   LPCONFLICT*      lpconflict,         /**< LP conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< LP data */
   CONFLICT*        conflict,           /**< conflict analysis data */
   int              maxsize,            /**< maximal size of conflict set */
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
   int iterations;

   assert(lpconflict != NULL);
   assert(success != NULL);
  
   *success = FALSE;
  
   ncols = SCIPlpGetNCols(lp);
   nrows = SCIPlpGetNRows(lp);

   /* get temporary memory for alternative LP information */
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &alt_rowvars, 2*nrows + 1) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &alt_colvars, 2*ncols) );
   CHECK_OKAY( SCIPsetAllocBufferArray(set, &alt_islower, 2*ncols) );

   /* create alternative LP */
   CHECK_OKAY( SCIPlpiClear(lpconflict->lpi) );
   CHECK_OKAY( lpGenerateAltLP(lp, set, lpconflict->lpi, 
                  &alt_nrowvars, &alt_ncolvars, alt_rowvars, alt_colvars, alt_islower) );

   /* solve alternative LP with primal simplex, because the LP is primal feasible but not necessary dual feasible */
   CHECK_OKAY( SCIPlpiSolvePrimal(lpconflict->lpi) );

   /* update the LP iteration count */
   CHECK_OKAY( SCIPlpiGetIntpar(lpconflict->lpi, SCIP_LPPAR_LPITER, &iterations) );
   lpconflict->nlpiterations += iterations;

   /* check, if the LP is stable and we have a primal feasible solution */
   if( SCIPlpiIsStable(lpconflict->lpi)
      && (SCIPlpiIsOptimal(lpconflict->lpi) || SCIPlpiIsPrimalUnbounded(lpconflict->lpi)) )
   {
      COL* col;
      Real* psol;
      Bool valid;
      int i;
      int j;

      debugMessage("found IIS\n");
     
      /* get primal solution of alternative LP: the non-zero variables define an IIS in the original LP */
      CHECK_OKAY( SCIPsetAllocBufferArray(set, &psol, alt_nrowvars + alt_ncolvars) );
      CHECK_OKAY( SCIPlpiGetSol(lpconflict->lpi, NULL, psol, NULL, NULL, NULL) ); 

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
                     CHECK_OKAY( SCIPconflictAddVar(conflict, memhdr, set, stat, col->var) );
                     debugMessage(" -> adding <%s> to IIS\n", SCIPvarGetName(col->var));
                  }
                  else
                  {
                     /* we cannot put a non-binary variable into the conflict set: we have to abort */
                     debugMessage(" -> IIS includes changed bound of non-binary variable <%s>\n", SCIPvarGetName(col->var));
                     valid = FALSE;
                     break;
                  }
               }      
            }
         }

         if( valid )
         {
            /* analyze the conflict set, and create a conflict constraint on success */
            CHECK_OKAY( conflictAnalyze(conflict, set, maxsize, success) );
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

/** analyzes an infeasible LP to find out the bound changes on binary variables that were responsible for the infeasibility;
 *  on success, calls standard conflict analysis with the responsible variables as starting conflict set, thus creating
 *  a conflict constraint out of the resulting conflict set;
 *  updates statistics for infeasible LP conflict analysis
 */
RETCODE SCIPlpconflictAnalyze(
   LPCONFLICT*      lpconflict,         /**< LP conflict analysis data */
   MEMHDR*          memhdr,             /**< block memory of transformed problem */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< LP data */
   CONFLICT*        conflict,           /**< conflict analysis data */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   Bool valid;
   int maxsize;

   assert(lpconflict != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if infeasible LP conflict analysis is enabled */
   if( !set->uselpconflict )
      return SCIP_OKAY;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
      return SCIP_OKAY;

   /* if the root LP was infeasible, nothing has to be done */
   if( stat->nnodes == 1 )
      return SCIP_OKAY;

   /* calculate the maximal size of the conflict set */
   maxsize = (int)(set->maxconfvarsfac * prob->nbin);
   maxsize = MAX(maxsize, set->minmaxconfvars);
   if( maxsize < 2 )
      return SCIP_OKAY;

   /* start timing */
   SCIPclockStart(lpconflict->analyzetime, set);

   lpconflict->ncalls++;

   /* analyze conflict */
   CHECK_OKAY( lpconflictAnalyze(lpconflict, memhdr, set, stat, lp, conflict, maxsize, &valid) );

   /* check, if a valid conflict have been found */
   if( valid )
   {
      lpconflict->nconflicts++;
      if( success != NULL )
         *success = TRUE;
   }

   /* stop timing */
   SCIPclockStop(lpconflict->analyzetime, set);

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing infeasible LP conflicts */
Real SCIPlpconflictGetTime(
   LPCONFLICT*      lpconflict          /**< LP conflict analysis data */
   )
{
   assert(lpconflict != NULL);

   return SCIPclockGetTime(lpconflict->analyzetime);
}

/** gets number of calls to infeasible LP conflict analysis */
Longint SCIPlpconflictGetNCalls(
   LPCONFLICT*      lpconflict          /**< LP conflict analysis data */
   )
{
   assert(lpconflict != NULL);

   return lpconflict->ncalls;
}

/** gets number of valid conflicts detected in infeasible LP conflict analysis */
Longint SCIPlpconflictGetNConflicts(
   LPCONFLICT*      lpconflict          /**< LP conflict analysis data */
   )
{
   assert(lpconflict != NULL);

   return lpconflict->nconflicts;
}

/** gets number of LP iterations in infeasible LP conflict analysis */
Longint SCIPlpconflictGetNLPIterations(
   LPCONFLICT*      lpconflict          /**< LP conflict analysis data */
   )
{
   assert(lpconflict != NULL);

   return lpconflict->nlpiterations;
}
