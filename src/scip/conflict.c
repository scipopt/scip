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
#pragma ident "@(#) $Id: conflict.c,v 1.17 2003/12/01 16:14:27 bzfpfend Exp $"

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
#include "clock.h"
#include "lpi.h"
#include "misc.h"
#include "paramset.h"
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

/** resizes conflictvars array to be able to store at least num constraints */
static
RETCODE lpconflictEnsureConflictvarsMem(
   LPCONFLICT*      lpconflict,         /**< LP conflict analysis data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(lpconflict != NULL);

   if( num > lpconflict->conflictvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&lpconflict->conflictvars, newsize) );
      lpconflict->conflictvarssize = newsize;
   }
   assert(num <= lpconflict->conflictvarssize);

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
   else if( SCIPvarGetInferNum(var1) > SCIPvarGetInferNum(var2) )
      return -1;
   else if( SCIPvarGetInferNum(var1) < SCIPvarGetInferNum(var2) )
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
      SCIPvarGetInferDepth(var), SCIPvarGetInferNum(var), SCIPvarGetInferCons(var));

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

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and creates a conflict set in the
 *  conflict analysis data structure
 */
static
RETCODE conflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   const SET*       set,                /**< global SCIP settings */
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

   debugMessage("analyzing conflict with %d conflict candidates\n", SCIPpqueueNElems(conflict->varqueue));

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
         printf("processing conflict var <%s>\n", SCIPvarGetName(var));
         printf(" - conflict set   :");
         for( v = 0; v < conflict->nconflictvars; ++v )
            printf(" <%s>", SCIPvarGetName(conflict->conflictvars[v]));
         printf("\n");
         printf(" - candidate queue:");
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
               SCIPvarGetInferDepth(var), SCIPvarGetInferNum(var));
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

   *success = (var == NULL);

   return SCIP_OKAY;
}

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set
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

   /* analyze conflict */
   CHECK_OKAY( conflictAnalyze(conflict, set, maxsize, &valid) );

   /* if a valid conflict set was found, call the conflict handlers */
   if( valid )
   {
      RESULT result;
      Bool resolved;
      int h;

      assert(conflict->conflictvars != NULL);
      assert(conflict->nconflictvars > 0);

      /* sort conflict handlers by priority */
      SCIPsetSortConflicthdlrs(set);

      /* call conflict handlers to create a conflict constraint */
      resolved = FALSE;
      for( h = 0; h < set->nconflicthdlrs; ++h )
      {
         CHECK_OKAY( SCIPconflicthdlrExec(set->conflicthdlrs[h], set->scip, 
                        conflict->conflictvars, conflict->nconflictvars, resolved, &result) );
         resolved = resolved || (result == SCIP_CONSADDED);
      }

      /* check, if a conflict constraint was created */
      if( resolved )
      {
         if( success != NULL )
            *success = TRUE;
         conflict->nconflicts++;
      }
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

/* ?????? MARC: Hier sollten Deine Methoden hinkommen; bitte halte Dich ein bisschen an das Sourcecode-Format,
 *              das Du oben siehst, insbesondere alle externen Methoden mit SCIP... beginnen und ins conflict.h
 *              file eintragen, alle internen Methoden static machen und den Namen mit nem Kleinbuchstaben
 *              anfangen lassen
 */

/** creates conflict analysis data for infeasible LP conflicts */
RETCODE SCIPlpconflictCreate(
   LPCONFLICT**     lpconflict          /**< pointer to LP conflict analysis data */
   )
{
   assert(lpconflict != NULL);

   ALLOC_OKAY( allocMemory(lpconflict) );

   CHECK_OKAY( SCIPclockCreate(&(*lpconflict)->analyzetime, SCIP_CLOCKTYPE_DEFAULT) );
   CHECK_OKAY( SCIPlpiCreate(&(*lpconflict)->lpi, "LPconflict") );
   (*lpconflict)->conflictvars = NULL;
   (*lpconflict)->conflictvarssize = 0;
   (*lpconflict)->nconflictvars = 0;
   (*lpconflict)->ncalls = 0;
   (*lpconflict)->nconflicts = 0;
   
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
   freeMemoryArrayNull(&(*lpconflict)->conflictvars);
   freeMemory(lpconflict);

   return SCIP_OKAY;
}

/** analyzes an infeasible LP trying to create a conflict set in the LP conflict analysis data structure */
static
RETCODE lpconflictAnalyze(
   LPCONFLICT*      lpconflict,         /**< LP conflict analysis data */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   int              maxsize,            /**< maximal size of the conflict set or -1 for no restriction */
   Bool*            success             /**< pointer to store whether the conflict set is valid */
   )
{
   assert(lpconflict != NULL);
   assert(success != NULL);

   *success = FALSE;

   /*lint --e{715}
    * ????? MARC: Hier kommen Deine wesentlichen Methoden!
    *             Diese sollten das LP analysieren (Zugriff ueber Methoden aus lp.h) und
    *             als Output sollten Sie das lpconflict->conflictvars array fuellen, und zwar
    *             mit der Bedeutung, dass das Setzen aller Variablen in dem Array auf FALSE zu
    *             einem unzulaessigen LP fuehrt (das Constraint "mindestens eine dieser Variablen
    *             muss TRUE sein" wird dann automatisch erzeugt).
    *             Die Methode der Wahl, um das Array in die richtige Groesse zu bringen, ist
    *               CHECK_OKAY( lpconflictEnsureConflictvarsMem(lpconflict, set, num) );
    *             wobei num die gewuenschte Mindestgroesse ist.
    *             Das LPI interface kennst Du ja schon ein bisschen. Fuer den Anfang wuerde ich
    *             die Analyze erst mal mit nem 
    *               CHECK_OKAY( SCIPlpiClear(lpi) );
    *             beginnen. Wie man das macht, dass das LP nur upgedated werden muss, koennen wir
    *             uns ja spaeter noch ueberlegen.
    */
   /*printf("LP infeasible: analyze LP conflict\n");*/ /*??????????????????????*/

   return SCIP_OKAY;
}

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and on success, calls the
 *  conflict handlers to create a conflict constraint out of the resulting conflict set
 */
RETCODE SCIPlpconflictAnalyze(
   LPCONFLICT*      lpconflict,         /**< LP conflict analysis data */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   LP*              lp,                 /**< LP data */
   Bool*            success             /**< pointer to store whether a conflict constraint was created, or NULL */
   )
{
   Bool valid;
   int maxsize;

   assert(lpconflict != NULL);
   assert(set != NULL);

   if( success != NULL )
      *success = FALSE;

   /* check, if there are any conflict handlers to use a conflict set */
   if( set->nconflicthdlrs == 0 )
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
   CHECK_OKAY( lpconflictAnalyze(lpconflict, set, lp, maxsize, &valid) );

   /* if a valid conflict set was found, call the conflict handlers */
   if( valid )
   {
      RESULT result;
      Bool resolved;
      int h;

      assert(lpconflict->conflictvars != NULL);
      assert(lpconflict->nconflictvars > 0);

      /* call conflict handlers to create a conflict constraint */
      resolved = FALSE;
      for( h = 0; h < set->nconflicthdlrs; ++h )
      {
         CHECK_OKAY( SCIPconflicthdlrExec(set->conflicthdlrs[h], set->scip, 
                        lpconflict->conflictvars, lpconflict->nconflictvars, resolved, &result) );
         resolved = resolved || (result == SCIP_CONSADDED);
      }

      /* check, if a conflict constraint was created */
      if( resolved )
      {
         if( success != NULL )
            *success = TRUE;
         lpconflict->nconflicts++;
      }
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
