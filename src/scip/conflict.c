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
 *     is at the top of the queue. Create an empty conflict set.
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
 *     (b) Otherwise, the assignment to variable v was a branching decision or
 *         a deduction with missing inference reason.
 *          - Put v in the conflict set.
 *            Note that if v was a branching variable, all deduced variables remaining
 *            in the priority queue have smaller inference depth level than v, since
 *            deductions are allways applied after the branching decisions. However,
 *            there is the possibility, that v was a deduced variable, where the
 *            inference reason was not given. With this lack of information, we must
 *            treat the deduced variable like a branching variable, and there may
 *            exist other deduced variables of the same inference depth level in
 *            the priority queue.
 *  4. If priority queue is non-empty, goto step 2.
 *  5. The conflict set represents the conflict clause saying that at least one
 *     of the conflict variables must be set to TRUE. The caller of the conflict
 *     analysis can form a corresponding constraint (e.g. a logicor constraint)
 *     out of these conflict variables and add it to the problem.
 *
 * If all deduced variables come with inference information, the resulting
 * conflict set has the property, that for each depth level at most one variable
 * assigned at that level is member of the conflict set. This conflict variable
 * is the first unique implication point of its depth level (1UIP).
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
 *  - If an (with the current assignment) infeasible constraint is detected,
 *    the constraint handler should
 *     1. call SCIPaddConflictVar() for each variable in the infeasible
 *        constraint,
 *     2. call SCIPanalyzeConflict() to get a set of variables that lead
 *        to the conflicing assignment when all set to FALSE,
 *     3. use the conflict set to create an appropriate constraint.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <math.h>

#include "misc.h"
#include "message.h"
#include "conflict.h"



struct Conflict
{
   CLOCK*           analyzetime;        /**< time used for conflict analysis */
   PQUEUE*          varqueue;           /**< unprocessed conflict variables */
   VAR**            conflictvars;       /**< variables resembling the conflict clause */
   int              conflictvarssize;   /**< size of conflictvars array */
   int              nconflictvars;      /**< number of variables in the conflict set (used slots of conflictvars array) */
   Longint          ncalls;             /**< number of calls to conflict analysis */
   Longint          nconflicts;         /**< number of valid conflicts detected in conflict analysis */
};




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




/** compares two binary variables w.r.t. their inference information, such that variables infered later are
 *  ordered prior to variables infered earlier
 */
static
DECL_SORTPTRCOMP(conflictVarCmp)
{
   VAR* var1;
   VAR* var2;
   
   var1 = (VAR*)elem1;
   var2 = (VAR*)elem2;
   assert(var1 != NULL);
   assert(var1->vartype == SCIP_VARTYPE_BINARY);
   assert(var2 != NULL);
   assert(var2->vartype == SCIP_VARTYPE_BINARY);

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

/** creates conflict analysis data */
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

/** frees conflict analysis data */
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

/** initializes the conflict analysis by clearing the conflict variable candidate queue */
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
   assert(var->vartype == SCIP_VARTYPE_BINARY);
   assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

   debugMessage("adding variable <%s>[%g,%g] [status:%d, depth:%d, num:%d, cons:%p] to conflict candidates\n",
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), SCIPvarGetStatus(var),
      SCIPvarGetInferDepth(var), SCIPvarGetInferNum(var), SCIPvarGetInferCons(var));

   /* get active problem variable */
   var = SCIPvarGetProbvar(var);

   /* we can ignore fixed variables */
   if( var != NULL )
   {
      assert(var->vartype == SCIP_VARTYPE_BINARY);
      assert(SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
      
      /* choose between variable or its negation, such that the literal is fixed to FALSE */
      if( SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), 1.0) )
      {
         /* variable is fixed to TRUE -> use the negation */
         CHECK_OKAY( SCIPvarGetNegated(var, memhdr, set, stat, &var) );
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
   int              maxsize,            /**< maximal size of the conflict set or -1 for no restriction */
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

   /* check, if there is something to analyze */
   if( SCIPpqueueNElems(conflict->varqueue) == 0 )
   {
      errorMessage("no conflict variables to analyze");
      return SCIP_INVALIDDATA;
   }

   if( maxsize == -1 )
      maxsize = INT_MAX;

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

/** analyzes conflict variables that were added with calls to SCIPconflictAddVar(), and returns a conflict set, that
 *  can be used to create a conflict constraint
 */
RETCODE SCIPconflictAnalyze(
   CONFLICT*        conflict,           /**< conflict analysis data */
   const SET*       set,                /**< global SCIP settings */
   int              maxsize,            /**< maximal size of the conflict set or -1 for no restriction */
   VAR***           conflictvars,       /**< pointer to store the conflict set (user must not change this array) */
   int*             nconflictvars,      /**< pointer to store the number of conflict variables */
   Bool*            success             /**< pointer to store whether the conflict set is valid */
   )
{
   assert(conflict != NULL);
   assert(conflictvars != NULL);
   assert(nconflictvars != NULL);

   conflict->ncalls++;

   /* start timing */
   SCIPclockStart(conflict->analyzetime, set);

   /* analyze conflict */
   CHECK_OKAY( conflictAnalyze(conflict, set, maxsize, success) );
   assert(conflict->conflictvars != NULL);
   assert(conflict->nconflictvars > 0);

   /* stop timing */
   SCIPclockStop(conflict->analyzetime, set);

   *conflictvars = conflict->conflictvars;
   *nconflictvars = conflict->nconflictvars;

   if( *success )
      conflict->nconflicts++;

   return SCIP_OKAY;
}

/** gets time in seconds used for analyzing conflicts */
Real SCIPconflictGetTime(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return SCIPclockGetTime(conflict->analyzetime);
}

/** gets number of calls to conflict analysis */
Longint SCIPconflictGetNCalls(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->ncalls;
}

/** gets number of valid conflicts detected in conflict analysis */
Longint SCIPconflictGetNConflicts(
   CONFLICT*        conflict            /**< conflict analysis data */
   )
{
   assert(conflict != NULL);

   return conflict->nconflicts;
}

