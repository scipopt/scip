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
#pragma ident "@(#) $Id: primal.c,v 1.50 2004/11/26 14:22:12 bzfpfend Exp $"

/**@file   primal.c
 * @brief  methods for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "stat.h"
#include "vbc.h"
#include "event.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "sol.h"
#include "primal.h"
#include "tree.h"
#include "disp.h"



/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that sols array can store at least num entries */
static
RETCODE ensureSolsSize(
   PRIMAL*          primal,             /**< primal data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(primal->nsols <= primal->solssize);
   
   if( num > primal->solssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&primal->sols, newsize) );
      primal->solssize = newsize;
   }
   assert(num <= primal->solssize);

   return SCIP_OKAY;
}

/** ensures, that existingsols array can store at least num entries */
static
RETCODE ensureExistingsolsSize(
   PRIMAL*          primal,             /**< primal data */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(primal->nexistingsols <= primal->existingsolssize);
   
   if( num > primal->existingsolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&primal->existingsols, newsize) );
      primal->existingsolssize = newsize;
   }
   assert(num <= primal->existingsolssize);

   return SCIP_OKAY;
}




/** creates primal data */
RETCODE SCIPprimalCreate(
   PRIMAL**         primal              /**< pointer to primal data */
   )
{
   assert(primal != NULL);

   ALLOC_OKAY( allocMemory(primal) );
   (*primal)->sols = NULL;
   (*primal)->existingsols = NULL;
   (*primal)->currentsol = NULL;
   (*primal)->solssize = 0;
   (*primal)->nsols = 0;
   (*primal)->existingsolssize = 0;
   (*primal)->nexistingsols = 0;
   (*primal)->nsolsfound = 0;
   (*primal)->nbestsolsfound = 0;
   (*primal)->upperbound = SCIP_INVALID;
   (*primal)->cutoffbound = SCIP_INVALID;

   return SCIP_OKAY;
}

/** frees primal data */
RETCODE SCIPprimalFree(
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr              /**< block memory */
   )
{
   int s;

   assert(primal != NULL);
   assert(*primal != NULL);

   /* free temporary solution for storing current solution */
   if( (*primal)->currentsol != NULL )
   {
      CHECK_OKAY( SCIPsolFree(&(*primal)->currentsol, memhdr, *primal) );
   }

   /* free feasible primal CIP solutions */
   for( s = 0; s < (*primal)->nsols; ++s )
   {
      CHECK_OKAY( SCIPsolFree(&(*primal)->sols[s], memhdr, *primal) );
   }
   assert((*primal)->nexistingsols == 0);

   freeMemoryArrayNull(&(*primal)->sols);
   freeMemoryArrayNull(&(*primal)->existingsols);
   freeMemory(primal);

   return SCIP_OKAY;
}

/**@todo implement an objective function domain propagator instead of the following code */
#if 0
/** applies domain propagation of the objective function after a new upper bound was found */
static
RETCODE primalPropagateGlobalObj(
   PRIMAL*          primal,             /**< primal data */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree                /**< branch and bound tree */
   )
{
   VAR* var;
   Real minobj;
   Real maxobj;
   Real obj;
   Real maxdelta;
   Real lb;
   Real ub;
   Real upperbound;
   Real lowerbound;
   int v;

   assert(primal != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   /* if there are unknown variables, we cannot propagate on the objective function */
   if( set->nactivepricers > 0 )
      return SCIP_OKAY;

   /* calculate minimal and maximal objective value in current global bounds */
   minobj = 0.0;
   maxobj = 0.0;
   for( v = 0; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      obj = SCIPvarGetObj(var);
      if( obj > 0.0 )
      {
         minobj += SCIPvarGetLbGlobal(var) * obj;
         maxobj += SCIPvarGetUbGlobal(var) * obj;
      }
      else if( obj < 0.0 )
      {
         minobj += SCIPvarGetUbGlobal(var) * obj;
         maxobj += SCIPvarGetLbGlobal(var) * obj;
      }
   }

   /* update global bounds w.r.t. the objective function */
   lowerbound = SCIPtreeGetLowerbound(tree, set);
   upperbound = primal->cutoffbound;
   for( v = 0; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      obj = SCIPvarGetObj(var);
      if( SCIPsetIsPositive(set, obj) )
      {
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);

         maxdelta = (lowerbound - maxobj) / obj;
         if( SCIPsetIsLbBetter(set, ub + maxdelta, lb) )
         {
            debugMessage("tightening global lower bound of <%s> due to dual bound: [%g,%g] -> [%g,%g]\n",
               SCIPvarGetName(var), lb, ub, ub+maxdelta, ub);
            CHECK_OKAY( SCIPvarSetLbGlobal(var, set, ub+maxdelta) );
            lb = ub + maxdelta;
            stat->nrootboundchgs++;
            stat->nrootboundchgsrun++;
         }

         maxdelta = (upperbound - minobj) / obj;
         if( SCIPsetIsUbBetter(set, lb + maxdelta, ub) )
         {
            debugMessage("tightening global upper bound of <%s> due to primal bound: [%g,%g] -> [%g,%g]\n",
               SCIPvarGetName(var), lb, ub, lb, lb+maxdelta);
            CHECK_OKAY( SCIPvarSetUbGlobal(var, set, lb+maxdelta) );
            ub = lb + maxdelta;
            stat->nrootboundchgs++;
            stat->nrootboundchgsrun++;
         }
      }
      else if( SCIPsetIsNegative(set, obj) )
      {
         lb = SCIPvarGetLbGlobal(var);
         ub = SCIPvarGetUbGlobal(var);

         maxdelta = (lowerbound - maxobj) / obj;
         if( SCIPsetIsUbBetter(set, lb + maxdelta, ub) )
         {
            debugMessage("tightening global upper bound of <%s> due to dual bound: [%g,%g] -> [%g,%g]\n",
               SCIPvarGetName(var), lb, ub, lb, lb+maxdelta);
            CHECK_OKAY( SCIPvarSetUbGlobal(var, set, lb+maxdelta) );
            ub = lb + maxdelta;
            stat->nrootboundchgs++;
            stat->nrootboundchgsrun++;
         }

         maxdelta = (upperbound - minobj) / obj;
         if( SCIPsetIsLbBetter(set, ub + maxdelta, lb) )
         {
            debugMessage("tightening global lower bound of <%s> due to primal bound: [%g,%g] -> [%g,%g]\n",
               SCIPvarGetName(var), lb, ub, ub+maxdelta, ub);
            CHECK_OKAY( SCIPvarSetLbGlobal(var, set, ub+maxdelta) );
            lb = ub + maxdelta;
            stat->nrootboundchgs++;
            stat->nrootboundchgsrun++;
         }
      }
   }

   return SCIP_OKAY;
}
#endif

/** sets upper bound in primal data and in LP solver */
static
RETCODE primalSetUpperbound(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   Real             upperbound          /**< new upper bound */
   )
{
   assert(primal != NULL);
   assert(set != NULL);
   assert(stat != NULL);
   assert(upperbound <= SCIPsetInfinity(set));
   assert(upperbound <= primal->upperbound || tree == NULL);

   debugMessage("changing upper bound from %g to %g\n", primal->upperbound, upperbound);

   primal->upperbound = upperbound;
   
   /* if objective value is always integral, the cutoff bound can be reduced to nearly the previous integer number */
   if( SCIPprobIsObjIntegral(prob) )
   {
      primal->cutoffbound = SCIPsetFeasCeil(set, upperbound) - (1.0 - 10.0*SCIPsetFeastol(set));
   }
   else
      primal->cutoffbound = upperbound;

   /* set cut off value in LP solver */
   CHECK_OKAY( SCIPlpSetCutoffbound(lp, set, primal->cutoffbound) );

   if( tree != NULL )
   {
      /* cut off leaves of the tree */
      CHECK_OKAY( SCIPtreeCutoff(tree, memhdr, set, stat, lp, primal->cutoffbound) );

      /* update upper bound in VBC output */
      SCIPvbcUpperbound(stat->vbc, stat, primal->upperbound);

#if 0
      /* propagate global bounds on variables */
      CHECK_OKAY( primalPropagateGlobalObj(primal, set, stat, prob, tree) );
#endif
   }

   return SCIP_OKAY;
}

/** sets upper bound in primal data and in LP solver */
RETCODE SCIPprimalSetUpperbound(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   Real             upperbound          /**< new upper bound */
   )
{
   assert(primal != NULL);
   assert(upperbound <= SCIPsetInfinity(set));

   if( upperbound < primal->upperbound )
   {
      /* update primal bound */
      CHECK_OKAY( primalSetUpperbound(primal, memhdr, set, stat, prob, tree, lp, upperbound) );
   }
   else if( upperbound > primal->upperbound )
   {
      errorMessage("invalid increase in upper bound\n");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** recalculates upper bound in primal data after a change of the problem's objective offset */
RETCODE SCIPprimalUpdateUpperbound(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp                  /**< current LP data */
   )
{
   SOL* sol;
   Real upperbound;
   Real objval;
   int i;
   int j;

   assert(set != NULL);

   /* recalculate internal objective limit */
   upperbound = SCIPprobGetInternObjlim(prob, set);
   upperbound = MIN(upperbound, SCIPsetInfinity(set));

   /* resort current primal solutions */
   for( i = 1; i < primal->nsols; ++i )
   {
      sol = primal->sols[i];
      objval = SCIPsolGetObj(sol);
      for( j = i; j > 0 && objval < SCIPsolGetObj(primal->sols[j-1]); --j )
         primal->sols[j] = primal->sols[j-1];
      primal->sols[j] = sol;
   }

   /* delete all solutions worse than the current objective limit */
   for( ; primal->nsols > 0 && SCIPsolGetObj(primal->sols[primal->nsols-1]) > upperbound; primal->nsols-- )
   {
      CHECK_OKAY( SCIPsolFree(&primal->sols[primal->nsols-1], memhdr, primal) );
   }

   /* compare objective limit to currently best solution */
   if( primal->nsols > 0 )
      upperbound = MIN(upperbound, SCIPsolGetObj(primal->sols[0]));

   /* set new upper bound */
   if( upperbound != primal->upperbound )
   {
      CHECK_OKAY( primalSetUpperbound(primal, memhdr, set, stat, prob, tree, lp, upperbound) );
   }

   return SCIP_OKAY;
}

/** returns whether the current primal bound is justified with a feasible primal solution; if not, the primal bound
 *  was set from the user as objective limit
 */
Bool SCIPprimalUpperboundIsSol(
   PRIMAL*          primal              /**< primal data */
   )
{
   assert(primal != NULL);

   return (primal->nsols > 0 && primal->upperbound == SCIPsolGetObj(primal->sols[0]));
}

/** adds primal solution to solution storage at given position */
static
RETCODE primalAddSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol,                /**< primal CIP solution */
   int              insertpos           /**< position in solution storage to add solution to */
   )
{
   EVENT event;
   Real obj;
   int pos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(0 <= insertpos && insertpos < set->limit_maxsol);

   debugMessage("insert primal solution at position %d:", insertpos);
   debug( SCIPsolPrint(sol, set, stat, prob, NULL, NULL) );

#if 0
#ifndef NDEBUG
   /* check solution again completely
    * (this may fail, because in the LP solver, the feasibility tolerance is a relative measure against the row's norm
    */
   {
      Bool feasible;
      CHECK_OKAY( SCIPsolCheck(sol, memhdr, set, stat, prob, TRUE, TRUE, &feasible) );
      if( !feasible )
      {
         errorMessage("infeasible solution accepted:\n");
         CHECK_OKAY( SCIPsolPrint(sol, set, stat, prob, NULL, NULL) );
      }
      assert(feasible);
   }
#endif
#endif

   /* completely fill the solution's own value array to unlink it from the LP or pseudo solution */
   CHECK_OKAY( SCIPsolUnlink(sol, set, prob) );

   /* allocate memory for solution storage */
   CHECK_OKAY( ensureSolsSize(primal, set, set->limit_maxsol) );
   
   /* if the solution storage is full, free the last solution(s)
    * more than one solution may be freed, if set->limit_maxsol was decreased in the meantime
    */
   for( pos = set->limit_maxsol-1; pos < primal->nsols; ++pos )
   {
      CHECK_OKAY( SCIPsolFree(&primal->sols[pos], memhdr, primal) );
   }

   /* insert solution at correct position */
   primal->nsols = MIN(primal->nsols+1, set->limit_maxsol);
   for( pos = primal->nsols-1; pos > insertpos; --pos )
      primal->sols[pos] = primal->sols[pos-1];

   assert(0 <= insertpos && insertpos < primal->nsols);
   primal->sols[insertpos] = sol;
   primal->nsolsfound++;
   debugMessage(" -> stored at position %d of %d solutions, found %lld solutions\n", 
      insertpos, primal->nsols, primal->nsolsfound);

   /* update the solution value sums in variables */
   SCIPsolUpdateVarsum(sol, set, stat, prob, (Real)(primal->nsols - insertpos)/(Real)(2.0*primal->nsols - 1.0));

   /* change color of node in VBC output */
   SCIPvbcFoundSolution(stat->vbc, stat, SCIPtreeGetCurrentNode(tree));

   /* check, if the global upper bound has to be updated */
   obj = SCIPsolGetObj(sol);
   if( obj < primal->upperbound )
   {
      /* update the upper bound */
      CHECK_OKAY( SCIPprimalSetUpperbound(primal, memhdr, set, stat, prob, tree, lp, obj) );
      
      /* display node information line */
      CHECK_OKAY( SCIPdispPrintLine(set, stat, TRUE) );

      /* issue BESTLPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_BESTSOLFOUND) );
      primal->nbestsolsfound++;
   }
   else
   {
      /* issue POORLPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_POORSOLFOUND) );
   }
   CHECK_OKAY( SCIPeventChgSol(&event, sol) );
   CHECK_OKAY( SCIPeventProcess(&event, memhdr, set, NULL, NULL, NULL, eventfilter) );

   return SCIP_OKAY;
}

/** uses binary search to find position in solution storage */
static
int primalSearchSolPos(
   PRIMAL*          primal,             /**< primal data */
   Real             obj                 /**< objective value of solution to search position for */
   )
{
   Real middleobj;
   int left;
   int right;
   int middle;

   assert(primal != NULL);

   left = -1;
   right = primal->nsols;
   while( left < right-1 )
   {
      middle = (left+right)/2;
      assert(left < middle && middle < right);
      assert(0 <= middle && middle < primal->nsols);
      middleobj = SCIPsolGetObj(primal->sols[middle]);
      if( obj < middleobj )
         right = middle;
      else
         left = middle;
   }
   assert(left == right-1);

   return right;
}

/** adds primal solution to solution storage by copying it */
RETCODE SCIPprimalAddSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol,                /**< primal CIP solution */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   int insertpos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(stored != NULL);

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(sol));

   if( insertpos < set->limit_maxsol )
   {
      SOL* solcopy;

      /* create a copy of the solution */
      CHECK_OKAY( SCIPsolCopy(&solcopy, memhdr, set, primal, sol) );
      
      /* insert copied solution into solution storage */
      CHECK_OKAY( primalAddSol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, solcopy, insertpos) );

      *stored = TRUE;
   }
   else
      *stored = FALSE;

   return SCIP_OKAY;
}

/** adds primal solution to solution storage, frees the solution afterwards */
RETCODE SCIPprimalAddSolFree(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   int insertpos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(stored != NULL);

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(*sol));

   if( insertpos < set->limit_maxsol )
   {
      /* insert solution into solution storage */
      CHECK_OKAY( primalAddSol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, *sol, insertpos) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;

      *stored = TRUE;
   }
   else
   {
      /* the solution is too bad -> free it immediately */
      CHECK_OKAY( SCIPsolFree(sol, memhdr, primal) );

      *stored = FALSE;
   }
   assert(*sol == NULL);

   return SCIP_OKAY;
}

/** links temporary solution of primal data to current solution */
static
RETCODE primalLinkCurrentSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's from the tree) */
   )
{
   assert(primal != NULL);

   if( primal->currentsol == NULL )
   {
      CHECK_OKAY( SCIPsolCreateCurrentSol(&primal->currentsol, memhdr, set, stat, primal, tree, lp, heur) );
   }
   else
   {
      CHECK_OKAY( SCIPsolLinkCurrentSol(primal->currentsol, memhdr, set, stat, tree, lp) );
      SCIPsolSetHeur(primal->currentsol, heur);
   }

   return SCIP_OKAY;
}

/** adds current LP/pseudo solution to solution storage */
RETCODE SCIPprimalAddCurrentSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   HEUR*            heur,               /**< heuristic that found the solution (or NULL if it's from the tree) */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   assert(primal != NULL);

   /* link temporary solution to current solution */
   CHECK_OKAY( primalLinkCurrentSol(primal, memhdr, set, stat, tree, lp, heur) );

   /* add solution to solution storage */
   CHECK_OKAY( SCIPprimalAddSol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, primal->currentsol, stored) );

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage by copying it */
RETCODE SCIPprimalTrySol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol,                /**< primal CIP solution */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   Bool feasible;
   int insertpos;

   assert(primal != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(sol != NULL);
   assert(stored != NULL);

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || set->misc_exactsolve;

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(sol));

   if( insertpos < set->limit_maxsol )
   {
      /* check solution for feasibility */
      CHECK_OKAY( SCIPsolCheck(sol, memhdr, set, stat, prob, checkintegrality, checklprows, &feasible) );
   }
   else
      feasible = FALSE;

   if( feasible )
   {
      SOL* solcopy;
      
      /* create a copy of the solution */
      CHECK_OKAY( SCIPsolCopy(&solcopy, memhdr, set, primal, sol) );
      
      /* insert copied solution into solution storage */
      CHECK_OKAY( primalAddSol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, solcopy, insertpos) );
      
      *stored = TRUE;
   }
   else
      *stored = FALSE;

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage; solution is freed afterwards */
RETCODE SCIPprimalTrySolFree(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether solution was feasible and good enough to keep */
   )
{
   Bool feasible;
   int insertpos;

   assert(primal != NULL);
   assert(tree != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(stored != NULL);

   *stored = FALSE;

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || set->misc_exactsolve;

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(*sol));

   if( insertpos < set->limit_maxsol )
   {
      /* check solution for feasibility */
      CHECK_OKAY( SCIPsolCheck(*sol, memhdr, set, stat, prob, checkintegrality, checklprows, &feasible) );
   }
   else
      feasible = FALSE;

   if( feasible )
   {
      /* insert solution into solution storage */
      CHECK_OKAY( primalAddSol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, *sol, insertpos) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;
      *stored = TRUE;
   }
   else
   {
      /* the solution is too bad or infeasible -> free it immediately */
      CHECK_OKAY( SCIPsolFree(sol, memhdr, primal) );
      *stored = FALSE;
   }
   assert(*sol == NULL);

   return SCIP_OKAY;
}

/** checks current LP/pseudo solution; if feasible, adds it to storage */
RETCODE SCIPprimalTryCurrentSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   HEUR*            heur,               /**< heuristic that found the solution (or NULL if it's from the tree) */
   Bool             checkintegrality,   /**< has integrality to be checked? */
   Bool             checklprows,        /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   assert(primal != NULL);

   /* link temporary solution to current solution */
   CHECK_OKAY( primalLinkCurrentSol(primal, memhdr, set, stat, tree, lp, heur) );

   /* add solution to solution storage */
   CHECK_OKAY( SCIPprimalTrySol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, primal->currentsol,
         checkintegrality, checklprows, stored) );

   return SCIP_OKAY;
}

/** inserts solution into the global array of all existing primal solutions */
RETCODE SCIPprimalSolCreated(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   SET*             set,                /**< global SCIP settings */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   assert(primal != NULL);
   assert(sol != NULL);
   assert(SCIPsolGetPrimalIndex(sol) == -1);

   /* allocate memory for solution storage */
   CHECK_OKAY( ensureExistingsolsSize(primal, set, primal->nexistingsols+1) );

   /* append solution */
   SCIPsolSetPrimalIndex(sol, primal->nexistingsols);
   primal->existingsols[primal->nexistingsols] = sol;
   primal->nexistingsols++;

   return SCIP_OKAY;
}

/** removes solution from the global array of all existing primal solutions */
void SCIPprimalSolFreed(
   PRIMAL*          primal,             /**< primal data */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   int idx;

   assert(primal != NULL);
   assert(sol != NULL);

   /* remove solution */
   idx = SCIPsolGetPrimalIndex(sol);
   assert(0 <= idx && idx < primal->nexistingsols);
   if( idx < primal->nexistingsols-1 )
   {
      primal->existingsols[idx] = primal->existingsols[primal->nexistingsols-1];
      SCIPsolSetPrimalIndex(primal->existingsols[idx], idx);
   }
   primal->nexistingsols--;
}

/** updates all existing primal solutions after a change in a variable's objective value */
void SCIPprimalUpdateVarObj(
   PRIMAL*          primal,             /**< primal data */
   VAR*             var,                /**< problem variable */
   Real             oldobj,             /**< old objective value */
   Real             newobj              /**< new objective value */
   )
{
   int i;

   assert(primal != NULL);

   for( i = 0; i < primal->nexistingsols; ++i )
   {
      SCIPsolUpdateVarObj(primal->existingsols[i], var, oldobj, newobj);
   }
}
