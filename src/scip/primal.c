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
#pragma ident "@(#) $Id: primal.c,v 1.37 2004/04/30 11:16:25 bzfpfend Exp $"

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
   (*primal)->solssize = 0;
   (*primal)->nsols = 0;
   (*primal)->existingsolssize = 0;
   (*primal)->nexistingsols = 0;
   (*primal)->nsolsfound = 0;
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
   assert(upperbound <= set->infinity);
   assert(upperbound <= primal->upperbound || tree == NULL);

   debugMessage("changing upper bound from %g to %g\n", primal->upperbound, upperbound);

   primal->upperbound = upperbound;
   
   /* if objective value is always integral, the cutoff bound can be reduced to nearly the previous integer number */
   if( SCIPprobIsObjIntegral(prob) )
   {
      primal->cutoffbound = SCIPsetCeil(set, upperbound) - (1.0 - 10.0*set->feastol);
   }
   else
      primal->cutoffbound = upperbound;

   /* set cut off value in LP solver */
   CHECK_OKAY( SCIPlpSetCutoffbound(lp, set, primal->cutoffbound) );

   if( tree != NULL )
   {
      /* cut off leaves of the tree */
      CHECK_OKAY( SCIPtreeCutoff(tree, memhdr, set, lp, primal->cutoffbound) );

      /* update upper bound in VBC output */
      SCIPvbcUpperbound(stat->vbc, stat, primal->upperbound);
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
   assert(upperbound <= set->infinity);

   if( upperbound < primal->upperbound )
   {
      /* update primal bound */
      CHECK_OKAY( primalSetUpperbound(primal, memhdr, set, stat, prob, tree, lp, upperbound) );
   }
   else if( upperbound > primal->upperbound )
   {
      errorMessage("Invalid increase in upper bound\n");
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
   upperbound = MIN(upperbound, set->infinity);

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
   assert(0 <= insertpos && insertpos < set->maxsol);

   debugMessage("insert primal solution at position %d:", insertpos);
   debug( SCIPsolPrint(sol, set, stat, prob, NULL) );

#if 1
#ifndef NDEBUG
   /* check solution again completely
    * (this may fail, because in the LP solver, the feasibility tolerance is a relative measure against the row's norm
    */
   {
      Bool feasible;
      CHECK_OKAY( SCIPsolCheck(sol, memhdr, set, prob, TRUE, TRUE, &feasible) );
      if( !feasible )
      {
         errorMessage("infeasible solution accepted:\n");
         CHECK_OKAY( SCIPsolPrint(sol, set, stat, prob, NULL) );
      }
      assert(feasible);
   }
#endif
#endif

   /* completely fill the solution's own value array to unlink it from the LP or pseudo solution */
   CHECK_OKAY( SCIPsolUnlink(sol, set, prob) );

   /* allocate memory for solution storage */
   CHECK_OKAY( ensureSolsSize(primal, set, set->maxsol) );
   
   /* if the solution storage is full, free the last solution(s)
    * more than one solution may be freed, if set->maxsol was decreased in the meantime
    */
   for( pos = set->maxsol-1; pos < primal->nsols; ++pos )
   {
      CHECK_OKAY( SCIPsolFree(&primal->sols[pos], memhdr, primal) );
   }

   /* insert solution at correct position */
   primal->nsols = MIN(primal->nsols+1, set->maxsol);
   for( pos = primal->nsols-1; pos > insertpos; --pos )
      primal->sols[pos] = primal->sols[pos-1];

   assert(0 <= insertpos && insertpos < primal->nsols);
   primal->sols[insertpos] = sol;
   primal->nsolsfound++;
   debugMessage(" -> stored at position %d of %d solutions, found %lld solutions\n", 
      insertpos, primal->nsols, primal->nsolsfound);
   
   /* issue POORLPSOLVED or BESTLPSOLVED event */
   if( insertpos == 0 )
   {
      /* issue BESTLPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_BESTSOLFOUND) );
   }
   else
   {
      /* issue POORLPSOLVED event */
      CHECK_OKAY( SCIPeventChgType(&event, SCIP_EVENTTYPE_POORSOLFOUND) );
   }
   CHECK_OKAY( SCIPeventChgSol(&event, sol) );
   CHECK_OKAY( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

   /* change color of node in VBC output */
   SCIPvbcFoundSolution(stat->vbc, stat, tree->actnode);

   /* check, if the global upper bound has to be updated */
   obj = SCIPsolGetObj(sol);
   if( obj < primal->upperbound )
   {
      /* update the upper bound */
      CHECK_OKAY( SCIPprimalSetUpperbound(primal, memhdr, set, stat, prob, tree, lp, obj) );
      
      /* display node information line */
      CHECK_OKAY( SCIPdispPrintLine(set, stat, TRUE) );
   }

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

   if( insertpos < set->maxsol )
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

   if( insertpos < set->maxsol )
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
   assert(sol != NULL);
   assert(stored != NULL);

   /* if we want to solve exactly, the constraint handlers cannot rely on the LP's feasibility */
   checklprows = checklprows || set->exactsolve;

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(sol));

   if( insertpos < set->maxsol )
   {
      /* check solution for feasibility */
      CHECK_OKAY( SCIPsolCheck(sol, memhdr, set, prob, checkintegrality, checklprows, &feasible) );
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
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(stored != NULL);

   *stored = FALSE;

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(*sol));

   if( insertpos < set->maxsol )
   {
      /* check solution for feasibility */
      CHECK_OKAY( SCIPsolCheck(*sol, memhdr, set, prob, checkintegrality, checklprows, &feasible) );
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
