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

/**@file   primal.c
 * @brief  datastructures and methods for collecting primal CIP solutions and primal informations
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "primal.h"
#include "disp.h"



/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that sols array can store at least num entries */
static
RETCODE ensureSolsSize(
   PRIMAL*          primal,             /**< primal data */
   const SET*       set,                /**< global SCIP settings */
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




/** sets upper bound in primal data and in LP solver */
static
RETCODE primalSetUpperbound(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   Real             upperbound          /**< new upper bound */
   )
{
   assert(primal != NULL);
   assert(lp != NULL);
   assert(upperbound <= set->infinity);

   debugMessage("set upper bound from %g to %g\n", primal->upperbound, upperbound);
   if( upperbound < primal->upperbound )
   {
      primal->upperbound = upperbound;
      CHECK_OKAY( SCIPlpSetUpperbound(lp, upperbound) );
      if( tree != NULL )
      {
         CHECK_OKAY( SCIPtreeCutoff(tree, memhdr, set, lp, upperbound) );
      }
   }
   else
   {
      errorMessage("Invalid increase in upper bound");
      return SCIP_INVALIDDATA;
   }

   return SCIP_OKAY;
}

/** creates primal data */
RETCODE SCIPprimalCreate(
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(primal != NULL);

   ALLOC_OKAY( allocMemory(primal) );
   (*primal)->sols = NULL;
   (*primal)->solssize = 0;
   (*primal)->nsols = 0;
   (*primal)->nsolsfound = 0;
   (*primal)->upperbound = SCIP_INVALID;

   CHECK_OKAY( primalSetUpperbound(*primal, memhdr, set, NULL, lp, MIN(prob->objlim, set->infinity)) );

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

   /* free primal CIP solutions */
   for( s = 0; s < (*primal)->nsols; ++s )
   {
      CHECK_OKAY( SCIPsolFree(&(*primal)->sols[s], memhdr) );
   }
   freeMemoryArrayNull(&(*primal)->sols);
   freeMemory(primal);

   return SCIP_OKAY;
}

/** adds primal solution to solution storage at given position */
static
RETCODE primalAddSol(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol,                /**< primal CIP solution */
   int              insertpos           /**< position in solution storage to add solution to */
   )
{
   EVENT event;
   int pos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(0 <= insertpos && insertpos < set->maxsol);

   debugMessage("insert primal solution at position %d:", insertpos);
   debug( SCIPsolPrint(sol, memhdr, set, stat, prob, NULL) );

   /* completely fill the solution's own value array to unlink it from the LP or pseudo solution */
   CHECK_OKAY( SCIPsolUnlink(sol, memhdr, set, prob) );

   /* allocate memory for solution storage */
   CHECK_OKAY( ensureSolsSize(primal, set, set->maxsol) );
   
   /* if the solution storage is full, free the last solution(s)
    * more than one solution may be freed, if set->maxsol was decreased in the meantime
    */
   for( pos = set->maxsol-1; pos < primal->nsols; ++pos )
   {
      CHECK_OKAY( SCIPsolFree(&primal->sols[pos], memhdr) );
   }

   /* insert solution at correct position */
   primal->nsols = MIN(primal->nsols+1, set->maxsol);
   for( pos = primal->nsols-1; pos > insertpos; --pos )
      primal->sols[pos] = primal->sols[pos-1];

   primal->sols[insertpos] = sol;
   primal->nsolsfound++;
   debugMessage(" -> stored at position %d of %d solutions, found %d solutions\n", 
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
   CHECK_OKAY( SCIPeventProcess(&event, memhdr, set, NULL, NULL, NULL, eventfilter) );

   /* check, if the global upper bound has to be updated */
   if( SCIPsolGetObj(sol) < primal->upperbound )
   {
      /* update the upper bound */
      CHECK_OKAY( primalSetUpperbound(primal, memhdr, set, tree, lp, SCIPsolGetObj(sol)) );
      
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

/** adds primal solution to solution storage by moving it */
RETCODE SCIPprimalAddSolMove(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL**            sol                 /**< pointer to primal CIP solution; is cleared in function call */
   )
{
   int insertpos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(*sol));

   if( insertpos < set->maxsol )
   {
      /* insert solution into solution storage */
      CHECK_OKAY( primalAddSol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, *sol, insertpos) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;
   }
   else
   {
      /* the solution is too bad -> free it immediately */
      CHECK_OKAY( SCIPsolFree(sol, memhdr) );
   }
   assert(*sol == NULL);

   return SCIP_OKAY;
}

/** adds primal solution to solution storage by copying it */
RETCODE SCIPprimalAddSolCopy(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol                 /**< primal CIP solution */
   )
{
   int insertpos;

   assert(primal != NULL);
   assert(sol != NULL);

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(sol));

   if( insertpos < set->maxsol )
   {
      SOL* solcopy;

      /* create a copy of the solution */
      CHECK_OKAY( SCIPsolCopy(&solcopy, memhdr, sol) );
      
      /* insert copied solution into solution storage */
      CHECK_OKAY( primalAddSol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, solcopy, insertpos) );
   }

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage by moving it */
RETCODE SCIPprimalTrySolMove(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL**            sol,                /**< pointer to primal CIP solution; is cleared in function call */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   Bool feasible;
   int insertpos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(stored != NULL);

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(*sol));

   if( insertpos < set->maxsol )
   {
      /* check solution for feasibility */
      CHECK_OKAY( SCIPsolCheck(*sol, memhdr, set, prob, chckintegrality, chcklprows, &feasible) );
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
      CHECK_OKAY( SCIPsolFree(sol, memhdr) );
      *stored = FALSE;
   }
   assert(*sol == NULL);

   return SCIP_OKAY;
}

/** checks primal solution; if feasible, adds it to storage by copying it */
RETCODE SCIPprimalTrySolCopy(
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SOL*             sol,                /**< primal CIP solution */
   Bool             chckintegrality,    /**< has integrality to be checked? */
   Bool             chcklprows,         /**< have current LP rows to be checked? */
   Bool*            stored              /**< stores whether given solution was feasible and good enough to keep */
   )
{
   Bool feasible;
   int insertpos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(stored != NULL);

   /* search the position to insert solution in storage */
   insertpos = primalSearchSolPos(primal, SCIPsolGetObj(sol));

   if( insertpos < set->maxsol )
   {
      /* check solution for feasibility */
      CHECK_OKAY( SCIPsolCheck(sol, memhdr, set, prob, chckintegrality, chcklprows, &feasible) );
   }
   else
      feasible = FALSE;

   if( feasible )
   {
      SOL* solcopy;

      /* create a copy of the solution */
      CHECK_OKAY( SCIPsolCopy(&solcopy, memhdr, sol) );
      
      /* insert copied solution into solution storage */
      CHECK_OKAY( primalAddSol(primal, memhdr, set, stat, prob, tree, lp, eventfilter, solcopy, insertpos) );

      *stored = TRUE;
   }
   else
      *stored = FALSE;

   return SCIP_OKAY;
}
