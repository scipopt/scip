/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
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

static
RETCODE ensureSolsSize(                 /**< ensures, that sols array can store at least num entries */
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
      ALLOC_OKAY( reallocMemoryArray(primal->sols, newsize) );
      primal->solssize = newsize;
   }
   assert(num <= primal->solssize);

   return SCIP_OKAY;
}




static
RETCODE primalSetUpperbound(            /**< sets upper bound in primal data and in LP solver */
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

RETCODE SCIPprimalCreate(               /**< creates primal data */
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(primal != NULL);

   ALLOC_OKAY( allocMemory(*primal) );
   (*primal)->sols = NULL;
   (*primal)->solssize = 0;
   (*primal)->nsols = 0;
   (*primal)->nsolsfound = 0;
   (*primal)->upperbound = SCIP_INVALID;

   CHECK_OKAY( primalSetUpperbound(*primal, memhdr, set, NULL, lp, MIN(prob->objlim, set->infinity)) );

   return SCIP_OKAY;
}

RETCODE SCIPprimalFree(                 /**< frees primal data */
   PRIMAL**         primal,             /**< pointer to primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   int s;

   assert(primal != NULL);

   /* release primal CIP solutions */
   for( s = 0; s < (*primal)->nsols; ++s )
   {
      SCIPsolRelease(&(*primal)->sols[s], memhdr, set, lp);
   }
   freeMemoryArrayNull((*primal)->sols);
   freeMemory(*primal);

   return SCIP_OKAY;
}

RETCODE SCIPprimalAddSol(               /**< adds solution to primal solution storage */
   PRIMAL*          primal,             /**< primal data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< actual LP data */
   SOL**            sol                 /**< pointer to primal CIP solution */
   )
{
   int insertpos;
   int pos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);

   debugMessage("added primal solution:");
   debug( SCIPsolPrint(*sol, set, NULL) );

   SCIPsolCapture(*sol);

   /* search the position to insert solution in storage */
   for( insertpos = 0; insertpos < primal->nsols && (*sol)->obj >= primal->sols[insertpos]->obj; ++insertpos )
   {
   }

   if( insertpos < set->maxsol )
   {
      CHECK_OKAY( ensureSolsSize(primal, set, set->maxsol) );
      
      /* insert solution at correct position */
      primal->nsols = MIN(primal->nsols+1, set->maxsol);
      for( pos = primal->nsols-1; pos > insertpos; --pos )
         primal->sols[pos] = primal->sols[pos-1];

      primal->sols[insertpos] = *sol;
      primal->nsolsfound++;
      debugMessage(" -> stored at position %d of %d solutions, found %d solutions\n", 
         insertpos, primal->nsols, primal->nsolsfound);

      /* check, if the global upper bound has to be updated */
      if( (*sol)->obj < primal->upperbound )
      {
         /* update the upper bound */
         CHECK_OKAY( primalSetUpperbound(primal, memhdr, set, tree, lp, (*sol)->obj) );
         
         /* display node information line */
         CHECK_OKAY( SCIPdispPrintLine(set, stat, TRUE) );
      }
   }
   else
   {
      /* we don't need the solution -> release it */
      SCIPsolRelease(sol, memhdr, set, lp);
   }

   return SCIP_OKAY;
}
