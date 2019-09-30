/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   primalex.c
 * @brief  methods for collecting exact primal CIP solutions and primal informations
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/set.h"
#include "scip/stat.h"
#include "scip/visual.h"
#include "scip/event.h"
#include "scip/lp.h"
#include "scip/lpex.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/sol.h"
#include "scip/solex.h"
#include "scip/primal.h"
#include "scip/tree.h"
#include "scip/reopt.h"
#include "scip/disp.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/rational.h"
#include "scip/struct_solex.h"

#ifdef SCIP_WITH_GMP
#include "gmp.h"
#endif

/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that existingsols array can store at least num entries */
static
SCIP_RETCODE ensureExistingsolexsSize(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(primal->nexistingsols <= primal->existingsolssize);

   if( num > primal->existingsolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&primal->existingsols, newsize) );
      primal->existingsolssize = newsize;
   }
   assert(num <= primal->existingsolssize);

   return SCIP_OKAY;
}

/** uses binary search to find position in solution storage */
static
int primalexSearchSolPos(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            transprob,          /**< tranformed problem data */
   SCIP_PROB*            origprob,           /**< original problem data */
   SCIP_SOLEX*           sol                 /**< primal solution to search position for */
   )
{
   SCIP_SOLEX** sols;
   SCIP_Rational* obj;
   SCIP_Rational* middleobj;
   int left;
   int right;
   int middle;

   assert(primal != NULL);

   obj = SCIPsolexGetObj(sol, set, transprob, origprob);
   sols = primal->sols;

   left = -1;
   right = primal->nsols;
   while( left < right-1 )
   {
      middle = (left+right)/2;
      assert(left < middle && middle < right);
      assert(0 <= middle && middle < primal->nsols);

      middleobj = SCIPsolexGetObj(sols[middle], set, transprob, origprob);

      if( RisLT(obj, middleobj) )
         right = middle;
      else
         left = middle;
   }
   assert(left == right-1);

   /* prefer solutions that live in the transformed space */
   if( sol->solorigin != SCIP_SOLORIGIN_ORIGINAL )
   {
      while( right > 0 && sol->solorigin == SCIP_SOLORIGIN_ORIGINAL 
         && RisEqual(SCIPsolexGetObj(sols[right-1], set, transprob, origprob), obj) )
         --right;
   }

   return right;
}

/** returns whether the given primal solution is already existent in the solution storage */
static
SCIP_Bool primalexExistsSol(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_SOLEX*           sol,                /**< primal solution to search position for */
   int*                  insertpos,          /**< pointer to insertion position returned by primalSearchSolPos(); the
                                              *   position might be changed if an existing solution should be replaced */
   SCIP_Bool*            replace             /**< pointer to store whether the solution at insertpos should be replaced */
   )
{
   SCIP_Rational* obj;
   int i;

   assert(primal != NULL);
   assert(insertpos != NULL);
   assert(replace != NULL);
   assert(0 <= (*insertpos) && (*insertpos) <= primal->nsols);

   obj = SCIPsolexGetObj(sol, set, transprob, origprob);

   assert(primal->sols != NULL || primal->nsols == 0);
   assert(primal->sols != NULL || (*insertpos) == 0);

   /* search in the better solutions */
   for( i = (*insertpos)-1; i >= 0; --i )
   {
      SCIP_Rational* solobj;

      solobj = SCIPsolexGetObj(primal->sols[i], set, transprob, origprob);

      assert(RisLE(solobj, obj));

      if( RisLT(solobj, obj) )
         break;

      if( SCIPsolexsAreEqual(sol, primal->sols[i], set, stat, origprob, transprob) )
      {
         if( SCIPsolexIsOriginal(primal->sols[i]) && !SCIPsolexIsOriginal(sol) )
         {
            (*insertpos) = i;
            (*replace) = TRUE;
         }
         return TRUE;
      }
   }

   /* search in the worse solutions */
   for( i = (*insertpos); i < primal->nsols; ++i )
   {
      SCIP_Rational* solobj;

      solobj = SCIPsolexGetObj(primal->sols[i], set, transprob, origprob);

      /* due to transferring the objective value of transformed solutions to the original space, small numerical errors might occur
       * which can lead to SCIPsetIsLE() failing in case of high absolute numbers
       */
      assert( RisGE(solobj, obj));

      if( RisGT(solobj, obj) )
         break;

      if( SCIPsolexsAreEqual(sol, primal->sols[i], set, stat, origprob, transprob) )
      {
         if( SCIPsolexIsOriginal(primal->sols[i]) && !SCIPsolexIsOriginal(sol) )
         {
            (*insertpos) = i;
            (*replace) = TRUE;
         }
         return TRUE;
      }
   }

   return FALSE;
}

/** check if we are willing to check the solution for feasibility */
static
SCIP_Bool solOfInterest(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_SOLEX*           sol,                /**< primal CIP solution */
   int*                  insertpos,          /**< pointer to store the insert position of that solution */
   SCIP_Bool*            replace             /**< pointer to store whether the solution at insertpos should be replaced
                                              *   (e.g., because it lives in the original space) */
   )
{
   SCIP_Rational* obj;

   obj = SCIPsolexGetObj(sol, set, transprob, origprob);

   /* check if we are willing to check worse solutions; a solution is better if the objective is smaller than the
    * current cutoff bound; solutions with infinite objective value are never accepted
    */
   if( (!set->misc_improvingsols || RisLT(obj, primal->cutoffbound)) && !RisInfinity(obj) )
   {
      /* find insert position for the solution */
      (*insertpos) = primalexSearchSolPos(primal, set, transprob, origprob, sol);
      (*replace) = FALSE;

      /* the solution should be added, if the insertpos is smaller than the maximum number of solutions to be stored
       * and it does not already exist or it does exist, but the existing solution should be replaced by the new one
       */
      if( (*insertpos) < set->limit_maxsol && (!primalexExistsSol(primal, set, stat, origprob, transprob, sol, insertpos, replace) || (*replace)) )
         return TRUE;
   }

   return FALSE;
}

/** ensures, that sols array can store at least num entries */
static
SCIP_RETCODE ensureSolsSize(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(primal->nsols <= primal->solssize);

   if( num > primal->solssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&primal->sols, newsize) );
      primal->solssize = newsize;
   }
   assert(num <= primal->solssize);

   return SCIP_OKAY;
}

/** ensures, that existingsols array can store at least num entries */
static
SCIP_RETCODE ensureExistingsolsSize(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   num                 /**< minimum number of entries to store */
   )
{
   assert(primal->nexistingsols <= primal->existingsolssize);

   if( num > primal->existingsolssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      SCIP_ALLOC( BMSreallocMemoryArray(&primal->existingsols, newsize) );
      primal->existingsolssize = newsize;
   }
   assert(num <= primal->existingsolssize);

   return SCIP_OKAY;
}

/** adds exact primal solution to solution storage at given position */
static
SCIP_RETCODE primalexAddSol(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOLEX**          solptr,             /**< pointer to primal CIP solution */
   int                   insertpos,          /**< position in solution storage to add solution to */
   SCIP_Bool             replace             /**< should the solution at insertpos be replaced by the new solution? */
   )
{
   SCIP_SOLEX* sol;
   SCIP_Bool stored;
   /* cppcheck-suppress unassignedVariable */
   SCIP_EVENT event;
   SCIP_Rational* obj;
   SCIP_Real fpobj;
   int pos;

   assert(primal != NULL);
   assert(set != NULL);
   assert(solptr != NULL);
   assert(stat != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(0 <= insertpos && insertpos < set->limit_maxsol);
   assert(tree == NULL || !SCIPtreeInRepropagation(tree));

   sol = *solptr;
   assert(sol != NULL);
   obj = SCIPsolexGetObj(sol, set, transprob, origprob);
   fpobj = RgetRealRelax(obj, SCIP_ROUND_UPWARDS);

   SCIPsetDebugMsg(set, "insert exact primal solution %p with obj %g at position %d (replace=%u):\n",
      (void*)sol, RgetRealApprox(obj), insertpos, replace);

   SCIPdebug( SCIP_CALL( SCIPsolexPrint(sol, set, messagehdlr, stat, transprob, NULL, NULL, FALSE, FALSE) ) );

#if 0 /* this is not a valid debug check, but can be used to track down numerical troubles */
#ifndef NDEBUG
   /* check solution again completely
    * it fail for different reasons:
    * - in the LP solver, the feasibility tolerance is a relative measure against the row's norm
    * - in SCIP, the feasibility tolerance is a relative measure against the row's rhs/lhs
    * - the rhs/lhs of a row might drastically change during presolving when variables are fixed or (multi-)aggregated
    */
   if( !SCIPsolIsOriginal(sol) )
   {
      SCIP_Bool feasible;

      SCIP_CALL( SCIPsolCheck(sol, set, messagehdlr, blkmem, stat, transprob, TRUE, TRUE, TRUE, TRUE, &feasible) );

      if( !feasible )
      {
         SCIPerrorMessage("infeasible solution accepted:\n");
         SCIP_CALL( SCIPsolPrint(sol, set, messagehdlr, stat, origprob, transprob, NULL, FALSE, FALSE) );
      }
      assert(feasible);
   }
#endif
#endif

   /* completely fill the solution's own value array to unlink it from the LP or pseudo solution */
   SCIP_CALL( SCIPsolexUnlink(sol, set, transprob) );

#if 0
   /* allocate memory for solution storage */
   SCIP_CALL( ensureSolsSize(primal, set, set->limit_maxsol) );

   /* if set->limit_maxsol was decreased in the meantime, free all solutions exceeding the limit */
   for( pos = set->limit_maxsol; pos < primal->nsols; ++pos )
   {
      SCIP_CALL( SCIPsolexFree(&primal->sols[pos], blkmem, primal) );
   }
   primal->nsols = MIN(primal->nsols, set->limit_maxsol);

   /* if the solution should replace an existing one, free this solution, otherwise,
    * free the last solution if the solution storage is full;
    */
   if( replace )
   {
      //SCIP_CALL( SCIPsolexTransform(primal->sols[insertpos], solptr, blkmem, set, primal) );
      sol = primal->sols[insertpos];
   }
   else
   {
      if( primal->nsols == set->limit_maxsol )
      {
         SCIP_CALL( SCIPsolexFree(&primal->sols[set->limit_maxsol - 1], blkmem, primal) );
      }
      else
      {
         primal->nsols = primal->nsols + 1;
         assert(primal->nsols <= set->limit_maxsol);
      }

      /* move all solutions with worse objective value than the new solution */
      for( pos = primal->nsols-1; pos > insertpos; --pos )
         primal->sols[pos] = primal->sols[pos-1];

      /* insert solution at correct position */
      assert(0 <= insertpos && insertpos < primal->nsols);
      primal->sols[insertpos] = sol;
      primal->nsolsfound++;

      /** @todo: exip, is this correct? */
      /* check if solution is better than objective limit */
      if( SCIPsetIsFeasLE(set, fpobj, SCIPprobInternObjval(transprob, origprob, set, SCIPprobGetObjlim(origprob, set))) )
         primal->nlimsolsfound++;
   }

   SCIPsetDebugMsg(set, " -> stored at position %d of %d solutions, found %" SCIP_LONGINT_FORMAT " solutions\n",
      insertpos, primal->nsols, primal->nsolsfound);

   /* update the solution value sums in variables */
   /* if( !SCIPsolIsOriginal(sol) )
   {
      SCIPsolUpdateVarsum(sol, set, stat, transprob,
         (SCIP_Real)(primal->nsols - insertpos)/(SCIP_Real)(2.0*primal->nsols - 1.0));
   }
 */
   /* change color of node in visualization output */
   // SCIPvisualFoundSolution(stat->visual, set, stat, SCIPtreeGetCurrentNode(tree), insertpos == 0 ? TRUE : FALSE, sol);
   #endif
   SCIP_CALL( SCIPsolexOverwriteFPSol(sol->fpsol, sol, set, stat, origprob, transprob, tree) );
   /* note: we copy the solution so to not destroy the double-link between sol and fpsol */
   SCIP_CALL( SCIPprimalAddSolFree(primal->fpstorage, blkmem, set, messagehdlr, stat,
            origprob, transprob, tree, reopt,
            lp->fplp, eventqueue, eventfilter, &sol->fpsol, &stored) );

   /* check, if the global upper bound has to be updated */
   if( RisLT(obj, primal->cutoffbound) && insertpos == 0 )
   {
      Rset(primal->upperbound, obj);
      Rset(primal->cutoffbound, obj);

      primal->nbestsolsfound++;
      stat->bestsolnode = stat->nnodes;
   }
   // else
   // {
   //    /* issue POORSOLFOUND event */
   //    SCIP_CALL( SCIPeventChgType(&event, SCIP_EVENTTYPE_POORSOLFOUND) );
   // }
   // SCIP_CALL( SCIPeventChgSol(&event, sol->fpsol) );
   // SCIP_CALL( SCIPeventProcess(&event, set, NULL, NULL, NULL, eventfilter) );

   /* if an original solution was added during solving, try to transfer it to the transformed space */
   /* if( SCIPsolIsOriginal(sol) && SCIPsetGetStage(set) == SCIP_STAGE_SOLVING && set->misc_transorigsols )
   {
      SCIP_Bool added;

      SCIP_CALL( SCIPprimalTransformSol(primal, sol, blkmem, set, messagehdlr, stat, origprob, transprob, tree, reopt,
            lp, eventqueue, eventfilter, NULL, NULL, 0, &added) );

      SCIPsetDebugMsg(set, "original solution %p was successfully transferred to the transformed problem space\n",
         (void*)sol);
   } */

   return SCIP_OKAY;
}

/** creates exact primal data */
SCIP_RETCODE SCIPprimalexCreate(
   SCIP_PRIMALEX**       primal,             /**< pointer to exact primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_PRIMAL*          fpstorage           /**< the fp primal storage */
   )
{
   assert(primal != NULL);

   SCIP_ALLOC( BMSallocMemory(primal) );
   (*primal)->sols = NULL;
   (*primal)->existingsols = NULL;
   (*primal)->solssize = 0;
   (*primal)->nsolsfound = 0;
   (*primal)->nlimsolsfound = 0;
   (*primal)->nbestsolsfound = 0;
   (*primal)->nlimbestsolsfound = 0;
   (*primal)->solssize = 0;
   (*primal)->partialsolssize = 0;
   (*primal)->nsols = 0;
   (*primal)->npartialsols = 0;
   (*primal)->existingsolssize = 0;
   (*primal)->nexistingsols = 0;
   (*primal)->updateviolations = TRUE;
   /* double-link the two storages */
   (*primal)->fpstorage = fpstorage;
   fpstorage->primalex = (*primal);
   SCIP_CALL( RcreateString(blkmem, &(*primal)->upperbound, "inf") );
   SCIP_CALL( RcreateString(blkmem, &(*primal)->cutoffbound, "inf") );

   return SCIP_OKAY;
}

/** frees exact primal data */
SCIP_RETCODE SCIPprimalexFree(
   SCIP_PRIMALEX**       primal,             /**< pointer to exact primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int s;

   assert(primal != NULL);
   assert(*primal != NULL);

   /* free feasible primal CIP solutions */
   for( s = 0; s < (*primal)->nsols; ++s )
   {
      SCIP_CALL( SCIPsolexFree(&(*primal)->sols[s], blkmem, *primal) );
   }

   Rdelete(blkmem, &(*primal)->cutoffbound);
   Rdelete(blkmem, &(*primal)->upperbound);

   BMSfreeMemoryArrayNull(&(*primal)->existingsols);
   BMSfreeMemoryArrayNull(&(*primal)->sols);
   BMSfreeMemory(primal);

   return SCIP_OKAY;
}

/** adds exact primal solution to solution storage, frees the solution afterwards */
SCIP_RETCODE SCIPprimalexTrySolFree(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   SCIP_STAT*            stat,               /**< problem statistics data */
   SCIP_PROB*            origprob,           /**< original problem */
   SCIP_PROB*            transprob,          /**< transformed problem after presolve */
   SCIP_TREE*            tree,               /**< branch and bound tree */
   SCIP_REOPT*           reopt,              /**< reoptimization data structure */
   SCIP_LPEX*            lp,                 /**< current LP data */
   SCIP_EVENTQUEUE*      eventqueue,         /**< event queue */
   SCIP_EVENTFILTER*     eventfilter,        /**< event filter for global (not variable dependent) events */
   SCIP_SOLEX**          sol,                /**< pointer to primal CIP solution; is cleared in function call */
   SCIP_Bool             printreason,        /**< Should all the reasons of violations be printed? */
   SCIP_Bool             completely,         /**< Should all violations be checked? */
   SCIP_Bool             checkbounds,        /**< Should the bounds of the variables be checked? */
   SCIP_Bool             checkintegrality,   /**< Has integrality to be checked? */
   SCIP_Bool             checklprows,        /**< Do constraints represented by rows in the current LP have to be checked? */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   SCIP_Bool replace = FALSE;
   SCIP_Bool feasible;
   int insertpos;

   assert(primal != NULL);
   assert(transprob != NULL);
   assert(origprob != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(stored != NULL);

   insertpos = -1;

   /* insert solution into solution storage */

   if( solOfInterest(primal, set, stat, origprob, transprob, *sol, &insertpos, &replace) )
   {
      /* check solution for feasibility */
      SCIP_CALL( SCIPsolexCheck(*sol, set, messagehdlr, blkmem, stat, transprob, printreason, completely, checkbounds,
            checkintegrality, checklprows, &feasible) );
   }
   else
      feasible  = FALSE;

   if( feasible )
   {
      SCIP_CALL( primalexAddSol(primal, blkmem, set, messagehdlr, stat, origprob, transprob,
            tree, reopt, lp, eventqueue, eventfilter, sol, insertpos, replace) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;
      *stored = TRUE;
   }
   else
   {
      SCIP_SOL* fpsol;
      /* the solution is too bad or infeasible -> free it immediately */
      fpsol = (*sol)->fpsol;
      SCIP_CALL( SCIPsolexFree(sol, blkmem, primal) );
      SCIP_CALL( SCIPsolFree(&fpsol, blkmem, primal->fpstorage) );
   
      *stored = FALSE;
   }

   assert(*sol == NULL);

   return SCIP_OKAY;
}

/** inserts solution into the global array of all existing primal solutions */
SCIP_RETCODE SCIPprimalexSolexCreated(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_SOLEX*           sol                 /**< primal CIP solution */
   )
{
   assert(primal != NULL);
   assert(sol != NULL);
   assert(SCIPsolexGetPrimalexIndex(sol) == -1);

   /* allocate memory for solution storage */
   SCIP_CALL( ensureExistingsolexsSize(primal, set, primal->nexistingsols+1) );

   /* append solution */
   SCIPsolexSetPrimalexIndex(sol, primal->nexistingsols);
   primal->existingsols[primal->nexistingsols] = sol;
   primal->nexistingsols++;

   return SCIP_OKAY;
}

/** removes solution from the global array of all existing primal solutions */
void SCIPprimalexSolexFreed(
   SCIP_PRIMALEX*        primal,             /**< primal data */
   SCIP_SOLEX*           sol                 /**< primal CIP solution */
   )
{
   int idx;

   assert(primal != NULL);
   assert(sol != NULL);

#ifndef NDEBUG
   for( idx = 0; idx < primal->nexistingsols; ++idx )
   {
      assert(idx == SCIPsolexGetPrimalexIndex(primal->existingsols[idx]));
   }
#endif

   /* remove solution */
   idx = SCIPsolexGetPrimalexIndex(sol);
   assert(0 <= idx && idx < primal->nexistingsols);
   assert(sol == primal->existingsols[idx]);
   if( idx < primal->nexistingsols-1 )
   {
      primal->existingsols[idx] = primal->existingsols[primal->nexistingsols-1];
      SCIPsolexSetPrimalexIndex(primal->existingsols[idx], idx);
   }
   primal->nexistingsols--;
}

/** clears primal data */
SCIP_RETCODE SCIPprimalexClear(
   SCIP_PRIMALEX**       primal,             /**< pointer to primal data */
   BMS_BLKMEM*           blkmem              /**< block memory */
   )
{
   int s;

   assert(primal != NULL);
   assert(*primal != NULL);

   /* free temporary solution for storing current solution */
   if( (*primal)->currentsol != NULL )
   {
      SCIP_CALL( SCIPsolexFree(&(*primal)->currentsol, blkmem, *primal) );
   }

   /* free solution for storing primal ray */
   if( (*primal)->primalray != NULL )
   {
      SCIP_CALL( SCIPsolexFree(&(*primal)->primalray, blkmem, *primal) );
   }

   /* free feasible primal CIP solutions */
   for( s = 0; s < (*primal)->nsols; ++s )
   {
      SCIP_CALL( SCIPsolexFree(&(*primal)->sols[s], blkmem, *primal) );
   }

   (*primal)->currentsol = NULL;
   (*primal)->primalray = NULL;
   (*primal)->nsols = 0;
   (*primal)->nsolsfound = 0;
   (*primal)->nlimsolsfound = 0;
   (*primal)->nbestsolsfound = 0;
   (*primal)->nlimbestsolsfound = 0;
   (*primal)->updateviolations = TRUE;

   return SCIP_OKAY;
}