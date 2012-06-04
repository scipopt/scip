/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   primalex.c
 * @brief  methods for collecting exact primal CIP solutions and exact primal information
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/def.h"
#include "scip/message.h"
#include "scip/set.h"
#include "scip/var.h"
#include "scip/prob.h"
#include "scip/solex.h"
#include "scip/primalex.h"


#ifdef WITH_EXACTSOLVE
#include "gmp.h"

/*
 * memory growing methods for dynamically allocated arrays
 */

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

/** uses binary search to find position in solution storage */
static
int primalexSearchSolPos(
   SCIP_PRIMALEX*        primal,             /**< exact primal data */
   SCIP_SOLEX*           sol                 /**< exact primal solution to search position for */
   )
{
   mpq_t obj;
   mpq_t middleobj;
   int left;
   int right;
   int middle;

   assert(primal != NULL);

   mpq_init(obj);
   mpq_init(middleobj);

   SCIPsolexGetObj(sol, obj);

   left = -1;
   right = primal->nsols;
   while( left < right-1 )
   {
      middle = (left+right)/2;
      assert(left < middle && middle < right);
      assert(0 <= middle && middle < primal->nsols);
      SCIPsolexGetObj(primal->sols[middle], middleobj);

      if( mpq_cmp(obj, middleobj) < 0 )
         right = middle;
      else
         left = middle;
   }
   assert(left == right-1);

   return right;

   mpq_clear(middleobj);
   mpq_clear(obj);
}

/** returns whether the given exact primal solution is already existant in the solution storage */
static
SCIP_Bool primalexExistsSol(
   SCIP_PRIMALEX*        primal,             /**< exact primal data */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_SOLEX*           sol,                /**< exact primal solution to search position for */
   int                   insertpos           /**< insertion position returned by primalexSearchSolPos() */
   )
{
   SCIP_Bool exists;
   mpq_t obj;
   mpq_t solobj;
   int i;

   assert(primal != NULL);
   assert(0 <= insertpos && insertpos <= primal->nsols);

   exists = FALSE;

   /* original solutions are always accepted */
   if( SCIPsolexGetOrigin(sol) == SCIP_SOLORIGIN_ORIGINAL )
      return FALSE;

   mpq_init(obj);
   mpq_init(solobj);

   SCIPsolexGetObj(sol, obj);

   /* search in the better solutions */
   for( i = insertpos-1; i >= 0 && !exists; --i )
   {
      SCIPsolexGetObj(primal->sols[i], solobj);
      assert(mpq_cmp(solobj, obj) <= 0);
      if( mpq_cmp(solobj, obj) < 0 )
         break;
      if( SCIPsolexGetOrigin(primal->sols[i]) != SCIP_SOLORIGIN_ORIGINAL && SCIPsolexsAreEqual(sol, primal->sols[i], prob) )
         exists = TRUE;
   }

   /* search in the worse solutions */
   for( i = insertpos; i < primal->nsols; ++i )
   {
      SCIPsolexGetObj(primal->sols[i], solobj);
      assert(mpq_cmp(solobj, obj) >= 0);
      if( mpq_cmp(solobj, obj) > 0 )
         break;
      if( SCIPsolexGetOrigin(primal->sols[i]) != SCIP_SOLORIGIN_ORIGINAL && SCIPsolexsAreEqual(sol, primal->sols[i], prob) )
         exists = TRUE;
   }

   mpq_clear(solobj);
   mpq_clear(obj);

   return exists;
}


/** adds exact primal solution to solution storage at given position */
static
SCIP_RETCODE primalexAddSol(
   SCIP_PRIMALEX*        primal,             /**< exact primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_SOLEX*           sol,                /**< exact primal CIP solution */
   int                   insertpos           /**< position in solution storage to add solution to */
   )
{
   int pos;
#ifdef SCIP_DEBUG
   char s[SCIP_MAXSTRLEN];
   mpq_t obj;
#endif

   assert(primal != NULL);
   assert(sol != NULL);
   assert(0 <= insertpos && insertpos < set->limit_maxsol);

#ifdef SCIP_DEBUG
   mpq_init(obj);
   SCIPsolexGetObj(sol, obj);
   gmp_snprintf(s, SCIP_MAXSTRLEN, "insert primal solution %p with obj %Qd at position %d:\n", (void*)sol, obj, insertpos);
   SCIPdebugMessage(s);
   SCIPdebug( SCIPsolexPrint(sol, prob, NULL, NULL, FALSE) );
   mpq_clear(obj);
#endif

   /* allocate memory for solution storage */
   SCIP_CALL( ensureSolsSize(primal, set, set->limit_maxsol) );

   /* if the solution storage is full, free the last solution(s)
    * more than one solution may be freed, if set->limit_maxsol was decreased in the meantime
    */
   for( pos = set->limit_maxsol-1; pos < primal->nsols; ++pos )
   {
      SCIP_CALL( SCIPsolexFree(&primal->sols[pos], blkmem) );
   }

   /* insert solution at correct position */
   primal->nsols = MIN(primal->nsols+1, set->limit_maxsol);
   for( pos = primal->nsols-1; pos > insertpos; --pos )
      primal->sols[pos] = primal->sols[pos-1];

   assert(0 <= insertpos && insertpos < primal->nsols);
   primal->sols[insertpos] = sol;
   SCIPdebugMessage(" -> stored at position %d of %d solutions\n", insertpos, primal->nsols);

   return SCIP_OKAY;
}

/** creates exact primal data */
SCIP_RETCODE SCIPprimalexCreate(
   SCIP_PRIMALEX**       primal              /**< pointer to exact primal data */
   )
{
   assert(primal != NULL);

   SCIP_ALLOC( BMSallocMemory(primal) );
   (*primal)->sols = NULL;
   (*primal)->solssize = 0;
   (*primal)->nsols = 0;

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
      SCIP_CALL( SCIPsolexFree(&(*primal)->sols[s], blkmem) );
   }

   BMSfreeMemoryArrayNull(&(*primal)->sols);
   BMSfreeMemory(primal);

   return SCIP_OKAY;
}

/** adds exact primal solution to solution storage, frees the solution afterwards */
SCIP_RETCODE SCIPprimalexAddSolFree(
   SCIP_PRIMALEX*        primal,             /**< exact primal data */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PROB*            prob,               /**< transformed problem after presolve */
   SCIP_SOLEX**          sol,                /**< pointer to exact primal CIP solution; is cleared in function call */
   SCIP_Bool*            stored              /**< stores whether given solution was good enough to keep */
   )
{
   int insertpos;

   assert(primal != NULL);
   assert(sol != NULL);
   assert(*sol != NULL);
   assert(stored != NULL);

   /* search the position to insert solution in storage */
   insertpos = primalexSearchSolPos(primal, *sol);

   if( insertpos < set->limit_maxsol && !primalexExistsSol(primal, prob, *sol, insertpos) )
   {
      /* insert solution into solution storage */
      SCIP_CALL( primalexAddSol(primal, blkmem, set, prob, *sol, insertpos) );

      /* clear the pointer, such that the user cannot access the solution anymore */
      *sol = NULL;

      *stored = TRUE;
   }
   else
   {
      /* the solution is too bad -> free it immediately */
      SCIP_CALL( SCIPsolexFree(sol, blkmem) );

      *stored = FALSE;
   }
   assert(*sol == NULL);

   return SCIP_OKAY;
}

#endif
