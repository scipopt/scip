/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch.c,v 1.54 2004/11/19 17:27:23 bzfpfend Exp $"

/**@file   branch.c
 * @brief  methods for branching rules and branching candidate storage
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "memory.h"
#include "set.h"
#include "stat.h"
#include "clock.h"
#include "paramset.h"
#include "event.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "tree.h"
#include "sepastore.h"
#include "scip.h"
#include "branch.h"

#include "struct_branch.h"



/*
 * memory growing methods for dynamically allocated arrays
 */

/** ensures, that lpcands array can store at least num entries */
static
RETCODE ensureLpcandsSize(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(branchcand->nlpcands <= branchcand->lpcandssize);
   
   if( num > branchcand->lpcandssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&branchcand->lpcands, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&branchcand->lpcandssol, newsize) );
      ALLOC_OKAY( reallocMemoryArray(&branchcand->lpcandsfrac, newsize) );
      branchcand->lpcandssize = newsize;
   }
   assert(num <= branchcand->lpcandssize);

   return SCIP_OKAY;
}

/** ensures, that pseudocands array can store at least num entries */
static
RETCODE ensurePseudocandsSize(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   )
{
   assert(branchcand->npseudocands <= branchcand->pseudocandssize);
   
   if( num > branchcand->pseudocandssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&branchcand->pseudocands, newsize) );
      branchcand->pseudocandssize = newsize;
   }
   assert(num <= branchcand->pseudocandssize);

   return SCIP_OKAY;
}



/*
 * branching candidate storage methods
 */

/** creates a branching candidate storage */
RETCODE SCIPbranchcandCreate(
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   )
{
   assert(branchcand != NULL);

   ALLOC_OKAY( allocMemory(branchcand) );
   (*branchcand)->lpcands = NULL;
   (*branchcand)->lpcandssol = NULL;
   (*branchcand)->lpcandsfrac = NULL;
   (*branchcand)->pseudocands = NULL;
   (*branchcand)->lpcandssize = 0;
   (*branchcand)->nlpcands = 0;
   (*branchcand)->npriolpcands = 0;
   (*branchcand)->npriolpbins = 0;
   (*branchcand)->lpmaxpriority = INT_MIN;
   (*branchcand)->pseudocandssize = 0;
   (*branchcand)->npseudocands = 0;
   (*branchcand)->npriopseudocands = 0;
   (*branchcand)->npriopseudobins = 0;
   (*branchcand)->npriopseudoints = 0;
   (*branchcand)->pseudomaxpriority = INT_MIN;
   (*branchcand)->validlpcandslp = -1;
   
   return SCIP_OKAY;
}

/** frees branching candidate storage */
RETCODE SCIPbranchcandFree(
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   )
{
   assert(branchcand != NULL);

   freeMemoryArrayNull(&(*branchcand)->lpcands);
   freeMemoryArrayNull(&(*branchcand)->lpcandssol);
   freeMemoryArrayNull(&(*branchcand)->lpcandsfrac);
   freeMemoryArrayNull(&(*branchcand)->pseudocands);
   freeMemory(branchcand);

   return SCIP_OKAY;
}

/** calculates branching candidates for LP solution branching (fractional variables) */
static
RETCODE branchcandCalcLPCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp                  /**< current LP data */
   )
{
   assert(branchcand != NULL);
   assert(stat != NULL);
   assert(branchcand->validlpcandslp <= stat->lpcount);
   assert(lp != NULL);
   assert(lp->solved);
   assert(SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_OPTIMAL || SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_UNBOUNDEDRAY);

   debugMessage("calculating LP branching candidates: validlp=%d, lpcount=%d\n",
      branchcand->validlpcandslp, stat->lpcount);

   /* check, if the current LP branching candidate array is invalid */
   if( branchcand->validlpcandslp < stat->lpcount )
   {
      COL** cols;
      VAR* var;
      COL* col;
      Real primsol;
      Real frac;
      VARTYPE vartype;
      int branchpriority;
      int ncols;
      int c;
      int insertpos;

      debugMessage(" -> recalculating LP branching candidates\n");

      cols = SCIPlpGetCols(lp);
      ncols = SCIPlpGetNCols(lp);

      /* construct the LP branching candidate set, moving the candidates with maximal priority to the front */
      CHECK_OKAY( ensureLpcandsSize(branchcand, set, ncols) );

      branchcand->lpmaxpriority = INT_MIN;
      branchcand->nlpcands = 0;
      branchcand->npriolpcands = 0;
      branchcand->npriolpbins = 0;
      for( c = 0; c < ncols; ++c )
      {
         col = cols[c];
         assert(col != NULL);
         assert(col->lppos == c);
         assert(col->lpipos >= 0);

         primsol = SCIPcolGetPrimsol(col);
         assert(primsol < SCIP_INVALID);
         assert(SCIPsetIsFeasGE(set, primsol, col->lb));
         assert(SCIPsetIsFeasLE(set, primsol, col->ub));

         var = col->var;
         assert(var != NULL);
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(var) == col);
         
         /* LP branching candidates are fractional binary and integer variables */
         vartype = SCIPvarGetType(var);
         if( vartype == SCIP_VARTYPE_BINARY || vartype == SCIP_VARTYPE_INTEGER )
         {
            frac = SCIPsetFrac(set, primsol);
            if( !SCIPsetIsFracIntegral(set, frac) )
            {
               assert(branchcand->nlpcands < branchcand->lpcandssize);

               /* insert candidate in candidate list */
               branchpriority = SCIPvarGetBranchPriority(var);
               insertpos = branchcand->nlpcands;
               branchcand->nlpcands++;
               if( branchpriority > branchcand->lpmaxpriority )
               {
                  /* candidate has higher priority than the current maximum:
                   * move it to the front and declare it to be the single best candidate
                   */
                  if( insertpos != 0 )
                  {
                     branchcand->lpcands[insertpos] = branchcand->lpcands[0];
                     branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[0];
                     branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[0];
                     insertpos = 0;
                  }
                  branchcand->npriolpcands = 1;
                  branchcand->npriolpbins = (vartype == SCIP_VARTYPE_BINARY ? 1 : 0);
                  branchcand->lpmaxpriority = branchpriority;
               }
               else if( branchpriority == branchcand->lpmaxpriority )
               {
                  /* candidate has equal priority as the current maximum:
                   * move away the first non-maximal priority candidate, move the current candidate to the correct
                   * slot (binaries first) and increase the number of maximal priority candidates
                   */
                  if( insertpos != branchcand->npriolpcands )
                  {
                     branchcand->lpcands[insertpos] = branchcand->lpcands[branchcand->npriolpcands];
                     branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[branchcand->npriolpcands];
                     branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[branchcand->npriolpcands];
                     insertpos = branchcand->npriolpcands;
                  }
                  branchcand->npriolpcands++;
                  if( vartype == SCIP_VARTYPE_BINARY )
                  {
                     if( insertpos != branchcand->npriolpbins )
                     {
                        branchcand->lpcands[insertpos] = branchcand->lpcands[branchcand->npriolpbins];
                        branchcand->lpcandssol[insertpos] = branchcand->lpcandssol[branchcand->npriolpbins];
                        branchcand->lpcandsfrac[insertpos] = branchcand->lpcandsfrac[branchcand->npriolpbins];
                        insertpos = branchcand->npriolpbins;
                     }
                     branchcand->npriolpbins++;
                  }
               }
               branchcand->lpcands[insertpos] = var;
               branchcand->lpcandssol[insertpos] = primsol;
               branchcand->lpcandsfrac[insertpos] = frac;

               debugMessage(" -> candidate %d: var=<%s>, sol=%g, frac=%g, prio=%d (max: %d) -> pos %d\n", 
                  branchcand->nlpcands, SCIPvarGetName(var), primsol, frac, branchpriority, branchcand->lpmaxpriority,
                  insertpos);
            }
         }
      }

      branchcand->validlpcandslp = stat->lpcount;
   }
   assert(0 <= branchcand->npriolpcands && branchcand->npriolpcands <= branchcand->nlpcands);

   debugMessage(" -> %d fractional variables (%d of maximal priority)\n", branchcand->nlpcands, branchcand->npriolpcands);

   return SCIP_OKAY;
}

/** gets branching candidates for LP solution branching (fractional variables) */
RETCODE SCIPbranchcandGetLPCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< current LP data */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands,           /**< pointer to store the number of LP branching candidates, or NULL */
   int*             npriolpcands        /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   /* calculate branching candidates */
   CHECK_OKAY( branchcandCalcLPCands(branchcand, set, stat, lp) );

   /* assign return values */
   if( lpcands != NULL )
      *lpcands = branchcand->lpcands;
   if( lpcandssol != NULL )
      *lpcandssol = branchcand->lpcandssol;
   if( lpcandsfrac != NULL )
      *lpcandsfrac = branchcand->lpcandsfrac;
   if( nlpcands != NULL )
      *nlpcands = branchcand->nlpcands;
   if( npriolpcands != NULL )
      *npriolpcands = (set->branch_preferbinary && branchcand->npriolpbins > 0 ? branchcand->npriolpbins
         : branchcand->npriolpcands);

   return SCIP_OKAY;
}

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
RETCODE SCIPbranchcandGetPseudoCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands,       /**< pointer to store the number of pseudo branching candidates, or NULL */
   int*             npriopseudocands    /**< pointer to store the number of candidates with maximal priority, or NULL */
   )
{
   assert(branchcand != NULL);

#ifndef NDEBUG
   /* check, if the current pseudo branching candidate array is correct */
   {
      VAR* var;
      int npcs;
      int v;
      
      assert(prob != NULL);
      
      /* pseudo branching candidates are non-fixed binary, integer, and implicit integer variables */
      npcs = 0;
      for( v = 0; v < prob->nbinvars + prob->nintvars + prob->nimplvars; ++v )
      {
         var = prob->vars[v];
         assert(var != NULL);
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY
            || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER
            || SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT);
         assert(SCIPsetIsIntegral(set, SCIPvarGetLbLocal(var)));
         assert(SCIPsetIsIntegral(set, SCIPvarGetUbLocal(var)));
         assert(SCIPsetIsLE(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

         if( SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
         {
            assert(0 <= var->pseudocandindex && var->pseudocandindex < branchcand->npseudocands);
            assert(branchcand->pseudocands[var->pseudocandindex] == var);
            npcs++;
         }
         else
         {
            assert(var->pseudocandindex == -1);
         }
      }
      assert(branchcand->npseudocands == npcs);
   }
#endif

   /* assign return values */
   if( pseudocands != NULL )
      *pseudocands = branchcand->pseudocands;
   if( npseudocands != NULL )
      *npseudocands = branchcand->npseudocands;
   if( npriopseudocands != NULL )
      *npriopseudocands = (set->branch_preferbinary && branchcand->npriopseudobins > 0 ? branchcand->npriopseudobins
         : branchcand->npriopseudocands);

   return SCIP_OKAY;
}

/** gets number of branching candidates for pseudo solution branching (nonfixed variables) */
int SCIPbranchcandGetNPseudoCands(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npseudocands;
}

/** gets number of branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoCands(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriopseudocands;
}

/** gets number of binary branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoBins(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriopseudobins;
}

/** gets number of integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoInts(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriopseudoints;
}

/** gets number of implicit integer branching candidates with maximal branch priority for pseudo solution branching */
int SCIPbranchcandGetNPrioPseudoImpls(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   assert(branchcand != NULL);

   return branchcand->npriopseudocands - branchcand->npriopseudobins - branchcand->npriopseudoints;
}

/** insert pseudocand at given position, or to the first positions of the maximal priority candidates, using the
 *  given position as free slot for the other candidates
 */
static
void branchcandInsertPseudoCand(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var,                /**< variable to insert */
   int              insertpos           /**< free position to insert the variable */
   )
{
   VARTYPE vartype;
   int branchpriority;

   assert(branchcand != NULL);
   assert(var != NULL);
   assert(branchcand->npriopseudocands <= insertpos && insertpos < branchcand->npseudocands);
   assert(branchcand->npseudocands <= branchcand->pseudocandssize);

   vartype = SCIPvarGetType(var);
   branchpriority = SCIPvarGetBranchPriority(var);

   debugMessage("inserting pseudo candidate <%s> of type %d and priority %d into candidate set at position %d (maxprio: %d)\n",
      SCIPvarGetName(var), vartype, branchpriority, insertpos, branchcand->pseudomaxpriority);

   /* insert the variable into pseudocands, making sure, that the highest priority candidates are at the front
    * and ordered binaries, integers, implicit integers
    */
   if( branchpriority > branchcand->pseudomaxpriority )
   {
      /* candidate has higher priority than the current maximum:
       * move it to the front and declare it to be the single best candidate
       */
      if( insertpos != 0 )
      {
         branchcand->pseudocands[insertpos] = branchcand->pseudocands[0];
         branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
         insertpos = 0;
      }
      branchcand->npriopseudocands = 1;
      branchcand->npriopseudobins = (vartype == SCIP_VARTYPE_BINARY ? 1 : 0);
      branchcand->npriopseudoints = (vartype == SCIP_VARTYPE_INTEGER ? 1 : 0);
      branchcand->pseudomaxpriority = branchpriority;
   }
   else if( branchpriority == branchcand->pseudomaxpriority )
   {
      /* candidate has equal priority as the current maximum:
       * move away the first non-maximal priority candidate, move the current candidate to the correct
       * slot (binaries first, integers next, implicits last) and increase the number of maximal priority candidates
       */
      if( insertpos != branchcand->npriopseudocands )
      {
         branchcand->pseudocands[insertpos] = branchcand->pseudocands[branchcand->npriopseudocands];
         branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
         insertpos = branchcand->npriopseudocands;
      }
      branchcand->npriopseudocands++;
      if( vartype == SCIP_VARTYPE_BINARY || vartype == SCIP_VARTYPE_INTEGER )
      {
         if( insertpos != branchcand->npriopseudobins + branchcand->npriopseudoints )
         {
            branchcand->pseudocands[insertpos] =
               branchcand->pseudocands[branchcand->npriopseudobins + branchcand->npriopseudoints];
            branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
            insertpos = branchcand->npriopseudobins + branchcand->npriopseudoints;
         }
         branchcand->npriopseudoints++;

         if( vartype == SCIP_VARTYPE_BINARY )
         {
            if( insertpos != branchcand->npriopseudobins )
            {
               branchcand->pseudocands[insertpos] = branchcand->pseudocands[branchcand->npriopseudobins];
               branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
               insertpos = branchcand->npriopseudobins;
            }
            branchcand->npriopseudobins++;
            branchcand->npriopseudoints--;
         }
      }
   }
   branchcand->pseudocands[insertpos] = var;
   var->pseudocandindex = insertpos;

   debugMessage(" -> inserted at position %d (npriopseudocands=%d)\n", insertpos, branchcand->npriopseudocands);

   assert(0 <= branchcand->npriopseudocands && branchcand->npriopseudocands <= branchcand->npseudocands);
   assert(0 <= branchcand->npriopseudobins && branchcand->npriopseudobins <= branchcand->npriopseudocands);
   assert(0 <= branchcand->npriopseudoints && branchcand->npriopseudoints <= branchcand->npriopseudocands);
}

/** sorts the pseudo branching candidates, such that the candidates of maximal priority are at the front,
 *  ordered by binaries, integers, implicit integers
 */
static
void branchcandSortPseudoCands(
   BRANCHCAND*      branchcand          /**< branching candidate storage */
   )
{
   VAR* var;
   int i;

   assert(branchcand != NULL);
   assert(branchcand->npriopseudocands == 0); /* is only be called after removal of last maximal candidate */
   assert(branchcand->npriopseudobins == 0);
   assert(branchcand->npriopseudoints == 0);

   debugMessage("resorting pseudo candidates\n");

   branchcand->pseudomaxpriority = INT_MIN;
   
   for( i = 0; i < branchcand->npseudocands; ++i )
   {
      var = branchcand->pseudocands[i];
      assert(var->pseudocandindex == i);

      if( SCIPvarGetBranchPriority(var) >= branchcand->pseudomaxpriority )
         branchcandInsertPseudoCand(branchcand, var, i);
   }

   assert(0 <= branchcand->npriopseudocands && branchcand->npriopseudocands <= branchcand->npseudocands);
   assert(0 <= branchcand->npriopseudobins && branchcand->npriopseudobins <= branchcand->npriopseudocands);
   assert(0 <= branchcand->npriopseudoints && branchcand->npriopseudoints <= branchcand->npriopseudocands);
}

/** removes pseudo candidate from pseudocands array
 */
static
void branchcandRemovePseudoCand(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< variable to remove */
   )
{
   VARTYPE vartype;
   int branchpriority;
   int freepos;

   assert(branchcand != NULL);
   assert(var != NULL);
   assert(var->pseudocandindex < branchcand->npseudocands);
   assert(branchcand->pseudocands[var->pseudocandindex] == var);
   assert(branchcand->pseudocands[branchcand->npseudocands-1] != NULL);

   vartype = SCIPvarGetType(var);
   branchpriority = SCIPvarGetBranchPriority(var);

   debugMessage("removing pseudo candidate <%s> of type %d and priority %d at %d from candidate set (maxprio: %d)\n",
      SCIPvarGetName(var), vartype, branchpriority, var->pseudocandindex, branchcand->pseudomaxpriority);

   /* delete the variable from pseudocands, making sure, that the highest priority candidates are at the front
    * and ordered binaries, integers, implicit integers
    */
   freepos = var->pseudocandindex;
   var->pseudocandindex = -1;
   assert(0 <= freepos && freepos < branchcand->npseudocands);

   if( freepos < branchcand->npriopseudobins )
   {
      /* a binary candidate of maximal priority was removed */
      assert(vartype == SCIP_VARTYPE_BINARY);
      assert(branchpriority == branchcand->pseudomaxpriority);
      if( freepos != branchcand->npriopseudobins - 1 )
      {
         branchcand->pseudocands[freepos] = branchcand->pseudocands[branchcand->npriopseudobins - 1];
         branchcand->pseudocands[freepos]->pseudocandindex = freepos;
         freepos = branchcand->npriopseudobins - 1;
      }
      branchcand->npriopseudobins--;
      branchcand->npriopseudoints++;
   }
   if( freepos < branchcand->npriopseudobins + branchcand->npriopseudoints )
   {
      /* a binary or integer candidate of maximal priority was removed */
      assert(vartype == SCIP_VARTYPE_BINARY || vartype == SCIP_VARTYPE_INTEGER);
      assert(branchpriority == branchcand->pseudomaxpriority);
      if( freepos != branchcand->npriopseudobins + branchcand->npriopseudoints - 1 )
      {
         branchcand->pseudocands[freepos] =
            branchcand->pseudocands[branchcand->npriopseudobins + branchcand->npriopseudoints - 1];
         branchcand->pseudocands[freepos]->pseudocandindex = freepos;
         freepos = branchcand->npriopseudobins + branchcand->npriopseudoints - 1;
      }
      branchcand->npriopseudoints--;
   }
   if( freepos < branchcand->npriopseudocands )
   {
      /* a candidate of maximal priority was removed */
      assert(branchpriority == branchcand->pseudomaxpriority);
      if( freepos != branchcand->npriopseudocands - 1 )
      {
         branchcand->pseudocands[freepos] = branchcand->pseudocands[branchcand->npriopseudocands - 1];
         branchcand->pseudocands[freepos]->pseudocandindex = freepos;
         freepos = branchcand->npriopseudocands - 1;
      }
      branchcand->npriopseudocands--;
   }
   if( freepos != branchcand->npseudocands - 1 )
   {
      branchcand->pseudocands[freepos] = branchcand->pseudocands[branchcand->npseudocands - 1];
      branchcand->pseudocands[freepos]->pseudocandindex = freepos;
   }
   branchcand->npseudocands--;

   assert(0 <= branchcand->npriopseudocands && branchcand->npriopseudocands <= branchcand->npseudocands);
   assert(0 <= branchcand->npriopseudobins && branchcand->npriopseudobins <= branchcand->npriopseudocands);
   assert(0 <= branchcand->npriopseudoints && branchcand->npriopseudoints <= branchcand->npriopseudocands);

   /* if all maximal priority candidates were removed, resort the array s.t. the new maximal priority candidates
    * are at the front
    */
   if( branchcand->npriopseudocands == 0 )
      branchcandSortPseudoCands(branchcand);
}

/** updates branching candidate list for a given variable */
RETCODE SCIPbranchcandUpdateVar(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   SET*             set,                /**< global SCIP settings */
   VAR*             var                 /**< variable that changed its bounds */
   )
{
   assert(branchcand != NULL);
   assert(var != NULL);
   assert(SCIPsetIsFeasLE(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
   
   if( (SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN)
      && SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS
      && SCIPsetIsLT(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      /* variable is neither continuous nor fixed: make sure it is member of the pseudo branching candidate list */
      if( var->pseudocandindex == -1 )
      {
         CHECK_OKAY( ensurePseudocandsSize(branchcand, set, branchcand->npseudocands+1) );

         branchcand->npseudocands++;
         branchcandInsertPseudoCand(branchcand, var, branchcand->npseudocands-1);
      }
   }
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR
         || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED
         || SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
         || SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

      /* variable is continuous or fixed: make sure it is not member of the pseudo branching candidate list */
      if( var->pseudocandindex >= 0 )
      {
         branchcandRemovePseudoCand(branchcand, var);
      }
   }

   return SCIP_OKAY;
}




/*
 * branching rule methods
 */

/** compares two branching rules w. r. to their priority */
DECL_SORTPTRCOMP(SCIPbranchruleComp)
{  /*lint --e{715}*/
   return ((BRANCHRULE*)elem2)->priority - ((BRANCHRULE*)elem1)->priority;
}

/** method to call, when the priority of a branching rule was changed */
static
DECL_PARAMCHGD(paramChgdBranchrulePriority)
{  /*lint --e{715}*/
   PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetBranchrulePriority() to mark the branchrules unsorted */
   CHECK_OKAY( SCIPsetBranchrulePriority(scip, (BRANCHRULE*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a branching rule */
RETCODE SCIPbranchruleCreate(
   BRANCHRULE**     branchrule,         /**< pointer to store branching rule */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   SET*             set,                /**< global SCIP settings */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   int              maxdepth,           /**< maximal depth level, up to which this branching rule should be used (or -1) */
   Real             maxbounddist,       /**< maximal relative distance from current node's dual bound to primal bound
                                         *   compared to best node's dual bound for applying branching rule
                                         *   (0.0: only on current best node, 1.0: on all nodes) */
   DECL_BRANCHFREE  ((*branchfree)),    /**< destructor of branching rule */
   DECL_BRANCHINIT  ((*branchinit)),    /**< initialize branching rule */
   DECL_BRANCHEXIT  ((*branchexit)),    /**< deinitialize branching rule */
   DECL_BRANCHEXECLP((*branchexeclp)),  /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXECPS((*branchexecps)),  /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
   char paramname[MAXSTRLEN];
   char paramdesc[MAXSTRLEN];

   assert(branchrule != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   ALLOC_OKAY( allocMemory(branchrule) );
   ALLOC_OKAY( duplicateMemoryArray(&(*branchrule)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*branchrule)->desc, desc, strlen(desc)+1) );
   (*branchrule)->priority = priority;
   (*branchrule)->maxdepth = maxdepth;
   (*branchrule)->maxbounddist = maxbounddist;
   (*branchrule)->branchfree = branchfree;
   (*branchrule)->branchinit = branchinit;
   (*branchrule)->branchexit = branchexit;
   (*branchrule)->branchexeclp = branchexeclp;
   (*branchrule)->branchexecps = branchexecps;
   (*branchrule)->branchruledata = branchruledata;
   CHECK_OKAY( SCIPclockCreate(&(*branchrule)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*branchrule)->nlpcalls = 0;
   (*branchrule)->npseudocalls = 0;
   (*branchrule)->ncutoffs = 0;
   (*branchrule)->ncutsfound = 0;
   (*branchrule)->nconssfound = 0;
   (*branchrule)->ndomredsfound = 0;
   (*branchrule)->nchildren = 0;
   (*branchrule)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "branching/%s/priority", name);
   sprintf(paramdesc, "priority of branching rule <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
         &(*branchrule)->priority, priority, INT_MIN, INT_MAX, 
         paramChgdBranchrulePriority, (PARAMDATA*)(*branchrule)) ); /*lint !e740*/
   sprintf(paramname, "branching/%s/maxdepth", name);
   sprintf(paramdesc, "maximal depth level, up to which branching rule <%s> should be used (-1 for no limit)", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
         &(*branchrule)->maxdepth, maxdepth, -1, INT_MAX, 
         NULL, NULL) ); /*lint !e740*/
   sprintf(paramname, "branching/%s/maxbounddist", name);
   sprintf(paramdesc, "maximal relative distance from current node's dual bound to primal bound compared to best node's dual bound for applying branching rule (0.0: only on current best node, 1.0: on all nodes)");
   CHECK_OKAY( SCIPsetAddRealParam(set, memhdr, paramname, paramdesc,
         &(*branchrule)->maxbounddist, maxbounddist, 0.0, 1.0, 
         NULL, NULL) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** frees memory of branching rule */   
RETCODE SCIPbranchruleFree(
   BRANCHRULE**     branchrule,         /**< pointer to branching rule data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(branchrule != NULL);
   assert(*branchrule != NULL);
   assert(!(*branchrule)->initialized);

   /* call destructor of branching rule */
   if( (*branchrule)->branchfree != NULL )
   {
      CHECK_OKAY( (*branchrule)->branchfree(scip, *branchrule) );
   }

   SCIPclockFree(&(*branchrule)->clock);
   freeMemoryArray(&(*branchrule)->name);
   freeMemoryArray(&(*branchrule)->desc);
   freeMemory(branchrule);

   return SCIP_OKAY;
}

/** initializes branching rule */
RETCODE SCIPbranchruleInit(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(branchrule != NULL);
   assert(scip != NULL);

   if( branchrule->initialized )
   {
      errorMessage("branching rule <%s> already initialized\n", branchrule->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(branchrule->clock);

   branchrule->nlpcalls = 0;
   branchrule->npseudocalls = 0;
   branchrule->ncutoffs = 0;
   branchrule->ncutsfound = 0;
   branchrule->nconssfound = 0;
   branchrule->ndomredsfound = 0;
   branchrule->nchildren = 0;

   if( branchrule->branchinit != NULL )
   {
      CHECK_OKAY( branchrule->branchinit(scip, branchrule) );
   }
   branchrule->initialized = TRUE;

   return SCIP_OKAY;
}

/** deinitializes branching rule */
RETCODE SCIPbranchruleExit(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(branchrule != NULL);
   assert(scip != NULL);

   if( !branchrule->initialized )
   {
      errorMessage("branching rule <%s> not initialized\n", branchrule->name);
      return SCIP_INVALIDCALL;
   }

   if( branchrule->branchexit != NULL )
   {
      CHECK_OKAY( branchrule->branchexit(scip, branchrule) );
   }
   branchrule->initialized = FALSE;

   return SCIP_OKAY;
}

/** executes branching rule for fractional LP solution */
RETCODE SCIPbranchruleExecLPSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   SEPASTORE*       sepastore,          /**< separation storage */
   Real             upperbound,         /**< global upper bound */
   Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->focusnode != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexeclp != NULL
      && (branchrule->maxdepth == -1 || branchrule->maxdepth >= SCIPtreeGetCurrentDepth(tree)) )
   {
      Real loclowerbound;
      Real glblowerbound;

      loclowerbound = SCIPnodeGetLowerbound(tree->focusnode);
      glblowerbound = SCIPtreeGetLowerbound(tree, set);
      if( SCIPsetIsLE(set, loclowerbound - glblowerbound, branchrule->maxbounddist * (upperbound - glblowerbound)) )
      {
         Longint oldndomchgs;
         int oldncutsstored;
         int oldnactiveconss;

         debugMessage("executing LP branching rule <%s>\n", branchrule->name);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldncutsstored = SCIPsepastoreGetNCutsStored(sepastore);
         oldnactiveconss = stat->nactiveconss;

         /* start timing */
         SCIPclockStart(branchrule->clock, set);
   
         /* call external method */
         CHECK_OKAY( branchrule->branchexeclp(set->scip, branchrule, allowaddcons, result) );

         /* stop timing */
         SCIPclockStop(branchrule->clock, set);
      
         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_SEPARATED
            && *result != SCIP_BRANCHED
            && *result != SCIP_DIDNOTRUN )
         {
            errorMessage("branching rule <%s> returned invalid result code <%d> from LP solution branching\n",
               branchrule->name, *result);
            return SCIP_INVALIDRESULT;
         }
         if( *result == SCIP_CONSADDED && !allowaddcons )
         {
            errorMessage("branching rule <%s> added a constraint in LP solution branching without permission\n",
               branchrule->name);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            branchrule->nlpcalls++;
         if( *result == SCIP_CUTOFF )
            branchrule->ncutoffs++;
         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);
            branchrule->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            branchrule->ncutsfound += SCIPsepastoreGetNCutsStored(sepastore) - oldncutsstored;
            branchrule->nconssfound += stat->nactiveconss - oldnactiveconss;
         }
         else
            branchrule->nchildren += tree->nchildren;
      }
   }

   return SCIP_OKAY;
}

/** executes branching rule for not completely fixed pseudo solution */
RETCODE SCIPbranchruleExecPseudoSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   Real             upperbound,         /**< global upper bound */
   Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   assert(tree != NULL);
   assert(tree->nchildren == 0);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexecps != NULL
      && (branchrule->maxdepth == -1 || branchrule->maxdepth >= SCIPtreeGetCurrentDepth(tree)) )
   {
      Real loclowerbound;
      Real glblowerbound;

      loclowerbound = SCIPnodeGetLowerbound(tree->focusnode);
      glblowerbound = SCIPtreeGetLowerbound(tree, set);
      if( SCIPsetIsLE(set, loclowerbound - glblowerbound, branchrule->maxbounddist * (upperbound - glblowerbound)) )
      {
         Longint oldndomchgs;
         Longint oldnactiveconss;

         debugMessage("executing pseudo branching rule <%s>\n", branchrule->name);

         oldndomchgs = stat->nboundchgs + stat->nholechgs;
         oldnactiveconss = stat->nactiveconss;

         /* start timing */
         SCIPclockStart(branchrule->clock, set);
   
         /* call external method */
         CHECK_OKAY( branchrule->branchexecps(set->scip, branchrule, allowaddcons, result) );

         /* stop timing */
         SCIPclockStop(branchrule->clock, set);
      
         /* evaluate result */
         if( *result != SCIP_CUTOFF
            && *result != SCIP_CONSADDED
            && *result != SCIP_REDUCEDDOM
            && *result != SCIP_BRANCHED
            && *result != SCIP_DIDNOTRUN )
         {
            errorMessage("branching rule <%s> returned invalid result code <%d> from pseudo solution branching\n",
               branchrule->name, *result);
            return SCIP_INVALIDRESULT;
         }
         if( *result == SCIP_CONSADDED && !allowaddcons )
         {
            errorMessage("branching rule <%s> added a constraint in pseudo solution branching without permission\n",
               branchrule->name);
            return SCIP_INVALIDRESULT;
         }

         /* update statistics */
         if( *result != SCIP_DIDNOTRUN )
            branchrule->npseudocalls++;
         if( *result == SCIP_CUTOFF )
            branchrule->ncutoffs++;
         if( *result != SCIP_BRANCHED )
         {
            assert(tree->nchildren == 0);
            branchrule->ndomredsfound += stat->nboundchgs + stat->nholechgs - oldndomchgs;
            branchrule->nconssfound += stat->nactiveconss - oldnactiveconss;
         }
         else
            branchrule->nchildren += tree->nchildren;
      }
   }

   return SCIP_OKAY;
}

/** gets user data of branching rule */
BRANCHRULEDATA* SCIPbranchruleGetData(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->branchruledata;
}

/** sets user data of branching rule; user has to free old data in advance! */
void SCIPbranchruleSetData(
   BRANCHRULE*      branchrule,         /**< branching rule */
   BRANCHRULEDATA*  branchruledata      /**< new branching rule user data */
   )
{
   assert(branchrule != NULL);

   branchrule->branchruledata = branchruledata;
}

/** gets name of branching rule */
const char* SCIPbranchruleGetName(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->name;
}

/** gets description of branching rule */
const char* SCIPbranchruleGetDesc(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->desc;
}

/** gets priority of branching rule */
int SCIPbranchruleGetPriority(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->priority;
}

/** sets priority of branching rule */
void SCIPbranchruleSetPriority(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the branching rule */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   
   branchrule->priority = priority;
   set->branchrulessorted = FALSE;
}

/** gets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
int SCIPbranchruleGetMaxdepth(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->maxdepth;
}

/** sets maximal depth level, up to which this branching rule should be used (-1 for no limit) */
void SCIPbranchruleSetMaxdepth(
   BRANCHRULE*      branchrule,         /**< branching rule */
   int              maxdepth            /**< new maxdepth of the branching rule */
   )
{
   assert(branchrule != NULL);
   assert(maxdepth >= -1);

   branchrule->maxdepth = maxdepth;
}

/** gets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
Real SCIPbranchruleGetMaxbounddist(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->maxbounddist;
}

/** sets maximal relative distance from current node's dual bound to primal bound for applying branching rule */
void SCIPbranchruleSetMaxbounddist(
   BRANCHRULE*      branchrule,         /**< branching rule */
   Real             maxbounddist        /**< new maxbounddist of the branching rule */
   )
{
   assert(branchrule != NULL);
   assert(maxbounddist >= -1);

   branchrule->maxbounddist = maxbounddist;
}

/** gets time in seconds used in this branching rule */
Real SCIPbranchruleGetTime(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return SCIPclockGetTime(branchrule->clock);
}

/** gets the total number of times, the branching rule was called on an LP solution */
Longint SCIPbranchruleGetNLPCalls(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->nlpcalls;
}

/** gets the total number of times, the branching rule was called on a pseudo solution */
Longint SCIPbranchruleGetNPseudoCalls(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->npseudocalls;
}

/** gets the total number of times, the branching rule detected a cutoff */
Longint SCIPbranchruleGetNCutoffs(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->ncutoffs;
}

/** gets the total number of cuts, the branching rule separated */
Longint SCIPbranchruleGetNCutsFound(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->ncutsfound;
}

/** gets the total number of constraints, the branching rule added to the respective local nodes (not counting constraints
 *  that were added to the child nodes as branching decisions)
 */
Longint SCIPbranchruleGetNConssFound(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->nconssfound;
}

/** gets the total number of domain reductions, the branching rule found */
Longint SCIPbranchruleGetNDomredsFound(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->ndomredsfound;
}

/** gets the total number of children, the branching rule created */
Longint SCIPbranchruleGetNChildren(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->nchildren;
}

/** is branching rule initialized? */
Bool SCIPbranchruleIsInitialized(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->initialized;
}




/*
 * branching methods
 */

/** calculates the branching score out of the gain predictions for a binary branching */
Real SCIPbranchGetScore(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   Real             downgain,           /**< prediction of objective gain for rounding downwards */
   Real             upgain              /**< prediction of objective gain for rounding upwards */
   )
{
   Real score;

   assert(set != NULL);

   /* weigh the two child nodes with branchscorefac and 1-branchscorefac */
   if( downgain > upgain )
      score = set->branch_scorefac * downgain + (1.0-set->branch_scorefac) * upgain;
   else
      score = set->branch_scorefac * upgain + (1.0-set->branch_scorefac) * downgain;

   /* slightly increase gains, such that for zero gains, the branch factor comes into account */
   score += SCIPsetSumepsilon(set);

   /* apply the branch factor of the variable */
   if( var != NULL )
      score *= SCIPvarGetBranchFactor(var);

   return score;
}

/** calculates the branching score out of the gain predictions for a branching with arbitrary many children */
Real SCIPbranchGetScoreMultiple(
   SET*             set,                /**< global SCIP settings */
   VAR*             var,                /**< variable, of which the branching factor should be applied, or NULL */
   int              nchildren,          /**< number of children that the branching will create */
   Real*            gains               /**< prediction of objective gain for each child */
   )
{
   Real min1;
   Real min2;
   int c;

   assert(nchildren == 0 || gains != NULL);

   /* search for the two minimal gains in the child list and use these to calculate the branching score */
   min1 = SCIPsetInfinity(set);
   min2 = SCIPsetInfinity(set);
   for( c = 0; c < nchildren; ++c )
   {
      if( gains[c] < min1 )
      {
         min2 = min1;
         min1 = gains[c];
      }
      else if( gains[c] < min2 )
         min2 = gains[c];
   }

   return SCIPbranchGetScore(set, var, min1, min2);
}

/** calls branching rules to branch on an LP solution; if no fractional variables exist, the result is SCIP_DIDNOTRUN;
 *  if the branch priority of an unfixed variable is larger than the maximal branch priority of the fractional
 *  variables, pseudo solution branching is applied on the unfixed variables with maximal branch priority
 */
RETCODE SCIPbranchExecLP(
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   SEPASTORE*       sepastore,          /**< separation storage */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             upperbound,         /**< global upper bound */
   Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   int i;

   assert(branchcand != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   /* calculate branching candidates */
   CHECK_OKAY( branchcandCalcLPCands(branchcand, set, stat, lp) );
   assert(0 <= branchcand->npriolpcands && branchcand->npriolpcands <= branchcand->nlpcands);
   assert((branchcand->npriolpcands == 0) == (branchcand->nlpcands == 0));

   debugMessage("branching on LP solution with %d fractional variables (%d of maximal priority)\n",
      branchcand->nlpcands, branchcand->npriolpcands);

   /* do nothing, if no fractional variables exist */
   if( branchcand->nlpcands == 0 )
      return SCIP_OKAY;

   /* if there is a non-fixed variable with higher priority than the maximal priority of the fractional candidates,
    * use pseudo solution branching instead
    */
   if( branchcand->pseudomaxpriority > branchcand->lpmaxpriority )
   {
      CHECK_OKAY( SCIPbranchExecPseudo(memhdr, set, stat, tree, lp, branchcand, eventqueue, upperbound, allowaddcons,
            result) );
      assert(*result != SCIP_DIDNOTRUN);
      return SCIP_OKAY;
   }

   /* sort the branching rules by priority */
   SCIPsetSortBranchrules(set);

   /* try all branching rules until one succeeded to branch */
   for( i = 0; i < set->nbranchrules && *result == SCIP_DIDNOTRUN; ++i )
   {
      CHECK_OKAY( SCIPbranchruleExecLPSol(set->branchrules[i], set, stat, tree, sepastore, upperbound, allowaddcons, 
            result) );
   }

   if( *result == SCIP_DIDNOTRUN )
   {
      VAR* var;
      Real factor;
      Real bestfactor;
      int priority;
      int bestpriority;
      int bestcand;

      /* no branching method succeeded in choosing a branching: just branch on the first fractional variable with maximal
       * priority, and out of these on the one with maximal branch factor
       */
      bestcand = -1;
      bestpriority = INT_MIN;
      bestfactor = REAL_MIN;
      for( i = 0; i < branchcand->nlpcands; ++i )
      {
         priority = SCIPvarGetBranchPriority(branchcand->lpcands[i]);
         factor = SCIPvarGetBranchFactor(branchcand->lpcands[i]);
         if( priority > bestpriority || (priority == bestpriority && factor > bestfactor) )
         {
            bestcand = i;
            bestpriority = priority;
            bestfactor = factor;
         }
      }
      assert(0 <= bestcand && bestcand < branchcand->nlpcands);

      var = branchcand->lpcands[bestcand];
      assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
      assert(!SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

      CHECK_OKAY( SCIPtreeBranchVar(tree, memhdr, set, stat, lp, branchcand, eventqueue, var) );

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

/** calls branching rules to branch on a pseudo solution; if no unfixed variables exist, the result is SCIP_DIDNOTRUN */
RETCODE SCIPbranchExecPseudo(
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch and bound tree */
   LP*              lp,                 /**< current LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   Real             upperbound,         /**< global upper bound */
   Bool             allowaddcons,       /**< should adding constraints be allowed to avoid a branching? */
   RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   )
{
   int i;
   
   assert(branchcand != NULL);
   assert(result != NULL);

   debugMessage("branching on pseudo solution with %d unfixed variables\n", branchcand->npseudocands);

   *result = SCIP_DIDNOTRUN;

   /* do nothing, if no unfixed variables exist */
   if( branchcand->npseudocands == 0 )
      return SCIP_OKAY;

   /* sort the branching rules by priority */
   SCIPsetSortBranchrules(set);

   /* try all branching rules until one succeeded to branch */
   for( i = 0; i < set->nbranchrules && *result == SCIP_DIDNOTRUN; ++i )
   {
      CHECK_OKAY( SCIPbranchruleExecPseudoSol(set->branchrules[i], set, stat, tree, upperbound, allowaddcons, result) );
   }

   if( *result == SCIP_DIDNOTRUN )
   {
      VAR* var;
      Real factor;
      Real bestfactor;
      int priority;
      int bestpriority;
      int bestcand;

      /* no branching method succeeded in choosing a branching: just branch on the first unfixed variable with maximal
       * priority, and out of these on the one with maximal branch factor
       */
      bestcand = -1;
      bestpriority = INT_MIN;
      bestfactor = REAL_MIN;
      for( i = 0; i < branchcand->npseudocands; ++i )
      {
         priority = SCIPvarGetBranchPriority(branchcand->pseudocands[i]);
         factor = SCIPvarGetBranchFactor(branchcand->pseudocands[i]);
         if( priority > bestpriority || (priority == bestpriority && factor > bestfactor) )
         {
            bestcand = i;
            bestpriority = priority;
            bestfactor = factor;
         }
      }
      assert(0 <= bestcand && bestcand < branchcand->npseudocands);

      var = branchcand->pseudocands[bestcand];
      assert(SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS);
      assert(!SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));

      CHECK_OKAY( SCIPtreeBranchVar(tree, memhdr, set, stat, lp, branchcand, eventqueue, var) );

      *result = SCIP_BRANCHED;
   }

   return SCIP_OKAY;
}

