/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: branch.c,v 1.31 2003/12/03 18:08:12 bzfpfend Exp $"

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
#include "paramset.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
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
   const SET*       set,                /**< global SCIP settings */
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
   const SET*       set,                /**< global SCIP settings */
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
   (*branchcand)->pseudocandssize = 0;
   (*branchcand)->npseudocands = 0;
   (*branchcand)->npseudobins = 0;
   (*branchcand)->npseudoints = 0;
   (*branchcand)->npseudoimpls = 0;
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

/** gets branching candidates for LP solution branching (fractional variables) */
RETCODE SCIPbranchcandGetLPCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandssol,         /**< pointer to store the array of LP candidate solution values, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands            /**< pointer to store the number of LP branching candidates, or NULL */
   )
{
   assert(branchcand != NULL);
   assert(stat != NULL);
   assert(branchcand->validlpcandslp <= stat->lpcount);
   assert(lp != NULL);
   assert(lp->solved);
   assert(lp->flushed);
   assert(lp->lpsolstat == SCIP_LPSOLSTAT_OPTIMAL || lp->lpsolstat == SCIP_LPSOLSTAT_UNBOUNDED);

   debugMessage("getting LP branching candidates: validlp=%d, lpcount=%d\n", branchcand->validlpcandslp, stat->lpcount);

   /* check, if the actual LP branching candidate array is invalid */
   if( branchcand->validlpcandslp < stat->lpcount )
   {
      COL** cols;
      VAR* var;
      COL* col;
      Real frac;
      int ncols;
      int c;

      debugMessage(" -> recalculating LP branching candidates\n");

      cols = SCIPlpGetCols(lp);
      ncols = SCIPlpGetNCols(lp);

      /* construct the LP branching candidate set */
      CHECK_OKAY( ensureLpcandsSize(branchcand, set, ncols) );

      branchcand->nlpcands = 0;
      for( c = 0; c < ncols; ++c )
      {
         col = cols[c];
         assert(col != NULL);
         assert(col->primsol < SCIP_INVALID);
         assert(col->lppos == c);
         assert(col->lpipos >= 0);

         var = col->var;
         assert(var != NULL);
         assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);
         assert(SCIPvarGetCol(var) == col);
         
         /* LP branching candidates are fractional binary and integer variables */
         if( SCIPvarGetType(var) == SCIP_VARTYPE_BINARY || SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
         {
            frac = SCIPsetFrac(set, col->primsol);
            if( !SCIPsetIsFracIntegral(set, frac) )
            {
               debugMessage(" -> candidate %d: var=<%s>, sol=%g, frac=%g\n", 
                  branchcand->nlpcands, SCIPvarGetName(var), col->primsol, frac);

               assert(branchcand->nlpcands < branchcand->lpcandssize);
               branchcand->lpcands[branchcand->nlpcands] = var;
               branchcand->lpcandssol[branchcand->nlpcands] = col->primsol;
               branchcand->lpcandsfrac[branchcand->nlpcands] = frac;
               branchcand->nlpcands++;
            }
         }
      }

      branchcand->validlpcandslp = stat->lpcount;
   }

   debugMessage(" -> %d fractional variables\n", branchcand->nlpcands);

   /* assign return values */
   if( lpcands != NULL )
      *lpcands = branchcand->lpcands;
   if( lpcandssol != NULL )
      *lpcandssol = branchcand->lpcandssol;
   if( lpcandsfrac != NULL )
      *lpcandsfrac = branchcand->lpcandsfrac;
   if( nlpcands != NULL )
      *nlpcands = branchcand->nlpcands;

   return SCIP_OKAY;
}

/** gets branching candidates for pseudo solution branching (nonfixed variables) */
RETCODE SCIPbranchcandGetPseudoCands(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< problem data */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands        /**< pointer to store the number of pseudo branching candidates, or NULL */
   )
{
   assert(branchcand != NULL);

#ifndef NDEBUG
   /* check, if the actual pseudo branching candidate array is correct */
   {
      VAR* var;
      int npcs;
      int v;
      
      assert(prob != NULL);
      
      /* pseudo branching candidates are non-fixed binary, integer, and implicit integer variables */
      npcs = 0;
      for( v = 0; v < prob->nbin + prob->nint + prob->nimpl; ++v )
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

   return SCIP_OKAY;
}

/** updates branching candidate list for a given variable */
RETCODE SCIPbranchcandUpdateVar(
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< variable that changed its bounds */
   )
{
   assert(branchcand != NULL);
   assert(var != NULL);
   assert(SCIPsetIsFeasLE(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)));
   
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_FIXED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_AGGREGATED
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_MULTAGGR
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_NEGATED
      || SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
      || SCIPsetIsEQ(set, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) )
   {
      /* variable is continuous or fixed: make sure it is not member of the pseudo branching candidate list */
      if( var->pseudocandindex >= 0 )
      {
         int freepos;
         int intstart;
         int implstart;

         assert(var->pseudocandindex < branchcand->npseudocands);
         assert(branchcand->pseudocands[branchcand->npseudocands-1] != NULL);

         debugMessage("deleting pseudo candidate <%s> of type %d at %d from candidate set (%d/%d/%d)\n",
            SCIPvarGetName(var), SCIPvarGetType(var), var->pseudocandindex, 
            branchcand->npseudobins, branchcand->npseudoints, branchcand->npseudoimpls);

         /* delete the variable from pseudocands, retaining the ordering binaries, integers, implicit integers */
         intstart = branchcand->npseudobins;
         implstart = intstart + branchcand->npseudoints;
         freepos = var->pseudocandindex;
         assert(0 <= freepos && freepos < branchcand->npseudocands);

         if( freepos < intstart )
            branchcand->npseudobins--;
         else if( freepos < implstart )
            branchcand->npseudoints--;
         else
            branchcand->npseudoimpls--;

         if( freepos < intstart-1 )
         {
            /* move last binary to free slot */
            branchcand->pseudocands[freepos] = branchcand->pseudocands[intstart-1];
            branchcand->pseudocands[freepos]->pseudocandindex = freepos;
            freepos = intstart-1;
         }
         if( freepos < implstart-1 )
         {
            /* move last integer to free slot */
            branchcand->pseudocands[freepos] = branchcand->pseudocands[implstart-1];
            branchcand->pseudocands[freepos]->pseudocandindex = freepos;
            freepos = implstart-1;
         }
         if( freepos < branchcand->npseudocands-1 )
         {
            /* move last integer to free slot */
            branchcand->pseudocands[freepos] = branchcand->pseudocands[branchcand->npseudocands-1];
            branchcand->pseudocands[freepos]->pseudocandindex = freepos;
            freepos = branchcand->npseudocands-1;
         }
         assert(freepos == branchcand->npseudocands-1);

         branchcand->npseudocands--;
         var->pseudocandindex = -1;

         assert(branchcand->npseudocands == branchcand->npseudobins + branchcand->npseudoints + branchcand->npseudoimpls);
      }
   }
   else
   {
      assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

      /* variable is not fixed: make sure it is member of the pseudo branching candidate list */
      if( var->pseudocandindex == -1 )
      {
         int insertpos;
         int intstart;
         int implstart;

         debugMessage("adding pseudo candidate <%s> of type %d to candidate set (%d/%d/%d)\n",
            SCIPvarGetName(var), SCIPvarGetType(var),
            branchcand->npseudobins, branchcand->npseudoints, branchcand->npseudoimpls);

         CHECK_OKAY( ensurePseudocandsSize(branchcand, set, branchcand->npseudocands+1) );

         /* insert the variable into pseudocands, retaining the ordering binaries, integers, implicit integers */
         intstart = branchcand->npseudobins;
         implstart = intstart + branchcand->npseudoints;
         insertpos = branchcand->npseudocands;
         if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
            branchcand->npseudoimpls++;
         else
         {
            if( insertpos > implstart )
            {
               branchcand->pseudocands[insertpos] = branchcand->pseudocands[implstart];
               branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
               insertpos = implstart;
            }
            assert(insertpos == implstart);

            if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
               branchcand->npseudoints++;
            else
            {
               assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
               if( insertpos > intstart )
               {
                  branchcand->pseudocands[insertpos] = branchcand->pseudocands[intstart];
                  branchcand->pseudocands[insertpos]->pseudocandindex = insertpos;
                  insertpos = intstart;
               }
               assert(insertpos == intstart);

               branchcand->npseudobins++;
            }
         }
         branchcand->npseudocands++;

         assert(branchcand->npseudocands == branchcand->npseudobins + branchcand->npseudoints + branchcand->npseudoimpls);
         assert((SCIPvarGetType(var) == SCIP_VARTYPE_BINARY && insertpos == branchcand->npseudobins - 1)
            || (SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER
               && insertpos == branchcand->npseudobins + branchcand->npseudoints - 1)
            || (SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT && insertpos == branchcand->npseudocands - 1));
         
         branchcand->pseudocands[insertpos] = var;
         var->pseudocandindex = insertpos;
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
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
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
   (*branchrule)->branchfree = branchfree;
   (*branchrule)->branchinit = branchinit;
   (*branchrule)->branchexit = branchexit;
   (*branchrule)->branchexeclp = branchexeclp;
   (*branchrule)->branchexecps = branchexecps;
   (*branchrule)->branchruledata = branchruledata;
   (*branchrule)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "branching/%s/priority", name);
   sprintf(paramdesc, "priority of branching rule <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
                  &(*branchrule)->priority, priority, INT_MIN, INT_MAX, 
                  paramChgdBranchrulePriority, (PARAMDATA*)(*branchrule)) ); /*lint !e740*/

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
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexeclp != NULL )
   {
      CHECK_OKAY( branchrule->branchexeclp(set->scip, branchrule, result) );
      if( *result != SCIP_CUTOFF
         && *result != SCIP_BRANCHED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_DIDNOTRUN )
      {
         errorMessage("branching rule <%s> returned invalid result code <%d> from LP solution branching\n",
            branchrule->name, *result);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** executes branching rule for not completely fixed pseudo solution */
RETCODE SCIPbranchruleExecPseudoSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   const SET*       set,                /**< global SCIP settings */
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(set != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexecps != NULL )
   {
      CHECK_OKAY( branchrule->branchexecps(set->scip, branchrule, result) );
      if( *result != SCIP_CUTOFF
         && *result != SCIP_BRANCHED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_DIDNOTRUN )
      {
         errorMessage("branching rule <%s> returned invalid result code <%d> from LP solution branching\n",
            branchrule->name, *result);
         return SCIP_INVALIDRESULT;
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

/** calculates the branching score out of the downward and upward gain prediction */
Real SCIPbranchGetScore(
   const SET*       set,                /**< global SCIP settings */
   Real             downgain,           /**< prediction of objective gain for branching downwards */
   Real             upgain              /**< prediction of objective gain for branching upwards */
   )
{
   assert(set != NULL);

   if( downgain < upgain )
      return set->branchscorefac * downgain + (1.0-set->branchscorefac) * upgain;
   else
      return set->branchscorefac * upgain + (1.0-set->branchscorefac) * downgain;
}
