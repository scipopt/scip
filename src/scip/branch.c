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

/**@file   branch.c
 * @brief  datastructures and methods for branching methods
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "branch.h"


/** branching candidate storage */
struct BranchCand
{
   VAR**            lpcands;            /**< candidates for branching on LP solution (fractional integer variables) */
   Real*            lpcandsfrac;        /**< fractionalities of LP candidates */
   VAR**            pseudocands;        /**< candidates for branching on pseudo solution (non-fixed integer variables) */
   int              lpcandssize;        /**< number of available slots in lpcands array */
   int              nlpcands;           /**< number of candidates for branching on LP solution */
   int              pseudocandssize;    /**< number of available slots in pseudocands array */
   int              npseudocands;       /**< number of candidates for branching on pseudo solution */
   int              validlpcandslp;     /**< lp number for which lpcands are valid */
   int              validpseudocandsbc; /**< bound change number for which pseudocands are valid */
};


/** branching method data */
struct Branch
{
   char*            name;               /**< name of branching method */
   char*            desc;               /**< description of branching method */
   DECL_BRANCHFREE((*branchfree));      /**< destructor of branching method */
   DECL_BRANCHINIT((*branchinit));      /**< initialise branching method */
   DECL_BRANCHEXIT((*branchexit));      /**< deinitialise branching method */
   DECL_BRANCHEXLP((*branchexlp));      /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXPS((*branchexps));      /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHDATA*      branchdata;         /**< branching method data */
   unsigned int     initialized:1;      /**< is branching method initialized? */
};




/*
 * memory growing methods for dynamically allocated arrays
 */

static
RETCODE ensureLpcandsSize(              /**< ensures, that lpcands array can store at least num entries */
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
      ALLOC_OKAY( reallocMemoryArray(branchcand->lpcands, newsize) );
      ALLOC_OKAY( reallocMemoryArray(branchcand->lpcandsfrac, newsize) );
      branchcand->lpcandssize = newsize;
   }
   assert(num <= branchcand->lpcandssize);

   return SCIP_OKAY;
}

static
RETCODE ensurePseudocandsSize(          /**< ensures, that pseudocands array can store at least num entries */
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
      ALLOC_OKAY( reallocMemoryArray(branchcand->pseudocands, newsize) );
      branchcand->pseudocandssize = newsize;
   }
   assert(num <= branchcand->pseudocandssize);

   return SCIP_OKAY;
}



/*
 * branching candidate storage methods
 */

RETCODE SCIPbranchcandCreate(           /**< creates a branching candidate storage */
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   )
{
   assert(branchcand != NULL);

   ALLOC_OKAY( allocMemory(*branchcand) );
   (*branchcand)->lpcands = NULL;
   (*branchcand)->lpcandsfrac = NULL;
   (*branchcand)->pseudocands = NULL;
   (*branchcand)->lpcandssize = 0;
   (*branchcand)->nlpcands = 0;
   (*branchcand)->pseudocandssize = 0;
   (*branchcand)->npseudocands = 0;
   (*branchcand)->validlpcandslp = -1;
   (*branchcand)->validpseudocandsbc = -1;

   return SCIP_OKAY;
}

RETCODE SCIPbranchcandFree(             /**< frees branching candidate storage */
   BRANCHCAND**     branchcand          /**< pointer to store branching candidate storage */
   )
{
   assert(branchcand != NULL);

   freeMemoryArrayNull((*branchcand)->lpcands);
   freeMemoryArrayNull((*branchcand)->lpcandsfrac);
   freeMemoryArrayNull((*branchcand)->pseudocands);
   freeMemory(*branchcand);

   return SCIP_OKAY;
}

RETCODE SCIPbranchcandGetLPCands(       /**< gets branching candidates for LP solution branching (fractional variables) */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   VAR***           lpcands,            /**< pointer to store the array of LP branching candidates, or NULL */
   Real**           lpcandsfrac,        /**< pointer to store the array of LP candidate fractionalities, or NULL */
   int*             nlpcands            /**< pointer to store the number of LP branching candidates, or NULL */
   )
{
   assert(branchcand != NULL);
   assert(stat != NULL);
   assert(branchcand->validlpcandslp <= stat->nlp);
   assert(lp != NULL);
   assert(lp->solved);
   assert(lp->flushed);

   debugMessage("getting LP branching candidates: validlp=%d, nlp=%d\n", branchcand->validlpcandslp, stat->nlp);

   /* check, if the actual LP branching candidate array is invalid */
   if( branchcand->validlpcandslp < stat->nlp )
   {
      VAR* var;
      COL* col;
      Real frac;
      int c;

      debugMessage(" -> recalculating LP branching candidates\n");

      /* construct the LP branching candidate set */
      CHECK_OKAY( ensureLpcandsSize(branchcand, set, lp->ncols) );

      branchcand->nlpcands = 0;
      for( c = 0; c < lp->ncols; ++c )
      {
         col = lp->cols[c];
         assert(col != NULL);
         assert(col->primsol < SCIP_INVALID);
         assert(col->inlp);
         assert(col->lpipos >= 0);

         var = col->var;
         assert(var != NULL);
         assert(var->varstatus == SCIP_VARSTATUS_COLUMN);
         assert(var->data.col == col);
         
         /* LP branching candidates are fractional binary and integer variables */
         if( var->vartype == SCIP_VARTYPE_BINARY || var->vartype == SCIP_VARTYPE_INTEGER )
         {
            frac = SCIPsetFrac(set, col->primsol);
            if( !SCIPsetIsFracIntegral(set, frac) )
            {
               debugMessage(" -> candidate %d: var=<%s>, sol=%g, frac=%g\n", 
                  branchcand->nlpcands, var->name, col->primsol, frac);

               assert(branchcand->nlpcands < branchcand->lpcandssize);
               branchcand->lpcands[branchcand->nlpcands] = var;
               branchcand->lpcandsfrac[branchcand->nlpcands] = frac;
               branchcand->nlpcands++;
            }
         }
      }

      branchcand->validlpcandslp = stat->nlp;
   }

   debugMessage(" -> %d fractional variables\n", branchcand->nlpcands);

   /* assign return values */
   if( lpcands != NULL )
      *lpcands = branchcand->lpcands;
   if( lpcandsfrac != NULL )
      *lpcandsfrac = branchcand->lpcandsfrac;
   if( nlpcands != NULL )
      *nlpcands = branchcand->nlpcands;

   return SCIP_OKAY;
}

RETCODE SCIPbranchcandGetPseudoCands(   /**< gets branching candidates for pseudo solution branching (nonfixed variables) */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB*            prob,               /**< problem data */
   VAR***           pseudocands,        /**< pointer to store the array of pseudo branching candidates, or NULL */
   int*             npseudocands        /**< pointer to store the number of pseudo branching candidates, or NULL */
   )
{
   assert(branchcand != NULL);
   assert(stat != NULL);
   assert(prob != NULL);

   /* check, if the actual pseudo branching candidate array is invalid */
   if( branchcand->validpseudocandsbc != stat->nboundchanges )
   {
      VAR* var;
      int v;

      /* construct the pseudo branching candidate set */
      CHECK_OKAY( ensurePseudocandsSize(branchcand, set, prob->nbin + prob->nint + prob->nimpl) );

      branchcand->npseudocands = 0;
      for( v = 0; v < prob->nvars; ++v )
      {
         var = prob->vars[v];
         assert(var != NULL);
         assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);

         /* pseudo branching candidates are non-fixed binary, integer, and implicit integer variables */
         if( var->vartype == SCIP_VARTYPE_BINARY
            || var->vartype == SCIP_VARTYPE_INTEGER
            || var->vartype == SCIP_VARTYPE_IMPLINT )
         {
            assert(SCIPsetIsIntegral(set, var->dom.lb));
            assert(SCIPsetIsIntegral(set, var->dom.ub));
            
            if( !SCIPsetIsFixed(set, var->dom.lb, var->dom.ub) )
            {
               assert(branchcand->npseudocands < branchcand->pseudocandssize);
               branchcand->pseudocands[branchcand->npseudocands] = var;
               branchcand->npseudocands++;
            }
         }
      }

      branchcand->validpseudocandsbc = stat->nboundchanges;
   }

   /* assign return values */
   if( pseudocands != NULL )
      *pseudocands = branchcand->pseudocands;
   if( npseudocands != NULL )
      *npseudocands = branchcand->npseudocands;

   return SCIP_OKAY;
}



/*
 * branching methods
 */

RETCODE SCIPbranchCreate(               /**< creates a branching method */
   BRANCH**         branch,             /**< pointer to store branching method */
   const char*      name,               /**< name of branching method */
   const char*      desc,               /**< description of branching method */
   DECL_BRANCHFREE((*branchfree)),      /**< destructor of branching method */
   DECL_BRANCHINIT((*branchinit)),      /**< initialise branching method */
   DECL_BRANCHEXIT((*branchexit)),      /**< deinitialise branching method */
   DECL_BRANCHEXLP((*branchexlp)),      /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXPS((*branchexps)),      /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHDATA*      branchdata          /**< branching method data */
   )
{
   assert(branch != NULL);
   assert(name != NULL);
   assert(desc != NULL);

   ALLOC_OKAY( allocMemory(*branch) );
   ALLOC_OKAY( duplicateMemoryArray((*branch)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray((*branch)->desc, desc, strlen(desc)+1) );
   (*branch)->branchfree = branchfree;
   (*branch)->branchinit = branchinit;
   (*branch)->branchexit = branchexit;
   (*branch)->branchexlp = branchexlp;
   (*branch)->branchexps = branchexps;
   (*branch)->branchdata = branchdata;
   (*branch)->initialized = FALSE;

   return SCIP_OKAY;
}
   
RETCODE SCIPbranchFree(                 /**< frees memory of branching method */
   BRANCH**         branch,             /**< pointer to branching method data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(branch != NULL);
   assert(*branch != NULL);
   assert(!(*branch)->initialized);

   /* call destructor of branching method */
   if( (*branch)->branchfree != NULL )
   {
      CHECK_OKAY( (*branch)->branchfree(*branch, scip) );
   }

   freeMemoryArray((*branch)->name);
   freeMemoryArray((*branch)->desc);
   freeMemory(*branch);

   return SCIP_OKAY;
}

RETCODE SCIPbranchInit(                 /**< initializes branching method */
   BRANCH*          branch,             /**< branching method */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(branch != NULL);
   assert(scip != NULL);

   if( branch->initialized )
   {
      char s[255];
      sprintf(s, "branching method <%s> already initialized", branch->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( branch->branchinit != NULL )
   {
      CHECK_OKAY( branch->branchinit(branch, scip) );
   }
   branch->initialized = TRUE;

   return SCIP_OKAY;
}

RETCODE SCIPbranchExit(                 /**< deinitializes branching method */
   BRANCH*          branch,             /**< branching method */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(branch != NULL);
   assert(scip != NULL);

   if( !branch->initialized )
   {
      char s[255];
      sprintf(s, "branching method <%s> not initialized", branch->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( branch->branchexit != NULL )
   {
      CHECK_OKAY( branch->branchexit(branch, scip) );
   }
   branch->initialized = FALSE;

   return SCIP_OKAY;
}

RETCODE SCIPbranchExecLPSol(            /**< executes branching method for fractional LP solution */
   BRANCH*          branch,             /**< branching method */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branch != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;
   if( branch->branchexlp != NULL )
   {
      CHECK_OKAY( branch->branchexlp(branch, scip, result) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPbranchExecPseudoSol(        /**< executes branching method for not completely fixed pseudo solution */
   BRANCH*          branch,             /**< branching method */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branch != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_FEASIBLE;
   if( branch->branchexps != NULL )
   {
      CHECK_OKAY( branch->branchexps(branch, scip, result) );
   }

   return SCIP_OKAY;
}

const char* SCIPbranchGetName(          /**< gets name of branching method */
   BRANCH*          branch              /**< branching method */
   )
{
   assert(branch != NULL);

   return branch->name;
}

BRANCHDATA* SCIPbranchGetData(          /**< gets user data of branching method */
   BRANCH*          branch              /**< branching method */
   )
{
   assert(branch != NULL);

   return branch->branchdata;
}

void SCIPbranchSetData(                 /**< sets user data of branching method; user has to free old data in advance! */
   BRANCH*          branch,             /**< branching method */
   BRANCHDATA*      branchdata          /**< new branching method user data */
   )
{
   assert(branch != NULL);

   branch->branchdata = branchdata;
}

Bool SCIPbranchIsInitialized(           /**< is branching method initialized? */
   BRANCH*          branch              /**< branching method */
   )
{
   assert(branch != NULL);

   return branch->initialized;
}

