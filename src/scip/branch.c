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
   Real*            lpcandssol;         /**< solution values of LP candidates */
   Real*            lpcandsfrac;        /**< fractionalities of LP candidates */
   VAR**            pseudocands;        /**< candidates for branching on pseudo solution (non-fixed integer variables) */
   int              lpcandssize;        /**< number of available slots in lpcands array */
   int              nlpcands;           /**< number of candidates for branching on LP solution */
   int              pseudocandssize;    /**< number of available slots in pseudocands array */
   int              npseudocands;       /**< number of candidates for branching on pseudo solution */
   int              validlpcandslp;     /**< lp number for which lpcands are valid */
};


/** branching rule */
struct BranchRule
{
   char*            name;               /**< name of branching rule */
   char*            desc;               /**< description of branching rule */
   int              priority;           /**< priority of the branching rule */
   DECL_BRANCHFREE((*branchfree));      /**< destructor of branching rule */
   DECL_BRANCHINIT((*branchinit));      /**< initialise branching rule */
   DECL_BRANCHEXIT((*branchexit));      /**< deinitialise branching rule */
   DECL_BRANCHEXLP((*branchexlp));      /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXPS((*branchexps));      /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata;     /**< branching rule data */
   unsigned int     initialized:1;      /**< is branching rule initialized? */
};




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
   BRANCHCAND**     branchcand,         /**< pointer to store branching candidate storage */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< problem data */
   )
{
   int v;

   assert(branchcand != NULL);
   assert(prob != NULL);

   ALLOC_OKAY( allocMemory(branchcand) );
   (*branchcand)->lpcands = NULL;
   (*branchcand)->lpcandssol = NULL;
   (*branchcand)->lpcandsfrac = NULL;
   (*branchcand)->pseudocands = NULL;
   (*branchcand)->lpcandssize = 0;
   (*branchcand)->nlpcands = 0;
   (*branchcand)->pseudocandssize = 0;
   (*branchcand)->npseudocands = 0;
   (*branchcand)->validlpcandslp = -1;

   /* init pseudo branching candidate list */
   CHECK_OKAY( ensurePseudocandsSize(*branchcand, set, prob->nbin + prob->nint + prob->nimpl) );
   for( v = 0; v < prob->nbin + prob->nint + prob->nimpl; ++v )
   {
      CHECK_OKAY( SCIPbranchcandUpdateVar(*branchcand, set, prob->vars[v]) );
   }
   
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
         assert(col->lppos == c);
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
               branchcand->lpcandssol[branchcand->nlpcands] = col->primsol;
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
      int npseudocands;
      int v;
      
      assert(prob != NULL);
      
      /* pseudo branching candidates are non-fixed binary, integer, and implicit integer variables */
      npseudocands = 0;
      for( v = 0; v < prob->nbin + prob->nint + prob->nimpl; ++v )
      {
         var = prob->vars[v];
         assert(var != NULL);
         assert(var->varstatus == SCIP_VARSTATUS_LOOSE || var->varstatus == SCIP_VARSTATUS_COLUMN);
         assert(var->vartype == SCIP_VARTYPE_BINARY
            || var->vartype == SCIP_VARTYPE_INTEGER
            || var->vartype == SCIP_VARTYPE_IMPLINT);
         assert(SCIPsetIsIntegral(set, var->dom.lb));
         assert(SCIPsetIsIntegral(set, var->dom.ub));
         
         if( !SCIPsetIsFixed(set, var->dom.lb, var->dom.ub) )
         {
            assert(0 <= var->pseudocandindex && var->pseudocandindex < branchcand->npseudocands);
            assert(branchcand->pseudocands[var->pseudocandindex] == var);
            npseudocands++;
         }
         else
         {
            assert(var->pseudocandindex == -1);
         }
      }
      assert(branchcand->npseudocands == npseudocands);
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

   if( var->vartype == SCIP_VARTYPE_BINARY
      || var->vartype == SCIP_VARTYPE_INTEGER
      || var->vartype == SCIP_VARTYPE_IMPLINT )
   {
      if( SCIPsetIsFixed(set, var->dom.lb, var->dom.ub) )
      {
         /* variable is fixed: make sure it is not member of the pseudo branching candidate list */
         if( var->pseudocandindex >= 0 )
         {
            assert(var->pseudocandindex < branchcand->npseudocands);
            assert(branchcand->pseudocands[branchcand->npseudocands-1] != NULL);
            branchcand->pseudocands[var->pseudocandindex] = branchcand->pseudocands[branchcand->npseudocands-1];
            branchcand->pseudocands[var->pseudocandindex]->pseudocandindex = var->pseudocandindex;
            branchcand->npseudocands--;
            var->pseudocandindex = -1;
         }
      }
      else
      {
         /* variable is not fixed: make sure it is member of the pseudo branching candidate list */
         if( var->pseudocandindex == -1 )
         {
            CHECK_OKAY( ensurePseudocandsSize(branchcand, set, branchcand->npseudocands+1) );
            branchcand->pseudocands[branchcand->npseudocands] = var;
            var->pseudocandindex = branchcand->npseudocands;
            branchcand->npseudocands++;
         }
      }
   }

   return SCIP_OKAY;
}



/*
 * branching rule methods
 */

/** creates a branching rule */
RETCODE SCIPbranchruleCreate(
   BRANCHRULE**     branchrule,         /**< pointer to store branching rule */
   const char*      name,               /**< name of branching rule */
   const char*      desc,               /**< description of branching rule */
   int              priority,           /**< priority of the branching rule */
   DECL_BRANCHFREE((*branchfree)),      /**< destructor of branching rule */
   DECL_BRANCHINIT((*branchinit)),      /**< initialise branching rule */
   DECL_BRANCHEXIT((*branchexit)),      /**< deinitialise branching rule */
   DECL_BRANCHEXLP((*branchexlp)),      /**< branching execution method for fractional LP solutions */
   DECL_BRANCHEXPS((*branchexps)),      /**< branching execution method for not completely fixed pseudo solutions */
   BRANCHRULEDATA*  branchruledata      /**< branching rule data */
   )
{
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
   (*branchrule)->branchexlp = branchexlp;
   (*branchrule)->branchexps = branchexps;
   (*branchrule)->branchruledata = branchruledata;
   (*branchrule)->initialized = FALSE;

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
      CHECK_OKAY( (*branchrule)->branchfree(*branchrule, scip) );
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
      char s[255];
      sprintf(s, "branching rule <%s> already initialized", branchrule->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( branchrule->branchinit != NULL )
   {
      CHECK_OKAY( branchrule->branchinit(branchrule, scip) );
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
      char s[255];
      sprintf(s, "branching rule <%s> not initialized", branchrule->name);
      errorMessage(s);
      return SCIP_INVALIDCALL;
   }

   if( branchrule->branchexit != NULL )
   {
      CHECK_OKAY( branchrule->branchexit(branchrule, scip) );
   }
   branchrule->initialized = FALSE;

   return SCIP_OKAY;
}

/** executes branching rule for fractional LP solution */
RETCODE SCIPbranchruleExecLPSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexlp != NULL )
   {
      CHECK_OKAY( branchrule->branchexlp(branchrule, scip, result) );
      if( *result != SCIP_CUTOFF
         && *result != SCIP_BRANCHED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_SEPARATED
         && *result != SCIP_DIDNOTRUN )
      {
         char s[255];
         sprintf(s, "branching rule <%s> returned invalid result code <%d> from LP solution branching",
            branchrule->name, *result);
         errorMessage(s);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** executes branching rule for not completely fixed pseudo solution */
RETCODE SCIPbranchruleExecPseudoSol(
   BRANCHRULE*      branchrule,         /**< branching rule */
   SCIP*            scip,               /**< SCIP data structure */   
   RESULT*          result              /**< pointer to store the result of the callback method */
   )
{
   assert(branchrule != NULL);
   assert(scip != NULL);
   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;
   if( branchrule->branchexps != NULL )
   {
      CHECK_OKAY( branchrule->branchexps(branchrule, scip, result) );
      if( *result != SCIP_CUTOFF
         && *result != SCIP_BRANCHED
         && *result != SCIP_REDUCEDDOM
         && *result != SCIP_DIDNOTRUN )
      {
         char s[255];
         sprintf(s, "branching rule <%s> returned invalid result code <%d> from LP solution branching",
            branchrule->name, *result);
         errorMessage(s);
         return SCIP_INVALIDRESULT;
      }
   }

   return SCIP_OKAY;
}

/** gets name of branching rule */
const char* SCIPbranchruleGetName(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->name;
}

/** gets priority of branching rule */
int SCIPbranchruleGetPriority(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->priority;
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

/** is branching rule initialized? */
Bool SCIPbranchruleIsInitialized(
   BRANCHRULE*      branchrule          /**< branching rule */
   )
{
   assert(branchrule != NULL);

   return branchrule->initialized;
}

