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

/**@file   prob.c
 * @brief  Methods and datastructures for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "prob.h"



/*
 * dymanic memory arrays
 */

static
RETCODE probEnsureVarsMem(              /**< resizes vars array to be able to store at least num entries */
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   if( num > prob->varssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&prob->vars, newsize) );
      prob->varssize = newsize;
   }
   assert(num <= prob->varssize);

   return SCIP_OKAY;
}

static
RETCODE probEnsureConssMem(             /**< resizes conss array to be able to store at least num entries */
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   if( num > prob->consssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&prob->conss, newsize) );
      prob->consssize = newsize;
   }
   assert(num <= prob->consssize);

   return SCIP_OKAY;
}



/*
 * problem creation
 */

RETCODE SCIPprobCreate(                 /**< creates problem data structure */
   PROB**           prob,               /**< pointer to problem data structure */
   const char*      name                /**< problem name */
   )
{
   assert(prob != NULL);

   ALLOC_OKAY( allocMemory(prob) );
   ALLOC_OKAY( duplicateMemoryArray(&(*prob)->name, name, strlen(name)+1) );

   (*prob)->fixedvars = NULL;
   (*prob)->vars = NULL;
   CHECK_OKAY( SCIPhashtableCreate(&(*prob)->varnames, SCIP_HASHSIZE_NAMES,
                  SCIPhashGetKeyVar, SCIPhashKeyEqString, SCIPhashKeyValString) );
   (*prob)->conss = NULL;
   CHECK_OKAY( SCIPhashtableCreate(&(*prob)->consnames, SCIP_HASHSIZE_NAMES,
                  SCIPhashGetKeyCons, SCIPhashKeyEqString, SCIPhashKeyValString) );
   (*prob)->objsense = SCIP_OBJSENSE_MINIMIZE;
   (*prob)->objoffset = 0.0;
   (*prob)->objlim = SCIP_INVALID;
   (*prob)->fixedvarssize = 0;
   (*prob)->nfixedvars = 0;
   (*prob)->varssize = 0;
   (*prob)->nvars = 0;
   (*prob)->nbin = 0;
   (*prob)->nint = 0;
   (*prob)->nimpl = 0;
   (*prob)->ncont = 0;
   (*prob)->consssize = 0;
   (*prob)->nconss = 0;
   (*prob)->nmodelconss = 0;

   return SCIP_OKAY;
}

RETCODE SCIPprobFree(                   /**< frees problem data structure */
   PROB**           prob,               /**< pointer to problem data structure */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's the original problem) */
   )
{
   int c;
   int v;

   assert(prob != NULL);
   assert(*prob != NULL);
   
   freeMemoryArray(&(*prob)->name);

   /* release constraints and free constraint array */
   for( c = 0; c < (*prob)->nconss; ++c )
   {
      CHECK_OKAY( SCIPconsRelease(&(*prob)->conss[c], memhdr, set) );
   }
   freeMemoryArrayNull(&(*prob)->conss);
   
   /* release problem variables */
   for( v = 0; v < (*prob)->nvars; ++v )
   {
      assert((*prob)->vars[v]->probindex >= 0);
      (*prob)->vars[v]->probindex = -1;
      CHECK_OKAY( SCIPvarRelease(&(*prob)->vars[v], memhdr, set, lp) );
   }
   freeMemoryArrayNull(&(*prob)->vars);

   /* free hash tables for names */
   SCIPhashtableFree(&(*prob)->varnames, memhdr);
   SCIPhashtableFree(&(*prob)->consnames, memhdr);

   freeMemory(prob);
   
   return SCIP_OKAY;
}

RETCODE SCIPprobTransform(              /**< transform problem data into normalized form */
   PROB*            source,             /**< problem to transform */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   PROB**           target              /**< pointer to target problem data structure */
   )
{
   VAR* targetvar;
   CONS* targetcons;
   int v;
   int c;

   assert(source != NULL);
   assert(memhdr != NULL);
   assert(target != NULL);

   debugMessage("transform problem: original has %d variables\n", source->nvars);

   /* create target problem data */
   CHECK_OKAY( SCIPprobCreate(target, source->name) );

   /* transform and copy all variables to target problem */
   CHECK_OKAY( probEnsureVarsMem(*target, set, source->nvars) );
   for( v = 0; v < source->nvars; ++v )
   {
      CHECK_OKAY( SCIPvarTransform(&targetvar, memhdr, set, stat, source->objsense, source->vars[v]) );
      CHECK_OKAY( SCIPprobAddVar(*target, memhdr, set, targetvar) );
      CHECK_OKAY( SCIPvarRelease(&targetvar, memhdr, set, NULL) );
   }
   assert((*target)->nvars == source->nvars);

   /* transform and copy all constraints to target problem */
   for( c = 0; c < source->nconss; ++c )
   {
      CHECK_OKAY( SCIPconsTransform(&targetcons, memhdr, set, source->conss[c]) );
      CHECK_OKAY( SCIPprobAddCons(*target, memhdr, set, targetcons) );
      CHECK_OKAY( SCIPconsRelease(&targetcons, memhdr, set) );
   }

   /* add transformed model constraints to the constraint handlers' probconss arrays */
   for( c = 0; c < (*target)->nmodelconss; ++c )
   {
      targetcons = (*target)->conss[c];
      assert(targetcons != NULL);
      CHECK_OKAY( SCIPconshdlrAddProbCons(SCIPconsGetConsHdlr(targetcons), set, targetcons) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPprobActivate(               /**< activates constraints in the problem */
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   int c;

   assert(prob != NULL);
   assert(prob->nconss == 0 || prob->conss != NULL);
   assert(set != NULL);

   for( c = 0; c < prob->nconss; ++c )
   {
      CHECK_OKAY( SCIPconsActivate(prob->conss[c], set) );
   }

   return SCIP_OKAY;
}

RETCODE SCIPprobDeactivate(             /**< deactivates constraints in the problem */
   PROB*            prob                /**< problem data */
   )
{
   int c;

   assert(prob != NULL);
   assert(prob->nconss == 0 || prob->conss != NULL);

   for( c = 0; c < prob->nconss; ++c )
   {
      CHECK_OKAY( SCIPconsDeactivate(prob->conss[c]) );
   }

   return SCIP_OKAY;
}



/*
 * problem modification
 */

static
void probInsertVar(                     /**< insert var at the correct position in vars array, depending on its type */
   PROB*            prob,               /**< problem data */
   VAR*             var                 /**< variable to insert */
   )
{
   int insertpos;
   int intstart;
   int implstart;
   int contstart;

   assert(prob != NULL);
   assert(prob->vars != NULL);
   assert(prob->nvars < prob->varssize);
   assert(var != NULL);
   assert(var->probindex == -1);

   /* insert variable in array */
   insertpos = prob->nvars;
   intstart = prob->nbin;
   implstart = intstart + prob->nint;
   contstart = implstart + prob->nimpl;

   if( var->vartype == SCIP_VARTYPE_CONTINOUS )
      prob->ncont++;
   else
   {
      if( insertpos > contstart )
      {
         prob->vars[insertpos] = prob->vars[contstart];
         prob->vars[insertpos]->probindex = insertpos;
         insertpos = contstart;
      }
      assert(insertpos == contstart);

      if( var->vartype == SCIP_VARTYPE_IMPLINT )
         prob->nimpl++;
      else
      {
         if( insertpos > implstart )
         {
            prob->vars[insertpos] = prob->vars[implstart];
            prob->vars[insertpos]->probindex = insertpos;
            insertpos = implstart;
         }
         assert(insertpos == implstart);

         if( var->vartype == SCIP_VARTYPE_INTEGER )
            prob->nint++;
         else
         {
            assert(var->vartype == SCIP_VARTYPE_BINARY);
            if( insertpos > intstart )
            {
               prob->vars[insertpos] = prob->vars[intstart];
               prob->vars[insertpos]->probindex = insertpos;
               insertpos = intstart;
            }
            assert(insertpos == intstart);

            prob->nbin++;
         }
      }
   }
   prob->nvars++;

   assert(prob->nvars == prob->nbin + prob->nint + prob->nimpl + prob->ncont);
   assert((var->vartype == SCIP_VARTYPE_BINARY && insertpos == prob->nbin - 1)
      || (var->vartype == SCIP_VARTYPE_INTEGER && insertpos == prob->nbin + prob->nint - 1)
      || (var->vartype == SCIP_VARTYPE_IMPLINT && insertpos == prob->nbin + prob->nint + prob->nimpl - 1)
      || (var->vartype == SCIP_VARTYPE_CONTINOUS && insertpos == prob->nbin + prob->nint + prob->nimpl + prob->ncont - 1));

   prob->vars[insertpos] = var;
   var->probindex = insertpos;
}

static
void probRemoveVar(                     /**< removes variable from vars array */
   PROB*            prob,               /**< problem data */
   VAR*             var                 /**< variable to remove */
   )
{
   int freepos;
   int intstart;
   int implstart;
   int contstart;

   assert(prob != NULL);
   assert(var != NULL);
   assert(var->probindex >= 0);
   assert(prob->vars != NULL);
   assert(prob->vars[var->probindex] == var);

   intstart = prob->nbin;
   implstart = intstart + prob->nint;
   contstart = implstart + prob->nimpl;

   switch( var->vartype )
   {
   case SCIP_VARTYPE_BINARY:
      assert(0 <= var->probindex && var->probindex < intstart);
      prob->nbin--;
      break;
   case SCIP_VARTYPE_INTEGER:
      assert(intstart <= var->probindex && var->probindex < implstart);
      prob->nint--;
      break;
   case SCIP_VARTYPE_IMPLINT:
      assert(implstart <= var->probindex && var->probindex < contstart);
      prob->nimpl--;
      break;
   case SCIP_VARTYPE_CONTINOUS:
      assert(contstart <= var->probindex && var->probindex < prob->nvars);
      prob->ncont--;
      break;
   default:
      errorMessage("unknown variable type");
      abort();
   }

   /* move last binary, last integer, last implicit, and last continous variable forward to fill the free slot */
   freepos = var->probindex;
   if( freepos < intstart-1 )
   {
      /* move last binary variable to free slot */
      prob->vars[freepos] = prob->vars[intstart-1];
      prob->vars[freepos]->probindex = freepos;
      freepos = intstart-1;
   }
   if( freepos < implstart-1 )
   {
      /* move last integer variable to free slot */
      prob->vars[freepos] = prob->vars[implstart-1];
      prob->vars[freepos]->probindex = freepos;
      freepos = implstart-1;
   }
   if( freepos < contstart-1 )
   {
      /* move last implicit integer variable to free slot */
      prob->vars[freepos] = prob->vars[contstart-1];
      prob->vars[freepos]->probindex = freepos;
      freepos = contstart-1;
   }
   if( freepos < prob->nvars-1 )
   {
      /* move last implicit integer variable to free slot */
      prob->vars[freepos] = prob->vars[prob->nvars-1];
      prob->vars[freepos]->probindex = freepos;
      freepos = prob->nvars-1;
   }
   assert(freepos == prob->nvars-1);

   prob->nvars--;
   var->probindex = -1;

   assert(prob->nvars == prob->nbin + prob->nint + prob->nimpl + prob->ncont);
}

RETCODE SCIPprobAddVar(                 /**< adds variable to the problem and captures it */
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< variable to add */
   )
{
   assert(prob != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(var->probindex == -1);

   /* allocate additional memory */
   CHECK_OKAY( probEnsureVarsMem(prob, set, prob->nvars+1) );
   
   /* insert variable in vars array and mark it to be in problem */
   probInsertVar(prob, var);

   /* capture variable */
   SCIPvarCapture(var);

   /* add variable's name to the namespace */
   CHECK_OKAY( SCIPhashtableInsert(prob->varnames, memhdr, (void*)var) );

   debugMessage("added variable <%s> to problem (%d variables: %d binary, %d integer, %d implicit, %d continous)\n",
      var->name, prob->nvars, prob->nbin, prob->nint, prob->nimpl, prob->ncont);

   return SCIP_OKAY;
}

RETCODE SCIPprobChgVarType(             /**< changes the type of a variable in the problem */
   PROB*            prob,               /**< problem data */
   VAR*             var,                /**< variable to add */
   VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(prob != NULL);
   assert(var != NULL);
   assert(var->probindex >= 0);

   if( var->vartype == vartype )
      return SCIP_OKAY;

   /* temporarily remove variable from problem */
   probRemoveVar(prob, var);

   /* change the type of the variable */
   CHECK_OKAY( SCIPvarChgType(var, vartype) );

   /* reinsert variable into problem */
   probInsertVar(prob, var);

   return SCIP_OKAY;
}

static
void probInsertCons(                    /**< inserts constraint into the problems constraint array */
   PROB*            prob,               /**< problem data */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(prob != NULL);
   assert(prob->nmodelconss <= prob->nconss);
   assert(prob->nconss < prob->consssize);
   assert(cons != NULL);

   if( SCIPconsIsModel(cons) )
   {
      prob->conss[prob->nconss] = prob->conss[prob->nmodelconss];
      prob->conss[prob->nmodelconss] = cons;
      prob->nmodelconss++;
      prob->nconss++;
   }
   else
   {
      prob->conss[prob->nconss] = cons;
      prob->nconss++;
   }
}

RETCODE SCIPprobAddCons(                /**< adds constraint to the problem and captures it */
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(prob != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);

   /* add the constraint to the problem's constraint array */
   CHECK_OKAY( probEnsureConssMem(prob, set, prob->nconss+1) );
   probInsertCons(prob, cons);

   /* capture constraint */
   SCIPconsCapture(cons);

   /* add constraint's name to the namespace */
   CHECK_OKAY( SCIPhashtableInsert(prob->consnames, memhdr, (void*)cons) );

   return SCIP_OKAY;
}

void SCIPprobSetObjsense(               /**< sets objective sense: minimization or maximization */
   PROB*            prob,               /**< problem data */
   OBJSENSE         objsense            /**< new objective sense */
   )
{
   assert(prob != NULL);
   assert(prob->objsense == SCIP_OBJSENSE_MAXIMIZE || prob->objsense == SCIP_OBJSENSE_MINIMIZE);
   assert(objsense == SCIP_OBJSENSE_MAXIMIZE || objsense == SCIP_OBJSENSE_MINIMIZE);

   if( prob->objlim < SCIP_INVALID && objsense != prob->objsense )
      prob->objlim *= -1;
   prob->objsense = objsense;
}

void SCIPprobSetObjlim(                 /**< sets limit on objective function, such that only solutions better than this
                                           limit are accepted */
   PROB*            prob,               /**< problem data */
   Real             objlim              /**< objective limit */
   )
{
   assert(prob != NULL);

   prob->objlim = objlim;
}

Real SCIPprobExternObjval(              /**< returns the external value of the given internal objective value */
   PROB*            prob,               /**< problem data */
   Real             objval              /**< internal objective value */
   )
{
   assert(prob != NULL);

   return prob->objsense * (objval + prob->objoffset);
}


/*
 * problem information
 */

const char* SCIPprobGetName(            /**< gets problem name */
   PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);
   return prob->name;
}

VAR* SCIPprobFindVar(                   /**< returns variable of the problem with given name */
   PROB*            prob,               /**< problem data */
   const char*      name                /**< name of variable to find */
   )
{
   assert(prob != NULL);
   assert(name != NULL);

   return (VAR*)(SCIPhashtableRetrieve(prob->varnames, (void*)name));
}

CONS* SCIPprobFindCons(                 /**< returns constraint of the problem with given name */
   PROB*            prob,               /**< problem data */
   const char*      name                /**< name of variable to find */
   )
{
   assert(prob != NULL);
   assert(name != NULL);

   return (CONS*)(SCIPhashtableRetrieve(prob->consnames, (void*)name));
}

void SCIPprobPrintPseudoSol(            /**< displays actual pseudo solution */
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   VAR* var;
   Real solval;
   int v;
   
   for( v = 0; v < prob->nvars; ++v )
   {
      var = prob->vars[v];
      assert(var != NULL);
      solval = SCIPvarGetPseudoSol(var);
      if( !SCIPsetIsZero(set, solval) )
         printf(" <%s>=%g", SCIPvarGetName(var), solval);
   }
   printf("\n");
}
