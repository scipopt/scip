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
#pragma ident "@(#) $Id: prob.c,v 1.40 2004/02/04 17:27:35 bzfpfend Exp $"

/**@file   prob.c
 * @brief  Methods and datastructures for storing and manipulating the main problem
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "misc.h"
#include "lp.h"
#include "var.h"
#include "prob.h"
#include "tree.h"
#include "branch.h"
#include "cons.h"



/*
 * dymanic memory arrays
 */

/** resizes fixedvars array to be able to store at least num entries */
static
RETCODE probEnsureFixedvarsMem(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimal number of slots in array */
   )
{
   assert(prob != NULL);
   assert(set != NULL);

   if( num > prob->fixedvarssize )
   {
      int newsize;

      newsize = SCIPsetCalcMemGrowSize(set, num);
      ALLOC_OKAY( reallocMemoryArray(&prob->fixedvars, newsize) );
      prob->fixedvarssize = newsize;
   }
   assert(num <= prob->fixedvarssize);

   return SCIP_OKAY;
}

/** resizes vars array to be able to store at least num entries */
static
RETCODE probEnsureVarsMem(
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

/** resizes conss array to be able to store at least num entries */
static
RETCODE probEnsureConssMem(
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

/** creates problem data structure
 *  If the problem type requires the use of variable pricers, these pricers should be activated with calls
 *  to SCIPactivatePricer(). These pricers are automatically deactivated, when the problem is freed.
 */
RETCODE SCIPprobCreate(
   PROB**           prob,               /**< pointer to problem data structure */
   const char*      name,               /**< problem name */
   DECL_PROBDELORIG ((*probdelorig)),   /**< frees user data of original problem */
   DECL_PROBTRANS   ((*probtrans)),     /**< creates user data of transformed problem by transforming original user data */
   DECL_PROBDELTRANS((*probdeltrans)),  /**< frees user data of transformed problem */
   PROBDATA*        probdata,           /**< user problem data set by the reader */
   Bool             transformed         /**< is this the transformed problem? */
   )
{
   assert(prob != NULL);

   ALLOC_OKAY( allocMemory(prob) );
   ALLOC_OKAY( duplicateMemoryArray(&(*prob)->name, name, strlen(name)+1) );

   (*prob)->probdata = probdata;
   (*prob)->probdelorig = probdelorig;
   (*prob)->probtrans = probtrans;
   (*prob)->probdeltrans = probdeltrans;
   CHECK_OKAY( SCIPhashtableCreate(&(*prob)->varnames, SCIP_HASHSIZE_NAMES,
                  SCIPhashGetKeyVar, SCIPhashKeyEqString, SCIPhashKeyValString) );
   (*prob)->fixedvars = NULL;
   (*prob)->fixedvarssize = 0;
   (*prob)->nfixedvars = 0;
   (*prob)->vars = NULL;
   (*prob)->varssize = 0;
   (*prob)->nvars = 0;
   (*prob)->nbinvars = 0;
   (*prob)->nintvars = 0;
   (*prob)->nimplvars = 0;
   (*prob)->ncontvars = 0;
   (*prob)->ncolvars = 0;
   CHECK_OKAY( SCIPhashtableCreate(&(*prob)->consnames, SCIP_HASHSIZE_NAMES,
                  SCIPhashGetKeyCons, SCIPhashKeyEqString, SCIPhashKeyValString) );
   (*prob)->conss = NULL;
   (*prob)->consssize = 0;
   (*prob)->nconss = 0;
   (*prob)->maxnconss = 0;
   (*prob)->startnconss = 0;
   (*prob)->objsense = SCIP_OBJSENSE_MINIMIZE;
   (*prob)->objoffset = 0.0;
   (*prob)->objlim = SCIP_INVALID;
   (*prob)->objisintegral = FALSE;
   (*prob)->transformed = transformed;

   return SCIP_OKAY;
}

/** frees problem data structure */
RETCODE SCIPprobFree(
   PROB**           prob,               /**< pointer to problem data structure */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's the original problem) */
   )
{
   int v;

   assert(prob != NULL);
   assert(*prob != NULL);
   assert(set != NULL);
   
   /* free user problem data */
   if( (*prob)->transformed )
   {
      if( (*prob)->probdeltrans != NULL )
      {
         CHECK_OKAY( (*prob)->probdeltrans(set->scip, &(*prob)->probdata) );
      }
   }
   else
   {
      if( (*prob)->probdelorig != NULL )
      {
         CHECK_OKAY( (*prob)->probdelorig(set->scip, &(*prob)->probdata) );
      }
   }

   /* remove all constraints from the problem */
   while( (*prob)->nconss > 0 )
   {
      assert((*prob)->conss != NULL);
      CHECK_OKAY( SCIPprobDelCons(*prob, memhdr, set, (*prob)->conss[0]) );
   }

   freeMemoryArray(&(*prob)->name);

   /* free constraint array */
   freeMemoryArrayNull(&(*prob)->conss);
   
   /* release problem variables */
   for( v = 0; v < (*prob)->nvars; ++v )
   {
      assert((*prob)->vars[v]->probindex >= 0);
      (*prob)->vars[v]->probindex = -1;
      CHECK_OKAY( SCIPvarRelease(&(*prob)->vars[v], memhdr, set, lp) );
   }
   freeMemoryArrayNull(&(*prob)->vars);

   /* release fixed problem variables */
   for( v = 0; v < (*prob)->nfixedvars; ++v )
   {
      assert((*prob)->fixedvars[v]->probindex == -1);
      CHECK_OKAY( SCIPvarRelease(&(*prob)->fixedvars[v], memhdr, set, lp) );
   }
   freeMemoryArrayNull(&(*prob)->fixedvars);

   /* free hash tables for names */
   SCIPhashtableFree(&(*prob)->varnames, memhdr);
   SCIPhashtableFree(&(*prob)->consnames, memhdr);

   freeMemory(prob);
   
   return SCIP_OKAY;
}

/** transform problem data into normalized form */
RETCODE SCIPprobTransform(
   PROB*            source,             /**< problem to transform */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   PROB**           target              /**< pointer to target problem data structure */
   )
{
   VAR* targetvar;
   CONS* targetcons;
   char transname[MAXSTRLEN];
   int v;
   int c;

   assert(set != NULL);
   assert(source != NULL);
   assert(memhdr != NULL);
   assert(target != NULL);

   debugMessage("transform problem: original has %d variables\n", source->nvars);

   /* create target problem data (probdelorig and probtrans are not needed, probdata is set later) */
   sprintf(transname, "t_%s", source->name);
   CHECK_OKAY( SCIPprobCreate(target, transname, NULL, NULL, source->probdeltrans, NULL, TRUE) );

   /* transform objective limit */
   if( source->objlim < SCIP_INVALID )
      SCIPprobSetExternObjlim(*target, SCIPprobInternObjval(source, set, source->objlim));

   /* transform and copy all variables to target problem */
   CHECK_OKAY( probEnsureVarsMem(*target, set, source->nvars) );
   for( v = 0; v < source->nvars; ++v )
   {
      CHECK_OKAY( SCIPvarTransform(source->vars[v], memhdr, set, stat, source->objsense, &targetvar) );
      CHECK_OKAY( SCIPprobAddVar(*target, memhdr, set, lp, branchcand, targetvar) );
      CHECK_OKAY( SCIPvarRelease(&targetvar, memhdr, set, NULL) );
   }
   assert((*target)->nvars == source->nvars);

   /* transform and copy all constraints to target problem */
   for( c = 0; c < source->nconss; ++c )
   {
      CHECK_OKAY( SCIPconsTransform(source->conss[c], memhdr, set, &targetcons) );
      CHECK_OKAY( SCIPprobAddCons(*target, memhdr, set, targetcons) );
      CHECK_OKAY( SCIPconsRelease(&targetcons, memhdr, set) );
   }

   /* objective value is always integral, iff original objective value is always integral and shift is integral */
   (*target)->objisintegral = source->objisintegral && SCIPsetIsIntegral(set, (*target)->objoffset);

   /* call user data transformation */
   if( source->probtrans != NULL )
   {
      CHECK_OKAY( source->probtrans(set->scip, source->probdata, &(*target)->probdata) );
   }
   else
      (*target)->probdata = source->probdata;

   return SCIP_OKAY;
}




/*
 * problem modification
 */

/** sets user problem data */
void SCIPprobSetData(
   PROB*            prob,               /**< problem */
   PROBDATA*        probdata            /**< user problem data to use */
   )
{
   assert(prob != NULL);

   prob->probdata = probdata;
}

/** inserts variable at the correct position in vars array, depending on its type */
static
void probInsertVar(
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
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   /* insert variable in array */
   insertpos = prob->nvars;
   intstart = prob->nbinvars;
   implstart = intstart + prob->nintvars;
   contstart = implstart + prob->nimplvars;

   if( SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS )
      prob->ncontvars++;
   else
   {
      if( insertpos > contstart )
      {
         prob->vars[insertpos] = prob->vars[contstart];
         prob->vars[insertpos]->probindex = insertpos;
         insertpos = contstart;
      }
      assert(insertpos == contstart);

      if( SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT )
         prob->nimplvars++;
      else
      {
         if( insertpos > implstart )
         {
            prob->vars[insertpos] = prob->vars[implstart];
            prob->vars[insertpos]->probindex = insertpos;
            insertpos = implstart;
         }
         assert(insertpos == implstart);

         if( SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER )
            prob->nintvars++;
         else
         {
            assert(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY);
            if( insertpos > intstart )
            {
               prob->vars[insertpos] = prob->vars[intstart];
               prob->vars[insertpos]->probindex = insertpos;
               insertpos = intstart;
            }
            assert(insertpos == intstart);

            prob->nbinvars++;
         }
      }
   }
   prob->nvars++;

   assert(prob->nvars == prob->nbinvars + prob->nintvars + prob->nimplvars + prob->ncontvars);
   assert((SCIPvarGetType(var) == SCIP_VARTYPE_BINARY && insertpos == prob->nbinvars - 1)
      || (SCIPvarGetType(var) == SCIP_VARTYPE_INTEGER && insertpos == prob->nbinvars + prob->nintvars - 1)
      || (SCIPvarGetType(var) == SCIP_VARTYPE_IMPLINT && insertpos == prob->nbinvars + prob->nintvars + prob->nimplvars - 1)
      || (SCIPvarGetType(var) == SCIP_VARTYPE_CONTINUOUS
         && insertpos == prob->nbinvars + prob->nintvars + prob->nimplvars + prob->ncontvars - 1));

   prob->vars[insertpos] = var;
   var->probindex = insertpos;

   /* update number of column variables in problem */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      prob->ncolvars++;
   assert(0 <= prob->ncolvars && prob->ncolvars <= prob->nvars);
}

/** removes variable from vars array */
static
void probRemoveVar(
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

   intstart = prob->nbinvars;
   implstart = intstart + prob->nintvars;
   contstart = implstart + prob->nimplvars;

   switch( SCIPvarGetType(var) )
   {
   case SCIP_VARTYPE_BINARY:
      assert(0 <= var->probindex && var->probindex < intstart);
      prob->nbinvars--;
      break;
   case SCIP_VARTYPE_INTEGER:
      assert(intstart <= var->probindex && var->probindex < implstart);
      prob->nintvars--;
      break;
   case SCIP_VARTYPE_IMPLINT:
      assert(implstart <= var->probindex && var->probindex < contstart);
      prob->nimplvars--;
      break;
   case SCIP_VARTYPE_CONTINUOUS:
      assert(contstart <= var->probindex && var->probindex < prob->nvars);
      prob->ncontvars--;
      break;
   default:
      errorMessage("unknown variable type\n");
      abort();
   }

   /* move last binary, last integer, last implicit, and last continuous variable forward to fill the free slot */
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

   assert(prob->nvars == prob->nbinvars + prob->nintvars + prob->nimplvars + prob->ncontvars);

   /* update number of column variables in problem */
   if( SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN )
      prob->ncolvars--;
   assert(0 <= prob->ncolvars && prob->ncolvars <= prob->nvars);
}

/** adds variable to the problem and captures it */
RETCODE SCIPprobAddVar(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< variable to add */
   )
{
   assert(prob != NULL);
   assert(set != NULL);
   assert(var != NULL);
   assert(var->probindex == -1);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   /* allocate additional memory */
   CHECK_OKAY( probEnsureVarsMem(prob, set, prob->nvars+1) );
   
   /* insert variable in vars array and mark it to be in problem */
   probInsertVar(prob, var);

   /* capture variable */
   SCIPvarCapture(var);

   /* add variable's name to the namespace */
   CHECK_OKAY( SCIPhashtableInsert(prob->varnames, memhdr, (void*)var) );

   /* update branching candidates and pseudo and loose objective value in the LP */
   if( SCIPvarGetStatus(var) != SCIP_VARSTATUS_ORIGINAL )
   {
      CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
      CHECK_OKAY( SCIPlpUpdateAddVar(lp, set, var) );
   }

   debugMessage("added variable <%s> to problem (%d variables: %d binary, %d integer, %d implicit, %d continuous)\n",
      SCIPvarGetName(var), prob->nvars, prob->nbinvars, prob->nintvars, prob->nimplvars, prob->ncontvars);

   return SCIP_OKAY;
}

/** changes the type of a variable in the problem */
RETCODE SCIPprobChgVarType(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var,                /**< variable to add */
   VARTYPE          vartype             /**< new type of variable */
   )
{
   assert(prob != NULL);
   assert(var != NULL);
   assert(var->probindex >= 0);
   assert(SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_LOOSE
      || SCIPvarGetStatus(var) == SCIP_VARSTATUS_COLUMN);

   if( SCIPvarGetType(var) == vartype )
      return SCIP_OKAY;

   /* temporarily remove variable from problem */
   probRemoveVar(prob, var);

   /* change the type of the variable */
   CHECK_OKAY( SCIPvarChgType(var, vartype) );

   /* reinsert variable into problem */
   probInsertVar(prob, var);

   /* update branching candidates */
   assert(branchcand != NULL || SCIPvarGetStatus(var) == SCIP_VARSTATUS_ORIGINAL);
   if( branchcand != NULL )
   {
      CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
   }

   return SCIP_OKAY;
}

/** informs problem, that the given loose problem variable changed its status */
RETCODE SCIPprobVarChangedStatus(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   VAR*             var                 /**< problem variable */
   )
{
   assert(prob != NULL);
   assert(var != NULL);
   assert(SCIPvarGetProbindex(var) != -1);

   /* get current status of variable */
   switch( SCIPvarGetStatus(var) )
   {
   case SCIP_VARSTATUS_ORIGINAL:
   case SCIP_VARSTATUS_LOOSE:
      errorMessage("variables cannot switch to ORIGINAL or LOOSE status\n");
      return SCIP_INVALIDDATA;

   case SCIP_VARSTATUS_COLUMN:
      /* variable switched from non-column to column */
      prob->ncolvars++;
      break;

   case SCIP_VARSTATUS_FIXED:
   case SCIP_VARSTATUS_AGGREGATED:
   case SCIP_VARSTATUS_MULTAGGR:
   case SCIP_VARSTATUS_NEGATED:
      /* variable switched from unfixed to fixed (if it was fixed before, probindex would have been -1) */

      /* remove variable from problem */
      probRemoveVar(prob, var);
      
      /* insert variable in fixedvars array */
      CHECK_OKAY( probEnsureFixedvarsMem(prob, set, prob->nfixedvars+1) );
      prob->fixedvars[prob->nfixedvars] = var;
      prob->nfixedvars++;

      /* update branching candidates */
      CHECK_OKAY( SCIPbranchcandUpdateVar(branchcand, set, var) );
      break;

   default:
      errorMessage("invalid variable status <%d>\n", SCIPvarGetStatus(var));
      return SCIP_INVALIDDATA;
   }
   assert(0 <= prob->ncolvars && prob->ncolvars <= prob->nvars);

   return SCIP_OKAY;
}

/** adds constraint to the problem and captures it; a local constraint is automatically upgraded into a global constraint */
RETCODE SCIPprobAddCons(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(prob != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);
   assert(cons->addconssetchg == NULL);
   assert(cons->addarraypos == -1);

   /* mark the constraint as problem constraint, and remember the constraint's position */
   cons->addconssetchg = NULL;
   cons->addarraypos = prob->nconss;

   /* add the constraint to the problem's constraint array */
   CHECK_OKAY( probEnsureConssMem(prob, set, prob->nconss+1) );
   prob->conss[prob->nconss] = cons;
   prob->nconss++;
   prob->maxnconss = MAX(prob->maxnconss, prob->nconss);

   /* undelete constraint, if it was globally deleted in the past */
   cons->deleted = FALSE;

   /* mark constraint to be globally valid */
   cons->local = FALSE;

   /* capture constraint */
   SCIPconsCapture(cons);

   /* add constraint's name to the namespace */
   CHECK_OKAY( SCIPhashtableInsert(prob->consnames, memhdr, (void*)cons) );

   /* if the problem is the transformed problem, activate and lock constraint */
   if( prob->transformed )
   {
      /* activate constraint */
      CHECK_OKAY( SCIPconsActivate(cons, set) );

      /* if constraint is a check-constraint, lock roundings of constraint's variables */
      if( SCIPconsIsChecked(cons) )
      {
         CHECK_OKAY( SCIPconsLockVars(cons, set, 1, 0) );
      }
   }

   return SCIP_OKAY;
}

/** releases and removes constraint from the problem; if the user has not captured the constraint for his own use, the
 *  constraint may be invalid after the call
 */
RETCODE SCIPprobDelCons(
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   CONS*            cons                /**< constraint to remove */
   )
{
   int arraypos;

   assert(prob != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);
   assert(cons->addconssetchg == NULL);
   assert(0 <= cons->addarraypos && cons->addarraypos < prob->nconss);
   assert(prob->conss != NULL);
   assert(prob->conss[cons->addarraypos] == cons);

   /* if the problem is the transformed problem, deactivate and unlock constraint */
   if( prob->transformed )
   {
      /* if constraint is a check-constraint, unlock roundings of constraint's variables */
      if( SCIPconsIsChecked(cons) )
      {
         CHECK_OKAY( SCIPconsUnlockVars(cons, set, 1, 0) );
      }

      /* deactivate constraint, if it is currently active */
      if( cons->active && !cons->updatedeactivate )
      {
         CHECK_OKAY( SCIPconsDeactivate(cons, set) );
      }
   }
   assert(!cons->active || cons->updatedeactivate);
   assert(!cons->enabled || cons->updatedeactivate);

   /* remove constraint's name from the namespace */
   CHECK_OKAY( SCIPhashtableRemove(prob->consnames, memhdr, (void*)cons) );

   /* remove the constraint from the problem's constraint array */
   arraypos = cons->addarraypos;
   prob->conss[arraypos] = prob->conss[prob->nconss-1];
   assert(prob->conss[arraypos] != NULL);
   assert(prob->conss[arraypos]->addconssetchg == NULL);
   prob->conss[arraypos]->addarraypos = arraypos;
   prob->nconss--;

   /* mark the constraint to be no longer in the problem */
   cons->addarraypos = -1;

   /* release constraint */
   CHECK_OKAY( SCIPconsRelease(&cons, memhdr, set) );

   return SCIP_OKAY;
}

/** informs problem, that the solution process is about to beginn;
 *  - resets maximum number of constraints to current number of constraints
 *  - remembers actual number of constraints as starting number of constraints
 */
void SCIPprobSolvingStarts(
   PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);

   /* remember number of constraints for statistic */
   prob->maxnconss = prob->nconss;
   prob->startnconss = prob->nconss;
}

/** sets objective sense: minimization or maximization */
void SCIPprobSetObjsense(
   PROB*            prob,               /**< problem data */
   OBJSENSE         objsense            /**< new objective sense */
   )
{
   assert(prob != NULL);
   assert(prob->objsense == SCIP_OBJSENSE_MAXIMIZE || prob->objsense == SCIP_OBJSENSE_MINIMIZE);
   assert(objsense == SCIP_OBJSENSE_MAXIMIZE || objsense == SCIP_OBJSENSE_MINIMIZE);

   prob->objsense = objsense;
}

/** increases objective offset */
void SCIPprobIncObjoffset(
   PROB*            prob,               /**< problem data */
   Real             incval              /**< value to add to objective offset */
   )
{
   assert(prob != NULL);

   prob->objoffset += incval;
}

/** sets limit on objective function, such that only solutions better than this limit are accepted */
void SCIPprobSetExternObjlim(
   PROB*            prob,               /**< problem data */
   Real             objlim              /**< external objective limit */
   )
{
   assert(prob != NULL);

   prob->objlim = objlim;
}

/** sets limit on objective function as transformed internal objective value */
void SCIPprobSetInternObjlim(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             objlim              /**< transformed internal objective limit */
   )
{
   assert(prob != NULL);

   prob->objlim = SCIPprobExternObjval(prob, set, objlim);
}

/** informs the problem, that its objective value is always integral in every feasible solution */
void SCIPprobSetObjIntegral(
   PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);
   
   prob->objisintegral = TRUE;
}

/** sets integral objective value flag, if all variables with non-zero objective values are integral and have 
 *  integral objective value
 */
void SCIPprobCheckObjIntegral(
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   Real obj;
   int v;

   assert(prob != NULL);
   
   /* if we know already, that the objective value is integral, nothing has to be done */
   if( prob->objisintegral )
      return;

   /* if there exist unknown variables, we cannot conclude that the objective value is always integral */
   if( set->nactivepricers != 0 )
      return;

   /* if the objective value offset is fractional, the value itself is possibly fractional */
   if( !SCIPsetIsIntegral(set, prob->objoffset) )
      return;

   /* scan through the variables */
   for( v = 0; v < prob->nvars; ++v )
   {
      /* get objective value of variable */
      obj = SCIPvarGetObj(prob->vars[v]);

      /* check, if objective value is non-zero */
      if( !SCIPsetIsZero(set, obj) )
      {
         /* if variable's objective value is fractional, the problem's objective value may also be fractional */
         if( !SCIPsetIsIntegral(set, obj) )
            break;
         
         /* if variable with non-zero objective value is continuous, the problem's objective value may be fractional */
         if( SCIPvarGetType(prob->vars[v]) == SCIP_VARTYPE_CONTINUOUS )
            break;
      }
   }

   /* objective value is integral, if the variable loop scanned all variables */
   prob->objisintegral = (v == prob->nvars);
}




/*
 * problem information
 */

/** gets problem name */
const char* SCIPprobGetName(
   PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);
   return prob->name;
}

/** gets user problem data */
PROBDATA* SCIPprobGetData(
   PROB*            prob                /**< problem */
   )
{
   assert(prob != NULL);

   return prob->probdata;
}

/** returns the external value of the given internal objective value */
Real SCIPprobExternObjval(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             objval              /**< internal objective value */
   )
{
   assert(prob != NULL);

   if( SCIPsetIsInfinity(set, objval) )
      return (Real)prob->objsense * set->infinity;
   else if( SCIPsetIsInfinity(set, -objval) )
      return -(Real)prob->objsense * set->infinity;
   else
      return (Real)prob->objsense * (objval + prob->objoffset);
}

/** returns the internal value of the given external objective value */
Real SCIPprobInternObjval(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   Real             objval              /**< external objective value */
   )
{
   assert(prob != NULL);

   if( SCIPsetIsInfinity(set, objval) )
      return (Real)prob->objsense * set->infinity;
   else if( SCIPsetIsInfinity(set, -objval) )
      return -(Real)prob->objsense * set->infinity;
   else
      return (Real)prob->objsense * objval - prob->objoffset;
}

/** gets limit on objective function in external space */
Real SCIPprobGetExternObjlim(
   PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);

   return prob->objlim;
}

/** gets limit on objective function as transformed internal objective value */
Real SCIPprobGetInternObjlim(
   PROB*            prob,               /**< problem data */
   const SET*       set                 /**< global SCIP settings */
   )
{
   assert(prob != NULL);

   return SCIPprobInternObjval(prob, set, prob->objlim);
}

/** returns whether the objective value is known to be integral in every feasible solution */
Bool SCIPprobIsObjIntegral(
   PROB*            prob                /**< problem data */
   )
{
   assert(prob != NULL);

   return prob->objisintegral;
}

/** returns variable of the problem with given name */
VAR* SCIPprobFindVar(
   PROB*            prob,               /**< problem data */
   const char*      name                /**< name of variable to find */
   )
{
   assert(prob != NULL);
   assert(name != NULL);

   return (VAR*)(SCIPhashtableRetrieve(prob->varnames, (void*)name));
}

/** returns constraint of the problem with given name */
CONS* SCIPprobFindCons(
   PROB*            prob,               /**< problem data */
   const char*      name                /**< name of variable to find */
   )
{
   assert(prob != NULL);
   assert(name != NULL);

   return (CONS*)(SCIPhashtableRetrieve(prob->consnames, (void*)name));
}

/** returns TRUE iff all columns, i.e. every variable with non-empty column w.r.t. all ever created rows, are present
 *  in the LP, and FALSE, if there are additional already existing columns, that may be added to the LP in pricing
 */
Bool SCIPprobAllColsInLP(
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   )
{
   assert(SCIPlpGetNCols(lp) <= prob->ncolvars && prob->ncolvars <= prob->nvars);

   return (SCIPlpGetNCols(lp) == prob->ncolvars && set->nactivepricers == 0);
}

/** displays actual pseudo solution */
void SCIPprobPrintPseudoSol(
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

/** outputs problem statistics */
void SCIPprobPrintStatistics(
   PROB*            prob,               /**< problem data */
   FILE*            file                /**< output file (or NULL for standard output) */
   )
{
   assert(prob != NULL);

   if( file == NULL )
      file = stdout;

   fprintf(file, "  Problem name     : %s\n", prob->name);
   fprintf(file, "  Variables        : %d (%d binary, %d integer, %d implicit integer, %d continuous)\n",
      prob->nvars, prob->nbinvars, prob->nintvars, prob->nimplvars, prob->ncontvars);
   fprintf(file, "  Constraints      : %d initial, %d maximal\n", prob->startnconss, prob->maxnconss);
}

