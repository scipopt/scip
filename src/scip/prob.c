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
      ALLOC_OKAY( reallocMemoryArray(prob->vars, newsize) );
      prob->varssize = newsize;
   }
   assert(num <= prob->varssize);

   return SCIP_OKAY;
}



/*
 * problem information
 */

const char* SCIPprobGetName(            /**< gets problem name */
   const PROB*      prob                /**< problem data */
   )
{
   assert(prob != NULL);
   return prob->name;
}




/*
 * problem modification
 */

RETCODE SCIPprobAddVar(                 /**< adds variable to the problem */
   PROB*            prob,               /**< problem data */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var                 /**< variable to add */
   )
{
   int insertpos;
   int intstart;
   int implstart;
   int contstart;

   assert(prob != NULL);
   assert(set != NULL);
   assert(var != NULL);

   CHECK_OKAY( probEnsureVarsMem(prob, set, prob->nvars+1) );

   /* insert variable at the correct position, depending on its type */
   insertpos = prob->nvars;
   intstart = prob->nbin;
   implstart = intstart + prob->nint;
   contstart = implstart + prob->nimpl;
   if( var->vartype == SCIP_VARTYPE_CONTINOUS )
      prob->ncont++;
   else
   {
      prob->vars[insertpos] = prob->vars[contstart];
      insertpos = contstart;
      if( var->vartype == SCIP_VARTYPE_IMPLINT )
         prob->nimpl++;
      else
      {
         prob->vars[insertpos] = prob->vars[implstart];
         insertpos = implstart;
         if( var->vartype == SCIP_VARTYPE_INTEGER )
            prob->nint++;
         else
         {
            assert(var->vartype == SCIP_VARTYPE_BINARY);
            prob->vars[insertpos] = prob->vars[intstart];
            insertpos = intstart;
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

   return SCIP_OKAY;
}

RETCODE SCIPprobAddCons(                /**< adds constraint to the problem */
   PROB*            prob,               /**< problem data */
   MEMHDR*          memhdr,             /**< block memory */
   CONS*            cons                /**< constraint to add */
   )
{
   assert(prob != NULL);
   assert(memhdr != NULL);
   assert(cons != NULL);

   CHECK_OKAY( SCIPconslistAdd(&(prob->conslist), memhdr, cons) );

   return SCIP_OKAY;
}



RETCODE SCIPprobCreate(                 /**< creates problem data structure */
   PROB**           prob,               /**< pointer to problem data structure */
   const char*      name                /**< problem name */
   )
{
   assert(prob != NULL);

   ALLOC_OKAY( allocMemory(*prob) );
   ALLOC_OKAY( allocMemoryArray((*prob)->name, strlen(name)+1) );
   copyMemoryArray((*prob)->name, name, strlen(name)+1);

   (*prob)->vars = NULL;
   (*prob)->conslist = NULL;
   (*prob)->varssize = 0;
   (*prob)->nvars = 0;
   (*prob)->nbin = 0;
   (*prob)->nint = 0;
   (*prob)->nimpl = 0;
   (*prob)->ncont = 0;

   return SCIP_OKAY;
}

RETCODE SCIPprobFree(                   /**< frees problem data structure */
   PROB**           prob,               /**< pointer to problem data structure */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's the original problem) */
   )
{
   int v;

   assert(prob != NULL);
   assert(*prob != NULL);
   
   freeMemoryArray((*prob)->name);

   /* free problem variables */
   for( v = 0; v < (*prob)->nvars; ++v )
   {
      SCIPvarFree(&(*prob)->vars[v], memhdr, set, lp);
   }
   freeMemoryArrayNull((*prob)->vars);

   /* free constraints */
   SCIPconslistFree(&(*prob)->conslist, memhdr, set);

   freeMemory(*prob);
   
   return SCIP_OKAY;
}

RETCODE SCIPprobTransform(              /**< transform problem data into normalized form */
   PROB*            source,             /**< problem to transform */
   MEMHDR*          memhdr,             /**< block memory buffer */
   const SET*       set,                /**< global SCIP settings */
   PROB**           target              /**< pointer to target problem data structure */
   )
{
   VAR* targetvar;
   CONS* targetcons;
   CONSLIST* conslist;
   int v;

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
      CHECK_OKAY( SCIPvarTransform(source->vars[v], memhdr, set, &targetvar) );
      CHECK_OKAY( SCIPprobAddVar(*target, set, targetvar) );
   }
   assert((*target)->nvars == source->nvars);

   /* transform and copy all constraints to target problem */
   conslist = source->conslist;
   while( conslist != NULL )
   {
      CHECK_OKAY( SCIPconsTransform(conslist->cons, memhdr, set, &targetcons) );
      CHECK_OKAY( SCIPprobAddCons(*target, memhdr, targetcons) );
      conslist = conslist->next;
   }

   return SCIP_OKAY;
}
