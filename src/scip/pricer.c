/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pricer.c,v 1.9 2005/01/18 09:26:51 bzfpfend Exp $"

/**@file   pricer.c
 * @brief  methods for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "def.h"
#include "message.h"
#include "set.h"
#include "clock.h"
#include "paramset.h"
#include "lp.h"
#include "prob.h"
#include "scip.h"
#include "pricer.h"

#include "struct_pricer.h"



/** compares two pricers w. r. to their activity and their priority */
DECL_SORTPTRCOMP(SCIPpricerComp)
{  /*lint --e{715}*/
   if( ((PRICER*)elem1)->active != ((PRICER*)elem2)->active )
      return ((PRICER*)elem1)->active ? -1 : +1;
   else
      return ((PRICER*)elem2)->priority - ((PRICER*)elem1)->priority;
}

/** method to call, when the priority of a pricer was changed */
static
DECL_PARAMCHGD(paramChgdPricerPriority)
{  /*lint --e{715}*/
   PARAMDATA* paramdata;

   paramdata = SCIPparamGetData(param);
   assert(paramdata != NULL);

   /* use SCIPsetPricerPriority() to mark the pricers unsorted */
   CHECK_OKAY( SCIPsetPricerPriority(scip, (PRICER*)paramdata, SCIPparamGetInt(param)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** creates a variable pricer
 *  To use the variable pricer for solving a problem, it first has to be activated with a call to SCIPactivatePricer().
 */
RETCODE SCIPpricerCreate(
   PRICER**         pricer,             /**< pointer to variable pricer data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of variable pricer */
   const char*      desc,               /**< description of variable pricer */
   int              priority,           /**< priority of the variable pricer */
   DECL_PRICERFREE  ((*pricerfree)),    /**< destructor of variable pricer */
   DECL_PRICERINIT  ((*pricerinit)),    /**< initialize variable pricer */
   DECL_PRICEREXIT  ((*pricerexit)),    /**< deinitialize variable pricer */
   DECL_PRICERREDCOST((*pricerredcost)),/**< reduced cost pricing method of variable pricer for feasible LPs */
   DECL_PRICERFARKAS((*pricerfarkas)),  /**< farkas pricing method of variable pricer for infeasible LPs */
   PRICERDATA*      pricerdata          /**< variable pricer data */
   )
{
   char paramname[MAXSTRLEN];
   char paramdesc[MAXSTRLEN];

   assert(pricer != NULL);
   assert(name != NULL);
   assert(desc != NULL);
   assert(pricerredcost != NULL);

   ALLOC_OKAY( allocMemory(pricer) );
   ALLOC_OKAY( duplicateMemoryArray(&(*pricer)->name, name, strlen(name)+1) );
   ALLOC_OKAY( duplicateMemoryArray(&(*pricer)->desc, desc, strlen(desc)+1) );
   (*pricer)->priority = priority;
   (*pricer)->pricerfree = pricerfree;
   (*pricer)->pricerinit = pricerinit;
   (*pricer)->pricerexit = pricerexit;
   (*pricer)->pricerredcost = pricerredcost;
   (*pricer)->pricerfarkas = pricerfarkas;
   (*pricer)->pricerdata = pricerdata;
   CHECK_OKAY( SCIPclockCreate(&(*pricer)->clock, SCIP_CLOCKTYPE_DEFAULT) );
   (*pricer)->ncalls = 0;
   (*pricer)->nvarsfound = 0;
   (*pricer)->active = FALSE;
   (*pricer)->initialized = FALSE;

   /* add parameters */
   sprintf(paramname, "pricers/%s/priority", name);
   sprintf(paramdesc, "priority of pricer <%s>", name);
   CHECK_OKAY( SCIPsetAddIntParam(set, memhdr, paramname, paramdesc,
                  &(*pricer)->priority, priority, INT_MIN, INT_MAX, 
                  paramChgdPricerPriority, (PARAMDATA*)(*pricer)) ); /*lint !e740*/

   return SCIP_OKAY;
}

/** calls destructor and frees memory of variable pricer */
RETCODE SCIPpricerFree(
   PRICER**         pricer,             /**< pointer to variable pricer data structure */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(pricer != NULL);
   assert(*pricer != NULL);
   assert(!(*pricer)->initialized);
   assert(scip != NULL);

   /* call destructor of variable pricer */
   if( (*pricer)->pricerfree != NULL )
   {
      CHECK_OKAY( (*pricer)->pricerfree(scip, *pricer) );
   }

   SCIPclockFree(&(*pricer)->clock);
   freeMemoryArray(&(*pricer)->name);
   freeMemoryArray(&(*pricer)->desc);
   freeMemory(pricer);

   return SCIP_OKAY;
}

/** initializes variable pricer */
RETCODE SCIPpricerInit(
   PRICER*          pricer,             /**< variable pricer */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(pricer != NULL);
   assert(pricer->active);
   assert(scip != NULL);

   if( pricer->initialized )
   {
      errorMessage("variable pricer <%s> already initialized\n", pricer->name);
      return SCIP_INVALIDCALL;
   }

   SCIPclockReset(pricer->clock);

   pricer->ncalls = 0;
   pricer->nvarsfound = 0;

   if( pricer->pricerinit != NULL )
   {
      CHECK_OKAY( pricer->pricerinit(scip, pricer) );
   }
   pricer->initialized = TRUE;

   return SCIP_OKAY;
}

/** calls exit method of variable pricer */
RETCODE SCIPpricerExit(
   PRICER*          pricer,             /**< variable pricer */
   SCIP*            scip                /**< SCIP data structure */   
   )
{
   assert(pricer != NULL);
   assert(pricer->active);
   assert(scip != NULL);

   if( !pricer->initialized )
   {
      errorMessage("variable pricer <%s> not initialized\n", pricer->name);
      return SCIP_INVALIDCALL;
   }

   if( pricer->pricerexit != NULL )
   {
      CHECK_OKAY( pricer->pricerexit(scip, pricer) );
   }
   pricer->initialized = FALSE;

   return SCIP_OKAY;
}

/** activates pricer such that it is called in LP solving loop */
RETCODE SCIPpricerActivate(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);
   assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);

   if( !pricer->active )
   {
      pricer->active = TRUE;
      set->nactivepricers++;
      set->pricerssorted = FALSE;
   }

   return SCIP_OKAY;
}

/** deactivates pricer such that it is no longer called in LP solving loop */
RETCODE SCIPpricerDeactivate(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set                 /**< global SCIP settings */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);
   assert(SCIPgetStage(set->scip) == SCIP_STAGE_PROBLEM);

   if( pricer->active )
   {
      pricer->active = FALSE;
      set->nactivepricers--;
      set->pricerssorted = FALSE;
   }

   return SCIP_OKAY;
}

/** returns whether the given pricer is in use in the current problem */
Bool SCIPpricerIsActive(
   PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->active;
}

/** calls reduced cost pricing method of variable pricer */
RETCODE SCIPpricerRedcost(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob                /**< transformed problem */
   )
{
   int oldnvars;

   assert(pricer != NULL);
   assert(pricer->active);
   assert(pricer->pricerredcost != NULL);
   assert(set != NULL);
   assert(prob != NULL);

   debugMessage("executing reduced cost pricing of variable pricer <%s>\n", pricer->name);
   
   oldnvars = prob->nvars;
   
   /* start timing */
   SCIPclockStart(pricer->clock, set);
   
   /* call external method */
   CHECK_OKAY( pricer->pricerredcost(set->scip, pricer) );
   
   /* stop timing */
   SCIPclockStop(pricer->clock, set);
   
   /* evaluate result */
   pricer->ncalls++;
   pricer->nvarsfound += prob->nvars - oldnvars;
   
   return SCIP_OKAY;
}

/** calls farkas pricing method of variable pricer */
RETCODE SCIPpricerFarkas(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob                /**< transformed problem */
   )
{
   int oldnvars;

   assert(pricer != NULL);
   assert(pricer->active);
   assert(set != NULL);
   assert(prob != NULL);

   /* check, if pricer implemented a farkas pricing algorithm */
   if( pricer->pricerfarkas == NULL )
      return SCIP_OKAY;

   debugMessage("executing farkas pricing of variable pricer <%s>\n", pricer->name);
   
   oldnvars = prob->nvars;
   
   /* start timing */
   SCIPclockStart(pricer->clock, set);
   
   /* call external method */
   CHECK_OKAY( pricer->pricerfarkas(set->scip, pricer) );
   
   /* stop timing */
   SCIPclockStop(pricer->clock, set);
   
   /* evaluate result */
   pricer->ncalls++;
   pricer->nvarsfound += prob->nvars - oldnvars;
   
   return SCIP_OKAY;
}

/** depending on the LP's solution status, calls reduced cost or farkas pricing method of variable pricer */
RETCODE SCIPpricerExec(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< transformed problem */
   LP*              lp                  /**< LP data */
   )
{
   if( SCIPlpGetSolstat(lp) == SCIP_LPSOLSTAT_INFEASIBLE )
   {
      CHECK_OKAY( SCIPpricerFarkas(pricer, set, prob) );
   }
   else
   {
      CHECK_OKAY( SCIPpricerRedcost(pricer, set, prob) );
   }

   return SCIP_OKAY;
}

/** gets user data of variable pricer */
PRICERDATA* SCIPpricerGetData(
   PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->pricerdata;
}

/** sets user data of variable pricer; user has to free old data in advance! */
void SCIPpricerSetData(
   PRICER*          pricer,             /**< variable pricer */
   PRICERDATA*      pricerdata          /**< new variable pricer user data */
   )
{
   assert(pricer != NULL);

   pricer->pricerdata = pricerdata;
}

/** gets name of variable pricer */
const char* SCIPpricerGetName(
   PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->name;
}

/** gets description of variable pricer */
const char* SCIPpricerGetDesc(
   PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->desc;
}

/** gets priority of variable pricer */
int SCIPpricerGetPriority(
   PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->priority;
}

/** sets priority of variable pricer */
void SCIPpricerSetPriority(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the variable pricer */
   )
{
   assert(pricer != NULL);
   assert(set != NULL);
   
   pricer->priority = priority;
   set->pricerssorted = FALSE;
}

/** gets the number of times, the pricer was called and tried to find a variable with negative reduced costs */
int SCIPpricerGetNCalls(
   PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->ncalls;
}

/** gets the number of variables with negative reduced costs found by this pricer */
int SCIPpricerGetNVarsFound(
   PRICER*          pricer              /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->nvarsfound;
}

/** is variable pricer initialized? */
Bool SCIPpricerIsInitialized(
   PRICER*            pricer                /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return pricer->initialized;
}

/** gets time in seconds used in this pricer */
Real SCIPpricerGetTime(
   PRICER*            pricer                /**< variable pricer */
   )
{
   assert(pricer != NULL);

   return SCIPclockGetTime(pricer->clock);
}

