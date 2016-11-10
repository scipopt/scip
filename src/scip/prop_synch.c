/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_synch.c
 * @brief  propagator for applying global bound changes that were communicated by other
 *         concurrent solvers
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include "scip/prop_synch.h"
#include "scip/concurrent.h"
#include "tpi/tpi.h"

/* fundamental propagator properties */
#define PROP_NAME              "synch"
#define PROP_DESC              "synchronization propagator"
#define PROP_PRIORITY                     (INT_MAX/4) /**< propagator priority */
#define PROP_FREQ                                  -1 /**< propagator frequency */
#define PROP_DELAY                              FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING            SCIP_PROPTIMING_ALWAYS /**< propagation timing mask */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_VAR**       bndvar;        /**< array of variables with a bound change */
   SCIP_Real*       bndval;        /**< array of new bound values */
   SCIP_BOUNDTYPE*  bndtype;       /**< array of bound types */
   int              nbnds;         /**< number of boundchanges */
   int              bndsize;       /**< current size of bound change array */
   SCIP_Longint     ntightened;    /**< number of tightened bounds */
   SCIP_Longint     ntightenedint; /**< number of tightened bounds of integer variables */
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of propagator
 */

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeSynch)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}

/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitSynch)
{  /*lint --e{715}*/
   SCIP_PROPDATA* data;

   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   data->bndsize = 0;
   data->nbnds = 0;
   data->bndvar = NULL;
   data->bndval = NULL;
   data->bndtype = NULL;
   data->ntightened = 0;
   data->ntightenedint = 0;

   return SCIP_OKAY;
}

/** deinitialization method of propagator (called before transformed problem is freed) */
static
SCIP_DECL_PROPEXIT(propExitSynch)
{  /*lint --e{715}*/
   SCIP_PROPDATA* data;

   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   SCIPfreeBlockMemoryArrayNull(scip, &data->bndvar, data->bndsize);
   SCIPfreeBlockMemoryArrayNull(scip, &data->bndval, data->bndsize);
   SCIPfreeBlockMemoryArrayNull(scip, &data->bndtype, data->bndsize);

   return SCIP_OKAY;
}

/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecSynch)
{  /*lint --e{715}*/
   SCIP_PROPDATA*  data;
   int             i;
   int             ntightened;
   int             ntightenedint;

   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   *result = SCIP_DIDNOTRUN;

   if( SCIPinProbing(scip) )
      return SCIP_OKAY;

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   SCIPdisableConcurrentBoundStorage(scip);
   *result = SCIP_DIDNOTFIND;
   ntightened = 0;
   ntightenedint = 0;

   for( i = 0; i < data->nbnds; ++i )
   {
      SCIP_Bool infeas, tightened;
      SCIP_CALL( SCIPvarGetProbvarBound(&data->bndvar[i], &data->bndval[i], &data->bndtype[i]) );

      /* cannot change bounds of multi-aggregated variables so skip this bound-change */
      if( SCIPvarGetStatus(data->bndvar[i]) == SCIP_VARSTATUS_MULTAGGR )
         continue;

      if( data->bndtype[i] == SCIP_BOUNDTYPE_LOWER )
      {
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, data->bndvar[i], data->bndval[i], FALSE, &infeas, &tightened) );
      }
      else
      {
         assert(data->bndtype[i] == SCIP_BOUNDTYPE_UPPER);
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, data->bndvar[i], data->bndval[i], FALSE, &infeas, &tightened) );
      }
      if( tightened )
      {
         ++ntightened;
         if( SCIPvarGetType(data->bndvar[i]) <= SCIP_VARTYPE_INTEGER )
            ++ntightenedint;
      }
      if( infeas )
      {
#ifndef NDEBUG
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "synch propagator found cutoff in thread %i\n", SCIPtpiGetThreadNum());
#endif
         data->ntightened += ntightened;
         data->ntightenedint += ntightenedint;
         *result = SCIP_CUTOFF;
         return SCIP_OKAY;
      }
   }

   data->nbnds = 0;
   SCIPenableConcurrentBoundStorage(scip);
   SCIPpropSetFreq(prop, -1);
   if( ntightened > 0 )
   {
      data->ntightened += ntightened;
      data->ntightenedint += ntightenedint;
      *result = SCIP_REDUCEDDOM;
      SCIPdebugMsg(scip, "tightened %i var bounds in thread %i of which %i where integral vars\n", ntightened, SCIPtpiGetThreadNum(), ntightenedint);
   }
   return SCIP_OKAY;
}

/*
 * propagator specific interface methods
 */

/** creates the synch propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropSynch(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create xyz propagator data */
   propdata = NULL;
   /* create propagator specific data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );

   prop = NULL;

   /* include propagator */

   /* use SCIPincludePropBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecSynch, propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeSynch) );
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitSynch) );
   SCIP_CALL( SCIPsetPropExit(scip, prop, propExitSynch) );

   return SCIP_OKAY;
}


/** adds a boundchange to the synch propagator */
EXTERN
SCIP_RETCODE SCIPpropSynchAddBndchg(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< synch propagator */
   SCIP_VAR*             var,                /**< variable for bound */
   SCIP_Real             val,                /**< value of bound */
   SCIP_BOUNDTYPE        bndtype             /**< type of bound */
   )
{
   SCIP_PROPDATA* data;

   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   if( data->nbnds + 1 > data->bndsize )
   {
      int newsize;
      newsize = SCIPcalcMemGrowSize(scip, data->nbnds+1);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &data->bndvar, data->bndsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &data->bndval, data->bndsize, newsize) );
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &data->bndtype, data->bndsize, newsize) );
      data->bndsize = newsize;
   }

   data->bndvar[data->nbnds] = var;
   data->bndval[data->nbnds] = val;
   data->bndtype[data->nbnds] = bndtype;

   if( data->nbnds == 0 )
   {
      SCIPpropSetFreq(prop, 1);
   }
   ++data->nbnds;

   return SCIP_OKAY;
}

/** gives the total number of tightened bounds found by the synch propagator */
SCIP_Longint SCIPpropSynchGetNTightenedBnds(
   SCIP_PROP*            prop                /**< synch propagator */
   )
{
   SCIP_PROPDATA* data;

   assert(prop != NULL);

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   return data->ntightened;
}

/** gives the total number of tightened bounds for integer variables found by the synch propagator */
SCIP_Longint SCIPpropSynchGetNTightenedIntBnds(
   SCIP_PROP*            prop                /**< synch propagator */
   )
{
   SCIP_PROPDATA* data;

   assert(prop != NULL);

   data = SCIPpropGetData(prop);
   assert(data != NULL);

   return data->ntightenedint;
}
