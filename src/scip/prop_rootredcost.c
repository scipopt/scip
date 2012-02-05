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

/**@file   prop_rootredcost.c
 * @brief  reduced cost strengthening at the root node
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_rootredcost.h"


#define PROP_NAME              "rootredcost"
#define PROP_DESC              "reduced cost strengthening at the root node"
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP
#define PROP_PRIORITY         +10000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_PRESOL_PRIORITY          0 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOL_DELAY          TRUE /**< should presolving be delay, if other presolvers found reductions?  */
#define PROP_PRESOL_MAXROUNDS         0 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */




/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_Real             lastcutoffbound;    /**< cutoff bound for which the root reduced costs were already processed */
};




/*
 * Callback methods of propagator
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyRootredcost)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropRootredcost(scip) );
 
   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeRootredcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}










/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolRootredcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   propdata->lastcutoffbound = SCIPinfinity(scip);

   return SCIP_OKAY;
}






/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecRootredcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_Real cutoffbound;
   SCIP_Real lpobjval;
   int nvars;
   int v;

   *result = SCIP_DIDNOTRUN;

   /* propagator can only be applied during solving stage */
   if( SCIPgetStage(scip) != SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* propagator can only be applied if the root lp was a valid relaxation */
   if( !SCIPisRootLPRelax(scip) )
      return SCIP_OKAY;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* get current cutoff bound */
   cutoffbound = SCIPgetCutoffbound(scip);
   assert(cutoffbound <= propdata->lastcutoffbound);
   if( cutoffbound == propdata->lastcutoffbound ) /*lint !e777*/
      return SCIP_OKAY;

   propdata->lastcutoffbound = cutoffbound;

   /* reduced cost strengthening can only be applied, if we have a finite upper bound on the LP value */
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   /* get root LP objective value */
   lpobjval = SCIPgetLPRootObjval(scip);
   if( lpobjval == SCIP_INVALID ) /*lint !e777*/
      return SCIP_OKAY;

   /* get variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   if( nvars == 0 )
      return SCIP_OKAY;

   SCIPdebugMessage("searching for root node reduced cost fixings due to new cutoffbound %g\n", cutoffbound);

   *result = SCIP_DIDNOTFIND;

   /* check reduced costs for variables that were columns of the root LP */
   for( v = 0; v < nvars; ++v )
   {
      SCIP_VAR* var;
      SCIP_Real rootsol;
      SCIP_Real redcost;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      var = vars[v];
      redcost = SCIPvarGetBestRootRedcost(var);
      if( redcost == SCIP_INVALID ) /*lint !e777*/
         continue;

      rootsol = SCIPvarGetBestRootSol(var);
      lpobjval = SCIPvarGetBestRootLPobjval(var);

      if( SCIPisFeasPositive(scip, redcost) )
      {
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newub;

         oldlb = SCIPvarGetLbGlobal(var);
         oldub = SCIPvarGetUbGlobal(var);
         assert(SCIPisFeasLE(scip, rootsol, SCIPvarGetLbGlobal(var))); /* lb might have been increased in the meantime */

         /* calculate reduced cost based bound */
         newub = (cutoffbound - lpobjval) / redcost + rootsol;

         /* strengthen upper bound */
         SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newub, FALSE, &infeasible, &tightened) );

         if( infeasible )
         {
            /* we are done with solving: cutoff root node */
            SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
            *result = SCIP_CUTOFF;
            break;
         }

         if( tightened )
         {
            SCIPdebugMessage("root redcost strengthening upper bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
               SCIPvarGetName(var), oldlb, oldub, oldlb, newub, cutoffbound, lpobjval, redcost);

            *result = SCIP_REDUCEDDOM;
         }
      }
      else if( SCIPisFeasNegative(scip, redcost) )
      {
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newlb;

         oldlb = SCIPvarGetLbGlobal(var);
         oldub = SCIPvarGetUbGlobal(var);
         assert(SCIPisGE(scip, rootsol, SCIPvarGetUbGlobal(var))); /* ub might have been decreased in the meantime */

         /* calculate reduced cost based bound */
         newlb = (cutoffbound - lpobjval) / redcost + rootsol;

         /* strengthen lower bound */
         SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newlb, FALSE, &infeasible, &tightened) );

         if( infeasible )
         {
            /* we are done with solving: cutoff root node */
            SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
            *result = SCIP_CUTOFF;
            break;
         }

         if( tightened )
         {
            SCIPdebugMessage("root redcost strengthening lower bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
               SCIPvarGetName(var), oldlb, oldub, newlb, oldub, cutoffbound, lpobjval, redcost);

            *result = SCIP_REDUCEDDOM;
         }
      }
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropRootredcost)
{  /*lint --e{715}*/
   SCIPerrorMessage("cannot resolve root node reduced cost fixings -- should not happen!\n");

   return SCIP_INVALIDCALL;
}




/*
 * propagator specific interface methods
 */

/** creates the root node reduced cost strengthening propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropRootredcost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;
   SCIP_PROP* prop;

   /* create rootredcost propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->lastcutoffbound = SCIPinfinity(scip);

   /* include propagator */
   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecRootredcost, propRespropRootredcost,
         propdata) );

   assert(prop != NULL);

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropCopy(scip, prop, propCopyRootredcost) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeRootredcost) );
   SCIP_CALL( SCIPsetPropInitsol(scip, prop, propInitsolRootredcost) );

   return SCIP_OKAY;
}
