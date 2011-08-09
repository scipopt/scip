/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   prop_rootredcost.c
 * @ingroup PROPAGATORS
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
#define PROP_PRIORITY          +1000000 /**< propagator priority */ 
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


/** initialization method of propagator (called after problem was transformed) */
#define propInitRootredcost NULL


/** deinitialization method of propagator (called before transformed problem is freed) */
#define propExitRootredcost NULL


/** presolving initialization method of propagator (called when presolving is about to begin) */
#define propInitpreRootredcost NULL


/** presolving deinitialization method of propagator (called after presolving has been finished) */
#define propExitpreRootredcost NULL


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


/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
#define propExitsolRootredcost NULL


/** presolving method of propagator */
#define propPresolRootredcost NULL


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
   if ( !SCIPisRootLPRelax(scip) )
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
      SCIP_Real redcost;

      var = vars[v];
      redcost = SCIPvarGetRootRedcost(var);
      if( redcost == SCIP_INVALID ) /*lint !e777*/
         continue;

      if( SCIPisFeasPositive(scip, redcost) )
      {
         SCIP_Real rootsol;
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newub;
         SCIP_Bool strengthen;

         rootsol = SCIPvarGetRootSol(var);
         oldlb = SCIPvarGetLbGlobal(var);
         oldub = SCIPvarGetUbGlobal(var);
         assert(SCIPisLE(scip, rootsol, SCIPvarGetLbGlobal(var))); /* lb might have been increased in the meantime */

         /* calculate reduced cost based bound */
         newub = (cutoffbound - lpobjval) / redcost + rootsol;

         /* check, if new bound is good enough */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            newub = SCIPadjustedVarUb(scip, var, newub);
            strengthen = (newub < oldub - 0.5);
         }
         else
            strengthen = SCIPisUbBetter(scip, newub, oldlb, oldub);

         if( strengthen )
         {
            SCIP_Bool infeasible;
            SCIP_Bool tightened;

            /* strengthen upper bound */
            SCIPdebugMessage("root redcost strengthening upper bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
               SCIPvarGetName(var), oldlb, oldub, oldlb, newub, cutoffbound, lpobjval, redcost);
            SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newub, FALSE, &infeasible, &tightened) );
            if( infeasible )
            {
               /* we are done with solving: cutoff root node */
               SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
               *result = SCIP_CUTOFF;
               break;
            }
            else if( tightened )
               *result = SCIP_REDUCEDDOM;
         }
      }
      else if( SCIPisFeasNegative(scip, redcost) )
      {
         SCIP_Real rootsol;
         SCIP_Real oldlb;
         SCIP_Real oldub;
         SCIP_Real newlb;
         SCIP_Bool strengthen;

         rootsol = SCIPvarGetRootSol(var);
         oldlb = SCIPvarGetLbGlobal(var);
         oldub = SCIPvarGetUbGlobal(var);
         assert(SCIPisGE(scip, rootsol, SCIPvarGetUbGlobal(var))); /* ub might have been decreased in the meantime */

         /* calculate reduced cost based bound */
         newlb = (cutoffbound - lpobjval) / redcost + rootsol;

         /* check, if new bound is good enough */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            newlb = SCIPadjustedVarLb(scip, var, newlb);
            strengthen = (newlb > oldlb + 0.5);
         }
         else
            strengthen = SCIPisLbBetter(scip, newlb, oldlb, oldub);

         if( strengthen )
         {
            SCIP_Bool infeasible;
            SCIP_Bool tightened;

            /* strengthen lower bound */
            SCIPdebugMessage("root redcost strengthening lower bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
               SCIPvarGetName(var), oldlb, oldub, newlb, oldub, cutoffbound, lpobjval, redcost);
            SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newlb, FALSE, &infeasible, &tightened) );
            if( infeasible )
            {
               /* we are done with solving: cutoff root node */
               SCIP_CALL( SCIPcutoffNode(scip, SCIPgetRootNode(scip)) );
               *result = SCIP_CUTOFF;
               break;
            }
            else if( tightened )
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

   /* create rootredcost propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->lastcutoffbound = SCIPinfinity(scip);

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOL_DELAY,
         propCopyRootredcost,
         propFreeRootredcost, propInitRootredcost, propExitRootredcost, propInitpreRootredcost, propExitpreRootredcost,
         propInitsolRootredcost, propExitsolRootredcost, propPresolRootredcost, propExecRootredcost, propRespropRootredcost,
         propdata) );

   return SCIP_OKAY;
}
