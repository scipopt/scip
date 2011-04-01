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

/**@file   prop_pseudoobj.c
 * @ingroup PROPAGATORS
 * @brief  pseudoobj propagator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_pseudoobj.h"


#define PROP_NAME              "pseudoobj"
#define PROP_DESC              "pseudo objective function propagator"
#define PROP_PRIORITY                 0
#define PROP_FREQ                     1
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

#define EVENTHDLR_NAME         "pseudoobj"
#define EVENTHDLR_DESC         "bound change event handler for pseudo objective function propagator"

#define DEFAULT_MAXCANDS            100 /**< maximal number of variables to look at in a single propagation round
                                         *   (-1: process all variables) */




/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   SCIP_EVENTHDLR*       eventhdlr;          /**< event handler for global bound change events */
   SCIP_VAR**            objvars;            /**< variables with non-zero objective */
   SCIP_Real             maxpseudoobjact;    /**< maximal global pseudo objective activity */
   int                   maxpseudoobjactinf; /**< number of coefficients contributing with infinite value to maxpseudoobjact */
   int                   nobjvars;           /**< number of variables with non-zero objective */
   int                   maxcands;           /**< maximal number of variables to look at in a single propagation round */
   int                   lastvarnum;         /**< last variable number that was looked at */
};




/*
 * Local methods
 */

/** resolves a propagation by supplying the variables whose bound changes increased the pseudo objective value */
static
SCIP_RETCODE resolvePropagation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL for conflict analysis initialization */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_VAR** objvars;
   SCIP_VAR* var;
   SCIP_Real obj;
   int nobjvars;
   int v;

   assert(propdata != NULL);

   /**@todo improve pseudo objective propagator conflict resolving method:
    *       only add bound changes up to the point, the primal bound is reached
    */

   /* the variables responsible for the propagation are the ones with
    *  - obj > 0 and local lb > global lb
    *  - obj < 0 and local ub < global ub
    */
   objvars = propdata->objvars;
   nobjvars = propdata->nobjvars;
   assert(nobjvars == 0 || objvars != NULL);

   for( v = 0; v < nobjvars; ++v )
   {
      var = objvars[v];
      if( var == infervar )
         continue;

      obj = SCIPvarGetObj(var);
      assert(!SCIPisZero(scip, obj));

      if( obj > 0.0 )
      {
         SCIP_Real loclb;
         SCIP_Real glblb;

         glblb = SCIPvarGetLbGlobal(var);
         loclb = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);
         if( SCIPisGT(scip, loclb, glblb) )
         {
            SCIP_CALL( SCIPaddConflictLb(scip, var, bdchgidx) );
         }
      }
      else
      {
         SCIP_Real locub;
         SCIP_Real glbub;

         glbub = SCIPvarGetUbGlobal(var);
         locub = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);
         if( SCIPisLT(scip, locub, glbub) )
         {
            SCIP_CALL( SCIPaddConflictUb(scip, var, bdchgidx) );
         }
      }
   }

   return SCIP_OKAY;
}

/** propagates the cutoff bound for the given variables */
static
SCIP_RETCODE propagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_VAR*             var,                /**< variables to propagate */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   int*                  nchgbds             /**< pointer to store the number of changed bounds */
   )
{
   SCIP_Real lb;
   SCIP_Real ub;
   SCIP_Real obj;
   
   assert(!SCIPisInfinity(scip, -pseudoobjval));
   assert(!SCIPisInfinity(scip, cutoffbound));
   assert(SCIPisLT(scip, pseudoobjval, cutoffbound) );
   
   lb = SCIPvarGetLbLocal(var);
   ub = SCIPvarGetUbLocal(var);
   
   if( SCIPisFeasEQ(scip, lb, ub) )
      return SCIP_OKAY;
   
   obj = SCIPvarGetObj(var);
   assert(!SCIPisZero(scip, obj));

   if( obj > 0.0 )
   {
      SCIP_Real newub;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      newub = lb + (cutoffbound - pseudoobjval)/obj;
      if( SCIPisUbBetter(scip, newub, lb, ub) )
      {
         SCIPdebugMessage(" -> new upper bound of variable <%s>[%.10f,%.10f]: %.10f\n",
            SCIPvarGetName(var), lb, ub, newub);
         SCIP_CALL( SCIPinferVarUbProp(scip, var, newub, prop, 0, FALSE, &infeasible, &tightened) );
         assert(!infeasible);
         if( tightened ) /* might not be tightened due to numerical reasons */
            (*nchgbds)++;
      }
   }
   else
   {
      SCIP_Real newlb;
      SCIP_Bool infeasible;
      SCIP_Bool tightened;

      newlb = ub + (cutoffbound - pseudoobjval)/obj;
      if( SCIPisLbBetter(scip, newlb, lb, ub) )
      {
         SCIPdebugMessage(" -> new lower bound of variable <%s>[%g,%g]: %g\n", SCIPvarGetName(var), lb, ub, newlb);
         SCIP_CALL( SCIPinferVarLbProp(scip, var, newlb, prop, 0, FALSE, &infeasible, &tightened) );
         assert(!infeasible);
         if( tightened ) /* might not be tightened due to numerical reasons */
            (*nchgbds)++;
      }
   }
   
   return SCIP_OKAY;
}

/** propagates the cutoff bound c*x <= cutoff */
static
SCIP_RETCODE propagateCutoffbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )

{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** objvars;
   SCIP_Real pseudoobjval;
   SCIP_Real cutoffbound;
   int nchgbds;
   int ncands;
   int nobjvars;
   int c;
   int v;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   objvars = propdata->objvars;
   nobjvars = propdata->nobjvars;
   assert(nobjvars == 0 || objvars != NULL);

   /* nothing to do for empty objective */
   if( nobjvars == 0 )
      return SCIP_OKAY;

   /* get current pseudo objective value and cutoff bound */
   pseudoobjval = SCIPgetPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;
   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;
   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      SCIPdebugMessage("pseudo objective value %g exceeds cutoff bound %g\n", pseudoobjval, cutoffbound);

      /* initialize conflict analysis, and add all variables of infeasible constraint to conflict candidate queue */
      SCIP_CALL( SCIPinitConflictAnalysis(scip) );
      SCIP_CALL( resolvePropagation(scip, propdata, NULL, NULL) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
      
      *result = SCIP_CUTOFF;

      return SCIP_OKAY;
   }

   SCIPdebugMessage("propagating pseudo objective function (pseudoobj: %g, cutoffbound: %g)\n", pseudoobjval, cutoffbound);

   *result = SCIP_DIDNOTFIND;
   nchgbds = 0;
   
   /* tighten domains, if they would increase the pseudo objective value above the upper bound */
   ncands = (propdata->maxcands >= 0 ? MIN(propdata->maxcands, nobjvars) : nobjvars);
   v = propdata->lastvarnum;
   for( c = 0; c < ncands; ++c )
   {
      v++;
      if( v >= nobjvars )
         v = 0;

      SCIP_CALL( propagateCutoffboundVar(scip, prop, objvars[v], cutoffbound, pseudoobjval, &nchgbds) );
   }
   propdata->lastvarnum = v;

   if( nchgbds > 0 )
      *result = SCIP_REDUCEDDOM;
   
   return SCIP_OKAY;
}

/** recalculates the maximum objective pseudoactivity */
static
void calcMaxObjPseudoactivity(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata            /**< propagator data */
   )
{
   SCIP_VAR** vars;
   int nvars;
   int v;

   assert(propdata != NULL);

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* calculate current max pseudo activity and largest contribution */
   propdata->maxpseudoobjact = 0.0;
   propdata->maxpseudoobjactinf = 0;
   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real obj;
      SCIP_Real contrib;

      obj = SCIPvarGetObj(vars[v]);
      if( SCIPisPositive(scip, obj) )
      {
         contrib = SCIPvarGetUbGlobal(vars[v]);
         if( !SCIPisInfinity(scip, contrib) )
            contrib *= obj;
      }
      else if( SCIPisNegative(scip, obj) )
      {
         contrib = SCIPvarGetLbGlobal(vars[v]);
         if( !SCIPisInfinity(scip, contrib) )
            contrib *= obj;
         else
            contrib *= -1.0;
      }
      else
         continue;

      if( SCIPisInfinity(scip, contrib) )
         propdata->maxpseudoobjactinf++;
      else
         propdata->maxpseudoobjact += contrib;
   }
}

/** returns the residual pseudo objective activity without the given value */
static
SCIP_Real getMaxObjPseudoactivityResidualValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_Real             contrib             /**< value to eliminate from pseudo objective activity */
   )
{
   SCIP_Real residual;

   assert(propdata != NULL);

   /* if necessary, calculate the maximum pseudo objective activity */
   if( propdata->maxpseudoobjact == SCIP_INVALID ) /*lint !e777*/
      calcMaxObjPseudoactivity(scip, propdata);
   assert(propdata->maxpseudoobjact != SCIP_INVALID); /*lint !e777*/

   if( SCIPisInfinity(scip, contrib) )
   {
      assert(propdata->maxpseudoobjactinf >= 1);
      /* check if this variable yields the only infinite contribution */
      if( propdata->maxpseudoobjactinf == 1 )
         residual = propdata->maxpseudoobjact;
      else
         residual = SCIPinfinity(scip);
   }
   else
   {
      /* check if there is an infinite contribution */
      if( propdata->maxpseudoobjactinf >= 1 )
         residual = SCIPinfinity(scip);
      else
         residual = propdata->maxpseudoobjact - contrib;
   }

   return residual;
}

/** returns the residual pseudo objective activity */
static
SCIP_Real getMaxObjPseudoactivityResidual(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< propagator data */
   SCIP_VAR*             var                 /**< variable to get residual activity for */
   )
{
   SCIP_Real obj;
   SCIP_Real contrib;

   assert(propdata != NULL);

   contrib = 0.0;
   obj = SCIPvarGetObj(var);
   if( SCIPisPositive(scip, obj) )
   {
      contrib = SCIPvarGetUbGlobal(var);
      if( !SCIPisInfinity(scip, contrib) )
         contrib *= obj;
   }
   else if( SCIPisNegative(scip, obj) )
   {
      contrib = SCIPvarGetLbGlobal(var);
      if( !SCIPisInfinity(scip, -contrib) )
         contrib *= obj;
      else
         contrib *= -1.0;
   }
   
   return getMaxObjPseudoactivityResidualValue(scip, propdata, contrib);
}

/** propagates the global lower (dual) bound c*x >= lowerbound */
static
SCIP_RETCODE propagateLowerbound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   )

{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Real lowerbound;
   SCIP_Real residual;
   SCIP_VAR** objvars;
   int nobjvars;
   int k;

   assert(result != NULL);

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   objvars = propdata->objvars;
   nobjvars = propdata->nobjvars;
   assert(nobjvars == 0 || objvars != NULL);

   /* nothing to do for empty objective */
   if( nobjvars == 0 )
      return SCIP_OKAY;

   /* check if there is a chance to find a reduction */
   /**@todo This must be done in a clever way: Store the smallest global lower bound at which
    *       a propagation is triggered and skip the propagation if lowerbound < triggerlowerbound.
    *       The same can be done for the cutoff propagation, which most probably makes the lastvarnum approach irrelevant.
    */
   lowerbound = SCIPgetLowerbound(scip);

   if( SCIPisInfinity(scip, -lowerbound) )
      return SCIP_OKAY;

   /* propagate c*x >= lowerbound */
   for( k = 0; k < nobjvars; k++ )
   {
      SCIP_VAR* var;
      SCIP_Real obj;
      SCIP_Real lb;
      SCIP_Real ub;

      var = objvars[k];
      obj = SCIPvarGetObj(var);
      assert(!SCIPisZero(scip, obj));

      lb = SCIPvarGetLbGlobal(var);
      ub = SCIPvarGetUbGlobal(var);
      residual = getMaxObjPseudoactivityResidual(scip, propdata, var);
      if( SCIPisInfinity(scip, residual) )
         continue;

      if( obj > 0.0 )
      {
         SCIP_Real newlb;
         SCIP_Bool infeasible;
         SCIP_Bool tightened;

         newlb = (lowerbound - residual)/obj;

         if( SCIPisLbBetter(scip, newlb, lb, ub) )
         {
            SCIPdebugMessage(" -> new global lower bound of variable <%s>[%.10f,%.10f]: %.10f (obj=%g, residual=%g)\n",
               SCIPvarGetName(var), lb, ub, newlb, obj, residual);
            SCIP_CALL( SCIPtightenVarLbGlobal(scip, var, newlb, FALSE, &infeasible, &tightened) );
            if( infeasible )
            {
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            if( tightened ) /* might not be tightened due to numerical reasons */
               *result = SCIP_REDUCEDDOM;
         }
      }
      else
      {
         SCIP_Real newub;
         SCIP_Bool infeasible;
         SCIP_Bool tightened;

         newub = (lowerbound - residual)/obj;
         if( SCIPisUbBetter(scip, newub, lb, ub) )
         {
            SCIPdebugMessage(" -> new global upper bound of variable <%s>[%.10f,%.10f]: %.10f (obj=%g, residual=%g)\n",
               SCIPvarGetName(var), lb, ub, newub, obj, residual);
            SCIP_CALL( SCIPtightenVarUbGlobal(scip, var, newub, FALSE, &infeasible, &tightened) );
            if( infeasible )
            {
               *result = SCIP_CUTOFF;
               return SCIP_OKAY;
            }
            if( tightened ) /* might not be tightened due to numerical reasons */
               *result = SCIP_REDUCEDDOM;
         }
      }
   }

   return SCIP_OKAY;
}



/*
 * Callback methods of propagator
 */

/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyPseudoobj)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of propagator */
   SCIP_CALL( SCIPincludePropPseudoobj(scip) );
 
   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreePseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   SCIPfreeMemory(scip, &propdata);
   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
#define propInitPseudoobj NULL


/** deinitialization method of propagator (called before transformed problem is freed) */
#define propExitPseudoobj NULL


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
static
SCIP_DECL_PROPINITSOL(propInitsolPseudoobj)
{
   SCIP_PROPDATA* propdata;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_VAR** vars;
   int nvars;
   int nobjvars;
   int v;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   eventhdlr = propdata->eventhdlr;

   /* get problem variables */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* count variables with non-zero objective value */
   nobjvars = 0;
   for( v = 0; v < nvars; v++ )
   {
      SCIP_Real obj;

      obj = SCIPvarGetObj(vars[v]);
      if( !SCIPisZero(scip, obj) )
         nobjvars++;
   }
   SCIPdebugMessage("There are %d variables in the objective function.\n",nobjvars);

   /* collect the variables with non-zero objective and catch global bound tighten events that decrease the
    * maximal pseudo objective activity
    */
   if( nobjvars > 0 )
   {
      int k;

      /* allocate memory for non-zero objective variables */
      SCIP_CALL( SCIPallocMemoryArray(scip, &propdata->objvars, nobjvars) );
      SCIPdebugMessage("Store objective variables at address <%p>.\n", (void*)propdata->objvars);

      k = 0;
      for( v = 0; v < nvars; v++ )
      {
         SCIP_Real obj;
         SCIP_VAR* var;

         var = vars[v];
         obj = SCIPvarGetObj(var);
         if( SCIPisZero(scip, obj) )
            continue;

         /* store and capture variable */
         assert(k < nobjvars);
         propdata->objvars[k] = var;
         SCIP_CALL( SCIPcaptureVar(scip, var) ) ;
         k++;

         /* catch bound change events */
         if( obj > 0.0 )
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, var,
                  SCIP_EVENTTYPE_GUBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
         }
         else
         {
            SCIP_CALL( SCIPcatchVarEvent(scip, var,
                  SCIP_EVENTTYPE_GLBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, NULL) );
         }
      }
      assert(k == nobjvars);
   }

   propdata->nobjvars = nobjvars;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;

   return SCIP_OKAY;
}

/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
static
SCIP_DECL_PROPEXITSOL(propExitsolPseudoobj)
{
   SCIP_PROPDATA* propdata;
   SCIP_EVENTHDLR* eventhdlr;
   SCIP_VAR** objvars;
   int nobjvars;
   int k;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   eventhdlr = propdata->eventhdlr;
   objvars = propdata->objvars;
   nobjvars = propdata->nobjvars;
   assert(nobjvars == 0 || objvars != NULL);

   /* drop global bound tighten events that decrease the maximal pseudo objective activity and release variables */
   for( k = 0; k < nobjvars; k++ )
   {
      SCIP_VAR* var;
      SCIP_Real obj;

      var = objvars[k];
      obj = SCIPvarGetObj(var);
      assert(!SCIPisZero(scip, obj));

      /* drop bound change event */
      if( obj > 0.0 )
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, 
               SCIP_EVENTTYPE_GUBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }
      else
      {
         SCIP_CALL( SCIPdropVarEvent(scip, var, 
               SCIP_EVENTTYPE_GLBCHANGED, eventhdlr, (SCIP_EVENTDATA*)propdata, -1) );
      }

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &objvars[k]) );
   }

   propdata->nobjvars = 0;
   /* free memory for non-zero objective variables */
   SCIPfreeMemoryArrayNull(scip, &propdata->objvars);

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecPseudoobj)
{  /*lint --e{715}*/

   *result = SCIP_DIDNOTRUN;
   
   /* propagate c*x <= cutoff */
   SCIP_CALL( propagateCutoffbound(scip, prop, result) );

   if( *result != SCIP_CUTOFF && SCIPgetStage(scip) == SCIP_STAGE_SOLVING )
   {
      SCIP_RESULT dualresult;

      /* propagate c*x >= lowerbound */
      SCIP_CALL( propagateLowerbound(scip, prop, &dualresult) );

      if( dualresult == SCIP_REDUCEDDOM || dualresult == SCIP_CUTOFF )
         *result = dualresult;
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIP_CALL( resolvePropagation(scip, propdata, infervar, bdchgidx) );

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}




/*
 * Event handler
 */

/** execution methode of bound change event handler */
static
SCIP_DECL_EVENTEXEC(eventExecPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   propdata = (SCIP_PROPDATA*)eventdata;
   assert(propdata != NULL);

   /* invalidate the maximum pseudo objective activity */
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;

   return SCIP_OKAY;
}




/*
 * propagator specific interface methods
 */

/** creates the pseudo objective function propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropPseudoobj(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;

   /* include event handler for bound change events */
   SCIP_CALL( SCIPincludeEventhdlr(scip, EVENTHDLR_NAME, EVENTHDLR_DESC,
         NULL,
         NULL, NULL, NULL, NULL, NULL, NULL, eventExecPseudoobj, NULL) );

   /* create pseudoobj propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->objvars = NULL;
   propdata->nobjvars = 0;
   propdata->maxpseudoobjact = SCIP_INVALID;
   propdata->maxpseudoobjactinf = 0;
   propdata->lastvarnum = -1;

   /* get event handler for bound change events */
   propdata->eventhdlr = SCIPfindEventhdlr(scip, EVENTHDLR_NAME);
   if( propdata->eventhdlr == NULL )
   {
      SCIPerrorMessage("event handler for pseudo objective propagator not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY,
         propCopyPseudoobj,
         propFreePseudoobj, propInitPseudoobj, propExitPseudoobj, 
         propInitsolPseudoobj, propExitsolPseudoobj, propExecPseudoobj, propRespropPseudoobj,
         propdata) );

   /* add pseudoobj propagator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/pseudoobj/maxcands", 
         "maximal number of variables to look at in a single propagation round (-1: process all variables)",
         &propdata->maxcands, TRUE, DEFAULT_MAXCANDS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}

/** propagates the cutoff bound for the given variables */
SCIP_RETCODE SCIPpropagateCutoffboundVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROP*            prop,               /**< propagator, or NULL */
   SCIP_VAR*             var,                /**< variables to propagate */
   SCIP_Real             cutoffbound,        /**< cutoff bound to use */
   SCIP_Real             pseudoobjval,       /**< pseudo objective value to use */
   int*                  nchgbds             /**< pointer to store the number of changed bounds */
   )
{
   SCIP_CALL( propagateCutoffboundVar(scip, prop, var, cutoffbound, pseudoobjval, nchgbds) );

   return SCIP_OKAY;
}
