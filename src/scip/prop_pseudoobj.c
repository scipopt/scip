/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: prop_pseudoobj.c,v 1.17 2007/03/12 14:50:20 bzfpfend Exp $"

/**@file   prop_pseudoobj.c
 * @brief  pseudoobj propagator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/prop_pseudoobj.h"


#define PROP_NAME              "pseudoobj"
#define PROP_DESC              "pseudo objective function propagator"
#define PROP_PRIORITY                 0
#define PROP_FREQ                     1
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */

#define DEFAULT_MAXCANDS            100 /**< maximal number of variables to look at in a single propagation round
                                         *   (-1: process all variables) */




/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
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
   SCIP_VAR*             infervar,           /**< variable that was deduced, or NULL for conflict analysis initialization */
   SCIP_BDCHGIDX*        bdchgidx            /**< bound change index (time stamp of bound change), or NULL for current time */
   )
{
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real obj;
   int nvars;
   int v;

   /**@todo improve pseudo objective propagator conflict resolving method:
    *       only add bound changes up to the point, the primal bound is reached
    */

   /* the variables responsible for the propagation are the ones with
    *  - obj > 0 and local lb > global lb
    *  - obj < 0 and local ub < global ub
    */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      if( var == infervar )
         continue;

      obj = SCIPvarGetObj(var);
      if( SCIPisPositive(scip, obj) )
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
      else if( SCIPisNegative(scip, obj) )
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




/*
 * Callback methods of propagator
 */

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
#define propInitsolPseudoobj NULL


/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
#define propExitsolPseudoobj NULL


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecPseudoobj)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_VAR** vars;
   SCIP_VAR* var;
   SCIP_Real pseudoobjval;
   SCIP_Real cutoffbound;
   SCIP_Real obj;
   SCIP_Real lb;
   SCIP_Real ub;
   int ncands;
   int nvars;
   int c;
   int v;

   *result = SCIP_DIDNOTRUN;

   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

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
      SCIP_CALL( resolvePropagation(scip, NULL, NULL) );

      /* analyze the conflict */
      SCIP_CALL( SCIPanalyzeConflict(scip, 0, NULL) );
      
      *result = SCIP_CUTOFF;

      return SCIP_OKAY;
   }

   SCIPdebugMessage("propagating pseudo objective function (pseudoobj: %g, cutoffbound: %g)\n", pseudoobjval, cutoffbound);

   *result = SCIP_DIDNOTFIND;

   /* tighten domains, if they would increase the pseudo objective value above the upper bound */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   ncands = (propdata->maxcands >= 0 ? MIN(propdata->maxcands, nvars) : nvars);
   v = propdata->lastvarnum;
   for( c = 0; c < ncands; ++c )
   {
      v++;
      if( v >= nvars )
         v = 0;
      var = vars[v];
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      if( SCIPisFeasEQ(scip, lb, ub) )
         continue;
      obj = SCIPvarGetObj(var);

      if( SCIPisPositive(scip, obj) )
      {
         SCIP_Real newub;
         SCIP_Bool infeasible;
         SCIP_Bool tightened;

         newub = lb + (cutoffbound - pseudoobjval)/obj;
         if( SCIPisUbBetter(scip, newub, lb, ub) )
         {
            SCIPdebugMessage(" -> new upper bound of variable <%s>[%.10f,%.10f]: %.10f\n",
               SCIPvarGetName(var), lb, ub, newub);
            SCIP_CALL( SCIPinferVarUbProp(scip, var, newub, prop, 0, &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened ) /* might not be tightened due to numerical reasons */
               *result = SCIP_REDUCEDDOM;
         }
      }
      else if( SCIPisNegative(scip, obj) )
      {
         SCIP_Real newlb;
         SCIP_Bool infeasible;
         SCIP_Bool tightened;

         newlb = ub + (cutoffbound - pseudoobjval)/obj;
         if( SCIPisLbBetter(scip, newlb, lb, ub) )
         {
            SCIPdebugMessage(" -> new lower bound of variable <%s>[%g,%g]: %g\n", SCIPvarGetName(var), lb, ub, newlb);
            SCIP_CALL( SCIPinferVarLbProp(scip, var, newlb, prop, 0, &infeasible, &tightened) );
            assert(!infeasible);
            if( tightened ) /* might not be tightened due to numerical reasons */
               *result = SCIP_REDUCEDDOM;
         }
      }
   }
   propdata->lastvarnum = v;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropPseudoobj)
{  /*lint --e{715}*/
   SCIP_CALL( resolvePropagation(scip, infervar, bdchgidx) );

   *result = SCIP_SUCCESS;

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

   /* create pseudoobj propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );
   propdata->lastvarnum = -1;

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY,
         propFreePseudoobj, propInitPseudoobj, propExitPseudoobj, 
         propInitsolPseudoobj, propExitsolPseudoobj, propExecPseudoobj, propRespropPseudoobj,
         propdata) );

   /* add pseudoobj propagator parameters */
   SCIP_CALL( SCIPaddIntParam(scip,
         "propagating/pseudoobj/maxcands", 
         "maximal number of variables to look at in a single propagation round (-1: process all variables)",
         &propdata->maxcands, DEFAULT_MAXCANDS, -1, INT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
