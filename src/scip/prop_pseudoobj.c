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
#pragma ident "@(#) $Id: prop_pseudoobj.c,v 1.1 2004/09/23 15:46:31 bzfpfend Exp $"

/**@file   prop_pseudoobj.c
 * @brief  pseudoobj propagator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "prop_pseudoobj.h"


#define PROP_NAME              "pseudoobj"
#define PROP_DESC              "pseudo objective function propagator"
#define PROP_PRIORITY                 0
#define PROP_FREQ                     1




/*
 * Data structures
 */

/* TODO: fill in the necessary propagator data */

/** propagator data */
struct PropData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of propagator
 */

/** destructor of propagator to free user data (called when SCIP is exiting) */
#define propFreePseudoobj NULL


/** initialization method of propagator (called after problem was transformed) */
#define propInitPseudoobj NULL


/** deinitialization method of propagator (called before transformed problem is freed) */
#define propExitPseudoobj NULL


/** execution method of propagator */
static
DECL_PROPEXEC(propExecPseudoobj)
{  /*lint --e{715}*/
   VAR** vars;
   VAR* var;
   Real pseudoobjval;
   Real cutoffbound;
   Real obj;
   Real lb;
   Real ub;
   int nvars;
   int v;

   *result = SCIP_DIDNOTRUN;

   /* get current pseudo objective value and cutoff bound */
   pseudoobjval = SCIPgetPseudoObjval(scip);
   if( SCIPisInfinity(scip, -pseudoobjval) )
      return SCIP_OKAY;
   cutoffbound = SCIPgetCutoffbound(scip);
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;
   if( SCIPisGE(scip, pseudoobjval, cutoffbound) )
   {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
   }

   debugMessage("propagating pseudo objective function\n");

   *result = SCIP_DIDNOTFIND;

   /* tighten domains, if they would increase the pseudo objective value above the upper bound */
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   for( v = 0; v < nvars; ++v )
   {
      var = vars[v];
      lb = SCIPvarGetLbLocal(var);
      ub = SCIPvarGetUbLocal(var);
      if( SCIPisFeasEQ(scip, lb, ub) )
         continue;
      obj = SCIPvarGetObj(var);

      if( SCIPisPositive(scip, obj) )
      {
         Real newub;
         Bool infeasible;
         Bool tightened;

         newub = lb + (cutoffbound - pseudoobjval)/obj;
         if( SCIPisUbBetter(scip, newub, ub) )
         {
            debugMessage(" -> new upper bound of variable <%s>[%g,%g]: %g\n", SCIPvarGetName(var), lb, ub, newub);
            CHECK_OKAY( SCIPinferVarUbProp(scip, var, newub, prop, 0, &infeasible, &tightened) );
            assert(!infeasible);
            assert(tightened);
            *result = SCIP_REDUCEDDOM;
         }
      }
      else if( SCIPisNegative(scip, obj) )
      {
         Real newlb;
         Bool infeasible;
         Bool tightened;

         newlb = ub + (cutoffbound - pseudoobjval)/obj;
         if( SCIPisLbBetter(scip, newlb, lb) )
         {
            debugMessage(" -> new lower bound of variable <%s>[%g,%g]: %g\n", SCIPvarGetName(var), lb, ub, newlb);
            CHECK_OKAY( SCIPinferVarLbProp(scip, var, newlb, prop, 0, &infeasible, &tightened) );
            assert(!infeasible);
            assert(tightened);
            *result = SCIP_REDUCEDDOM;
         }
      }
   }

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
DECL_PROPRESPROP(propRespropPseudoobj)
{  /*lint --e{715}*/
   VAR** vars;
   VAR* var;
   Real obj;
   int nvars;
   int v;

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
         Real loclb;
         Real glblb;

         glblb = SCIPvarGetLbGlobal(var);
         loclb = SCIPvarGetLbAtIndex(var, bdchgidx, FALSE);
         if( SCIPisGT(scip, loclb, glblb) )
         {
            CHECK_OKAY( SCIPaddConflictLb(scip, var, bdchgidx) );
         }
      }
      else if( SCIPisNegative(scip, obj) )
      {
         Real locub;
         Real glbub;

         glbub = SCIPvarGetUbGlobal(var);
         locub = SCIPvarGetUbAtIndex(var, bdchgidx, FALSE);
         if( SCIPisLT(scip, locub, glbub) )
         {
            CHECK_OKAY( SCIPaddConflictUb(scip, var, bdchgidx) );
         }
      }
   }

   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}




/*
 * propagator specific interface methods
 */

/** creates the pseudo objective function propagator and includes it in SCIP */
RETCODE SCIPincludePropPseudoobj(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   PROPDATA* propdata;

   /* create pseudoobj propagator data */
   propdata = NULL;
   /* TODO: (optional) create propagator specific data here */

   /* include propagator */
   CHECK_OKAY( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ,
         propFreePseudoobj, propInitPseudoobj, propExitPseudoobj, propExecPseudoobj, propRespropPseudoobj,
         propdata) );

   /* add pseudoobj propagator parameters */
   /* TODO: (optional) add propagator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
