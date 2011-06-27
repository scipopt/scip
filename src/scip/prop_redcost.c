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

/**@file   prop_redcost.c
 * @ingroup PROPAGATORS
 * @brief  redcost propagator
 * @author Tobias Achterberg
 * @author Matthias Miltenberger
 * @author Michael Winkler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/prop_redcost.h"


#define PROP_NAME              "redcost"
#define PROP_DESC              "reduced cost strengthening propagator"
#define PROP_TIMING             SCIP_PROPTIMING_DURINGLPLOOP
#define PROP_PRIORITY         +10000000 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_PRESOL_PRIORITY          0 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOL_DELAY          TRUE /**< should presolving be delay, if other presolvers found reductions?  */
#define PROP_PRESOL_MAXROUNDS         0 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */
#define MAXBOUNDDIST                1.0
#define DEFAULT_CONTINUOUS        FALSE /**< should reduced cost fixing be also applied to continuous variables? */




/*
 * Data structures
 */


/** propagator data */
struct SCIP_PropData
{
   SCIP_Bool             continuous;         /**< should reduced cost fixing be also applied to continuous variables? */
};


/** copy method for propagator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_PROPCOPY(propCopyRedcost)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(prop != NULL);
   assert(strcmp(SCIPpropGetName(prop), PROP_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludePropRedcost(scip) );

   return SCIP_OKAY;
}

/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeRedcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   /* free propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   SCIPfreeMemory(scip, &propdata);

   SCIPpropSetData(prop, NULL);

   return SCIP_OKAY;
}


/** initialization method of propagator (called after problem was transformed) */
#define propInitRedcost NULL


/** deinitialization method of propagator (called before transformed problem is freed) */
#define propExitRedcost NULL


/** presolving initialization method of propagator (called when presolving is about to begin) */
#define propInitpreRedcost NULL


/** presolving deinitialization method of propagator (called after presolving has been finished) */
#define propExitpreRedcost NULL


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
#define propInitsolRedcost NULL


/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
#define propExitsolRedcost NULL


/** presolving method of propagator */
#define propPresolRedcost NULL


/** reduced cost propagation method for an LP solution */
static
SCIP_DECL_PROPEXEC(propExecRedcost)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_COL** cols;
   SCIP_Real cutoffbound;
   SCIP_Real lpobjval;
   int ncols;
   int c;

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* we cannot apply reduced cost fixing, if we want to solve exactly */
   /**@todo implement reduced cost fixing with interval arithmetics */
   if( SCIPisExactSolve(scip) )
      return SCIP_OKAY;

   /* only call propagator, if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* we cannot apply reduced cost strengthening, if no simplex basis is available */
   if( !SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* get current cutoff bound */
   cutoffbound = SCIPgetCutoffbound(scip);

   /* reduced cost strengthening can only be applied, if we have a finite upper bound on the LP value */
   if( SCIPisInfinity(scip, cutoffbound) )
      return SCIP_OKAY;

   /* only call propagator, if the current LP is a valid relaxation */
   if( !SCIPisLPRelax(scip) )
      return SCIP_OKAY;

   /* get LP columns */
   cols = SCIPgetLPCols(scip);
   ncols = SCIPgetNLPCols(scip);
   if( ncols == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get propagator data */
   propdata = SCIPpropGetData(prop);
   assert(propdata != NULL);

   /* get LP objective value */
   lpobjval = SCIPgetLPObjval(scip);

   /* check reduced costs for non-basic columns */
   for( c = 0; c < ncols; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real redcost;

      var = SCIPcolGetVar(cols[c]);
      if( !propdata->continuous && !SCIPvarIsIntegral(var) )
         continue;

      switch( SCIPcolGetBasisStatus(cols[c]) )
      {
      case SCIP_BASESTAT_LOWER:
         redcost = SCIPgetColRedcost(scip, cols[c]);
         assert( !SCIPisFeasNegative(scip, redcost) || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );
         if( SCIPisFeasPositive(scip, redcost) )
         {
            SCIP_Real oldlb;
            SCIP_Real oldub;

            oldlb = SCIPvarGetLbLocal(var);
            oldub = SCIPvarGetUbLocal(var);
            assert(SCIPisEQ(scip, oldlb, SCIPcolGetLb(cols[c])));
            assert(SCIPisEQ(scip, oldub, SCIPcolGetUb(cols[c])));
            if( SCIPisFeasLT(scip, oldlb, oldub) )
            {
               SCIP_Real newub;
               SCIP_Bool strengthen;

               /* calculate reduced cost based bound */
               newub = (cutoffbound - lpobjval) / redcost + oldlb;

               /* check, if new bound is good enough:
                *  - integer variables: take all possible strengthenings
                *  - continuous variables: strengthening must cut part of the variable's dynamic range, and
                *                          at least 20% of the current domain
                */
               if( SCIPvarIsIntegral(var) )
               {
                  newub = SCIPadjustedVarUb(scip, var, newub);
                  strengthen = (newub < oldub - 0.5);
               }
               else
                  strengthen = (newub < SCIPcolGetMaxPrimsol(cols[c]) && newub <= 0.2 * oldlb + 0.8 * oldub);

               if( strengthen )
               {
                  /* strengthen upper bound */
                  SCIPdebugMessage("redcost strengthening upper bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                     SCIPvarGetName(var), oldlb, oldub, oldlb, newub, cutoffbound, lpobjval, redcost);
                  SCIP_CALL( SCIPchgVarUb(scip, var, newub) );
                  *result = SCIP_REDUCEDDOM;
               }
            }
         }
         break;

      case SCIP_BASESTAT_BASIC:
         break;

      case SCIP_BASESTAT_UPPER:
         redcost = SCIPgetColRedcost(scip, cols[c]);
         assert( !SCIPisFeasPositive(scip, redcost) || SCIPisFeasEQ(scip, SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var)) );
         if( SCIPisFeasNegative(scip, redcost) )
         {
            SCIP_Real oldlb;
            SCIP_Real oldub;

            oldlb = SCIPvarGetLbLocal(var);
            oldub = SCIPvarGetUbLocal(var);
            assert(SCIPisEQ(scip, oldlb, SCIPcolGetLb(cols[c])));
            assert(SCIPisEQ(scip, oldub, SCIPcolGetUb(cols[c])));
            if( SCIPisFeasLT(scip, oldlb, oldub) )
            {
               SCIP_Real newlb;
               SCIP_Bool strengthen;

               /* calculate reduced cost based bound */
               newlb = (cutoffbound - lpobjval) / redcost + oldub;

               /* check, if new bound is good enough:
                *  - integer variables: take all possible strengthenings
                *  - continuous variables: strengthening must cut part of the variable's dynamic range, and
                *                          at least 20% of the current domain
                */
               if( SCIPvarIsIntegral(var) )
               {
                  newlb = SCIPadjustedVarLb(scip, var, newlb);
                  strengthen = (newlb > oldlb + 0.5);
               }
               else
                  strengthen = (newlb > SCIPcolGetMinPrimsol(cols[c]) && newlb >= 0.8 * oldlb + 0.2 * oldub);

               /* check, if new bound is good enough: at least 20% strengthening for continuous variables */
               if( strengthen )
               {
                  /* strengthen lower bound */
                  SCIPdebugMessage("redcost strengthening lower bound: <%s> [%g,%g] -> [%g,%g] (ub=%g, lb=%g, redcost=%g)\n",
                     SCIPvarGetName(var), oldlb, oldub, newlb, oldub, cutoffbound, lpobjval, redcost);
                  SCIP_CALL( SCIPchgVarLb(scip, var, newlb) );
                  *result = SCIP_REDUCEDDOM;
               }
            }
         }
         break;

      case SCIP_BASESTAT_ZERO:
         assert(SCIPisFeasZero(scip, SCIPgetColRedcost(scip, cols[c])));
         break;

      default:
         SCIPerrorMessage("invalid basis state\n");
         return SCIP_INVALIDDATA;
      }
   }

   return SCIP_OKAY;
}



/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropRedcost)
{  /*lint --e{715}*/
   *result = SCIP_DIDNOTRUN;

   return SCIP_OKAY;
}




/*
 * propagator specific interface methods
 */

/** creates the redcost propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropRedcost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;

   /* create redcost propagator data */
   SCIP_CALL( SCIPallocMemory(scip, &propdata) );


   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING, PROP_PRESOL_PRIORITY, PROP_PRESOL_MAXROUNDS, PROP_PRESOL_DELAY,
         propCopyRedcost,
         propFreeRedcost, propInitRedcost, propExitRedcost, propInitpreRedcost, propExitpreRedcost,
         propInitsolRedcost, propExitsolRedcost, propPresolRedcost, propExecRedcost, propRespropRedcost,
         propdata) );

   /* add redcost propagator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "propagating/redcost/continuous",
         "should reduced cost fixing be also applied to continuous variables?",
         &propdata->continuous, FALSE, DEFAULT_CONTINUOUS, NULL, NULL) );
   return SCIP_OKAY;
}
