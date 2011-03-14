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

/**@file   sepa_redcost.c
 * @ingroup SEPARATORS
 * @brief  reduced cost strengthening separator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/sepa_redcost.h"


#define SEPA_NAME              "redcost"
#define SEPA_DESC              "reduced cost strengthening separator"
#define SEPA_PRIORITY         +10000000
#define SEPA_FREQ                     1
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */

#define DEFAULT_CONTINUOUS        FALSE /**< should reduced cost fixing be also applied to continuous variables? */


/** separator data */
struct SCIP_SepaData
{
   SCIP_Bool             continuous;         /**< should reduced cost fixing be also applied to continuous variables? */
};



/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyRedcost)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaRedcost(scip) );
 
   return SCIP_OKAY;
}

/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeRedcost)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   /* free separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeMemory(scip, &sepadata);

   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
#define sepaInitRedcost NULL


/** deinitialization method of separator (called before transformed problem is freed) */
#define sepaExitRedcost NULL


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#define sepaInitsolRedcost NULL


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#define sepaExitsolRedcost NULL


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpRedcost)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;
   SCIP_COL** cols;
   SCIP_Real cutoffbound;
   SCIP_Real lpobjval;
   int ncols;
   int c;

   *result = SCIP_DIDNOTRUN;

   /* we cannot apply reduced cost fixing, if we want to solve exactly */
   /**@todo implement reduced cost fixing with interval arithmetics */
   if( SCIPisExactSolve(scip) )
      return SCIP_OKAY;

   /* only call separator, if an optimal LP solution is at hand */
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

   /* only call separator, if the current LP is a valid relaxation */
   if( !SCIPisLPRelax(scip) )
      return SCIP_OKAY;

   /* get LP columns */
   cols = SCIPgetLPCols(scip);
   ncols = SCIPgetNLPCols(scip);
   if( ncols == 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* get separator data */
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   /* get LP objective value */
   lpobjval = SCIPgetLPObjval(scip);

   /* check reduced costs for non-basic columns */
   for( c = 0; c < ncols; ++c )
   {
      SCIP_VAR* var;
      SCIP_Real redcost;

      var = SCIPcolGetVar(cols[c]);
      if( !sepadata->continuous && !SCIPvarIsIntegral(var) )
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


/** arbitrary primal solution separation method of separator */
#define sepaExecsolRedcost NULL




/*
 * separator specific interface methods
 */

/** creates the reduced cost strengthening separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaRedcost(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create redcost separator data */
   SCIP_CALL( SCIPallocMemory(scip, &sepadata) );

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY, sepaCopyRedcost,
         sepaFreeRedcost, sepaInitRedcost, sepaExitRedcost, 
         sepaInitsolRedcost, sepaExitsolRedcost,
         sepaExeclpRedcost, sepaExecsolRedcost,
         sepadata) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "separating/redcost/continuous",
         "should reduced cost fixing be also applied to continuous variables?",
         &sepadata->continuous, FALSE, DEFAULT_CONTINUOUS, NULL, NULL) );

   return SCIP_OKAY;
}
