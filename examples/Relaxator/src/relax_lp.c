/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_lp.c
 * @brief  lp relaxator
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "relax_lp.h"

#define RELAX_NAME             "lp"
#define RELAX_DESC             "relaxator solving LP relaxation"
#define RELAX_PRIORITY         0
#define RELAX_FREQ             0


/*
 * Data structures
 */


/*
 * Local methods
 */


/*
 * Callback methods of relaxator
 */

/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecLp)
{  /*lint --e{715}*/
   SCIP* relaxscip;
   SCIP_HASHMAP* varmap;
   SCIP_CONS** conss;
   SCIP_Real relaxval;
   SCIP_Bool valid;
   int nconss;
   int i;
   int c;

   *lowerbound = -SCIPinfinity(scip);
   *result = SCIP_DIDNOTRUN;

   /* we can only run if none of the present constraints expect their variables to be binary or integer during transformation */
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   for( c = 0; c < nconss; ++c )
   {
      const char* conshdlrname;

      conshdlrname = SCIPconshdlrGetName(SCIPconsGetHdlr(conss[c]));

      /* skip if there are any "and", "linking", or", "orbitope", "pseudoboolean", "superindicator", "xor" or new/unknown constraints */
      if( strcmp(conshdlrname, "SOS1") != 0 && strcmp(conshdlrname, "SOS2") != 0 && strcmp(conshdlrname, "abspower") != 0
            && strcmp(conshdlrname, "bivariate") != 0 && strcmp(conshdlrname, "bounddisjunction") != 0
            && strcmp(conshdlrname, "cardinality") != 0 && strcmp(conshdlrname, "components") != 0
            && strcmp(conshdlrname, "conjunction") != 0 && strcmp(conshdlrname, "countsols") != 0
            && strcmp(conshdlrname, "cumulative") != 0 && strcmp(conshdlrname, "disjunction") != 0
            && strcmp(conshdlrname, "indicator") != 0 && strcmp(conshdlrname, "integral") != 0
            && strcmp(conshdlrname, "knapsack") != 0 && strcmp(conshdlrname, "linear") != 0
            && strcmp(conshdlrname, "logicor") != 0 && strcmp(conshdlrname, "nonlinear") != 0
            && strcmp(conshdlrname, "orbisack") != 0 && strcmp(conshdlrname, "quadratic") != 0
            && strcmp(conshdlrname, "setppc") != 0 && strcmp(conshdlrname, "soc") != 0
            && strcmp(conshdlrname, "symresack") != 0 && strcmp(conshdlrname, "varbound") != 0 )
         return SCIP_OKAY;
   }

   /* create the variable mapping hash map */
   SCIP_CALL( SCIPcreate(&relaxscip) );
   SCIP_CALL( SCIPhashmapCreate(&varmap, SCIPblkmem(relaxscip), SCIPgetNVars(scip)) );
   valid = FALSE;
   SCIP_CALL( SCIPcopy(scip, relaxscip, varmap, NULL, "relaxscip", FALSE, FALSE, FALSE, &valid) );

   /* change variable types */
   for( i = 0; i < SCIPgetNVars(relaxscip); ++i )
   {
      SCIP_VAR* var;
      SCIP_Bool infeasible;

      var = SCIPgetVars(relaxscip)[i];
      assert(var != NULL);

      SCIP_CALL( SCIPchgVarType(relaxscip, var, SCIP_VARTYPE_CONTINUOUS, &infeasible) );
      assert(!infeasible);
   }

   SCIPsetMessagehdlrQuiet(relaxscip, TRUE);
   SCIP_CALL( SCIPtransformProb(relaxscip) );
   SCIP_CALL( SCIPsolve(relaxscip) );
   relaxval = SCIPgetPrimalbound(relaxscip);
   SCIPdebugMessage("relaxation bound = %e status = %d\n", relaxval, SCIPgetStatus(relaxscip));

   if( SCIPgetStatus(relaxscip) == SCIP_STATUS_OPTIMAL )
   {
      /* store relaxation solution in original SCIP if it improves the best relaxation solution thus far */
      if( (! SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)) )
      {
         SCIPdebugMsg(scip, "Setting LP relaxation solution, which improved upon earlier solution\n");
         SCIP_CALL( SCIPclearRelaxSolVals(scip) );

         for( i = 0; i < SCIPgetNVars(scip); ++i )
         {
            SCIP_VAR* relaxvar;
            SCIP_Real solval;

            relaxvar = SCIPhashmapGetImage(varmap, SCIPgetVars(scip)[i]);
            assert(relaxvar != NULL);

            solval = SCIPgetSolVal(relaxscip, SCIPgetBestSol(relaxscip), relaxvar);

            SCIP_CALL( SCIPsetRelaxSolVal(scip, SCIPgetVars(scip)[i], solval) );
         }

         /* mark relaxation solution to be valid and inform SCIP that the relaxation included all LP rows */
         SCIP_CALL( SCIPmarkRelaxSolValid(scip, TRUE) );
      }

      SCIPdebugMsg(scip, "LP lower bound = %g\n", relaxval);
      *lowerbound = relaxval;
      *result = SCIP_SUCCESS;
   }

   /* free memory */
   SCIPhashmapFree(&varmap);
   SCIP_CALL( SCIPfree(&relaxscip) );

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods
 */

/** creates the lp relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxLp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   /* create lp relaxator data */
   relaxdata = NULL;
   relax = NULL;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecLp, relaxdata) );
   assert(relax != NULL);

   return SCIP_OKAY;
}
