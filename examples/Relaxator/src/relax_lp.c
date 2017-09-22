/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
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

#include "relax_lp.h"

#define RELAX_NAME             "lp"
#define RELAX_DESC             "relaxator solving LP relaxation"
#define RELAX_PRIORITY         0
#define RELAX_FREQ             0
#define RELAX_FULLLPINFO       TRUE


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
   SCIP_Real relaxval;
   SCIP_Bool valid;
   int i;

   *lowerbound = -SCIPinfinity(scip);
   *result = SCIP_DIDNOTRUN;

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
      *lowerbound = relaxval;
      *result = SCIP_SUCCESS;

      /* store relaxation solution in original SCIP */
      for( i = 0; i < SCIPgetNVars(scip); ++i )
      {
         SCIP_VAR* relaxvar;
         SCIP_Real solval;

         relaxvar = SCIPhashmapGetImage(varmap, SCIPgetVars(scip)[i]);
         assert(relaxvar != NULL);

         solval = SCIPgetSolVal(relaxscip, SCIPgetBestSol(relaxscip), relaxvar);

         SCIP_CALL( SCIPsetRelaxSolVal(scip, SCIPgetVars(scip)[i], solval) );
      }

      /* mark relaxation solution to be valid */
      SCIP_CALL( SCIPmarkRelaxSolValid(scip) );
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
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ, RELAX_FULLLPINFO,
         relaxExecLp, relaxdata) );
   assert(relax != NULL);

   return SCIP_OKAY;
}
