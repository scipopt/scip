/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   relax_nlp.c
 * @brief  nlp relaxator
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "relax_nlp.h"
#include "scip/scip_nlpi.h"


#define RELAX_NAME             "nlp"
#define RELAX_DESC             "relaxator solving a convex NLP relaxation"
#define RELAX_PRIORITY         10
#define RELAX_FREQ             1

#define NLPITERLIMIT           500       /**< iteration limit of NLP solver */
#define FEASTOLFAC             0.01      /**< factor for NLP feasibility tolerance */
#define RELOBJTOLFAC           0.01      /**< factor for NLP relative objective tolerance */

/*
 * Data structures
 */


/*
 * Local methods
 */


/*
 * Callback methods of relaxator
 */


/** solving process initialization method of relaxator (called when branch and bound process is about to begin) */
static
SCIP_DECL_RELAXINITSOL(relaxInitsolNlp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** solving process deinitialization method of relaxator (called before branch and bound process data is freed) */
static
SCIP_DECL_RELAXEXITSOL(relaxExitsolNlp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** execution method of relaxator */
static
SCIP_DECL_RELAXEXEC(relaxExecNlp)
{  /*lint --e{715}*/
   SCIP_NLROW** nlrows;
   SCIP_NLPIPROBLEM* nlpiprob;
   SCIP_HASHMAP* var2idx;
   SCIP_NLPI* nlpi;
   SCIP_NLPPARAM nlpparam = SCIP_NLPPARAM_DEFAULT(scip);
   int nnlrows;

   *result = SCIP_DIDNOTRUN;
   *lowerbound = -SCIPinfinity(scip);

   /* check if it is not possible to run the relaxator */
   if( !SCIPisNLPConstructed(scip) || SCIPinProbing(scip) || SCIPinDive(scip) || !SCIPallColsInLP(scip) || (SCIPgetNNlpis(scip) == 0) )
      return SCIP_OKAY;

   nlrows = SCIPgetNLPNlRows(scip);
   nnlrows = SCIPgetNNLPNlRows(scip);

   /* create a convex NLP relaxation */
   nlpi = SCIPgetNlpis(scip)[0];
   assert(nlpi != NULL);

   SCIP_CALL( SCIPhashmapCreate(&var2idx, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   SCIP_CALL( SCIPcreateNlpiProblemFromNlRows(scip, nlpi, &nlpiprob, "relax-NLP", nlrows, nnlrows, var2idx, NULL, NULL, SCIPgetCutoffbound(scip),
         TRUE, TRUE) );
   SCIP_CALL( SCIPaddNlpiProblemRows(scip, nlpi, nlpiprob, var2idx, SCIPgetLPRows(scip), SCIPgetNLPRows(scip)) );

   nlpparam.iterlimit = NLPITERLIMIT;
   nlpparam.feastol = SCIPfeastol(scip) * FEASTOLFAC;
   nlpparam.opttol = SCIPfeastol(scip) * RELOBJTOLFAC;

   /* solve NLP */
   SCIP_CALL( SCIPsolveNlpiParam(scip, nlpi, nlpiprob, nlpparam) );

   /* forward solution if we solved to optimality; local optimality is enough since the NLP is convex */
   if( SCIPgetNlpiSolstat(scip, nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT )
   {
      SCIP_VAR** vars;
      SCIP_Real* primal;
      SCIP_Real relaxval;
      int nvars;
      int i;

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      SCIP_CALL( SCIPgetNlpiSolution(scip, nlpi, nlpiprob, &primal, NULL, NULL, NULL, &relaxval) );

      /* store relaxation solution in original SCIP if it improves the best relaxation solution thus far */
      if( (! SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)) )
      {
         SCIPdebugMsg(scip, "Setting NLP relaxation solution, which improved upon earlier solution\n");
         SCIP_CALL( SCIPclearRelaxSolVals(scip, relax) );

         for( i = 0; i < nvars; ++i )
         {
   #ifndef NDEBUG
            SCIP_Real lb;
            SCIP_Real ub;

            lb = SCIPvarGetLbLocal(vars[i]);
            ub = SCIPvarGetUbLocal(vars[i]);
            assert(SCIPisInfinity(scip, -lb) || SCIPisFeasLE(scip, lb, primal[i]));
            assert(SCIPisInfinity(scip, ub) || SCIPisFeasLE(scip, primal[i], ub));
            SCIPdebugMsg(scip, "relax value of %s = %g in [%g,%g]\n", SCIPvarGetName(vars[i]), primal[i], lb, ub);
   #endif

            SCIP_CALL( SCIPsetRelaxSolVal(scip, relax, vars[i], primal[i]) );
         }

         /* mark relaxation solution to be valid */
         SCIP_CALL( SCIPmarkRelaxSolValid(scip, relax, TRUE) );
      }

      SCIPdebugMsg(scip, "NLP lower bound = %g\n", relaxval);
      *lowerbound = relaxval;
      *result = SCIP_SUCCESS;
   }

   /* free memory */
   SCIPhashmapFree(&var2idx);
   SCIP_CALL( SCIPfreeNlpiProblem(scip, nlpi, &nlpiprob) );

   return SCIP_OKAY;
}


/*
 * relaxator specific interface methods
 */

/** creates the nlp relaxator and includes it in SCIP */
SCIP_RETCODE SCIPincludeRelaxNlp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_RELAXDATA* relaxdata;
   SCIP_RELAX* relax;

   /* create nlp relaxator data */
   relaxdata = NULL;
   relax = NULL;

   /* include relaxator */
   SCIP_CALL( SCIPincludeRelaxBasic(scip, &relax, RELAX_NAME, RELAX_DESC, RELAX_PRIORITY, RELAX_FREQ,
         relaxExecNlp, relaxdata) );

   assert(relax != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetRelaxInitsol(scip, relax, relaxInitsolNlp) );
   SCIP_CALL( SCIPsetRelaxExitsol(scip, relax, relaxExitsolNlp) );

   return SCIP_OKAY;
}
