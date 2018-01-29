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

/**@file   relax_nlp.c
 * @brief  nlp relaxator
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "nlpi/nlpi.h"
#include "relax_nlp.h"


#define RELAX_NAME             "nlp"
#define RELAX_DESC             "relaxator solving a convex NLP relaxation"
#define RELAX_PRIORITY         10
#define RELAX_FREQ             1

#define NLPITERLIMIT           500       /**< iteration limit of NLP solver */
#define NLPVERLEVEL            0         /**< verbosity level of NLP solver */
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
   SCIP_Real timelimit;
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

   SCIP_CALL( SCIPnlpiCreateProblem(nlpi, &nlpiprob, "relax-NLP") );
   SCIP_CALL( SCIPhashmapCreate(&var2idx, SCIPblkmem(scip), SCIPgetNVars(scip)) );

   SCIP_CALL( SCIPcreateNlpiProb(scip, nlpi, nlrows, nnlrows, nlpiprob, var2idx, NULL, SCIPgetCutoffbound(scip),
         TRUE, TRUE) );
   SCIP_CALL( SCIPaddNlpiProbRows(scip, nlpi, nlpiprob, var2idx, SCIPgetLPRows(scip), SCIPgetNLPRows(scip)) );

   /* set working limits */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit <= 1.0 )
      {
         SCIPdebugMsg(scip, "skip NLP solve; no time left\n");
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_TILIM, timelimit) );
   SCIP_CALL( SCIPnlpiSetIntPar(nlpi, nlpiprob, SCIP_NLPPAR_ITLIM, NLPITERLIMIT) );
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_FEASTOL, SCIPfeastol(scip) * FEASTOLFAC) );
   SCIP_CALL( SCIPnlpiSetRealPar(nlpi, nlpiprob, SCIP_NLPPAR_RELOBJTOL, SCIPfeastol(scip) * RELOBJTOLFAC) );
   SCIP_CALL( SCIPnlpiSetIntPar(nlpi, nlpiprob, SCIP_NLPPAR_VERBLEVEL, NLPVERLEVEL) );

   /* solve NLP */
   SCIP_CALL( SCIPnlpiSolve(nlpi, nlpiprob) );

   /* forward solution if we solved to optimality; local optimality is enough since the NLP is convex */
   if( SCIPnlpiGetSolstat(nlpi, nlpiprob) <= SCIP_NLPSOLSTAT_LOCOPT )
   {
      SCIP_VAR** vars;
      SCIP_Real* primal;
      SCIP_Real relaxval;
      int nvars;
      int i;

      vars = SCIPgetVars(scip);
      nvars = SCIPgetNVars(scip);

      SCIP_CALL( SCIPnlpiGetSolution(nlpi, nlpiprob, &primal, NULL, NULL, NULL, &relaxval) );

      /* store relaxation solution in original SCIP if it improves the best relaxation solution thus far */
      if( (! SCIPisRelaxSolValid(scip)) || SCIPisGT(scip, relaxval, SCIPgetRelaxSolObj(scip)) )
      {
         SCIPdebugMsg(scip, "Setting NLP relaxation solution, which improved upon earlier solution\n");
         SCIP_CALL( SCIPclearRelaxSolVals(scip) );

         for( i = 0; i < nvars; ++i )
         {
   #ifndef NDEBUG
            SCIP_Real lb;
            SCIP_Real ub;

            lb = SCIPvarGetLbLocal(vars[i]);
            ub = SCIPvarGetUbLocal(vars[i]);
            assert(SCIPisInfinity(scip, -lb) || SCIPisLE(scip, lb, primal[i]));
            assert(SCIPisInfinity(scip, ub) || SCIPisLE(scip, primal[i], ub));
            SCIPdebugMsg(scip, "relax value of %s = %g in [%g,%g]\n", SCIPvarGetName(vars[i]), primal[i], lb, ub);
   #endif

            SCIP_CALL( SCIPsetRelaxSolVal(scip, vars[i], primal[i]) );
         }

         /* mark relaxation solution to be valid */
         SCIP_CALL( SCIPmarkRelaxSolValid(scip, TRUE) );
      }

      SCIPdebugMsg(scip, "NLP lower bound = %g\n", relaxval);
      *lowerbound = relaxval;
      *result = SCIP_SUCCESS;
   }

   /* free memory */
   SCIPhashmapFree(&var2idx);
   SCIP_CALL( SCIPnlpiFreeProblem(nlpi, &nlpiprob) );

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
