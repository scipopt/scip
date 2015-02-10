/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_trivialnegation.c
 * @brief  trivialnegation primal heuristic
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/heur_trivialnegation.h"
#include "scip/reopt.h"
#include "scip/struct_scip.h"

#define HEUR_NAME             "trivialnegation"
#define HEUR_DESC             "negate solution entries if a objective coefficient changes the sign, enters or leaves the objective."
#define HEUR_DISPCHAR         'j'
#define HEUR_PRIORITY         30000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         0
#define HEUR_TIMING           SCIP_HEURTIMING_BEFOREPRESOL
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
   SCIP_Bool             reopt_enabled;           /**< reoptimization enabled */
};


/*
 * Local methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTrivialnegation)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTrivialnegation(scip) );

   return SCIP_OKAY;
}

/* free data of the heuristic */
static
SCIP_DECL_HEURFREE(heurFreeTrivialnegation)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL );
   assert(heur != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL );

   SCIPfreeMemory(scip, &heurdata);
   SCIPheurSetData(heur, NULL);

   return SCIP_OKAY;
}

/* initialize the heuristic */
static
SCIP_DECL_HEURINIT(heurInitTrivialnegation)
{
   SCIP_HEURDATA* heurdata;

   assert(scip != NULL );
   assert(heur != NULL );

   heurdata = SCIPheurGetData(heur);
   assert(heurdata != NULL );

   heurdata->reopt_enabled = SCIPisReoptEnabled(scip);

   return SCIP_OKAY;
}

#define heurExitTrivialnegation NULL
#define heurInitsolTrivialnegation NULL
#define heurExitsolTrivialnegation NULL


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTrivialnegation)
{
   /*lint --e{715}*/
   SCIP_SOL* lastbestsol;         /** best solution from last run */
   SCIP_SOL* allchanged;          /** solution with all entries negated */
   SCIP_SOL* feasiblechanged;     /** solution with all feasible entries negated */
   SCIP_SOL* singlenegatedsol;    /** solution with exactly one negated entry */
   SCIP_VAR** vars;

   int nvars;
   int i;

   SCIP_Real solval;
   SCIP_Bool success;
   SCIP_Bool reopt;

   SCIP_CALL( SCIPgetBoolParam(scip, "reoptimization/enable", &reopt) );

   if( !reopt )
   {
      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   vars = SCIPgetOrigVars(scip);
   nvars = SCIPgetNOrigVars(scip);

   *result = SCIP_DIDNOTFIND;

   /* get best solution from the run */
   lastbestsol = SCIPreoptGetLastBestSol(scip->reopt);

   if( lastbestsol == NULL )
      return SCIP_OKAY;

   /* initialize data structure */
   SCIP_CALL( SCIPcreateSol(scip, &allchanged, heur) );
   SCIP_CALL( SCIPcreateSol(scip, &feasiblechanged, heur) );
   SCIP_CALL( SCIPcreateSol(scip, &singlenegatedsol, heur) );

   /* copy the solutions */
   for(i = 0; i < nvars; i++)
   {
      solval = SCIPgetSolVal(scip, lastbestsol, vars[i]);
      SCIP_CALL( SCIPsetSolVal(scip, allchanged, vars[i], solval) );
      SCIP_CALL( SCIPsetSolVal(scip, feasiblechanged, vars[i], solval) );
      SCIP_CALL( SCIPsetSolVal(scip, singlenegatedsol, vars[i], solval) );
   }

   assert(SCIPsolGetHeur(allchanged) == heur);
   assert(SCIPsolGetHeur(feasiblechanged) == heur);
   assert(SCIPsolGetHeur(singlenegatedsol) == heur);

   /* change the entries */
   for(i = 0; i < nvars; i++)
   {
      SCIP_Bool entering;
      SCIP_Bool leaving;
      SCIP_Bool negated;
      SCIP_Real obj;

      if( SCIPvarGetType(vars[i]) == SCIP_VARTYPE_BINARY
       && SCIPvarGetStatus(SCIPvarGetTransVar(vars[i])) != SCIP_VARSTATUS_FIXED
       && SCIPvarGetStatus(SCIPvarGetTransVar(vars[i])) != SCIP_VARSTATUS_AGGREGATED
       && SCIPvarGetStatus(SCIPvarGetTransVar(vars[i])) != SCIP_VARSTATUS_MULTAGGR )
      {
         /* check the changes of the variable */
         SCIPgetVarCoefChg(scip, SCIPvarGetIndex(vars[i]), &negated, &entering, &leaving);

         if( negated || entering || leaving )
         {
            solval = SCIPgetSolVal(scip, lastbestsol, vars[i]);

            /* change solution value */
            SCIP_CALL( SCIPsetSolVal(scip, allchanged, vars[i], 1 - solval) );
            SCIP_CALL( SCIPsetSolVal(scip, feasiblechanged, vars[i], 1 - solval) );
            SCIP_CALL( SCIPsetSolVal(scip, singlenegatedsol, vars[i], 1 - solval) );

            /* try solution with all changes */
            success = FALSE;
            obj = SCIPgetSolTransObj(scip, allchanged);
            if( SCIPisFeasLT(scip, obj, SCIPgetCutoffbound(scip)) )
            {
               SCIPdebugMessage("try solution with all negations\n");
               SCIP_CALL( SCIPtrySol(scip, allchanged, FALSE, FALSE, TRUE, TRUE, &success) );

               if( success )
               {
                  SCIPdebugMessage("found feasible solution solution:\n");
                  SCIPdebug( SCIP_CALL( SCIPprintSol(scip, allchanged, NULL, FALSE) ) );

                  *result = SCIP_FOUNDSOL;
               }
            }

            /* try solution with feasible changes */
            success = FALSE;
            obj = SCIPgetSolTransObj(scip, feasiblechanged);
            if( SCIPisFeasLT(scip, obj, SCIPgetCutoffbound(scip)) )
            {
               SCIPdebugMessage("try solution with feasible negations\n");
               SCIP_CALL( SCIPtrySol(scip, feasiblechanged, FALSE, FALSE, TRUE, TRUE, &success) );

               if( success )
               {
                  SCIPdebugMessage("found feasible solution solution:\n");
                  SCIPdebug( SCIP_CALL( SCIPprintSol(scip, feasiblechanged, NULL, FALSE) ) );

                  *result = SCIP_FOUNDSOL;
               }
            }

            if( !success )
            {
               /* reset solution with feasible changes */
               SCIP_CALL( SCIPsetSolVal(scip, feasiblechanged, vars[i], solval) );
            }

            /* try solution with exactly one changed value */
            obj = SCIPgetSolTransObj(scip, singlenegatedsol);
            if( SCIPisFeasLT(scip, obj, SCIPgetCutoffbound(scip)) )
            {
               success = FALSE;
               SCIPdebugMessage("try solution with a single negation\n");
               SCIP_CALL( SCIPtrySol(scip, singlenegatedsol, FALSE, FALSE, TRUE, TRUE, &success) );

               if( success )
               {
                  SCIPdebugMessage("found feasible solution:\n");
                  SCIPdebug( SCIP_CALL( SCIPprintSol(scip, singlenegatedsol, NULL, FALSE) ) );

                  *result = SCIP_FOUNDSOL;
               }
            }

            /* reset solution with exactly one changed value */
            SCIP_CALL( SCIPsetSolVal(scip, singlenegatedsol, vars[i], solval) );
         }
      }
   }

   /* free solutions */
   SCIP_CALL( SCIPfreeSol(scip, &allchanged) );
   SCIP_CALL( SCIPfreeSol(scip, &feasiblechanged) );
   SCIP_CALL( SCIPfreeSol(scip, &singlenegatedsol) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the trivialnegation primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTrivialnegation(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;
   SCIP_HEUR* heur;

   /* create trivialnegation primal heuristic data */
   SCIP_CALL( SCIPallocMemory(scip, &heurdata) );

   /* use SCIPincludeHeurBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTrivialnegation, heurdata) );

   assert(heur != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTrivialnegation) );
   SCIP_CALL( SCIPsetHeurInit(scip, heur, heurInitTrivialnegation) );
   SCIP_CALL( SCIPsetHeurFree(scip, heur, heurFreeTrivialnegation) );


   return SCIP_OKAY;
}
