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
#pragma ident "@(#) $Id: heur_fixandinfer.c,v 1.2 2004/10/12 14:06:06 bzfpfend Exp $"

/**@file   heur_fixandinfer.c
 * @brief  fix-and-infer primal heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "heur_fixandinfer.h"


#define HEUR_NAME             "fixandinfer"
#define HEUR_DESC             "iteratively fixes variables and propagates inferences"
#define HEUR_DISPCHAR         'i'
#define HEUR_PRIORITY         -500000
#define HEUR_FREQ             10
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_PSEUDONODES      TRUE      /* call heuristic at nodes where only a pseudo solution exist? */
#define HEUR_DURINGPLUNGING   FALSE     /* call heuristic during plunging? (should be FALSE for diving heuristics!) */

#define MAXDIVEDEPTH          100



/*
 * Local methods
 */

/** selects a variable and fixes it to its current pseudo solution value */
static
RETCODE fixVariable(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            pseudocands,        /**< array of unfixed variables */
   int              npseudocands        /**< number of unfixed variables */
   )
{
   VAR* var;
   Real bestscore;
   Real score;
   Real solval;
   int bestcand;
   int ncands;
   int c;

   assert(pseudocands != NULL);
   assert(npseudocands > 0);

   /* if existing, choose one of the highest priority binary variables; if no high priority binary variables
    * exist, choose a variable among all unfixed integral variables
    */
   ncands = SCIPgetNPrioPseudoBranchBins(scip);
   if( ncands == 0 )
      ncands = npseudocands;

   /* select variable to tighten the domain for */
   bestscore = -SCIPinfinity(scip);
   bestcand = -1;
   for( c = 0; c < ncands; ++c )
   {
      score = SCIPgetVarAvgInferenceScore(scip, pseudocands[c]);
      if( score > bestscore )
      {
         bestscore = score;
         bestcand = c;
      }
   }
   assert(bestcand != -1);
   
   /* fix variable to its current pseudo solution value */
   var = pseudocands[bestcand];
   solval = SCIPgetVarSol(scip, var);
   assert(SCIPisIntegral(scip, solval)); /* in probing, we always have the pseudo solution */
   debugMessage(" -> fixed variable <%s>[%g,%g] = %g\n", 
      SCIPvarGetName(var), SCIPvarGetLbLocal(var), SCIPvarGetUbLocal(var), solval);
   CHECK_OKAY( SCIPfixVarProbing(scip, var, solval) );

   return SCIP_OKAY;
}




/*
 * Callback methods of primal heuristic
 */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#define heurFreeFixandinfer NULL


/** initialization method of primal heuristic (called after problem was transformed) */
#define heurInitFixandinfer NULL


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#define heurExitFixandinfer NULL


/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecFixandinfer)
{  /*lint --e{715}*/
   VAR** cands;
   int ncands;
   int startncands;
   int divedepth;
   Bool cutoff;

   *result = SCIP_DIDNOTRUN;

   /* we cannot run on problems with continuous variables */
   if( SCIPgetNContVars(scip) > 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* start probing */
   CHECK_OKAY( SCIPstartProbing(scip) );

   /* get unfixed variables */
   CHECK_OKAY( SCIPgetPseudoBranchCands(scip, &cands, &ncands, NULL) );

   debugMessage("starting fix-and-infer heuristic with %d unfixed integral variables\n", ncands);

   /* fix variables and propagate inferences as long as the problem is still feasible and there are 
    * unfixed integral variables
    */
   cutoff = FALSE;
   divedepth = 0;
   startncands = ncands;
   while( !cutoff && ncands > 0
      && (divedepth < MAXDIVEDEPTH || (startncands - ncands) * 2 * MAXDIVEDEPTH >= startncands * divedepth) )
   {
      divedepth++;

      /* fix next variable */
      CHECK_OKAY( fixVariable(scip, cands, ncands) );

      /* propagate the fixing */
      CHECK_OKAY( SCIPpropagateProbing(scip, &cutoff) );

      /* get remaining unfixed variables */
      if( !cutoff )
      {
         CHECK_OKAY( SCIPgetPseudoBranchCands(scip, &cands, &ncands, NULL) );
      }
   }

   /* check, if we are still feasible */
   if( !cutoff )
   {
      Bool success;

      assert(ncands == 0);
      success = FALSE;

      /* try to add solution to SCIP */
      CHECK_OKAY( SCIPtryCurrentSol(scip, heur, FALSE, TRUE, &success) );

      if( success )
      {
         debugMessage("found primal feasible solution\n");
         *result = SCIP_FOUNDSOL;
      }
      else
      {
         debugMessage("primal solution was rejected\n");
      }
   }
   else
   {
      debugMessage("propagation detected a cutoff\n");
   }

   /* end probing */
   CHECK_OKAY( SCIPendProbing(scip) );

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the fix-and-infer primal heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurFixandinfer(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create fixandinfer primal heuristic data */
   heurdata = NULL;

   /* include primal heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING,
         heurFreeFixandinfer, heurInitFixandinfer, heurExitFixandinfer, heurExecFixandinfer,
         heurdata) );

   return SCIP_OKAY;
}
