/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   heur_trivial.c
 * @ingroup DEFPLUGINS_HEUR
 * @brief  trivial primal heuristic
 * @author Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/heur_trivial.h"
#include "scip/pub_heur.h"
#include "scip/pub_message.h"
#include "scip/pub_var.h"
#include "scip/scip_heur.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include <string.h>

#define HEUR_NAME             "trivial"
#define HEUR_DESC             "start heuristic which tries some trivial solutions"
#define HEUR_DISPCHAR         SCIP_HEURDISPCHAR_TRIVIAL
#define HEUR_PRIORITY         10000
#define HEUR_FREQ             0
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_BEFOREPRESOL | SCIP_HEURTIMING_BEFORENODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */

/*
 * Local methods
 */

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
static
SCIP_DECL_HEURCOPY(heurCopyTrivial)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(heur != NULL);
   assert(strcmp(SCIPheurGetName(heur), HEUR_NAME) == 0);

   /* call inclusion method of primal heuristic */
   SCIP_CALL( SCIPincludeHeurTrivial(scip) );

   return SCIP_OKAY;
}


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecTrivial)
{  /*lint --e{715}*/
   SCIP_VAR** vars;
   SCIP_SOL* zerosol;                   /* solution where all variables are set next to zero within bounds */
   SCIP_SOL* lbsol;                     /* solution where all variables are set to their lower bounds */
   SCIP_SOL* ubsol;                     /* solution where all variables are set to their upper bounds */
   SCIP_SOL* locksol;                   /* solution where all variables are set to the bound with the fewer locks */
   SCIP_Real large;
   SCIP_Bool difflb;
   SCIP_Bool diffub;
   SCIP_Bool difflock;
   SCIP_Bool success;
   int nvars;
   int i;

   *result = SCIP_DIDNOTFIND;

   /* initialize data structure */
   SCIP_CALL( SCIPcreateSol(scip, &zerosol, heur) );
   SCIP_CALL( SCIPcreateSol(scip, &lbsol, heur) );
   SCIP_CALL( SCIPcreateSol(scip, &ubsol, heur) );
   SCIP_CALL( SCIPcreateSol(scip, &locksol, heur) );

   /* determine large value to set variables to */
   large = SCIPround(scip, MIN(1.0 / SCIPfeastol(scip), SCIPgetHugeValue(scip)) / 10.0); /*lint !e666 */

   /* check zero solution once */
   difflb = FALSE;
   diffub = FALSE;
   difflock = FALSE;

   SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL) );
   assert(vars != NULL || nvars == 0);

   for( i = 0; i < nvars; ++i )
   {
      SCIP_Real lb;
      SCIP_Real ub;
      SCIP_Real zeroval;
      SCIP_Real solval;

      assert(vars != NULL); /* this assert is needed for flexelint */

      lb = SCIPvarGetLbLocal(vars[i]);
      ub = SCIPvarGetUbLocal(vars[i]);

      /* if problem is obviously infeasible due to empty domain, stop */
      if( SCIPisFeasGT(scip, lb, ub) )
         goto TERMINATE;

      /* set bounds to sufficient large value */
      if( SCIPisInfinity(scip, -lb) )
         lb = MIN(-large, ub);
      if( SCIPisInfinity(scip, ub) )
         ub = MAX(large, lb);

      /* set value next to zero within bounds */
      zeroval = MAX(MIN(0.0, ub), lb);

      /* set value to the bound with fewer locks, if tie choose an average value */
      if( SCIPvarGetNLocksDownType(vars[i], SCIP_LOCKTYPE_MODEL) < SCIPvarGetNLocksUpType(vars[i], SCIP_LOCKTYPE_MODEL) )
         solval = lb;
      else if( SCIPvarGetNLocksDownType(vars[i], SCIP_LOCKTYPE_MODEL) > SCIPvarGetNLocksUpType(vars[i], SCIP_LOCKTYPE_MODEL) )
         solval = ub;
      else
      {
         solval = (lb+ub)/2.0;

         /* if a tie occurs, roughly every third integer variable will be rounded up */
         if( SCIPvarGetType(vars[i]) != SCIP_VARTYPE_CONTINUOUS )
            solval = i % 3 == 0 ? SCIPceil(scip,solval) : SCIPfloor(scip,solval);

         assert(SCIPisFeasLE(scip,SCIPvarGetLbLocal(vars[i]),solval) && SCIPisFeasLE(scip,solval,SCIPvarGetUbLocal(vars[i])));
      }

      if( !SCIPisEQ(scip, lb, zeroval) )
         difflb = TRUE;

      if( !SCIPisEQ(scip, ub, zeroval) )
         diffub = TRUE;

      if( !SCIPisEQ(scip, solval, zeroval) )
         difflock = TRUE;

      /* set variable to values */
      SCIP_CALL( SCIPsetSolVal(scip, zerosol, vars[i], zeroval) );
      SCIP_CALL( SCIPsetSolVal(scip, lbsol, vars[i], lb) );
      SCIP_CALL( SCIPsetSolVal(scip, ubsol, vars[i], ub) );
      SCIP_CALL( SCIPsetSolVal(scip, locksol, vars[i], solval) );
   }

   /* try zero solution */
   SCIPdebugMsg(scip, "try zero solution\n");
   SCIP_CALL( SCIPtrySol(scip, zerosol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );

   if( success )
   {
      SCIPdebugMsg(scip, "found feasible zero solution:\n");
      SCIPdebug( SCIP_CALL( SCIPprintSol(scip, zerosol, NULL, FALSE) ) );

      *result = SCIP_FOUNDSOL;
   }

   /* try lower bound solution */
   if( difflb )
   {
      SCIPdebugMsg(scip, "try lower bound solution\n");
      SCIP_CALL( SCIPtrySol(scip, lbsol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );

      if( success )
      {
         SCIPdebugMsg(scip, "found feasible lower bound solution:\n");
         SCIPdebug( SCIP_CALL( SCIPprintSol(scip, lbsol, NULL, FALSE) ) );

         *result = SCIP_FOUNDSOL;
      }
   }

   /* try upper bound solution */
   if( diffub )
   {
      SCIPdebugMsg(scip, "try upper bound solution\n");
      SCIP_CALL( SCIPtrySol(scip, ubsol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );

      if( success )
      {
         SCIPdebugMsg(scip, "found feasible upper bound solution:\n");
         SCIPdebug( SCIP_CALL( SCIPprintSol(scip, ubsol, NULL, FALSE) ) );

         *result = SCIP_FOUNDSOL;
      }
   }

   /* try lock solution */
   if( difflock )
   {
      SCIPdebugMsg(scip, "try lock solution\n");
      SCIP_CALL( SCIPtrySol(scip, locksol, FALSE, FALSE, FALSE, FALSE, TRUE, &success) );

      if( success )
      {
         SCIPdebugMsg(scip, "found feasible lock solution:\n");
         SCIPdebug( SCIP_CALL( SCIPprintSol(scip, locksol, NULL, FALSE) ) );

         *result = SCIP_FOUNDSOL;
      }
   }

TERMINATE:
   /* free solutions */
   SCIP_CALL( SCIPfreeSol(scip, &locksol) );
   SCIP_CALL( SCIPfreeSol(scip, &ubsol) );
   SCIP_CALL( SCIPfreeSol(scip, &lbsol) );
   SCIP_CALL( SCIPfreeSol(scip, &zerosol) );

   return SCIP_OKAY;
}


/*
 * primal heuristic specific interface methods
 */

/** creates the trivial primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurTrivial(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEUR* heur;

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeurBasic(scip, &heur,
         HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP, heurExecTrivial, NULL) );

   assert(heur != NULL);

   /* set non-NULL pointers to callback methods */
   SCIP_CALL( SCIPsetHeurCopy(scip, heur, heurCopyTrivial) );

   return SCIP_OKAY;
}
