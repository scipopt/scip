/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_xxx.c,v 1.13 2005/03/02 19:04:56 bzfpfend Exp $"

/**@file   heur_xxx.c
 * @brief  xxx primal heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/heur_xxx.h"


#define HEUR_NAME             "xxx"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         '?'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_PSEUDONODES      TRUE      /* call heuristic at nodes where only a pseudo solution exist? */
#define HEUR_DURINGPLUNGING   TRUE      /* call heuristic during plunging? (should be FALSE for diving heuristics!) */
#define HEUR_AFTERNODE        TRUE      /* call heuristic after or before the current node was solved? */




/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct HeurData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of primal heuristic
 */

/* TODO: Implement all necessary primal heuristic methods. The methods with an #if 0 ... #else #define ... are optional */

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
DECL_HEURFREE(heurFreeXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx primal heuristic not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreeXxx NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
DECL_HEURINIT(heurInitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx primal heuristic not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitXxx NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
DECL_HEUREXIT(heurExitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx primal heuristic not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitXxx NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
DECL_HEURINITSOL(heurInitsolXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx primal heuristic not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolXxx NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
DECL_HEUREXITSOL(heurExitsolXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx primal heuristic not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolXxx NULL
#endif


/** execution method of primal heuristic */
static
DECL_HEUREXEC(heurExecXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx primal heuristic not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the xxx primal heuristic and includes it in SCIP */
RETCODE SCIPincludeHeurXxx(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   HEURDATA* heurdata;

   /* create xxx primal heuristic data */
   heurdata = NULL;
   /* TODO: (optional) create primal heuristic specific data here */

   /* include primal heuristic */
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_PSEUDONODES, HEUR_DURINGPLUNGING, HEUR_AFTERNODE,
         heurFreeXxx, heurInitXxx, heurExitXxx, 
         heurInitsolXxx, heurExitsolXxx, heurExecXxx,
         heurdata) );

   /* add xxx primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
