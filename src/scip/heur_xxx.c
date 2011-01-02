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
#pragma ident "@(#) $Id: heur_xxx.c,v 1.28 2011/01/02 11:10:46 bzfheinz Exp $"

/**@file   heur_xxx.c
 * @ingroup PRIMALHEURISTICS
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
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE
#define HEUR_USESSUBSCIP      FALSE  /**< does the heuristic use a secondary SCIP instance? */




/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
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

/** copy method for primal heuristic plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_HEURCOPY(heurCopyXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurCopyXxx NULL
#endif

/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_HEURFREE(heurFreeXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurFreeXxx NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 0
static
SCIP_DECL_HEURINIT(heurInitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitXxx NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitXxx NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_HEURINITSOL(heurInitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurInitsolXxx NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define heurExitsolXxx NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the xxx primal heuristic and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurXxx(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create xxx primal heuristic data */
   heurdata = NULL;
   /* TODO: (optional) create primal heuristic specific data here */

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
         HEUR_MAXDEPTH, HEUR_TIMING, HEUR_USESSUBSCIP,
         heurCopyXxx,
         heurFreeXxx, heurInitXxx, heurExitXxx, 
         heurInitsolXxx, heurExitsolXxx, heurExecXxx,
         heurdata) );

   /* add xxx primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
