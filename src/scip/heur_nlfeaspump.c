/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2009 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_nlfeaspump.c,v 1.2 2009/07/10 08:33:10 bzfbelot Exp $"

/**@file   heur_nlfeaspump.c
 * @ingroup PRIMALHEURISTICS
 * @brief  MINLP feasibility pump primal heuristic
 * @author Pietro Belotti and Timo Berthold
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/heur_nlfeaspump.h"


#define HEUR_NAME             "MINLP Feasibility Pump"
#define HEUR_DESC             "primal heuristic template"
#define HEUR_DISPCHAR         'F'
#define HEUR_PRIORITY         0
#define HEUR_FREQ             1
#define HEUR_FREQOFS          0
#define HEUR_MAXDEPTH         -1
#define HEUR_TIMING           SCIP_HEURTIMING_AFTERNODE

#include "nlFeasPump/solution_pool.h"
#include "nlFeasPump/minlp_problem.h"
#include "nlFeasPump/minlp_cutgenerator.h"


/*
 * Data structures
 */

/* TODO: fill in the necessary primal heuristic data */

/** primal heuristic data */
struct SCIP_HeurData
{
  /*solution_pool *pool_;
  CouenneProblem *problem_;
  CouenneCutGenerator *cutgen_;  */
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
#if 1
static
SCIP_DECL_HEURFREE(heurFreeNlFeasPump)
{  /*lint --e{715}*/
  /*SCIPerrorMessage("method of the MINLP feasibility pump not implemented yet\n");
    SCIPABORT();*/ /*lint --e{527}*/

  printf ("Destroy pools, CouenneProblem, NL problem, LP problem\n");

  return SCIP_OKAY;
}
#else
#define heurFreeNlFeasPump NULL
#endif


/** initialization method of primal heuristic (called after problem was transformed) */
#if 1
static
SCIP_DECL_HEURINIT(heurInitNlFeasPump)
{  /*lint --e{715}*/
  /*SCIPerrorMessage("method of the MINLP feasibility pump not implemented yet\n");*/
  /*SCIPABORT();*/ /*lint --e{527}*/

  printf ("preprocess problem (FBBT)\n");
  printf ("generate initial convexification Ax <= b\n");

  return SCIP_OKAY;
}
#else
#define heurInitNlFeasPump NULL
#endif


/** deinitialization method of primal heuristic (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_HEUREXIT(heurExitNlFeasPump)
{  /*lint --e{715}*/
  /*SCIPerrorMessage("method of the MINLP feasibility pump not implemented yet\n");*/
  /*SCIPABORT();*/ /*lint --e{527}*/

  return SCIP_OKAY;
}
#else
#define heurExitNlFeasPump NULL
#endif


/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin) */
#if 1
static
SCIP_DECL_HEURINITSOL(heurInitsolNlFeasPump)
{  /*lint --e{715}*/
  /*SCIPerrorMessage("method of the MINLP feasibility pump not implemented yet\n");*/
  /*SCIPABORT();*/ /*lint --e{527}*/

  printf ("Generate nl_x [0]\n");

  return SCIP_OKAY;
}
#else
#define heurInitsolNlFeasPump NULL
#endif


/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_HEUREXITSOL(heurExitsolNlFeasPump)
{  /*lint --e{715}*/
  /*SCIPerrorMessage("method of the MINLP feasibility pump not implemented yet\n");
    SCIPABORT();*/ /*lint --e{527}*/

  return SCIP_OKAY;
}
#else
#define heurExitsolNlFeasPump NULL
#endif


/** execution method of primal heuristic */
static
SCIP_DECL_HEUREXEC(heurExecNlFeasPump)
{  /*lint --e{715}*/
  /*SCIPerrorMessage("method of the MINLP feasibility pump not implemented yet\n");
    SCIPABORT();*/ /*lint --e{527}*/

  printf ("Call heuristic\n");

  /*NlFeasPumpWrapper (data);*/

  return SCIP_OKAY;
}





/*
 * primal heuristic specific interface methods
 */

/** creates the MINLP feasibility pump and includes it in SCIP */
SCIP_RETCODE SCIPincludeHeurNlFeasPump(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HEURDATA* heurdata;

   /* create MINLP feasibility pump data */
   heurdata = NULL;
   /* TODO: (optional) create primal heuristic specific data here */

   SCIP_HEURDATA *data = (SCIP_HEURDATA *) malloc (sizeof (SCIP_HEURDATA));
   /*fillNlFeasPumpData (data);*/

   /* include primal heuristic */
   SCIP_CALL( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, 
			      HEUR_PRIORITY, HEUR_FREQ, HEUR_FREQOFS,
			      HEUR_MAXDEPTH, HEUR_TIMING,
			      heurFreeNlFeasPump, heurInitNlFeasPump, heurExitNlFeasPump, 
			      heurInitsolNlFeasPump, heurExitsolNlFeasPump, heurExecNlFeasPump,
			      heurdata) );

   /* add MINLP feasibility pump parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
