/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: heur_xxx.c,v 1.1 2003/11/26 16:08:59 bzfpfend Exp $"

/**@file   heur_xxx.c
 * @brief  xxx primal heuristic
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "heur_xxx.h"


#define HEUR_NAME            "xxx"
#define HEUR_DESC            "primal heuristic template"
#define HEUR_DISPCHAR        '?'
#define HEUR_PRIORITY        0
#define HEUR_FREQ            1
#define HEUR_PSEUDONODES     TRUE       /** call heuristic at nodes where only a pseudo solution exist? */




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


/** initialization method of primal heuristic (called when problem solving starts) */
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


/** deinitialization method of primal heuristic (called when problem solving exits) */
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
   CHECK_OKAY( SCIPincludeHeur(scip, HEUR_NAME, HEUR_DESC, HEUR_DISPCHAR, HEUR_PRIORITY, HEUR_FREQ, HEUR_PSEUDONODES,
                  heurFreeXxx, heurInitXxx, heurExitXxx, heurExecXxx,
                  heurdata) );

   /* add xxx primal heuristic parameters */
   /* TODO: (optional) add primal heuristic specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
