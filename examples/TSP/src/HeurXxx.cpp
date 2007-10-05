/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_xxx.c
 * @brief  xxx primal heuristic
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "HeurXxx.h"
#include "ProbDataTSP.h"

using namespace tsp;
using namespace std;


/*
 * Local methods
 */


/*
 * Callback methods of primal heuristic
 */


/** destructor of primal heuristic to free user data (called when SCIP is exiting) */
SCIP_RETCODE HeurXxx::scip_free(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** initialization method of primal heuristic (called after problem was transformed) */
SCIP_RETCODE HeurXxx::scip_init(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** deinitialization method of primal heuristic (called before transformed problem is freed) */
SCIP_RETCODE HeurXxx::scip_exit(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 */
SCIP_RETCODE HeurXxx::scip_initsol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 */
SCIP_RETCODE HeurXxx::scip_exitsol(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur                /**< the primal heuristic itself */
   )
{
   return SCIP_OKAY;
}
   
/** execution method of primal heuristic */
SCIP_RETCODE HeurXxx::scip_exec(
   SCIP*              scip,               /**< SCIP data structure */
   SCIP_HEUR*         heur,               /**< the primal heuristic itself */
   SCIP_RESULT*       result              /**< pointer to store the result of the heuristic call */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx primal heuristic not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;

}
