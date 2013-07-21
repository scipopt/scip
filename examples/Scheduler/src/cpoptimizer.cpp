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

/**@file   cpoptimizer.cpp
 * @brief  contains method to solve a single cumulative condition via IBM ILOG CP Optimiter
 * @author Stefan Heinz
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifdef WITH_CPOPTIMIZER

#include <ilcp/cp.h>
#include "scip/cons_cumulative.h"

/** solve single cumulative condition using CP Optimizer */
SCIP_DECL_SOLVECUMULATIVE(cpoptimizer)
{
   IloEnv env;

   (*solved) = FALSE;
   (*infeasible) = FALSE;
   (*unbounded) = FALSE;
   (*error) = FALSE;

   try
   {
      int totaldemand;
      int v;

      IloModel model(env);
      IloIntervalVarArray jobs(env, njobs);

      IloCumulFunctionExpr cumulative = IloCumulFunctionExpr(env);
      IloNumExpr costs(env);

      totaldemand = 0;

      for( v = 0; v < njobs; ++v )
      {
         jobs[v] = IloIntervalVar(env);

         /* set bounds */
         jobs[v].setStartMin(ests[v]);
         jobs[v].setEndMax(lsts[v] + durations[v]);

         /* set job duration */
         jobs[v].setSizeMin(durations[v]);
         jobs[v].setSizeMax(durations[v]);

         /* add job to cumulative constraint  with corresponding demand */
         cumulative += IloPulse(jobs[v], demands[v]);

         costs += IloStartOf(jobs[v]) * objvals[v];

         totaldemand += demands[v];
      }

      if( totaldemand <= capacity )
      {
         for( v = 0; v < njobs; ++v )
         {
            lsts[v] = ests[v];
         }

         return SCIP_OKAY;
      }

      /* add objective */
      IloObjective objective(env);
      objective.setExpr(costs);
      objective.setSense(IloObjective::Minimize);

      /* add cumulative constraint to the model */

      IloIntervalVar horizon(env);
      horizon.setStartMin(hmin);
      horizon.setEndMax(hmax);
      horizon.setSizeMin(hmax - hmin);
      horizon.setSizeMax(hmax - hmin);

      cumulative += IloPulse(horizon, totaldemand - capacity);
      model.add(cumulative <= totaldemand);

      IloCP cp(model);

      cp.setParameter(IloCP::TimeLimit, timelimit);
      cp.setParameter(IloCP::ChoicePointLimit, maxnodes);
      cp.setParameter(IloCP::LogVerbosity, IloCP::Quiet);
      // cp.setParameter(IloCP::DefaultInferenceLevel, "extended");
      // cp.out() << env.getVersion() << std::endl;
      // cp.out() << cp.getVersion() << std::endl;

      // cp.out() << "Instance \t: " << filename << std::endl;
      cp.solve();

      switch( cp.getStatus() )
      {
      case IloAlgorithm::Feasible:
      case IloAlgorithm::Optimal:
         /* collect optimal solution */
         for( v = 0; v < njobs; ++v )
         {
            ests[v] = cp.getStart(jobs[v]);
            lsts[v] = cp.getStart(jobs[v]);
         }
         break;
      case IloAlgorithm::InfeasibleOrUnbounded:
         (*infeasible) = TRUE;
         (*unbounded) = TRUE;
         abort();
         break;
      case IloAlgorithm::Infeasible:
         (*infeasible) = TRUE;
         break;
      case IloAlgorithm::Unbounded:
         (*unbounded) = TRUE;
         break;
      case IloAlgorithm::Unknown:
      case IloAlgorithm::Error:
         (*error) = TRUE;
         break;
      }
   }
   catch( IloException& e )
   {
      (*error) = FALSE;
   }

   env.end();

   return SCIP_OKAY;
}

#endif
