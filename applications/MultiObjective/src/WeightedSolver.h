/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   WeightedSolver.h
 * @brief  Class providing all methods for using the algorithm
 * @author Timo Strunk
 *
 * Abstract superclass for a multi objective solver using weighted objective functions.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef WEIGHTED_SOLVER
#define WEIGHTED_SOLVER

#include <vector>
#include <string>
#include <map>

#include "scip/scip.h"
#include "scip/message_default.h"

/** generic weight based solver */
class WeightedSolver
{
 public:
   /** SCIP style constructor */
   WeightedSolver(
      const char*        paramfilename       /**< name of file with SCIP parameters */
      );

   /** destructor */
   virtual ~WeightedSolver();

   /** reads problem data from file */
   SCIP_RETCODE readProblem(
      const char*        filename            /**< name of instance file */
      );

   /** returns true if there is a weight left to check */
   virtual bool hasNext() const = 0;

   /** solves the next weighted problem */
   virtual SCIP_RETCODE solveNext() = 0;

   /** returns true if the last weighted run found a new pareto optimum */
   virtual bool foundNewOptimum() const;

   /** gets the last weight loaded into the solver */
   virtual const std::vector<SCIP_Real>* getWeight() const;

   /** gets cost vector of last found pareto optimum */
   virtual const std::vector<SCIP_Real>* getCost() const;

   /** true if the last solved weighted problem is unbounded */
   virtual bool isWeightedUnbounded() const = 0;

   /** returns the last found pareto optimal solution */
   virtual SCIP_SOL* getSolution() const;

   /** returns the number of branch and bound nodes in the last weighted run */
   virtual SCIP_Longint getNNodesLastRun() const;

   /** returns the number of LP iterations used in the last run */
   virtual SCIP_Longint getNLPIterationsLastRun() const;

   /** returns the time needed for the last iteration in seconds */
   SCIP_Real getDurationLastRun() const;

   /** returns the number of objective functions */
   int getNObjs() const;

   /** returns the SCIP problem status of the multiobjective problem */
   virtual SCIP_Status getStatus() const;

   /** returns the number of weighted runs so far */
   int getNRuns() const;

   /** returns the number of found pareto optima so far */
   virtual int getNSolutions() const;

   /** get total time for algorithm */
   virtual SCIP_Real getTotalDuration() const=0;

   /** get number of new vertices in the 1-skeleton added in last step*/
   virtual int getNNewVertices() const=0;

   /** get number of vertices in the 1-skeleton processed in last step*/
   virtual int getNProcessedVertices() const=0;

   /** delete non extremal solutions */
   SCIP_RETCODE checkAndWriteSolutions();

   /** return verblevel parameter set in SCIP */
   int getVerbosity() const;

   /** return a list of names for objective functions */
   const std::vector<std::string>* getObjNames() const;

   /** return a map from cost vectors to solution file names */
   const std::map< const std::vector<SCIP_Real>*, const char* >* getCostToFilename() const;

   /** return all cost vectors of unbounded primal rays */
   const std::vector< const std::vector<SCIP_Real>* >* getCostRays() const;

 protected:
   SCIP*                 scip_;                   /**< SCIP solver */
   SCIP_Real             timelimit_;              /**< remaining time in seconds */
   int                   verbosity_;
   SCIP_Status           multiopt_status_;        /**< multiobjective problem status */
   int                   nruns_;                  /**< number of solveNext() calls */

   bool                  found_new_optimum_;      /**< true if last call of solveNext() found a new optimum */
   SCIP_Longint          nnodes_last_run_;        /**< number of branch and bound nodes in last call of solveNext() */
   SCIP_Longint          niterations_last_run_;   /**< number of lp iterations in last call of solveNext() */
   SCIP_Real             duration_last_run_;      /**< duration of last call of solveNext() in seconds */
   SCIP_SOL*             solution_;               /**< last solution found by solveNext() */

   const std::vector<SCIP_Real>*                  weight_;            /**< weight used in last solveNext() call */
   const std::vector<SCIP_Real>*                  cost_vector_;       /**< cost vector of last found solution */

   std::vector< SCIP_SOL* >                       solutions_;         /**< list of all pareto optimal SCIP solutions */
   std::vector< const std::vector< SCIP_Real>* >  nondom_points_;     /**< list of found non dominated points*/
   std::vector< const std::vector<SCIP_Real>* >   cost_rays_;         /**< cost vectors of unbounded primal rays*/

   std::map< const std::vector<SCIP_Real>*,
             SCIP_SOL* >                          cost_to_sol_;       /**< maps cost vectors to solutions */

 private:
   std::string           filename_;               /**< name of problem file */
   std::string           outfilestump_;           /**< beginning of outfile names */
   SCIP_LPI*             extremality_lpi_;        /**< lp interfaced for determining extremality */
   SCIP_MESSAGEHDLR*     extremality_msg_;        /**< message handler dealing with extremality lp messages */
   bool                  candidate_is_extremal_;  /**< true if nondom point checked by lp is extremal */
   int                   n_written_sols_;         /**< number of solutions written to files */

   std::map< const std::vector<SCIP_Real>*,
             const char* >                        cost_to_filename_;  /**< maps confirmed extremal points to solution files */

   /** prepare the LP for the extremality check */
   SCIP_RETCODE createExtremalityLP();

   /** solve the extremality lp for one particular nondom point */
   SCIP_RETCODE solveExtremalityLP(const std::vector<SCIP_Real>* nondom_point, int point_index);

   /** writes a solution to a file in folder solutions/ with name <instance name>-<solution number>.sol*/
   SCIP_RETCODE writeSolution(
      const std::vector<SCIP_Real>*     cost_vector         /**< cost vector of solution */
      );

};

#endif
