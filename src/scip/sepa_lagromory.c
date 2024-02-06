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

/**@file   sepa_lagromory.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  Lagromory separator
 * @author Suresh Bolusani
 * @author Mark Ruben Turner
 * @author Mathieu Besançon
 *
 * This separator is based on the following article that discusses Lagromory separation using the relax-and-cut
 * framework. Multiple enhancements have been implemented on top of the basic algorithm described in the article.
 *
 * Fischetti M. and Salvagnin D. (2011).@n
 * A relax-and-cut framework for Gomory mixed-integer cuts.@n
 * Mathematical Programming Computation, 3, 79-102.
 *
 * Consider the following linear relaxation at a node:
 *
 * \f[
 *    \begin{array}{rrl}
 *       \min & c^T x &\\
 *            & x & \in P,
 *    \end{array}
 * \f]
 *
 * where \f$P\f$ is the feasible region of the relaxation. Let the following be the cuts generated so far in the current
 * separation round.
 *
 * \f[
 *    {\alpha^i}^T x \leq \alpha^i_0, i = 1, 2, \hdots, M
 * \f]
 *
 * Then, the following is the Lagrangian dual problem considered in the relax-and-cut framework used in the separator.
 *
 * \f[
 *    z_D := \max\limits_{u \geq 0} \left\{L(u) := \min \left\{c^T x + \sum\limits_{i = 1}^{M} \left(u_i
 *    \left({\alpha^i}^T x - \alpha^i_0\right) \right) \mid x \in P\right\} \right\},
 * \f]
 *
 * where \f$u\f$ are the Lagrangian multipliers (referred to as \a dualvector in this separator) used for penalizing the
 * violation of the generated cuts, and \f$z_D\f$ is the optimal objective value (which is approximated via \a ubparam in this separator).
 * Then, the following are the steps of the relax-and-cut algorithm implemented in this separator.
 *
 * \begin{itemize}
 *    \item Generate an initial pool of cuts to build the initial Lagrangian dual problem.
 *    \item Select initial values for Lagrangian multipliers \f$u^0\f$ (e.g., all zeroes vector).
 *    \item In the outer main loop \f$i\f$ of the algorithm:
 *       \begin{enumerate}
 *          \item Solve the Lagrangian dual problem until certain termination criterion is met. This results in an inner
 *          subgradient loop, whose iteration \f$j\f$ is described below.
 *             \begin{enumerate}
 *                \item Fix \f$u^j\f$, and solve the LP corresponding to the Lagrangian dual with fixed multipliers.
 *                Gather its optimal simplex tableau and optimal objective value (i.e., the Lagrangian value)
 *                \f$L(u^j)\f$.
 *                \item Update \f$u^j\f$ to \f$u^{j+1}\f$ as follows.
 *                   \f[
 *                      u^{j+1} = \left(u^j + \lambda_j s^k\right)_+,
 *                   \f]
 *                   where \f$\lambda_j\f$ is the step length:
 *                   \f[
 *                      \lambda_j = \frac{\mu_j (UB - L(u^j))}{\|s^j\|^2_2},
 *                   \f]
 *                   where \f$mu_j\f$ is a factor (i.e., \a muparam) such that \f$0 < \mu_j \leq 2\f$, UB is \p ubparam,
 *                   and \f$s^j\f$ is the subgradient vector defined as:
 *                   \f[
 *                      s^j_k = \left({\alpha^k}^T x - \alpha^k_0\right), k = 1, 2, \hdots, M.
 *                   \f]
 *                   The factor \f$mu_j\f$ is updated as below.
 *                   \f[
 *                      mu_j = \begin{cases}
 *                               \p mubacktrackfactor * mu_j & \text{if } L(u^j) < bestLB - \delta\\
 *                               \begin{cases}
 *                                  \p muslab1factor * mu_j & \text{if } bestLB - avgLB < \p deltaslab1ub * delta\\
 *                                  \p muslab2factor * mu_j & \text{if } \p deltaslab1ub * \delta \leq bestLB - avgLB < \p deltaslab2ub * delta\\
 *                                  \p muslab3factor * mu_j & \text{otherwise}
 *                               \end{cases} & \text{otherwise},
 *                             \end{cases}
 *                   \f]
 *                   where \f$bestLB\f$ and \f$avgLB\f$ are best and average Lagrangian values found so far, and
 *                   \f$\delta = UB - bestLB\f$.
 *                \item Stabilize \f$u^{j+1}\f$ by projecting onto a norm ball followed by taking a convex combination
 *                with a core vector of Lagrangian multipliers.
 *                \item Generate GMI cuts based on the optimal simplex tableau.
 *                \item Relax the newly generated cuts by penalizing and adding them to the objective function.
 *                \item Go to the next iteration \f$j+1\f$.
 *             \end{enumerate}
 *          \item Gather all the generated cuts and build an LP by adding all these cuts to the node relaxation.
 *          \item Solve this LP to obtain its optimal primal and dual solutions.
 *          \item If this primal solution is MIP primal feasible, then add this solution to the solution pool, add all
 *          the generated cuts to the cutpool or sepastore as needed, and exit the separator.
 *          \item Otherwise, update the Lagrangian multipliers based on this optimal dual solution, and go to the next
 *          iteration \f$i+1\f$.
 *       \end{enumerate}
 * \end{itemize}
 *
 * @todo store all LP sols in a data structure, and send them to fix-and-propagate at the end.
 *
 * @todo test heuristics such as feasibility pump with multiple input solutions.
 *
 * @todo find dual degenerate problems and test the separator on these problems.
 *
 * @todo identify instance classes where these cuts work better.
 *
 * @todo add termination criteria based on failed efforts.
 *
 * @todo for warm starting, if there are additional rows/cols, set their basis status to non-basic and then set WS info.
 *
 * @todo filter cuts using multiple explored LP solutions.
 *
 * @todo for bases on optimal face only, aggregate to get a new basis and separate it.
 *
 * @todo generate other separators in addition to GMI cuts (0-1/2)
 *
 * @todo: convert iters from int to SCIP_Longint
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/heur_trysol.h"
#include "scip/sepa_lagromory.h"

#define SEPA_NAME                              "lagromory"
#define SEPA_DESC                              "separator for Lagromory cuts for MIP relaxations"
#define SEPA_PRIORITY                                -8000
#define SEPA_FREQ                                       -1
#define SEPA_MAXBOUNDDIST                              1.0
#define SEPA_USESSUBSCIP                             FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                                   FALSE /**< should separation method be delayed, if other separators found cuts? */

/* generic parameters */
#define DEFAULT_AWAY                                  0.01 /**< minimal integrality violation of a basis variable to try separation */
#define DEFAULT_DELAYEDCUTS                          FALSE /**< should cuts be added to the delayed cut pool? */
#define DEFAULT_SEPARATEROWS                          TRUE /**< separate rows with integral slack? */
#define DEFAULT_SORTCUTOFFSOL                         TRUE /**< sort fractional integer columns based on fractionality? */
#define DEFAULT_SIDETYPEBASIS                         TRUE /**< choose side types of row (lhs/rhs) based on basis information? */
#define DEFAULT_DYNAMICCUTS                           TRUE /**< should generated cuts be removed from the LP if they are no longer tight? */
#define DEFAULT_MAKEINTEGRAL                         FALSE /**< try to scale all cuts to integral coefficients? */
#define DEFAULT_FORCECUTS                            FALSE /**< force cuts to be added to the LP? */
#define DEFAULT_ALLOWLOCAL                           FALSE /**< should locally valid cuts be generated? */

/* parameters related to the separator's termination check */
#define DEFAULT_MAXROUNDSROOT                            1 /**< maximal number of separation rounds in the root node (-1: unlimited) */
#define DEFAULT_MAXROUNDS                                1 /**< maximal number of separation rounds per node (-1: unlimited) */
#define DEFAULT_DUALDEGENERACYRATETHRESHOLD            0.5 /**< minimum dual degeneracy rate for separator execution */
#define DEFAULT_VARCONSRATIOTHRESHOLD                  1.0 /**< minimum variable-constraint ratio on optimal face for separator execution */
#define DEFAULT_MINRESTART                               1 /**< minimum restart round for separator execution (0: from beginning of the
                                                             instance solving, >= n with n >= 1: from restart round n) */
#define DEFAULT_PERLPMAXCUTSROOT                        50 /**< maximal number of cuts separated per Lagromory LP in the root node */
#define DEFAULT_PERLPMAXCUTS                            10 /**< maximal number of cuts separated per Lagromory LP in the non-root node */
#define DEFAULT_PERROUNDLPITERLIMITFACTOR             -1.0 /**< factor w.r.t. root node LP iterations for maximal separating LP iterations
                                                             per separation round (negative for no limit) */
#define DEFAULT_ROOTLPITERLIMITFACTOR                 -1.0 /**< factor w.r.t. root node LP iterations for maximal separating LP iterations
                                                             in the root node (negative for no limit) */
#define DEFAULT_TOTALLPITERLIMITFACTOR                -1.0 /**< factor w.r.t. root node LP iterations for maximal separating LP iterations
                                                             in the tree (negative for no limit) */
#define DEFAULT_PERROUNDMAXLPITERS                   50000 /**< maximal number of separating LP iterations per separation round (-1: unlimited) */
#define DEFAULT_PERROUNDCUTSFACTORROOT                 1.0 /**< factor w.r.t. number of integer columns for number of cuts separated per
                                                             separation round in root node */
#define DEFAULT_PERROUNDCUTSFACTOR                     0.5 /**< factor w.r.t. number of integer columns for number of cuts separated per
                                                             separation round at a non-root node */
#define DEFAULT_TOTALCUTSFACTOR                       50.0 /**< factor w.r.t. number of integer columns for total number of cuts separated */
#define DEFAULT_MAXMAINITERS                             4 /**< maximal number of main loop iterations of the relax-and-cut algorithm */
#define DEFAULT_MAXSUBGRADIENTITERS                      6 /**< maximal number of subgradient loop iterations of the relax-and-cut algorithm */

/* parameters related to the relax-and-cut algorithm */
#define DEFAULT_MUPARAMCONST                          TRUE /**< is the mu parameter (factor for step length) constant? */
#define DEFAULT_MUPARAMINIT                           0.01 /**< initial value of the mu parameter (factor for step length) */
#define DEFAULT_MUPARAMLB                              0.0 /**< lower bound for the mu parameter (factor for step length) */
#define DEFAULT_MUPARAMUB                              2.0 /**< upper bound for the mu parameter (factor for step length) */
#define DEFAULT_MUBACKTRACKFACTOR                      0.5 /**< factor of mu while backtracking the mu parameter (factor for step length) -
                                                             see updateMuSteplengthParam() */
#define DEFAULT_MUSLAB1FACTOR                         10.0 /**< factor of mu parameter (factor for step length) for larger increment - see
                                                             updateMuSteplengthParam() */
#define DEFAULT_MUSLAB2FACTOR                          2.0 /**< factor of mu parameter (factor for step length) for smaller increment - see
                                                             updateMuSteplengthParam() */
#define DEFAULT_MUSLAB3FACTOR                          0.5 /**< factor of mu parameter (factor for step length) for reduction - see
                                                             updateMuSteplengthParam() */
#define DEFAULT_DELTASLAB1UB                         0.001 /**< factor of delta deciding larger increment of mu parameter (factor for step
                                                             length) - see updateMuSteplengthParam() */
#define DEFAULT_DELTASLAB2UB                          0.01 /**< factor of delta deciding smaller increment of mu parameter (factor for step
                                                             length) - see updateMuSteplengthParam() */
#define DEFAULT_UBPARAMPOSFACTOR                       2.0 /**< factor for positive upper bound used as an estimate for the optimal
                                                             Lagrangian dual value */
#define DEFAULT_UBPARAMNEGFACTOR                       0.5 /**< factor for negative upper bound used as an estimate for the optimal
                                                             Lagrangian dual value */
#define DEFAULT_MAXLAGRANGIANVALSFORAVG                  2 /**< maximal number of iterations for rolling average of Lagrangian value */
#define DEFAULT_MAXCONSECITERSFORMUUPDATE               10 /**< consecutive number of iterations used to determine if mu needs to be backtracked */
#define DEFAULT_PERROOTLPITERFACTOR                    0.2 /**< factor w.r.t. root node LP iterations for iteration limit of each separating
                                                             LP (negative for no limit) */
#define DEFAULT_PERLPITERFACTOR                        0.1 /**< factor w.r.t. non-root node LP iterations for iteration limit of each
                                                             separating LP (negative for no limit) */
#define DEFAULT_CUTGENFREQ                               1 /**< frequency of subgradient iterations for generating cuts */
#define DEFAULT_CUTADDFREQ                               1 /**< frequency of subgradient iterations for adding cuts to objective function */
#define DEFAULT_CUTSFILTERFACTOR                       1.0 /**< fraction of generated cuts per explored basis to accept from separator */
#define DEFAULT_OPTIMALFACEPRIORITY                      2 /**< priority of the optimal face for separator execution (0: low priority, 1:
                                                             medium priority, 2: high priority) */
#define DEFAULT_AGGREGATECUTS                         TRUE /**< aggregate all generated cuts using the Lagrangian multipliers? */
/* parameters for stabilization of the Lagrangian multipliers */
#define DEFAULT_PROJECTIONTYPE                           2 /**< the ball into which the Lagrangian multipliers are projected for
                                                             stabilization (0: no projection, 1: L1-norm ball projection, 2: L2-norm ball
                                                             projection, 3: L_inf-norm ball projection) */
#define DEFAULT_STABILITYCENTERTYPE                      1 /**< type of stability center for taking weighted average of Lagrangian multipliers for
                                                             stabilization (0: no weighted stabilization, 1: best Lagrangian multipliers) */
#define DEFAULT_RADIUSINIT                             0.5 /**< initial radius of the ball used in stabilization of Lagrangian multipliers */
#define DEFAULT_RADIUSMAX                             20.0 /**< maximum radius of the ball used in stabilization of Lagrangian multipliers */
#define DEFAULT_RADIUSMIN                             1e-6 /**< minimum radius of the ball used in stabilization of Lagrangian multipliers */
#define DEFAULT_CONST                                  2.0 /**< a constant for stablity center based stabilization of Lagrangian multipliers */
#define DEFAULT_RADIUSUPDATEWEIGHT                    0.98 /**< multiplier to evaluate cut violation score used for updating ball radius */

/* macros that are used directly */
#define RANDSEED                                        42 /**< random seed */
#define MAKECONTINTEGRAL                             FALSE /**< convert continuous variable to integral variables in SCIPmakeRowIntegral()? */
#define POSTPROCESS                                   TRUE /**< apply postprocessing after MIR calculation? - see SCIPcalcMIR() */
#define BOUNDSWITCH                                 0.9999 /**< threshold for bound switching - see SCIPcalcMIR() */
#define USEVBDS                                       TRUE /**< use variable bounds? - see SCIPcalcMIR() */
#define FIXINTEGRALRHS                               FALSE /**< try to generate an integral rhs? - see SCIPcalcMIR() */
#define MAXAGGRLEN(ncols)               (0.1*(ncols)+1000) /**< maximal length of base inequality */

/*
 * Data structures
 */

/** separator data */
struct SCIP_SepaData
{
   /* generic variables */
   SCIP_Real             away;                        /**< minimal integrality violation of a basis variable to try separation */
   SCIP_Bool             delayedcuts;                 /**< should cuts be added to the delayed cut pool? */
   SCIP_Bool             separaterows;                /**< separate rows with integral slack? */
   SCIP_Bool             sortcutoffsol;               /**< sort fractional integer columns based on fractionality? */
   SCIP_Bool             sidetypebasis;               /**< choose side types of row (lhs/rhs) based on basis information? */
   SCIP_Bool             dynamiccuts;                 /**< should generated cuts be removed from the LP if they are no longer tight? */
   SCIP_Bool             makeintegral;                /**< try to scale all cuts to integral coefficients? */
   SCIP_Bool             forcecuts;                   /**< force cuts to be added to the LP? */
   SCIP_Bool             allowlocal;                  /**< should locally valid cuts be generated? */
   SCIP_RANDNUMGEN*      randnumgen;                  /**< random number generator */
   SCIP_HEUR*            heurtrysol;                  /**< a pointer to the trysol heuristic, if available */

   /* used to define separating LPs */
   SCIP_LPI*             lpiwithsoftcuts;             /**< pointer to the lpi interface of Lagrangian dual with fixed multipliers */
   SCIP_LPI*             lpiwithhardcuts;             /**< pointer to the lpi interface of LP with all generated cuts */
   int                   nrowsinhardcutslp;           /**< nrows of \a lpiwithhardcuts */
   int                   nrunsinsoftcutslp;           /**< number of branch-and-bound runs on current instance */

   /* used for termination checks */
   SCIP_Longint          ncalls;                      /**< number of calls to the separator */
   int                   maxroundsroot;               /**< maximal number of separation rounds in the root node (-1: unlimited) */
   int                   maxrounds;                   /**< maximal number of separation rounds per node (-1: unlimited) */
   SCIP_Real             dualdegeneracyratethreshold; /**< minimum dual degeneracy rate for separator execution */
   SCIP_Real             varconsratiothreshold;       /**< minimum variable-constraint ratio on optimal face for separator execution */
   int                   minrestart;                  /**< minimum restart round for separator execution (0: from beginning of the instance
                                                        solving, >= n with n >= 1: from restart round n) */
   int                   nmaxcutsperlproot;           /**< maximal number of cuts separated per Lagromory LP in the root node */
   int                   nmaxcutsperlp;               /**< maximal number of cuts separated per Lagromory LP in the non-root node */
   SCIP_Real             perroundlpiterlimitfactor;   /**< factor w.r.t. root node LP iterations for maximal separating LP iterations per
                                                        separation round (negative for no limit) */
   SCIP_Real             rootlpiterlimitfactor;       /**< factor w.r.t. root node LP iterations for maximal separating LP iterations in the
                                                        root node (negative for no limit) */
   SCIP_Real             totallpiterlimitfactor;      /**< factor w.r.t. root node LP iterations for maximal separating LP iterations in the
                                                        tree (negative for no limit) */
   int                   perroundnmaxlpiters;         /**< maximal number of separating LP iterations per separation round (-1: unlimited) */
   SCIP_Real             perroundcutsfactorroot;      /**< factor w.r.t. number of integer columns for number of cuts separated per
                                                        separation round in root node */
   SCIP_Real             perroundcutsfactor;          /**< factor w.r.t. number of integer columns for number of cuts separated per
                                                        separation round at a non-root node */
   SCIP_Real             totalcutsfactor;             /**< factor w.r.t. number of integer columns for total number of cuts separated */
   int                   nmaxmainiters;               /**< maximal number of main loop iterations of the relax-and-cut algorithm */
   int                   nmaxsubgradientiters;        /**< maximal number of subgradient loop iterations of the relax-and-cut algorithm */
   int                   nmaxperroundlpiters;         /**< maximal number of separating LP iterations per separation round */
   int                   nmaxrootlpiters;             /**< maximal number of separating LP iterations in the root node */
   int                   nrootlpiters;                /**< number of separating LP iterations in the root node */
   int                   nmaxtotallpiters;            /**< maximal number of separating LP iterations in the tree */
   int                   ntotallpiters;               /**< number of separating LP iterations in the tree */
   int                   nmaxperroundcutsroot;        /**< maximal number of cuts separated per separation round in root node */
   int                   nmaxperroundcuts;            /**< maximal number of cuts separated per separation round */
   int                   nmaxtotalcuts;               /**< maximal number of cuts separated in the tree */
   int                   ntotalcuts;                  /**< number of cuts separated in the tree */

   /* used for the relax-and-cut algorithm */
   SCIP_Bool             muparamconst;                /**< is the mu parameter (factor for step length) constant? */
   SCIP_Real             muparaminit;                 /**< initial value of the mu parameter (factor for step length) */
   SCIP_Real             muparamlb;                   /**< lower bound for the mu parameter (factor for step length) */
   SCIP_Real             muparamub;                   /**< upper bound for the mu parameter (factor for step length) */
   SCIP_Real             mubacktrackfactor;           /**< factor of mu while backtracking the mu parameter (factor for step length) - see
                                                        updateMuSteplengthParam() */
   SCIP_Real             muslab1factor;               /**< factor of mu parameter (factor for step length) for larger increment - see
                                                        updateMuSteplengthParam() */
   SCIP_Real             muslab2factor;               /**< factor of mu parameter (factor for step length) for smaller increment - see
                                                        updateMuSteplengthParam() */
   SCIP_Real             muslab3factor;               /**< factor of mu parameter (factor for step length) for reduction - see updateMuSteplengthParam() */
   SCIP_Real             deltaslab1ub;                /**< factor of delta deciding larger increment of mu parameter (factor for step
                                                        length) - see updateMuSteplengthParam() */
   SCIP_Real             deltaslab2ub;                /**< factor of delta deciding smaller increment of mu parameter (factor for step
                                                        length) - see updateMuSteplengthParam() */
   SCIP_Real             ubparamposfactor;            /**< factor for positive upper bound used as an estimate for the optimal Lagrangian
                                                        dual value */
   SCIP_Real             ubparamnegfactor;            /**< factor for negative upper bound used as an estimate for the optimal Lagrangian
                                                        dual value */
   int                   nmaxlagrangianvalsforavg;    /**< maximal number of iterations for rolling average of Lagrangian value */
   int                   nmaxconsecitersformuupdate;  /**< consecutive number of iterations used to determine if mu needs to be backtracked */
   SCIP_Real             perrootlpiterfactor;         /**< factor w.r.t. root node LP iterations for iteration limit of each separating LP
                                                        (negative for no limit) */
   SCIP_Real             perlpiterfactor;             /**< factor w.r.t. non-root node LP iterations for iteration limit of each separating
                                                        LP (negative for no limit) */
   int                   cutgenfreq;                  /**< frequency of subgradient iterations for generating cuts */
   int                   cutaddfreq;                  /**< frequency of subgradient iterations for adding cuts to objective function */
   SCIP_Real             cutsfilterfactor;            /**< fraction of generated cuts per explored basis to accept from separator */
   int                   optimalfacepriority;         /**< priority of the optimal face for separator execution (0: low priority, 1: medium
                                                        priority, 2: high priority) */
   SCIP_Bool             aggregatecuts;               /**< aggregate all generated cuts using the Lagrangian multipliers? */

   /* for stabilization of Lagrangian multipliers */
   int                   projectiontype;              /**< the ball into which the Lagrangian multipliers are projected for stabilization
                                                        (0: no projection, 1: L1-norm ball projection, 2: L2-norm ball projection, 3:
                                                        L_inf-norm ball projection) */
   int                   stabilitycentertype;         /**< type of stability center for taking weighted average of Lagrangian multipliers for
                                                        stabilization (0: no weighted stabilization, 1: best Lagrangian multipliers) */
   SCIP_Real             radiusinit;                  /**< initial radius of the ball used in stabilization of Lagrangian multipliers */
   SCIP_Real             radiusmax;                   /**< maximum radius of the ball used in stabilization of Lagrangian multipliers */
   SCIP_Real             radiusmin;                   /**< minimum radius of the ball used in stabilization of Lagrangian multipliers */
   SCIP_Real             constant;                    /**< a constant for stablity center based stabilization of Lagrangian multipliers */
   SCIP_Real             radiusupdateweight;          /**< multiplier to evaluate cut violation score used for updating ball radius */
};


/*
 * Local methods
 */

/** start the diving mode for solving LPs corresponding to the Lagrangian dual with fixed multipliers in the subgradient
 * loop of the separator, and update some sepadata values */
static
SCIP_RETCODE createLPWithSoftCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data structure */
   )
{
   SCIP_ROW** rows;
   SCIP_COL** cols;
   int nruns;
   int runnum;
   int nrows;
   int ncols;
   unsigned int nintcols;
   SCIP_Longint nrootlpiters;

   assert(scip != NULL);
   assert(sepadata != NULL);

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   runnum = SCIPgetNRuns(scip);  /* current run number of SCIP (starts from 1) indicating restarts */
   nruns = sepadata->nrunsinsoftcutslp; /* previous run number of SCIP in which the diving LP was created */
   nintcols = 0;

   /* start diving mode so that the diving LP can be used for all subgradient iterations */
   SCIP_CALL( SCIPstartDive(scip) );

   /* store LPI pointer to be able to use LPI functions directly (e.g., setting time limit) */
   SCIP_CALL( SCIPgetLPI(scip, &(sepadata->lpiwithsoftcuts)) );

   /* if called for the first time in a restart (including the initial run), set certain sepadata values */
   if( nruns != runnum )
   {
      /* get number of LP iterations of root node's first LP solving */
      nrootlpiters = SCIPgetNRootFirstLPIterations(scip);

      /* calculate maximum number of LP iterations allowed for all separation calls in the root node */
      if( (sepadata->rootlpiterlimitfactor >= 0.0) && !SCIPisInfinity(scip, sepadata->rootlpiterlimitfactor) )
      {
         sepadata->nmaxrootlpiters = (int)(sepadata->rootlpiterlimitfactor * nrootlpiters);
      }
      else
      {
         sepadata->nmaxrootlpiters = -1; /* no finite limit */
      }

      /* calculate maximum number of LP iterations allowed for all separation calls in the entire tree */
      if( (sepadata->totallpiterlimitfactor >= 0.0) && !SCIPisInfinity(scip, sepadata->totallpiterlimitfactor) )
      {
         sepadata->nmaxtotallpiters = (int)(sepadata->totallpiterlimitfactor * nrootlpiters);
      }
      else
      {
         sepadata->nmaxtotallpiters = -1; /* no finite limit */
      }

      /* calculate maximum number of LP iterations allowed per separation call */
      if( (sepadata->perroundlpiterlimitfactor >= 0.0) && !SCIPisInfinity(scip, sepadata->perroundlpiterlimitfactor) )
      {
         sepadata->nmaxperroundlpiters = (int)(sepadata->perroundlpiterlimitfactor * nrootlpiters);
      }
      else
      {
         sepadata->nmaxperroundlpiters = -1; /* no finite limit */
      }

      /* update maximum number of LP iterations allowed per separation call using absolute limits */
      if( sepadata->perroundnmaxlpiters > 0 )
      {
         sepadata->nmaxperroundlpiters = ((sepadata->nmaxperroundlpiters >= 0) ? MIN(sepadata->nmaxperroundlpiters,
                  sepadata->perroundnmaxlpiters) : sepadata->perroundnmaxlpiters);
      }

      /* set maximum number of cuts allowed to generate per round in root and non-root nodes as well as the total tree */
      for( int i = 0; i < ncols; ++i )
      {
         nintcols += SCIPcolIsIntegral(cols[i]);
      }
      sepadata->nmaxperroundcutsroot = (int)(sepadata->perroundcutsfactorroot * nintcols);
      sepadata->nmaxperroundcuts = (int)(sepadata->perroundcutsfactor * nintcols);

      if( sepadata->ncalls == 0 )
      {
         sepadata->nmaxtotalcuts = (int)(sepadata->totalcutsfactor * nintcols);
         sepadata->ntotalcuts = 0;
      }

      /* update the run number of solving to represent the restart number in which the above limits were set */
      sepadata->nrunsinsoftcutslp = runnum;
   }

   return SCIP_OKAY;
}

/** end the diving mode that was used for solving LPs corresponding to the Lagrangian dual with fixed multipliers */
static
SCIP_RETCODE deleteLPWithSoftCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata            /**< separator data structure */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);

   SCIP_CALL( SCIPendDive(scip) );

   return SCIP_OKAY;
}

/** set up LP interface to solve LPs in the (outer) main loop of the relax-and-cut algorithm; these LPs are built by
 * adding all the generated cuts to the node relaxation */
/* @todo add lpi iters to global statistics */
static
SCIP_RETCODE createLPWithHardCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_ROW**            cuts,               /**< generated cuts to be added to the LP */
   int                   ncuts               /**< number of generated cuts to be added to the LP */
   )
{
   SCIP_ROW** rows;
   SCIP_COL** cols;
   SCIP_Real* collb;
   SCIP_Real* colub;
   SCIP_Real* colobj;
   SCIP_Real* rowlhs;
   SCIP_Real* rowrhs;
   SCIP_Real rowconst;
   SCIP_Real* rowvals;
   SCIP_Real* rowval;
   SCIP_COL** rowcols;
   SCIP_LPI* wslpi;
   SCIP_LPISTATE* lpistate;
   BMS_BLKMEM* blkmem;
   SCIP_Real pinf;
   SCIP_Real ninf;
   int nrows;
   int ncols;
   int nrownonz;
   int collppos;
   int* rowcolinds;
   int* rowbegs;

   assert(scip != NULL);
   assert(sepadata != NULL);

   SCIP_LPI** lpi = &(sepadata->lpiwithhardcuts);
   blkmem = SCIPblkmem(scip);

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );

   /* if this function is called for the first time in this separation round, then create an LPI and add cols & rows */
   if( ncuts == 0 )
   {
      if( *lpi != NULL )
      {
         SCIP_CALL( SCIPlpiFree(lpi) );
         *lpi = NULL;
      }
      assert(*lpi == NULL);

      /* create an LPI with appropriate objective sense */
      if( SCIPgetObjsense(scip) == SCIP_OBJSENSE_MAXIMIZE )
      {
         SCIP_CALL( SCIPlpiCreate(lpi, SCIPgetMessagehdlr(scip), "node LP with generated cuts", SCIP_OBJSEN_MAXIMIZE) );
      }
      else
      {
         SCIP_CALL( SCIPlpiCreate(lpi, SCIPgetMessagehdlr(scip), "node LP with generated cuts", SCIP_OBJSEN_MINIMIZE) );
      }

      /* add cols to the LP interface */
      SCIP_CALL( SCIPallocBufferArray(scip, &colobj, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &collb, ncols) );
      SCIP_CALL( SCIPallocBufferArray(scip, &colub, ncols) );
      /* gather required column information */
      for( int i = 0; i < ncols; ++i )
      {
         colobj[i] = SCIPcolGetObj(cols[i]);
         collb[i] = SCIPcolGetLb(cols[i]);
         colub[i] = SCIPcolGetUb(cols[i]);
      }
      /* add cols */
      SCIP_CALL( SCIPlpiAddCols(*lpi, ncols, colobj, collb, colub, NULL, 0, NULL, NULL, NULL) );
      SCIPfreeBufferArray(scip, &colub);
      SCIPfreeBufferArray(scip, &collb);
      SCIPfreeBufferArray(scip, &colobj);

      /* add rows to the LP interface */
      /* find number of nonzeroes in rows */
      nrownonz = 0;
      for( int i = 0; i < nrows; ++i )
      {
         assert(!(SCIPisInfinity(scip, -SCIProwGetLhs(rows[i])) && SCIPisInfinity(scip, SCIProwGetRhs(rows[i]))));
         nrownonz += SCIProwGetNLPNonz(rows[i]);
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &rowcolinds, nrownonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowvals, nrownonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowbegs, nrows + 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowlhs, nrows) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowrhs, nrows) );
      /* gather required row information */
      rowbegs[0] = 0;
      pinf = SCIPlpiInfinity(*lpi);
      ninf = -SCIPlpiInfinity(*lpi);
      for( int i = 0; i < nrows; ++i )
      {
         nrownonz = SCIProwGetNLPNonz(rows[i]);
         assert(nrownonz <= ncols);
         rowval = SCIProwGetVals(rows[i]);
         rowcols = SCIProwGetCols(rows[i]);

         rowbegs[i + 1] = rowbegs[i] + nrownonz;
         rowconst = SCIProwGetConstant(rows[i]);
         rowlhs[i] = SCIPisInfinity(scip, -SCIProwGetLhs(rows[i])) ? ninf : SCIProwGetLhs(rows[i]) - rowconst;
         rowrhs[i] = SCIPisInfinity(scip, SCIProwGetRhs(rows[i])) ? pinf : SCIProwGetRhs(rows[i]) - rowconst;

         for( int j = 0; j < nrownonz; ++j )
         {
            collppos = SCIPcolGetLPPos(rowcols[j]);
            assert(collppos >= 0);
            assert(collppos <= ncols);

            rowcolinds[rowbegs[i] + j] = collppos;
            rowvals[rowbegs[i] + j] = rowval[j];
         }
      }
      /* add rows */
      SCIP_CALL( SCIPlpiAddRows(*lpi, nrows, rowlhs, rowrhs, NULL, rowbegs[nrows], rowbegs, rowcolinds, rowvals) );

      /* get warm starting info */
      SCIP_CALL( SCIPgetLPI(scip, &wslpi) );
      SCIP_CALL( SCIPlpiGetState(wslpi, blkmem, &lpistate) );
   }
   /* if there are any cuts, then add the cuts that were not added earlier to the LPI */
   else
   {
      assert(nrows + ncuts >= sepadata->nrowsinhardcutslp);

      /* get warm starting info */
      wslpi = *lpi;
      SCIP_CALL( SCIPlpiGetState(wslpi, blkmem, &lpistate) );

      /* find number of nonzeros in cuts and allocate memory */
      nrownonz = 0;
      pinf = SCIPlpiInfinity(*lpi);
      ninf = -SCIPlpiInfinity(*lpi);
      for( int i = sepadata->nrowsinhardcutslp - nrows; i < ncuts; ++i )
      {
         assert(!(SCIPisInfinity(scip, -SCIProwGetLhs(cuts[i])) && SCIPisInfinity(scip, SCIProwGetRhs(cuts[i]))));
         assert(SCIPisInfinity(scip, -SCIProwGetLhs(cuts[i])));
         assert(!SCIPisInfinity(scip, SCIProwGetRhs(cuts[i])));
         nrownonz += SCIProwGetNNonz(cuts[i]);
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &rowcolinds, nrownonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowvals, nrownonz) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowbegs, (ncuts - sepadata->nrowsinhardcutslp + nrows) + 1) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowlhs, (ncuts - sepadata->nrowsinhardcutslp + nrows)) );
      SCIP_CALL( SCIPallocBufferArray(scip, &rowrhs, (ncuts - sepadata->nrowsinhardcutslp + nrows)) );

      /* gather required cut information */
      rowbegs[0] = 0;
      for( int i = sepadata->nrowsinhardcutslp - nrows; i < ncuts; ++i )
      {
         nrownonz = SCIProwGetNNonz(cuts[i]);
         assert(nrownonz <= ncols);
         rowval = SCIProwGetVals(cuts[i]);
         rowcols = SCIProwGetCols(cuts[i]);

         rowbegs[i - sepadata->nrowsinhardcutslp + nrows + 1] = rowbegs[i - sepadata->nrowsinhardcutslp + nrows] +
            nrownonz;
         rowconst = SCIProwGetConstant(cuts[i]);
         rowlhs[i - sepadata->nrowsinhardcutslp + nrows] = SCIPisInfinity(scip, -SCIProwGetLhs(cuts[i])) ? ninf :
            SCIProwGetLhs(cuts[i]) - rowconst;
         rowrhs[i - sepadata->nrowsinhardcutslp + nrows] = SCIPisInfinity(scip, SCIProwGetRhs(cuts[i])) ? pinf :
            SCIProwGetRhs(cuts[i]) - rowconst;

         for( int j = 0; j < nrownonz; ++j )
         {
            collppos = SCIPcolGetLPPos(rowcols[j]);
            assert(collppos >= 0);
            assert(collppos <= ncols);

            rowcolinds[rowbegs[i - sepadata->nrowsinhardcutslp + nrows] + j] = collppos;
            rowvals[rowbegs[i - sepadata->nrowsinhardcutslp + nrows] + j] = rowval[j];
         }
      }

      /* add cuts */
      SCIP_CALL( SCIPlpiAddRows(*lpi, (ncuts - sepadata->nrowsinhardcutslp + nrows), rowlhs, rowrhs, NULL,
               rowbegs[(ncuts - sepadata->nrowsinhardcutslp + nrows)], rowbegs, rowcolinds, rowvals) );
   }

   /* set warm starting basis */
   SCIP_CALL( SCIPlpiSetState(*lpi, blkmem, lpistate) );

   /* reset remaining sepadata values */
   sepadata->nrowsinhardcutslp = nrows + ncuts;

   /* free memory */
   SCIP_CALL( SCIPlpiFreeState(*lpi, blkmem, &lpistate) );
   SCIPfreeBufferArray(scip, &rowrhs);
   SCIPfreeBufferArray(scip, &rowlhs);
   SCIPfreeBufferArray(scip, &rowbegs);
   SCIPfreeBufferArray(scip, &rowvals);
   SCIPfreeBufferArray(scip, &rowcolinds);

   return SCIP_OKAY;
}

/** free separator data */
static
SCIP_RETCODE sepadataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA**       sepadata            /**< separator data structure */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(*sepadata != NULL);

   if( (*sepadata)->lpiwithhardcuts != NULL )
   {
      SCIP_CALL( SCIPlpiFree(&((*sepadata)->lpiwithhardcuts)) );
   }

   (*sepadata)->nrowsinhardcutslp = 0;
   (*sepadata)->nrunsinsoftcutslp = 0;
   (*sepadata)->ncalls = 0;
   (*sepadata)->nmaxperroundlpiters = 0;
   (*sepadata)->nmaxrootlpiters = 0;
   (*sepadata)->nrootlpiters = 0;
   (*sepadata)->nmaxtotallpiters = 0;
   (*sepadata)->ntotallpiters = 0;
   (*sepadata)->nmaxperroundcutsroot = 0;
   (*sepadata)->nmaxperroundcuts = 0;
   (*sepadata)->nmaxtotalcuts = 0;
   (*sepadata)->ntotalcuts = 0;

   SCIPfreeBlockMemory(scip, sepadata);

   return SCIP_OKAY;
}

/** update mu parameter which is used as a factor in the step length calculation; refer to the top of the file for a
 * description of the formula.
 */
/* @todo some adaptive strategy like constant after certain changes? */
static
SCIP_RETCODE updateMuSteplengthParam(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   int                   subgradientiternum, /**< subgradient iteration number */
   SCIP_Real             ubparam,            /**< estimate of the optimal Lagrangian dual value */
   SCIP_Real*            lagrangianvals,     /**< vector of Lagrangian values found so far */
   SCIP_Real             bestlagrangianval,  /**< best Lagrangian value found so far */
   SCIP_Real             avglagrangianval,   /**< rolling average of the Lagrangian values found so far */
   SCIP_Real*            muparam,            /**< mu parameter to be updated */
   SCIP_Bool*            backtrack           /**< whether mu parameter has been backtracked */
   )
{
   SCIP_Real delta;
   SCIP_Real deltaslab1ub;
   SCIP_Real deltaslab2ub;
   SCIP_Real muslab1factor;
   SCIP_Real muslab2factor;
   SCIP_Real muslab3factor;
   int maxiters;
   int i;

   *backtrack = FALSE;

   /* update the mu parameter only if it is not set to be a constant value */
   if( !sepadata->muparamconst )
   {
      delta = ubparam - bestlagrangianval;
      deltaslab1ub = MIN(sepadata->deltaslab1ub, sepadata->deltaslab2ub);
      deltaslab2ub = MAX(sepadata->deltaslab1ub, sepadata->deltaslab2ub);
      /* ensure that the ordering of different user input parameters is as expected */
      if( SCIPisPositive(scip, sepadata->muslab1factor - sepadata->muslab2factor) )
      {
         if( SCIPisPositive(scip, sepadata->muslab2factor - sepadata->muslab3factor) )
         {
            muslab1factor = sepadata->muslab1factor;
            muslab2factor = sepadata->muslab2factor;
            muslab3factor = sepadata->muslab3factor;
         }
         else
         {
            if( SCIPisPositive(scip, sepadata->muslab1factor - sepadata->muslab3factor) )
            {
               muslab1factor = sepadata->muslab1factor;
               muslab2factor = sepadata->muslab3factor;
               muslab3factor = sepadata->muslab2factor;
            }
            else
            {
               muslab1factor = sepadata->muslab3factor;
               muslab2factor = sepadata->muslab1factor;
               muslab3factor = sepadata->muslab2factor;
            }
         }
      }
      else
      {
         if( SCIPisPositive(scip, sepadata->muslab1factor - sepadata->muslab3factor) )
         {
            muslab1factor = sepadata->muslab2factor;
            muslab2factor = sepadata->muslab1factor;
            muslab3factor = sepadata->muslab3factor;
         }
         else
         {
            if( SCIPisPositive(scip, sepadata->muslab2factor - sepadata->muslab3factor) )
            {
               muslab1factor = sepadata->muslab2factor;
               muslab2factor = sepadata->muslab3factor;
               muslab3factor = sepadata->muslab1factor;
            }
            else
            {
               muslab1factor = sepadata->muslab3factor;
               muslab2factor = sepadata->muslab2factor;
               muslab3factor = sepadata->muslab1factor;
            }
         }
      }

      maxiters = MIN(sepadata->nmaxconsecitersformuupdate, sepadata->nmaxlagrangianvalsforavg);
      i = -1;

      /* if certain number of iterations are done, then check for a possibility of backtracking and apply accordingly */
      if( subgradientiternum >= maxiters )
      {
         for( i = subgradientiternum - maxiters; i < subgradientiternum; i++ )
         {
            if( SCIPisGE(scip, lagrangianvals[i], bestlagrangianval - delta) )
               break;
         }

         if( i == subgradientiternum )
         {
            *muparam *= sepadata->mubacktrackfactor;
            *backtrack = TRUE;
         }
      }

      /* update mu parameter based on the different between best and average Lagrangian values */
      if( (subgradientiternum < maxiters) || (i >= 0 && i < subgradientiternum) )
      {
         if( bestlagrangianval - avglagrangianval < deltaslab1ub * delta )
            *muparam *= muslab1factor;
         else if( bestlagrangianval - avglagrangianval < deltaslab2ub * delta )
            *muparam *= muslab2factor;
         else
            *muparam *= muslab3factor;
      }

      /* reset the mu parameter to within its bounds */
      *muparam = MAX(*muparam, sepadata->muparamlb);
      *muparam = MIN(*muparam, sepadata->muparamub);
   }

   return SCIP_OKAY;
}

/** update subgradient, i.e., residuals of generated cuts */
/* @note: assumed that \f$i^{th}\f$ cut is of the form \f${\alpha^i}^T x \leq {\alpha^i_0}\f$ */
static
void updateSubgradient(
   SCIP*                 scip,                     /**< SCIP data structure */
   SCIP_SOL*             sol,                      /**< LP solution used in updating subgradient vector */
   SCIP_ROW**            cuts,                     /**< cuts generated so far */
   int                   ncuts,                    /**< number of cuts generated so far */
   SCIP_Real*            subgradient,              /**< vector of subgradients to be updated */
   SCIP_Real*            dualvector,               /**< Lagrangian multipliers */
   SCIP_Bool*            subgradientzero,          /**< whether the subgradient vector is all zero */
   int*                  ncutviols,                /**< number of violations of generated cuts */
   SCIP_Real*            maxcutviol,               /**< maximum violation of generated cuts */
   int*                  nnzsubgradientdualprod,   /**< number of nonzero products of subgradient vector and Lagrangian multipliers (i.e.,
                                                     number of complementarity slackness violations) */
   SCIP_Real*            maxnzsubgradientdualprod  /**< maximum value of nonzero products of subgradient vector and Lagrangian multipliers
                                                     (i.e., maximum value of complementarity slackness violations) */
   )
{
   int nzerosubgradient;
   SCIP_Real term;

   assert(subgradientzero != NULL);
   assert(ncutviols != NULL);
   assert(maxcutviol != NULL);
   assert(nnzsubgradientdualprod != NULL);
   assert(maxnzsubgradientdualprod != NULL);

   *ncutviols = 0;
   *maxcutviol = 0.0;
   *nnzsubgradientdualprod = 0;
   *maxnzsubgradientdualprod = 0.0;
   nzerosubgradient = 0;
   *subgradientzero = FALSE;

   /* for each cut, calculate the residual along with various violation metrics */
   for( int i = 0; i < ncuts; i++ )
   {
      assert(SCIPisInfinity(scip, -SCIProwGetLhs(cuts[i])));
      assert(!SCIPisInfinity(scip, SCIProwGetRhs(cuts[i])));
      subgradient[i] = SCIPgetRowSolActivity(scip, cuts[i], sol) + SCIProwGetConstant(cuts[i]) - SCIProwGetRhs(cuts[i]);
      if( SCIPisFeasZero(scip, subgradient[i]) )
      {
         subgradient[i] = 0.0;
         nzerosubgradient++;
      }
      else
      {
         /* check for cut violation */
         if( SCIPisFeasPositive(scip, subgradient[i]) )
         {
            (*ncutviols)++;
            *maxcutviol = MAX(*maxcutviol, subgradient[i]);
         }

         /* check for violation of complementarity slackness associated with the cut */
         if( !SCIPisZero(scip, subgradient[i] * dualvector[i]) )
         {
            (*nnzsubgradientdualprod)++;
            term = REALABS(subgradient[i] * dualvector[i]);
            *maxnzsubgradientdualprod = MAX(*maxnzsubgradientdualprod, term);
         }
      }
   }

   /* indicator for all zero subgradient vector */
   if( nzerosubgradient == ncuts )
   {
      *subgradientzero = TRUE;
   }
}

/** update Lagrangian value, i.e., optimal value of the Lagrangian dual with fixed multipliers */
static
void updateLagrangianValue(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             objval,             /**< objective value of the Lagrangian dual with fixed multipliers */
   SCIP_Real*            dualvector,         /**< Lagrangian multipliers */
   SCIP_ROW**            cuts,               /**< cuts generated so far */
   int                   ncuts,              /**< number of cuts generated so far */
   SCIP_Real*            lagrangianval       /**< Lagrangian value to be updated */
   )
{
   *lagrangianval = objval;

   for( int i = 0; i < ncuts; i++ )
   {
      assert(SCIPisInfinity(scip, -SCIProwGetLhs(cuts[i])));
      assert(!SCIPisInfinity(scip, SCIProwGetRhs(cuts[i])));
      *lagrangianval += dualvector[i] * (SCIProwGetConstant(cuts[i]) - SCIProwGetRhs(cuts[i]));
   }
}

/** update step length based on various input arguments; refer to the top of the file for an expression */
static
void updateStepLength(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             muparam,            /**< mu parameter used as a factor for step length */
   SCIP_Real             ubparam,            /**< estimate of the optimal Lagrangian dual value */
   SCIP_Real             lagrangianval,      /**< Lagrangian value */
   SCIP_Real*            subgradient,        /**< subgradient vector */
   int                   ncuts,              /**< number of cuts generated so far */
   SCIP_Real*            steplength          /**< step length to be updated */
   )
{
   SCIP_Real normsquared = 0.0;

   for( int i = 0; i < ncuts; i++ )
   {
      normsquared += SQR(subgradient[i]);
   }

   if( !SCIPisFeasZero(scip, normsquared) )
   {
      *steplength = (muparam * (ubparam - lagrangianval))/(normsquared); /*lint !e795*/
   }
}

/** update the ball radius (based on various violation metrics) that is used for stabilization of Lagrangian multipliers */
static
void updateBallRadius(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_Real             maxviolscore,       /**< weighted average of maximum value of generated cut violations and maximum value of
                                               complementarity slackness violations, in the current iteration */
   SCIP_Real             maxviolscoreold,    /**< weighted average of maximum value of generated cut violations and maximum value of
                                               complementarity slackness violations, in the previous iteration */
   SCIP_Real             nviolscore,         /**< weighted average of number of generated cut violations and number of complementarity
                                               slackness violations, in the current iteration */
   SCIP_Real             nviolscoreold,      /**< weighted average of number of generated cut violations and number of complementarity
                                               slackness violations, in the previous iteration */
   int                   nlpiters,           /**< number of LP iterations taken for solving the Lagrangian dual with fixed multipliers in
                                               current iteration */
   SCIP_Real*            ballradius          /**< norm ball radius to be updated */
   )
{
   SCIP_Bool maxviolscoreimproved;
   SCIP_Bool nviolscoreimproved;

   assert(ballradius != NULL);

   maxviolscoreimproved = !SCIPisNegative(scip, maxviolscoreold - maxviolscore);
   nviolscoreimproved = !SCIPisNegative(scip, nviolscoreold - nviolscore);

   if( maxviolscoreimproved && nviolscoreimproved )
   {
      /* both the maximum violation and number of violations scores have become better, so, increase the radius */
      if( sepadata->optimalfacepriority <= 1 )
      {
         *ballradius *= 2.0;
         *ballradius = MIN(*ballradius, sepadata->radiusmax);
      }
      else
      {
         *ballradius *= 1.5;
         *ballradius = MIN(*ballradius, sepadata->radiusmax/2.0);
      }
   }
   else if( !maxviolscoreimproved && !nviolscoreimproved )
   {
      /* both the maximum violation and number of violations scores have become worse, so, decrease the radius */
      *ballradius *= 0.5;
      *ballradius = MAX(*ballradius, sepadata->radiusmin);
   }
   else if( nlpiters == 0 )
   {
      /* only one among the maximum violation and number of violations scores has become better, and the LP basis did
       * not change (i.e., nlpters = 0), so, increase the radius slightly */
      if( sepadata->optimalfacepriority <= 1 )
      {
         *ballradius *= 1.5;
         *ballradius = MIN(*ballradius, sepadata->radiusmax);
      }
      else
      {
         *ballradius *= 1.2;
         *ballradius = MIN(*ballradius, sepadata->radiusmax/2.0);
      }
   }
}

/** projection of Lagrangian multipliers onto L1-norm ball. This algorithm is based on the following article.
 *
 * Condat L. (2016).@n
 * Fast projection onto the simplex and the \f$l_1\f$ ball.@n
 * Mathematical Programming, 158, 1-2, 575-585.
 *
 */
static
SCIP_RETCODE l1BallProjection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            dualvector,         /**< Lagrangian multipliers to be projected onto L1-norm vall */
   int                   dualvectorlen,      /**< length of the Lagrangian multipliers vector */
   SCIP_Real             radius              /**< radius of the L1-norm ball */
   )
{
   SCIP_Real* temp1vals;
   SCIP_Real* temp2vals;
   SCIP_Real pivotparam;
   SCIP_Real val;
   SCIP_Real term;
   SCIP_Bool temp1changed;
   int ntemp1removed;
   int* temp1inds;
   int* temp2inds;
   int temp1len;
   int temp2len;

   assert(!SCIPisNegative(scip, radius));
   val = REALABS(dualvector[0]);
   /* calculate the L1-norm of the Lagrangian multipliers */
   for( int i = 1; i < dualvectorlen; i++ )
   {
      val += REALABS(dualvector[i]);
   }

   /* if the vector of Lagrangian multipliers lies outside the L1-norm ball, then do the projection */
   if( SCIPisGT(scip, val, radius) )
   {
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &temp1vals, dualvectorlen) );
      SCIP_CALL( SCIPallocCleanBufferArray(scip, &temp2vals, dualvectorlen) );
      SCIP_CALL( SCIPallocBufferArray(scip, &temp1inds, dualvectorlen) );
      SCIP_CALL( SCIPallocBufferArray(scip, &temp2inds, dualvectorlen) );
      for( int i = 0; i < dualvectorlen; i++ )
      {
         temp1inds[i] = -1;
         temp2inds[i] = -1;
      }
      temp2len = 0;

      temp1vals[0] = REALABS(dualvector[0]);
      temp1inds[0] = 0;
      temp1len = 1;
      pivotparam = REALABS(dualvector[0]) - radius;

      for( int i = 1; i < dualvectorlen; i++ )
      {
         if( SCIPisGT(scip, REALABS(dualvector[i]), pivotparam) )
         {
            pivotparam += ((REALABS(dualvector[i]) - pivotparam) / (temp1len + 1));
            if( SCIPisGT(scip, pivotparam, REALABS(dualvector[i]) - radius) )
            {
               temp1vals[temp1len] = REALABS(dualvector[i]);
               temp1inds[temp1len] = i;
               temp1len++;
            }
            else
            {
               for( int j = 0; j < temp1len; j++ )
               {
                  temp2vals[temp2len + j] = temp1vals[j];
                  temp2inds[temp2len + j] = temp1inds[j];
               }
               temp2len += temp1len;
               temp1vals[0] = REALABS(dualvector[i]);
               temp1inds[0] = i;
               temp1len = 1;
               pivotparam = REALABS(dualvector[i]) - radius;
            }
         }
      }

      for( int i = 0; i < temp2len; i++ )
      {
         if( SCIPisGT(scip, temp2vals[i], pivotparam) )
         {
            temp1vals[temp1len] = temp2vals[i];
            temp1inds[temp1len] = temp2inds[i];
            temp1len++;
            pivotparam += ((temp2vals[i] - pivotparam) / temp1len);
         }
      }

      temp1changed = TRUE;
      ntemp1removed = 0;
      while( temp1changed )
      {
         temp1changed = FALSE;

         for( int i = 0; i < temp1len; i++ )
         {
            /* @note: the third condition (temp1len - ntemp1removed > 0) is true as long as the first condition
             * (temp1inds[i] >= 0) is true.
             */
            if( (temp1inds[i] >= 0) && SCIPisLE(scip, temp1vals[i], pivotparam) )
            {
               temp1inds[i] = -1;
               temp1changed = TRUE;
               ntemp1removed++;
               assert(temp1len - ntemp1removed > 0);
               /* coverity[divide_by_zero] */
               pivotparam += ((pivotparam - temp1vals[i]) / (temp1len - ntemp1removed));
            }
         }
      }

      for( int i = 0; i < dualvectorlen; i++ )
      {
         term = REALABS(dualvector[i]);
         val = MAX(term - pivotparam, 0.0);

         if( SCIPisPositive(scip, dualvector[i]) )
         {
            dualvector[i] = val;
         }
         else if( SCIPisNegative(scip, dualvector[i]) )
         {
            dualvector[i] = -val;
         }
      }

      /* free memory */
      for( int i = 0; i < dualvectorlen; i++ )
      {
         temp2vals[i] = 0.0;
         temp1vals[i] = 0.0;
      }
      SCIPfreeBufferArray(scip, &temp2inds);
      SCIPfreeBufferArray(scip, &temp1inds);
      SCIPfreeCleanBufferArray(scip, &temp2vals);
      SCIPfreeCleanBufferArray(scip, &temp1vals);
   }

   return SCIP_OKAY;
}

/** projection of Lagrangian multipliers onto L2-norm ball */
static
void l2BallProjection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            dualvector,         /**< Lagrangian multipliers to be projected onto L2-norm vall */
   int                   dualvectorlen,      /**< length of the Lagrangian multipliers vector */
   SCIP_Real             radius              /**< radius of the L2-norm ball */
   )
{
   SCIP_Real l2norm;
   SCIP_Real factor;

   assert(!SCIPisNegative(scip, radius));

   l2norm = 0.0;
   /* calculate the L2-norm of the Lagrangian multipliers */
   for( int i = 0; i < dualvectorlen; i++ )
   {
      l2norm += SQR(dualvector[i]);
   }
   l2norm = sqrt(l2norm);
   factor = radius/(1.0 + l2norm);

   /* if the vector of Lagrangian multipliers is outside the L2-norm ball, then do the projection */
   if( SCIPisGT(scip, l2norm, radius) && SCIPisLT(scip, factor, 1.0) )
   {
      for( int i = 0; i < dualvectorlen; i++ )
      {
         dualvector[i] *= factor;
      }
   }
}

/** projection of Lagrangian multipliers onto L_infinity-norm ball */
static
void linfBallProjection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            dualvector,         /**< Lagrangian multipliers to be projected onto L_infinity-norm vall */
   int                   dualvectorlen,      /**< length of the Lagrangian multipliers vector */
   SCIP_Real             radius              /**< radius of the L_infinity-norm ball */
   )
{
   assert(!SCIPisNegative(scip, radius));

   /* if the vector of Lagrangian multipliers is outside the L_infinity-norm ball, then do the projection */
   for( int i = 0; i < dualvectorlen; i++ )
   {
      if( SCIPisLT(scip, dualvector[i], -radius) )
      {
         dualvector[i] = -radius;
      }
      else if( SCIPisGT(scip, dualvector[i], radius) )
      {
         dualvector[i] = radius;
      }
   }
}

/** weighted Lagrangian multipliers based on a given vector as stability center */
/* @todo calculate weight outside this function and pass it (so that this function becomes generic and independent of
 * the terminology related to best Lagrangian multipliers)
 */
static
SCIP_RETCODE weightedDualVector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_Real*            dualvector,         /**< Lagrangian multipliers */
   int                   dualvectorlen,      /**< length of the Lagrangian multipliers vector */
   SCIP_Real*            stabilitycenter,    /**< stability center (i.e., core vector of Lagrangian multipliers) */
   int                   stabilitycenterlen, /**< length of the stability center */
   int                   nbestdualupdates,   /**< number of best Lagrangian values found so far */
   int                   totaliternum        /**< total number of iterations of the relax-and-cut algorithm performed so far */
   )
{
   SCIP_Real constant;
   SCIP_Real weight;
   SCIP_Real alpha;

   constant = MAX(2.0, sepadata->constant);
   /* weight factor from the literature on Dantzig-Wolfe decomposition stabilization schemes */
   weight = MIN(constant, (totaliternum + 1 + nbestdualupdates) / 2.0);
   alpha = 1.0 / weight;

   assert(dualvectorlen >= stabilitycenterlen);

   /* weighted Lagrangian multipliers */
   for( int i = 0; i < stabilitycenterlen; i++ )
   {
      dualvector[i] = alpha * dualvector[i] + (1 - alpha) * stabilitycenter[i];
   }
   for( int i = stabilitycenterlen; i < dualvectorlen; i++ )
   {
      dualvector[i] = alpha * dualvector[i];
   }

   return SCIP_OKAY;
}

/** stabilize Lagrangian multipliers */
static
SCIP_RETCODE stabilizeDualVector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_Real*            dualvector,         /**< Lagrangian multipliers */
   int                   dualvectorlen,      /**< length of the Lagrangian multipliers vector */
   SCIP_Real*            bestdualvector,     /**< best Lagrangian multipliers found so far */
   int                   bestdualvectorlen,  /**< length of the best Lagrangian multipliers vector */
   int                   nbestdualupdates,   /**< number of best Lagrangian values found so far */
   int                   subgradientiternum, /**< iteration number of the subgradient algorithm */
   int                   totaliternum,       /**< total number of iterations of the relax-and-cut algorithm performed so far */
   SCIP_Real             maxviolscore,       /**< weighted average of maximum value of generated cut violations and maximum value of
                                               complementarity slackness violations, in the current iteration */
   SCIP_Real             maxviolscoreold,    /**< weighted average of maximum value of generated cut violations and maximum value of
                                               complementarity slackness violations, in the previous iteration */
   SCIP_Real             nviolscore,         /**< weighted average of number of generated cut violations and number of complementarity
                                               slackness violations, in the current iteration */
   SCIP_Real             nviolscoreold,      /**< weighted average of number of generated cut violations and number of complementarity
                                               slackness violations, in the previous iteration */
   int                   nlpiters,           /**< number of LP iterations taken for solving the Lagrangian dual with fixed multipliers in
                                               current iteration */
   SCIP_Real*            ballradius          /**< norm ball radius */
   )
{
   if( sepadata->projectiontype > 0 )
   {
      if( subgradientiternum >= 1 )
      {
         /* update the ball radius */
         updateBallRadius(scip, sepadata, maxviolscore, maxviolscoreold, nviolscore, nviolscoreold, nlpiters,
               ballradius);
      }

      if( sepadata->projectiontype == 1 )
      {
         /* projection of Lagrangian multipliers onto L1-norm ball */
         SCIP_CALL( l1BallProjection(scip, dualvector, dualvectorlen, *ballradius) );
      }
      else if( sepadata->projectiontype == 2 )
      {
         /* projection of Lagrangian multipliers onto L2-norm ball */
         l2BallProjection(scip, dualvector, dualvectorlen, *ballradius);
      }
      else if( sepadata->projectiontype == 3 )
      {
         /* projection of Lagrangian multipliers onto L_inf-norm ball */
         linfBallProjection(scip, dualvector, dualvectorlen, *ballradius);
      }
   }

   if( sepadata->stabilitycentertype == 1 )
   {
      /* weighted Lagrangian multipliers based on best Langrangian multipliers as stability center */
      SCIP_CALL( weightedDualVector(scip, sepadata, dualvector, dualvectorlen, bestdualvector,
               bestdualvectorlen, nbestdualupdates, totaliternum) );
   }

   return SCIP_OKAY;
}

/** update Lagrangian multipliers */
static
SCIP_RETCODE updateDualVector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_Real*            dualvector1,        /**< Lagrangian multipliers vector to be updated */
   SCIP_Real*            dualvector2,        /**< Lagrangian multipliers vector used for backtracking */
   int                   dualvector2len,     /**< length of the Lagrangian multipliers vector used for backtracking */
   int                   ndualvector2updates,/**< number of best Lagrangian values found so far */
   int                   subgradientiternum, /**< iteration number of the subgradient algorithm */
   int                   totaliternum,       /**< total number of iterations of the relax-and-cut algorithm performed so far */
   SCIP_Real             steplength,         /**< step length used for updating Lagrangian multipliers */
   SCIP_Real*            subgradient,        /**< subgradient vector */
   int                   ncuts,              /**< number of generated cuts so far */
   SCIP_Bool             backtrack,          /**< whether the Lagrangian multipliers need to be backtracked */
   SCIP_Real             maxviolscore,       /**< weighted average of maximum value of generated cut violations and maximum value of
                                               complementarity slackness violations, in the current iteration */
   SCIP_Real             maxviolscoreold,    /**< weighted average of maximum value of generated cut violations and maximum value of
                                               complementarity slackness violations, in the previous iteration */
   SCIP_Real             nviolscore,         /**< weighted average of number of generated cut violations and number of complementarity
                                               slackness violations, in the current iteration */
   SCIP_Real             nviolscoreold,      /**< weighted average of number of generated cut violations and number of complementarity
                                               slackness violations, in the previous iteration */
   int                   nlpiters,           /**< number of LP iterations taken for solving the Lagrangian dual with fixed multipliers in
                                               current iteration */
   SCIP_Bool*            dualvecsdiffer,     /**< whether the updated Lagrangian multipliers differ from the old one */
   SCIP_Real*            ballradius          /**< norm ball radius */
   )
{
   SCIP_Real* dualvector1copy;

   assert(dualvector2len <= ncuts);
   assert(dualvecsdiffer != NULL);

   *dualvecsdiffer = FALSE;
   /* @todo do allocation and free operations outside only once instead of every time this function is called? */
   /* copy of the Lagrangian multipliers to be used to check if the updated vector is different than this */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &dualvector1copy, ncuts) );
   for( int i = 0; i < ncuts; i++ )
   {
      dualvector1copy[i] = dualvector1[i];
   }

   /* if backtracking was not identified at the time of the mu parameter update, then update the Lagrangian multipliers
    * based on the given subgradient vector
    */
   if( !backtrack )
   {
      assert((subgradient != NULL) || (ncuts == 0));
      assert(subgradientiternum >= 0);

      /* update Lagrangian multipliers */
      for( int i = 0; i < ncuts; i++ )
      {
         dualvector1[i] += steplength * subgradient[i];
      }

      /* projection onto non-negative orthant */
      for( int i = 0; i < ncuts; i++ )
      {
         dualvector1[i] = MAX(dualvector1[i], 0.0);
      }

      /* stabilization of Lagrangian multipliers */
      SCIP_CALL( stabilizeDualVector(scip, sepadata, dualvector1, ncuts, dualvector2, dualvector2len,
               ndualvector2updates, subgradientiternum, totaliternum, maxviolscore, maxviolscoreold, nviolscore,
               nviolscoreold, nlpiters, ballradius) );

      /* projection onto non-negative orthant again in case stabilization changed some components negative*/
      for( int i = 0; i < ncuts; i++ )
      {
         dualvector1[i] = MAX(dualvector1[i], 0.0);
      }
   }
   /* if backtracking was identified at the time of the mu parameter update, then backtrack the Lagrangian multipliers
    * based on the given backtracking multipliers
    */
   else
   {
      for( int i = 0; i < dualvector2len; i++ )
      {
         dualvector1[i] = dualvector2[i];
      }

      for( int i = dualvector2len; i < ncuts; i++ )
      {
         dualvector1[i] = 0.0;
      }
   }

   /* identify if the vector of Lagrangian multipliers is indeed different after updating */
   for( int i = 0; i < ncuts; i++ )
   {
      if( !SCIPisEQ(scip, dualvector1[i], dualvector1copy[i]) )
      {
         *dualvecsdiffer = TRUE;
         break;
      }
   }

   /* free memory */
   for( int i = 0; i < ncuts; i++ )
   {
      dualvector1copy[i] = 0.0;
   }
   SCIPfreeCleanBufferArray(scip, &dualvector1copy);

   return SCIP_OKAY;
}

/** check different termination criteria */
/* @note: the criterion based on objvecsdiffer assumes deterministic solving process (i.e., we would get same LP solution
 * for "Lagrangian dual with fixed Lagrangian multipliers" when the objective vector remains the same across iterations).
 */
/* @todo nlpssolved criterion? */
static
SCIP_RETCODE checkLagrangianDualTermination(
   SCIP_SEPADATA*        sepadata,                 /**< separator data structure */
   int                   nnewaddedsoftcuts,        /**< number of cuts that were recently penalized and added to the Lagrangian dual's
                                                     objective function */
   int                   nyettoaddsoftcuts,        /**< number of cuts that are yet to be penalized and added to the Lagrangian dual's
                                                     objective function */
   SCIP_Bool             objvecsdiffer,            /**< whether the Lagrangian dual's objective function has changed */
   int                   ngeneratedcurrroundcuts,  /**< number of cuts generated in the current separation round */
   int                   nmaxgeneratedperroundcuts,/**< maximal number of cuts allowed to generate per separation round */
   int                   ncurrroundlpiters,        /**< number of separating LP iterations in the current separation round */
   int                   depth,                    /**< depth of the current node */
   SCIP_Bool*            terminate                 /**< whether to terminate the subgradient algorithm loop */
   )
{
   *terminate = FALSE;

   /* check if no new cuts were added to the Lagrangian dual, no cuts are remaining to be added, and the objective
    * function of the Lagrangian dual with fixed multipliers had not changed from the previous iteration
    */
   if( (nnewaddedsoftcuts == 0) && (nyettoaddsoftcuts == 0) && !objvecsdiffer )
      *terminate = TRUE;

   /* check if allowed number of cuts in this separation round have already been generated */
   if( ngeneratedcurrroundcuts >= nmaxgeneratedperroundcuts )
      *terminate = TRUE;

   /* check if allowed number of cuts in the tree have already been generated */
   if( sepadata->ntotalcuts >= sepadata->nmaxtotalcuts )
      *terminate = TRUE;

   /* check if allowed number of simplex iterations in this separation round have already been used up */
   if( (sepadata->nmaxperroundlpiters >= 0) && (ncurrroundlpiters >= sepadata->nmaxperroundlpiters) )
      *terminate = TRUE;

   /* check if allowed number of simplex iterations in the root node have already been used up */
   if( (depth == 0) && (sepadata->nmaxrootlpiters >= 0) && (sepadata->nrootlpiters >= sepadata->nmaxrootlpiters) )
      *terminate = TRUE;

   /* check if allowed number of simplex iterations in the tree have already been used up */
   if( (sepadata->nmaxtotallpiters >= 0) && (sepadata->ntotallpiters >= sepadata->nmaxtotallpiters) )
      *terminate = TRUE;

   return SCIP_OKAY;
}

/** solve the LP corresponding to the Lagrangian dual with fixed Lagrangian multipliers */
static
SCIP_RETCODE solveLagromoryLP(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   int                   depth,              /**< depth of the current node in the tree */
   SCIP_Real             origobjoffset,      /**< objective offset in the current node's relaxation */
   SCIP_Bool*            solfound,           /**< whether an LP optimal solution has been found */
   SCIP_SOL*             sol,                /**< data structure to store LP optimal solution, if found */
   SCIP_Real*            solvals,            /**< values of the LP optimal solution, if found */
   SCIP_Real*            objval,             /**< optimal objective value of the LP optimal solution, if found */
   int*                  ncurrroundlpiters   /**< number of LP iterations taken for solving Lagrangian dual problems with fixed multipliers
                                               in the current separator round */
   )
{
   SCIP_Real timelimit;
   SCIP_COL** cols;
   SCIP_COL* col;
   SCIP_VAR* var;
   SCIP_LPI* lpi;
   SCIP_Bool lperror;
   SCIP_Bool cutoff;
   SCIP_LPSOLSTAT stat;
   SCIP_Longint ntotallpiters;
   SCIP_Longint nlpiters;
   int ncols;
   int iterlimit;

   assert(solfound != NULL);
   assert(sol != NULL);
   assert(solvals != NULL);
   assert(ncurrroundlpiters != NULL);
   assert(objval != NULL);

   *solfound = FALSE;
   lperror = FALSE;
   cutoff = FALSE;
   iterlimit = -1;
   lpi = sepadata->lpiwithsoftcuts;

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* set time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit <= 0.0 )
      {
         SCIPdebugMsg(scip, "skip Lagromory cut generation since no time left\n");
         goto TERMINATE;
      }
      /* @note: the following direct LPI call is being used because of the lack of an equivalent function call in
       * scip_lp.c (lpSetRealpar exists in lp.c though)
       */
      SCIP_CALL( SCIPlpiSetRealpar(lpi, SCIP_LPPAR_LPTILIM, timelimit) );
   }

   /* find iteration limit */
   if( (depth == 0) &&
         (sepadata->perrootlpiterfactor >= 0.0 && !SCIPisInfinity(scip, sepadata->perrootlpiterfactor)) )
   {
      iterlimit = (int)(sepadata->perrootlpiterfactor * SCIPgetNRootFirstLPIterations(scip));
   }
   else if( (depth > 0) &&
         (sepadata->perlpiterfactor >= 0.0 && !SCIPisInfinity(scip, sepadata->perlpiterfactor)) )
   {
      iterlimit = (int)(sepadata->perlpiterfactor * SCIPgetNNodeInitLPIterations(scip));
   }
   if( sepadata->nmaxperroundlpiters >= 0 )
   {
      if( sepadata->nmaxperroundlpiters - *ncurrroundlpiters >= 0 )
      {
         if( iterlimit >= 0 )
         {
            iterlimit = MIN(iterlimit, sepadata->nmaxperroundlpiters - *ncurrroundlpiters);
         }
         else
         {
            iterlimit = sepadata->nmaxperroundlpiters - *ncurrroundlpiters;
         }
      }
      else
      {
         iterlimit = 0;
      }
   }
   /* @todo impose a finite iteration limit only when the dualvector changes from zero to non-zero for the first time because
    * many simplex pivots are performed in this case even with warm starting (compared to the case when the
    * dualvector changes from non-zero to non-zero).
    */

   /* solve the LP with an iteration limit and get number of simplex iterations taken */
   ntotallpiters = SCIPgetNLPIterations(scip);

   SCIP_CALL( SCIPsolveDiveLP(scip, iterlimit, &lperror, &cutoff) );

   nlpiters = SCIPgetNLPIterations(scip) - ntotallpiters;

   /* get the solution and objective value if optimal */
   stat = SCIPgetLPSolstat(scip);
   /* @todo is there any way to accept terminations due to iterlimit and timelimit as well? It is not possible
    * currently because primal sol is not saved in these cases.
    */
   /* @note: ideally, only primal feasibility is sufficient. But, there is no such option with SCIPgetLPSolstat. */
   if( stat == SCIP_LPSOLSTAT_OPTIMAL )
   {
      if( SCIPisLPSolBasic(scip) )
      {
         *solfound = TRUE;

         /* update sol */
         for( int i = 0; i < ncols; ++i )
         {
            col = cols[i];
            assert(col != NULL);

            var = SCIPcolGetVar(col);

            solvals[i] = SCIPcolGetPrimsol(col);
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, solvals[i]) );
         }

         *objval = SCIPgetLPObjval(scip);
         *objval = *objval + origobjoffset;
      }
   }

   /* update some statistics */
   if( depth == 0 )
   {
      sepadata->nrootlpiters += (int)nlpiters;
   }
   sepadata->ntotallpiters += (int)nlpiters;
   *ncurrroundlpiters += (int)nlpiters;

TERMINATE:
   return SCIP_OKAY;
}

/** solve the LP corresponding to the node relaxation upon adding all the generated cuts */
static
SCIP_RETCODE solveLPWithHardCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_Bool*            solfound,           /**< whether an LP optimal solution has been found */
   SCIP_SOL*             sol,                /**< data structure to store LP optimal solution, if found */
   SCIP_Real*            solvals             /**< values of the LP optimal solution, if found */
   )
{
   SCIP_Real timelimit;
   SCIP_COL** cols;
   SCIP_COL* col;
   SCIP_VAR* var;
   int ncols;

   assert(solfound != NULL);
   assert(sol != NULL);
   assert(solvals != NULL);

   *solfound = FALSE;

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* set time limit */
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
   if( !SCIPisInfinity(scip, timelimit) )
   {
      timelimit -= SCIPgetSolvingTime(scip);
      if( timelimit <= 0.0 )
      {
         SCIPdebugMsg(scip, "skip Lagromory cut generation since no time left\n");
         goto TERMINATE;
      }
      SCIP_CALL( SCIPlpiSetRealpar(sepadata->lpiwithhardcuts, SCIP_LPPAR_LPTILIM, timelimit) );
   }

   /* solve the LP */
   SCIP_CALL( SCIPlpiSolvePrimal(sepadata->lpiwithhardcuts) );

   /* get the solution if primal feasible */
   if( SCIPlpiIsPrimalFeasible(sepadata->lpiwithhardcuts) )
   {
      *solfound = TRUE;
      SCIP_CALL( SCIPlpiGetSol(sepadata->lpiwithhardcuts, NULL, solvals, NULL, NULL, NULL) );

      /* update sol */
      for( int i = 0; i < ncols; ++i )
      {
         col = cols[i];
         assert(col != NULL);

         var = SCIPcolGetVar(col);

         SCIP_CALL( SCIPsetSolVal(scip, sol, var, solvals[i]) );
      }
   }

TERMINATE:
   return SCIP_OKAY;
}

/** construct a cut based on the input cut coefficients, sides, etc */
static
SCIP_RETCODE constructCutRow(
   SCIP*                 scip,                     /**< SCIP data structure */
   SCIP_SEPA*            sepa,                     /**< pointer to the separator */
   SCIP_SEPADATA*        sepadata,                 /**< separator data structure */
   int                   mainiternum,              /**< iteration number of the outer loop of the relax-and-cut algorithm */
   int                   subgradientiternum,       /**< iteration number of the subgradient algorithm */
   int                   cutnnz,                   /**< number of nonzeros in cut */
   int*                  cutinds,                  /**< column indices in cut */
   SCIP_Real*            cutcoefs,                 /**< cut cofficients */
   SCIP_Real             cutefficacy,              /**< cut efficacy */
   SCIP_Real             cutrhs,                   /**< RHS of cut */
   SCIP_Bool             cutislocal,               /**< whether cut is local */
   int                   cutrank,                  /**< rank of cut */
   SCIP_ROW**            generatedcuts,            /**< array of generated cuts */
   SCIP_Real*            generatedcutefficacies,   /**< array of generated cut efficacies w.r.t. the respective LP bases used for cut
                                                     generations */
   int                   ngeneratedcurrroundcuts,  /**< number of cuts generated until the previous basis in the current separation round */
   int*                  ngeneratednewcuts,        /**< number of new cuts generated using the current basis */
   SCIP_Bool*            cutoff                    /**< should the current node be cutoff? */
   )
{
   SCIP_COL** cols;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real lhs;
   SCIP_Real rhs;

   assert(scip != NULL);
   assert(mainiternum >= 0);
   assert(ngeneratednewcuts != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;

   if( cutnnz == 0 && SCIPisFeasNegative(scip, cutrhs) ) /*lint !e644*/
   {
      SCIPdebugMsg(scip, " -> Lagromory cut detected node infeasibility with cut 0 <= %g.\n", cutrhs);
      *cutoff = TRUE;
      return SCIP_OKAY;
   }

   /* only take efficient cuts */
   if( SCIPisEfficacious(scip, cutefficacy) )
   {
      SCIP_ROW* cut;
      SCIP_VAR* var;
      char cutname[SCIP_MAXSTRLEN];
      int v;

      /* construct cut name */
      if( subgradientiternum >= 0 )
      {
         (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_%" SCIP_LONGINT_FORMAT "_%d" "_%d" "_%d", SCIPsepaGetName(sepa),
               sepadata->ncalls, mainiternum, subgradientiternum, ngeneratedcurrroundcuts + *ngeneratednewcuts);
      }
      else
      {
         (void) SCIPsnprintf(cutname, SCIP_MAXSTRLEN, "%s_%" SCIP_LONGINT_FORMAT "_%d" "_%d", SCIPsepaGetName(sepa),
               sepadata->ncalls, mainiternum, ngeneratedcurrroundcuts + *ngeneratednewcuts);
      }

      /* create empty cut */
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &cut, sepa, cutname, -SCIPinfinity(scip), cutrhs, cutislocal, FALSE,
               sepadata->dynamiccuts) );

      /* set cut rank */
      SCIProwChgRank(cut, cutrank); /*lint !e644*/

      /* cache the row extension and only flush them if the cut gets added */
      SCIP_CALL( SCIPcacheRowExtensions(scip, cut) );

      cols = SCIPgetLPCols(scip);

      /* collect all non-zero coefficients */
      for( v = 0; v < cutnnz; ++v )
      {
         var = SCIPcolGetVar(cols[cutinds[v]]);
         SCIP_CALL( SCIPaddVarToRow(scip, cut, var, cutcoefs[v]) );
      }

      /* flush all changes before adding the cut */
      SCIP_CALL( SCIPflushRowExtensions(scip, cut) );

      if( SCIProwGetNNonz(cut) == 0 )
      {
         assert( SCIPisFeasNegative(scip, cutrhs) );
         SCIPdebugMsg(scip, " -> Lagromory cut detected node infeasibility with cut 0 <= %g.\n", cutrhs);
         *cutoff = TRUE;
         return SCIP_OKAY;
      }
      else
      {
         /* gathering lhs and rhs again in case the separator is extended later to add other cuts/constraints that may
          * have non-inf lhs or inf rhs */
         lhs = SCIProwGetLhs(cut);
         rhs = SCIProwGetRhs(cut);
         assert(SCIPisInfinity(scip, -lhs));
         assert(!SCIPisInfinity(scip, rhs));

         SCIPdebugMsg(scip, " -> %s cut <%s>: rhs=%f, eff=%f\n", "lagromory", cutname, cutrhs, cutefficacy);

         /* check if the cut leads to infeasibility (i.e., *cutoff = TRUE) */
         {
            /* modifiable cuts cannot be declared infeasible, since we don't know all coefficients */
            if( SCIProwIsModifiable(cut) )
               *cutoff = FALSE;

            /* check for activity infeasibility */
            minactivity = SCIPgetRowMinActivity(scip, cut);
            maxactivity = SCIPgetRowMaxActivity(scip, cut);

            if( (!SCIPisInfinity(scip,  rhs) && SCIPisFeasPositive(scip, minactivity - rhs)) ||
                  (!SCIPisInfinity(scip, -lhs) && SCIPisFeasNegative(scip, maxactivity - lhs)) )
            {
               SCIPdebugMsg(scip, "cut <%s> is infeasible (sides=[%g,%g], act=[%g,%g])\n",
                     SCIProwGetName(cut), lhs, rhs, minactivity, maxactivity);
               *cutoff = TRUE;
            }
         }

         /* store the newly generated cut in an array and update some statistics */
         generatedcuts[ngeneratedcurrroundcuts + *ngeneratednewcuts] = cut;
         generatedcutefficacies[ngeneratedcurrroundcuts + *ngeneratednewcuts] = cutefficacy;
         ++(*ngeneratednewcuts);
      }
   }

   return SCIP_OKAY;
}

/** aggregated generated cuts based on the best Lagrangian multipliers */
static
SCIP_RETCODE aggregateGeneratedCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< pointer to the separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_ROW**            generatedcuts,      /**< cuts generated in the current separation round */
   SCIP_Real*            bestdualvector,     /**< best Lagrangian multipliers vector */
   int                   bestdualvectorlen,  /**< length of the best Lagrangian multipliers vector */
   SCIP_ROW**            aggrcuts,           /**< aggregated cuts generated so far in the current separation round */
   int*                  naggrcuts,          /**< number of aggregated cuts generated so far in the current separation round */
   SCIP_Bool*            cutoff              /**< should the current node be cutoff? */
   )
{
   SCIP_Real* cutvals;           /**< cut cofficients */
   SCIP_COL** cutcols;
   SCIP_COL** cols;
   SCIP_VAR* var;
   SCIP_ROW* cut;
   SCIP_ROW* aggrcut;
   SCIP_Real minactivity;
   SCIP_Real maxactivity;
   SCIP_Real cutlhs;
   SCIP_Real cutrhs;
   SCIP_Real cutconst;
   SCIP_Real aggrcutlhs;
   SCIP_Real aggrcutrhs;
   SCIP_Real aggrcutconst;
   SCIP_Real* aggrcutvals;
   SCIP_Real* aggrcutcoefs;
   SCIP_Real multiplier;
   SCIP_Bool aggrcutislocal;
   SCIP_Bool aggrindicator;
   SCIP_Real* tmpcutvals;
   SCIP_Real QUAD(quadterm);
   SCIP_Real QUAD(tmpcutconst);
   SCIP_Real QUAD(tmpcutrhs);
   SCIP_Real QUAD(quadprod);
   char aggrcutname[SCIP_MAXSTRLEN];
   int cutnnz;             /**< number of nonzeros in cut */
   int aggrcutnnz;
   int* aggrcutinds;            /**< column indices in cut */
   int aggrcutrank;            /**< rank of cut */
   int cutrank;
   int ncols;
   int collppos;
   int nlocalcuts;

   assert(scip != NULL);
   assert(naggrcuts != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;
   aggrcutlhs = -SCIPinfinity(scip);
   aggrcutnnz = 0;
   aggrcutrank = -1;
   nlocalcuts = 0;
   aggrindicator = FALSE;

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* allocate memory */
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &tmpcutvals, QUAD_ARRAY_SIZE(ncols)) );
   QUAD_ASSIGN(tmpcutconst, 0.0);
   QUAD_ASSIGN(tmpcutrhs, 0.0);
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &aggrcutvals, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrcutcoefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggrcutinds, ncols) );

   /* aggregate cuts based on the input Lagrangian multipliers */
   for( int i = 0; i < bestdualvectorlen; i++ )
   {
      multiplier = bestdualvector[i];
      if( SCIPisGE(scip, multiplier, 1e-4) )
      {
         cut = generatedcuts[i];
         cutnnz = SCIProwGetNNonz(cut);
         cutlhs = SCIProwGetLhs(cut);
         cutrhs = SCIProwGetRhs(cut);
         assert(SCIPisInfinity(scip, -cutlhs));
         assert(!SCIPisInfinity(scip, cutrhs));
         cutvals = SCIProwGetVals(cut);
         cutcols = SCIProwGetCols(cut);
         cutconst = SCIProwGetConstant(cut);

         for( int j = 0; j < cutnnz; j++ )
         {
            collppos = SCIPcolGetLPPos(cutcols[j]);
            assert(collppos >= 0);
            assert(collppos <= ncols);

            QUAD_ARRAY_LOAD(quadterm, tmpcutvals, collppos);
            SCIPquadprecProdDD(quadprod, multiplier, cutvals[j]);
            SCIPquadprecSumQQ(quadterm, quadterm, quadprod);
            QUAD_ARRAY_STORE(tmpcutvals, collppos, quadterm);
         }

         SCIPquadprecProdDD(quadprod, multiplier, cutconst);
         SCIPquadprecSumQQ(tmpcutconst, tmpcutconst, quadprod);
         SCIPquadprecProdDD(quadprod, multiplier, cutrhs);
         SCIPquadprecSumQQ(tmpcutrhs, tmpcutrhs, quadprod);

         cutrank = SCIProwGetRank(cut);
         aggrcutrank = MAX(aggrcutrank, cutrank);
         nlocalcuts += (SCIProwIsLocal(cut) ? 1 : 0);
         aggrindicator = TRUE;
      }
   }
   /* if at least one cut is local, then the aggregated cut is local */
   aggrcutislocal = (nlocalcuts > 0 ? TRUE : FALSE);

   /* if the aggregation was successful, then create an empty row and build a cut */
   if( aggrindicator )
   {
      aggrcutconst = QUAD_TO_DBL(tmpcutconst);
      aggrcutrhs = QUAD_TO_DBL(tmpcutrhs);

      /* build sparse representation of the aggregated cut */
      for( int i = 0; i < ncols; i++ )
      {
         QUAD_ARRAY_LOAD(quadterm, tmpcutvals, i);
         aggrcutvals[i] = QUAD_TO_DBL(quadterm);
         if( !SCIPisZero(scip, aggrcutvals[i]) )
         {
            aggrcutcoefs[aggrcutnnz] = aggrcutvals[i];
            aggrcutinds[aggrcutnnz] = i;
            aggrcutnnz++;
         }
      }

      if( aggrcutnnz == 0 && SCIPisFeasNegative(scip, aggrcutrhs) ) /*lint !e644*/
      {
         SCIPdebugMsg(scip, " -> Lagromory cut detected node infeasibility with cut 0 <= %g.\n", aggrcutrhs);
         *cutoff = TRUE;
         goto TERMINATE;
      }

      /* construct aggregated cut name */
      (void) SCIPsnprintf(aggrcutname, SCIP_MAXSTRLEN, "%s_%" SCIP_LONGINT_FORMAT "_aggregated", SCIPsepaGetName(sepa),
            sepadata->ncalls);

      /* create empty cut */
      SCIP_CALL( SCIPcreateEmptyRowSepa(scip, &aggrcut, sepa, aggrcutname, aggrcutlhs - aggrcutconst, aggrcutrhs -
               aggrcutconst, aggrcutislocal, FALSE, sepadata->dynamiccuts) );

      /* set cut rank */
      SCIProwChgRank(aggrcut, aggrcutrank); /*lint !e644*/

      /* cache the row extension and only flush them if the cut gets added */
      SCIP_CALL( SCIPcacheRowExtensions(scip, aggrcut) );

      /* collect all non-zero coefficients */
      for( int i = 0; i < aggrcutnnz; i++ )
      {
         var = SCIPcolGetVar(cols[aggrcutinds[i]]);
         SCIP_CALL( SCIPaddVarToRow(scip, aggrcut, var, aggrcutcoefs[i]) );
      }

      /* flush all changes before adding the cut */
      SCIP_CALL( SCIPflushRowExtensions(scip, aggrcut) );

      if( SCIProwGetNNonz(aggrcut) == 0 )
      {
         assert( SCIPisFeasNegative(scip, aggrcutrhs) );
         SCIPdebugMsg(scip, " -> Lagromory cut detected node infeasibility with cut 0 <= %g.\n", aggrcutrhs);
         *cutoff = TRUE;
         goto TERMINATE;
      }
      else
      {
         /* gathering lhs and rhs again in case the separator is extended later to add other cuts/constraints that may
          * have non-inf lhs or inf rhs */
         cutlhs = SCIProwGetLhs(aggrcut);
         cutrhs = SCIProwGetRhs(aggrcut);
         assert(SCIPisInfinity(scip, -cutlhs));
         assert(!SCIPisInfinity(scip, cutrhs));

         SCIPdebugMsg(scip, " -> %s cut <%s>: rhs=%f\n", "lagromory", aggrcutname, aggrcutrhs);

         /* check if the cut leads to infeasibility (i.e., *cutoff = TRUE) */
         {
            /* modifiable cuts cannot be declared infeasible, since we don't know all coefficients */
            if( SCIProwIsModifiable(aggrcut) )
               *cutoff = FALSE;

            /* check for activity infeasibility */
            minactivity = SCIPgetRowMinActivity(scip, aggrcut);
            maxactivity = SCIPgetRowMaxActivity(scip, aggrcut);

            if( (!SCIPisInfinity(scip, cutrhs) && SCIPisFeasPositive(scip, minactivity - cutrhs)) ||
                  (!SCIPisInfinity(scip, -cutlhs) && SCIPisFeasNegative(scip, maxactivity - cutlhs)) )
            {
               SCIPdebugMsg(scip, "cut <%s> is infeasible (sides=[%g,%g], act=[%g,%g])\n",
                     SCIProwGetName(aggrcut), cutlhs, cutrhs, minactivity, maxactivity);
               *cutoff = TRUE;
            }
         }

         /* add the aggregated cut to a separate data structure */
         aggrcuts[*naggrcuts] = aggrcut;
         (*naggrcuts)++;
      }

      QUAD_ASSIGN(quadterm, 0.0);
      for( int i = 0; i < ncols; i++ )
      {
         aggrcutvals[i] = 0.0;
         QUAD_ARRAY_STORE(tmpcutvals, i, quadterm);
      }
   }

TERMINATE:
   /* free memory */
   SCIPfreeBufferArray(scip, &aggrcutinds);
   SCIPfreeBufferArray(scip, &aggrcutcoefs);
   SCIPfreeCleanBufferArray(scip, &aggrcutvals);
   SCIPfreeCleanBufferArray(scip, &tmpcutvals);

   return SCIP_OKAY;
}

/** main method: LP solution separation method of separator */
static
SCIP_RETCODE generateGMICuts(
   SCIP*                 scip,                     /**< SCIP data structure */
   SCIP_SEPA*            sepa,                     /**< pointer to the separator */
   SCIP_SEPADATA*        sepadata,                 /**< separator data structure */
   int                   mainiternum,              /**< iteration number of the outer loop of the relax-and-cut algorithm */
   int                   subgradientiternum,       /**< iteration number of the subgradient algorithm */
   SCIP_SOL*             sol,                      /**< LP solution to be used for cut generation */
   SCIP_Real*            solvals,                  /**< values of the LP solution to be used for cut generation */
   int                   nmaxgeneratedperroundcuts,/**< maximal number of cuts allowed to generate per separation round */
   SCIP_Bool             allowlocal,               /**< should locally valid cuts be generated? */
   SCIP_ROW**            generatedcurrroundcuts,   /**< cuts generated in the current separation round */
   SCIP_Real*            generatedcutefficacies,   /**< array of generated cut efficacies w.r.t. the respective LP bases used for cut
                                                     generations */
   int                   ngeneratedcurrroundcuts,  /**< number of cuts generated until the previous basis in the current separation round */
   int*                  ngeneratednewcuts,        /**< number of new cuts generated using the current basis */
   int                   depth,                    /**< depth of the current node in the tree */
   SCIP_Bool*            cutoff                    /**< should the current node be cutoff? */
   )
{
   SCIP_Real minfrac;
   SCIP_Real maxfrac;
   SCIP_Real frac;
   SCIP_Real* basisfrac;
   SCIP_Real* binvrow;
   SCIP_AGGRROW* aggrrow;
   SCIP_Bool success;
   SCIP_Real* cutcoefs;
   SCIP_Real cutrhs;
   SCIP_Real cutefficacy;
   SCIP_Bool cutislocal;
   SCIP_ROW** rows;
   SCIP_COL** cols;
   SCIP_VAR* var;
   SCIP_ROW* row;
   int cutrank;
   int cutnnz;
   int* cutinds;
   int* basisind;
   int* inds;
   int nrows;
   int ncols;
   int* basisperm;
   int k;
   int c;
   int ninds;
   int nmaxcutsperlp;

   assert(ngeneratednewcuts != NULL);

   minfrac = sepadata->away;
   maxfrac = 1.0 - sepadata->away;
   SCIP_CALL( SCIPgetLPRowsData(scip, &rows, &nrows) );
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   *ngeneratednewcuts = 0;
   nmaxcutsperlp = ((depth == 0) ? sepadata->nmaxcutsperlproot : sepadata->nmaxcutsperlp);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &basisperm, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisfrac, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &binvrow, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inds, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &basisind, nrows) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutcoefs, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutinds, ncols) );
   SCIP_CALL( SCIPaggrRowCreate(scip, &aggrrow) );

   SCIP_CALL( SCIPgetLPBasisInd(scip, basisind) );

   /* check if the rows of the simplex tableau are suitable for cut generation and build an array of fractions */
   for( int i = 0; i < nrows; ++i )
   {
      frac = 0.0;

      c = basisind[i];

      basisperm[i] = i;

      /* if the simplex tableau row corresponds to an LP column */
      if( c >= 0 )
      {
         assert(c < ncols);

         var = SCIPcolGetVar(cols[c]);
         /* if the column is non-continuous one */
         if( SCIPvarGetType(var) != SCIP_VARTYPE_CONTINUOUS )
         {
            frac = SCIPfeasFrac(scip, solvals[c]);
            frac = MIN(frac, 1.0 - frac);
         }
      }
      /* if the simplex tableau row corresponds to an LP row and separation on rows is allowed */
      else if( sepadata->separaterows )
      {
         assert(0 <= -c-1);
         assert(-c-1 < nrows);

         row = rows[-c-1];
         /* if the row is suitable for cut generation */
         if( SCIProwIsIntegral(row) && !SCIProwIsModifiable(row) )
         {
            frac = SCIPfeasFrac(scip, SCIPgetRowActivity(scip, row));
            frac = MIN(frac, 1.0 - frac);
         }
      }

      if( frac >= minfrac )
      {
         /* slightly change fractionality to have random order for equal fractions */
         basisfrac[i] = frac + SCIPrandomGetReal(sepadata->randnumgen, -1e-6, 1e-6);
      }
      else
      {
         basisfrac[i] = 0.0;
      }
   }

   /* if there is a need to sort the fractionalities */
   if( sepadata->sortcutoffsol )
   {
      /* sort basis indices by fractionality */
      SCIPsortDownRealInt(basisfrac, basisperm, nrows);
   }

   /* for all basic columns belonging to integer variables, try to generate a GMI cut */
   for( int i = 0; i < nrows && !SCIPisStopped(scip) && !*cutoff; ++i )
   {
      if( (ngeneratedcurrroundcuts + *ngeneratednewcuts >= nmaxgeneratedperroundcuts) ||
            (sepadata->ntotalcuts + *ngeneratednewcuts >= sepadata->nmaxtotalcuts) ||
            (*ngeneratednewcuts >= nmaxcutsperlp) )
         break;

      ninds = -1;
      cutefficacy = 0.0;

      /* either break the loop or proceed to the next iteration if the fractionality is zero */
      if( basisfrac[i] == 0.0 )
      {
         if( sepadata->sortcutoffsol )
            break;
         else
            continue;
      }

      k = basisperm[i];

      /* get the row of B^-1 for this basic integer variable with fractional solution value and call aggregate function */
      SCIP_CALL( SCIPgetLPBInvRow(scip, k, binvrow, inds, &ninds) );

      SCIP_CALL( SCIPaggrRowSumRows(scip, aggrrow, binvrow, inds, ninds, sepadata->sidetypebasis, allowlocal, 2,
               (int) MAXAGGRLEN(ncols), &success) );

      if( !success )
         continue;

      /* @todo Currently we are using the SCIPcalcMIR() function to compute the coefficients of the Gomory
       *       cut. Alternatively, we could use the direct version (see thesis of Achterberg formula (8.4)) which
       *       leads to cut a of the form \sum a_i x_i \geq 1. Rumor has it that these cuts are better.
       */

      /* try to create GMI cut out of the aggregation row */
      SCIP_CALL( SCIPcalcMIR(scip, sol, POSTPROCESS, BOUNDSWITCH, USEVBDS, allowlocal, FIXINTEGRALRHS, NULL,
               NULL, minfrac, maxfrac, 1.0, aggrrow, cutcoefs, &cutrhs, cutinds, &cutnnz, &cutefficacy, &cutrank,
               &cutislocal, &success) );

      if( success )
      {
         assert(allowlocal || !cutislocal); /*lint !e644*/

         SCIP_CALL( constructCutRow(scip, sepa, sepadata, mainiternum, subgradientiternum, cutnnz, cutinds, cutcoefs,
                  cutefficacy, cutrhs, cutislocal, cutrank, generatedcurrroundcuts, generatedcutefficacies,
                  ngeneratedcurrroundcuts, ngeneratednewcuts, cutoff));
      }
   }

   /* free temporary memory */
   SCIPfreeBufferArray(scip, &cutinds);
   SCIPfreeBufferArray(scip, &cutcoefs);
   SCIPfreeBufferArray(scip, &basisind);
   SCIPfreeBufferArray(scip, &inds);
   SCIPfreeBufferArray(scip, &binvrow);
   SCIPfreeBufferArray(scip, &basisfrac);
   SCIPfreeBufferArray(scip, &basisperm);
   SCIPaggrRowFree(scip, &aggrrow);

   return SCIP_OKAY;
}

/** update objective vector w.r.t. the fixed Lagrangian multipliers */
static
SCIP_RETCODE updateObjectiveVector(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            dualvector,         /**< Lagrangian multipliers vector */
   SCIP_ROW**            cuts,               /**< cuts generated so far in the current separation round */
   int                   ncuts,              /**< number of cuts generated so far in the current separation round */
   SCIP_Real*            origobjcoefs,       /**< original objective function coefficients of the node linear relaxation */
   SCIP_Bool*            objvecsdiffer       /**< whether the updated objective function coefficients differ from the old ones */
   )
{
   SCIP_Real* objvals;
   SCIP_Real* prod;
   SCIP_Real* oldobjvals;
   SCIP_Real* cutvals;
   SCIP_COL** cutcols;
   SCIP_COL** cols;
   SCIP_VAR* var;
   int cutnnz;
   int collppos;
   int ncols;

   assert(objvecsdiffer != NULL);
   assert(ncuts > 0);

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &objvals, ncols) );
   SCIP_CALL( SCIPallocBufferArray(scip, &oldobjvals, ncols) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &prod, ncols) );
   *objvecsdiffer = FALSE;

   /* find the product of Lagrangian multipliers and cut coefficients */
   for( int i = 0; i < ncuts; i++ )
   {
      if( !SCIPisZero(scip, dualvector[i]) )
      {
         assert(!(SCIPisInfinity(scip, -SCIProwGetLhs(cuts[i])) && SCIPisInfinity(scip, SCIProwGetRhs(cuts[i]))));
         assert(SCIPisInfinity(scip, -SCIProwGetLhs(cuts[i])));
         assert(!SCIPisInfinity(scip, SCIProwGetRhs(cuts[i])));

         cutnnz = SCIProwGetNNonz(cuts[i]);
         assert(cutnnz <= ncols);
         cutvals = SCIProwGetVals(cuts[i]);
         cutcols = SCIProwGetCols(cuts[i]);

         for( int j = 0; j < cutnnz; ++j )
         {
            collppos = SCIPcolGetLPPos(cutcols[j]);
            assert(collppos >= 0);
            assert(collppos <= ncols);

            prod[collppos] += dualvector[i] * cutvals[j];
         }
      }
   }

   /* change objective coefficients */
   for( int i = 0; i < ncols; i++ )
   {
      var = SCIPcolGetVar(cols[i]);
      oldobjvals[i] = SCIPgetVarObjDive(scip, var);
      objvals[i] = origobjcoefs[i] + prod[i];

      SCIP_CALL( SCIPchgVarObjDive(scip, var, objvals[i]) );

      /* identify if the updated objective vector is indeed different from the previous one */
      if( !(*objvecsdiffer) && !SCIPisEQ(scip, oldobjvals[i], objvals[i]) )
         *objvecsdiffer = TRUE;
   }

   for( int i = 0; i < ncols; i++)
   {
      prod[i] = 0.0;
   }

   /* free memory  */
   SCIPfreeCleanBufferArray(scip, &prod);
   SCIPfreeBufferArray(scip, &oldobjvals);
   SCIPfreeBufferArray(scip, &objvals);

   return SCIP_OKAY;
}

/** add GMI cuts to the objective function of the Lagrangian dual problem by introducing new Lagrangian multipliers */
static
SCIP_RETCODE addGMICutsAsSoftConss(
   SCIP_Real*            dualvector,         /**< Lagrangian multipliers vector */
   int                   ngeneratedcuts,     /**< number of cuts generated so far in the current separation round */
   int*                  naddedcuts,         /**< number of cuts added so far in the current separation round to the Lagrangian dual problem
                                               upon penalization */
   int*                  nnewaddedsoftcuts   /**< number of cuts added newly to the Lagrangian dual problem upon penalization */
   )
{
   assert(*naddedcuts <= ngeneratedcuts);

   /* set the initial penalty of the newly penalized cuts as zero */
   for( int i = *naddedcuts; i < ngeneratedcuts; i++ )
      dualvector[i] = 0.0;

   *nnewaddedsoftcuts = ngeneratedcuts - *naddedcuts;
   *naddedcuts = ngeneratedcuts;

   return SCIP_OKAY;
}

/** solve the Lagrangian dual problem */
static
SCIP_RETCODE solveLagrangianDual(
   SCIP*                 scip,                     /**< SCIP data structure */
   SCIP_SEPA*            sepa,                     /**< pointer to the separator */
   SCIP_SEPADATA*        sepadata,                 /**< separator data structure */
   SCIP_SOL*             sol,                      /**< data structure to store an LP solution upon solving a Lagrangian dual problem with
                                                     fixed Lagrangian multipliers */
   SCIP_Real*            solvals,                  /**< values of the LP solution obtained upon solving a Lagrangian dual problem with fixed
                                                     Lagrangian multipliers */
   int                   mainiternum,              /**< iteration number of the outer loop of the relax-and-cut algorithm */
   SCIP_Real             ubparam,                  /**< estimate of the optimal Lagrangian dual value */
   int                   depth,                    /**< depth of the current node in the tree */
   SCIP_Bool             allowlocal,               /**< should locally valid cuts be generated? */
   int                   nmaxgeneratedperroundcuts,/**< maximal number of cuts allowed to generate per separation round */
   SCIP_Real*            origobjcoefs,             /**< original objective function coefficients of the node linear relaxation */
   SCIP_Real             origobjoffset,            /**< original objective function offset of the node linear relaxation */
   SCIP_Real*            dualvector,               /**< Lagrangian multipliers vector */
   int*                  nsoftcuts,                /**< number of generated cuts that were penalized and added to the Lagrangian dual problem */
   SCIP_ROW**            generatedcurrroundcuts,   /**< cuts generated in the current separation round */
   SCIP_Real*            generatedcutefficacies,   /**< array of generated cut efficacies w.r.t. the respective LP bases used for cut
                                                     generations */
   int*                  ngeneratedcutsperiter,    /**< number of cuts generated per subgradient iteration in the current separation round */
   int*                  ngeneratedcurrroundcuts,  /**< number of cuts generated so far in the current separation round */
   int*                  ncurrroundlpiters,        /**< number of LP iterations taken for solving Lagrangian dual problems with fixed
                                                     multipliers in the current separator round */
   SCIP_Bool*            cutoff,                   /**< should the current node be cutoff? */
   SCIP_Real*            bestlagrangianval,        /**< best Lagrangian value found so far */
   SCIP_Real*            bestdualvector,           /**< Lagrangian multipliers corresponding to the best Lagrangian value found so far */
   int*                  bestdualvectorlen,        /**< length of the Lagrangian multipliers corresponding to the best Lagrangian value
                                                     found so far */
   int*                  nbestdualupdates,         /**< number of best Lagrangian values found so far */
   int*                  totaliternum              /**< total number of iterations of the relax-and-cut algorithm performed so far */
   )
{
   SCIP_Real* subgradient;
   SCIP_Real muparam;
   SCIP_Real steplength;
   SCIP_Real objval;
   SCIP_Real lagrangianval;
   SCIP_Real* lagrangianvals;
   SCIP_Real avglagrangianval;
   SCIP_Real maxsoftcutviol;
   SCIP_Real maxnzsubgradientdualprod;
   SCIP_Real maxviolscore;
   SCIP_Real maxviolscoreold;
   SCIP_Real nviolscore;
   SCIP_Real nviolscoreold;
   SCIP_Real scoreweight;
   SCIP_Real ballradius;
   SCIP_Bool solfound;
   SCIP_Bool backtrack;
   SCIP_Bool terminate;
   SCIP_Bool subgradientzero;
   SCIP_Bool objvecsdiffer;
   SCIP_Bool dualvecsdiffer;
   SCIP_Bool solvelp;
   int ncurrroundlpiterslast;
   int nlpiters;
   int ngeneratednewcuts;
   int nnewaddedsoftcuts;
   int nsoftcutviols;
   int nnzsubgradientdualprod;

   SCIP_CALL( SCIPallocBufferArray(scip, &lagrangianvals, sepadata->nmaxsubgradientiters) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &subgradient, nmaxgeneratedperroundcuts) );

   muparam = sepadata->muparaminit;
   steplength = 0.0;
   objval = 0.0;
   avglagrangianval = 0.0;
   maxsoftcutviol = 0.0;
   maxnzsubgradientdualprod = 0.0;
   maxviolscore = 0.0;
   nviolscore = 0.0;
   scoreweight = 1.0;
   ballradius = sepadata->radiusinit;
   ngeneratednewcuts = 0;
   nsoftcutviols = 0;
   nnzsubgradientdualprod = 0;
   terminate = FALSE;
   subgradientzero = FALSE;
   objvecsdiffer = FALSE;
   dualvecsdiffer = FALSE;
   solvelp = TRUE;

   /* update objective vector based on input Lagrangian multipliers */
   if( *nsoftcuts > 0 )
   {
      SCIP_CALL( updateObjectiveVector(scip, dualvector, generatedcurrroundcuts, *nsoftcuts, origobjcoefs, &objvecsdiffer) );
   }

   /* termination check */
   SCIP_CALL( checkLagrangianDualTermination(sepadata, -1, -1, FALSE, *ngeneratedcurrroundcuts,
            nmaxgeneratedperroundcuts, *ncurrroundlpiters, depth, &terminate) );

   /* the subgradient algorithm loop */
   for( int i = 0; i < sepadata->nmaxsubgradientiters && !SCIPisStopped(scip) && !terminate; i++ )
   {
      solfound = FALSE;
      subgradientzero = FALSE;
      objvecsdiffer = FALSE;
      dualvecsdiffer = FALSE;
      nnewaddedsoftcuts = 0;
      scoreweight *= sepadata->radiusupdateweight;

      ncurrroundlpiterslast = *ncurrroundlpiters;
      if( solvelp )
      {
         /* solve Lagrangian dual for fixed Lagrangian multipliers */
         SCIP_CALL( solveLagromoryLP(scip, sepadata, depth, origobjoffset, &solfound, sol, solvals, &objval,
                  ncurrroundlpiters) );
      }
      nlpiters = *ncurrroundlpiters - ncurrroundlpiterslast;

      /* if an optimal solution has been found, then generate cuts and do other operations */
      if( solfound )
      {
         /* generate GMI cuts if a new basis solution is found */
         if( (nlpiters >= 1) && (i % sepadata->cutgenfreq == 0) )
         {
            ngeneratednewcuts = 0;
            SCIP_CALL( generateGMICuts(scip, sepa, sepadata, mainiternum, i, sol, solvals,
                     nmaxgeneratedperroundcuts, allowlocal, generatedcurrroundcuts, generatedcutefficacies,
                     *ngeneratedcurrroundcuts, &ngeneratednewcuts, depth, cutoff));
            sepadata->ntotalcuts += ngeneratednewcuts;
            *ngeneratedcurrroundcuts += ngeneratednewcuts;
            ngeneratedcutsperiter[mainiternum * sepadata->nmaxsubgradientiters + i + 1] = ngeneratednewcuts;
         }

         /* update subgradient, i.e., find the residuals of the penalized cuts, and determine various violations */
         updateSubgradient(scip, sol, generatedcurrroundcuts, *nsoftcuts, subgradient, dualvector, &subgradientzero,
               &nsoftcutviols, &maxsoftcutviol, &nnzsubgradientdualprod, &maxnzsubgradientdualprod);

         /* calculate Lagrangian value for the fixed Lagrangian multipliers, and update best and avg values */
         updateLagrangianValue(scip, objval, dualvector, generatedcurrroundcuts, *nsoftcuts, &lagrangianval);
         if( SCIPisPositive(scip, lagrangianval - *bestlagrangianval) )
         {
            *bestlagrangianval = lagrangianval;
            for( int j = 0; j < *nsoftcuts; j++ )
            {
               bestdualvector[j] = dualvector[j];
            }
            *bestdualvectorlen = *nsoftcuts;
            (*nbestdualupdates)++;
         }
         lagrangianvals[i] = lagrangianval;
         if( i < sepadata->nmaxlagrangianvalsforavg )
         {
            avglagrangianval = (avglagrangianval * i + lagrangianval)/(i+1);
         }
         else
         {
            avglagrangianval = (avglagrangianval * sepadata->nmaxlagrangianvalsforavg -
                  lagrangianvals[i - sepadata->nmaxlagrangianvalsforavg] +
                  lagrangianval)/(sepadata->nmaxlagrangianvalsforavg);
         }

         /* if the subgradient vector is non-zero, then update the mu parameter and the Lagrangian multipliers */
         if( !subgradientzero )
         {
            /* update mu param */
            SCIP_CALL( updateMuSteplengthParam(scip, sepadata, i, ubparam, lagrangianvals, *bestlagrangianval, avglagrangianval,
                     &muparam, &backtrack) );

            /* update step length */
            updateStepLength(scip, muparam, ubparam, lagrangianval, subgradient, *nsoftcuts, &steplength);

            /* update scores to determine how to update the stabilization ball radius */
            maxviolscoreold = maxviolscore;
            nviolscoreold = nviolscore;
            maxviolscore = (1.0 - scoreweight) * maxsoftcutviol + scoreweight * maxnzsubgradientdualprod;
            nviolscore = (1.0 - scoreweight) * nsoftcutviols + scoreweight * nnzsubgradientdualprod;

            /* update Lagrangian multipliers */
            SCIP_CALL( updateDualVector(scip, sepadata, dualvector, bestdualvector, *bestdualvectorlen,
                     *nbestdualupdates, i, *totaliternum, steplength, subgradient, *nsoftcuts, backtrack, maxviolscore,
                     maxviolscoreold, nviolscore, nviolscoreold, nlpiters, &dualvecsdiffer, &ballradius) );

            /* update objective vector based on updated Lagrangian multipliers */
            if( dualvecsdiffer )
            {
               SCIP_CALL( updateObjectiveVector(scip, dualvector, generatedcurrroundcuts, *nsoftcuts, origobjcoefs, &objvecsdiffer) );
            }
         }
         /* if the subgradient vector if zero, then simply mark that the Lagrangian multipliers and the objective
          * function of the Lagrangian dual did not change */
         else
         {
            dualvecsdiffer = FALSE;
            objvecsdiffer = FALSE;
         }


         /* add generated GMI cuts to the objective function of the Lagrangian dual problem by introducing new
          * Lagrangian multipliers */
         if( (i % sepadata->cutaddfreq == 0) || (!dualvecsdiffer && !objvecsdiffer &&
                  (*ngeneratedcurrroundcuts - *nsoftcuts > 0)) )
         {
            SCIP_CALL( addGMICutsAsSoftConss(dualvector, *ngeneratedcurrroundcuts, nsoftcuts, &nnewaddedsoftcuts) );
         }
      }
      else
      {
         /* add any remaining generated GMI cuts to the objective function of the Lagrangian dual problem by introducing
          * new Lagrangian multipliers */
         if( (*ngeneratedcurrroundcuts - *nsoftcuts) > 0 )
         {
            SCIP_CALL( addGMICutsAsSoftConss(dualvector, *ngeneratedcurrroundcuts, nsoftcuts, &nnewaddedsoftcuts) );
         }

         solvelp = FALSE;
      }

      /* termination check */
      SCIP_CALL( checkLagrangianDualTermination(sepadata, nnewaddedsoftcuts, *ngeneratedcurrroundcuts - *nsoftcuts,
               objvecsdiffer, *ngeneratedcurrroundcuts, nmaxgeneratedperroundcuts, *ncurrroundlpiters, depth,
               &terminate) );

      (*totaliternum)++;
   }

   /* add any remaining generated GMI cuts to the objective function of the Lagrangian dual problem by introducing new
    * Lagrangian multipliers */
   if( (*ngeneratedcurrroundcuts - *nsoftcuts) > 0 )
   {
      SCIP_CALL( addGMICutsAsSoftConss(dualvector, *ngeneratedcurrroundcuts, nsoftcuts, &nnewaddedsoftcuts) );
   }

   /* free memory */
   for( int i = 0; i < nmaxgeneratedperroundcuts; i++)
   {
      subgradient[i] = 0.0;
   }
   SCIPfreeCleanBufferArray(scip, &subgradient);
   SCIPfreeBufferArray(scip, &lagrangianvals);

   return SCIP_OKAY;
}

/** generates initial cut pool before solving the Lagrangian dual */
static
SCIP_RETCODE generateInitCutPool(
   SCIP*                 scip,                     /**< SCIP data structure */
   SCIP_SEPA*            sepa,                     /**< separator */
   SCIP_SEPADATA*        sepadata,                 /**< separator data structure */
   int                   mainiternum,              /**< iteration number of the outer loop of the relax-and-cut algorithm */
   SCIP_SOL*             sol,                      /**< LP solution to be used for cut generation */
   SCIP_Real*            solvals,                  /**< values of the LP solution to be used for cut generation */
   int                   nmaxgeneratedperroundcuts,/**< maximal number of cuts allowed to generate per separation round */
   SCIP_Bool             allowlocal,               /**< should locally valid cuts be generated? */
   SCIP_ROW**            generatedcurrroundcuts,   /**< cuts generated in the current separation round */
   SCIP_Real*            generatedcutefficacies,   /**< array of generated cut efficacies w.r.t. the respective LP bases used for cut
                                                     generations */
   int*                  ngeneratedcutsperiter,    /**< number of cuts generated per subgradient iteration in the current separation round */
   int*                  ngeneratedcurrroundcuts,  /**< number of cuts generated so far in the current separation round */
   int                   depth,                    /**< depth of the current node in the tree */
   SCIP_Bool*            cutoff                    /**< should the current node be cutoff? */
   )
{
   int ngeneratednewcuts;

   ngeneratednewcuts = 0;

   /* generate initial set of cuts */
   SCIP_CALL( generateGMICuts(scip, sepa, sepadata, mainiternum, -1, sol, solvals, nmaxgeneratedperroundcuts,
            allowlocal, generatedcurrroundcuts, generatedcutefficacies, *ngeneratedcurrroundcuts, &ngeneratednewcuts,
            depth, cutoff) );

   /* update certain statistics */
   sepadata->ntotalcuts += ngeneratednewcuts;
   *ngeneratedcurrroundcuts += ngeneratednewcuts;
   ngeneratedcutsperiter[sepadata->nmaxsubgradientiters * mainiternum] = ngeneratednewcuts;

   return SCIP_OKAY;
}

/** add cuts to SCIP */
static
SCIP_RETCODE addCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_ROW**            cuts,               /**< cuts generated so far in the current separation round */
   int                   ncuts,              /**< number of cuts generated so far in the current separation round */
   SCIP_Longint          maxdnom,            /**< maximum denominator in the rational representation of cuts */
   SCIP_Real             maxscale,           /**< maximal scale factor to scale the cuts to integral values */
   int*                  naddedcuts,         /**< number of cuts added to either global cutpool or sepastore */
   SCIP_Bool*            cutoff              /**< should the current node be cutoff? */
   )
{
   SCIP_ROW* cut;
   SCIP_Bool madeintegral;

   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(naddedcuts != NULL);
   assert(cutoff != NULL);

   *cutoff = FALSE;
   madeintegral = FALSE;

   for( int i = 0; i < ncuts && !*cutoff; i++ )
   {
      cut = cuts[i];

      if( SCIProwGetNNonz(cut) == 1 )
      {
         /* Add the bound change as cut to avoid that the LP gets modified. This would mean that the LP is not flushed
          * and the method SCIPgetLPBInvRow() fails; SCIP internally will apply this bound change automatically. */
         SCIP_CALL( SCIPaddRow(scip, cut, TRUE, cutoff) );
         ++(*naddedcuts);
      }
      else
      {
         if( sepadata->makeintegral && SCIPgetRowNumIntCols(scip, cut) == SCIProwGetNNonz(cut) )
         {
            /* try to scale the cut to integral values */
            SCIP_CALL( SCIPmakeRowIntegral(scip, cut, -SCIPepsilon(scip), SCIPsumepsilon(scip),
                     maxdnom, maxscale, MAKECONTINTEGRAL, &madeintegral) );

            /* if RHS = plus infinity (due to scaling), the cut is useless, so we are not adding it */
            if( madeintegral && SCIPisInfinity(scip, SCIProwGetRhs(cut)) )
               return SCIP_OKAY;
         }

         if( SCIPisCutNew(scip, cut) )
         {
            /* add global cuts which are not implicit bound changes to the cut pool */
            if( !SCIProwIsLocal(cut) )
            {
               if( sepadata->delayedcuts )
               {
                  SCIP_CALL( SCIPaddDelayedPoolCut(scip, cut) );
               }
               else
               {
                  SCIP_CALL( SCIPaddPoolCut(scip, cut) );
               }
            }
            else
            {
               /* local cuts we add to the sepastore */
               SCIP_CALL( SCIPaddRow(scip, cut, sepadata->forcecuts, cutoff) );
            }
            ++(*naddedcuts);
         }
      }
   }

   return SCIP_OKAY;
}

/** check different termination criteria */
/* @todo nlpssolved criterion? */
static
SCIP_RETCODE checkMainLoopTermination(
   SCIP_SEPADATA*        sepadata,                 /**< separator data structure */
   SCIP_Bool             cutoff,                   /**< should the current node be cutoff? */
   SCIP_Bool             dualvecsdiffer,           /**< whether the updated Lagrangian multipliers differ from the old one */
   int                   ngeneratedcurrroundcuts,  /**< number of cuts generated in the current separation round */
   int                   nsoftcuts,                /**< number of generated cuts that were penalized and added to the Lagrangian dual problem */
   int                   nmaxgeneratedperroundcuts,/**< maximal number of cuts allowed to generate per separation round */
   int                   ncurrroundlpiters,        /**< number of LP iterations taken for solving Lagrangian dual problems with fixed
                                                     multipliers in the current separator round */
   int                   depth,                    /**< depth of the current node in the tree */
   SCIP_Bool*            terminate                 /**< whether to terminate the relax-and-cut algorithm */
   )
{
   *terminate = FALSE;

   /* check if the node has been identified to be cutoff */
   if( cutoff )
      *terminate = TRUE;

   /* check if the Lagrangian multipliers do not differ from the previous iteration and no new cuts exist for penalizing */
   if( !dualvecsdiffer && (ngeneratedcurrroundcuts == nsoftcuts) )
      *terminate = TRUE;

   /* check if allowed number of cuts in this separation round have already been generated */
   if( ngeneratedcurrroundcuts >= nmaxgeneratedperroundcuts )
      *terminate = TRUE;

   /* check if allowed number of cuts in the tree have already been generated */
   if( sepadata->ntotalcuts >= sepadata->nmaxtotalcuts )
      *terminate = TRUE;

   /* check if allowed number of simplex iterations in this separation round have already been used up */
   if( (sepadata->nmaxperroundlpiters >= 0) && (ncurrroundlpiters >= sepadata->nmaxperroundlpiters) )
      *terminate = TRUE;

   /* check if allowed number of simplex iterations in the root node have already been used up */
   if( (depth == 0) && (sepadata->nmaxrootlpiters >= 0) && (sepadata->nrootlpiters >= sepadata->nmaxrootlpiters) )
      *terminate = TRUE;

   /* check if allowed number of simplex iterations in the tree have already been used up */
   if( (sepadata->nmaxtotallpiters >= 0) && (sepadata->ntotallpiters >= sepadata->nmaxtotallpiters) )
      *terminate = TRUE;

   return SCIP_OKAY;
}

/** Searches and tries to add Lagromory cuts */
static
SCIP_RETCODE separateCuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   SCIP_Real             ubparam,            /**< estimate of the optimal Lagrangian dual value */
   int                   depth,              /**< depth of the current node in the tree */
   SCIP_Bool             allowlocal,         /**< should locally valid cuts be generated? */
   SCIP_RESULT*          result              /**< final result of the separation round */
   )
{
   SCIP_ROW** generatedcurrroundcuts;
   SCIP_ROW** aggregatedcurrroundcuts;
   SCIP_Real* generatedcutefficacies;
   SCIP_Bool solfound;
   SCIP_SOL* softcutslpsol;
   SCIP_Real* softcutslpsolvals;
   SCIP_SOL* hardcutslpsol;
   SCIP_Real* hardcutslpsolvals;
   SCIP_Real* dualsol;
   SCIP_Real* dualvector;
   SCIP_Real* bestdualvector;
   SCIP_Real bestlagrangianval;
   SCIP_Real* origobjcoefs;
   SCIP_Real origobjoffset;
   SCIP_Real objval;
   SCIP_Real maxscale;
   SCIP_Longint maxdnom;
   SCIP_Bool cutoff;
   SCIP_Bool cutoff2;
   SCIP_Bool dualvecsdiffer;
   SCIP_Bool terminate;
   SCIP_COL** cols;

   int* ngeneratedcutsperiter;
   int bestdualvectorlen;
   int nbestdualupdates;
   int totaliternum;
   int* cutindsperm;
   int nprocessedcuts;
   int ncols;
   int nrows;
   int ngeneratedcurrroundcuts;
   int nselectedcurrroundcuts;
   int nselectedcuts;
   int naddedcurrroundcuts;
   int naggregatedcurrroundcuts;
   int nmaxgeneratedperroundcuts;
   int ncurrroundlpiters;
   int nsoftcuts;
   int nsoftcutsold;
   int maxdepth;

   assert(*result == SCIP_DIDNOTRUN);
   assert(sepadata != NULL);
   sepadata->ncalls = SCIPsepaGetNCalls(sepa);
   objval = SCIPinfinity(scip);
   bestlagrangianval = -SCIPinfinity(scip);
   bestdualvectorlen = 0;
   nbestdualupdates = 0;
   totaliternum = 0;
   ngeneratedcurrroundcuts = 0;
   naddedcurrroundcuts = 0;
   naggregatedcurrroundcuts = 0;
   ncurrroundlpiters = 0;
   nsoftcuts = 0;
   solfound = FALSE;
   cutoff = FALSE;
   cutoff2 = FALSE;
   dualvecsdiffer = FALSE;
   terminate = FALSE;

   SCIPdebugMsg(scip, "Separating cuts...\n");

   /* initialize the LP that will be used to solve the Lagrangian dual with fixed Lagrangian multipliers */
   SCIP_CALL( createLPWithSoftCuts(scip, sepadata) );
   /* create the LP that represents the node relaxation including all the generated cuts in this separator */
   SCIP_CALL( createLPWithHardCuts(scip, sepadata, NULL, 0) );

   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   nrows = SCIPgetNLPRows(scip);

   /* get the maximal number of cuts allowed in a separation round */
   nmaxgeneratedperroundcuts = ((depth == 0) ? sepadata->nmaxperroundcutsroot : sepadata->nmaxperroundcuts);

   /* allocate memory */
   SCIP_CALL( SCIPallocBufferArray(scip, &generatedcurrroundcuts, nmaxgeneratedperroundcuts) );
   SCIP_CALL( SCIPallocBufferArray(scip, &aggregatedcurrroundcuts, nmaxgeneratedperroundcuts) );
   SCIP_CALL( SCIPallocBufferArray(scip, &generatedcutefficacies, nmaxgeneratedperroundcuts) );
   SCIP_CALL( SCIPallocBufferArray(scip, &cutindsperm, nmaxgeneratedperroundcuts) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &ngeneratedcutsperiter, sepadata->nmaxmainiters *
            (sepadata->nmaxsubgradientiters + 1)) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &dualsol, nrows + nmaxgeneratedperroundcuts) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &dualvector, nmaxgeneratedperroundcuts) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &bestdualvector, nmaxgeneratedperroundcuts) );
   SCIP_CALL( SCIPallocBufferArray(scip, &softcutslpsolvals, ncols) );
   SCIP_CALL( SCIPcreateSol(scip, &softcutslpsol, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &hardcutslpsolvals, ncols) );
   SCIP_CALL( SCIPcreateSol(scip, &hardcutslpsol, NULL) );
   SCIP_CALL( SCIPallocBufferArray(scip, &origobjcoefs, ncols) );

   /* store current objective function */
   for( int i = 0; i < ncols; i++ )
   {
      origobjcoefs[i] = SCIPcolGetObj(cols[i]);
   }
   origobjoffset = SCIPgetTransObjoffset(scip);

   /* solve node LP relaxation to have an initial simplex tableau */
   SCIP_CALL( solveLagromoryLP(scip, sepadata, depth, origobjoffset, &solfound, softcutslpsol, softcutslpsolvals, &objval,
            &ncurrroundlpiters));

   /* generate initial cut pool */
   SCIP_CALL( generateInitCutPool(scip, sepa, sepadata, 0, softcutslpsol, softcutslpsolvals, nmaxgeneratedperroundcuts, allowlocal,
            generatedcurrroundcuts, generatedcutefficacies, ngeneratedcutsperiter, &ngeneratedcurrroundcuts, depth,
            &cutoff) );

   /* termination check */
   SCIP_CALL( checkMainLoopTermination(sepadata, cutoff, TRUE, ngeneratedcurrroundcuts, nsoftcuts,
            nmaxgeneratedperroundcuts, ncurrroundlpiters, depth, &terminate) );

   /* compute cuts for each integer col with fractional val */
   for( int i = 0; i < sepadata->nmaxmainiters && !SCIPisStopped(scip) && !terminate; ++i )
   {
      nsoftcutsold = nsoftcuts;
      nsoftcuts = ngeneratedcurrroundcuts;
      dualvecsdiffer = FALSE;

      /* solve Lagrangain dual */
      SCIP_CALL( solveLagrangianDual(scip, sepa, sepadata, softcutslpsol, softcutslpsolvals, i, ubparam, depth, allowlocal,
               nmaxgeneratedperroundcuts, origobjcoefs, origobjoffset, dualvector, &nsoftcuts, generatedcurrroundcuts,
               generatedcutefficacies, ngeneratedcutsperiter, &ngeneratedcurrroundcuts, &ncurrroundlpiters, &cutoff,
               &bestlagrangianval, bestdualvector, &bestdualvectorlen, &nbestdualupdates, &totaliternum) );

      /* @todo filter cuts before adding them to the new LP that was created based on the node relaxation? */

      /* update the LP representing the node relaxation by adding newly generated cuts */
      if( !cutoff && (ngeneratedcurrroundcuts - nsoftcutsold > 0) )
      {
         SCIP_CALL( createLPWithHardCuts(scip, sepadata, generatedcurrroundcuts, ngeneratedcurrroundcuts) );

         /* solve the LP and get relevant information */
         SCIP_CALL( solveLPWithHardCuts(scip, sepadata, &solfound, hardcutslpsol, hardcutslpsolvals) );

         /* if primal solution is found, then pass it to trysol heuristic */
         /* @note if trysol heuristic is was not present, then the solution found above, which can potentially be a MIP
          * primal feasible solution, will go to waste.
          */
         if( solfound && sepadata->heurtrysol != NULL )
         {
            SCIP_CALL( SCIPheurPassSolTrySol(scip, sepadata->heurtrysol, hardcutslpsol) );
         }

         /* if dual feasible, then fetch dual solution and reset Lagrangian multipliers based on it. otherwise, retain the
          * Lagrangian multipliers and simply initialize the new multipliers to zeroes. */
         if( SCIPlpiIsDualFeasible(sepadata->lpiwithhardcuts) )
         {
            SCIP_CALL( SCIPlpiGetSol(sepadata->lpiwithhardcuts, NULL, NULL, dualsol, NULL, NULL) );
            SCIP_CALL( updateDualVector(scip, sepadata, dualvector, &(dualsol[nrows]),
                     ngeneratedcurrroundcuts, 0, -1, -1, 0.0, NULL, ngeneratedcurrroundcuts, TRUE, 0.0, 0.0, 0.0, 0.0, -1,
                     &dualvecsdiffer, NULL) );
         }
         else
         {
            SCIP_CALL( updateDualVector(scip, sepadata, dualvector, dualvector, nsoftcuts, 0, -1, -1, 0.0, NULL,
                     ngeneratedcurrroundcuts, TRUE, 0.0, 0.0, 0.0, 0.0, -1, &dualvecsdiffer, NULL) );
         }
      }

      /* termination check */
      SCIP_CALL( checkMainLoopTermination(sepadata, cutoff, dualvecsdiffer, ngeneratedcurrroundcuts, nsoftcuts,
               nmaxgeneratedperroundcuts, ncurrroundlpiters, depth, &terminate) );
   }

   /* set the maximal denominator in rational representation of gomory cut and the maximal scale factor to
    * scale resulting cut to integral values to avoid numerical instabilities
    */
   /**@todo find better but still stable gomory cut settings: look at dcmulti, gesa3, khb0525, misc06, p2756 */
   /* @note: above todo was copied from sepa_gomory.c. So, if gomory code is changed, same changes can be done here. */
   maxdepth = SCIPgetMaxDepth(scip);
   if( depth == 0 )
   {
      maxdnom = 1000;
      maxscale = 1000.0;
   }
   else if( depth <= maxdepth/4 )
   {
      maxdnom = 1000;
      maxscale = 1000.0;
   }
   else if( depth <= maxdepth/2 )
   {
      maxdnom = 100;
      maxscale = 100.0;
   }
   else
   {
      maxdnom = 10;
      maxscale = 10.0;
   }

   /* filter cuts by calling cut selection algorithms and add cuts accordingly */
   /* @todo an idea is to filter cuts after every main iter */
   /* @todo we can remove !cutoff criterion by adding forcedcuts */
   if( !cutoff && (ngeneratedcurrroundcuts > 0) )
   {
      if( SCIPisGE(scip, sepadata->cutsfilterfactor, 1.0) )
      {
         nselectedcurrroundcuts = ngeneratedcurrroundcuts;
         SCIP_CALL( addCuts(scip, sepadata, generatedcurrroundcuts, nselectedcurrroundcuts, maxdnom, maxscale,
                  &naddedcurrroundcuts, &cutoff2) );
         cutoff = cutoff2;
      }
      else if( SCIPisPositive(scip, sepadata->cutsfilterfactor) )
      {
         nprocessedcuts = 0;
         for( int i = 0; i < sepadata->nmaxmainiters * (sepadata->nmaxsubgradientiters + 1); i++ )
         {
            if( ngeneratedcutsperiter[i] != 0 )
            {
               for( int j = 0; j < ngeneratedcutsperiter[i]; j++ )
                  cutindsperm[j] = j + nprocessedcuts;

               /* sort cut efficacies by fractionality */
               SCIPsortDownRealInt(&generatedcutefficacies[nprocessedcuts], cutindsperm, ngeneratedcutsperiter[i]);

               nselectedcuts = (int)SCIPceil(scip, sepadata->cutsfilterfactor * ngeneratedcutsperiter[i]);

               SCIP_CALL( addCuts(scip, sepadata, &generatedcurrroundcuts[nprocessedcuts], nselectedcuts, maxdnom,
                        maxscale, &naddedcurrroundcuts, &cutoff2) );
               cutoff = cutoff2;

               nprocessedcuts += ngeneratedcutsperiter[i];
            }
         }
      }
   }
   else if( ngeneratedcurrroundcuts > 0 )
   {
      nselectedcurrroundcuts = ngeneratedcurrroundcuts;
      SCIP_CALL( addCuts(scip, sepadata, generatedcurrroundcuts, nselectedcurrroundcuts, maxdnom, maxscale,
               &naddedcurrroundcuts, &cutoff2) );
   }

   /* add an aggregated cut based on best Lagrangian multipliers */
   if( sepadata->aggregatecuts && (ngeneratedcurrroundcuts > 0) && (bestdualvectorlen > 0) )
   {
      assert(bestdualvectorlen <= ngeneratedcurrroundcuts);
      SCIP_CALL( aggregateGeneratedCuts(scip, sepa, sepadata, generatedcurrroundcuts, bestdualvector, bestdualvectorlen,
               aggregatedcurrroundcuts, &naggregatedcurrroundcuts, &cutoff2) );
      cutoff = (!cutoff ? cutoff2 : cutoff);
      if( naggregatedcurrroundcuts > 0 )
      {
         SCIP_CALL( addCuts(scip, sepadata, aggregatedcurrroundcuts, naggregatedcurrroundcuts, maxdnom, maxscale,
                  &naddedcurrroundcuts, &cutoff2) );
         cutoff = (!cutoff ? cutoff2 : cutoff);
      }
   }

   if( cutoff )
   {
      *result = SCIP_CUTOFF;
   }
   else if( naddedcurrroundcuts > 0 )
   {
      *result = SCIP_SEPARATED;
   }
   else
   {
      *result = SCIP_DIDNOTFIND;
   }

   /* delete the LP representing the Lagrangian dual problem with fixed Lagrangian multipliers */
   SCIP_CALL( deleteLPWithSoftCuts(scip, sepadata) );

   /* release the rows in aggregatedcurrroundcuts */
   for( int i = 0; i < naggregatedcurrroundcuts; ++i )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(aggregatedcurrroundcuts[i])) );
   }
   /* release the rows in generatedcurrroundcuts */
   for( int i = 0; i < ngeneratedcurrroundcuts; ++i )
   {
      SCIP_CALL( SCIPreleaseRow(scip, &(generatedcurrroundcuts[i])) );
   }

   for( int i = 0; i < sepadata->nmaxmainiters * (sepadata->nmaxsubgradientiters + 1); i++ )
   {
      ngeneratedcutsperiter[i] = 0;
   }
   for( int i = 0; i < nrows; i++ )
   {
      dualsol[i] = 0.0;
   }
   for( int i = 0; i < nmaxgeneratedperroundcuts; i++ )
   {
      dualsol[nrows + i] = 0.0;
      dualvector[i] = 0.0;
      bestdualvector[i] = 0.0;
   }
   /* free memory */
   SCIPfreeBufferArray(scip, &origobjcoefs);
   SCIPfreeBufferArray(scip, &hardcutslpsolvals);
   SCIP_CALL( SCIPfreeSol(scip, &hardcutslpsol) );
   SCIPfreeBufferArray(scip, &softcutslpsolvals);
   SCIP_CALL( SCIPfreeSol(scip, &softcutslpsol) );
   SCIPfreeCleanBufferArray(scip, &ngeneratedcutsperiter);
   SCIPfreeBufferArray(scip, &cutindsperm);
   SCIPfreeBufferArray(scip, &generatedcutefficacies);
   SCIPfreeBufferArray(scip, &aggregatedcurrroundcuts);
   SCIPfreeBufferArray(scip, &generatedcurrroundcuts);
   SCIPfreeCleanBufferArray(scip, &dualvector);
   SCIPfreeCleanBufferArray(scip, &bestdualvector);
   SCIPfreeCleanBufferArray(scip, &dualsol);

   return SCIP_OKAY;
}

/** creates separator data */
static
SCIP_RETCODE sepadataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA**       sepadata            /**< separator data structure */
   )
{
   assert(scip != NULL);
   assert(sepadata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, sepadata) );
   BMSclearMemory(*sepadata);

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/** copy method for separator plugins (called when SCIP copies plugins) */
static
SCIP_DECL_SEPACOPY(sepaCopyLagromory)
{  /*lint --e{715}*/
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);

   /* call inclusion method of constraint handler */
   SCIP_CALL( SCIPincludeSepaLagromory(scip) );

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeLagromory)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIP_CALL( sepadataFree(scip, &sepadata) );
   SCIPsepaSetData(sepa, NULL);

   return SCIP_OKAY;
}


/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitLagromory)
{
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(scip != NULL);
   assert(sepadata != NULL);

   /* create and initialize random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &(sepadata->randnumgen), RANDSEED, TRUE) );

   /* find trysol heuristic */
   if ( sepadata->heurtrysol == NULL )
   {
      sepadata->heurtrysol = SCIPfindHeur(scip, "trysol");
   }

   return SCIP_OKAY;
}

/** deinitialization method of separator (called before transformed problem is freed) */
static
SCIP_DECL_SEPAEXIT(sepaExitLagromory)
{  /*lint --e{715}*/
   SCIP_SEPADATA* sepadata;

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPfreeRandom(scip, &(sepadata->randnumgen));

   return SCIP_OKAY;
}

/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpLagromory)
{
   /*lint --e{715}*/

   SCIP_SEPADATA* sepadata;
   SCIP_SOL* bestsol;
   SCIP_COL** cols;
   SCIP_COL* col;
   SCIP_VAR* var;
   SCIP_Real ubparam;
   SCIP_Real dualdegeneracyrate;
   SCIP_Real varconsratio;
   SCIP_Real threshold1;
   SCIP_Real threshold2;
   int nrows;
   int ncols;
   int ncalls;
   int runnum;

   assert(sepa != NULL);
   assert(strcmp(SCIPsepaGetName(sepa), SEPA_NAME) == 0);
   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   assert(scip != NULL);

   ncalls = SCIPsepaGetNCallsAtNode(sepa);
   runnum = SCIPgetNRuns(scip);
   dualdegeneracyrate = 0.0;
   varconsratio = 0.0;

   /* only generate Lagromory cuts if we are not close to terminating */
   if( SCIPisStopped(scip) )
      return SCIP_OKAY;

   /* only call the separator starting sepadata->minrestart runs */
   if( (sepadata->minrestart >= 1) && (runnum < sepadata->minrestart + 1) )
      return SCIP_OKAY;

   /* only call the separator a given number of times at each node */
   if( (depth == 0 && sepadata->maxroundsroot >= 0 && ncalls >= sepadata->maxroundsroot)
      || (depth > 0 && sepadata->maxrounds >= 0 && ncalls >= sepadata->maxrounds) )
      return SCIP_OKAY;

   /* only generate cuts if an optimal LP solution is at hand */
   if( SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL )
      return SCIP_OKAY;

   /* only generate cuts if the LP solution is basic */
   if( ! SCIPisLPSolBasic(scip) )
      return SCIP_OKAY;

   /* get LP data */
   SCIP_CALL( SCIPgetLPColsData(scip, &cols, &ncols) );
   assert(cols != NULL);

   nrows = SCIPgetNLPRows(scip);

   /* return if LP has no columns or no rows */
   if( ncols == 0 || nrows == 0 )
      return SCIP_OKAY;

   /* return if dual degeneracy metrics are below threshold values */
   threshold1 = sepadata->dualdegeneracyratethreshold;
   threshold2 = sepadata->varconsratiothreshold;
   SCIP_CALL( SCIPgetLPDualDegeneracy(scip, &dualdegeneracyrate, &varconsratio) );
   if( (dualdegeneracyrate < threshold1) && (varconsratio < threshold2) )
      return SCIP_OKAY;

   bestsol = SCIPgetBestSol(scip);
   ubparam = 0.0;

   /* if a MIP primal solution exists, use it to estimate the optimal value of the Lagrangian dual problem */
   if( bestsol != NULL )
   {
      for( int i = 0; i < ncols; ++i )
      {
         col = cols[i];
         assert(col != NULL);

         var = SCIPcolGetVar(col);

         ubparam += SCIPgetSolVal(scip, bestsol, var) * SCIPcolGetObj(col);
      }

      ubparam += SCIPgetTransObjoffset(scip);
   }
   /* if a MIP primal solution does not exist, then use the node relaxation's LP solutin to estimate the optimal value
    * of the Lagrangian dual problem
    */
   else
   {
      for( int i = 0; i < ncols; ++i )
      {
         col = cols[i];
         assert(col != NULL);

         ubparam += SCIPcolGetPrimsol(col) * SCIPcolGetObj(col);
      }

      ubparam += SCIPgetTransObjoffset(scip);
      ubparam *= SCIPisPositive(scip, ubparam) ? sepadata->ubparamposfactor : sepadata->ubparamnegfactor;
   }

   /* the main separation function of the separator */
   SCIP_CALL( separateCuts(scip, sepa, sepadata, ubparam, depth, (allowlocal && sepadata->allowlocal), result) );

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the Lagromory separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaLagromory(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create separator data */
   SCIP_CALL( sepadataCreate(scip, &sepadata) );

   sepadata->heurtrysol = NULL;

   /* include separator */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
            SEPA_USESSUBSCIP, SEPA_DELAY, sepaExeclpLagromory, NULL, sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyLagromory) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeLagromory) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitLagromory) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitLagromory) );

   /* add separator parameters */
   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/away", "minimal integrality violation of a basis "
            "variable to try separation", &sepadata->away, FALSE, DEFAULT_AWAY, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/rootlpiterlimitfactor", "factor w.r.t. root node LP "
            "iterations for maximal separating LP iterations in the root node (negative for no limit)",
            &sepadata->rootlpiterlimitfactor, TRUE, DEFAULT_ROOTLPITERLIMITFACTOR, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/totallpiterlimitfactor", "factor w.r.t. root node LP "
            "iterations for maximal separating LP iterations in the tree (negative for no limit)",
            &sepadata->totallpiterlimitfactor, TRUE, DEFAULT_TOTALLPITERLIMITFACTOR, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/perroundlpiterlimitfactor", "factor w.r.t. root node LP "
            "iterations for maximal separating LP iterations per separation round (negative for no limit)",
            &sepadata->perroundlpiterlimitfactor, TRUE, DEFAULT_PERROUNDLPITERLIMITFACTOR, -1.0, SCIP_REAL_MAX, NULL,
            NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/perroundcutsfactorroot", "factor w.r.t. number of integer "
            "columns for number of cuts separated per separation round in root node", &sepadata->perroundcutsfactorroot,
            TRUE, DEFAULT_PERROUNDCUTSFACTORROOT, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/perroundcutsfactor", "factor w.r.t. number of integer "
            "columns for number of cuts separated per separation round at a non-root node",
            &sepadata->perroundcutsfactor, TRUE, DEFAULT_PERROUNDCUTSFACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/totalcutsfactor", "factor w.r.t. number of integer "
            "columns for total number of cuts separated", &sepadata->totalcutsfactor, TRUE, DEFAULT_TOTALCUTSFACTOR, 0.0,
            SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/muparaminit", "initial value of the mu parameter (factor "
      "for step length)", &sepadata->muparaminit, TRUE, DEFAULT_MUPARAMINIT, 0.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/muparamlb", "lower bound of the mu parameter (factor for "
      "step length)", &sepadata->muparamlb, TRUE, DEFAULT_MUPARAMLB, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/muparamub", "upper bound of the mu parameter (factor for "
      "step length)", &sepadata->muparamub, TRUE, DEFAULT_MUPARAMUB, 1.0, 10.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/mubacktrackfactor", "factor of mu while backtracking the "
            "mu parameter (factor for step length)", &sepadata->mubacktrackfactor, TRUE, DEFAULT_MUBACKTRACKFACTOR,
            0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/muslab1factor", "factor of mu parameter (factor for step "
      "length) for larger increment" , &sepadata->muslab1factor, TRUE, DEFAULT_MUSLAB1FACTOR, 0.0, SCIP_REAL_MAX, NULL,
            NULL));

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/muslab2factor", "factor of mu parameter (factor for step "
      "length) for smaller increment", &sepadata->muslab2factor, TRUE, DEFAULT_MUSLAB2FACTOR, 0.0, SCIP_REAL_MAX, NULL,
            NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/muslab3factor", "factor of mu parameter (factor for step "
      "length) for reduction", &sepadata->muslab3factor, TRUE, DEFAULT_MUSLAB3FACTOR, 0.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/deltaslab1ub", "factor of delta deciding larger increment "
            "of mu parameter (factor for step length)", &sepadata->deltaslab1ub, TRUE, DEFAULT_DELTASLAB1UB, 0.0, 1.0,
            NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/deltaslab2ub", "factor of delta deciding smaller "
            "increment of mu parameter (factor for step length)", &sepadata->deltaslab2ub, TRUE, DEFAULT_DELTASLAB2UB,
            0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/ubparamposfactor", "factor for positive upper bound used "
            "as an estimate for the optimal Lagrangian dual value", &sepadata->ubparamposfactor, TRUE, DEFAULT_UBPARAMPOSFACTOR,
            1.0, 100.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/ubparamnegfactor", "factor for negative upper bound used "
            "as an estimate for the optimal Lagrangian dual value", &sepadata->ubparamnegfactor, TRUE, DEFAULT_UBPARAMNEGFACTOR,
            0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/perrootlpiterfactor", "factor w.r.t. root node LP "
            "iterations for iteration limit of each separating LP (negative for no limit)",
            &sepadata->perrootlpiterfactor, TRUE, DEFAULT_PERROOTLPITERFACTOR, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/perlpiterfactor", "factor w.r.t. non-root node LP "
            "iterations for iteration limit of each separating LP (negative for no limit)", &sepadata->perlpiterfactor, TRUE,
            DEFAULT_PERLPITERFACTOR, -1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/cutsfilterfactor", "fraction of generated cuts per "
            "explored basis to accept from separator", &sepadata->cutsfilterfactor, TRUE, DEFAULT_CUTSFILTERFACTOR, 0.0,
            1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/radiusinit", "initial radius of the ball used in "
            "stabilization of Lagrangian multipliers", &sepadata->radiusinit, TRUE, DEFAULT_RADIUSINIT, 0.0, 1.0, NULL,
            NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/radiusmax", "maximum radius of the ball used in "
            "stabilization of Lagrangian multipliers", &sepadata->radiusmax, TRUE, DEFAULT_RADIUSMAX, 0.0, SCIP_REAL_MAX,
            NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/radiusmin", "minimum radius of the ball used in "
            "stabilization of Lagrangian multipliers", &sepadata->radiusmin, TRUE, DEFAULT_RADIUSMIN, 0.0, SCIP_REAL_MAX,
            NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/constant", "a constant for stablity center based "
            "stabilization of Lagrangian multipliers", &sepadata->constant, TRUE, DEFAULT_CONST, 2.0, SCIP_REAL_MAX,
            NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/radiusupdateweight", "multiplier to evaluate cut "
            "violation score used for updating ball radius", &sepadata->radiusupdateweight, TRUE,
            DEFAULT_RADIUSUPDATEWEIGHT, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/dualdegeneracyratethreshold", "minimum dual degeneracy "
            "rate for separator execution", &sepadata->dualdegeneracyratethreshold, FALSE,
            DEFAULT_DUALDEGENERACYRATETHRESHOLD, 0.0, 1.0, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "separating/" SEPA_NAME "/varconsratiothreshold", "minimum variable-constraint "
            "ratio on optimal face for separator execution", &sepadata->varconsratiothreshold, FALSE,
            DEFAULT_VARCONSRATIOTHRESHOLD, 1.0, SCIP_REAL_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/muparamconst", "is the mu parameter (factor for step "
      "length) constant?" , &sepadata->muparamconst, TRUE, DEFAULT_MUPARAMCONST, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/separaterows", "separate rows with integral slack?",
            &sepadata->separaterows, TRUE, DEFAULT_SEPARATEROWS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/sortcutoffsol", "sort fractional integer columns"
            "based on fractionality?" , &sepadata->sortcutoffsol, TRUE, DEFAULT_SORTCUTOFFSOL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/sidetypebasis", "choose side types of row (lhs/rhs) "
            "based on basis information?", &sepadata->sidetypebasis, TRUE, DEFAULT_SIDETYPEBASIS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/dynamiccuts", "should generated cuts be removed from "
            "LP if they are no longer tight?", &sepadata->dynamiccuts, FALSE, DEFAULT_DYNAMICCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/makeintegral", "try to scale all cuts to integral "
            "coefficients?", &sepadata->makeintegral, TRUE, DEFAULT_MAKEINTEGRAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/forcecuts", "force cuts to be added to the LP?",
            &sepadata->forcecuts, TRUE, DEFAULT_FORCECUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/delayedcuts", "should cuts be added to the delayed cut "
            "pool", &sepadata->delayedcuts, TRUE, DEFAULT_DELAYEDCUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/allowlocal", "should locally valid cuts be generated?",
            &sepadata->allowlocal, TRUE, DEFAULT_ALLOWLOCAL, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "separating/" SEPA_NAME "/aggregatecuts", "aggregate all generated cuts using the "
            "Lagrangian multipliers?", &sepadata->aggregatecuts, TRUE, DEFAULT_AGGREGATECUTS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxrounds", "maximal number of separation rounds per node "
            "(-1: unlimited)", &sepadata->maxrounds, FALSE, DEFAULT_MAXROUNDS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/maxroundsroot", "maximal number of separation rounds in "
            "the root node (-1: unlimited)", &sepadata->maxroundsroot, FALSE, DEFAULT_MAXROUNDSROOT, -1, INT_MAX, NULL,
            NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/perroundnmaxlpiters", "maximal number of separating LP "
            "iterations per separation round (-1: unlimited)", &sepadata->perroundnmaxlpiters, FALSE,
            DEFAULT_PERROUNDMAXLPITERS, -1, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nmaxcutsperlp", "maximal number of cuts separated per "
            "Lagromory LP in the non-root node", &sepadata->nmaxcutsperlp, FALSE, DEFAULT_PERLPMAXCUTS, 0, INT_MAX, NULL,
            NULL));

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nmaxcutsperlproot", "maximal number of cuts separated per "
            "Lagromory LP in the root node", &sepadata->nmaxcutsperlproot, FALSE, DEFAULT_PERLPMAXCUTSROOT, 0, INT_MAX,
            NULL, NULL));

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nmaxmainiters", "maximal number of main loop iterations of "
            "the relax-and-cut algorithm", &sepadata->nmaxmainiters, TRUE, DEFAULT_MAXMAINITERS, 0, INT_MAX, NULL, NULL)
         );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nmaxsubgradientiters", "maximal number of subgradient loop "
            "iterations of the relax-and-cut algorithm", &sepadata->nmaxsubgradientiters, TRUE,
            DEFAULT_MAXSUBGRADIENTITERS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/cutgenfreq", "frequency of subgradient iterations for "
            "generating cuts", &sepadata->cutgenfreq, TRUE, DEFAULT_CUTGENFREQ, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/cutaddfreq", "frequency of subgradient iterations for "
            "adding cuts to objective function", &sepadata->cutaddfreq, TRUE, DEFAULT_CUTADDFREQ, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nmaxlagrangianvalsforavg", "maximal number of iterations "
            "for rolling average of Lagrangian value", &sepadata->nmaxlagrangianvalsforavg, TRUE,
            DEFAULT_MAXLAGRANGIANVALSFORAVG, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/nmaxconsecitersformuupdate", "consecutive number of "
            "iterations used to determine if mu needs to be backtracked", &sepadata->nmaxconsecitersformuupdate, TRUE,
            DEFAULT_MAXCONSECITERSFORMUUPDATE, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/projectiontype", "the ball into which the Lagrangian multipliers "
            "are projected for stabilization (0: no projection, 1: L1-norm ball projection, 2: L2-norm ball projection, 3: "
            "L_inf-norm ball projection)", &sepadata->projectiontype, TRUE, DEFAULT_PROJECTIONTYPE, 0, 3, NULL, NULL));

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/stabilitycentertype", "type of stability center for "
            "taking weighted average of Lagrangian multipliers for stabilization (0: no weighted stabilization, 1: best "
            "Lagrangian multipliers)", &sepadata->stabilitycentertype, TRUE, DEFAULT_STABILITYCENTERTYPE, 0, 1, NULL,
            NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/optimalfacepriority", "priority of the optimal face for "
            "separator execution (0: low priority, 1: medium priority, 2: high priority)",
            &sepadata->optimalfacepriority, TRUE, DEFAULT_OPTIMALFACEPRIORITY, 0, 2, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "separating/" SEPA_NAME "/minrestart", "minimum restart round for separator "
            "execution (0: from beginning of the instance solving, >= n with n >= 1: from restart round n)",
            &sepadata->minrestart, TRUE, DEFAULT_MINRESTART, 0, INT_MAX, NULL, NULL) );

  return SCIP_OKAY;
}
