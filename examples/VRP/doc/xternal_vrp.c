/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal_vrp.c
 * @brief  main document page of VRP example
 * @author Andreas Bley
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page VRP_MAIN Vehicle Routing
 * @version  0.1
 * @author   Andreas Bley
 *
 *
 * We want to solve the vehicle routing problem (VRP) on a graph \f$G = (V,E)\f$ with \f$V = J \cup {d}\f$, where
 * \f$d\f$ is the depot and the distances are given by the length function \f$l_e: E \to R_{\ge 0}\f$.
 *
 * Consider the MIP formulation
 *
 * \f[
 *  \begin{array}[t]{rll}
 *    \min &  \displaystyle \sum_{e \in E} l_e y_e \\
 *         & & \\
 *   s.t.  & -y_e + \sum_{t \in T_k} a^t_e x_t  \leq 0, &  \forall e \in E\\
 *         &  \displaystyle \sum_{t \in T_k} a^t_j x_t = 1, &  \forall j \in J \\
 *         &  y(\delta(j)) = 2, &  \forall j \in J \\
 *         &  y_e \in \{0,1,2\},  & \forall e \in E \\
 *         &  x_t  \in [0,1], & \forall t \in T_k
 *  \end{array}
 * \f]
 *
 * where \f$T_k\f$ is the set of tours visiting at most \f$k\f$ customers with repetitions of customers allowed and
 * \f$a^t_e (a^t_j)\f$ counts how often edge e (node j) is traversed in \f$t \in T_k\f$. The model contains two types of
 * variables, namely \f$ x \f$ which selects tours fractionally and \f$ y \f$ which indicates which edges of the graph
 * are in at least one selected tour. Note that it is possible to use an edge as a forward and backward edge in a tour.
 * This is necessary to ensure that a customer \f$ j \f$ with \f$ |\delta(j)| = 1 \f$ can be served.
 *
 * Since the number of tours can be exponential in the size of the graph, the algorithm starts with some subset \f$
 * \bar{T} \subseteq T_k \f$ and adds further tours during the solution process. This way it tries to improve the
 * current LP solution.
 *
 * Let \f$ \lambda_e \f$ and \f$ \gamma_i \f$ be the dual multipliers for the first and seconds constraint of the
 * MIP and we define the costs of a tour \f$ T \in T_k \f$ as:
 * \f[
 *   C(T) := \sum_{e \in E(T)} \lambda_e - \sum_{j \in V(T)} \gamma_j
 * \f]
 *
 * The resulting pricing problem \f$ \min_{T \in T_k} C(T) \f$ can be solved with dynamic programming. The algorithm is
 * similar to Dijkstra's shortest path algorithm if we shift the the costs \f$ \gamma_j \f$ from the nodes to the edges
 * of \f$ G \f$.
 *
 * Branching decisions on the variables \f$ y \f$ modify the pricing problem only slightly. The branch \f$ y_e = 0\f$
 * forbids \f$ e \f$ to be contained in a tour which can be easily realized if we remove \f$ e \f$ from \f$ E \f$. The
 * branch \f$ y_e \ge 1 \f$ does not have an impact on the pricing problem.
 *
 * Further information about the pricing routine and the dynamic program can be found in the documentation of the
 * corresponding files.
 *
 * The pricer pricer_vrp.cpp shows how to perform column generation in SCIP and how to solve the above described pricing
 * problem which uses an implementation of a priority queue implemented in pqueue.h. In main_vrp.cpp we read the
 * instance, create all necessary data and set up SCIP.
 */
