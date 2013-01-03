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

/**@file   xternal.c
 * @brief  main document page of VRP example
 * @author Andreas Bley
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 * @version  0.1
 * @author   Andreas Bley
 *
 *
 * We want to solve the vehicle routing problem (VRP) on a graph \f$G = (V,E)\f$ with
 * \f$V = J \cup {d}\f$, where \f$d\f$ is the depot and the distances are given by the
 * length function \f$l_e: E \to R_{\ge 0}\f$.
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
 * where \f$T_k\f$ is the set of tours visiting at most \f$k\f$ customers
 * with repetitions of customers allowed and \f$a^t_e (a^t_j)\f$ counts how often
 * edge e (node j) is traversed in \f$t \in T_k\f$.
 */
