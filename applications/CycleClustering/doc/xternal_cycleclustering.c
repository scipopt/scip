/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   xternal_cycleclustering.c
 * @brief  main document page
 * @author Leon Eifler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page CYCLECLUSTERING_MAIN CycleClustering
 * @version  0.1
 * @author   Leon Eifler
 *
 * This application can be used to solve the cycle clustering problem as described in
 *
 * "Mixed-Integer Programming for Cycle Detection in Non-reversible Markov Processes"@n
 * by Isabel Beckenbach, Leon Eifler, Konstantin Fackeldey, Ambros Gleixner, Andreas Grever, Marcus Weber, and Jakob Witzig, @n
 * <em> Multiscale Modeling and Simulation </em> , 2016 (accepted for publication,  preprint available as <a href="https://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/6035">ZIB-Report 16-39</a> ).
 *
 * The input format is an \f$n \times n\f$ - matrix \f$Q\f$ of unconditional transition probabilities with a header of the form
 * "# p nstates ncluster"; nstates is the size of the matrix, ncluster the desired number of clusters; the name of the file must end with ".spa".
 *
 * The cycle clustering problem is the following:
 *
 * Consider a set of states \f$ \mathcal B = \{1,\ldots,n\}\f$ and a set of clusters \f$\mathcal{C}=\{1,\ldots,m\}\f$.
 * Let \f$Q \in \mathbb{R}^{n \times n}\f$ with entries \f$ q_{ij}\f$. Then the problem is given by the MINLP
 *
 * 	\f{align*}{
 *	\max \ \ \ \ \ \sum_{t \in \mathcal{K}}f_t \ + \ &\alpha \cdot \sum_{t \in \mathcal{K}} g_t  \notag\\
 *	\text{s.t.} \quad \sum_{t \in \mathcal{K}} x_{it} &= 1 &&  \text{ for all } i \in \mathcal{S}  \\
 *	\sum_{i \in \mathcal{S}} x_{it} &\ge 1 &&  \text{ for all } t \in \mathcal{K} \label{eq:setcover} \\
 *	g_t &=  \sum_{\substack{i,j \in \mathcal{S}\\  i < j}} (q_{ij} + q_{ji}) x_{it} x_{jt} &&  \text{ for all } t \in \mathcal{K}\\
 *	f_t &= {\sum_{\substack{i,j \in \mathcal{S},\\  i \neq j}} (q_{ij}-q_{ji}) x_{it} x_{j \phi(t)}} &&  \text{ for all } t \in \mathcal{K}  \\
 *	x_{it} &\in \{0,1\} && \text{ for all } t \in \mathcal{K},  i \in \mathcal{S} \notag \\
 *	f_t, g_t &\in \mathbb{R}_{\geq 0} && \text{ for all } t \in \mathcal{K}. \notag
 *	\f}
 *
 * Further information about particular modules like heuristics and separation routines
 * can be found in the documentation of the corresponding files.
 *
 * Installation
 * ------------
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 *
 *
 */
