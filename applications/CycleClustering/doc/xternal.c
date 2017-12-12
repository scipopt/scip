/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author Leon Eifler
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 * @version  0.1
 * @author   Leon Eifler
 *
 * This application can be used to solve the cycle clustering problem as described in
 *
 * "Mixed-Integer Programming for Cycle Detection in Non-reversible Markov Processes"@n
 * by Isabel Beckenbach, Leon Eifler, Konstantin Fackeldey, Ambros Gleixner, Andreas Grever, Marcus Weber, and Jakob Witzig, @n
 * Multiscale Modling and Simulation, 2016 (accepted for publication,  preprint available as <a href="https://opus4.kobv.de/opus4-zib/frontdoor/index/index/docId/6035">ZIB-Report 16-39</a> ).
 *
 * The input format is an \f$n \times n\f$ - matrix \f$Q\f$ of unconditional transition probabilities with a header of the form
 * "# p nstates ncluster"; nstates is the size of the matrix, ncluster the desired number of clusters; the name of the file must end with ".spa".
 *
 * The cycle clustering problem is the following:
 *
 * Consider a set of states \f$ \mathcal B = \{1,\ldots,n\}\f$ and a set of clusters \f$\mathcal{C}=\{1,\ldots,m\}\f$.
 * Let \f$Q \in \mathbb{R}^{n \times n}\f$ with entries \f$ q_{ij}\f$. Then the problem is given by the MINLP
 *
 * @image html model.jpg
 *
 * Further information about particular modules like heuristics and separation routines
 * can be found in the documentation of the corresponding files.
 *
 * Compiling the CycleClustering application
 * -----------------------------------------
 *
 * See the @ref INSTALL "Install file"
 *
 *
 */
