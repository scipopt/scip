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

/**@file   xternal_lop.c
 * @brief  main document page
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page LOP_MAIN Linear Ordering
 * @version  1.0
 * @author   Marc Pfetsch
 *
 * The linear ordering gives another example for setting up a
 * constraint handler.
 *
 * The linear ordering problem is the following:
 *
 * Given a positive integer \f$ n \f$  and an \f$ n \times n  \f$ matrix \f$ W \f$  the goal is to
 * find a linear order of \f$ \{1, \dots, n\}  \f$  such that the sum of weights \f$
 * w_{ij}\f$ for all pairs in which \f$ x \f$  comes before \f$ j \f$  in the order is
 * maximized.
 *
 * We use the integer programming following model: We have binary
 * variables \f$ x_{ij}\f$ for all pairs \f$ (i,j)\f$ with \f$ i \neq
 * j\f$, where \f$ x_{ij} = 1\f$ if and only if \f$ i \f$  comes before \f$ j \f$  in the
 * encoded order. The basic model is then:
 * \f[
 *     \max \{ \sum_{i,j} w_{ij} x_{ij}\;:\; x_{ij} + x_{ji} = 1 \mbox{ for all } j \neq i\}.
 * \f]
 * To ensure that x encodes a linear order one has to add the
 * following @em triangle @em inequalities:
 * \f[
 *     x_{ij} + x_{jk} + x_{ki} \leq 2 \quad\mbox{for all }i,j,k.
 * \f]
 * Using the equations above, one can of course eliminate half of the
 * variables (and change the triangle inequalies accordingly), but we
 * do not do this explicitly in order to keep a simpler
 * formulation. In fact, SCIP will do some of the eliminations
 * automatically.
 *
 * The following files provide the example code:
 * - cmain.c: Here the main function is located. It sets up SCIP, the
 * linear order project, and solves the problem.
 * - reader_lop.c: this file provides code for reading the corresponding weight matrix
 * and setting up the above model.
 * - cons_lop.c: contains the constraint handler that takes care of the
 * equations and the triangle inequalities.
 * - genRandomLOPInstance.c: problem generator (see \ref LOP_PROBLEMGENERATORUSEIT "below")
 *
 *
 * To use the problem generator you have do two things. First
 * \ref LOP_PROBLEMGENERATORCOMPILE "compile the generator" and second \ref LOP_PROBLEMGENERATORUSEIT "use it".
 *
 * @section LOP_PROBLEMGENERATORCOMPILE Compile the Problem Generator
 *
 * Call the command
 *
 * <code>make genRandomLOPInstance</code>
 *
 * in main directory of the example. This will create a binary in the <code>bin/</code> directory
 * with the name <code>genRandomLOPInstance</code>.
 *
 * @section LOP_PROBLEMGENERATORUSEIT Use the Problem Generator
 *
 * The problem generator needs three parameter:
 * -# the name of the file to create
 * -# matrix dimension
 * -# the range of the integer values
 *
 * For example the call (in the main directory of the example)
 *
 * <code>bin/genRandomLOPInstance instance 10 6</code>
 *
 * produces a file named "instance" containing a matrix of dimension 10x10 with entries between 0 and 6.
 *
 * Installation
 * ============
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 *
 */
