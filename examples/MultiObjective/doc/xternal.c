/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
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
 * @author Timo Strunk
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Multi Objective Mixed Integer Programming using a Lifted Weight Space Algorithm
 * @version  0.1
 * @author   Timo Strunk
 *
 * This example is a method for solving mixed integer programs with multiple objective functions (MOMIP).
 *
 * There are different concepts of optimality in multi criterial optimization and in general there is more than one optimal solution.
 *
 * Given a problem of the form
 * @TODO integrality restrictions
 * \f$ \min Cx \f$ 
 * \f$ Ax \leq b \f$
 * 
 * the algorithm computes all feasible solutions which minimize \f$ w^TCx \f$ for some \f$ w \geq 0 \f$.
 * 
 * @TODO
 * - specify path to a SCIP installation
 * - install the lemon library headers (link?)
 * - specify lemon include path and library path in lib (lemoninc, lemonlib.a)
 *
 * After the sources were sucessfully built, you can call the program from the example directory root via
 *
 *
 *     bin/multiopt <instance>.mop [<SCIP-settings-file>.set]
 *
 *
 * The optional settings file argument can be used to pass customized settings for the solving process inside SCIP
 * as the underlying MIP solver. Since solving times and other statistics are accumulated over during the sequence
 * of MIP optimization processes, a time limit for SCIP will be used as a threshold for the overall solving process
 * of the Lifted Weight Space algorithm.
 *
 *
 * The problem file has to be in mps format, except that more than one row is marked "N" indicating a cost vector.
 * The program expects a minimization problem.  The file ending has to be ".mop".
 *
 *  @TODO exact specifications as list
 *
 * Tip: If you want to use ZIMPL to generate your instance, model the extra objectives as equality constraints with a right hand side of 0.
 * Then use ZIMPL to translate your instance to MPS. Now open the output file with an editor, find the names of the extra objectives and
 * change the mark "E" to "N" manually.
 * 
 * 
 * 
 * 
 * The program implements the Lifted Weightspace Algorithm.  Basically, in each iteration, the algorithm chooses a weight 
 * \f$ w \in W = \{ w \geq 0 | \sum w = 1 \} \f$ and uses
 * SCIP to solve the single objective problem with \f$ c = w^TC \f$ .
 * To calculate the weights, the algorithm keeps track of the lifted weight space polyhedron 
 * \f$ P = \{ (a,w) \in R \times W | a \leq w^TCx,  x \in X \} \f$
 * where \f$ X \f$ is the set of solutions found so far. \f$ P \f$ is a convex polyhedron.
 * In each iteration, the algorithm chooses a vertex \f$ (a,w) \f$ of \f$ P \f$ such that \f$ w \f$ has not yet been used in weighted optimization.
 * Every time a new solution is found, \f$ P \f$ is updated.
 * 
 * To understand the algorithm, consider this: The problem of finding the minimal weighted objective value \f$ a_w \f$ for every weight \f$ w \f$
 * is equivalent to computing the set of all pairs \f$ (a,w) \f$ such that there exists a solution \f$ x \f$ with \f$ a \geq w^TCx\f$. \f$ P \f$ 
 * is exactely the complement of that set.
 * @TODO add one or two pictures
 *
 * @TODO explanation of classes with reference reader_mop.h
 *
 * @TODO add eigene referenz
 *
 *
 * Further information about multi objective mixed integer programs can be found in
 *
 * Özgür Özpeynirci and Murat Köksalan,
 * An Exact Algorithm for Finding Extreme Supported Nondominated Points of Multiobjective Mixed Integer Programs.
 * _Management Science_, 56(12):2302--2315, 2010.
 */
