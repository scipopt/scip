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

/**@mainpage multi objective  example
 * @version  0.1
 * @author   Timo Strunk
 *
 * This example is a method for solving mixed integer programs with multiple objective functions (MOMIP).
 *
 * There are different concepts of optimality in multi criterial optimization and in general there is more than one optimal solution.
 *
 * Given a problem of the form
 *
 * \f$ \min Cx \f$ 
 * \f$ Ax \leq b \f$
 * 
 * the algorithm computes all feasible solutions which minimize \f$ w^TCx \f$ for some \f$ w \geq 0 \f$.
 * From the example directory root you can call the program with the command
 * 
 * bin/multiopt <instance>.mop [-v] [-t <time in seconds>]
 *
 * The problem file has to be in mps format, except that more than one row is marked "N" indicating a cost vector.
 * The program expects a minimization problem.  The file ending has to be ".mop".
 * 
 * Tip: If you want to use ZIMPL to generate your instance, model the extra objectives as equality constraints with a right hand side of 0.
 * Then use ZIMPL to translate your instance to MPS. Now open the output file with an editor, find the names of the extra objectives and
 * change the mark "E" to "N" manually.
 * 
 * Options:\\
 * -v	output SCIP messages as well as computation results\\
 * -t	set a time limit.
 * 
 * 
 * The example implements a version of the weight space decomposition algorithm specified in 
 * "An Exact Algorithm for Finding Extreme Supported Nondominated Points of Multiobjective Mixed Integer Programs"
 * by Özgür Özpeynirci and Murat Köksalan
 * Manage. Sci., 56(12):2302{2315, December 2010}
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
 * To understand the algorithm, consider this: The problem of finding the minimal weighted objective value a for every weight \f$ w \f$
 * is equivalent to computing the set of all pairs \f$ (a,w) \f$ such that there exists a solution \f$ x \f$ with \f$ a \geq w^TCx\f$. \f$ P \f$ 
 * is exactely the complement of that set. 
 */
