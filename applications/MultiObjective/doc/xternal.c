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
 * A multi objective mixed integer problem is given by a cost matrix \f$C \in Q^{p \times (n+m)} \f$,
 * a constraint matrix \f$ A \in Q^{\ell \times (n+m)} \f$ and a right hand vector \f$ b \in Q^\ell \f$
 * for some integers \f$ n, m, \ell \f$.
 * It has the form
 *
 * \f{eqnarray*}{
 *    \min Cx \\
 *    Ax \leq b \\
 *    x \in Z^n \times R^m
 * \f}
 *
 * The algorithm computes all feasible solutions which minimize \f$ w^TCx \f$ for some \f$ w \geq 0 \f$.
 * Those are called extremal supported non-dominated solutions.
 *
 * @section INSTALL Quick installation guide
 *
 * The program can be built following these steps:
 * - Install the <a href="http://lemon.cs.elte.hu/trac/lemon/">lemon</a> graph library for <code>c++</code>
 * - Install SCIP
 * - In the Makefile, make sure the variable <code>SCIPDIR</code> correctly points to the SCIP root directory
 * - Build the binary using <code>make</code>; follow the instructions for creating the necessary link to lemon
 * - If building the softlinks fails, go to the <code>lib</code> directory.
 *   Try manually creating a softlink named <code>lemon</code> to the directory containing the lemon header files,
 *   see also the <code>INSTALL</code>-file for further information
 *
 * @section RUNNING Running The Executable
 *
 * After the sources were sucessfully built, you can call the program from the example directory root via
 *
 *
 *     bin/multiopt <instance>.mop [<SCIP-settings-file>.set]
 *
 * The optional settings file argument can be used to pass customized settings for the solving process inside SCIP
 * as the underlying MIP solver. Some settings are handled differently from usual SCIP behaviour.
 * For example, if a time limit is given, it is applied to the whole MOMIP solving process, not every SCIP call individually.
 *
 * The input file has to follow some simple conventions
 * - It has to contain a problem in <a href="http://en.wikipedia.org/wiki/MPS_%28format%29">MPS</a> format
 * - The file extension must be <code>.mop</code>
 * - Every row marked <code>N</code> is treated as one of the cost vectors
 *
 * Tip: If you want to use ZIMPL to generate your instance, model the extra objectives as equality constraints with a right hand side of 0.
 * Then use ZIMPL to translate your instance to MPS. Now open the output file with an editor, find the names of the extra objectives and
 * change the mark <code>E</code> to <code>N</code> manually.
 *
 * @section ALGORITHM Algorithmic description
 *
 * The program implements the Lifted Weightspace Algorithm.  Basically, in each iteration, the algorithm chooses a weight
 * \f$ w \in W = \{ w \geq 0 | \sum w = 1 \} \f$ and uses
 * SCIP to solve the single objective problem with \f$ c = w^TC \f$ .
 * To calculate the weights, the algorithm keeps track of the lifted weight space polyhedron
 * \f$ P = \{ (a,w) \in R \times W | a \leq w^TCx,  x \in X \} \f$
 * where \f$ X \f$ is the set of solutions found so far. \f$ P \f$ is a convex polyhedron.
 * In each iteration, the algorithm chooses a vertex \f$ (a,w) \f$ of \f$ P \f$ such that \f$ w \f$ has not yet been used in weighted optimization.
 * Every time a new solution is found, \f$ P \f$ is updated.

 * \image html before_update.png "before update"
 * \image html after_update.png "after update"

 * To understand the algorithm, consider this: The problem of finding the minimal weighted objective value \f$ a_w \f$ for every weight \f$ w \f$
 * is equivalent to computing the set of all pairs \f$ (a,w) \f$ such that there exists a solution \f$ x \f$ with \f$ a \geq w^TCx\f$. \f$ P \f$
 * is exactly the complement of that set.
 *
 * @section LIBRARY The design of this example
 *
 * The functionality is distributed amongst several classes as follows:
 * - main.h reads the arguments, calls all other functions and produces text output.
 * - reader_mop.h parses a problem instance from a file; this is a modified version of the MPS reader in SCIP.
 * - Objectives.h contains the data structures for storing the objective function and methods for manipulating the weighted objective function.
 * - WeightedSolver.h is an abstract weight based iterative MOMIP solver, contains functions yielding statistical data about the solving process.
 * - LiftedWeightSpaceSolver.h is the implementation of the lifted weight space algorithm, inheriting from WeightedSolver.h.
 * - Skeleton.h contains the 1-skeleton of the lifted weight space polyhedron.
 * - WeightSpaceVertex.h represents individual vertices of the weight space polyhedron.
 *
 * The output from the program should look something like this:
 *
 *     reading parameter file <scipmip.set>
 *     reading problem file <data/tenfelde-podehl.mop>
 *     solving
 *                                       Weight                                    Cost    Time/s     Nodes   LP Iter     V_new    V_proc
 *      [     0.33333     0.33333     0.33333 ] [          11          11          14 ]      0.01         1         6         3         0
 *      [           1           0           0 ]                                       -      0.02         1         5         0         1
 *      [           0           1           0 ] [          15           9          17 ]      0.01         1         5         2         1
 *      [           0           0           1 ] [          19          14          10 ]      0.01         1         8         2         1
 *      [           0     0.57143     0.42857 ]                                       -      0.02         1         8         0         1
 *      [           0         0.6         0.4 ]                                       -      0.01         1         8         0         1
 *      [     0.33333     0.66667           0 ]                                       -      0.02         1         5         0         1
 *      [     0.33333           0     0.66667 ] [          13          16          11 ]      0.01         1         5         3         1
 *      [     0.14286           0     0.85714 ]                                       -      0.02         1         6         0         1
 *      [     0.18033      0.2623     0.55738 ]                                       -      0.01         1         8         0         1
 *      [         0.6           0         0.4 ]                                       -      0.02         1         5         0         1
 *     extremal solutions
 *         OBJECTIV     cost2_1     cost3_1   solution file
 *            11.00       11.00       14.00   solutions/tenfelde-podehl-1.sol
 *            19.00       14.00       10.00   solutions/tenfelde-podehl-3.sol
 *            15.00        9.00       17.00   solutions/tenfelde-podehl-2.sol
 *            13.00       16.00       11.00   solutions/tenfelde-podehl-4.sol
 *
 *                                         Runs                               Solutions    Time/s     Nodes   LP Iter     V_new    V_proc
 *
 *                                           11                                       4      0.16        11        69        10        10
 *
 * As you can see, the program generates two tables of output. In the first table, each row represents a weighted optimization call.
 * The used weight is displayed as well as the cost vector of the optimal solution, if one was found.
 * Furthermore, the time (in seconds) as well as the number of lp iterations and nodes for solving the MIP is given.
 * The last two values pertain to the weight space polyhedron skeleton.
 * V_new is the number of vertices added to the skeleton due to a polyhedron update, if there was one.
 * V_proc is the number of vertices that have either been confirmed or removed during an update.
 *
 * A version of the lifted weight space algorithm was first mentioned in
 *
 * Matthias Ehrgott, Andreas Löhne, and Lizhen Shao:
 * A dual variant of Benson’s “outer approximation algorithm” for multiple objective linear programming,
 * _Journal of Global Optimization_, Volume 52, Issue 4, pp 757-778, 2012
 */
