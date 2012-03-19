/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
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
 * @author Timo Berthold
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage SCIP TSP example
 * @version  0.9
 * @author   Timo Berthold

 * This is an example of using SCIP to solve the TSP problem on undirected graphs.
 * Here is the CIP model that we use:
 *
 * Given: a graph G=(V,E) with edge weights c_e
 * Task: find hamiltonian cycle T in G with minimal length c(T)
 *
 * Variables:  \f$ x_e \in \{0,1\} \, \forall e \in E, x_e = 1 \Leftrightarrow e \in T \f$
 *
 * Constraints:
 * 1. \f$\sum_{e \in \delta(v)} x_e = 2 \, \forall v \in V \qquad \qquad \f$
 * 2. subtour(G,x)
 *
 * Semantics of constraints:
 * 1. usual linear constraints
 * 2. subtour(G,x) \f$\Leftrightarrow\f$ T defined by x does not contain any cycle of length < |V|
 *
 * A few remarks to the model and the implementation (references to code lines might
 * not be up to date):
 *
 * As one can see, the TSP-Model consists of |V| linear constraints (the
 * degree constraints) and one "subtour" constraint. The latter is a
 * complex, non-linear constraint for which one has to implement an own
 * constraint handler.
 * The variables are created in the TSP file reader ReaderTSP.cpp at line
 * 379. A pointer to each variable is stored in the data structure of the
 * corresponding edge (i.e., in edge->var and edge->back->var, since the
 * formally undirected graph is represented as a directed graph with
 * antiparallel arcs).
 * The degree constraints are created at line 397. The data for the
 * linear degree constraints are the coefficients (for each \f$e \in
 * \delta(v)\f$ the variable x_e has coefficient 1.0) which are generated at
 * line 405, and the left and right hand sides of the inequality, which
 * are both set to 2.0 at line 397, such that the linear constraint
 * becomes an equation.
 * The subtour constraint is created at line 417. The data for this
 * constraint is the graph and the variables (see above), but we only have
 * to store a pointer to the graph because the edges already have links to
 * their respective variables.
 *
 * Now the problem instance is defined, and the "only" thing left is to
 * implement the semantics of the "subtour" constraint. This is of
 * course done in the subtour constraint handler ConshdlrSubtour.cpp. The
 * main task of a constraint handler is to decide whether a given
 * solution is feasible for all constraints of the constraint handler's
 * type (i.e., for example, for all linear constraint in the instance,
 * for all knapsack constraints in the instance, for all subtour
 * constraints in the instance, ...). This is performed in the
 * scip_enfolp(), scip_enfops(), and scip_check() methods. To implement
 * these three methods and the scip_lock() method (called the
 * "fundamental callback methods" in the SCIP documentation) is already
 * enough to obtain a correct algorithm, which means that the solver will
 * find (after waiting long enough) the optimal solution. The remaining
 * methods are only needed to speed up the solving process (for example,
 * cutting plane separation and domain propagation).
 *
 * As there is only one subtour constraint in a TSP instance, all the
 * loops "for( int i = 0; i < nconss; ++i )" in ConshdlrSubtour.cpp are a
 * bit ridiculous, since nconss will always be equal to one. However,
 * nothing prevents a user from using the subtour constraint handler in a
 * different application where you have more than one graph and a
 * solution must not contain any subtour in each of the graphs.
 *
 * Additionally, this example contains two wellknown combinatorial heuristics for the TSP,
 * namely a farthest insert heuristic (see HeurFarthestInsert.cpp) and a 2-opt heuristic (see Heur2opt.cpp)
 * and a rounding heuristic (see HeurFrats.cpp) for TSPs.
 * The idea of the latter one is to take an LP solution, and construct a hamiltonian cycle in the following way:
 * If x_e is equal to one, add edge e to the cycle.
 * Iterate over the remaining variables in nonincreasing order of their LP value x_e and add the corresponding edge e, 
 * if it does not close a subtour.
 */


