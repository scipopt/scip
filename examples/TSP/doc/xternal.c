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
 * @brief  main document page
 * @author Timo Berthold
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage Overview
 * @version  0.9
 * @author   Timo Berthold

 * This is an example of using SCIP to solve the TSP problem on undirected graphs.
 * Here is the CIP model that we use:
 *
 * Given: a graph \f$ G=(V,E) \f$  with edge weights \f$ c_e \f$
 * Task: find hamiltonian cycle \f$ T \f$  in \f$ G \f$  with minimal length \f$ c(T) \f$
 *
 * Variables:  \f$ x_e \in \{0,1\} \, \forall e \in E, x_e = 1 \Leftrightarrow e \in T \f$
 *
 * Constraints:
 *  -# \f$\sum_{e \in \delta(v)} x_e = 2 \, \forall v \in V \qquad \qquad \f$
 *  -# subtour\f$(G,x) \f$
 *
 * Semantics of constraints:
 *  -# usual linear constraints
 *  -# subtour\f$(G,x)\Leftrightarrow\f$ \f$ T \f$  defined by \f$ x \f$  does not contain
 * any cycle of length \f$ < |V| \f$.
 *
 * A few remarks to the model and the implementation (references to code lines might
 * not be up to date):
 *
 * As one can see, the TSP-Model consists of \f$ |V| \f$  linear constraints (the
 * degree constraints) and one "subtour" constraint. The latter is a
 * complex, non-linear constraint for which one has to implement an own
 * constraint handler.
 * The variables are created in the TSP file reader ReaderTSP.cpp at line
 * 408:
 * \code
 *    SCIP_CALL( SCIPcreateVar(scip, &var, varname.str().c_str(), 0.0, 1.0, edge->length,
 *           SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
 *
 *    SCIP_CALL( SCIPaddVar(scip, var) );
 *    SCIP_CALL( addVarToEdges(scip, edge, var) );
 * \endcode
 *  A pointer to each variable is stored in the data structure of the
 * corresponding edge (i.e., in <code>edge->var</code> and <code>edge->back->var</code>,
 *  since the formally undirected graph is represented as a directed graph with
 * antiparallel arcs).
 * After that, the degree constraints are created at line 430.
 * \code
 * SCIP_CALL( SCIPcreateConsLinear(scip, &cons, consname.str().c_str(), 0, NULL, NULL, 2.0, 2.0,
 *    TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
 * \endcode
 * The data for the
 * linear degree constraints are the coefficients (for each \f$e \in
 * \delta(v)\f$ the variable \f$ x_e \f$  has coefficient 1.0) which are generated at
 * line 437, and the left and right hand sides of the inequality, which
 * are both set to 2.0 at line 430, such that the linear constraint
 * becomes an equation.
 * The subtour constraint is created at line 449. The data for this
 * constraint is the graph and the variables (see above), but we only have
 * to store a pointer to the graph because the edges already have links to
 * their respective variables.
 * \code
 *   SCIP_CALL( SCIPcreateConsSubtour(scip, &cons, "subtour", graph,
 *       FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE ) );
 * \endcode
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
 * loops
 * \code
 * for( int i = 0; i < nconss; ++i )
 * ...
 * \endcode in ConshdlrSubtour.cpp are a
 * bit ridiculous, since <code>nconss</code> will always be equal to one. However,
 * nothing prevents a user from using the subtour constraint handler in a
 * different application where you have more than one graph and a
 * solution must not contain any subtour in each of the graphs.
 *
 * Additionally, this example contains two well known combinatorial heuristics for the TSP,
 * namely a \ref HeurFarthestInsert.cpp "farthest insert heuristic"  and a \ref Heur2opt.cpp "2-opt heuristic",
 * and a \ref HeurFrats.cpp "rounding heuristic" for TSPs.
 * The idea of the latter one is to take an LP solution and construct a hamiltonian cycle in the following way:
 * If \f$ x_e \f$  is equal to one, add edge \f$ e \f$  to the cycle.
 * Iterate over the remaining variables in nonincreasing order of their LP value \f$ x_e \f$  and add
 * the corresponding edge \f$ e \f$ ,
 * if it does not close a subtour.
 */


