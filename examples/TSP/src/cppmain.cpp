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
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cppmain.cpp
 * @brief  main file for C++ TSP example using SCIP as a callable library
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 * This is an example of using SCIP to solve the TSP problem on undirected graphs.
 * Here is the CIP model that we use:
 *
 * Given: a graph G=(V,E) with edge weights c_e
 * Task: find hamiltonian cycle T subseteq E in G with minimal length c(T)
 *
 * Variables: x_e in {0,1} for all e in E, x_e = 1 <=> e in T
 * Constraints:
 * 1. sum_{e in delta(v)} x_e == 2 for all v in V
 * 2. subtour(G,x)
 *
 * Semantics of constraints:
 * 1. usual linear constraints
 * 2. subtour(G,x) <=> T defined by x does not contain any cycle of length < |V|
 *
 * A few remarks to the model and the implementation (references to code lines might
 * not be up to date):
 *
 * As one can see, the TSP-Model consists of |V| linear constraints (the
 * degree constraints) and one "subtour" constraint. The latter is a
 * complex, non-linear constraint for which one has to implement an own
 * constraint handler.
 * The variables are created in the TSP file reader ReaderTsp.cpp at line
 * 311. A pointer to each variable is stored in the data structure of the
 * corresponding edge (i.e., in edge->var and edge->back->var, since the
 * formally undirected graph is represented as a directed graph with
 * antiparallel arcs).
 * The degree constraints are created at line 330. The data for the
 * linear degree constraints are the coefficients (for each e in
 * delta(v) the variable x_e has coefficient 1.0) which are generated at
 * line 337, and the left and right hand sides of the inequality, which
 * are both set to 2.0 at line 330, such that the linear constraint
 * becomes an equation.
 * The subtour constraint is created at line 349. The data for this
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
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <iostream>

/* include SCIP components */
#include "objscip/objscip.h"
#include "objscip/objscipdefplugins.h"

/* include TSP specific components */
#include "ReaderTSP.h"
#include "ConshdlrSubtour.h"
#include "HeurFarthestInsert.h"
#include "Heur2opt.h"
#include "HeurFrats.h"
#include "EventhdlrNewSol.h"

using namespace scip;
using namespace tsp;
using namespace std;

/** creates and runs a SCIP instance with default and TSP plugins */
static
SCIP_RETCODE runSCIP(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{
   SCIP* scip = NULL;


   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include TSP specific plugins */
   SCIP_CALL( SCIPincludeObjReader(scip, new ReaderTSP(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjConshdlr(scip, new ConshdlrSubtour(scip), TRUE) ); 
   SCIP_CALL( SCIPincludeObjEventhdlr(scip, new EventhdlrNewSol(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjHeur(scip, new HeurFarthestInsert(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjHeur(scip, new Heur2opt(scip), TRUE) );
   SCIP_CALL( SCIPincludeObjHeur(scip, new HeurFrats(scip), TRUE) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );


   /**********************************
    * Process command line arguments *
    **********************************/

   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, "sciptsp.set") );


   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method starting TSP code */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode;

   retcode = runSCIP(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
