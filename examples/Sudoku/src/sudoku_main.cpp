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

/**@file   sudoku_main.cpp
 * @brief  Sudoku solver built using constrained integer programming
 * @author Naga V C Gudapati
*/

#include <iostream>
#include <sstream>

#include "sudoku_utils.h"

#include <scip/scip.h>
#include <scip/scipdefplugins.h>


int main(int argc, char *argv[])
{
   /* declaration of variables used in forloops later on */
   int i = 0;
   int j = 0;
   int k = 0;
   int p = 0;
   int q = 0;

   if( argc < 2 )
   {
      std::cerr << "call " << argv[0] << " <puzzle file> " << "\n";
      exit(1);
   }

   std::string puzzlefilepath = argv[1];
   auto puzzle = sudoku::getSudokuPuzzle(puzzlefilepath);
   std::cout << "The unsolved Sudoku Puzzle is: " << "\n";
   sudoku::printSudoku(puzzle);

   /*
    * Setting up the SCIP environment
   */
   SCIP* scip = nullptr; /* Declaring the scip environment*/
   SCIP_CALL( SCIPcreate(&scip) ); /*Creating the SCIP environment */
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) ); /* include default plugins */
   SCIP_CALL( SCIPcreateProbBasic(scip, "sudoku") ); /* creating the SCIP Problem. */

   /*
    * The Sudoku puzzle is a feasibility problem and the objsense can be anything
   */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /*
    * We have to define 9x9x9 variables. Let x_{ijk} where i = 1...9, j = 1...9 and k = 1..9 be those binary variables.
    * x_{ijk} is the the  binary variable representing the number k (1 or 2 or ... 9) in the ith row and jth column
   */
   std::vector<std::vector<std::vector< SCIP_VAR* >>> xvars(9, std::vector<std::vector< SCIP_VAR* >>(9, std::vector< SCIP_VAR* >(9)));
   std::ostringstream namebuf;

   for( i = 0; i < 9; ++i )
   {
      for( j = 0; j < 9; ++j )
      {
         for( k = 0; k < 9; ++k )
         {
            SCIP_VAR* var = nullptr;
            namebuf.str("");
            namebuf << "x[" << i << "," << j << "," << k << "]";
            SCIP_CALL( SCIPcreateVarBasic(
                                           scip,                  /* SCIP environment */
                                           &var,                  /* reference to the variable */
                                           namebuf.str().c_str(), /* name of the variable */
                                           0.0,                   /* lower bound of the variable */
                                           1.0,                   /* upper bound of the variable */
                                           1.0,                   /* obj. coefficient. */
                                           SCIP_VARTYPE_BINARY    /* variable is binary */
                                           ) );
            SCIP_CALL( SCIPaddVar(scip, var) );
            xvars[i][j][k] = var;
         }
      }
   }

   /* for each column j and each number k we must have that only one entry of column j is k; since x_{ijk} is 1
    * if and only if the i-th entry of the j-th column is k, we can model this requirement with the following linear
    * constraint:
    *             x_1jk + x_2jk + x_3jk + ... + x_9jk = 1 for each j and k
    */
   std::vector< SCIP_CONS* > columnconstraints;
   for( j = 0; j < 9; ++j )
   {
      for( k = 0; k < 9; ++k )
      {
         SCIP_CONS* cons = nullptr;
         namebuf.str("");
         namebuf << "col_" << j << "_" << k << "]";

         /* we first create an empty equality constraint (i.e. the lhs and rhs are equal to 1) and then add the
          * variables
          */
         SCIP_CALL( SCIPcreateConsBasicLinear(
                                               scip,
                                               &cons,                 /* pointer to hold the created constraint */
                                               namebuf.str().c_str(), /* name of constraint */
                                               0,                     /* number of nonzeros in the constraint */
                                               nullptr,               /* array with variables of constraint entries */
                                               nullptr,               /* array with coefficients of constraint entries */
                                               1.0,                   /* left hand side of constraint */
                                               1.0) );                /* right hand side of constraint */
         for( i = 0; i < 9; ++i )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, xvars[i][j][k], 1.0) );
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         columnconstraints.push_back(cons);
      }
   }

   /* for each row i and each number k we must have that only one entry of row i is k; we can model
    * this requirement with the following linear constraint:
    *             x_i1k + x_i2k + x_i3k + ... + x_i9k = 1 for each i and k
    */
   std::vector< SCIP_CONS* > rowconstraints;
   for( i = 0; i < 9; ++i )
   {
      for( k = 0; k < 9; ++k )
      {
         SCIP_CONS* cons = nullptr;

         namebuf.str("");
         namebuf << "row_" << i << "_" << k << "]";

         /* we first create an empty equality constraint (i.e. the lhs and rhs are equal to 1) and then add the
          * variables
          */
          SCIP_CALL( SCIPcreateConsBasicLinear(
                                                scip,
                                                &cons,                 /* pointer to hold the created constraint */
                                                namebuf.str().c_str(), /* name of constraint */
                                                0,                     /* number of nonzeros in the constraint */
                                                nullptr,               /* array with variables of constraint entries */
                                                nullptr,               /* array with coefficients of constraint entries */
                                                1.0,                   /* left hand side of constraint */
                                                1.0) );                /* right hand side of constraint */
         for( j = 0; j < 9; ++j )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, xvars[i][j][k], 1.0) );
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         rowconstraints.push_back(cons);
      }
   }

   /* for each 3x3 subgrid we must we must have that only one entry is k;
    * a subgrid is formed by the entries (i,j) such that i is in {p, p + 1, p + 2} and j is in {q, q + 1, q + 2} for
    * each (p, q) in {(1,1), (1,4), (1,7), (4,1), (4, 4), (4, 7), (7, 1), (7, 4), (7, 7)}
    * for the (p,q)-th subgrid we can model this requirement with the following linear constraint:
    *             sum_{i in [p, p + 2], j in [q, q + 2]} x_ijk = 1 for each k
    */
   std::vector< SCIP_CONS* > subgridconstraints;
   for( k = 0; k < 9; ++k )
   {
      for( p = 0; p < 3; ++p )
      {
         for( q = 0; q < 3; ++q )
         {
            SCIP_CONS* cons = nullptr;

            namebuf.str("");
            namebuf << "subgrid_" << k << "_" << p << "_" << q << "]";

            /* we first create an empty an equality constraint (i.e. the lhs and rhs are equal to 1) and then add the
             * variables
             */
            SCIP_CALL( SCIPcreateConsBasicLinear(
                                                  scip,
                                                  &cons,                 /* pointer to hold the created constraint */
                                                  namebuf.str().c_str(), /* name of constraint */
                                                  0,                     /* number of nonzeros in the constraint */
                                                  nullptr,               /* array with variables of constraint entries */
                                                  nullptr,               /* array with coefficients of constraint entries */
                                                  1.0,                   /* left hand side of constraint */
                                                  1.0) );                /* right hand side of constraint */

            /* since we are using 0 based indexing we should be careful with the loop indices. */
            for( j = 3 * (p + 1) - 3; j < 3 * (p + 1); ++j )
            {
               for( i = 3 * (q + 1) - 3; i < 3 * (q + 1); ++i )
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, cons, xvars[i][j][k], 1.0) );
               }
            }
            SCIP_CALL( SCIPaddCons(scip, cons) );
            subgridconstraints.push_back(cons);
         }
      }
   }


   /* so far we have required that for each column, row, and subgrid only one occurence of {1, ..., 9} appears.
    * However, we have not yet imposed that they must appear in different positions. So, for example, a solution
    * that says that the numbers 3 and 4 appear in entry (1,1) could be feasible.
    * However, that can only happen if some entry (i,j) has no number assigned to it.
    * Thus, by requiring that each entry (i,j) has a number assigned to it we obtain a correct model.
    * This requirement with the following linear constraint:
    *             x_ij1 + x_ij2k + x_ij3 + ... + x_ij9 = 1 for each i and j
    *
    */
   std::vector< SCIP_CONS* > fillgridconstraints;
   for( i = 0; i < 9; ++i )
   {
      for( j = 0; j < 9; ++j )
      {
         SCIP_CONS* cons = nullptr;

         namebuf.str("");
         namebuf << "fillgrid_" << i << "_" << j << "]";

         /* we first create an empty an equality constraint (i.e. the lhs and rhs are equal to 1) and then add the
          * variables
          */
         SCIP_CALL( SCIPcreateConsBasicLinear(
                                               scip,
                                               &cons,                 /* pointer to hold the created constraint */
                                               namebuf.str().c_str(), /* name of constraint */
                                               0,                     /* number of nonzeros in the constraint */
                                               nullptr,               /* array with variables of constraint entries */
                                               nullptr,               /* array with coefficients of constraint entries */
                                               1.0,                   /* left hand side of constraint */
                                               1.0) );                /* right hand side of constraint */

         for( k = 0; k < 9; ++k )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, xvars[i][j][k], 1.0) );
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         fillgridconstraints.push_back(cons);
      }
   }


   /* we use SCIPfixVar to fix the binary variables corresponding to the given value in the puzzle to 1.
    * see https://www.scipopt.org/doc-7.0.1/html/group__PublicVariableMethods.php#ga7965b16efcb2f8cdf7e289198c5cbe16
    * for more info
    */
   for( i = 0; i < 9; ++i )
   {
      for( j = 0; j < 9; ++j )
      {
         /* The unsolved puzzle where there are blanks are represented by -1 in the puzzle datastructure */
         if( puzzle[i][j] > 0 )
         {
            SCIP_Bool infeasible;
            SCIP_Bool fixed;

            SCIP_CALL( SCIPfixVar(scip, xvars[i][j][puzzle[i][j] - 1], 1.0, &infeasible, &fixed) );

            assert(fixed == TRUE);
            /* we are assuming that the puzzle is not instantly infeasible */
            assert(infeasible == FALSE);
         }
      }
   }


   /* in c++, we can set the SCIP parameters by <tt> sSCIPsetIntParam(SCIP* scip, "parameter", param_value) ) </tt> */
   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) ); /* turns off the SCIP output */
   SCIP_CALL( SCIPsolve(scip) );

   /* Some wrongly generated sudoku puzzles can be infeasible. So we use the solnstatus to display different types of
    * output.
   */
   SCIP_STATUS solutionstatus = SCIPgetStatus( scip );

   if( solutionstatus == SCIP_STATUS_OPTIMAL )
   {
      /* SCIP_STATUS_OPTIMAL status indicates that an optimal solution was found, hence we can print the final puzzle */
      SCIP_SOL* sol;
      sol = SCIPgetBestSol( scip );

      for( i = 0; i < 9; ++i )
      {
         for( j = 0; j < 9; ++j )
         {
            for( k = 0; k < 9; ++k )
            {
               if( SCIPgetSolVal(scip, sol, xvars[i][j][k]) > 0 )
               {
                  /* As we are using 0 based indices, to display the final puzzle, we should increment values by 1. */
                  puzzle[i][j] = k + 1;
               }
            }
         }
      }
      std::cout << "The solved puzzle is: "
                << "\n";
      sudoku::printSudoku( puzzle );
   }
   else if( solutionstatus == SCIP_STATUS_INFEASIBLE )
   {
      /* solutions status of SCIP_STATUS_INFEASIBLE indicates that the puzzle is infeasible. */
      std::cout << "Check the Input puzzle"
                << "\n";
   }
   else
   {
      std::cerr << "Something went wrong during the optimization." << "\n";
      exit(1);
   }

   /*freeing the variables */
   for( i = 0; i < 9; ++i )
   {
      for( j = 0; j < 9; ++j )
      {
         for( k = 0; k < 9; ++k )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &xvars[i][j][k]) );
         }
      }
   }
   xvars.clear();

   /* freeing the constraints
    * we use C++11's auto keyword to iterate over the constraints and free them.
    */
   for( auto &constr : columnconstraints )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &constr) );
   }
   columnconstraints.clear();

   for( auto &constr : rowconstraints )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &constr) );
   }
   rowconstraints.clear();

   for( auto &constr : subgridconstraints )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &constr) );
   }
   subgridconstraints.clear();

   for( auto &constr : fillgridconstraints )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &constr) );
   }
   fillgridconstraints.clear();

   /*freeing scip */

   SCIP_CALL( SCIPfree(&scip) );
}
