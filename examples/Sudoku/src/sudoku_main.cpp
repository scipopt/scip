/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sudoku_main.cpp
 * @brief  Sudoku solver built using constrained integer programming
 * @author Naga V C Gudapati


 *  To solve the given sudoku puzzle using integer programming, we shall use this handy tutorial;
 *  http://profs.sci.univr.it/~rrizzi/classes/PLS2015/sudoku/doc/497_Olszowy_Wiktor_Sudoku.pdf
 *
 *  An unsolved sudoku puzzle looks like below:<br>
 * <tt>
 *  +----------+-----------+-----------+ <br>
 *  |*   *   * | *   *   * | *   *   * | <br>
 *  |*   3   * | *   1   2 | *   *   8 | <br>
 *  |*   7   * | *   6   8 | *   2   * | <br>
 *  +----------+-----------+-----------+ <br>
 *  |*   *   * | *   *   9 | 8   7   * | <br>
 *  |1   2   * | 6   5   * | 4   *   * | <br>
 *  |*   *   * | *   *   * | *   *   6 | <br>
 *  +----------+-----------+-----------+ <br>
 *  |*   *   3 | 9   4   * | *   *   * | <br>
 *  |*   *   * | 2   *   * | *   6   * | <br>
 *  |4   *   * | *   *   * | *   3   1 | <br>
 *  +----------+-----------+-----------+ <br>
 * </tt>
 *
 *  The solved puzzle will have each of the nine rows and columns filled by numbers 1, ..., 9 each appearing
 *  exactly once. There are 9 subgrids and each subgrid also needs to be filled by numbers 1, ..., 9 each
 *  appearing exactly once.  As seen in the unsolved puzzles, some of the positions are already filled.
 *
 *  In this tutorial, we will learn how to model problems with many variables in SCIP, learn how to set parameters
 *  and learn how use the solution status to print custom output messages.
 *
*/



#include <iostream>
#include <sstream>

#include "sudoku_utils.h"

#include <scip/scip.h>
#include <scip/scipdefplugins.h>



int main( int args, char *argv[])
{
   if(args < 2)
   {
      std::cerr << "call " << argv[0] << " <puzzle file> " << "\n";
      exit(1);
   }

   std::string puzzle_file_path = argv[1];

   auto puzzle = sudoku::get_sudoku_puzzle(puzzle_file_path);

   std::cout << "The unsolved Sudoku Puzzle is: " << "\n";

   sudoku::print_sudoku(puzzle);

   /**
    * Setting up the SCIP environment
   */

   SCIP* scip = nullptr; /* Declaring the scip environment*/

   SCIP_CALL( SCIPcreate(&scip) ); /*Creating the SCIP environment */

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) ); /* include default plugins */


   SCIP_CALL( SCIPcreateProbBasic(scip, "sudoku") ); /* creating the SCIP Problem. */

   /**
    * The Sudoku puzzle is a feasibility problem and the objsense can be anything
   */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );


   /**
    * We have to define 9x9x9 variables. Let x_{ijk} where i = 1...9, j = 1...9 and k = 1..9 be those binary variables.
    * x_{ijk} is the the  binary variable representing the number k (1 or 2 or ... 9) in the ith row and jth column
   */

   std::vector<std::vector<std::vector< SCIP_VAR* >>> x_vars(9, std::vector<std::vector< SCIP_VAR* >>(9, std::vector< SCIP_VAR* >(9)));

   std::ostringstream namebuf;

   for( int i = 0; i < 9; ++i)
   {
      for( int j = 0; j < 9; ++j)
      {
         for( int k = 0; k < 9; ++k)
         {
            SCIP_VAR* var = nullptr;
            namebuf.str("");
            namebuf << "x[" << i << "," << j << "," << k << "]";
            SCIP_CALL( SCIPcreateVarBasic(
                                         scip,                  // SCIP environment
                                         &var,                  // reference to the variable
                                         namebuf.str().c_str(), // name of the variable
                                         0.0,                   // Lower bound of the variable
                                         1.0,                   // upper bound of the variable
                                         1.0,                   // Obj. coefficient.
                                         SCIP_VARTYPE_BINARY    // Binary variable
                                         ) );
            SCIP_CALL( SCIPaddVar(scip, var) );
            x_vars[i][j][k] = var;
         }
      }
   }

     /**
      * Since there is nothing advanced about our constraints,
      * we use the simple constraints provided by SCIPcreateConsBasicLinear()
     */


    /**
     * The column_constrs set of constraints will model that in each column, the numbers 1..9 should not repeat.
     * The constraint will look like x_1jk + x_2jk + x_3jk + ... + x_9jk = 1 for a given value of j and k
     */

   std::vector< SCIP_CONS* > column_constrs;

   for( int j = 0; j < 9; ++j)
   {
      for( int k = 0; k < 9; ++k)
      {
         SCIP_CONS* cons = nullptr;
         namebuf.str("");
         namebuf << "col_" << j << "_" << k << "]";
         /* This is an equality constraint and hence LHS and RHS are equal to 1 */
         SCIP_CALL( SCIPcreateConsBasicLinear(
                                                scip,
                                                &cons,                 // pointer to hold the created constraint
                                                namebuf.str().c_str(), // name of constraint
                                                0,                     // number of nonzeros in the constraint
                                                nullptr,               // array with variables of constraint entries
                                                nullptr,               // array with coefficients of constraint entries
                                                1.0,                   // left hand side of constraint
                                                1.0) );                // right hand side of constraint
         for( int i = 0; i < 9; ++i)
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, x_vars[i][j][k], 1.0));
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         column_constrs.push_back(cons);
      }
   }

   /**
    * The row_constrs set of constraints will model that in each row, the numbers 1..9 do not repeat.
    * The constraint will look like x_i1k + x_i2k + x_i3k + ... + x_i9k = 1 for a given value of i and k
    */
   std::vector< SCIP_CONS* > row_constrs;

   for( int i = 0; i < 9; ++i)
   {
      for( int k = 0; k < 9; ++k)
      {
         SCIP_CONS* cons = nullptr;

         namebuf.str("");
         namebuf << "row_" << i << "_" << k << "]";
          SCIP_CALL( SCIPcreateConsBasicLinear(
                                               scip,
                                               &cons,                 // pointer to hold the created constraint
                                               namebuf.str().c_str(), // name of constraint
                                               0,                     // number of nonzeros in the constraint
                                               nullptr,               // array with variables of constraint entries
                                               nullptr,               // array with coefficients of constraint entries
                                               1.0,                   // left hand side of constraint
                                               1.0) );                // right hand side of constraint
         for( int j = 0; j < 9; ++j)
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, x_vars[i][j][k], 1.0) );
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         row_constrs.push_back(cons);
      }
   }

   /**
    * The subgrid_constrs set of constraints model each of the 3x3 subgrids contains 1...9 without any repetition.
    */
   std::vector< SCIP_CONS* > subgrid_constrs; // These constraints will model

   for( int k = 0; k < 9; ++k)
   {
      for( int p = 0; p < 3; ++p)
      {
         for( int q = 0; q < 3; ++q)
         {
            SCIP_CONS* cons = nullptr;

            namebuf.str("");
            namebuf << "subgrid_" << k << "_" << p << "_" << q << "]";

            /* This is an equality constraint and hence LHS and RHS are equal to 1 */
            SCIP_CALL( SCIPcreateConsBasicLinear(
                                                scip,
                                                &cons,                 // pointer to hold the created constraint
                                                namebuf.str().c_str(), // name of constraint
                                                0,                     // number of nonzeros in the constraint
                                                nullptr,               // array with variables of constraint entries
                                                nullptr,               // array with coefficients of constraint entries
                                                1.0,                   // left hand side of constraint
                                                1.0) );                // right hand side of constraint

            /* since we are using 0 based indexing we should be careful with the loop indices. */
            for( int j = 3 * (p + 1) - 3; j < 3 * (p + 1); ++j)
            {
               for( int i = 3 * (q + 1) - 3; i < 3 * (q + 1); ++i)
               {
                  SCIP_CALL( SCIPaddCoefLinear(scip, cons, x_vars[i][j][k], 1.0) );
               }
            }
            SCIP_CALL( SCIPaddCons(scip, cons) );
            subgrid_constrs.push_back(cons);
         }
      }
   }


   /**
   * The fillgrid_constrs set of constraints ensure that every position in the whole 9x9 grid is
   * filled with exactly one number.
   */
   std::vector< SCIP_CONS* > fillgrid_constrs;

   for( int i = 0; i < 9; ++i )
   {
      for( int j = 0; j < 9; ++j )
      {
         SCIP_CONS* cons = nullptr;

         namebuf.str("");
         namebuf << "fillgrid_" << i << "_" << j << "]";
          SCIP_CALL( SCIPcreateConsBasicLinear(
                                               scip,
                                               &cons,                 // pointer to hold the created constraint
                                               namebuf.str().c_str(), // name of constraint
                                               0,                     // number of nonzeros in the constraint
                                               nullptr,               // array with variables of constraint entries
                                               nullptr,               // array with coefficients of constraint entries
                                               1.0,                   // left hand side of constraint
                                               1.0) );                // right hand side of constraint
         for( int k = 0; k < 9; ++k )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, cons, x_vars[i][j][k], 1.0) );
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         fillgrid_constrs.push_back(cons);
      }
   }


   /**
    * We use SCIPfixVar to fix the binary variables corresponding to the given value in the puzzle to 1.
    * see https://www.scipopt.org/doc-7.0.1/html/group__PublicVariableMethods.php#ga7965b16efcb2f8cdf7e289198c5cbe16
    * for more info
    */

   SCIP_Bool infeasible;
   SCIP_Bool fixed;
   for( int i = 0; i < 9; ++i )
   {
      for( int j = 0; j < 9; ++j )
      {
         /* The unsolved puzzle where there are blanks are represented by -1 in the puzzle datastructure */
         if( puzzle[i][j] > 0 )
         {
            SCIP_CALL( SCIPfixVar(scip, x_vars[i][j][(puzzle[i][j]) - 1], 1.0, &infeasible, &fixed) );
         }
      }
   }


   /**
    * In c++, we can set the SCIP parameters by <tt> sSCIPsetIntParam(SCIP* scip, "parameter", param_value) ) </tt>
   */

   SCIP_CALL( SCIPsetIntParam(scip, "display/verblevel", 0) ); // We use SCIPsetIntParams to turn off the logging.

   SCIP_CALL( SCIPsolve(scip) );

   /**
    * Some wrongly generated sudoku puzzles can be infeasible. So we use the solnstatus to display different types of output.
   */
   SCIP_STATUS soln_status = SCIPgetStatus( scip );

   if(soln_status == 11) // solution status of 11 indicates optimal solution was found. Hence we can print the final puzzle.
   {
      SCIP_SOL* sol;
      sol = SCIPgetBestSol( scip );

      for( int i = 0; i < 9; ++i )
      {
         for( int j = 0; j < 9; ++j )
         {
            for( int k = 0; k < 9; ++k )
            {
               if( SCIPgetSolVal(scip, sol, x_vars[i][j][k]) > 0 )
               {
                  /* As we are using 0 based indices, to display the final puzzle, we should increment values by 1. */
                  puzzle[i][j] = k + 1;
               }
            }
         }
      }
      std::cout << "The solved puzzle is: "
                << "\n";
      sudoku::print_sudoku( puzzle );
   }
   else if( soln_status == 12 ) // solutions status of 12 indicates that the puzzle is infeasible.
   {
      std::cout << "Check the Input puzzle"
                << "\n";
   }
   else
   {
      std::cerr << "Something went wrong during the optimization." << "\n";
      exit(1);
   }

   //Freeing the variables

   for( int i = 0; i < 9; ++i )
   {
      for( int j = 0; j < 9; ++j )
      {
         for( int k = 0; k < 9; ++k )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &x_vars[i][j][k]) );
         }
      }
   }
   x_vars.clear();

   // Freeing the constraints
   /**
    * We use C++11's auto keyword to iterate over the constraints and free them.
   */

   for( auto &constr : column_constrs)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &constr) );
   }
   column_constrs.clear();

   for( auto &constr : row_constrs)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &constr) );
   }
   row_constrs.clear();

   for( auto &constr : subgrid_constrs)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &constr) );
   }
   subgrid_constrs.clear();

   for( auto &constr : fillgrid_constrs)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &constr) );
   }
   fillgrid_constrs.clear();

   //freeing scip

   SCIP_CALL( SCIPfree(&scip) );
}
