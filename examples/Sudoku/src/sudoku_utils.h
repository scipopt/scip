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

/**@file   sudoku_utils.cpp
 * @brief  A set of utilities that are used to read the puzzle and display the puzzle
 * @author Naga V C Gudapati
 *
 *
 * A sudoku puzzle is in represented by a string of 81 charcaters. An already filled number in the puzzle is
 * represented by that number and if there is a blank, then it is represented by either '.' or '0' in the puzzle
 * string.
 * There are two functions:
 * 1) <tt>std::vector<std::vector<int>> get_sudoku_puzzle(std::string &file_path)</tt>
 *    which will read the textstring of sudoku puzzle into a 9x9 grid represented by a vector of vector of ints.
 *    The actual number is stores as itself and the blanks are stored as -1
 * 2) <tt>void print_sudoku(const std::vector<std::vector<int>> &sudoku_puzzle)</tt>
 *    which is a handy function to print the sudoku grid clearly.
 *
*/



#include <iostream>
#include <fstream>
#include <vector>
#include <string>


namespace sudoku
{
   std::vector<std::vector<int>> get_sudoku_puzzle( std::string &file_path )
   {
      /* setting up a 9x9 grid forstoring the sudoku puzzle. */
      std::vector<std::vector<int>> puzzle(9, std::vector<int>(9));

      /* Reading the puzzle into a stringstream */
      std::ifstream infile(file_path);

      std::string puzzle_data = "";

      if( infile.is_open() )
      {
         std::getline( infile, puzzle_data );
         if( puzzle_data.length() != 81 ) /* The puzzle should have 81 characters */
         {
            std::cerr << "Please check the puzzle file forinconsistencies"
                      << "\n";
            exit(1);
         }
      }

      int idx = 0; /* This variable will be used to access the numbers in the puzzle string */

      for( int i = 0; i < 9; ++i )
      {
         for( int j = 0; j < 9; ++j )
         {
            /* We will only convert the numeric string to an integer if it is not '.' or '0'. */
            if( (puzzle_data.substr(idx, 1) != ".") and (puzzle_data.substr(idx, 1) != "0") )
            {
               puzzle[i][j] = std::stoi( puzzle_data.substr(idx, 1) );
            }
            else
            {
               /* If we are currently reading a '.' or '0' make it -1. */
               puzzle[i][j] = -1;
            }
            idx++;
         }
      }

      return puzzle;
   }


   void print_sudoku( const std::vector<std::vector<int>> &sudoku_puzzle )
   {
      std::cout << "+----------+-----------+-----------+" << "\n";
      for( int i = 0; i < 9; ++i)
      {
         std::cout << "|";
         for( int j = 0; j < 9; ++j)
         {
            if( sudoku_puzzle[i][j] > 0)
            {

               if( j == 2 or j == 5 or j == 8 )
               {
                  std::cout << sudoku_puzzle[i][j] << " | ";
               }
               else
               {
                  std::cout << sudoku_puzzle[i][j] << "   ";
               }
            }
            else
            {
               if( j == 2 or j == 5 or j == 8 )
               {
                  std::cout << "*" << " | ";
               }
               else
               {
                  std::cout << "*" << "   ";
               }
            }
         }
         std::cout << "\n";

         if( i == 2 or i == 5 or i == 8 )
         {
            std::cout << "+----------+-----------+-----------+" << "\n";
         }
      }
   }
} // namespace sudoku
