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

/**@file   examples/Sudoku/doc/xternal_sudoku.c
 * @brief  main document page
 * @author Naga V C Gudapati
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page     SUDOKU_MAIN Sudoku Solver
 * @author Naga V C Gudapati
 *
 * Methodology
 * ============
 *
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
 *  In this example, we see
 *  -  how to model problems with many variables in SCIP,
 *  -  how to set parameters
 *  -  how use the solution status to print custom output messages.
 *
 * Data Format
 * ============
 *
 * A sudoku puzzle is in represented by a string of 81 charcaters. An already filled number in the
 * puzzle is represented by that number; a blank is represented by either '.' or '0' in the puzzle
 * string.
 * The input file is containing the configuration of a sudoku read rowwise as a string.
 * For example, the above puzzle is represented by
 * 000000000030012008070068020000009870120650400000000006003940000000200060400000031
 * or
 * ..........3..12..8.7..68.2......987.12.65.4..........6..394.......2...6.4......31
 *
 * Installation
 * ============
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 *
 */
