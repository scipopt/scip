/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal_relaxator.c
 * @brief  Main documentation page of the Relaxator example
 * @author Benjamin Mueller
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page RELAXATOR_MAIN Relaxator
 * @author   Benjamin Mueller
 *
 * This example illustrates how to write a relaxator for SCIP. It extends the default plugins of <a
 * href="http://scipopt.org">SCIP</a> by two additional relaxator plugins, one for solving the linear and another one
 * for solving the convex nonlinear relaxation. Both relaxators are called in each node of the branch-and-bound tree.
 *
 * Installation
 * ------------
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 */
