/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal_callablelibrary.c
 * @brief  main document page
 * @author Stefan Vigerske
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page CALLABLELIBRARY_MAIN Callable Library Example
 * @author   Stefan Vigerske
 *
 * This example illustrates how to setup nonlinear constraints when using SCIP as callable library.
 * The example implements several small models that use specific types of constraints.
 *
 * - string.c shows how to setup general nonlinear and quadratic constraints
 * - gastrans.c shows how to setup absolute power constraints
 * - circle.c shows how to setup second-order-cone constraints
 * - brachistochrone.c shows how to setup nonlinear constraints
 * - circlepacking.c shows how to setup quadratic constraints
 */
