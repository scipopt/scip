/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal_callablelibrary.c
 * @brief  main document page
 * @author Stefan Vigerske
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page CALLABLELIBRARY_MAIN Callable Library
 * @author   Stefan Vigerske
 *
 * This example implements several small models to illustrate how to set up different types of constraints when using SCIP as callable library.
 * At the moment, only examples that show how to set up nonlinear constraints are available:
 *
 * - \ref brachistochrone.c
 * - \ref circlepacking.c
 * - \ref gastrans.c
 * - \ref spring.c
 *
 * Installation
 * ------------
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 */
