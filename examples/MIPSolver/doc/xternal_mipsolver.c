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

/**@file   xternal_mipsolver.c
 * @brief  main documentation page of the MIP solver example
 * @author Timo Berthold
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page MIPSOLVER_MAIN SCIP as a MIP solver
 * @version  0.1
 * @author   Tobias Achterberg
 *
 * This very basic example illustrates how to integrate SCIP into your C++ source code.  The SCIP header files are
 * included from the directory which is denoted by SCIPDIR in the Makefile. If your SCIP headers are installed somewhere
 * else, just change this link. Since this example is written in C++, objscip/objscip.h is included at the beginning of
 * cppmain.cpp. It includes scip.h (which you would include when using SCIP as callable library in a C
 * program). The main function shows how the SCIP_RETCODE can be caught and handled. In runSCIP(), you see how a SCIP
 * instance is created and freed, plus a few more things. Here, the SCIP_RETCODEs are not checked explicitly but handled
 * by the SCIP_CALL macro. Please also note the CallableLibrary example for nonlinear problems.
 *
 * Installation
 * ------------
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 */

