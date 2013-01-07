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

/**@mainpage example project using SCIP as a MIP solver
 * @version  0.1
 * @author   Tobias Achterberg
 *
 * This very basic example illustrates how to integrate SCIP into your C++ source code.  The SCIP header files are
 * included from the directory which is denoted by SCIPDIR in the Makefile. If your SCIP headers are installed somewhere
 * else, just change this link. Since this example is written in C++, objscip/objscip.h is included at the beginning of
 * cppmain.cpp. It includes scip/scip.h (which you would include when using SCIP as callable library in a C
 * program). The main function shows how the SCIP_RETCODE can be caught and handled. In runSCIP(), you see how a SCIP
 * instance is created and freed, plus a few more things. Here, the SCIP_RETCODEs are not checked explicitly but handled
 * by the SCIP_CALL macro. Please also note the CallableLibrary example for nonlinear problems.
 */

