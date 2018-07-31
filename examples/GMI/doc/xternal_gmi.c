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

/**@file   xternal_gmi.c
 * @brief  main document page
 * @author Marc Pfetsch
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page GMI_MAIN Gomory Mixed Integer Cut Example
 * @version  1.0
 * @author   Giacomo Nannicini
 * @author   Marc Pfetsch
 *
 *
 * This example provides a textbook implementation of Gomory mixed integer (GMI) cuts.
 *
 * The default implementation in SCIP does not produce GMI cuts in the strict sense, since it applies the CMIR function
 * to the aggregated row. This function can, among other things, take variable bounds into account. Thus, the resulting
 * cuts cannot be used for comparison with standard GMI cuts. This example remedies this situation.
 *
 * The implementation has been used in the paper
 *
 * G. Cornuejols, F. Margot and G. Nannicini:@n
 * On the safety of Gomory cut generators.@n
 * Math. Program. Comput. 5(4), 2013.
 */
