/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   xternal_pbsolver.c
 * @brief  Pseudo-Boolean solver application
 * @author Alexander Hoen
 * @author Gioni Mexi
 * @author Dominik Kamp
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@page PBSOLVER_MAIN PBSolver
 * @author Alexander Hoen
 * @author Gioni Mexi
 * @author Dominik Kamp
 *
 *  A solver for pseudoboolean problems in OPB or WBO format. It complies by default with the technical regulations of
 *  the PB competition. Therefore, the following plugins are implemented:
 *
 * - a \ref message_pb.c "message handler" to produce a valid general log with accepted line initials
 * - a \ref event_bestsol.c "event handler" to signal achievements of best primal solutions (only for optimization)
 *
 **@todo Add problem specific parameter settings.
 *
 * Installation
 * ------------
 *
 * See the @ref INSTALL_APPLICATIONS_EXAMPLES "Install file"
 */
