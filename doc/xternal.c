/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SIP; see the file COPYING. If not email to scip@zib.de.       */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* $Id: xternal.c,v 1.1 2002/10/23 14:31:37 bzfpfend Exp $
/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage SCIP (Solving Constraint Integer Programs)
   @version  0.0.1
   @author   Tobias Achterberg
   @author   Thorsten Koch
   @author   Alexander Martin

   SCIP is a program and library to solve constraint mixed integer programs (CMIPs).

   - \ref CODE  "Coding style guidlines"
   - \ref SEPAR "How to add separators"
   - \ref HEUR  "How to add heuristics"
*/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CODE Coding style guidlines
 *
 * Here we explain how your code should look like, so it is easy for us to include it into the SCIP distribution.
 *
 * - Indentation is 3 spaces. No <tabs> anywhere in the code.
 * - Allways only one declaration in a line.
 * - Braces are on a new line and not indented.
 * - Spaces around all operators.
 * - Use assert() to show preconditions for the parameters, invariants and postconditions.
 * - All global functions start with "SCIP". In the usual naming scheme this is followed by a verb and a nown, like in
 *   SCIPaddSepar(). Functions return TRUE or FALSE should be named like SCIPisSeparEnabled().
 * - Make all functions that are not used outside the module 'static'. Naming should start with a lower case letter.
 * - Variable names should start with a lower case letter.
 * - For each structure there is a typedef with the name in all upper case.
 * - Defines should be named all upper case.
 * - Document functions, parameters and variables doxygen conform.
 *
 * As an example have a look at separate.c .
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page SEPAR How to add separators
 *
 * Here we explain how to add separation routines to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page HEUR How to add heuristics
 *
 * Here we explain how to add a new primal heuristic to SCIP.
 */




