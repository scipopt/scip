/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   xternal.c
 * @brief  main document page
 * @author Tobias Achterberg
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/**@mainpage SCIP (Solving Constraint Integer Programs)
 * @version  0.79a
 * @author   Tobias Achterberg
 * @author   Thorsten Koch
 * @author   Alexander Martin
 * @author   Kati Wolter
 *
 * SCIP is a program and library to solve constraint mixed integer programs (CMIPs).
 *
 * - \ref CODE  "Coding style guidlines"
 * - \ref CONS  "How to add constraint handlers"
 * - \ref SEPAR "How to add separators"
 * - \ref HEUR  "How to add heuristics"
 * - \ref OBJ   "Creating, capturing, releasing, and adding data objects"
 * - \ref PARAM "Adding additional user parameters"
*/

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CODE Coding style guidlines
 *
 * Here we explain how your code should look like, so it is easy for us to include it into the SCIP distribution.
 *
 * - Indentation is 3 spaces. No tabs anywhere in the code.
 * - Allways only one declaration in a line.
 * - Braces are on a new line and not indented.
 * - Spaces around all operators.
 * - Use assert() to show preconditions for the parameters, invariants and postconditions.
 * - All global functions start with "SCIP". In the usual naming scheme this is followed by the object and a method name
 *   like in SCIPlpAddRow(). Functions return TRUE or FALSE should be named like SCIPlpiIsOptimal().
 * - Make all functions that are not used outside the module 'static'. Naming should start with a lower case letter.
 * - Variable names should be all lower case.
 * - For each structure there is a typedef with the name in all upper case.
 * - Defines should be named all upper case.
 * - Document functions, parameters and variables doxygen conform.
 *
 * As an example have a look at tree.c .
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page CONS How to add constraint handlers
 *
 * This page is not yet written. Here we will explain how to add constraint handlers to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page SEPAR How to add separators
 *
 * This page is not yet written. Here we will explain how to add separation routines to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page HEUR How to add heuristics
 *
 * This page is not yet written. Here we will explain how to add a new primal heuristic to SCIP.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page OBJ Creating, capturing, releasing, and adding data objects
 *
 *  Data objects (variables, constraints, rows) are subject to reference counting
 *  to avoid expensive copying operations. Creating such an object will set the
 *  reference count to one. Capturing an object increases the reference counter,
 *  releasing it decreases the counter. If the reference counter gets zero, the
 *  object is destroyed.
 *
 *  Remember that a created data object is automatically captured. If the user
 *  doesn't need the object anymore, he has to call the object's release() method.
 *
 *  When a data object is added to SCIP, it is captured again, such that a
 *  release() call does not destroy the object. If SCIP doesn't need the object
 *  anymore, it is automatically relased.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/**@page PARAM Adding additional user parameters
 *
 *  The user may add own parameters to SCIP with a call to SCIPaddXxxParam(). Using
 *  these methods, he has two possibilities where to store the actual parameter value:
 *   - If the given valueptr is NULL, SCIP stores the parameter value internally, and
 *     the user can only access the value with the SCIPgetXxxParam() and
 *     SCIPsetXxxParam() calls.
 *   - If the given valueptr is not NULL, SCIP stores the parameter value at the given
 *     address, and the user can directly manipulate the value at this address.
 *     He has to be careful with memory management in string parameters: when the
 *     SCIPaddStringParam() method is called, the given address must hold a char*
 *     pointer with value NULL. The default value is then copied into this pointer,
 *     allocating memory with malloc(). If the parameter is changed, the old string
 *     is free()'d and the new one is copied to a new memory area allocated with
 *     malloc(). When the parameter is freed, the memory is freed with free().
 *     The user should not interfere with this internal memory management. Accessing
 *     the string parameter through the given valueptr is okay as long it does not
 *     involve reallocating memory for the string.
 */

