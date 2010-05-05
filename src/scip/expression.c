/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: expression.c,v 1.1 2010/05/05 17:53:12 bzfviger Exp $"

/**@file   expression.c
 * @brief  more methods for expressions and expression trees
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 *
 * This file contains methods for handling and manipulating expressions and expression trees
 * that are SCIP specific and thus not included in nlpi/expression.*
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "scip/expression.h"

#include "nlpi/struct_expression.h"

/** returns variables of expression tree */
SCIP_VAR** SCIPexprtreeGetVars(
   SCIP_EXPRTREE*        tree                /**< expression tree */
)
{
   assert(tree != NULL);

   return tree->vars;
}

/** stores array of variables in expression tree */
SCIP_RETCODE SCIPexprtreeSetVars(
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   int                   nvars,              /**< number of variables */
   SCIP_VAR**            vars                /**< variables */
)
{
   assert(tree != NULL);
   assert(vars != NULL || nvars == 0);

   if( tree->vars != NULL )
   {
      if( BMSreallocBlockMemoryArray(tree->blkmem, &tree->vars, tree->nvars, nvars) == NULL )
         return SCIP_NOMEMORY;
      BMScopyMemoryArray(tree->vars, vars, nvars);
   }
   else
   {
      if( BMSduplicateBlockMemoryArray(tree->blkmem, &tree->vars, vars, nvars) == NULL )
         return SCIP_NOMEMORY;
   }

   tree->nvars = nvars;

   return SCIP_OKAY;
}
