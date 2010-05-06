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
#pragma ident "@(#) $Id: expression.c,v 1.3 2010/05/06 18:30:23 bzfviger Exp $"

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
#include "scip/scip.h"

#include "nlpi/struct_expression.h"

/** translate from one value of infinity to another
 *
 *  if val is >= infty1, then give infty2, else give val
 */
#define infty2infty(infty1, infty2, val) (val >= infty1 ? infty2 : val)

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

/** evaluates an expression tree for a primal solution or LP solution */
SCIP_RETCODE SCIPexprtreeEvalSol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_SOL*             sol,                /**< a solution, or NULL for current LP solution */
   SCIP_Real*            val                 /**< buffer to store value */
)
{
   SCIP_Real* varvals;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   if( tree->nvars == 0 )
   {
      SCIP_CALL( SCIPexprEval(tree->root, NULL, tree->params, val) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, tree->nvars) );
   SCIP_CALL( SCIPgetSolVals(scip, sol, tree->nvars, tree->vars, varvals) );

   SCIP_CALL( SCIPexprEval(tree->root, varvals, tree->params, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. current global bounds */
SCIP_RETCODE SCIPexprtreeEvalIntGlobalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
)
{
   SCIP_INTERVAL* varvals;
   int i;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   if( tree->nvars == 0 )
   {
      SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, NULL, tree->params, val) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, tree->nvars) );
   for( i = 0; i < tree->nvars; ++i )
   {
      SCIPintervalSetBounds(&varvals[i],
         -infty2infty(SCIPinfinity(scip), infinity, -SCIPvarGetLbGlobal(tree->vars[i])),
          infty2infty(SCIPinfinity(scip), infinity,  SCIPvarGetUbGlobal(tree->vars[i])));
   }

   SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, varvals, tree->params, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}

/** evaluates an expression tree w.r.t. current local bounds */
SCIP_RETCODE SCIPexprtreeEvalIntLocalBounds(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value to use for infinity */
   SCIP_INTERVAL*        val                 /**< buffer to store result */
)
{
   SCIP_INTERVAL* varvals;
   int i;

   assert(scip != NULL);
   assert(tree != NULL);
   assert(val  != NULL);

   if( tree->nvars == 0 )
   {
      SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, NULL, tree->params, val) );
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &varvals, tree->nvars) );
   for( i = 0; i < tree->nvars; ++i )
   {
      SCIPintervalSetBounds(&varvals[i],
         -infty2infty(SCIPinfinity(scip), infinity, -SCIPvarGetLbLocal(tree->vars[i])),
          infty2infty(SCIPinfinity(scip), infinity,  SCIPvarGetUbLocal(tree->vars[i])));
   }

   SCIP_CALL( SCIPexprEvalInt(tree->root, infinity, varvals, tree->params, val) );

   SCIPfreeBufferArray(scip, &varvals);

   return SCIP_OKAY;
}
