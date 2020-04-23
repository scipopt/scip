/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_varex.c
 * @brief  public methods for exact SCIP variables
 * @author Leon Eifler
 *
 * @todo check all SCIP_STAGE_* switches, and include the new stages TRANSFORMED and INITSOLVE
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "blockmemshell/memory.h"
#include "lpi/lpi.h"
#include "scip/branch.h"
#include "scip/clock.h"
#include "scip/conflict.h"
#include "scip/debug.h"
#include "scip/history.h"
#include "scip/implics.h"
#include "scip/lp.h"
#include "scip/prob.h"
#include "scip/pub_cons.h"
#include "scip/pub_implics.h"
#include "scip/pub_lp.h"
#include "scip/pub_message.h"
#include "scip/pub_misc.h"
#include "scip/pub_tree.h"
#include "scip/pub_var.h"
#include "scip/relax.h"
#include "scip/rational.h"
#include "scip/scip_general.h"
#include "scip/scip_lp.h"
#include "scip/scip_mem.h"
#include "scip/scip_message.h"
#include "scip/scip_numerics.h"
#include "scip/scip_prob.h"
#include "scip/scip_probing.h"
#include "scip/scip_sol.h"
#include "scip/scip_solvingstats.h"
#include "scip/scip_tree.h"
#include "scip/scip_var.h"
#include "scip/set.h"
#include "scip/sol.h"
#include "scip/solve.h"
#include "scip/stat.h"
#include "scip/struct_lp.h"
#include "scip/struct_mem.h"
#include "scip/struct_primal.h"
#include "scip/struct_prob.h"
#include "scip/struct_scip.h"
#include "scip/struct_set.h"
#include "scip/struct_stat.h"
#include "scip/struct_tree.h"
#include "scip/struct_var.h"
#include "scip/tree.h"
#include "scip/var.h"
#include "scip/varex.h"
#include <ctype.h>

/** changes variable's objective value
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_PRESOLVED
 */
SCIP_RETCODE SCIPchgVarObjExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< variable to change the objective value for */
   SCIP_Rational*        newobj              /**< new objective value */
   )
{
   SCIP_CALL( SCIPcheckStage(scip, "SCIPchgVarObj", FALSE, TRUE, TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   assert( var->scip == scip );

   /* forbid infinite objective values */
   if( RatIsAbsInfinity(newobj) )
   {
      SCIPerrorMessage("invalid objective value: objective value is infinite\n");
      return SCIP_INVALIDDATA;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      assert(!SCIPvarIsTransformed(var));
      SCIP_CALL( SCIPvarChgObjExact(var, scip->mem->probmem, scip->set, scip->origprob, scip->primal, scip->lpex, scip->eventqueue, newobj) );
      return SCIP_OKAY;

   case SCIP_STAGE_TRANSFORMED:
   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_PRESOLVED:
      SCIP_CALL( SCIPvarChgObjExact(var, scip->mem->probmem, scip->set,  scip->transprob, scip->primal, scip->lpex, scip->eventqueue, newobj) );
      return SCIP_OKAY;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/
}

/** Add exact data to variable
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_RETCODE SCIPaddVarExactData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR*             var,                /**< pointer to variable */
   SCIP_Rational*        lb,                 /**< lower bounf of variable */
   SCIP_Rational*        ub,                 /**< upper bound of variable */
   SCIP_Rational*        obj                 /**< objective function value */
   )
{
   assert(var != NULL);
   assert(lb != NULL);
   assert(ub != NULL);

   assert(RatIsLE(lb, ub));

   SCIP_CALL( SCIPcheckStage(scip, "SCIPcreateVar", FALSE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE) );

   /* forbid infinite objective function values */
   if( obj != NULL && RatIsAbsInfinity(obj) )
   {
      SCIPerrorMessage("invalid objective function value: value is infinite\n");
      return SCIP_INVALIDDATA;
   }

   switch( scip->set->stage )
   {
   case SCIP_STAGE_PROBLEM:
      SCIP_CALL( SCIPvarAddExactData(var, scip->mem->probmem, lb, ub, obj) );
      break;

   case SCIP_STAGE_TRANSFORMING:
   case SCIP_STAGE_INITPRESOLVE:
   case SCIP_STAGE_PRESOLVING:
   case SCIP_STAGE_EXITPRESOLVE:
   case SCIP_STAGE_PRESOLVED:
   case SCIP_STAGE_SOLVING:
      SCIP_CALL( SCIPvarAddExactData(var, scip->mem->probmem, lb, ub, obj) );
      break;

   default:
      SCIPerrorMessage("invalid SCIP stage <%d>\n", scip->set->stage);
      return SCIP_INVALIDCALL;
   }  /*lint !e788*/

   return SCIP_OKAY;
}

/** print the given variables and rational coefficients as linear sum in the following form
 *  c1 \<x1\> + c2 \<x2\>   ... + cn \<xn\>
 *
 *  This string can be parsed by the method SCIPparseVarsLinearsum().
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_PROBLEM
 *       - \ref SCIP_STAGE_TRANSFORMING
 *       - \ref SCIP_STAGE_TRANSFORMED
 *       - \ref SCIP_STAGE_INITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVING
 *       - \ref SCIP_STAGE_EXITPRESOLVE
 *       - \ref SCIP_STAGE_PRESOLVED
 *       - \ref SCIP_STAGE_INITSOLVE
 *       - \ref SCIP_STAGE_SOLVING
 *       - \ref SCIP_STAGE_SOLVED
 *       - \ref SCIP_STAGE_EXITSOLVE
 *       - \ref SCIP_STAGE_FREETRANS
 *
 *  @note The printing process is done via the message handler system.
 */
SCIP_RETCODE SCIPwriteVarsExactLinearsum(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< output file, or NULL for stdout */
   SCIP_VAR**            vars,               /**< variable array to output */
   SCIP_Rational**       vals,               /**< array of coefficients or NULL if all coefficients are 1.0 */
   int                   nvars,              /**< number of variables */
   SCIP_Bool             type                /**< should the variable type be also posted */
   )
{
   int v;

   SCIP_CALL( SCIPcheckStage(scip, "SCIPwriteVarsLinearsum", FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE) );

   for( v = 0; v < nvars; ++v )
   {
      if( vals != NULL )
      {
         if( RatIsEqualReal(vals[v], 1.0) )
         {
            if( v > 0 )
               SCIPinfoMessage(scip, file, " +");
         }
         else if( RatIsEqualReal(vals[v], -1.0) )
            SCIPinfoMessage(scip, file, " -");
         else
         {
            char buf[SCIP_MAXSTRLEN];
            RatToString(vals[v], buf, SCIP_MAXSTRLEN);
            SCIPinfoMessage(scip, file, " %s", buf);
         }
      }
      else if( nvars > 0 )
         SCIPinfoMessage(scip, file, " +");

      /* print variable name */
      SCIP_CALL( SCIPwriteVarName(scip, file, vars[v], type) );
   }

   return SCIP_OKAY;
}
