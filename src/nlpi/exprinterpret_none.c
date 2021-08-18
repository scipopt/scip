/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    exprinterpret_none.c
 * @brief   function definitions for nonexisting expression interpreter to resolve linking references
 * @ingroup EXPRINTS
 * @author  Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_message.h"
#include "nlpi/exprinterpret.h"

/** gets name and version of expression interpreter */
const char* SCIPexprintGetName(
   void
   )
{
   return "NONE";
}

/** gets descriptive text of expression interpreter */
const char* SCIPexprintGetDesc(
   void
   )
{
   return "dummy expression interpreter which solely purpose it is to resolve linking symbols";
}

/** gets capabilities of expression interpreter (using bitflags) */
SCIP_EXPRINTCAPABILITY SCIPexprintGetCapability(
   void
   )
{
   return SCIP_EXPRINTCAPABILITY_NONE;
}

/** creates an expression interpreter object */ /*lint -e{715}*/
SCIP_RETCODE SCIPexprintCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
   )
{  /*lint --e{715}*/
   SCIPdebugMessage("SCIPexprintCreate()\n");
   SCIPdebugMessage("Note that there is no expression interpreter linked to the binary.\n");

   *exprint = (SCIP_EXPRINT*)1u;  /* some code checks that a non-NULL pointer is returned here, even though it may not point anywhere */

   return SCIP_OKAY;
}

/** frees an expression interpreter object */
SCIP_RETCODE SCIPexprintFree(
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
   )
{
   assert(exprint != NULL);
   *exprint = NULL;

   return SCIP_OKAY;
}

/** compiles an expression tree and stores compiled data in expression tree */ /*lint --e{715}*/
SCIP_RETCODE SCIPexprintCompile(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}


/** gives the capability to evaluate an expression by the expression interpreter
 *
 * In cases of user-given expressions, higher order derivatives may not be available for the user-expression,
 * even if the expression interpreter could handle these. This method allows to recognize that, e.g., the
 * Hessian for an expression is not available because it contains a user expression that does not provide
 * Hessians.
 */ /*lint -e{715}*/
SCIP_EXPRINTCAPABILITY SCIPexprintGetExprtreeCapability(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{  /*lint --e{715}*/
   return SCIP_EXPRINTCAPABILITY_NONE;
}

/** frees interpreter data */
SCIP_RETCODE SCIPexprintFreeData(
   SCIP_EXPRINTDATA**    interpreterdata     /**< interpreter data that should freed */
   )
{
   assert(interpreterdata  != NULL);
   assert(*interpreterdata == NULL);

   return SCIP_OKAY;
}

/** notify expression interpreter that a new parameterization is used
 * this probably causes retaping by AD algorithms
 */ /*lint -e{715}*/
SCIP_RETCODE SCIPexprintNewParametrization(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree                /**< expression tree */
   )
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/** evaluates an expression tree */ /*lint -e{715}*/
SCIP_RETCODE SCIPexprintEval(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Real*            val                 /**< buffer to store value */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}

/** evaluates an expression tree on intervals */ /*lint -e{715}*/
SCIP_RETCODE SCIPexprintEvalInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables */
   SCIP_INTERVAL*        val                 /**< buffer to store interval value of expression */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}

/** computes value and gradient of an expression tree */ /*lint -e{715}*/
SCIP_RETCODE SCIPexprintGrad(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to a point evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store expression value */
   SCIP_Real*            gradient            /**< buffer to store expression gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}

/** computes interval value and interval gradient of an expression tree */ /*lint -e{715}*/
SCIP_RETCODE SCIPexprintGradInt(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real             infinity,           /**< value for infinity */
   SCIP_INTERVAL*        varvals,            /**< interval values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable interval values changed since last call to an interval evaluation routine? */
   SCIP_INTERVAL*        val,                /**< buffer to store expression interval value */
   SCIP_INTERVAL*        gradient            /**< buffer to store expression interval gradient, need to have length at least SCIPexprtreeGetNVars(tree) */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}

/** gives sparsity pattern of hessian
 * NOTE: this function might be replaced later by something nicer 
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 */  /*lint -e{715}*/
SCIP_RETCODE SCIPexprintHessianSparsityDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables */
   SCIP_Bool*            sparsity            /**< buffer to store sparsity pattern of Hessian, sparsity[i+n*j] indicates whether entry (i,j) is nonzero in the hessian */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}

/** computes value and dense hessian of an expression tree
 * the full hessian is computed (lower left and upper right triangle)
 */ /*lint -e{715}*/
SCIP_RETCODE SCIPexprintHessianDense(
   SCIP_EXPRINT*         exprint,            /**< interpreter data structure */
   SCIP_EXPRTREE*        tree,               /**< expression tree */
   SCIP_Real*            varvals,            /**< values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /**< have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /**< buffer to store function value */
   SCIP_Real*            hessian             /**< buffer to store hessian values, need to have size at least n*n */
   )
{  /*lint --e{715}*/
   SCIPerrorMessage("No expression interpreter linked to SCIP, try recompiling with EXPRINT=cppad.\n");
   return SCIP_PLUGINNOTFOUND;
}
