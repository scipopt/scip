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
#pragma ident "@(#) $Id: exprinterpret_none.c,v 1.2 2010/05/10 19:03:33 bzfviger Exp $"

/**@file   exprinterpret_none.c
 * @brief  function definitions for nonexisting expression interpreter to resolve linking references 
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "nlpi/exprinterpret.h"

struct SCIP_ExprInt
{
};

const char* SCIPexprintGetName(void)
{
	return "NONE";
}

/** creates an expression interpreter object */
SCIP_RETCODE SCIPexprintCreate(
   BMS_BLKMEM*           blkmem,             /**< block memory data structure */
   SCIP_EXPRINT**        exprint             /**< buffer to store pointer to expression interpreter */
)
{
	SCIPdebugMessage("SCIPexprintCreate()\n");
	SCIPdebugMessage("Note that there is no expression interpreter linked to the binary.\n");
	
   if( BMSallocMemory(exprint) == NULL )
      return SCIP_NOMEMORY;
   
	return SCIP_OKAY;
}

/** frees an expression interpreter object */
SCIP_RETCODE SCIPexprintFree(
   SCIP_EXPRINT**        exprint             /**< expression interpreter that should be freed */
)
{
   BMSfreeMemory(exprint);
   
	return SCIP_OKAY;
}

/** compiles an expression tree and stores compiled data in expression tree */
SCIP_RETCODE SCIPexprintCompile(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree                /** expression tree */
)
{
	return SCIP_PLUGINNOTFOUND;
}

/** frees interpreter data */
SCIP_RETCODE SCIPexprintFreeData(
   SCIP_EXPRINTDATA**    interpreterdata     /** interpreter data that should freed */
)
{
	return SCIP_PLUGINNOTFOUND;
}

/** notify expression interpreter that a new parameterization is used
 * this probably causes retaping by AD algorithms
 */
SCIP_RETCODE SCIPexprintNewParametrization(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree                /** expression tree */
)
{
	return SCIP_PLUGINNOTFOUND;
}

/** evaluates an expression tree */
SCIP_RETCODE SCIPexprintEval(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables */
   SCIP_Real*            val                 /** buffer to store value */
)
{
	return SCIP_PLUGINNOTFOUND;
}

/** evaluates an expression tree on intervals */
SCIP_RETCODE SCIPexprintEvalInt(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real             infinity,           /** value for infinity */
   SCIP_INTERVAL*        varvals,            /** interval values of variables */
   SCIP_INTERVAL*        val                 /** buffer to store interval value of expression */
)
{
   return SCIP_PLUGINNOTFOUND;
}

/** gets number of nonzeros in gradient of expression tree */
SCIP_RETCODE SCIPexprintGetNGradPattern(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   int*                  nnz                 /** buffer to store number of nonzeros */
)
{
	return SCIP_PLUGINNOTFOUND;
}

/** gets sparsity pattern of expression trees gradient */
SCIP_RETCODE SCIPexprintGetGradPattern(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   int*                  gradidx             /** buffer to store gradient indices */
)
{
	return SCIP_PLUGINNOTFOUND;
}

/** computes value and gradient of an expression tree */
SCIP_RETCODE SCIPexprintGrad(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /** have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /** buffer to store value */
   SCIP_Real*            gradvals            /** buffer to store gradient values */
)
{
	return SCIP_PLUGINNOTFOUND;
}

/** computes value and dense gradient of an expression tree */
SCIP_RETCODE SCIPexprintGradDense(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /** have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /** buffer to store value */
   SCIP_Real*            gradient            /** buffer to store gradient */
)
{
   return SCIP_PLUGINNOTFOUND;
}

/** computes interval value and dense interval gradient of an expression tree */
SCIP_RETCODE SCIPexprintGradDenseInt(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real             infinity,           /** value for infinity */
   SCIP_INTERVAL*        varvals,            /** interval values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /** have variable values changed since last call to an interval evaluation routine? */
   SCIP_INTERVAL*        val,                /** buffer to store expression interval value */
   SCIP_INTERVAL*        gradient            /** buffer to store expression interval gradient */
)
{
   return SCIP_PLUGINNOTFOUND;
}

/** gives sparsity pattern of hessian
 * NOTE: this function might be replaced later by something nicer 
 * Since the AD code might need to do a forward sweep, you should pass variable values in here.
 */
SCIP_RETCODE SCIPexprintHessianSparsityDense(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables */
   SCIP_Bool*            sparsity            /** buffer to store sparsity pattern of Hessian, sparsity[i+n*j] indicates whether entry (i,j) is nonzero in the hessian */ 
)
{
   return SCIP_PLUGINNOTFOUND;
}

/** computes value and dense hessian of an expression tree
 * the full hessian is computed (lower left and upper right triangle)
 */
SCIP_RETCODE SCIPexprintHessianDense(
   SCIP_EXPRINT*         exprint,            /** interpreter data structure */
   SCIP_EXPRTREE*        tree,               /** expression tree */
   SCIP_Real*            varvals,            /** values of variables, can be NULL if new_varvals is FALSE */
   SCIP_Bool             new_varvals,        /** have variable values changed since last call to an evaluation routine? */
   SCIP_Real*            val,                /** buffer to store function value */
   SCIP_Real*            hessian             /** buffer to store hessian values, need to have size at least n*n */
)
{
	return SCIP_PLUGINNOTFOUND;
}
