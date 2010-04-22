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
#pragma ident "@(#) $Id: struct_expression.h,v 1.2 2010/04/22 19:15:12 bzfviger Exp $"

/**@file   struct_expression.h
 * @brief  data definitions for expressions and expression trees
 * @author Stefan Vigerske
 * @author Thorsten Gellermann
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_EXPRESSION_H__
#define __SCIP_STRUCT_EXPRESSION_H__

#include "scip/def.h"
#include "nlpi/type_expression.h"
#include "nlpi/type_exprinterpret.h"
#include "blockmemshell/memory.h"

#ifndef __SCIP_TYPE_VAR_H__
#define SCIP_VAR void
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** operator data of an expression */
union SCIP_ExprOpData {
   SCIP_VAR*             var;                /**< SCIP variable */
   int                   intval;             /**< index of a variable or parameter or a constant integer value */
   SCIP_Real             dbl;                /**< a constant double value */
   void*                 data;               /**< pointer to some data structure */
};

/** arithmetic expression node */
struct SCIP_Expr {
   SCIP_EXPROP           op;         /**< operator of the node */
   int                   nchildren;  /**< number of children */
   SCIP_EXPR**           children;   /**< children nodes */
   SCIP_EXPROPDATA       data;       /**< operator data */
};

/** expression tree */
struct SCIP_ExprTree {
   BMS_BLKMEM*           blkmem;           /**< block memory data structure */
   SCIP_EXPR*            root;             /**< root node expression of expression tree */
   int                   nvars;            /**< number of variables */
   SCIP_VAR**            vars;             /**< mapping of variable indices to SCIP variables, may be NULL if not used in context of SCIP */
   SCIP_Bool             varsasidx;        /**< are variables are stored as index (SCIP_EXPR_VARIDX)? */
   int                   nparams;          /**< number of parameters (modifiable constants) in expression */
   SCIP_Real*            params;           /**< current values for parameters, or NULL if no parameters */
   SCIP_EXPRINTDATA*     interpreterdata;  /**< data of expression interpreter (evaluator) */
};

#ifdef __cplusplus
}
#endif

#endif /* __SCIP_STRUCT_EXPRESSION_H__ */
