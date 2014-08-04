/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_reopt.h
 * @brief  datastructures for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_REOPT_H__
#define __SCIP_STRUCT_REOPT_H__


#include "scip/def.h"
#include "scip/type_reopt.h"

#ifdef __cplusplus
extern "C" {
#endif

/** reopt data and solution storage */
struct SCIP_Reopt
{
   SCIP_SOLNODE***       sols;               /**< solutions of the reoptimization runs */

   int                   run;                /**< current position in the sols array*/
   int                   runsize;            /**< allocated memory for runs */
   int*                  solssize;           /**< size of sols[x] arrays */
   int*                  nsols;              /**< number of solutions stored in sols[x] array */

   SCIP_Real**           objs;               /**< list of objective coefficients */
   SCIP_SOL*             lastbestsol;        /**< best solution of the last round */
   SCIP_SOLTREE*         soltree;            /**< tree to handle all saved solutions */

   SCIP_Real             simtolastobj;       /**< similarity to the last objective function */
   SCIP_Real             simtofirstobj;      /**< similarity to the first objective function */
};

/** nodes of SCIP_SolTree */
struct SCIP_SolNode
{
   SCIP_SOL*             sol;
   SCIP_SOLNODE*         father;
   SCIP_SOLNODE*         rchild;
   SCIP_SOLNODE*         lchild;
   SCIP_Bool             updated;
   SCIP_Bool             used;
   SCIP_Bool             infeasible;
   int                   val;
};

/** tree for solution */
struct SCIP_SolTree
{
   SCIP_SOLNODE*         root;
   int                   nsols;
};

#ifdef __cplusplus
}
#endif

#endif
