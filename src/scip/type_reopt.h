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

/**@file   type_reopt.h
 * @brief  type definitions for collecting reoptimization information
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_REOPT_H__
#define __SCIP_TYPE_REOPT_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Reopt SCIP_REOPT;           /**< primal data */

typedef struct SCIP_SolTree SCIP_SOLTREE;       /**< tree to check solutions */

typedef struct SCIP_SolNode SCIP_SOLNODE;       /**< nodes of SCIP_SOLTREE */

/* type of nodes during reoptimization */
enum SCIP_ReoptType
{
   SCIP_REOPTTYPE_NONE        = 0,
   SCIP_REOPTTYPE_TRANSIT     = 1,
   SCIP_REOPTTYPE_STRBRANCHED = 2,
   SCIP_REOPTTYPE_LOGICORNODE = 3,
   SCIP_REOPTTYPE_LEAF        = 4,
   SCIP_REOPTTYPE_INFEASIBLE  = 5,
   SCIP_REOPTTYPE_PRUNED      = 6,
   SCIP_REOPTTYPE_FEASIBLE    = 7
};
typedef enum SCIP_ReoptType SCIP_REOPTTYPE;     /**< type nodes during reoptimization */

#ifdef __cplusplus
}
#endif

#endif
