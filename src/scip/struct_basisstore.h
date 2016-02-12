/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_BASISSTORE.h
 * @brief  datastructures for storing starting basis
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BASISSTORE_H__
#define __SCIP_STRUCT_BASISSTORE_H__


#include "scip/def.h"
#include "scip/type_lp.h"
#include "scip/type_basisstore.h"

#ifdef __cplusplus
extern "C" {
#endif

/** storage for starting basis */
struct SCIP_BasisStore
{
   SCIP_BASIS**          bases;             /**< array of starting basis */
   int                   basessize;         /**< size of the basis array */
   int                   nbases;            /**< number of stored basis */
};

struct SCIP_Basis
{
   SCIP_VAR**            vars;               /**< array of variables (correspond to columns) */
   SCIP_CONS**           conss;              /**< array of constraints (correspond to rows) */
   int*                  varstat;            /**< array for storing variable/column basis status */
   int*                  consstat;           /**< array for storing constraint/row basis status */
   int                   nvars;              /**< number of LP variables/columns */
   int                   nconss;             /**< number of LP constraints/rows */
};

#ifdef __cplusplus
}
#endif

#endif
