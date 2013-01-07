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

/**@file   mem.h
 * @brief  methods for block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_MEM_H__
#define __SCIP_MEM_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_mem.h"
#include "scip/struct_mem.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates block memory structures */
extern
SCIP_RETCODE SCIPmemCreate(
   SCIP_MEM**            mem                 /**< pointer to block memory structure */
   );

/** frees block memory structures */
extern
SCIP_RETCODE SCIPmemFree(
   SCIP_MEM**            mem                 /**< pointer to block memory structure */
   );

/** returns the total number of bytes used in block memory */
extern
SCIP_Longint SCIPmemGetUsed(
   SCIP_MEM*             mem                 /**< pointer to block memory structure */
   );

#ifdef __cplusplus
}
#endif

#endif
