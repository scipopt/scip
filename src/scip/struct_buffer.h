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

/**@file   struct_buffer.h
 * @brief  datastructures for memory buffers for temporary objects
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_BUFFER_H__
#define __SCIP_STRUCT_BUFFER_H__


#include <assert.h>

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

/** memory buffer storage for temporary objects */
struct SCIP_Buffer
{
   void**                data;               /**< allocated memory chunks for arbitrary data */
   int*                  size;               /**< sizes of buffers in bytes */
   SCIP_Bool*            used;               /**< TRUE iff corresponding buffer is in use */
   int                   ndata;              /**< number of memory chunks */
   int                   firstfree;          /**< first unused memory chunk */
};

#ifdef __cplusplus
}
#endif

#endif
