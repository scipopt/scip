/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_random.h
 * @brief  data structures for random number generator
 * @author Jakob Witzig
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_RANDOM_H__
#define __SCIP_STRUCT_RANDOM_H__

#include "scip/def.h"
#include "scip/type_random.h"

#ifdef __cplusplus
extern "C" {
#endif

/** random number generator data */
struct SCIP_RandGen
{
   unsigned int          seed;               /**< start seed */
   unsigned int          xor_seed;           /**< Xorshift seed */
   unsigned int          mwc_seed;           /**< Multiply-with-carry seed */
   unsigned int          cst_seed;           /**< constant seed */
};

#ifdef __cplusplus
}
#endif

#endif
