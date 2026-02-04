/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   rationalgmp.h
 * @ingroup PUBLICCOREAPI
 * @brief  wrapper for rational number arithmetic that interacts with GMP
 * @author Leon Eifler
 * @author Dominik Kamp
 *
 * This header file needs to be explicitly included by source that needs to convert between SCIP_Rational and GMP.
 * It is intentionally not included by other header files that belong to the public SCIP API.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_RATIONALGMP_H__
#define __SCIP_RATIONALGMP_H__

#include <stdlib.h>
#include "scip/def.h"
#include "scip/type_rational.h"

#if defined(SCIP_WITH_BOOST) && defined(SCIP_WITH_GMP)
#include <gmp.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicRationalMethods
 *
 * @{
 */

#if defined(SCIP_WITH_BOOST) && defined(SCIP_WITH_GMP)
/** gets the underlying gmp rational pointer */
SCIP_EXPORT
mpq_t* SCIPrationalGetGMP(
   SCIP_RATIONAL*        rational            /**< rational to access */
   );

/** sets rational to gmp rational */
SCIP_EXPORT
void SCIPrationalSetGMP(
   SCIP_RATIONAL*        rational,           /**< rational to define */
   const mpq_t           numb                /**< gmp rational to set */
   );

/** creates rational from gmp rational */
SCIP_EXPORT
SCIP_RETCODE SCIPrationalCreateBlockGMP(
   BMS_BLKMEM*           mem,                /**< block memory */
   SCIP_RATIONAL**       rational,           /**< pointer to the rational to create */
   mpq_t                 numb                /**< gmp rational to set */
   );

/** sets gmp rational array to values of rational array */
SCIP_EXPORT
void SCIPrationalSetGMPArray(
   mpq_t*                mpqaaray,           /**< gmp rational array */
   SCIP_RATIONAL**       ratarrray,          /**< rational array */
   int                   len                 /**< array length */
   );

/** sets rational array to values of gmp rational array */
SCIP_EXPORT
void SCIPrationalSetArrayGMP(
   SCIP_RATIONAL**       ratarray,           /**< rational array */
   mpq_t*                mpqarray,           /**< gmp rational array */
   int                   len                 /**< array length */
   );

/** clears gmp rational array */
SCIP_EXPORT
void SCIPrationalClearArrayGMP(
   mpq_t*                mpqarray,           /**< gmp rational array */
   int                   len                 /**< array length */
   );
#endif

/** @} */

#ifdef __cplusplus
}
#endif

#endif
