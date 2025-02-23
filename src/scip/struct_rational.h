/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   struct_rational.h
 * @brief  datastructures for storing rational numbers
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_RATIONAL_H__
#define __SCIP_STRUCT_RATIONAL_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_rational.h"
#include "scip/multiprecision.hpp"

#ifdef __cplusplus
extern "C" {
#endif

/** rational wrapper struct */
struct SCIP_Rational
{
   scip_rational::Rational val;              /**< value of the rational */
   unsigned int isinf:1;                     /**< is the value infinite? sign is determined by val */
   unsigned int isfprepresentable:2;         /**< is the value exactly representable as floating point number?
                                              *   (0 - unknown, 1 - yes, 2 - no) */
};

/** rational array struct, essentially a std vector with all indices offset by firstidx*/
struct SCIP_RationalArray
{
   /* MP@LE Is this a good idea to have an STL vector here instead of an array? */
   std::vector<SCIP_Rational> vals;          /**< values of the array */
   int firstidx;                             /**< first used index */
};

#ifdef __cplusplus
}
#endif

#endif
