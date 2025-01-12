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

/**@file   type_rational.h
 * @brief  type definitions for rational numbers
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_RATIONAL_H__
#define __SCIP_TYPE_RATIONAL_H__

#ifdef __cplusplus
extern "C" {
#endif

/** type used for rational numbers */
typedef struct SCIP_Rational SCIP_Rational;

/** dynamic array for storing SCIP_Real values */
typedef struct SCIP_RationalArray SCIP_RATIONALARRAY;

/** information if a rational is exactly representable as a floating point number */
enum SCIP_IsFpRepresentable
{
   SCIP_ISFPREPRESENTABLE_UNKNOWN = 0,       /**< representability is unknown */
   SCIP_ISFPREPRESENTABLE_TRUE    = 1,       /**< is representable */
   SCIP_ISFPREPRESENTABLE_FALSE   = 2        /**< is not representable */
};
typedef enum SCIP_IsFpRepresentable SCIP_ISFPREPRESENTABLE;

/** defines the possible rounding direction for a rational number, when converting to a double */
enum SCIP_RoundModeRational
{
   SCIP_R_ROUND_DOWNWARDS = 0,               /**< always round to nearest smaller double */
   SCIP_R_ROUND_UPWARDS   = 1,               /**< always round to nearest larger double */
   SCIP_R_ROUND_NEAREST   = 2                /**< always round to nearest double */
};
typedef enum SCIP_RoundModeRational SCIP_ROUNDMODE_RAT;

#ifdef __cplusplus
}
#endif

#endif
