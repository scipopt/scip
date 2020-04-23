/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
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

/**< type used for rational numbers */
typedef struct SCIP_Rational SCIP_Rational;

/** dynamic array for storing SCIP_Real values */
typedef struct SCIP_RationalArray SCIP_RATIONALARRAY;

/** information if a rational is exactly representable as a floating point number */
enum SCIP_IsFpRepresentable
{
   SCIP_ISFPREPRESENTABLE_UNKNOWN = 0,
   SCIP_ISFPREPRESENTABLE_TRUE = 1,
   SCIP_ISFPREPRESENTABLE_FALSE = 2
};
typedef enum SCIP_IsFpRepresentable SCIP_ISFPREPRESENTABLE;

/** defines the possible rounding direction for a rational number, when converting to a double */
enum SCIP_RoundModeRational
{
   SCIP_ROUND_DOWNWARDS = 0,     /**< always round to nearest smaller double*/
   SCIP_ROUND_UPWARDS = 1,       /**< always round to nearest larger double*/
   SCIP_ROUND_NEAREST = 2        /**< always round to nearest double*/
};
typedef enum SCIP_RoundModeRational SCIP_ROUND;

#ifdef __cplusplus
}
#endif

#endif
