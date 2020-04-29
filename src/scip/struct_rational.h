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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
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
   Rational val;                             /**< value of the rational */
   unsigned int isinf:1;                     /**< is the value infinite? sign is determined by val */
   unsigned int isfprepresentable:2;         /**< is the value exactly representable as floating point number?
                                              *   (0 - unknown, 1 - yes, 2 - no) */
};

/** rational array struct, essentially a std vector with all indices offset by firstidx*/
struct SCIP_RationalArray
{
   std::vector<SCIP_Rational> vals;          /**< values of the array */
   int firstidx;                             /**< first used index */
};

#ifdef __cplusplus
}
#endif

#endif
