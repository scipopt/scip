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
 * @ingroup
 * @brief  datastructures for storing rational numbers
 * @author Leon Eifler
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_STRUCT_RATIONAL_H__
#define __SCIP_STRUCT_RATIONAL_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_rational.h"

#ifdef __cplusplus
extern "C" {
#endif

/** rational wrapper struct */
struct SCIP_Rational
{
   Rational val;
   unsigned int isinf:1;
   unsigned int isfprepresentable:2;
};

/** rational array struct */
struct SCIP_RationalArray
{
   std::vector<SCIP_Rational> vals;
   int firstidx;
};

#ifdef __cplusplus
}
#endif

#endif
