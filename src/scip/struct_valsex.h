/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_valsex.h
 * @ingroup
 * @brief  datastructures for storing exact primal CIP solutions
 * @author Kati Wolter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __SCIP_STRUCT_VALSEX_H__
#define __SCIP_STRUCT_VALSEX_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_valsex.h"

#ifdef __cplusplus
extern "C" {
#endif

/** exact primal CIP solution
 *
 *  Solutions with origin ORIGINAL contain the values for original variables. The stored objective value also
 *  corresponds to the original problem.
 */
struct SCIP_Valsex
{
   SCIP_Rational*        obj;                /**< objective value of solution */
   SCIP_RATIONALARRAY*   vals;               /**< solution values for variables */
   SCIP_BOOLARRAY*       valid;              /**< is value in vals array valid? otherwise it has to be retrieved from
                                              *   origin */
};

#ifdef __cplusplus
}
#endif

#endif
