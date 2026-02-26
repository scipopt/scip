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

/**@file   hash.c
 * @brief  unittest for the hash functions in misc.c
 * @author Rolf van der Hulst
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/pub_misc.h"

#include "include/scip_test.h"


TestSuite(hash);

Test(hash, setup_and_teardown, .description = "test that hashing numerics are correct")
{
   cr_expect_eq(SCIPrealHashCode(-1), SCIPrealHashCode(-0.9999999998));

   cr_expect_eq(SCIPrealHashCode(1.5), SCIPrealHashCode(1.4999999999));

   /* Things can go wrong around powers of two because the exponent of the floating point representation changes here */
   cr_expect_eq(SCIPrealHashCode(0.125), SCIPrealHashCode(0.12499999999999));

   cr_expect_eq(SCIPrealHashCode(0.125), SCIPrealHashCode(0.12500000005));

   cr_expect_neq(SCIPrealHashCode(0.125), SCIPrealHashCode(0.1255));

   cr_expect_neq(SCIPrealHashCode(0.125), SCIPrealHashCode(0.1245));

   /* Check if signs are correctly interpreted */
   cr_expect_eq(SCIPrealHashCode(1.00000000001), SCIPrealHashCode(0.9999999998));

   cr_expect_neq(SCIPrealHashCode(1.0), SCIPrealHashCode(-1.0));

   cr_expect_neq(SCIPrealHashCode(1.00000000001), SCIPrealHashCode(-0.9999999998));

   cr_expect_neq(SCIPrealHashCode(-1.00000000001), SCIPrealHashCode(0.9999999998));

   cr_expect_eq(SCIPrealHashCode(-1), SCIPrealHashCode(-0.9999999998));

   /* Check if sensible things happen at the numerical limits and at +- zero representation */
   cr_expect_neq(SCIPrealHashCode(DBL_MAX), SCIPrealHashCode(DBL_MIN));
   cr_expect_neq(SCIPrealHashCode(DBL_MAX), SCIPrealHashCode(0.0));
   cr_expect_neq(SCIPrealHashCode(DBL_MIN), SCIPrealHashCode(0.0));
   cr_expect_eq(SCIPrealHashCode(0.0), SCIPrealHashCode(-0.0));

}

