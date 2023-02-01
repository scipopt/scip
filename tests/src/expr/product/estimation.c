/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2023 Zuse Institute Berlin (ZIB)                      */
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

/**@file   estimation.c
 * @brief  tests estimation of products
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/expr_product.c"
#include "../estimation.h"

Test(estimation, bilinear, .init = setup, .fini = teardown,
   .description = "test estimation for a bilinear expression"
   )
{
   SCIP_EXPR* expr;
   SCIP_Real coefs[2];
   SCIP_Real constant;
   SCIP_INTERVAL bnds[2];
   SCIP_Real ref[2];
   SCIP_Bool islocal;
   SCIP_Bool branchcand = TRUE;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreateExprProduct(scip, &expr, 0, NULL, 1.5, NULL, NULL) );
   SCIP_CALL( SCIPappendExprChild(scip, expr, xexpr) );
   SCIP_CALL( SCIPappendExprChild(scip, expr, yexpr) );

   /*
    * compute an overestimator for 1.5*x*y with x* = 0, y* = -4
    * together with the bounds this should result in an estimator of the form
    *    -4.5x - 1.5y - 4.5
    */
   ref[0] = 0.0;
   ref[1] = -4.0;

   bnds[0].inf = SCIPvarGetLbLocal(x);
   bnds[0].sup = SCIPvarGetUbLocal(x);
   bnds[1].inf = SCIPvarGetLbLocal(y);
   bnds[1].sup = SCIPvarGetUbLocal(y);

   SCIP_CALL( estimateProduct(scip, expr, bnds, bnds, ref, TRUE, SCIPinfinity(scip), coefs, &constant, &islocal, &success, &branchcand) );

   cr_assert(success);
   cr_assert_float_eq(constant, -4.5, SCIPepsilon(scip));
   cr_assert_float_eq(coefs[0], -4.5, SCIPepsilon(scip));
   cr_assert_float_eq(coefs[1], -1.5, SCIPepsilon(scip));
   cr_assert(islocal);
   cr_assert(branchcand);


   /*
    * compute an underestimator for 1.5*x*y with x* = 0, y* = -4
    * together with the bounds this should result in en estimator of the form
    *    -9x - 1.5y - 9
    */
   SCIP_CALL( estimateProduct(scip, expr, bnds, bnds, ref, FALSE, -SCIPinfinity(scip), coefs, &constant, &islocal, &success, &branchcand) );

   cr_assert(success);
   cr_assert_float_eq(constant, -9.0, SCIPepsilon(scip));
   cr_assert_float_eq(coefs[0], -9.0, SCIPepsilon(scip));
   cr_assert_float_eq(coefs[1], -1.5, SCIPepsilon(scip));
   cr_assert(islocal);
   cr_assert(branchcand);

   /* release expression */
   SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
}


/* This next test separates a 4-linear function and then a 3-linear function:
 * The exact facets were computed using CDD and the following (ugly) Julia code:
 *
 * using Iterators, Polyhedra, CDDLib
 * #Compute (vertex, f(vertex)) where the domain is [-0.2, 0.7], [-10, 8], [1, 1.3], [0.09, 2.1] and f is -0.7 x*y*w*z
 * # choose only one of the following
 * domz = [1//1,1//1]; zstar = 1; tstar =-3.4;
 * domz = [9//100,21//10]; zstar = 0.18; tstar = 12.12;
 *
 * # get vertices in matrix notation for CDD
 * vertices = map(t->vcat(t[2:end]..., prod(t)), collect(product(-7//10,[-2//10, 7//10], [-10, 8], [1, 13//10], domz)));
 * vert = vcat([vertices[i]' for i in 1:length(vertices)]...);
 *
 * # compute polyhedron and get its H-representation
 * vrepre = SimpleVRepresentation(vert);
 * poly = polyhedron(vrepre, CDDLibrary(:exact));
 * hrepre = SimpleHRepresentation(gethrep(poly));
 *
 * #find the most violated facet
 * m = size(hrepre.A)[1] # number of facets
 * values = hrepre.A * [0.2, -4.0, 1.1, zstar, tstar]  -  hrepre.b
 * separates = values  .> 0 # this are the inequalities that separate
 * maxel = 0.0; argmax = -1;
 * for i in 1:m
 *     if separates[i]
 *         value = dot(hrepre.A[i,:], [0.2, -4.0, 1.1, zstar, tstar])  -  hrepre.b[i]
 *         value /= abs(dot(hrepre.A[i,:], [0,0,0,0, 1])) # we have to normalize
 *         @show dot(hrepre.A[i,:], [0,0,0,0, 1])
 *         if maxel <= value
 *             maxel = value
 *             argmax = i
 *         end
 *     end
 * end
 * println(maxel, " index: ", argmax)
 * violation = (hrepre.A[argmax,:]' * [0.2, -4.0, 1.1, zstar, tstar] - hrepre.b[argmax])/ hrepre.A[argmax,end]
 * */

Test(estimation, quadrilinear,
   .description = "test separation for a quadrilinear expression"
   )
{
   SCIP_EXPR* expr;
   SCIP_Bool islocal;
   SCIP_Bool branchcand = TRUE;
   SCIP_Bool success;
   SCIP_Real facetcoefs[4];
   SCIP_Real facetconstant;
   int i;
   int round;

   SCIP_VAR* vars[4];
   SCIP_EXPR* varexprs[4];
   SCIP_INTERVAL bnds[4];
   const char* names[4] = { "x", "y", "w", "z" };
   SCIP_Real lb[] = {-0.2, -10.0, 1.0, 0.09};
   SCIP_Real ub[] = { 0.7,   8.0, 1.3,  2.1};
   SCIP_Real solval[] = { 0.2, -4.0, 1.1, 0.18};
   SCIP_Real exact_facet[2][5] = {{63.0/100, 63.0/5000, 441.0/1000, 637.0/100, -8883.0/10000},{7.0, -49.0/100, -98.0/25, 0.0, -49.0/50}};

   /* first round: overestimation, expect exact_facet[1]
    * second round: underestimation, expect exact_facet[2]
    */
   for( round = 0; round < 2 ; ++round )
   {
      SCIP_CALL( SCIPcreate(&scip) );
      SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

      /* create problem */
      SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

      for( i = 0; i < 4; ++i )
      {
         SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], names[i], lb[i], ub[i], 1.0, SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, vars[i]) );
         SCIP_CALL( SCIPcreateExprVar(scip, &varexprs[i], vars[i], NULL, NULL) );
         bnds[i].inf = lb[i];
         bnds[i].sup = ub[i];
      }

      SCIP_CALL( SCIPcreateExprProduct(scip, &expr, 4, varexprs, -0.7, NULL, NULL) );

      /* round 0:
       * compute an overestimator for -0.7*x*y*w*z with x* = 0.2, y* = -4, w* = 1.1, z* = 0.18
       * together with the bounds x,y,w,z \in [-0.2, 0.7], [-10, 8], [1, 1.3], [0.09, 2.1]
       * -> 63/5000 * (50x + y + 35w + 4550/9 z - 141/2)
       *
       * round 1:
       * compute an underestimator for the same function, but now z is fixed to 1
       * so, we underestimate
       *   -0.7*x*y*w with x* = 0.2, y* = -4, w* = 1.1
       * together with the bounds x,y,w \in [-0.2, 0.7], [-10, 8], [1, 1.3]
       * -> -49/100 * (-100/7 x + y + 8w + 2)
       */
      if( round == 1 )
      {
         bnds[3].inf = 1.0;
         bnds[3].sup = 1.0;
         solval[3] = 1.0;
      }

      SCIP_CALL( estimateProduct(scip, expr, bnds, bnds, solval, round == 0, (round == 0 ? SCIPinfinity(scip) : -SCIPinfinity(scip)), facetcoefs, &facetconstant, &islocal, &success, &branchcand) );

      cr_assert(success);
      cr_assert(islocal);
      cr_assert(branchcand);
      for( i = 0; i < 4; ++i ) /* index 4 is the constant */
      {
         cr_expect_float_eq(facetcoefs[i], exact_facet[round][i], SCIPfeastol(scip), "coef %d: received %g instead of %g\n", i, facetcoefs[i], exact_facet[round][i]);
      }
      cr_expect_float_eq(facetconstant, exact_facet[round][i], SCIPfeastol(scip), "constant: received %g instead of %g\n", facetconstant, exact_facet[round][i]);

      /* release and free everything */
      SCIP_CALL( SCIPreleaseExpr(scip, &expr) );
      for( i = 0; i < 4; ++i )
      {
         SCIP_CALL( SCIPreleaseExpr(scip, &varexprs[i]) );
         SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
      }
      SCIP_CALL( SCIPfree(&scip) );

      cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
   }
}
