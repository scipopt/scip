/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2018 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   separation_prod.c
 * @brief  tests separation of products
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/cons_expr_product.c"
#include "separation.h"

Test(separation, bilinear, .init = setup, .fini = teardown,
   .description = "test separation for a bilinear expression"
   )
{
   SCIP_CONSEXPR_EXPR* expr;
   SCIP_Real coefs[2];
   SCIP_Real constant;
   SCIP_Bool islocal;
   SCIP_Bool success;

   SCIP_CALL( SCIPcreateConsExprExprProduct(scip, conshdlr, &expr, 0, NULL, 1.5) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, xexpr) );
   SCIP_CALL( SCIPappendConsExprExprProductExpr(scip, expr, yexpr) );

   /*
    * compute an overestimator for 1.5*x*y with x* = 0, y* = -4
    * together with the bounds this should result in an estimator of the form
    *    -4.5x - 1.5y - 4.5
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );

   SCIP_CALL( estimateProduct(scip, conshdlr, expr, sol, TRUE, SCIPinfinity(scip), coefs, &constant, &islocal, &success) );

   cr_assert(success);
   cr_assert_float_eq(constant, -4.5, SCIPepsilon(scip));
   cr_assert_float_eq(coefs[0], -4.5, SCIPepsilon(scip));
   cr_assert_float_eq(coefs[1], -1.5, SCIPepsilon(scip));
   cr_assert(islocal);


   /*
    * compute an underestimator for 1.5*x*y with x* = 0, y* = -4
    * together with the bounds this should result in en estimator of the form
    *    -9x - 1.5y - 9
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );

   SCIP_CALL( estimateProduct(scip, conshdlr, expr, sol, FALSE, -SCIPinfinity(scip), coefs, &constant, &islocal, &success) );

   cr_assert(success);
   cr_assert_float_eq(constant, -9.0, SCIPepsilon(scip));
   cr_assert_float_eq(coefs[0], -9.0, SCIPepsilon(scip));
   cr_assert_float_eq(coefs[1], -1.5, SCIPepsilon(scip));
   cr_assert(islocal);

   /* release expression */
   SCIP_CALL( SCIPreleaseConsExprExpr(scip, &expr) );
}

/*
 * the following test builds the matrices used to compute facets of convex/concave envelope from size 0 to 8 (otherwise
 * we get a overlength string warning) and checks it is correctly built by printing it to std output and comparing with
 * the matrices in matrices.sol
 */

#include "matrices.sol"

/* specify parameters of parameterized test */
ParameterizedTestParameters(separation, multilinearLP)
{
   static const int sizes[] = {1,2,3,4,5,6,7};

   /* type of the parameter; the parameter; number of parameters */
   return cr_make_param_array(const int, sizes, sizeof(sizes)/sizeof(int));
}

static
void printMatrix(int size)
{
   int i, j, nrows, ncols;
   SCIP_LPI* lp;

   SCIP_CALL( SCIPcreate(&scip) );

   cr_redirect_stdout();
   SCIP_CALL( buildMultilinearSeparationLP(scip, size, &lp) );

   SCIP_CALL( SCIPlpiGetNRows(lp, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lp, &ncols) );

   for( i = 0; i < nrows; ++i )
   {
      for( j = 0; j < ncols; ++j )
      {
         SCIP_Real val;
         SCIP_CALL( SCIPlpiGetCoef(lp, i, j, &val) );
         printf("%g", val);
      }
      printf("\n");
   }
   fflush(stdout);
   cr_assert_stdout_eq_str(matrices[size]);

   SCIP_CALL( SCIPlpiFree(&lp) );
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

/* generates matrix of size *size and prints it; checks it is the expected matrix */
ParameterizedTest(const int* size, separation, multilinearLP)
{
   printMatrix(*size);
}

/*
 * the following test checks that separation of multilinear terms is performed correctly, even when using matrices of
 * different sizes. It also test the re-use of the separation LP
 */

static void testBilinearLP(int lpsize)
{
   SCIP_VAR* vars[2];
   SCIP_Real facetcoefs[2];
   SCIP_Real facetconstant;
   SCIP_LPI* lp;
   SCIP_Bool success;

   /* build LP */
   SCIP_CALL( buildMultilinearSeparationLP(scip, lpsize, &lp) );

   /* set vars array */
   vars[0] = x;
   vars[1] = y;

   /*
    * compute a facet of the concave envelope of 1.5*x*y with x* = 0, y* = -4
    * together with the bounds (in separation.h) this should result in an estimator of the form
    *   -4.5x - 1.5y - 4.5, so the facet = -4.5, -1.5, -4.5.
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );

   SCIP_CALL( computeFacet(scip, randnumgen, lp, sol, vars, 2, 1.5, TRUE /* overestimate */,
      SCIPinfinity(scip), 1.0, (SCIP_INTERVAL){1.0, 1.0}, facetcoefs, &facetconstant, &success) );

   cr_assert(success);
   cr_expect_float_eq(facetcoefs[0], -4.5, SCIPepsilon(scip));
   cr_expect_float_eq(facetcoefs[1], -1.5, SCIPepsilon(scip));
   cr_expect_float_eq(facetconstant, -4.5, SCIPepsilon(scip));

   /*
    * compute a facet of the convex envelope for 1.5*x*y with x* = 0, y* = -4, re-using the separation lp
    * together with the bounds this should result in an estimator of the form
    *    -9x - 1.5y - 9
    */
   SCIP_CALL( SCIPsetSolVal(scip, sol, x, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, sol, y, -4.0) );

   SCIP_CALL( computeFacet(scip, randnumgen, lp, sol, vars, 2, 1.5, FALSE /* underestimate */,
      -SCIPinfinity(scip), 1.0, (SCIP_INTERVAL){1.0, 1.0}, facetcoefs, &facetconstant, &success) );

   cr_assert(success);
   cr_expect_float_eq(facetcoefs[0], -9.0, SCIPepsilon(scip));
   cr_expect_float_eq(facetcoefs[1], -1.5, SCIPepsilon(scip));
   cr_expect_float_eq(facetconstant, -9.0, SCIPepsilon(scip));

   SCIP_CALL( SCIPlpiFree(&lp) );
}

Test(separation, bilinear_with_LP, .init = setup, .fini = teardown,
   .description = "test separation for a bilinear expression, but using the LP"
   )
{
   for( int i = 2; i < 11; ++i)
      testBilinearLP(i);
}

/* This next test separates a 4-linear function and then a 3-linear function, for different sizes of the lp:
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

/* this test needs to create its own data, it doesn't use the fixtures! */
static
void testMultilinearLP(int lpsize)
{
   int i;
   SCIP_LPI* lp;
   SCIP_VAR* vars[4];
   SCIP_Real facetcoefs[4];
   SCIP_Real facetconstant;
   char const* names[] = {"x", "y", "w", "z"};
   SCIP_Real lb[] = {-0.2, -10.0, 1.0, 0.09};
   SCIP_Real ub[] = { 0.7,   8.0, 1.3,  2.1};
   SCIP_Real solval[] = { 0.2, -4.0, 1.1, 0.18};
   SCIP_Bool success;

   /* need a problem to create solutions */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );

   /* create random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randnumgen, 1, TRUE) );

   /* build separation LP */
   SCIP_CALL( buildMultilinearSeparationLP(scip, lpsize, &lp) );

   /* create and add vars and set solution value */
   SCIP_CALL( SCIPcreateSol(scip, &sol, NULL) );
   for( i = 0; i < 4; ++i )
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], names[i], lb[i], ub[i], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
      SCIP_CALL( SCIPsetSolVal(scip, sol, vars[i], solval[i]) );
   }

   /* compute an overestimator for -0.7*x*y*w*z with x* = 0.2, y* = -4, w* = 1.1, z* = 0.18
    * together with the bounds x,y,w,z \in [-0.2, 0.7], [-10, 8], [1, 1.3], [0.09, 2.1]
    * -> 63/5000 * (50x + y + 35w + 4550/9 z - 141/2)
    */
   SCIP_CALL( computeFacet(scip, randnumgen, lp, sol, vars, 4, -0.7, TRUE /* overestimate */,
            SCIPinfinity(scip), 1.0, (SCIP_INTERVAL){1.0, 1.0}, facetcoefs, &facetconstant, &success) );
   cr_assert(success);

   SCIP_Real exact_facet1[] = {63.0/100, 63.0/5000, 441.0/1000, 637.0/100, -8883.0/10000};
   for( i = 0; i < 4; ++i ) /* index 4 is the constant */
   {
      cr_expect_float_eq(facetcoefs[i], exact_facet1[i], SCIPfeastol(scip), "coef %d: received %g instead of %g\n", i, facetcoefs[i], exact_facet1[i]);
   }
   cr_expect_float_eq(facetconstant, exact_facet1[4], SCIPfeastol(scip), "constant: received %g instead of %g\n", facetconstant, exact_facet1[i]);

   /* compute an underestimator for the same function as before, but now z is fixed to 1 and we underestimate
    *   -0.7*x*y*w with x* = 0.2, y* = -4, w* = 1.1
    * together with the bounds x,y,w \in [-0.2, 0.7], [-10, 8], [1, 1.3]
    * -> -49/100 * (-100/7 x + y + 8w + 2)
    */

   /* compute cut */
   SCIP_CALL( computeFacet(scip, randnumgen, lp, sol, vars, 3, -0.7, FALSE /* underestimate */,
            -SCIPinfinity(scip), 1.0, (SCIP_INTERVAL){1.0, 1.0}, facetcoefs, &facetconstant, &success) );
   cr_assert(success);

   SCIP_Real exact_facet2[] = {7.0, -49.0/100, -98.0/25, -49.0/50};
   for( i = 0; i < 3; ++i ) /* index 3 is the constant */
   {
      cr_expect_float_eq(facetcoefs[i], exact_facet2[i], SCIPfeastol(scip), "coef %d: received %g instead of %g\n", i, facetcoefs[i], exact_facet2[i]);
   }
   cr_expect_float_eq(facetconstant, exact_facet2[3], SCIPfeastol(scip), "constant: received %g instead of %g\n", facetconstant, exact_facet2[i]);

   /* free all */
   SCIPfreeRandom(scip, &randnumgen);
   SCIP_CALL( SCIPfreeSol(scip, &sol) );
   SCIP_CALL( SCIPreleaseVar(scip, &vars[3]) );
   SCIP_CALL( SCIPreleaseVar(scip, &vars[2]) );
   SCIP_CALL( SCIPreleaseVar(scip, &vars[1]) );
   SCIP_CALL( SCIPreleaseVar(scip, &vars[0]) );
   SCIP_CALL( SCIPlpiFree(&lp) );
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(separation, multilinearseparation)
{
   for( int i = 4; i < 11; ++i)
      testMultilinearLP(i);
}

Test(separation, errorfacet)
{
   SCIP_VAR* vars[3];
   SCIP_Real maxfaceterror;
   char const* names[] = {"x", "y", "w"};
   SCIP_Real lb[] = {-0.2, -10.0, 1.0};
   SCIP_Real ub[] = { 0.7,   8.0, 1.3};
   int i;
   SCIP_Real* facet;
   SCIP_Real* funvals;

   /* a facet of the convex envelope of -0.7*x*y*w is 7 x - 133/250 y - 7/5 w - 98/25 */
   SCIP_Real facetunder[] = {7.0, -133.0/250, -7.0/5, -98.0/25};
   /* values of -0.7*x*y*w at the vertices of the domain */
   SCIP_Real funvalsunder[] = {-1.4, 4.9, 1.12, -3.92, -1.82, 6.37, 1.456, -5.096};
   /* a facet of the concave envelope of 0.7*x*y*w is -7 x + 133/250 y + 7/5 w + 98/25 */
   SCIP_Real facetover[] = {-7.0, 133.0/250, 7.0/5, 98.0/25};
   /* values of 0.7*x*y*w at the vertices of the domain */
   SCIP_Real funvalsover[] = {1.4, -4.9, -1.12, 3.92, 1.82, -6.37, -1.456, 5.096};

   /* create problem and vars */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPcreateProbBasic(scip, "test_problem") );
   for( i = 0; i < 3; ++i )
   {
      SCIP_CALL( SCIPcreateVarBasic(scip, &vars[i], names[i], lb[i], ub[i], 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, vars[i]) );
   }

   facet = facetunder;
   funvals = funvalsunder;

   /* compute the maximum error */
   printf("computing maximum error\n");
   maxfaceterror = computeMaxFacetError(scip, funvals, vars, 3, FALSE, 1.0, (SCIP_INTERVAL){1.0, 1.0}, facet, facet[3]);
   cr_expect_eq(maxfaceterror, 0.0);
   printf("done\n");

   /* perturb facet */
   facet[3] += 1.0;
   maxfaceterror = computeMaxFacetError(scip, funvals, vars, 3, FALSE, 1.0, (SCIP_INTERVAL){1.0, 1.0}, facet, facet[3]);
   cr_expect_eq(maxfaceterror, 1.0);
   facet[3] -= 1.0;

   /* now let us assume we got this function because we actually fixed a
    * variable whose bounds were -1, 3 to its middle point, 1.0.  at (0.7, -10,
    * 1.3) the function is 6.37 and the facet 4.48. However, if the fixed
    * variable would have been -1, the value of the function would be -6.37
    * giving an error of 6.37 + 4.48 = 10.85 */
   maxfaceterror = computeMaxFacetError(scip, funvals, vars, 3, FALSE, 1.0, (SCIP_INTERVAL){-1.0, 3.0}, facet, facet[3]);
   printf("maxfaceterror = %g\n", maxfaceterror);
   cr_expect_float_eq(maxfaceterror, 10.85, SCIPfeastol(scip), "difference is %g\n", maxfaceterror - 10.85);

   /* now we do the same, but overestimating */
   facet = facetover;
   funvals = funvalsover;

   /* compute the maximum error */
   printf("computing maximum error over\n");
   maxfaceterror = computeMaxFacetError(scip, funvals, vars, 3, TRUE, 1.0, (SCIP_INTERVAL){1.0, 1.0}, facet, facet[3]);
   cr_expect_eq(maxfaceterror, 0.0);
   printf("done\n");

   /* perturb facet */
   facet[3] -= 1.0;
   maxfaceterror = computeMaxFacetError(scip, funvals, vars, 3, TRUE, 1.0, (SCIP_INTERVAL){1.0, 1.0}, facet, facet[3]);
   cr_expect_eq(maxfaceterror, 1.0);
   facet[3] += 1.0;

   /* change interval */
   maxfaceterror = computeMaxFacetError(scip, funvals, vars, 3, TRUE, 1.0, (SCIP_INTERVAL){-1.0, 3.0}, facet, facet[3]);
   printf("maxfaceterror over = %g\n", maxfaceterror);
   cr_expect_float_eq(maxfaceterror, 10.85, SCIPfeastol(scip), "difference is %g\n", maxfaceterror - 10.85);

   /* free everything */
   for( i = 2; i >= 0; --i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &vars[i]) );
   }
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}
