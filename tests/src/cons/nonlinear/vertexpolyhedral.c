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

/**@file   vertexpolyhedral.c
 * @brief  tests estimation of vertexpolyhedral functions
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

/* ICC 2020.4.304 segfaults when compiling this file
 * __INTEL_COMPILER_BUILD_DATE=20200925
 * __INTEL_COMPILER=1910
 * __INTEL_COMPILER_UPDATE=3
 */
#if !defined(__INTEL_COMPILER) || (__INTEL_COMPILER != 1910)

#include "scip/nlhdlr.c"
#include "scip/cons_nonlinear.c"
#include "../../expr/estimation.h"

/*
 * the following test builds the matrices used to compute facets of convex/concave envelope from size 0 to 8 (otherwise
 * we get a overlength string warning) and checks it is correctly built by printing it to std output and comparing with
 * the matrices in matrices.sol
 */

#include "vertexpolyhedral_matrices.sol"

/* specify parameters of parameterized test */
ParameterizedTestParameters(separation, multilinearLP)
{
   static const int sizes[] = {1,2,3,4,5,6,7};

   /* type of the parameter; the parameter; number of parameters */
   return cr_make_param_array(const int, sizes, sizeof(sizes)/sizeof(int));
}

static
SCIP_RETCODE printMatrix(int size)
{
   int i, j, nrows, ncols;
   SCIP_LPI* lp;

   SCIP_CALL( SCIPcreate(&scip) );

   SCIP_CALL( buildVertexPolyhedralSeparationLP(scip, size, &lp) );

   SCIP_CALL( SCIPlpiGetNRows(lp, &nrows) );
   SCIP_CALL( SCIPlpiGetNCols(lp, &ncols) );

   cr_redirect_stdout();
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

   return SCIP_OKAY;
}

static
SCIP_DECL_VERTEXPOLYFUN(prodfunction)
{
   /* funcdata is a pointer to the double holding the coefficient */
   SCIP_Real ret = funcdata != NULL ? *(SCIP_Real*)funcdata : 1.0;
   int i;

   for( i = 0; i < nargs; ++i )
      ret *= args[i];

   /* these products can get fairly large when nargs is high, which is troublesome for new aspiring LP solvers;
    * let's take a sqrt() if we have more than 10 args
    */
   if( nargs > 10 )
      ret = sqrt(fabs(ret)) * (ret < 0.0 ? -1.0 : 1.0);

   return ret;
}


/* generates matrix of size *size and prints it; checks it is the expected matrix */
ParameterizedTest(const int* size, separation, multilinearLP)
{
   SCIP_CALL_ABORT( printMatrix(*size) );
}

/*
 * The following test checks that separation of multilinear terms is performed correctly.
 * It also test the re-use of the separation LP
 */
Test(separation, bilinear_with_LP, .init = setup, .fini = teardown,
   .description = "test separation for a bilinear expression, but using the code for general vertex-polyhedral functions"
   )
{
   SCIP_Real prodcoef;
   SCIP_Real xstar[2];
   SCIP_Real box[4];
   SCIP_Real facetcoefs[2];
   SCIP_Real facetconstant;
   SCIP_Bool success;

   /*
    * compute a facet of the concave envelope of 1.5*x*y with x* = 0, y* = -4
    * together with the bounds (in separation.h) this should result in an estimator of the form
    *   -4.5x - 1.5y - 4.5, so the facet = -4.5, -1.5, -4.5.
    */
   prodcoef = 1.5;

   box[0] = SCIPvarGetLbLocal(x);
   box[1] = SCIPvarGetUbLocal(x);
   box[2] = SCIPvarGetLbLocal(y);
   box[3] = SCIPvarGetUbLocal(y);

   xstar[0] = 0.0;
   xstar[1] = -4.0;

   SCIP_CALL( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, SCIPfindConshdlr(scip, "nonlinear"), TRUE /* overestimate */, prodfunction, &prodcoef, xstar, box, 2, SCIPinfinity(scip), &success, facetcoefs, &facetconstant) );

   cr_assert(success);
   cr_expect_float_eq(facetcoefs[0], -4.5, SCIPepsilon(scip));
   cr_expect_float_eq(facetcoefs[1], -1.5, SCIPepsilon(scip));
   cr_expect_float_eq(facetconstant, -4.5, SCIPepsilon(scip));

   /*
    * compute a facet of the convex envelope for 1.5*x*y with x* = 0, y* = -4, re-using the separation lp
    * together with the bounds this should result in an estimator of the form
    *    -9x - 1.5y - 9
    *
    * TODO reuse a previous LP
    */

   SCIP_CALL( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, SCIPfindConshdlr(scip, "nonlinear"), FALSE /* underestimate */, prodfunction, &prodcoef, xstar, box, 2, -SCIPinfinity(scip), &success, facetcoefs, &facetconstant) );

   cr_assert(success);
   cr_expect_float_eq(facetcoefs[0], -9.0, SCIPepsilon(scip));
   cr_expect_float_eq(facetcoefs[1], -1.5, SCIPepsilon(scip));
   cr_expect_float_eq(facetconstant, -9.0, SCIPepsilon(scip));
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

/* this test needs to create its own data, it doesn't use the fixtures! */
Test(separation, multilinearseparation)
{
   /* char const* names[] = {"x", "y", "w", "z"}; */
   SCIP_Real prodcoef = -0.7;
   SCIP_Real box[] = {-0.2, 0.7, -10.0, 8.0, 1.0, 1.3, 0.09, 2.1};
   SCIP_Real solval[] = { 0.2, -4.0, 1.1, 0.18};
   SCIP_Bool success;
   SCIP_Real facetcoefs[4];
   SCIP_Real facetconstant;
   int i;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );

   /* compute an overestimator for -0.7*x*y*w*z with x* = 0.2, y* = -4, w* = 1.1, z* = 0.18
    * together with the bounds x,y,w,z \in [-0.2, 0.7], [-10, 8], [1, 1.3], [0.09, 2.1]
    * -> 63/5000 * (50x + y + 35w + 4550/9 z - 141/2)
    */
   SCIP_CALL( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, SCIPfindConshdlr(scip, "nonlinear"), TRUE /* overestimate */, prodfunction, &prodcoef, solval, box, 4, SCIPinfinity(scip), &success, facetcoefs, &facetconstant) );

   cr_assert(success);

   SCIP_Real exact_facet1[] = {63.0/100, 63.0/5000, 441.0/1000, 637.0/100, -8883.0/10000};
   for( i = 0; i < 4; ++i ) /* index 4 is the constant */
   {
      cr_expect_float_eq(facetcoefs[i], exact_facet1[i], SCIPfeastol(scip), "coef %d: received %g instead of %g\n", i, facetcoefs[i], exact_facet1[i]);
   }
   cr_expect_float_eq(facetconstant, exact_facet1[4], SCIPfeastol(scip), "constant: received %g instead of %g\n", facetconstant, exact_facet1[i]);

   /* the code below assumes that we do the same permutations as before, so recreate scip to reset random number generator */
   SCIP_CALL( SCIPfree(&scip) );

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );

   /* compute an underestimator for the same function as before, but now z is fixed to 1 and we underestimate
    *   -0.7*x*y*w with x* = 0.2, y* = -4, w* = 1.1
    * together with the bounds x,y,w \in [-0.2, 0.7], [-10, 8], [1, 1.3]
    * -> -49/100 * (-100/7 x + y + 8w + 2)
    */
   SCIP_CALL( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, SCIPfindConshdlr(scip, "nonlinear"), FALSE /* underestimate */, prodfunction, &prodcoef, solval, box, 3, -SCIPinfinity(scip), &success, facetcoefs, &facetconstant) );

   cr_assert(success);

   SCIP_Real exact_facet2[] = {7.0, -49.0/100, -98.0/25, -49.0/50};
   for( i = 0; i < 3; ++i ) /* index 3 is the constant */
   {
      cr_expect_float_eq(facetcoefs[i], exact_facet2[i], SCIPfeastol(scip), "coef %d: received %g instead of %g\n", i, facetcoefs[i], exact_facet2[i]);
   }
   cr_expect_float_eq(facetconstant, exact_facet2[3], SCIPfeastol(scip), "constant: received %g instead of %g\n", facetconstant, exact_facet2[i]);

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

Test(separation, errorfacet)
{
   /* char const* names[] = {"x", "y", "w"}; */
   SCIP_Real box[] = {-0.2, 0.7, -10.0, 8.0, 1.0, 1.3};
   int nonfixedpos[] = { 0, 1, 2 };
   SCIP_Real* facet;
   SCIP_Real* funvals;
   SCIP_Real maxfaceterror;

   /* a facet of the convex envelope of -0.7*x*y*w is 7 x - 133/250 y - 7/5 w - 98/25 */
   SCIP_Real facetunder[] = {7.0, -133.0/250, -7.0/5, -98.0/25};
   /* values of -0.7*x*y*w at the vertices of the domain */
   SCIP_Real funvalsunder[] = {-1.4, 4.9, 1.12, -3.92, -1.82, 6.37, 1.456, -5.096};
   /* a facet of the concave envelope of 0.7*x*y*w is -7 x + 133/250 y + 7/5 w + 98/25 */
   SCIP_Real facetover[] = {-7.0, 133.0/250, 7.0/5, 98.0/25};
   /* values of 0.7*x*y*w at the vertices of the domain */
   SCIP_Real funvalsover[] = {1.4, -4.9, -1.12, 3.92, 1.82, -6.37, -1.456, 5.096};

   SCIP_CALL( SCIPcreate(&scip) );

   facet = facetunder;
   funvals = funvalsunder;

   /* compute the maximum error */
   printf("computing maximum error\n");
   maxfaceterror = computeVertexPolyhedralMaxFacetError(scip, FALSE, funvals, box, 3, 3, nonfixedpos, facet, facet[3]);
   cr_expect_eq(maxfaceterror, 0.0);
   printf("done\n");

   /* perturb facet */
   facet[3] += 1.0;
   maxfaceterror = computeVertexPolyhedralMaxFacetError(scip, FALSE, funvals, box, 3, 3, nonfixedpos, facet, facet[3]);
   cr_expect_eq(maxfaceterror, 1.0);
   facet[3] -= 1.0;

   /* now we do the same, but overestimating */
   facet = facetover;
   funvals = funvalsover;

   /* compute the maximum error */
   printf("computing maximum error over\n");
   maxfaceterror = computeVertexPolyhedralMaxFacetError(scip, TRUE, funvals, box, 3, 3, nonfixedpos, facet, facet[3]);
   cr_expect_eq(maxfaceterror, 0.0);
   printf("done\n");

   /* perturb facet */
   facet[3] -= 1.0;
   maxfaceterror = computeVertexPolyhedralMaxFacetError(scip, TRUE, funvals, box, 3, 3, nonfixedpos, facet, facet[3]);
   cr_expect_eq(maxfaceterror, 1.0);
   facet[3] += 1.0;

   /* TODO add some tests where nonfixedpos is not the trivial identity */

   /* free everything */
   SCIP_CALL( SCIPfree(&scip) );
   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

static
void test_vertexpolyhedral(
   int       dim,
   double*   box,
   double*   xstar,
   SCIP_DECL_VERTEXPOLYFUN((*function)),
   void*     functiondata,
   SCIP_Bool overestimate
)
{
   SCIP_Bool success = FALSE;
   SCIP_Real* corner;
   SCIP_Real* facetcoefs;
   SCIP_Real facetconstant;
   SCIP_Real funval;
   SCIP_Real facetval;
   SCIP_Real targetval;
   int ncorners;
   int ntight;
   int i, j;

   assert(scip != NULL);

   SCIP_ALLOC_ABORT( BMSallocMemoryArray(&facetcoefs, dim) );
   SCIP_ALLOC_ABORT( BMSallocMemoryArray(&corner, dim) );

   targetval = overestimate ? SCIPinfinity(scip) : -SCIPinfinity(scip);

   SCIP_CALL_ABORT( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, SCIPfindConshdlr(scip, "nonlinear"), overestimate,
      function, functiondata, xstar, box, dim, targetval, &success, facetcoefs, &facetconstant) );

   cr_assert(success);
   if( !success )
      return;

   ncorners = 1<<dim;
   ntight = 0;
   for( i = 0; i < ncorners; ++i )
   {
      for( j = 0; j < dim; ++j )
      {
         if( ((unsigned int)i >> j) & 0x1 )
            corner[j] = box[2 * j + 1]; /* ub of var */
         else
            corner[j] = box[2 * j]; /* lb of var */
      }

      funval = function(corner, dim, functiondata);
      facetval = facetconstant;
      for( j = 0; j < dim; ++j )
         facetval += facetcoefs[j] * corner[j];

      if( overestimate )
         cr_expect_geq(facetval, funval-SCIPdualfeastol(scip)*MAX(1.0,REALABS(funval)), "Facet value %.15g expected to be above function value %.15g", facetval, funval);
      else
         cr_expect_leq(facetval, funval+SCIPdualfeastol(scip)*MAX(1.0,REALABS(funval)), "Facet value %.15g expected to be below function value %.15g", facetval, funval);

      if( SCIPisFeasEQ(scip, facetval, funval) )
         ++ntight;
   }

   cr_assert_geq(ntight, dim+1);  /* for a facet of the envelope, the hyperplane should touch the function in at least dim+1 corner points */


   /* now set a target value */
   targetval = function(xstar, dim, functiondata);
   facetval = facetconstant;
   for( j = 0; j < dim; ++j )
      facetval += facetcoefs[j] * xstar[j];

   SCIP_CALL_ABORT( SCIPcomputeFacetVertexPolyhedralNonlinear(scip, SCIPfindConshdlr(scip, "nonlinear"), overestimate,
      function, functiondata, xstar, box, dim, targetval, &success, facetcoefs, &facetconstant) );

   /* if target couldn't be reached before, it should not have been reached now, so method should not have succeeded
    * (in principal it would also be allowed to succeed when missing the target, but current implementation doesn't) */
   if( !SCIPisFeasEQ(scip, targetval, facetval) )
      cr_assert(!success);

   BMSfreeMemoryArray(&corner);
   BMSfreeMemoryArray(&facetcoefs);
}

Test(separation, vertexpolyhedral,
   .description = "test facets of convex and concave envelopes for general vertex-polyhedral functions"
   )
{
   SCIP_Real box[2*SCIP_MAXVERTEXPOLYDIM];
   SCIP_Real xstar[SCIP_MAXVERTEXPOLYDIM];
   int dim;
   int nfixed;
   int i;

   SCIP_CALL_ABORT( SCIPcreate(&scip) );
   SCIP_CALL_ABORT( SCIPincludeConshdlrNonlinear(scip) );

   SCIP_CALL_ABORT( SCIPcreateRandom(scip, &randnumgen, 20181106, FALSE) );

   /* try every dimension */
   for( dim = 1; dim <= SCIP_MAXVERTEXPOLYDIM; ++dim )
   {
      /* make up a box, not too symmetric, and a reference point, not just in the middle */
      for( i = 0; i < dim; ++i )
      {
         box[2*i] = -i-1;
         box[2*i+1] = (i+1)/2.0;
         xstar[i] = sqrt(i) * (i%2 ? -1.0 : 1.0);
      }

      /* test that decent facets can be computed */
      printf("Dim: %d overestimate\n", dim);
      test_vertexpolyhedral(dim, box, xstar, prodfunction, NULL, TRUE);

      printf("Dim: %d underestimate\n", dim);
      test_vertexpolyhedral(dim, box, xstar, prodfunction, NULL, FALSE);

      /* fix ~50% of the variables and try again */
      nfixed = 0;
      for( i = 0; i < dim; ++i )
         if( SCIPrandomGetReal(randnumgen, 0.0, 1.0) < 0.5 )
         {
            box[2*i] = box[2*i+1];
            xstar[i] = box[2*i];
            ++nfixed;
         }

      if( nfixed == 0 || nfixed == dim )
         continue;

      printf("Dim: %d, Fixed: %d, overestimate\n", dim, nfixed);
      test_vertexpolyhedral(dim, box, xstar, prodfunction, NULL, TRUE);
   }

   SCIPfreeRandom(scip, &randnumgen);
   SCIP_CALL_ABORT( SCIPfree(&scip) );

   cr_assert_eq(BMSgetMemoryUsed(), 0, "Memory is leaking!!");
}

#endif  /* if !intel */
