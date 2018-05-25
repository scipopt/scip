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

/**@file   projection.c
 * @brief  unit test for checking the projection of points onto convex quadratic constraints
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/cons_quadratic.c"
#include "include/scip_test.h"

#define TOL 1e-6

/*
 * Gloval Variables
 */

static SCIP* scip;
static SCIP_SOL* point;
static SCIP_VAR* x;
static SCIP_VAR* y;
static SCIP_VAR* w;
static SCIP_VAR* z;
static SCIP_CONS* cons;
static SCIP_CONSDATA* consdata;
static SCIP_CONSHDLR* conshdlr;

/*
 * Local methods for tests
 */

/* stores in `matrix` the upper triangular part of the quadratic term matrix of cons */
static
void getMatrix(SCIP_Real* matrix)
{
   int i, n, row, col;
   SCIP_HASHMAP* var2index;

   n = consdata->nquadvars;

   BMSclearMemoryArray(matrix, n*n);
   SCIP_CALL( SCIPhashmapCreate(&var2index, SCIPblkmem(scip), n) );

   for( i = 0; i < n; ++i )
   {
      SCIP_CALL( SCIPhashmapInsert(var2index, consdata->quadvarterms[i].var, (void*)(size_t)i) );
      matrix[i*n + i] = consdata->quadvarterms[i].sqrcoef;
   }

   for( i = 0; i < consdata->nbilinterms; ++i )
   {
      cr_assert(SCIPhashmapExists(var2index, consdata->bilinterms[i].var1));
      cr_assert(SCIPhashmapExists(var2index, consdata->bilinterms[i].var2));

      row = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinterms[i].var1);
      col = (int)(size_t)SCIPhashmapGetImage(var2index, consdata->bilinterms[i].var2);
      if( row < col )
      {
         matrix[row * n + col] = consdata->bilinterms[i].coef/2;
      }
      else
      {
         cr_assert_neq(row, col, "there is a bilinear term with var1 == var2, this shouldn't happen by construction");
         matrix[col * n + row] = consdata->bilinterms[i].coef/2;
      }
   }

   SCIPhashmapFree(&var2index);
}

/* checks that the eigen decomposition of `cons` was computed correctly, ie, check that P is an orthonormal matrix and
 * that P D P^T = A;
 * in consdata->eigenvectors we store P^T row-wise, i.e., the first row of P^T is store in eigenvector[0..n-1], the second
 * row is stored in eigenvectors[n..2n-1], etc; equivalently, the first eigenvector is stored in eigenvector[0..n-1], the
 * second eigenvector is eigenvectors[n..2n-1], etc;
 */
static
void checkEigenDecomposition(void)
{
   int n, i, j, k;
   SCIP_Real* A;
   SCIP_Real sum;
   SCIP_Bool positive;
   SCIP_Bool firstnonzero; /* whether we have seen the first non zero eigenvalue */

   printf("SVD computed\n");

   cr_assert(consdata->isedavailable, "problem computing eigen decomposition, not available");

   n = consdata->nquadvars;

   /* check that all eigenvalues have the same sign */
   firstnonzero = FALSE;
   for( i = 0; i < n; i++ )
   {
      /* determine the sign of the first non-zero eigenvalue */
      if( SCIPisZero(scip, consdata->eigenvalues[i]) )
         continue;

      if( !firstnonzero )
      {
         positive = consdata->eigenvalues[i] > 0;
         firstnonzero = TRUE;
      }

      if( positive )
         cr_assert(SCIPisGE(scip, consdata->eigenvalues[i], 0.0),
               "eigenvalue %d is %g, expected it to be >= 0", i, consdata->eigenvalues[i]);
      else
         cr_assert(SCIPisLE(scip, consdata->eigenvalues[i], 0.0),
               "eigenvalue %d is %g, expected it to be <= 0", i, consdata->eigenvalues[i]);
   }

   /* check orthonormality of eigenvectors, ie, P^T */
   for( i = 0; i < n; i++ )
   {
      for( j = i; j < n; j++ )
      {
         sum = 0.0;

         /* compute dot product between eigenvector[i] and eigenvector[j] */
         for( k = 0; k < n; k++ )
            sum += consdata->eigenvectors[i*n + k] * consdata->eigenvectors[j*n +k];

         if( i == j )
            cr_assert_float_eq(sum, 1.0, TOL, "eigenvector %d has norm %f, expecting 1.0", i, sum);
         else
            cr_assert_float_eq(sum, 0.0, TOL, "eigenvectors %d, %d are not orthogonal, dot product is %f", i, j, sum);
      }
   }

   /* get matrix A */
   SCIP_CALL( SCIPallocBufferArray(scip, &A, n*n) );
   getMatrix(A);

   //printf("AAA:\n");
   //for( i = 0; i < n; i++ )
   //{
   //   for( j = 0; j < n; j++ )
   //      printf("%g ", A[i*n + j]);
   //   printf("\n");
   //}

   /* check that P * D * P^T = A */
   for( i = 0; i < n; i++ )
   {
      for( j = i; j < n; j++ )
      {
         sum = 0;

         /* compute (P * D * P^T)[i,j] */
         for( k = 0; k < n; k++ )
            sum += consdata->eigenvectors[k*n + i] * consdata->eigenvalues[k] * consdata->eigenvectors[k*n +j];

         cr_assert_float_eq(A[i*n + j], sum, TOL,
               "eigendecomposition is wrong, A(%d,%d) = %.10f != PDP^T(%d,%d) = %.10f\n", i, j, A[i*n+j], i, j, sum);
      }
   }

   SCIPfreeBufferArray(scip, &A);
}

/* computes and checks eigen decomposition of constraint; computes and checks projection
 * `n` is the length of the expected projection
 */
static
void checkEDandProjection(int n, SCIP_SOL* pointtoproject, SCIP_Real* expectedprojection)
{
   int i;
   SCIP_Real actualprojection[n];

   /* TODO: when we can skip test (cr_skip), remove this and just skip the test the setup.
    * this is basically a hack not to write the same code in every test */
   if( SCIPgetNNlpis(scip) == 0 )
   {
      printf("IPOPT not available, don't run test\n");
      return;
   }

   /* get consdata */
   consdata = SCIPconsGetData(cons);

   /* check coefficient of a sqr term to see whether constraint is convex or not
    * this is needed because of some asserts inside the methods we test
    * this is enough, since a matrix is positive semidefinite only if its diagonal is >= 0 */
   if( consdata->quadvarterms[0].sqrcoef > 0 )
      consdata->isconvex = TRUE;
   else
      consdata->isconcave = TRUE;

   /* add the constraints (why??) */
   //SCIP_CALL( SCIPaddCons(scip, cons) );

   /* compute eigen decomposition */
   SCIP_CALL( computeED(scip, conshdlr, cons) );

   /* check eigen decomposition */
   checkEigenDecomposition();

   /* check projection */
   SCIP_CALL( computeReferencePointProjection(scip, conshdlr, cons, point, actualprojection) );

   for( i = 0; i < n; i++ )
   {
      cr_expect_float_eq(actualprojection[i], expectedprojection[i], TOL,
            "got %.10f, expected %.10f", actualprojection[i], expectedprojection[i]);
   }
}

/** test projection */
static
void setup(void)
{
   SCIP_NLPI* nlpi;

   SCIP_CALL( SCIPcreate(&scip) );

   /* include quadratic conshdlr (need to include nonlinear) */
   SCIP_CALL( SCIPincludeConshdlrNonlinear(scip) );
   SCIP_CALL( SCIPincludeConshdlrQuadratic(scip) );

   /* include NLPI, if not available test will fail */
   SCIP_CALL( SCIPcreateNlpSolverIpopt(SCIPblkmem(scip), &nlpi) );
   if( nlpi != NULL )
   {
      SCIP_CALL( SCIPincludeNlpi(scip, nlpi) );
      SCIP_CALL( SCIPincludeExternalCodeInformation(scip, SCIPgetSolverNameIpopt(), SCIPgetSolverDescIpopt()) );
   }

   /* create a problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "problem") );

   /* create variables */
   SCIP_CALL( SCIPcreateVarBasic(scip, &x, "x", -100, 100, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &y, "y", -100, 100, 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", -100, 100, 1.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", -100, 100, 0.0, SCIP_VARTYPE_CONTINUOUS) );

   SCIP_CALL( SCIPaddVar(scip, x) );
   SCIP_CALL( SCIPaddVar(scip, y) );
   SCIP_CALL( SCIPaddVar(scip, w) );
   SCIP_CALL( SCIPaddVar(scip, z) );

   /* create solution (point to project) */
   SCIP_CALL( SCIPcreateSol(scip, &point, NULL) );

   /* find quadratic conshdlr, we need it to call computeED and computeReferencePointProjection */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME); /* we are including cons_quadratic.c */
   cr_assert_not_null(conshdlr);
}

static
void teardown(void)
{
   /* TODO: when we can skip test (cr_skip), remove this and just skip the test in setup. */
   if( SCIPgetNNlpis(scip) == 0 )
   {
      printf("IPOPT not available, don't run test\n");
      return;
   }

   /* release cons */
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* release vars */
   SCIP_CALL( SCIPreleaseVar(scip, &x) );
   SCIP_CALL( SCIPreleaseVar(scip, &y) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &z) );

   /* free sol */
   SCIP_CALL( SCIPfreeSol(scip, &point) );

   /* free SCIP */
   SCIP_CALL( SCIPfree(&scip) );

   /* check for memory leaks */
   cr_assert_eq(BMSgetMemoryUsed(), 0, "There are memory leaks!");
}

/*
 * TESTSUITE
 */
TestSuite(separation, .init = setup, .fini = teardown);

/*
 * TESTS
 */
/* test projection over: x^2 + y^2 <= 1 */
Test(separation, projection_simple)
{
   enum nquadterms {nquadterms = 2};

   SCIP_VAR* quadvars1[nquadterms];
   SCIP_VAR* quadvars2[nquadterms];
   SCIP_Real quadvals[nquadterms]  = {1.0, 1.0};

   SCIP_Real lhs = -SCIPinfinity(scip);
   SCIP_Real rhs = 1.0;

   quadvars1[0] = x; quadvars2[0] = x;
   quadvars1[1] = y; quadvars2[1] = y;

   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "cons",
            0, NULL, NULL, nquadterms, quadvars1, quadvars2, quadvals, lhs, rhs) );

   enum n {n = 2};
   SCIP_Real expectedprojection[n] = {1.0, 0.0};

   /* set point to project */
   SCIP_CALL( SCIPsetSolVal(scip, point, x, 5.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, y, 0.0) );

   /* check projection and eigen decomposition */
   checkEDandProjection(n, point, expectedprojection);
}

/* test projection over:
 * x^2 + 2829.0/2000*y^2  + 5266.0/5321*y*w + 4166.0/3319*y*z + 441.0/1033*w^2 + 2480.0/4413*w*z + 2093.0/1805*z^2 <= 1
 * @note: strictly convex quadratic function, no variable appears linearly.
 */
Test(separation, projection_complex_strictly_convex)
{
   enum nquadterms {nquadterms = 7};

   SCIP_VAR* quadvars1[nquadterms];
   SCIP_VAR* quadvars2[nquadterms];
   SCIP_Real quadvals[nquadterms]  = {2829.0/2000, 5266.0/5321, 4166.0/3319, 441.0/1033, 2480.0/4413, 2093.0/1805, 1.0};

   SCIP_Real lhs = -SCIPinfinity(scip);
   SCIP_Real rhs = 1.0;

   quadvars1[0] = y; quadvars2[0] = y;
   quadvars1[1] = y; quadvars2[1] = w;
   quadvars1[2] = y; quadvars2[2] = z;
   quadvars1[3] = w; quadvars2[3] = w;
   quadvars1[4] = w; quadvars2[4] = z;
   quadvars1[5] = z; quadvars2[5] = z;
   quadvars1[6] = x; quadvars2[6] = x;

   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "cons",
            0, NULL, NULL, nquadterms, quadvars1, quadvars2, quadvals, lhs, rhs) );


   enum n {n = 4};
   SCIP_Real expectedprojection[n] = {0.904237367218059, -0.375667504869449, -0.540802541465134, 0.529035034426650};

   /* set point to project */
   SCIP_CALL( SCIPsetSolVal(scip, point, x, 2.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, y, 3.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, z, -1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, w, 0.0) );

   /* check projection and eigen decomposition */
   checkEDandProjection(n, point, expectedprojection);
}

/* test projection over:
 * x + 2829.0/2000*y^2  + 5266.0/5321*y*w + 4166.0/3319*y*z + 441.0/1033*w^2 + 2480.0/4413*w*z + 2093.0/1805*z^2 <= 1
 * @note x doesn't appear quadratically
 */
Test(separation, projection_complex_with_linear_part)
{
   enum nquadterms {nquadterms = 6};
   enum nlinvars {nlinvars = 1};

   SCIP_VAR* quadvars1[nquadterms];
   SCIP_VAR* quadvars2[nquadterms];
   SCIP_Real quadvals[nquadterms]  = {2829.0/2000, 5266.0/5321, 4166.0/3319, 441.0/1033, 2480.0/4413, 2093.0/1805};

   SCIP_VAR* linvars[nlinvars];
   SCIP_Real linvals[nlinvars] = {1.0};

   SCIP_Real lhs = -SCIPinfinity(scip);
   SCIP_Real rhs = 1.0;

   quadvars1[0] = y; quadvars2[0] = y;
   quadvars1[1] = y; quadvars2[1] = w;
   quadvars1[2] = y; quadvars2[2] = z;
   quadvars1[3] = w; quadvars2[3] = w;
   quadvars1[4] = w; quadvars2[4] = z;
   quadvars1[5] = z; quadvars2[5] = z;

   linvars[0] = x;

   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "cons",
            nlinvars, linvars, linvals, nquadterms, quadvars1, quadvars2, quadvals, lhs, rhs) );


   enum n {n = 4};
   SCIP_Real expectedprojection[n-1] = {0.190585198986541, 0.543192501726410, 0.305475137359772};

   /* set point to project */
   SCIP_CALL( SCIPsetSolVal(scip, point, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, y, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, w, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, z, 1.0) );

   /* we only check the value for the non-linear variables. This is not ideal, but those values are
    * the only ones that matter when computing a gradient cut, so are the only ones that SCIP computes
    * for this reason we give n-1 as the number of variables (there are n-1 nonlinear vars) */
   checkEDandProjection(n-1, point, expectedprojection);
}

/* test projection over:
 * 9*x^2 + 10*x*y + 26*x*w + 20*x*z + 5*y^2 + 10*y*w + 20*y*z + 21*w^2 + 20*w*z + 20*z^2 - x + 2*y - w <= 30
 * @note: quadratic variables with linear terms (ie, x^2 + x) appear
 */
Test(separation, projection_complex_with_linear_terms)
{
   enum nquadterms {nquadterms = 10};
   enum nlinvars {nlinvars = 3};

   SCIP_VAR* quadvars1[nquadterms];
   SCIP_VAR* quadvars2[nquadterms];
   SCIP_Real quadvals[nquadterms]  = {9,10,26,20,5,10,20,21,20,20};

   SCIP_VAR* linvars[nlinvars];
   SCIP_Real linvals[nlinvars] = {-1.0,2.0,-1.0};

   SCIP_Real lhs = -SCIPinfinity(scip);
   SCIP_Real rhs = 30.0;

   quadvars1[0] = x; quadvars2[0] = x;
   quadvars1[1] = x; quadvars2[1] = y;
   quadvars1[2] = x; quadvars2[2] = w;
   quadvars1[3] = x; quadvars2[3] = z;
   quadvars1[4] = y; quadvars2[4] = y;
   quadvars1[5] = y; quadvars2[5] = w;
   quadvars1[6] = y; quadvars2[6] = z;
   quadvars1[7] = w; quadvars2[7] = w;
   quadvars1[8] = w; quadvars2[8] = z;
   quadvars1[9] = z; quadvars2[9] = z;

   linvars[0] = x;
   linvars[1] = y;
   linvars[2] = w;

   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "cons",
            nlinvars, linvars, linvals, nquadterms, quadvars1, quadvars2, quadvals, lhs, rhs) );


   enum n {n = 4};
   SCIP_Real expectedprojection[n] = {0.520865712781762, 0.630362637727899, 0.364426436908240, 0.323315076692312};

   /* set point to project */
   SCIP_CALL( SCIPsetSolVal(scip, point, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, y, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, w, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, z, 1.0) );

   /* check projection and eigen decomposition */
   checkEDandProjection(n, point, expectedprojection);
}

/* test projection over:
 * - x^2 - 2829.0/2000*y^2  - 5266.0/5321*y*w - 4166.0/3319*y*z - 441.0/1033*w^2 - 2480.0/4413*w*z - 2093.0/1805*z^2 >= -1
 * @note: convex quadratic constraint, but function is now concave
 */
Test(separation, projection_complex_strictly_concave)
{
   enum nquadterms {nquadterms = 7};

   SCIP_VAR* quadvars1[nquadterms];
   SCIP_VAR* quadvars2[nquadterms];
   SCIP_Real quadvals[nquadterms]  = {-2829.0/2000, -5266.0/5321, -4166.0/3319,
      -441.0/1033, -2480.0/4413, -2093.0/1805, -1.0};

   SCIP_Real lhs = -1.0;
   SCIP_Real rhs = SCIPinfinity(scip);

   quadvars1[0] = y; quadvars2[0] = y;
   quadvars1[1] = y; quadvars2[1] = w;
   quadvars1[2] = y; quadvars2[2] = z;
   quadvars1[3] = w; quadvars2[3] = w;
   quadvars1[4] = w; quadvars2[4] = z;
   quadvars1[5] = z; quadvars2[5] = z;
   quadvars1[6] = x; quadvars2[6] = x;

   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "cons",
            0, NULL, NULL, nquadterms, quadvars1, quadvars2, quadvals, lhs, rhs) );


   enum n {n = 4};
   SCIP_Real expectedprojection[n] = {0.904237367218059, -0.375667504869449, -0.540802541465134, 0.529035034426650};

   /* set point to project */
   SCIP_CALL( SCIPsetSolVal(scip, point, y, 3.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, w, 0.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, z, -1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, x, 2.0) );

   /* check projection and eigen decomposition */
   checkEDandProjection(n, point, expectedprojection);
}

/* test projection over:
 * 999.999001*x^2 -0.00199800000007144*x*y -1.99799800199401*x*w + 999.999001*y^2
 * -1.99799800199451*y*w + 1.00199799999999*w^2 -600*x + 200*w + 12*z <= -20
 * @note: quadratic with quadratic variables with linear terms (ie, x^2 + x) and also a linear part (12 *z)
 */
Test(separation, projection_complex_with_linear_part_and_terms)
{
   enum nquadterms {nquadterms = 6};
   enum nlinvars {nlinvars = 3};

   SCIP_VAR* quadvars1[nquadterms];
   SCIP_VAR* quadvars2[nquadterms];
   SCIP_Real quadvals[nquadterms]  = {999.999001, -0.00199800000007144, -1.99799800199401,
      999.999001, -1.99799800199451, 1.00199799999999};

   SCIP_VAR* linvars[nlinvars];
   SCIP_Real linvals[nlinvars] = {-600.0,200.0,12.0};

   SCIP_Real lhs = -SCIPinfinity(scip);
   SCIP_Real rhs = -20.0;

   quadvars1[0] = x; quadvars2[0] = x;
   quadvars1[1] = x; quadvars2[1] = y;
   quadvars1[2] = x; quadvars2[2] = w;
   quadvars1[3] = y; quadvars2[3] = y;
   quadvars1[4] = y; quadvars2[4] = w;
   quadvars1[5] = w; quadvars2[5] = w;

   linvars[0] = x;
   linvars[1] = w;
   linvars[2] = z;

   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "cons",
            nlinvars, linvars, linvals, nquadterms, quadvars1, quadvars2, quadvals, lhs, rhs) );


   enum n {n = 4};
   SCIP_Real expectedprojection[n-1] = {0.377889429089707, 0.111193729191078, 0.201511107483149};

   /* set point to project */
   SCIP_CALL( SCIPsetSolVal(scip, point, x, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, y, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, w, 1.0) );
   SCIP_CALL( SCIPsetSolVal(scip, point, z, 1.0) );

   /* we only check the value for the non-linear variables. This is not ideal, but those values are
    * the only ones that matter when computing a gradient cut, so are the only ones that SCIP computes
    * for this reason we give n-1 as the number of variables (there are n-1 nonlinear vars) */
   checkEDandProjection(n-1, point, expectedprojection);
}
