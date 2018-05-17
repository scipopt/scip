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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   circpacking.c
 * @brief  Packing circles in a rectangle of minimal size.
 * @author Jose Salmeron
 *
 * This example shows how to setup quadratic constraints in SCIP when using SCIP as callable library.
 * The example implements a model for the computation of a smallest rectangle that contains a number of
 * given circles in the plane.
 *
 * Given n circles with radii \f$(r_i)\f$, the task is to find a coordinates \f$(x_i, y_i)\f$ for the
 * circle midpoints and a minimal rectangle \f$W, H \geq 0\f$, such that every circle is places within
 * the rectangle (\f$r_i \leq x_i \leq W-r_i\f$, \f$r_i \leq y_i \leq H-r_i\f$) and circles are not
 * overlapping (\f$\sqrt{(x_i-x_j)^2 + (y_i-y_j)^2} \geq (r_i + r_j)^2\f$).
 */

#include <stdio.h>
#include <math.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/* Model parameters */

/** Number of possible circles **/
#define nCircles 6

/** Radius **/
static const SCIP_Real r[] = {0.25, 0.25, 0.4, 0.7, 0.1, 0.3};

/** sets up problem */
static SCIP_RETCODE setupProblem(SCIP* scip) {
	SCIP_VAR* x[nCircles]; /* x coordinate */
	SCIP_VAR* y[nCircles]; /* y coordinate */
	SCIP_VAR* z;           /* aux var */
	SCIP_VAR* w;           /* height */
	SCIP_VAR* h;           /* width */

	SCIP_CONS* xrzero[nCircles];
	SCIP_CONS* xrw[nCircles];
	SCIP_CONS* yrzero[nCircles];
	SCIP_CONS* yrh[nCircles];
	SCIP_CONS* whz;
	SCIP_CONS* quad[nCircles][nCircles];

	int i, j;
	char name[SCIP_MAXSTRLEN];

	/* create empty problem */
	SCIP_CALL( SCIPcreateProbBasic(scip, "Packing circles") );

	/* create variables */
	for (i = 0; i < nCircles; i++) {
		(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
	}
	for (i = 0; i < nCircles; i++) {
		(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d", i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &y[i], name, 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
	}
	SCIP_CALL( SCIPcreateVarBasic(scip, &z, "z", 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) ); // 1*z
	SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
	SCIP_CALL( SCIPcreateVarBasic(scip, &h, "h", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

	/* add variables to problem */
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPaddVar(scip, x[i]) );
	}
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPaddVar(scip, y[i]) );
	}
	SCIP_CALL( SCIPaddVar(scip, z) );
	SCIP_CALL( SCIPaddVar(scip, w) );
	SCIP_CALL( SCIPaddVar(scip, h) );

	/* linear constraint: x_i - r_i >= 0 --> x_i >= r_i */
	{
		for (i = 0; i < nCircles; i++) {
			(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d-r_%d>=0", i, i);
			SCIP_CALL( SCIPcreateConsBasicLinear(scip, &xrzero[i], name, 0, NULL, NULL, r[i], SCIPinfinity(scip)) );
			SCIP_CALL( SCIPaddCoefLinear(scip, xrzero[i], x[i], 1.0) );
		}
	}

	/* linear constraint: x_i + r_i <= w --> r_i <= w - x_i */
	{
		for (i = 0; i < nCircles; i++) {
			(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d+r_%d>=w", i, i);
			SCIP_CALL( SCIPcreateConsBasicLinear(scip, &xrw[i], name, 0, NULL, NULL, r[i], SCIPinfinity(scip)) );
			SCIP_CALL( SCIPaddCoefLinear(scip, xrw[i], w, 1.0) );
			SCIP_CALL( SCIPaddCoefLinear(scip, xrw[i], x[i], -1.0) );
		}
	}

	/* linear constraint: y_i - r_i >= 0 --> y_i >= r_i */
	{
		for (i = 0; i < nCircles; i++) {
			(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d-r_%d>=0", i, i);
			SCIP_CALL( SCIPcreateConsBasicLinear(scip, &yrzero[i], name, 0, NULL, NULL, r[i], SCIPinfinity(scip)) );
			SCIP_CALL( SCIPaddCoefLinear(scip, yrzero[i], y[i], 1.0) );
		}
	}

	/* linear constraint: y_i + r_i <= h --> r_i <= h - y_i */
	{
		for (i = 0; i < nCircles; i++) {
			(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d+r_%d<=h", i, i);
			SCIP_CALL( SCIPcreateConsBasicLinear(scip, &yrh[i], name, 0, NULL, NULL, r[i], SCIPinfinity(scip)) );
			SCIP_CALL( SCIPaddCoefLinear(scip, yrh[i], h, 1.0) );
			SCIP_CALL( SCIPaddCoefLinear(scip, yrh[i], y[i], -1.0) );
		}
	}

	/* nonlinear constraint: w * h <= z --> w * h - z <= 0 */
	{
		SCIP_EXPR* wexpr;
		SCIP_EXPR* hexpr;
		SCIP_EXPR* expr;
		SCIP_EXPRTREE* exprtree;
		SCIP_VAR* vars[2];

		SCIP_Real one;
		SCIP_Real minusone;

		one = 1.0;
		minusone = -1.0;

		/* expression for variables w (index 0) and h (index 1) */
		SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &wexpr, SCIP_EXPR_VARIDX, 0) );
		SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &hexpr, SCIP_EXPR_VARIDX, 1) );

		/* expression for w*h */
		SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_MUL, wexpr, hexpr) );

		/* expression tree from expr */
		SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 2, 0, NULL) );

		/* set variables in expression tree */
		vars[0] = w;
		vars[1] = h;
		SCIP_CALL( SCIPexprtreeSetVars(exprtree, 2, vars) );

		/* create nonlinear constraint for exprtree - z <= 0.0 */
		SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &whz, "w*h<=z", 1, &z, &minusone, 1, &exprtree, &one, -SCIPinfinity(scip), 0.0) );
		
		/* free expression tree, because it was copied by the constraint */
		SCIP_CALL( SCIPexprtreeFree(&exprtree) );
	}

	/* quadratic constraint: (x_i - x_j)^2 + (y_i - y_j)^2 >= (r_i + r_j)^2 */
	/* x_i^2 - 2 x_i x_j + x_j^2 + y_i^2 - 2 y_i y_j + y_j^2 >= (r_i + r_j)^2 */

	/* x_i^2 - 2 x_i x_j + x_j^2 + y_i^2 - 2 y_i y_j + y_j^2 >= r_i^2 - 2 r_i r_j + r_j^2 */
	/* x_i^2 - 2 x_i x_j + x_j^2 + y_i^2 - 2 y_i y_j + y_j^2 - r_i^2 + 2 r_i r_j - r_j^2 >= 0 */
	{
		for (i = 0; i < nCircles; i++) {
			for (j = i + 1; j < nCircles; j++) {
				/* create empty quadratic constraint with right-hand-side (r_i - r_j)^2 */
				(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "quad_%d,%d", i, j);				
				SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &quad[i][j], name, 0, NULL, NULL, 0, NULL, NULL, NULL, (r[i] + r[j])*(r[i] + r[j]), SCIPinfinity(scip)) );
				/*
				SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &quad[i][j], name, 0, NULL, NULL, 0, NULL, NULL, NULL, 0, SCIPinfinity(scip)) );				
				*/
				// x_i^2
				SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], x[i], 1.0) );
				// - 2 x_i x_j
				SCIP_CALL( SCIPaddBilinTermQuadratic(scip, quad[i][j], x[i], x[j], -2.0) );
				// x_j^2
				SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], x[j], 1.0) );

				// y_i^2
				SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], y[i], 1.0) );
				// - 2 y_i y_j
				SCIP_CALL( SCIPaddBilinTermQuadratic(scip, quad[i][j], y[i], y[j], -2.0) );
				// y_j^2
				SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], y[j], 1.0) );
				/*
				// - r_i^2
				SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], r[i], -1.0) );
				// - 2 r_i r_j
				SCIP_CALL( SCIPaddBilinTermQuadratic(scip, quad[i][j], r[i], r[j], -2.0) );
				// - r_j^2
				SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], r[j], -1.0) );
				*/
			}
		}
	}

	/* add constraints to problem */
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPaddCons(scip, xrzero[i]) );
	}
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPaddCons(scip, xrw[i]) );
	}
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPaddCons(scip, yrzero[i]) );
	}
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPaddCons(scip, yrh[i]) );
	}
	SCIP_CALL( SCIPaddCons(scip, whz) );
	for (i = 0; i < nCircles; i++) {
		for (j = i + 1; j < nCircles; j++) {
			SCIP_CALL( SCIPaddCons(scip, quad[i][j]) );
		}
	}

	/* release variables and constraints
	* the problem has them captured, and we do not require them anymore
	*/
	SCIP_CALL( SCIPreleaseVar(scip, &z) );
	SCIP_CALL( SCIPreleaseVar(scip, &w) );
	SCIP_CALL( SCIPreleaseVar(scip, &h) );
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPreleaseVar(scip, &x[i]) );
	}
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
	}

	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPreleaseCons(scip, &xrzero[i]) );
	}
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPreleaseCons(scip, &xrw[i]) );
	}
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPreleaseCons(scip, &yrzero[i]) );
	}
	for (i = 0; i < nCircles; i++) {
		SCIP_CALL( SCIPreleaseCons(scip, &yrh[i]) );
	}
	SCIP_CALL( SCIPreleaseCons(scip, &whz) );
	for (i = 0; i < nCircles; i++) {
		for (j = i + 1; j < nCircles; j++) {
			SCIP_CALL( SCIPreleaseCons(scip, &quad[i][j]) );
		}
	}

	return SCIP_OKAY;
}

/* runs packing circles */
static SCIP_RETCODE runPacking(void) {
	SCIP* scip;

	SCIP_CALL( SCIPcreate(&scip) );
	SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

	SCIPinfoMessage(scip, NULL, "\n");
	SCIPinfoMessage(scip, NULL, "***************************\n");
	SCIPinfoMessage(scip, NULL, "* Running Packing Circles *\n");
	SCIPinfoMessage(scip, NULL, "***************************\n");
	SCIPinfoMessage(scip, NULL, "\n");

	SCIP_CALL( setupProblem(scip) );

	SCIPinfoMessage(scip, NULL, "Original problem:\n");
	SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

	SCIPinfoMessage(scip, NULL, "\n");
	SCIP_CALL( SCIPpresolve(scip) );

	/* SCIPinfoMessage(scip, NULL, "Reformulated problem:\n");
	SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
	*/

	SCIPinfoMessage(scip, NULL, "\nSolving...\n");
	SCIP_CALL( SCIPsolve(scip) );

	SCIP_CALL( SCIPfreeTransform(scip) );

	if(SCIPgetNSols(scip) > 0) {
		SCIPinfoMessage(scip, NULL, "\nSolution:\n");
		SCIPinfoMessage(scip, NULL, "Name: Packing Circles\n");
		SCIPinfoMessage(scip, NULL, "N %d\n", nCircles);
		SCIPinfoMessage(scip, NULL, "r ");
		for(int i = 0; i < nCircles; i++) {
			SCIPinfoMessage(scip, NULL, "%f ", r[i]);
		}
		SCIPinfoMessage(scip, NULL, "\n");
		SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
	}

	SCIP_CALL( SCIPfree(&scip) );

	return SCIP_OKAY;
}

/** main method starting SCIP */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )  /*lint --e{715}*/
{
   SCIP_RETCODE retcode;

   retcode = runPacking();

   /* evaluate return code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
