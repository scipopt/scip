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
 * @author Stefan Vigerske
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
#define nCircles 5

/** Radius **/
static const SCIP_Real r[] = {0.25, 0.25, 0.4, 0.7, 0.1};

/** sets up problem */
static SCIP_RETCODE setupProblem(
   SCIP*                 scip                /**< SCIP data structure */
)
{
	SCIP_VAR* x[nCircles]; /* x coordinates */
	SCIP_VAR* y[nCircles]; /* y coordinates */
	SCIP_VAR* a;           /* area of rectangle */
	SCIP_VAR* w;           /* width of rectangle */
	SCIP_VAR* h;           /* height of rectangle */

	SCIP_CONS* xrw[nCircles];
	SCIP_CONS* yrh[nCircles];
	SCIP_CONS* wha;
	SCIP_CONS* quad[nCircles][nCircles];

	int i, j;
	char name[SCIP_MAXSTRLEN];
	SCIP_Real one = 1.0;

	/* create empty problem */
	SCIP_CALL( SCIPcreateProbBasic(scip, "Packing circles") );

	/* create variables */
	for( i = 0; i < nCircles; ++i )
	{
		(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, r[i], SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

		(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d", i);
		SCIP_CALL( SCIPcreateVarBasic(scip, &y[i], name, r[i], SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
	}
	SCIP_CALL( SCIPcreateVarBasic(scip, &a, "a", 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
	SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
	SCIP_CALL( SCIPcreateVarBasic(scip, &h, "h", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

	/* add variables to problem */
	for( i = 0; i < nCircles; ++i )
	{
		SCIP_CALL( SCIPaddVar(scip, x[i]) );
		SCIP_CALL( SCIPaddVar(scip, y[i]) );
	}
	SCIP_CALL( SCIPaddVar(scip, a) );
	SCIP_CALL( SCIPaddVar(scip, w) );
	SCIP_CALL( SCIPaddVar(scip, h) );

	/* circles must be within rectangle (left and bottom are passed in by variable bounds) */
	for( i = 0; i < nCircles; ++i )
	{
	   /* linear constraint: x_i + r_i <= w --> r_i <= w - x_i */
	   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "boundaryright_%d", i, i);
	   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &xrw[i], name, 0, NULL, NULL, r[i], SCIPinfinity(scip)) );
	   SCIP_CALL( SCIPaddCoefLinear(scip, xrw[i], w, 1.0) );
	   SCIP_CALL( SCIPaddCoefLinear(scip, xrw[i], x[i], -1.0) );

	   /* linear constraint: y_i + r_i <= h --> r_i <= h - y_i */
	   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "boundarytop_%d", i, i);
	   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &yrh[i], name, 0, NULL, NULL, r[i], SCIPinfinity(scip)) );
	   SCIP_CALL( SCIPaddCoefLinear(scip, yrh[i], h, 1.0) );
	   SCIP_CALL( SCIPaddCoefLinear(scip, yrh[i], y[i], -1.0) );
	}

	/* quadratic constraint: w * h <= a --> w * h - a <= 0 */
	SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &wha, "area", 0, NULL, NULL, 1, &w, &h, &one, -SCIPinfinity(scip), 0.0) );
	SCIP_CALL( SCIPaddLinearVarQuadratic(scip, wha, a, -1.0) );

	/* quadratic constraint: (x_i - x_j)^2 + (y_i - y_j)^2 >= (r_i + r_j)^2 */
	/* x_i^2 - 2 x_i x_j + x_j^2 + y_i^2 - 2 y_i y_j + y_j^2 >= (r_i + r_j)^2 */
	for( i = 0; i < nCircles; ++i )
	{
	   for( j = i + 1; j < nCircles; ++j )
	   {
	      /* create empty quadratic constraint with right-hand-side (r_i - r_j)^2 */
	      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nooverlap_%d,%d", i, j);
	      SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &quad[i][j], name, 0, NULL, NULL, 0, NULL, NULL, NULL, SQR(r[i] + r[j]), SCIPinfinity(scip)) );

	      SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], x[i], 1.0) ); /* x_i^2 */
         SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], x[j], 1.0) ); /* x_j^2 */
	      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, quad[i][j], x[i], x[j], -2.0) ); /* - 2 x_i x_j */

	      SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], y[i], 1.0) ); /* y_i^2 */
         SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, quad[i][j], y[j], 1.0) ); /* y_j^2 */
	      SCIP_CALL( SCIPaddBilinTermQuadratic(scip, quad[i][j], y[i], y[j], -2.0) ); /* - 2 y_i y_j */
	   }
	}

	/* add constraints to problem */
	for( i = 0; i < nCircles; ++i )
	{
		SCIP_CALL( SCIPaddCons(scip, xrw[i]) );
		SCIP_CALL( SCIPaddCons(scip, yrh[i]) );
	}

	SCIP_CALL( SCIPaddCons(scip, wha) );

	for( i = 0; i < nCircles; ++i )
	{
		for( j = i + 1; j < nCircles; ++j )
		{
			SCIP_CALL( SCIPaddCons(scip, quad[i][j]) );
		}
	}

	/* release variables and constraints
	 * the problem has them captured, and we do not require them anymore
	 */
	SCIP_CALL( SCIPreleaseVar(scip, &a) );
	SCIP_CALL( SCIPreleaseVar(scip, &w) );
	SCIP_CALL( SCIPreleaseVar(scip, &h) );
	for( i = 0; i < nCircles; ++i )
	{
		SCIP_CALL( SCIPreleaseVar(scip, &x[i]) );
		SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
		SCIP_CALL( SCIPreleaseCons(scip, &xrw[i]) );
		SCIP_CALL( SCIPreleaseCons(scip, &yrh[i]) );

      for( j = i + 1; j < nCircles; ++j )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &quad[i][j]) );
      }
	}
	SCIP_CALL( SCIPreleaseCons(scip, &wha) );

	return SCIP_OKAY;
}

/* runs packing circles */
static SCIP_RETCODE runPacking(void)
{
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

	/* closing the last bit of the gap can take very long */
	SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", 1e-4) );

	SCIPinfoMessage(scip, NULL, "\nSolving...\n");
	SCIP_CALL( SCIPsolve(scip) );

	if( SCIPgetNSols(scip) > 0 )
	{
		SCIPinfoMessage(scip, NULL, "\nSolution:\n");
		SCIPinfoMessage(scip, NULL, "Name: Packing Circles\n");
		SCIPinfoMessage(scip, NULL, "N %d\n", nCircles);
		SCIPinfoMessage(scip, NULL, "r ");
		for( int i = 0; i < nCircles; ++i )
		{
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
