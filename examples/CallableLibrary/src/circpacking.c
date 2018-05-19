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

/* variables */
SCIP_VAR* x[nCircles]; /**< x coordinates */
SCIP_VAR* y[nCircles]; /**< y coordinates */
SCIP_VAR* a;           /**< area of rectangle */
SCIP_VAR* w;           /**< width of rectangle */
SCIP_VAR* h;           /**< height of rectangle */

/** sets up problem */
static SCIP_RETCODE setupProblem(
   SCIP*                 scip                /**< SCIP data structure */
)
{
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

	/* release constraints
	 * the problem has them captured, and we do not require them anymore
	 */
	for( i = 0; i < nCircles; ++i )
	{
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

/** plots solution by use of Python/Matplotlib */
static
void visualizeSolutionMatplotlib(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution to plot */
)
{
#if _POSIX_C_SOURCE < 2
   SCIPinfoMessage(scip, NULL, "No POSIX version 2. Try http://distrowatch.com/.");
#else
   FILE* stream;
   int i;

   stream = popen("python", "w");
   if( stream == NULL )
   {
      SCIPerrorMessage("Could not open pipe to python.\n");
      return;
   }

   fputs("import numpy as np\n"
      "import matplotlib\n"
      "import matplotlib.pyplot as plt\n",
      stream);

   fputs("fig, ax = plt.subplots()\n"
      "patches = []\n",
      stream);

   for( i = 0; i < nCircles; ++i )
   {
      fprintf(stream, "patches.append(matplotlib.patches.Circle((%g, %g), %g))\n",
         SCIPgetSolVal(scip, sol, x[i]),
         SCIPgetSolVal(scip, sol, y[i]),
         r[i]);
   }

   fputs("colors = 100*np.random.rand(len(patches))\n"
       "p = matplotlib.collections.PatchCollection(patches, alpha=0.4)\n"
       "p.set_array(np.array(colors))\n"
       "ax.add_collection(p)\n",
       stream);

   fprintf(stream, "plt.xlim(xmax=%g)\n", SCIPgetSolVal(scip, sol, w));
   fprintf(stream, "plt.ylim(ymax=%g)\n", SCIPgetSolVal(scip, sol, h));
   fprintf(stream, "plt.title('Area = %.4f')\n", SCIPgetSolVal(scip, sol, a));
   fputs("plt.gca().set_aspect(1)\n", stream);
   fputs("plt.show()\n", stream);

   pclose(stream);
#endif
}

/** plots solution by use of gnuplot */
static
void visualizeSolutionGnuplot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution to plot */
)
{
#if _POSIX_C_SOURCE < 2
   SCIPinfoMessage(scip, NULL, "No POSIX version 2. Try http://distrowatch.com/.");
#else
   SCIP_Real wval;
   SCIP_Real hval;
   FILE* stream;
   int i;

   /* -p (persist) to keep the plot open after gnuplot terminates */
   stream = popen("gnuplot -p", "w");
   if( stream == NULL )
   {
      SCIPerrorMessage("Could not open pipe to gnuplot.\n");
      return;
   }

   fputs("unset xtics\n"
      "unset ytics\n"
      "unset border\n"
      "set size ratio 1\n",
      stream);

   wval = SCIPgetSolVal(scip, sol, w);
   hval = SCIPgetSolVal(scip, sol, h);

   fprintf(stream, "set xrange [0:%.2f]\n", MAX(wval, hval));
   fprintf(stream, "set yrange [0:%.2f]\n", MAX(wval, hval));
   fprintf(stream, "set object rectangle from 0,0 to %.2f,%.2f\n", wval, hval);
   fprintf(stream, "set xlabel \"Area = %.2f\"\n", wval * hval);

   fputs("plot '-' with circles notitle\n", stream);
   for( i = 0; i < nCircles; ++i )
   {
      fprintf(stream, "%g %g %g\n",
         SCIPgetSolVal(scip, sol, x[i]),
         SCIPgetSolVal(scip, sol, y[i]),
         r[i]);
   }
   fputs("e\n", stream);

   pclose(stream);
#endif
}

/* runs packing circles */
static SCIP_RETCODE runPacking(
   SCIP_Bool             dognuplot,          /**< whether to draw best solution with gnuplot */
   SCIP_Bool             domatplotlib        /**< whether to draw best solution with python/matplotlib */
   )
{
	SCIP* scip;
	int i;

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
		for( i = 0; i < nCircles; ++i )
		{
			SCIPinfoMessage(scip, NULL, "%f ", r[i]);
		}
		SCIPinfoMessage(scip, NULL, "\n");
		SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );

		if( dognuplot )
	      visualizeSolutionGnuplot(scip, SCIPgetBestSol(scip));

		if( domatplotlib )
		   visualizeSolutionMatplotlib(scip, SCIPgetBestSol(scip));
	}

   /* release variables */
   SCIP_CALL( SCIPreleaseVar(scip, &a) );
   SCIP_CALL( SCIPreleaseVar(scip, &w) );
   SCIP_CALL( SCIPreleaseVar(scip, &h) );
   for( i = 0; i < nCircles; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &x[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
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

   retcode = runPacking(FALSE, FALSE);

   /* evaluate return code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
