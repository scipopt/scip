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

/**@file   circlepacking.c
 * @brief  Packing circles in a rectangle of minimal size.
 * @author Jose Salmeron
 * @author Stefan Vigerske
 *
 * This example shows how to setup quadratic constraints in SCIP when using SCIP as callable library.
 * The example implements a model for the computation of a smallest rectangle that contains a number of
 * given circles in the plane or the computation of the maximal number of circles that can be placed
 * into a given rectangle.
 *
 * Suppose that n circles with radii \f$r_i\f$ are given.
 * The task is to find coordinates \f$(x_i, y_i)\f$ for the circle midpoints and a rectangle of
 * width \f$W \geq 0\f$ and height \f$H \geq 0\f$, such that
 * - every circle is placed within the rectangle (\f$r_i \leq x_i \leq W-r_i\f$, \f$r_i \leq y_i \leq H-r_i\f$),
 * - circles are not overlapping \f$\left((x_i-x_j)^2 + (y_i-y_j)^2 \geq (r_i + r_j)^2\right)\f$, and
 * - the area of the rectangle is minimal.
 *
 * Alternatively, one may fix the size of the rectangle and maximize the number of circles that
 * can be fit into the rectangle without being overlapping.
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/* Model parameters */

/** Number of possible circles **/
int ncircles = 0;

/** Radii **/
SCIP_Real* r = NULL;
int rsize = 0;

/* Variables */
SCIP_VAR** x;       /**< x coordinates */
SCIP_VAR** y;       /**< y coordinates */
SCIP_VAR** b;       /**< whether circle is placed into rectangle */
SCIP_VAR*  a;       /**< area of rectangle */
SCIP_VAR*  w;       /**< width of rectangle */
SCIP_VAR*  h;       /**< height of rectangle */
SCIP_Bool  minarea; /**< whether we minimize the area (TRUE) or maximize the number of circles in the rectangle (FALSE) */


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

   for( i = 0; i < ncircles; ++i )
   {
      /* The binary variable b[i] indicates which circles should be included in the current solution.
       * Here we want to skip circles that are not included, that is b[i] is zero (or close to zero due to tolerances).
       */
      if( !minarea && SCIPgetSolVal(scip, sol, b[i]) < 0.5 )
         continue;

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
   if( minarea )
      fprintf(stream, "plt.title('Area = %.4f')\n", SCIPgetSolVal(scip, sol, a));
   else
      fprintf(stream, "plt.title('Number of circles = %d')\n", (int)SCIPround(scip, SCIPgetSolOrigObj(scip, sol)));
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
   if( minarea )
      fprintf(stream, "set xlabel \"Area = %.4f\"\n", SCIPgetSolVal(scip, sol, a));
   else
      fprintf(stream, "set xlabel \"Number of circles = %d\"\n", (int)SCIPround(scip, SCIPgetSolOrigObj(scip, sol)));

   fputs("plot '-' with circles notitle\n", stream);
   for( i = 0; i < ncircles; ++i )
   {
      /* skip circles that are not included in current solution, that is b[i] is close to zero */
      if( !minarea && SCIPgetSolVal(scip, sol, b[i]) < 0.5 )
         continue;

      fprintf(stream, "%g %g %g\n",
         SCIPgetSolVal(scip, sol, x[i]),
         SCIPgetSolVal(scip, sol, y[i]),
         r[i]);
   }
   fputs("e\n", stream);

   pclose(stream);
#endif
}

/** plots solution by use of ascii graphics */
static
SCIP_RETCODE visualizeSolutionAscii(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< solution to plot */
)
{
   SCIP_Real wval;
   SCIP_Real hval;
   SCIP_Real xval;
   SCIP_Real yval;
   SCIP_Real radius;
   SCIP_Real scalex;
   SCIP_Real scaley;
   int dispwidth;
   int width;
   int height;
   char* picture;
   int i;

   wval = SCIPgetSolVal(scip, sol, w);
   hval = SCIPgetSolVal(scip, sol, h);

   /* scale so that picture is about as wide as SCIP B&B log */
   SCIP_CALL( SCIPgetIntParam(scip, "display/width", &dispwidth) );
   scalex = (dispwidth-3) / wval;
   scaley = scalex / 2.0;  /* this gives almost round looking circles on my (SV) terminal */

   width = SCIPceil(scip, scalex*wval)+3;  /* +2 for left and right border and +1 for \n */
   height = SCIPceil(scip, scaley*hval)+2; /* +2 for top and bottom border */

   SCIP_CALL( SCIPallocBufferArray(scip, &picture, width * height + 1) );

   /* initialize with white space and termination */
   memset(picture, ' ', width * height);
   picture[width*height] = '\0';

   /* add border and linebreaks */
   memset(picture, '*', width-1); /* top border */
   memset(picture + (height-1) * width, '*', width-1);  /* bottom border */
   for( i = 0; i < height; ++i )
   {
      picture[i*width] = '*';  /* left border */
      picture[i*width + width-2] = '*';  /* right border */
      picture[i*width + width-1] = '\n';
   }

   /* add circles */
   for( i = 0; i < ncircles; ++i )
   {
      SCIP_Real phi;
      int xcoord;
      int ycoord;

      /* skip circles that are not included in current solution, that is b[i] close to zero */
      if( !minarea && SCIPgetSolVal(scip, sol, b[i]) < 0.5 )
         continue;

      xval = SCIPgetSolVal(scip, sol, x[i]);
      yval = SCIPgetSolVal(scip, sol, y[i]);
      radius = r[i];

      for( phi = 0.0; phi < 6.283 /* 2*pi */; phi += 0.01 )
      {
         xcoord = SCIPround(scip, scalex * (xval + radius * cos(phi))) + 1; /* +1 for border */
         ycoord = SCIPround(scip, scaley * (yval + radius * sin(phi))) + 1; /* +1 for border */

         /* feasible solutions should be within box
          * due to rounding, they can be on the border
          */
         assert(xcoord >= 0);
         assert(ycoord >= 0);
         assert(xcoord < width);
         assert(ycoord < height);

         picture[ycoord * width + xcoord] = 'a' + i;
      }
   }

   /* print objective value inside top border */
   i = SCIPsnprintf(picture + width/2 - 8, width/2 + 8,
      minarea ? " Area = %g " : " #Circles = %.0f ", SCIPgetSolOrigObj(scip, sol));
   picture[width/2-8+i] = '*';

   /* show plot */
   SCIPinfoMessage(scip, NULL, picture);

   SCIPfreeBufferArray(scip, &picture);

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitDispsol)
{
   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitDispsol)
{
   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecDispsol)
{  /*lint --e{715}*/
   SCIP_SOL* sol;

   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND);

   sol = SCIPeventGetSol(event);
   assert(sol != NULL);

   SCIP_CALL( visualizeSolutionAscii(scip, sol) );

   return SCIP_OKAY;
}

/** creates event handler for dispsol event */
static
SCIP_RETCODE includeEventHdlrDispsol(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_EVENTHDLR* eventhdlr = NULL;

   /* include event handler into SCIP */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, "dispsol", "displays new solutions",
      eventExecDispsol, NULL) );
   assert(eventhdlr != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitDispsol) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitDispsol) );

   return SCIP_OKAY;
}

/** create problem in given SCIP and add all variables and constraints to model the circle packing problem */
static SCIP_RETCODE setupProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             fixwidth,           /**< a given fixed width for the rectangle, or SCIP_INVALID if not fixed */
   SCIP_Real             fixheight           /**< a given fixed height for the rectangle, or SCIP_INVALID if not fixed */
)
{
   char name[SCIP_MAXSTRLEN];
   SCIP_CONS* cons;
   SCIP_Real one = 1.0;
   int i, j;

   /* if both width and height are fixed, then we don't optimize the area anymore, but the number of circles */
   minarea = (fixwidth == SCIP_INVALID || fixheight == SCIP_INVALID);

   /* create empty problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "Packing circles") );

   /* change to maximization if optimizing number of circles instead of rectangle area */
   if( !minarea )
   {
      SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MAXIMIZE) );
   }

   /* Create variables, setup variable bounds, and setup objective function:
    * We add variables x[i], y[i] for the circle midpoints and, if optimizing the number of circles in the rectangle,
    * a binary variable b[i] to indicate whether circle i should be in the rectangle or not.
    * Additionally, we add variables for the rectangle area, width, and height.
    *
    * Since the rectangle lower-left corner is assumed to be at (0,0), we can set a lower bound
    * r[i] for both variables x[i] and y[i].
    *
    * In the objective function, we have 1.0*a, if minimizing the area of the rectangle,
    * otherwise the area does not contribute to the objective, but every b[i] will be present instead.
    *
    * Further below we fix the width and height of the rectangle, if fixwidth and fixheight are valid.
    * If both are valid, then we can also fix the area of the rectangle.
    */
   SCIP_CALL( SCIPallocMemoryArray(scip, &x, ncircles) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &y, ncircles) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &b, ncircles) );
   for( i = 0; i < ncircles; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, r[i], SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &y[i], name, r[i], SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, y[i]) );

      if( !minarea )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "b_%d", i);
         SCIP_CALL( SCIPcreateVarBasic(scip, &b[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
         SCIP_CALL( SCIPaddVar(scip, b[i]) );
      }
   }
   SCIP_CALL( SCIPcreateVarBasic(scip, &a, "a", 0.0, SCIPinfinity(scip), minarea ? 1.0 : 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, a) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &w, "w", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, w) );

   SCIP_CALL( SCIPcreateVarBasic(scip, &h, "h", 0.0, SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );
   SCIP_CALL( SCIPaddVar(scip, h) );

   /* fix width if a valid value for fixwidth is given */
   if( fixwidth != SCIP_INVALID )
   {
      SCIP_Bool infeas;
      SCIP_Bool fixed;

      SCIP_CALL( SCIPfixVar(scip, w, fixwidth, &infeas, &fixed) );

      assert(!infeas);
      assert(fixed);
   }

   /* fix height if a valid value for fixheight is given */
   if( fixheight != SCIP_INVALID )
   {
      SCIP_Bool infeas;
      SCIP_Bool fixed;

      SCIP_CALL( SCIPfixVar(scip, h, fixheight, &infeas, &fixed) );

      assert(!infeas);
      assert(fixed);
   }

   /* fix area if both width and height are fixed */
   if( !minarea )
   {
      SCIP_Bool infeas;
      SCIP_Bool fixed;

      SCIP_CALL( SCIPfixVar(scip, a, fixwidth * fixheight, &infeas, &fixed) );

      assert(!infeas);
      assert(fixed);
   }

   /* boundary constraints on circle coordinates
    *
    * If minimizing the area of the rectangle, then the coordinates of every circle must
    * satisfy the boundary conditions r_i <= x_i <= w - r_i and r_i <= y_i <= h - r_i.
    * The lower bounds r_i are already required by the variable bounds (see above).
    * For the upper bounds, we add the corresponding linear constraints.
    *
    * If the area of the rectangle is fixed, then it would be sufficient to place only
    * circles into the rectangle that have been decided to be put into the rectangle.
    * We could model this as a big-M constraint on x_i and y_i, but on the other hand,
    * we can also require that the circle coordinates always satisfy the boundary conditions,
    * even if the circle is not placed into the rectangle (b_i=0).
    * As the non-overlapping constraints do not apply for circles that are not placed into
    * the rectangle, satisfying these boundary conditions is always possible, unless the
    * circle itself is too big to be placed into the rectangle. In this case, though,
    * we can decide a-priori that the circle is not placed into the rectangle, i.e., fix b_i to 0.
    */
   for( i = 0; i < ncircles; ++i )
   {
      if( !minarea && SCIPisLT(scip, MIN(fixwidth, fixheight), 2*r[i]) )
      {
         SCIP_Bool infeas;
         SCIP_Bool fixed;

         SCIP_CALL( SCIPfixVar(scip, b[i], 0.0, &infeas, &fixed) );

         assert(!infeas);
         assert(fixed);

         continue;
      }

      /* linear constraint: x_i + r_i <= w --> r_i <= w - x_i */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "boundaryright_%d", i, i);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, r[i], SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, w, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i], -1.0) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );

      /* linear constraint: y_i + r_i <= h --> r_i <= h - y_i */
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "boundarytop_%d", i, i);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, r[i], SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, h, 1.0) );
      SCIP_CALL( SCIPaddCoefLinear(scip, cons, y[i], -1.0) );
      SCIP_CALL( SCIPaddCons(scip, cons) );
      SCIP_CALL( SCIPreleaseCons(scip, &cons) );
   }

   /* constraint that defines the area of the rectangle
    *
    * We could add the quadratic constraint w * h - a = 0.
    * But if we are minimizing a, then we can relax this constraint to w * h - a <= 0.
    * If the size the rectangle is fixed, then w, h, and a have been fixed above.
    * We could actually omit this constraint, but here SCIP presolve will take care of removing it.
    */
   SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, "area", 0, NULL, NULL, 1, &w, &h, &one, -SCIPinfinity(scip), 0.0) );
   SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, a, -1.0) );
   SCIP_CALL( SCIPaddCons(scip, cons) );
   SCIP_CALL( SCIPreleaseCons(scip, &cons) );

   /* quadratic constraints that require that circles are not overlapping.
    * For circles i and j, i<j, we require that the euclidean distance of the circle middle points
    * is at least the sum of the circle radii, i.e., || (x_i,y_i) - (x_j,y_j) || >= r_i + r_j.
    * Equivalently, (x_i - x_j)^2 + (y_i - y_j)^2 >= (r_i + r_j)^2, which can be expanded to
    *   x_i^2 - 2 x_i x_j + x_j^2 + y_i^2 - 2 y_i y_j + y_j^2 >= (r_i + r_j)^2
    *
    * When not minimizing the area of the circles, however, then this constraint only needs
    * to hold if both circles are placed into the rectangle, that is if b_i=1 and b_j=1.
    * We can achieve this by relaxing the right-hand-side to 0 or a negative value if b_i + b_j <= 1:
    *   (x_i - x_j)^2 + (y_i - y_j)^2 >= (r_i + r_j)^2 * (b_i+b_j-1), which can be expanded to
    *   x_i^2 - 2 x_i x_j + x_j^2 + y_i^2 - 2 y_i y_j + y_j^2 - (r_i+r_j)^2 b_i - (r_i+r_j)^2 b_j >= -(r_i + r_j)^2
    */
   for( i = 0; i < ncircles; ++i )
   {
      for( j = i + 1; j < ncircles; ++j )
      {
         /* create empty quadratic constraint with right-hand-side +/- (r_i - r_j)^2 */
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "nooverlap_%d,%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicQuadratic(scip, &cons, name, 0, NULL, NULL, 0, NULL, NULL, NULL, (minarea ? 1.0 : -1.0) * SQR(r[i] + r[j]), SCIPinfinity(scip)) );

         SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, x[i], 1.0) ); /* x_i^2 */
         SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, x[j], 1.0) ); /* x_j^2 */
         SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, x[i], x[j], -2.0) ); /* - 2 x_i x_j */

         SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, y[i], 1.0) ); /* y_i^2 */
         SCIP_CALL( SCIPaddSquareCoefQuadratic(scip, cons, y[j], 1.0) ); /* y_j^2 */
         SCIP_CALL( SCIPaddBilinTermQuadratic(scip, cons, y[i], y[j], -2.0) ); /* - 2 y_i y_j */

         if( !minarea )
         {
            /* add -(r_i+r_j)^2 (b_i + b_j) */
            SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, b[i], -SQR(r[i] + r[j])) );
            SCIP_CALL( SCIPaddLinearVarQuadratic(scip, cons, b[j], -SQR(r[i] + r[j])) );
         }

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   return SCIP_OKAY;
}

/* run circle packing example
 *
 * Sets up SCIP and the SCIP problem, solves the problem, and shows the solution.
 */
static SCIP_RETCODE runPacking(
   SCIP_Real             fixwidth,           /**< a given fixed width for the rectangle, or SCIP_INVALID if not fixed */
   SCIP_Real             fixheight,          /**< a given fixed height for the rectangle, or SCIP_INVALID if not fixed */
   SCIP_Bool             dognuplot,          /**< whether to draw best solution with gnuplot */
   SCIP_Bool             domatplotlib        /**< whether to draw best solution with python/matplotlib */
)
{
   SCIP* scip;
   int i;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( includeEventHdlrDispsol(scip) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "***************************\n");
   SCIPinfoMessage(scip, NULL, "* Running Packing Circles *\n");
   SCIPinfoMessage(scip, NULL, "***************************\n");
   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "%d circles given with radii", ncircles);
   for( i = 0; i < ncircles; ++i )
      SCIPinfoMessage(scip, NULL, " %.2f", r[i]);
   SCIPinfoMessage(scip, NULL, "\n\n");

   SCIP_CALL( setupProblem(scip, fixwidth, fixheight) );

   SCIPinfoMessage(scip, NULL, "Original problem:\n");
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

   /* By default, SCIP tries to close the gap between primal and dual bound completely.
    * This can take very long for this example, so we increase the gap tolerance to have
    * SCIP stop when the distance between primal and dual bound is already below 1e-4.
    */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", 1e-4) );

   SCIPinfoMessage(scip, NULL, "\nSolving...\n");
   SCIP_CALL( SCIPsolve(scip) );

   if( SCIPgetNSols(scip) > 0 )
   {
      SCIPinfoMessage(scip, NULL, "\nSolution:\n");
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
   for( i = 0; i < ncircles; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &x[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
      if( !minarea )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &b[i]) );
      }
   }

   /* free memory arrays */
   SCIPfreeMemoryArray(scip, &b);
   SCIPfreeMemoryArray(scip, &y);
   SCIPfreeMemoryArray(scip, &x);

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
   SCIP_Bool dognuplot = FALSE;
   SCIP_Bool domatplotlib = FALSE;
   SCIP_Real fixwidth = SCIP_INVALID;
   SCIP_Real fixheight = SCIP_INVALID;
   char* endptr;
   int i;

   ncircles = 0;
   rsize = 5;
   SCIP_ALLOC_ABORT( BMSallocMemoryArray(&r, rsize) );

   for( i = 1; i < argc; ++i )
   {
      if( strcmp(argv[i], "--help") == 0 )
      {
         printf("usage: %s [--help] [-w <width>] [-h <height>]", argv[0]);
#if _POSIX_C_SOURCE >= 2
         printf(" [-g] [-m]");
#endif
         puts(" { <radius> }");
         puts("  --help shows this help and exits");
         puts("  -w <width> fix rectangle width to given value");
         puts("  -h <height> fix rectangle height to given value");
#if _POSIX_C_SOURCE >= 2
         puts("  -g show final solution with gnuplot");
         puts("  -m show final solution with matplotlib");
#endif
         puts("If no radii are given, then 5 circles with radii 0.25, 0.25, 0.4, 0.7, and 0.1 used.");
         puts("If both width and height are fixed, then the number of circles that fit into the rectangle is maximized.");

         return EXIT_SUCCESS;
      }

#if _POSIX_C_SOURCE >= 2
      if( strcmp(argv[i], "-g") == 0 )
      {
         dognuplot = TRUE;
         continue;
      }

      if( strcmp(argv[i], "-m") == 0 )
      {
         domatplotlib = TRUE;
         continue;
      }
#endif

      if( strcmp(argv[i], "-w") == 0 )
      {
         if( i == argc-1 )
         {
            fprintf(stderr, "ERROR: Missing argument for -w.\n");
            return EXIT_FAILURE;
         }

         fixwidth = strtod(argv[i+1], &endptr);
         if( *endptr != '\0' )
         {
            fprintf(stderr, "ERROR: Could not parse argument %s into a number.\n", argv[i+1]);
            return EXIT_FAILURE;
         }

         ++i;
         continue;
      }

      if( strcmp(argv[i], "-h") == 0 )
      {
         if( i == argc-1 )
         {
            fprintf(stderr, "ERROR: Missing argument for -h.\n");
            return EXIT_FAILURE;
         }

         fixheight = strtod(argv[i+1], &endptr);
         if( *endptr != '\0' )
         {
            fprintf(stderr, "ERROR: Could not parse argument %s into a number.\n", argv[i+1]);
            return EXIT_FAILURE;
         }

         ++i;
         continue;
      }

      /* see whether the argument can be parsed as a positive real number */
      if( rsize <= ncircles )
      {
         rsize += 5;
         SCIP_ALLOC_ABORT( BMSreallocMemoryArray(&r, rsize) );
      }
      r[ncircles] = strtod(argv[i], &endptr);
      if( *endptr == '\0' && endptr != argv[i] && r[ncircles] > 0.0 )
      {
         ++ncircles;
         continue;
      }

      fprintf(stderr, "ERROR: Unknown option %s.\n", argv[i]);
      return EXIT_FAILURE;
   }

   if( ncircles == 0 )
   {
      assert(rsize >= 5);
      r[0] = 0.25;
      r[1] = 0.25;
      r[2] = 0.4;
      r[3] = 0.7;
      r[4] = 0.1;
      ncircles = 5;
   }

   /* run the actual circle packing example (setting up SCIP, solving the problem, showing the solution) */
   retcode = runPacking(fixwidth, fixheight, dognuplot, domatplotlib);

   /* evaluate return code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
