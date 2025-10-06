/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   brachistochrone.c
 * @brief  Computing a minimum-time trajectory for a particle to move from point A to B under gravity only
 * @author Anass Meskini
 * @author Stefan Vigerske
 *
 * This is an example that uses expressions to setup non-linear constraints in SCIP when used as
 * a callable library. This example implements a discretized model to obtain the trajectory associated with the shortest
 * time to go from point A to B for a particle under gravity only.
 *
 * The model:
 *
 * Given \f$N\f$ number of points for the discretisation of the trajectory, we can approximate the time to go from
 * \f$(x_0,y_0)\f$ to \f$ (x_N,y_N)\f$ for a given trajectory by \f[T = \sqrt{\frac{2}{g}}
 * \sum_{0}^{N-1} \frac{\sqrt{(y_{i+1} - y_i)^2 + (x_{i+1} - x_i)^2}}{\sqrt{1-y_{i+1}} + \sqrt{1 - y_i}}.\f]
 * We seek to minimize \f$T\f$.
 * A more detailed description of the model can be found in the brachistochrone directory of http://scipopt.org/workshop2018/pyscipopt-exercises.tgz
 *
 * Passing this equation as it is to SCIP does not lead to satisfying results, though, so we reformulate a bit.
 * Let \f$t_i \geq \frac{\sqrt{(y_{i+1} - y_i)^2 + (x_{i+1} - x_i)^2}}{\sqrt{1-y_{i+1}} + \sqrt{1 - y_i}}\f$.
 * To avoid a potential division by zero, we rewrite this as
 * \f$t_i (\sqrt{1-y_{i+1}} + \sqrt{1 - y_i}) \geq \sqrt{(y_{i+1} - y_i)^2 + (x_{i+1} - x_i)^2}, t_i\geq 0\f$.
 * Further, introduce \f$v_i \geq 0\f$ such that \f$v_i \geq \sqrt{(y_{i+1} - y_i)^2 + (x_{i+1} - x_i)^2}\f$.
 * Then we can state the optimization problem as
 * \f{align}{ \min \;& \sqrt{\frac{2}{g}} \sum_{i=0}^{N-1} t_i \\
 *   \text{s.t.} \; & t_i (\sqrt{1-y_{i+1}} + \sqrt{1 - y_i}) \geq v_i \\
 *     & v_i^2 \geq (y_{i+1} - y_i)^2 + (x_{i+1} - x_i)^2 \\
 *     & t_i \geq 0,\; v_i \geq 0 \\
 *     & i = 0, \ldots, N-1
 * \f}
 *
 * Further, we can require that the particle moves only in direction horizontally, that is
 * \f$x_i \leq x_{i+1}\f$ if \f$x_0 \leq x_N\f$, and \f$x_{i+1} \leq x_{i}\f$, otherwise,
 * and that it will not move higher than the start-coordinate.
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>

#include "scip/pub_misc.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/* default start and end points */
#define Y_START  1.0
#define Y_END    0.0
#define X_START  0.0
#define X_END   10.0

/** sets up problem */
static
SCIP_RETCODE setupProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   unsigned int          n,                  /**< number of points for discretization */
   SCIP_Real*            coord,              /**< array containing [y(0), y(N), x(0), x(N)] */
   SCIP_VAR***           xvars,              /**< buffer to store pointer to x variables array */
   SCIP_VAR***           yvars               /**< buffer to store pointer to y variables array */
   )
{
   /* variables:
    * t[i] i=0..N-1, such that: value function=sum t[i]
    * v[i] i=0..N-1, such that: v_i = ||(x_{i+1},y_{i+1})-(x_i,y_i)||_2
    * y[i] i=0..N, projection of trajectory on y-axis
    * x[i] i=0..N, projection of trajectory on x-axis
    */
   SCIP_VAR** t;
   SCIP_VAR** v;
   SCIP_VAR** x;
   SCIP_VAR** y;

   char namet[SCIP_MAXSTRLEN];
   char namev[SCIP_MAXSTRLEN];
   char namey[SCIP_MAXSTRLEN];
   char namex[SCIP_MAXSTRLEN];

   SCIP_Real ylb;
   SCIP_Real yub;
   SCIP_Real xlb;
   SCIP_Real xub;

   /* an upper bound for v */
   SCIP_Real maxdistance = 10.0 * sqrt(SQR(coord[1]-coord[0]) + SQR(coord[3]-coord[2]));
   unsigned int i;

   /* create empty problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "brachistochrone") );

   /* allocate arrays of SCIP_VAR* */
   SCIP_CALL( SCIPallocBufferArray(scip, &t, (size_t) n ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &v, (size_t) n ) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &y, (size_t) n + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &x, (size_t) n + 1) );
   *xvars = x;
   *yvars = y;

   /* create and add variables to the problem and set the initial and final point constraints through upper and lower
    * bounds
    */
   for( i = 0; i < n+1; ++i )
   {
      /* setting up the names of the variables */
      if( i < n )
      {
         (void)SCIPsnprintf(namet, SCIP_MAXSTRLEN, "t(%d)", i);
         (void)SCIPsnprintf(namev, SCIP_MAXSTRLEN, "v(%d)", i);
      }
      (void)SCIPsnprintf(namey, SCIP_MAXSTRLEN, "y(%d)", i);
      (void)SCIPsnprintf(namex, SCIP_MAXSTRLEN, "x(%d)", i);

      /* fixing y(0), y(N), x(0), x(N) through lower and upper bounds */
      if( i == 0 )
      {
         ylb = coord[0];
         yub = coord[0];
         xlb = coord[2];
         xub = coord[2];
      }
      else if( i == n )
      {
         ylb = coord[1];
         yub = coord[1];
         xlb = coord[3];
         xub = coord[3];
      }
      else
      {
         /* constraint the other variables to speed up solving */
         ylb = -SCIPinfinity(scip);
         yub = coord[0];
         xlb = MIN(coord[2], coord[3]);
         xub = MAX(coord[2], coord[3]);
      }

      if( i < n )
      {
         SCIP_CALL( SCIPcreateVarBasic(scip, &t[i], namet, 0.0, SCIPinfinity(scip), sqrt(2.0/9.80665), SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPcreateVarBasic(scip, &v[i], namev, 0.0, maxdistance, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      }
      SCIP_CALL( SCIPcreateVarBasic(scip, &y[i], namey, ylb, yub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], namex, xlb, xub, 0.0, SCIP_VARTYPE_CONTINUOUS) );

      if( i < n )
      {
         SCIP_CALL( SCIPaddVar(scip, t[i]) );
         SCIP_CALL( SCIPaddVar(scip, v[i]) );
      }
      SCIP_CALL( SCIPaddVar(scip, y[i]) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );
   }

   /* add constraints
    *
    * t(i) * sqrt(1 - y(i+1)) + t(i) * sqrt(1 - y(i)) >= v(i)
    * v(i)^2 >= (y(i+1) - y(i))^2 + (x(i+1) - x(i))^2
    * in the loop, create the i-th constraint
    */
   {
      SCIP_CONS* cons;
      char consname[SCIP_MAXSTRLEN];

      /* variable expressions:
       * yplusexpr: expression for y[i+1]
       * yexpr: expression for y[i]
       * texpr: expression for t[i]
       * vexpr: expression for v[i]
       */
      SCIP_EXPR* yplusexpr;
      SCIP_EXPR* yexpr;
      SCIP_EXPR* texpr;
      SCIP_EXPR* vexpr;

      /* intermediary expressions */
      SCIP_EXPR* expr1;
      SCIP_EXPR* expr2;
      SCIP_EXPR* expr3;
      SCIP_EXPR* expr4;
      SCIP_EXPR* expr5;
      SCIP_EXPR* expr6;
      SCIP_EXPR* expr7;

      SCIP_Real minusone = -1.0;

      for( i = 0; i < n; ++i )
      {
         SCIP_EXPR* exprs[3];
         SCIP_Real coefs[3];
         SCIP_VAR* quadvars1[7];
         SCIP_VAR* quadvars2[7];
         SCIP_Real quadcoefs[7];

         /* create expressions for variables that are used in nonlinear constraint */
         SCIP_CALL( SCIPcreateExprVar(scip, &yplusexpr, y[i+1], NULL, NULL) );
         SCIP_CALL( SCIPcreateExprVar(scip, &yexpr, y[i], NULL, NULL) );
         SCIP_CALL( SCIPcreateExprVar(scip, &texpr, t[i], NULL, NULL) );
         SCIP_CALL( SCIPcreateExprVar(scip, &vexpr, v[i], NULL, NULL) );

         /* set up the i-th constraint
          * expr1: 1 - y[i+1]
          * expr2: 1 - y[i]
          * expr3: sqrt(1 - y[i+1])
          * expr4: sqrt(1 - y[i])
          * expr5: t[i] * sqrt(1 - y[i+1])
          * expr6: t[i] * sqrt(1 - y[i])
          * expr7: t[i] * sqrt(1 - y[i+1]) + t[i] * sqrt(1 - y[i]) - v[i]
          */
         SCIP_CALL( SCIPcreateExprSum(scip, &expr1, 1, &yplusexpr, &minusone, 1.0, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprSum(scip, &expr2, 1, &yexpr, &minusone, 1.0, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &expr3, expr1, 0.5, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &expr4, expr2, 0.5, NULL, NULL) );
         exprs[0] = expr3;
         exprs[1] = texpr;
         SCIP_CALL( SCIPcreateExprProduct(scip, &expr5, 2, exprs, 1.0, NULL, NULL) );
         exprs[0] = expr4;
         SCIP_CALL( SCIPcreateExprProduct(scip, &expr6, 2, exprs, 1.0, NULL, NULL) );
         exprs[0] = expr5;  coefs[0] =  1.0;
         exprs[1] = expr6;  coefs[1] =  1.0;
         exprs[2] = vexpr;  coefs[2] = -1.0;
         SCIP_CALL( SCIPcreateExprSum(scip, &expr7, 3, exprs, coefs, 0.0, NULL, NULL) );

         /* create the constraint expr7 >= 0, add to the problem, and release it */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "timestep(%d)", i);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, consname, expr7, 0.0, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

         /* release exprs */
         SCIP_CALL( SCIPreleaseExpr(scip, &expr7) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr6) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr5) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr4) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr3) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr2) );
         SCIP_CALL( SCIPreleaseExpr(scip, &expr1) );

         SCIP_CALL( SCIPreleaseExpr(scip, &vexpr) );
         SCIP_CALL( SCIPreleaseExpr(scip, &texpr) );
         SCIP_CALL( SCIPreleaseExpr(scip, &yexpr) );
         SCIP_CALL( SCIPreleaseExpr(scip, &yplusexpr) );

         /* create constraint v_i^2 >= (y_{i+1}^2 - 2*y_{i+1}y_i + y_i^2) + (x_{i+1}^2 - 2*x_{i+1}x_i + x_i^2)
          * SCIP should recognize that this can be formulated as SOC
          */
         quadvars1[0] = y[i];   quadvars2[0] = y[i];   quadcoefs[0] =  1.0;
         quadvars1[1] = y[i+1]; quadvars2[1] = y[i+1]; quadcoefs[1] =  1.0;
         quadvars1[2] = x[i];   quadvars2[2] = x[i];   quadcoefs[2] =  1.0;
         quadvars1[3] = x[i+1]; quadvars2[3] = x[i+1]; quadcoefs[3] =  1.0;
         quadvars1[4] = y[i];   quadvars2[4] = y[i+1]; quadcoefs[4] = -2.0;
         quadvars1[5] = x[i];   quadvars2[5] = x[i+1]; quadcoefs[5] = -2.0;
         quadvars1[6] = v[i];   quadvars2[6] = v[i];   quadcoefs[6] = -1.0;

         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "steplength(%d)", i);
         SCIP_CALL( SCIPcreateConsQuadraticNonlinear(scip, &cons, consname, 0, NULL, NULL, 7, quadvars1, quadvars2, quadcoefs, -SCIPinfinity(scip), 0.0, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE) );

         /* add the constraint to the problem and forget it */
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

         /* add constraint x[i] <= x[i+1], if x[0] < x[N], otherwise add x[i+1] <= x[i] */
         (void)SCIPsnprintf(consname, SCIP_MAXSTRLEN, "xorder(%d)", i);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, consname, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i],   coord[2] < coord[3] ?  1.0 : -1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, cons, x[i+1], coord[2] < coord[3] ? -1.0 :  1.0) );

         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );
      }
   }

   /* release intermediate variables */
   for( i = 0; i < n; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &t[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &v[i]) );
   }

   /* free arrays allocated */
   SCIPfreeBufferArray(scip, &t);
   SCIPfreeBufferArray(scip, &v);

   return SCIP_OKAY;
}

/** plots solution by use of gnuplot */
static
void visualizeSolutionGnuplot(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< solution to plot */
   unsigned int          n,                  /**< number of points for discretization */
   SCIP_VAR**            x,                  /**< x coordinates */
   SCIP_VAR**            y                   /**< y coordinates */
   )
{
#if _POSIX_C_SOURCE < 2
   SCIPinfoMessage(scip, NULL, "No POSIX version 2. Try http://distrowatch.com/.");
#else
   FILE* stream;
   unsigned int i;

   /* -p (persist) to keep the plot open after gnuplot terminates (if terminal is not dumb) */
   stream = popen("gnuplot -p", "w");
   if( stream == NULL )
   {
      SCIPerrorMessage("Could not open pipe to gnuplot.\n");
      return;
   }
   /* take out this line to get a non-ascii plot */
   fputs("set terminal dumb\n", stream);

   fprintf(stream, "plot '-' smooth csplines title \"Time = %.4fs\"\n", SCIPgetSolOrigObj(scip, sol));
   for( i = 0; i < n+1; ++i )
      fprintf(stream, "%g %g\n", SCIPgetSolVal(scip, sol, x[i]), SCIPgetSolVal(scip, sol, y[i]));
   fputs("e\n", stream);

   (void)pclose(stream);
#endif
}

/** runs the brachistochrone example*/
static
SCIP_RETCODE runBrachistochrone(
   unsigned int          n,                  /**< number of points for discretization */
   SCIP_Real*            coord               /**< array containing [y(0), y(N), x(0), x(N)] */
   )
{
   SCIP* scip;
   SCIP_VAR** y;
   SCIP_VAR** x;
   unsigned int i;

   assert(n >= 2);

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "**********************************************\n");
   SCIPinfoMessage(scip, NULL, "* Running Brachistochrone Problem            *\n");
   SCIPinfoMessage(scip, NULL, "* between A=(%g,%g) and B=(%g,%g) with %u points *\n", coord[2], coord[0], coord[3], coord[1], n);
   SCIPinfoMessage(scip, NULL, "**********************************************\n");
   SCIPinfoMessage(scip, NULL, "\n");

   /* set gap at which SCIP will stop */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", 0.05) );

   SCIP_CALL( setupProblem(scip, n, coord, &x, &y) );

   SCIPinfoMessage(scip, NULL, "Original problem:\n");
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

   SCIPinfoMessage(scip, NULL, "\nSolving...\n");
   SCIP_CALL( SCIPsolve(scip) );

   if( SCIPgetNSols(scip) > 0 )
   {
      SCIPinfoMessage(scip, NULL, "\nSolution:\n");
      SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );

      visualizeSolutionGnuplot(scip, SCIPgetBestSol(scip), n, x, y);
   }

   /* release problem variables */
   for( i = 0; i < n+1; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &x[i]) );
   }

   /* free arrays allocated */
   SCIPfreeMemoryArray(scip, &x);
   SCIPfreeMemoryArray(scip, &y);

   SCIP_CALL( SCIPfree(&scip) );

   return SCIP_OKAY;
}

/** main method starting SCIP */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< arguments: number of points and end coordinates y(N), x(N)*/
   )
{
   SCIP_RETCODE retcode;

   /* setting up default problem parameters */
   unsigned int n = 3;
   SCIP_Real coord[4] = { Y_START, Y_END, X_START, X_END };

   /* change some parameters if given as arguments */
   if( argc == 4 || argc == 2  )
   {
       char *end1 = NULL;
       char *end2 = NULL;
       char *end3 = NULL;

       n = (unsigned int) strtol(argv[1], &end1, 10);
       if( argc == 4 )
       {
          coord[1] = strtof(argv[2], &end2);
          coord[3] = strtof(argv[3], &end3);
       }

       if( end1 == argv[1] || ( argc == 4 && ( end2 == argv[2] || end3 == argv[3] ) ) )
       {
          fprintf(stderr, "Error: expected real values as arguments.\n");
          return EXIT_FAILURE;
       }
   }
   else if( argc != 1 )
   {
       fprintf(stderr, "Usage: %s [<N> [<y(N)> <x(N)>]]\n", argv[0]);
       return EXIT_FAILURE;
   }

   /* check that y(0) > y(N) */
   if( coord[0] <= coord[1] )
   {
      fprintf(stderr, "Error: expected y(N) < 1.0\n");
      return EXIT_FAILURE;
   }

   retcode = runBrachistochrone(n, coord);

   /* evaluate return code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
