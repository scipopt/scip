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

/**@file   brachistochrone.c
 * @brief  Minimum time for a particle to go from point A to B under gravity only.
 * @author Anass Meskini
 * 
 * This is an example that uses expressions and expression trees to set up non-linear constraints in SCIP when used as 
 * a callable library. This example implements a discretized model to obtain the trajectory associated with the shortest
 * time to go from point A to B for a particle under gravity only.
 * 
 * The model:
 * 
 * Given \f$N\f$ number of points for the discretisation of the trajectory, we can approximate the time to go from
 * \f$(x_0,y_0)\f$ to \f$ (x_N,y_N)\f$ for a given trajectory by: \f$T = \sqrt{\frac{2}{g}} 
 * \sum_{0}^{N-1} \frac{\sqrt{(y_{i+1} - y_i)^2 + (x_{i+1} - x_i)^2}}{\sqrt{1-y_{i+1}} + \sqrt{1 - y_i}}\f$
 * 
 * The optimisation problem is \f$ \min\limits_{x_0,\dots x_N, y_0,\dots y_N } \sum_{i=0}^{N-1} t_i \f$, such that:
 * \f$t_i = \frac{\sqrt{(y_{i+1} - y_i)^2 + (x_{i+1} - x_i)^2}}{\sqrt{1-y_{i+1}} + \sqrt{1 - y_i}}\f$
 * 
 * A more detailed description of the model can be found in the brachistochrone directory in:
 *
 * http://scip.zib.de/workshop2018/pyscipopt-exercises.tgz
 * 
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <stdlib.h>

#include "scip/pub_misc.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/** default number of points for the discretization of trajectory */
#define DEFAULT_NPOINTS 4

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
   SCIP_Real*            coord               /**< array containing [y(0), y(N), x(0), x(N)] */
   )
{
   /* variables: 
    * t[i] i=0..N-1, such that: value function=sum t[i]
    * y[i] i=0..N, projection of trajectory on y-axis
    * x[i] i=0..N, projection of trajectory on x-axis
    */
   SCIP_VAR** t;
   SCIP_VAR** y;
   SCIP_VAR** x;
 
   char namet[SCIP_MAXSTRLEN];
   char namey[SCIP_MAXSTRLEN];
   char namex[SCIP_MAXSTRLEN];

   SCIP_Real ylb;
   SCIP_Real yub;
   SCIP_Real xlb;
   SCIP_Real xub;

   unsigned int i;

   /* create empty problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "brachistochrone") );

   /* allocate arrays of SCIP_VAR* */
   SCIP_CALL( SCIPallocBufferArray(scip, &t, (size_t) n ) );
   SCIP_CALL( SCIPallocBufferArray(scip, &y, (size_t) n + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &x, (size_t) n + 1) );

   /* create and add variables to the problem and set the initial and final point constraints through upper and lower
    * bounds
    */
   for( i = 0; i < n+1; ++i )
   {
      /* setting up the names of the variables */
      namet[0] = '\0';
      namey[0] = '\0';
      namex[0] = '\0';
      if( i < n )
         SCIPsnprintf(namet, SCIP_MAXSTRLEN, "t(%d)", i);
      SCIPsnprintf(namey, SCIP_MAXSTRLEN, "y(%d)", i);
      SCIPsnprintf(namex, SCIP_MAXSTRLEN, "x(%d)", i);

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
         xlb = coord[2];
         xub = coord[3];
      }

      if( i < n )
      {
         SCIP_CALL( SCIPcreateVarBasic(scip, &t[i], namet, 0.0, SCIPinfinity(scip), 1.0, SCIP_VARTYPE_CONTINUOUS) );
      }
      SCIP_CALL( SCIPcreateVarBasic(scip, &y[i], namey, ylb, yub, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], namex, xlb, xub, 0.0, SCIP_VARTYPE_CONTINUOUS) );

      if( i < n )
      {
         SCIP_CALL( SCIPaddVar(scip, t[i]) );
      }
      SCIP_CALL( SCIPaddVar(scip, y[i]) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );
   }

   /* add constraint related to the linear formulation of the value function 
    * f = sum{ t[i] }
    * t[i] =  sqrt( (y(i+1) - y(i))^2 +  (x(i+1) - x(i))^2 ) / ( sqrt(1 - y(i+1)) + sqrt(1 - y(i))  )
    * in the loop, create the i-th constraint
    */
   {
      SCIP_CONS* conslinfor;

      /* child expressions: 
       * yplusexpr: expression for y[i+1]
       * yplusexprsec: second expression for y[i+1] (it appears twice in the constraint)
       * yexpr: expression for y[i]
       * yexprsec: second expression for y[i]
       * xplusexpr: expression for x[i+1]
       * xexpr: expression for x[i]
       */
      SCIP_EXPR* yplusexpr;
      SCIP_EXPR* yplusexprsec;
      SCIP_EXPR* yexpr;
      SCIP_EXPR* yexprsec;
      SCIP_EXPR* xplusexpr;
      SCIP_EXPR* xexpr;

      /* intermediary expressions */
      SCIP_EXPR* expr;
      SCIP_EXPR* expr1;
      SCIP_EXPR* expr2;
      SCIP_EXPR* expr3;
      SCIP_EXPR* expr4;
      SCIP_EXPR* expr5;
      SCIP_EXPR* expr6;
      SCIP_EXPR* expr7;
      SCIP_EXPR* expr8;
      SCIP_EXPR* expr9;
      SCIP_EXPR* expr10;
      SCIP_EXPR* expr11;

      /* tree to hold the non-linear part of the constraint */ 
      SCIP_EXPRTREE* exprtree;

      SCIP_Real coef1[2] = {1.0, -1.0};
      SCIP_Real coef2 = -1.0;
      SCIP_Real coef3 = 1.0;

      /* at each iteration create an expression for the non-linear part of the i-th constraint and add it the problem */
      for( i = 0; i < n; ++i )
      {
         /* vars to be added to the exprtree in this step of the loop */
         SCIP_VAR* varstoadd[6] = { y[i+1], y[i+1], y[i], y[i], x[i+1], x[i] };

         /* create expressions meant to be child expressions in the tree. give different indexes to the expressions to 
          * assign the correct variables to them later
          */
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &yplusexpr, SCIP_EXPR_VARIDX, 0) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &yplusexprsec, SCIP_EXPR_VARIDX,  1) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &yexpr, SCIP_EXPR_VARIDX, 2) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &yexprsec, SCIP_EXPR_VARIDX, 3) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &xplusexpr, SCIP_EXPR_VARIDX, 4) );
         SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &xexpr, SCIP_EXPR_VARIDX, 5) );

         {
            /* set up the i-th constraint
             * expr1: y[i+1] - y[i]
             * expr2: x[i+1] - x[i]
             * expr6: sqrt( (y[i+1] - y[i])^2 + (x[i+1] - x[i])^2) [1]
             * expr7: 1 - y[i+1]
             * expr8: 1 - y[i]
             * expr11: sqrt(1 - y[i+1]) + sqrt(1 - y[i]) [2]
             * expr: [1] / [2]
             */
            SCIP_EXPR* children1[2] = {yplusexpr, yexpr};
            SCIP_EXPR* children2[2] = {xplusexpr, xexpr};
            SCIP_EXPR* children3 = yplusexprsec;
            SCIP_EXPR* children4 = yexprsec;

            SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &expr1, 2, children1, coef1, 0.0) );
            SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &expr2, 2, children2, coef1, 0.0) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr3, SCIP_EXPR_SQUARE, expr1) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr4, SCIP_EXPR_SQUARE, expr2) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr5, SCIP_EXPR_PLUS, expr3, expr4) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr6, SCIP_EXPR_SQRT, expr5) );
      
            SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &expr7, 1, &children3, &coef2, 1.0) );
            SCIP_CALL( SCIPexprCreateLinear(SCIPblkmem(scip), &expr8, 1, &children4, &coef2, 1.0) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr9, SCIP_EXPR_SQRT, expr7) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr10, SCIP_EXPR_SQRT, expr8) );
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr11, SCIP_EXPR_PLUS, expr9, expr10) );
      
            SCIP_CALL( SCIPexprCreate(SCIPblkmem(scip), &expr, SCIP_EXPR_DIV, expr6, expr11) );
         }
         /* create the expression tree with expr as root */
         SCIP_CALL( SCIPexprtreeCreate(SCIPblkmem(scip), &exprtree, expr, 0, 0, NULL) );
   
         SCIP_CALL( SCIPexprtreeAddVars(exprtree, 6, varstoadd) );
         /* use the tree and a linear term to add  the constraint */
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &conslinfor, "Linear Formulation", 1, &t[i],
                                                 &coef2, 1, &exprtree, &coef3, 0.0, 0.0) );
         /* add the constraint to the problem, release the constraint and free the tree */
         SCIP_CALL( SCIPaddCons(scip, conslinfor) );
         SCIP_CALL( SCIPexprtreeFree(&exprtree) );
         SCIP_CALL( SCIPreleaseCons(scip, &conslinfor) );
      }
   }
   
   /* release problem variables */
   for( i = 0; i < n+1; ++i )
   {
      if( i < n )
      {
         SCIP_CALL( SCIPreleaseVar(scip, &t[i]) );
      }
      SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &x[i]) );
   }

   /* free arrays allocated */
   SCIPfreeBufferArray(scip, &x);
   SCIPfreeBufferArray(scip, &y);
   SCIPfreeBufferArray(scip, &t);

   return SCIP_OKAY;
}

/** runs the brachistochrone example*/
static
SCIP_RETCODE runBrachistochrone(
   unsigned int          n,                  /**< number of points for discretization */
   SCIP_Real*            coord               /**< array containing [y(0), y(N), x(0), x(N)] */
)
{
   SCIP* scip;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "**********************************************\n");
   SCIPinfoMessage(scip, NULL, "* Running Brachistochrone Problem            *\n");
   SCIPinfoMessage(scip, NULL, "* between A=(%g,%g) and B=(%g,%g) with %d points *\n", coord[2], coord[0], coord[3], coord[1], n);
   SCIPinfoMessage(scip, NULL, "**********************************************\n");
   SCIPinfoMessage(scip, NULL, "\n");

   /* set gap at which SCIP will stop */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", 0.05) );

   /* Set larger constraint violation tolerence to speed up the solving process */
   SCIP_CALL( SCIPsetRealParam(scip, "numerics/feastol", 1e-4) );

   /* if no NLP solver, then stop at first solution */
   if( SCIPgetNNlpis(scip) == 0 )
   {
      SCIPwarningMessage(scip, "No NLP solver available. There is little hope to solve this problem. Stopping after first solution.\n");
      SCIP_CALL( SCIPsetIntParam(scip, "limits/solutions", 1) );
   }

   SCIP_CALL( setupProblem(scip, n, coord) );

   SCIPinfoMessage(scip, NULL, "Original problem:\n");
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPpresolve(scip) );

   SCIPinfoMessage(scip, NULL, "\nSolving...\n");
   SCIP_CALL( SCIPsolve(scip) );

   SCIP_CALL( SCIPfreeTransform(scip) );

   if( SCIPgetNSols(scip) > 0 )
   {
      SCIPinfoMessage(scip, NULL, "\nSolution:\n");
      SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
   }

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
   unsigned int n = DEFAULT_NPOINTS;
   SCIP_Real coord[4] = { Y_START, Y_END, X_START, X_END };

   /* change some parameters if given as arguments */
   if( argc == 4 || argc == 2  )
   {
       char *end1 = NULL;
       char *end2 = NULL;
       char *end3 = NULL;

       n = strtol(argv[1], &end1, 10);
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
