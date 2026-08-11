/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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

/**@file   lj.c
 * @brief  Lennard-Jones Cluster problem with custom nonlinear handler.
 * @author Stefan Vigerske
 *
 * https://www-wales.ch.cam.ac.uk/~wales/CCD/jon/structures/LJ.html
 * https://arxiv.org/abs/cond-mat/9803344
 * https://github.com/anatoliy-kuznetsov/EuclidLib/blob/main/gams/lj_10.gms
 */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "nlhdlr_lj.h"

#define BOX 9.0

/** create problem in given SCIP and add all variables and constraints to model the Lennard-Jones Cluster problem */
static SCIP_RETCODE setupProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   nparticles,         /**< number of particles */
   SCIP_Bool             distvars            /**< whether to introduce extra explicit variables and constraints for
                                              *   squared distance of particles */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR** x;
   SCIP_VAR** y;
   SCIP_VAR** z;
   SCIP_VAR** r = NULL;
   SCIP_VAR** p;
   int i, j;

   /* create empty problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "LennardJones") );

   /* create variables for position of particles */
   SCIP_CALL( SCIPallocMemoryArray(scip, &x, nparticles) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &y, nparticles) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &z, nparticles) );
   for( i = 0; i < nparticles; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "x_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &x[i], name, -BOX, BOX, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, x[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "y_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &y[i], name, -BOX, BOX, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, y[i]) );

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "z_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &z[i], name, -BOX, BOX, 0.0, SCIP_VARTYPE_CONTINUOUS) );
      SCIP_CALL( SCIPaddVar(scip, z[i]) );
   }

   /* create variables for squared distance between particles and their Lennard-Jones Potential (objcoef 4.0) */
   if( distvars )
   {
      SCIP_CALL( SCIPallocClearMemoryArray(scip, &r, nparticles * nparticles) );
   }
   SCIP_CALL( SCIPallocClearMemoryArray(scip, &p, nparticles * nparticles) );
   for( i = 0; i < nparticles; ++i )
   {
      for( j = i + 1; j < nparticles; ++j )
      {
         if( distvars )
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "r_%d_%d", i, j);
            SCIP_CALL( SCIPcreateVarBasic(scip, &r[i * nparticles + j], name, 0.0, SCIPinfinity(scip), 0.0,
               SCIP_VARTYPE_CONTINUOUS) );
            SCIP_CALL( SCIPaddVar(scip, r[i * nparticles + j]) );
         }

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "p_%d_%d", i, j);
         SCIP_CALL( SCIPcreateVarBasic(scip, &p[i * nparticles + j], name, -SCIPinfinity(scip), SCIPinfinity(scip), 4.0,
            SCIP_VARTYPE_CONTINUOUS) );
         SCIP_CALL( SCIPaddVar(scip, p[i * nparticles + j]) );
      }
   }

   /* create constraints that couple particle position (x,y,z) to squared particle distance (r)
    * and couple particle distance to LJ potential
    */
   for( i = 0; i < nparticles; ++i )
   {
      for( j = i + 1; j < nparticles; ++j )
      {
         SCIP_CONS* cons;
         SCIP_EXPR* distexpr;
         SCIP_EXPR* pexpr;
         SCIP_EXPR* rpow1;
         SCIP_EXPR* rpow2;
         SCIP_EXPR* sum;
         SCIP_EXPR* terms[3];
         SCIP_Real coefs[3];

         /* (x[i] - x[j])^2 + (y[i] - y[j])^2 + (z[i] - z[j])^2 - r[i,j] = 0 */
         SCIP_VAR* quadvars1[9];
         SCIP_VAR* quadvars2[9];
         SCIP_Real quadcoefs[9];

         /* x[i]^2 */
         quadvars1[0] = x[i];
         quadvars2[0] = x[i];
         quadcoefs[0] = 1.0;

         /* -2 x[i] x[j] */
         quadvars1[1] = x[i];
         quadvars2[1] = x[j];
         quadcoefs[1] = -2.0;

         /* x[i]^2 */
         quadvars1[2] = x[j];
         quadvars2[2] = x[j];
         quadcoefs[2] = 1.0;

         /* y[i]^2 */
         quadvars1[3] = y[i];
         quadvars2[3] = y[i];
         quadcoefs[3] = 1.0;

         /* -2 y[i] y[j] */
         quadvars1[4] = y[i];
         quadvars2[4] = y[j];
         quadcoefs[4] = -2.0;

         /* y[i]^2 */
         quadvars1[5] = y[j];
         quadvars2[5] = y[j];
         quadcoefs[5] = 1.0;

         /* z[i]^2 */
         quadvars1[6] = z[i];
         quadvars2[6] = z[i];
         quadcoefs[6] = 1.0;

         /* -2 z[i] z[j] */
         quadvars1[7] = z[i];
         quadvars2[7] = z[j];
         quadcoefs[7] = -2.0;

         /* z[i]^2 */
         quadvars1[8] = z[j];
         quadvars2[8] = z[j];
         quadcoefs[8] = 1.0;

         if( distvars )
         {
            SCIP_Real minusone = -1.0;

            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "dist_%d_%d", i, j);
            SCIP_CALL( SCIPcreateConsBasicQuadraticNonlinear(scip, &cons, name, 1, &r[i * nparticles + j], &minusone, 9,
               quadvars1, quadvars2, quadcoefs, 0.0, 0.0) );
            SCIP_CALL( SCIPaddCons(scip, cons) );
            SCIP_CALL( SCIPreleaseCons(scip, &cons) );

            SCIP_CALL( SCIPcreateExprVar(scip, &distexpr, r[i * nparticles + j], NULL, NULL) );
         }
         else
         {
            SCIP_CALL( SCIPcreateExprQuadratic(scip, &distexpr, 0, NULL, NULL, 9, quadvars1, quadvars2, quadcoefs, NULL,
               NULL) );
         }

         /* p[i,j] >= r[i,j]^-6 - r[i,j]^-3 */
         SCIP_CALL( SCIPcreateExprVar(scip, &pexpr, p[i * nparticles + j], NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &rpow1, distexpr, -6.0, NULL, NULL) );
         SCIP_CALL( SCIPcreateExprPow(scip, &rpow2, distexpr, -3.0, NULL, NULL) );
         terms[0] = pexpr;
         coefs[0] = -1.0;
         terms[1] = rpow1;
         coefs[1] = 1.0;
         terms[2] = rpow2;
         coefs[2] = -1.0;
         SCIP_CALL( SCIPcreateExprSum(scip, &sum, 3, terms, coefs, 0.0, NULL, NULL) );

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pot_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicNonlinear(scip, &cons, name, sum, -SCIPinfinity(scip), 0.0) );
         SCIP_CALL( SCIPaddCons(scip, cons) );
         SCIP_CALL( SCIPreleaseCons(scip, &cons) );

         SCIP_CALL( SCIPreleaseExpr(scip, &pexpr) );
         SCIP_CALL( SCIPreleaseExpr(scip, &distexpr) );
         SCIP_CALL( SCIPreleaseExpr(scip, &rpow1) );
         SCIP_CALL( SCIPreleaseExpr(scip, &rpow2) );
         SCIP_CALL( SCIPreleaseExpr(scip, &sum) );
      }
   }

   /* release variables */
   for( i = 0; i < nparticles; ++i )
   {
      for( j = i + 1; j < nparticles; ++j )
      {
         if( distvars )
         {
            SCIP_CALL( SCIPreleaseVar(scip, &r[i * nparticles + j]) );
         }
         SCIP_CALL( SCIPreleaseVar(scip, &p[i * nparticles + j]) );
      }
      SCIP_CALL( SCIPreleaseVar(scip, &x[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &y[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &z[i]) );
   }

   SCIPfreeMemoryArray(scip, &p);
   SCIPfreeMemoryArrayNull(scip, &r);
   SCIPfreeMemoryArray(scip, &z);
   SCIPfreeMemoryArray(scip, &y);
   SCIPfreeMemoryArray(scip, &x);

   return SCIP_OKAY;
}

/** run LJ Cluster example
 *
 *  Sets up SCIP and the SCIP problem, solves the problem, and shows the solution.
 */
static SCIP_RETCODE runLJ(
   int                   nparticles,         /**< number of particles */
   SCIP_Bool             distvars,           /**< whether to introduce extra explicit variables and constraints for
                                              *   squared distance of particles */
   const char*           setfile             /**< name of settings file to attempt reading */
   )
{
   SCIP* scip;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIP_CALL( SCIPincludeNlhdlrLJ(scip) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "**************************************************\n");
   SCIPinfoMessage(scip, NULL, "* Running Lennard-Jones Cluster for %d particles *\n", nparticles);
   SCIPinfoMessage(scip, NULL, "**************************************************\n");
   SCIPinfoMessage(scip, NULL, "\n");

   SCIPprintVersion(scip, NULL);
   SCIPprintExternalCodes(scip, NULL);

   SCIP_CALL( setupProblem(scip, nparticles, distvars) );

   if( nparticles <= 5 )
   {
      SCIPinfoMessage(scip, NULL, "\nOriginal problem:\n");
      SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );
   }

   /* By default, SCIP tries to close the gap between primal and dual bound completely.
    * This can take very long for this example, so we increase the gap tolerance to have
    * SCIP stop when the gap between primal and dual bound is already at most 0.01%.
    */
   SCIP_CALL( SCIPsetRealParam(scip, "limits/gap", 0.0001) );

   if( SCIPfileExists(setfile) )
   {
      SCIPinfoMessage(scip, NULL, "\nReading %s\n", setfile);
      SCIP_CALL( SCIPreadParams(scip, setfile) );
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "\nSettings file %s not found.\n", setfile);
   }

   SCIPinfoMessage(scip, NULL, "\nSolving...\n");
   SCIP_CALL( SCIPsolve(scip) );

   if( SCIPgetNSols(scip) > 0 )
   {
      SCIPinfoMessage(scip, NULL, "\nSolution:\n");
      SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
   }

   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   SCIP_CALL( SCIPfree(&scip) );

   return SCIP_OKAY;
}

/** main method starting SCIP */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{  /*lint --e{715}*/
   SCIP_RETCODE retcode;
   long nparticles = 0;
   char* endptr;

   if( argc <= 1 )
   {
      printf("usage: %s NPARTICLES <setfile> <\n", argv[0]);
      return EXIT_SUCCESS;
   }

   nparticles = strtol(argv[1], &endptr, 10);
   if( *endptr != '\0' )
   {
      fprintf(stderr, "ERROR: Could not parse argument %s into a number.\n", argv[1]);
      return EXIT_FAILURE;
   }

   if( nparticles < 2 )
   {
      fputs("ERROR: Number of particles needs to be >= 2.\n", stderr);
      return EXIT_FAILURE;
   }

   if( nparticles >= 100 )
   {
      fputs("ERROR: Number of particles needs to be < 100.\n", stderr);
      return EXIT_FAILURE;
   }

   /* run the example (setting up SCIP, solving the problem, showing the solution)
    * second argument specified whether to introduce explicit variables for particle distances (r_ij)
    */
   retcode = runLJ((int)nparticles, FALSE, argc >= 3 ? argv[2] : "scip.set");

   /* evaluate return code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
