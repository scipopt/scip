/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   gastrans.c
 * @brief  Simple Gas Transportation Model
 * @author Stefan Vigerske
 *
 * This example shows how to setup abspower constraints in SCIP when using SCIP as callable library.
 * The example implements a model for the distribution of gas through a network of pipelines, which
 * is formulated as a cost minimization subject to nonlinear flow-pressure relations, material balances,
 * and pressure bounds. The Belgian gas network is used as an example.
 *
 * The model is taken from the GAMS model library:
 * http://www.gams.com/modlib/libhtml/gastrans.htm
 *
 * Original model source:
 * @par
 *   D. de Wolf and Y. Smeers@n
 *   The Gas Transmission Problem Solved by and Extension of the Simplex Algorithm@n
 *   Management Science 46, 11 (2000), 1454-1465
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <math.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/* Model parameters */

/** node data structure */
typedef struct NodeData {
   const char*           name;
   SCIP_Real             supplylower;
   SCIP_Real             supplyupper;
   SCIP_Real             pressurelower;
   SCIP_Real             pressureupper;
   SCIP_Real             cost;
} NodeData;

/** arc data structure */
typedef struct ArcData {
   int                   node1;
   int                   node2;
   SCIP_Real             diameter;
   SCIP_Real             length;
   SCIP_Bool             active;
} ArcData;

/** number of nodes (towns) */
#define nnodes 20

/** number of arcs */
#define narcs  24

/** value we use to represent infinity */
#define infinity 1e+20

/** data of nodes */
static const NodeData nodedata[] =
{
   /* name         supplylo   supplyup pressurelo pressureup   cost */
   {"Anderlues",        0.0,       1.2,       0.0,      66.2,   0.0  },  /*  0 */
   {"Antwerpen",  -infinity,    -4.034,      30.0,      80.0,   0.0  },  /*  1 */
   {"Arlon",      -infinity,    -0.222,       0.0,      66.2,   0.0  },  /*  2 */
   {"Berneau",          0.0,       0.0,       0.0,      66.2,   0.0  },  /*  3 */
   {"Blaregnies", -infinity,   -15.616,      50.0,      66.2,   0.0  },  /*  4 */
   {"Brugge",     -infinity,    -3.918,      30.0,      80.0,   0.0  },  /*  5 */
   {"Dudzele",          0.0,       8.4,       0.0,      77.0,   2.28 },  /*  6 */
   {"Gent",       -infinity,    -5.256,      30.0,      80.0,   0.0  },  /*  7 */
   {"Liege",      -infinity,    -6.385,      30.0,      66.2,   0.0  },  /*  8 */
   {"Loenhout",         0.0,       4.8,       0.0,      77.0,   2.28 },  /*  9 */
   {"Mons",       -infinity,    -6.848,       0.0,      66.2,   0.0  },  /* 10 */
   {"Namur",      -infinity,    -2.120,       0.0,      66.2,   0.0  },  /* 11 */
   {"Petange",    -infinity,    -1.919,      25.0,      66.2,   0.0  },  /* 12 */
   {"Peronnes",         0.0,      0.96,       0.0,      66.2,   1.68 },  /* 13 */
   {"Sinsin",           0.0,       0.0,       0.0,      63.0,   0.0  },  /* 14 */
   {"Voeren",        20.344,    22.012,      50.0,      66.2,   1.68 },  /* 15 */
   {"Wanze",            0.0,       0.0,       0.0,      66.2,   0.0  },  /* 16 */
   {"Warnand",          0.0,       0.0,       0.0,      66.2,   0.0  },  /* 17 */
   {"Zeebrugge",       8.87,    11.594,       0.0,      77.0,   2.28 },  /* 18 */
   {"Zomergem",         0.0,       0.0,       0.0,      80.0,   0.0  }   /* 19 */
};

/** data of arcs */
static const ArcData arcdata[] =
{
  /* node1  node2  diameter length active */
   { 18,  6, 890.0,   4.0, FALSE },
   { 18,  6, 890.0,   4.0, FALSE },
   {  6,  5, 890.0,   6.0, FALSE },
   {  6,  5, 890.0,   6.0, FALSE },
   {  5, 19, 890.0,  26.0, FALSE },
   {  9,  1, 590.1,  43.0, FALSE },
   {  1,  7, 590.1,  29.0, FALSE },
   {  7, 19, 590.1,  19.0, FALSE },
   { 19, 13, 890.0,  55.0, FALSE },
   { 15,  3, 890.0,   5.0,  TRUE },
   { 15,  3, 395.0,   5.0,  TRUE },
   {  3,  8, 890.0,  20.0, FALSE },
   {  3,  8, 395.0,  20.0, FALSE },
   {  8, 17, 890.0,  25.0, FALSE },
   {  8, 17, 395.0,  25.0, FALSE },
   { 17, 11, 890.0,  42.0, FALSE },
   { 11,  0, 890.0,  40.0, FALSE },
   {  0, 13, 890.0,   5.0, FALSE },
   { 13, 10, 890.0,  10.0, FALSE },
   { 10,  4, 890.0,  25.0, FALSE },
   { 17, 16, 395.5,  10.5, FALSE },
   { 16, 14, 315.5,  26.0,  TRUE },
   { 14,  2, 315.5,  98.0, FALSE },
   {  2, 12, 315.5,   6.0, FALSE }
};

/** gas temperatur (K) */
static const SCIP_Real gastemp = 281.15;

/** absolute rugosity (mm) */
static const SCIP_Real rugosity = 0.05;

/** density of gas relative to air */
static const SCIP_Real density = 0.616;

/** compressibility factor */
static const SCIP_Real compressibility = 0.8;



/** sets up problem */
static
SCIP_RETCODE setupProblem(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_VAR* flow[narcs];
   SCIP_VAR* supply[nnodes];
   SCIP_VAR* pressure[nnodes];
   SCIP_VAR* pressurediff[narcs];

   SCIP_CONS* flowbalance[nnodes];
   SCIP_CONS* pressurediffcons[narcs];
   SCIP_CONS* pressureloss[narcs];

   char name[SCIP_MAXSTRLEN];
   int i;
   int j;

   /* create empty problem */
   SCIP_CALL( SCIPcreateProbBasic(scip, "gastrans") );

   /* create variables for flows and add to problem */
   for( i = 0; i < narcs; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "flow_%s_%s", nodedata[arcdata[i].node1].name, nodedata[arcdata[i].node2].name);
      SCIP_CALL( SCIPcreateVarBasic(scip, &flow[i], name, arcdata[i].active ? 0.0 : -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(scip, flow[i]) );
   }

   /* create variables for pressure difference and add to problem */
   for( i = 0; i < narcs; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pressurediff_%s_%s", nodedata[arcdata[i].node1].name, nodedata[arcdata[i].node2].name);
      SCIP_CALL( SCIPcreateVarBasic(scip, &pressurediff[i], name, -SCIPinfinity(scip), SCIPinfinity(scip), 0.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(scip, pressurediff[i]) );
   }

   /* create variables for supply at nodes and add to problem */
   for( i = 0; i < nnodes; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "supply_%s", nodedata[i].name);
      SCIP_CALL( SCIPcreateVarBasic(scip, &supply[i], name,
         nodedata[i].supplylower == -infinity ? -SCIPinfinity(scip) : nodedata[i].supplylower,
         nodedata[i].supplyupper ==  infinity ?  SCIPinfinity(scip) : nodedata[i].supplyupper,
         nodedata[i].cost, SCIP_VARTYPE_CONTINUOUS) ); /*lint !e777*/

      SCIP_CALL( SCIPaddVar(scip, supply[i]) );
   }

   /* create variables for squared pressure at nodes and add to problem */
   for( i = 0; i < nnodes; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pressure_%s", nodedata[i].name);
      SCIP_CALL( SCIPcreateVarBasic(scip, &pressure[i], name,
         nodedata[i].pressurelower*nodedata[i].pressurelower, nodedata[i].pressureupper*nodedata[i].pressureupper,
         0.0, SCIP_VARTYPE_CONTINUOUS) );

      SCIP_CALL( SCIPaddVar(scip, pressure[i]) );
   }

   /* create flow balance constraints and add to problem
    * for each node i: outflows - inflows = supply
    */
   for( i = 0; i < nnodes; ++i )
   {
      SCIP_Real minusone;

      minusone = -1.0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "flowbalance%s", nodedata[i].name);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &flowbalance[i], name, 1, &supply[i], &minusone, 0.0, 0.0) );

      for( j = 0; j < narcs; ++j )
      {
         /* check for outflows */
         if( arcdata[j].node1 == i )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, flowbalance[i], flow[j], +1.0) );
         }

         /* check for inflows */
         if( arcdata[j].node2 == i )
         {
            SCIP_CALL( SCIPaddCoefLinear(scip, flowbalance[i], flow[j], -1.0) );
         }
      }

      SCIP_CALL( SCIPaddCons(scip, flowbalance[i]) );
   }

   /* create pressure difference constraints and add to problem
    * pressurediff[node1 to node2] = pressure[node1] - pressure[2]
    */
   for( i = 0; i < narcs; ++i )
   {
      SCIP_VAR* vars[3]  = { pressurediff[i], pressure[arcdata[i].node1], pressure[arcdata[i].node2] };
      SCIP_Real coefs[3] = { 1.0, -1.0, 1.0 };

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pressurediff_%s_%s", nodedata[arcdata[i].node1].name, nodedata[arcdata[i].node2].name);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &pressurediffcons[i], name, 3, vars, coefs, 0.0, 0.0) );

      SCIP_CALL( SCIPaddCons(scip, pressurediffcons[i]) );
   }

   /* create pressure loss constraints and add to problem
    * pressure loss on   active arcs: flow(i) *     flow(i)  + coef * pressurediff(i) <= 0.0
    *               on regular pipes: flow(i) * abs(flow(i)) - coef * pressurediff(i) == 0.0,
    * where coef = 96.074830e-15*power(diameter(i)^5/lambda/compressibility/temperatur/length(i)/density
    * and lambda = (2*log10(3.7*diameter(i)/rugosity))^(-2);
    */
   for( i = 0; i < narcs; ++i )
   {
      SCIP_Real coef;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pressureloss_%s_%s", nodedata[arcdata[i].node1].name, nodedata[arcdata[i].node2].name);

      coef = 96.074830e-15 * pow(arcdata[i].diameter, 5.0) * pow(2.0*log10(3.7*arcdata[i].diameter / rugosity), 2.0)
         / compressibility / gastemp / arcdata[i].length / density;

      if( arcdata[i].active )
      {
         /* we can also use an abspower constraint here, because flow(i) is positive for active arcs */
         SCIP_CALL( SCIPcreateConsBasicAbspower(scip, &pressureloss[i], name, flow[i], pressurediff[i], 2.0, 0.0,  coef, -SCIPinfinity(scip), 0.0) );
      }
      else
      {
         SCIP_CALL( SCIPcreateConsBasicAbspower(scip, &pressureloss[i], name, flow[i], pressurediff[i], 2.0, 0.0, -coef, 0.0, 0.0) );
      }

      SCIP_CALL( SCIPaddCons(scip, pressureloss[i]) );
   }

   /* release variables and constraints
    * the problem has them captured, and we do not require them anymore
    */
   for( i = 0; i < nnodes; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &supply[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &pressure[i]) );

      SCIP_CALL( SCIPreleaseCons(scip, &flowbalance[i]) );
   }

   for( i = 0; i < narcs; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &flow[i]) );
      SCIP_CALL( SCIPreleaseVar(scip, &pressurediff[i]) );

      SCIP_CALL( SCIPreleaseCons(scip, &pressurediffcons[i]) );
      SCIP_CALL( SCIPreleaseCons(scip, &pressureloss[i]) );
   }

   return SCIP_OKAY;
}

/* runs gas transportation example */
static
SCIP_RETCODE runGastrans(void)
{
   SCIP* scip;

   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIPinfoMessage(scip, NULL, "*****************************************\n");
   SCIPinfoMessage(scip, NULL, "*    Running Gas Transportation Model   *\n");
   SCIPinfoMessage(scip, NULL, "*****************************************\n");
   SCIPinfoMessage(scip, NULL, "\n");

   SCIP_CALL( setupProblem(scip) );

   SCIPinfoMessage(scip, NULL, "Original problem:\n");
   SCIP_CALL( SCIPprintOrigProblem(scip, NULL, "cip", FALSE) );

   SCIPinfoMessage(scip, NULL, "\n");
   SCIP_CALL( SCIPpresolve(scip) );

   /* SCIPinfoMessage(scip, NULL, "Presolved problem:\n");
   SCIP_CALL( SCIPprintTransProblem(scip, NULL, "cip", FALSE) );
   */

   SCIPinfoMessage(scip, NULL, "\nSolving...\n");
   SCIP_CALL( SCIPsolve(scip) );

   SCIP_CALL( SCIPfreeTransform(scip) );

   /*
   if( SCIPgetNSols(scip) > 0 )
   {
      SCIPinfoMessage(scip, NULL, "\nSolution:\n");
      SCIP_CALL( SCIPprintSol(scip, SCIPgetBestSol(scip), NULL, FALSE) );
   }
   */

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

   retcode = runGastrans();

   /* evaluate return code of the SCIP process */
   if( retcode != SCIP_OKAY )
   {
      /* write error back trace */
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
