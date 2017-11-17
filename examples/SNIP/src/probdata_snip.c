/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   probdata_snip.c
 * @brief  Problem data for snip problem
 * @author Stephen J. Maher
 *
 * This file handles the main problem data used in that project. For more details see \ref SNIP_PROBLEMDATA page.
 *
 * @page SNIP_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the snip problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * \code
 *  ** @brief Problem data which is accessible in all places
 *  *
 *  *   This problem data is used to store the input of the snip instance, all variables which are created, and all
 *  *   constraints.
 *  *
 * struct SCIP_ProbData
 * {
 *    SCIP_VAR**            vars;         **< all exiting variables in the problem *
 *    SCIP_CONS**           conss;        **< set partitioning constraints for each item exactly one *
 *    SCIP_Longint*         weights;      **< array of item weights *
 *    int*                  ids;          **< array of item ids *
 *    int                   nvars;        **< number of generated variables *
 *    int                   varssize;     **< size of the variable array *
 *    int                   nitems;       **< number of items *
 *    SCIP_Longint          capacity;     **< bin capacity *
 * };
 * \endcode
 *
 * The function SCIPprobdataCreate(), which is called in the \ref reader_bpa.c "reader plugin" after the input file was
 * parsed, initializes the problem data structure and creates the problem in the SCIP environment. For this, it creates
 * for each item of the snip problem one set covering constraint and creates an initial set of variables for the
 * packings. Note that the set covering constraints have to have the <code>modifiable</code>-flag set to TRUE. This is
 * necessary to tell the solver that these constraints are not completed yet. This means, during the search new
 * variables/packings might be added.  The solver needs this information because certain reductions are not allowed.
 * See the body of the function SCIPprobdataCreate() for more details.
 *
 * A list of all interface methods can be found in probdata_snip.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "probdata_snip.h"

#include "scip/cons_linear.h"
#include "scip/benders_default.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the snip, all variables which are created, and all
 * constrsaints.
 */
struct SCIP_ProbData
{
   SCIP**                subproblems;        /**< the Benders' decomposition subproblems */
   int                   nscenarios;         /**< the number of scenarios */
   SCIP_Bool             usebenders;         /**< is Benders' decomposition used? */
};



/**@name Local methods
 *
 * @{
 */

/** releases the variables and constraints of the problem */
static
SCIP_RETCODE releaseVariablesAndConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   SCIP_VAR**            sensorvars,         /**< the sensor variable */
   SCIP_VAR***           subsensorvars,      /**< the sensor variables of the subproblems */
   SCIP_VAR***           scenarionodevars,   /**< the node variable per scenario */
   SCIP_CONS*            budgetcons,         /**< the budget constraint */
   SCIP_CONS**           sinkconss,          /**< the sink constraints per scenario */
   SCIP_CONS***          sensorarcconss,     /**< the sensor arc constraints */
   SCIP_CONS***          arcconss,           /**< the arc constraints */
   SCIP_CONS***          probarcconss,       /**< the probability arc constraints */
   int                   narcs,              /**< the number of arcs */
   int                   nnodes,             /**< the number of nodes */
   int                   nsensors,           /**< the number of sensors */
   int                   nscenarios,         /**< the number of scenarios */
   SCIP_Bool             usebenders          /**< is Benders' decomposition used? */
   )
{
   SCIP* varconsscip;
   int i;
   int j;

   assert(scip != NULL);

   varconsscip = scip;

   /* releasing the variables and constraints */
   for( j = 0; j < nscenarios; j++ )
   {
      if( usebenders )
         varconsscip = subproblems[j];

      for( i = 0; i < nnodes; i++ )
         SCIP_CALL( SCIPreleaseVar(varconsscip, &scenarionodevars[i][j]) );
   }

   for( i = 0; i < nsensors; i++ )
      SCIP_CALL( SCIPreleaseVar(scip, &sensorvars[i]) );

   if( usebenders )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         if( usebenders )
            varconsscip = subproblems[j];

         for( i = 0; i < nsensors; i++ )
            SCIP_CALL( SCIPreleaseVar(varconsscip, &subsensorvars[i][j]) );
      }
   }

   SCIP_CALL( SCIPreleaseCons(scip, &budgetcons) );
   for( j = 0; j < nscenarios; j++ )
   {
      if( usebenders )
         varconsscip = subproblems[j];

      for( i = 0; i < nsensors; i++ )
      {
         SCIP_CALL( SCIPreleaseCons(varconsscip, &sensorarcconss[i][j]) );
         SCIP_CALL( SCIPreleaseCons(varconsscip, &probarcconss[i][j]) );
      }

      for( i = 0; i < narcs - nsensors; i++ )
      {
         SCIP_CALL( SCIPreleaseCons(varconsscip, &arcconss[i][j]) );
      }

      SCIP_CALL( SCIPreleaseCons(varconsscip, &sinkconss[j]) );
   }

   return SCIP_OKAY;
}

/** creates the original problem */
static
SCIP_RETCODE createOriginalproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            sensorvars,         /**< the sensor variable */
   SCIP_VAR***           scenarionodevars,   /**< the node variable per scenario */
   SCIP_CONS**           budgetcons,         /**< the budget constraint */
   SCIP_CONS**           sinkconss,          /**< the sink constraints per scenario */
   SCIP_CONS***          sensorarcconss,     /**< the sensor arc constraints */
   SCIP_CONS***          arcconss,           /**< the arc constraints */
   SCIP_CONS***          probarcconss,       /**< the probability arc constraints */
   SCIP_Real*            scenariocost,       /**< the costs for the scenarios */
   SCIP_Real*            probwosensor,       /**< the probability of detection without a sensor */
   SCIP_Real*            intdictwosensor,    /**< the probability of detection without a sensor */
   SCIP_Real**           shortestpaths,      /**< the shortest paths for each scenario */
   int*                  scenarioarcids,     /**< the scenario arc ids */
   int*                  arcids,             /**< the arc ids */
   int*                  intdictarcids,      /**< the interdiction arc ids */
   int*                  nodemapping,        /**< mapping from the node ids to the node index */
   SCIP_Real             budget,             /**< the sensor budget */
   SCIP_Real             multiplier,         /**< the probability multiplier */
   int                   narcs,              /**< the number of arcs */
   int                   nnodes,             /**< the number of nodes */
   int                   nsensors,           /**< the number of sensors */
   int                   nscenarios          /**< the number of scenarios */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_Real coeff;
   int edgetail;
   int edgehead;
   int i;
   int j;
   char name[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "budget");
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, budgetcons, name, 0, NULL, NULL, -SCIPinfinity(scip), budget) );

   SCIP_CALL( SCIPaddCons(scip, (*budgetcons)) );

   /* adds the sink node constraints to the scenario */
   for( i = 0; i < nscenarios; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sinkcons_%d", i);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, 1.0, 1.0) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

      sinkconss[i] = cons;
   }

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < nsensors; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sensorarccons_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, 0.0, SCIPinfinity(scip)) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         sensorarcconss[i][j] = cons;
      }
   }

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < narcs - nsensors; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "arccons_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, 0.0, SCIPinfinity(scip)) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         arcconss[i][j] = cons;
      }
   }

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < nsensors; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "probarccons_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, 0.0, SCIPinfinity(scip)) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         probarcconss[i][j] = cons;
      }
   }

   for( i = 0; i < nsensors; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sensor_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0, 1, 0.0, SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, var) );

      /* storing the variable in the facility variable list */
      sensorvars[i] = var;

      /* adding the sensor vars to the probability arc constraints */
      for( j = 0; j < nscenarios; j++ )
      {
         edgehead = intdictarcids[i] % 1000;
         coeff = intdictwosensor[i]*(1 - multiplier)*shortestpaths[j][nodemapping[edgehead]];
         SCIP_CALL( SCIPaddCoefLinear(scip, probarcconss[i][j], var, coeff) );
      }

      /* adding the sensor vars to the budget constraint */
      SCIP_CALL( SCIPaddCoefLinear(scip, *budgetcons, var, 1.0) );
   }

   for( i = 0; i < nnodes; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         edgetail = (int) scenarioarcids[j]/1000;
         edgehead = scenarioarcids[j] % 1000;

         coeff = 0.0;
         if( i == nodemapping[edgetail] )
            coeff = scenariocost[j];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "node_%d_%d", i, j);
         SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0, SCIPinfinity(scip), coeff, SCIP_VARTYPE_CONTINUOUS) );

         SCIP_CALL( SCIPaddVar(scip, var) );

         /* storing the variable in the facility variable list */
         scenarionodevars[i][j] = var;

         /* adding the node variable to the scenario sink constraints */
         if( i == nodemapping[edgehead] )
            SCIP_CALL( SCIPaddCoefLinear(scip, sinkconss[j], var, 1.0) );
      }
   }

   /* adding the node variables to the constraints */

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < nsensors; i++ )
   {
      edgetail = (int) intdictarcids[i]/1000;
      edgehead = intdictarcids[i] % 1000;

      for( j = 0; j < nscenarios; j++ )
      {
         /* adding variable to interdiction arc constraints */
         SCIP_CALL( SCIPaddCoefLinear(scip, sensorarcconss[i][j], scenarionodevars[nodemapping[edgetail]][j], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, sensorarcconss[i][j], scenarionodevars[nodemapping[edgehead]][j],
               -intdictwosensor[i]*multiplier) );

         /* adding variable to probability arc constraints */
         SCIP_CALL( SCIPaddCoefLinear(scip, probarcconss[i][j], scenarionodevars[nodemapping[edgetail]][j], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, probarcconss[i][j], scenarionodevars[nodemapping[edgehead]][j],-intdictwosensor[i]) );
      }
   }

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < narcs - nsensors; i++ )
   {
      edgetail = (int) arcids[i]/1000;
      edgehead = arcids[i] % 1000;

      for( j = 0; j < nscenarios; j++ )
      {
         /* adding variable to arc constraints */
         SCIP_CALL( SCIPaddCoefLinear(scip, arcconss[i][j], scenarionodevars[nodemapping[edgetail]][j], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(scip, arcconss[i][j], scenarionodevars[nodemapping[edgehead]][j],-probwosensor[i]) );
      }
   }



   return SCIP_OKAY;
}

/** creates the Benders' decomposition master problem */
static
SCIP_RETCODE createMasterproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            sensorvars,         /**< the sensor variable */
   SCIP_CONS**           budgetcons,         /**< the budget constraint */
   SCIP_Real             budget,             /**< the sensor budget */
   int                   nsensors            /**< the number of sensors */
   )
{
   SCIP_VAR* var;
   int i;
   char name[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "budget");
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, budgetcons, name, 0, NULL, NULL, -SCIPinfinity(scip), budget) );

   SCIP_CALL( SCIPaddCons(scip, (*budgetcons)) );

   for( i = 0; i < nsensors; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sensor_%d", i);
      SCIP_CALL( SCIPcreateVarBasic(scip, &var, name, 0, 1, 0.0, SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, var) );

      /* storing the variable in the facility variable list */
      sensorvars[i] = var;

      /* adding the sensor vars to the budget constraint */
      SCIP_CALL( SCIPaddCoefLinear(scip, *budgetcons, var, 1.0) );
   }

   return SCIP_OKAY;
}

/** creates the scenario subproblems */
static
SCIP_RETCODE createSubproblems(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   SCIP_VAR***           subsensorvars,      /**< the sensor variable */
   SCIP_VAR***           scenarionodevars,   /**< the node variable per scenario */
   SCIP_CONS**           sinkconss,          /**< the sink constraints per scenario */
   SCIP_CONS***          sensorarcconss,     /**< the sensor arc constraints */
   SCIP_CONS***          arcconss,           /**< the arc constraints */
   SCIP_CONS***          probarcconss,       /**< the probability arc constraints */
   SCIP_Real*            scenariocost,       /**< the costs for the scenarios */
   SCIP_Real*            probwosensor,       /**< the probability of detection without a sensor */
   SCIP_Real*            intdictwosensor,    /**< the probability of detection without a sensor */
   SCIP_Real**           shortestpaths,      /**< the shortest paths for each scenario */
   int*                  scenarioarcids,     /**< the scenario arc ids */
   int*                  arcids,             /**< the arc ids */
   int*                  intdictarcids,      /**< the interdiction arc ids */
   int*                  nodemapping,        /**< mapping from the node ids to the node index */
   SCIP_Real             multiplier,         /**< the probability multiplier */
   int                   narcs,              /**< the number of arcs */
   int                   nnodes,             /**< the number of nodes */
   int                   nsensors,           /**< the number of sensors */
   int                   nscenarios          /**< the number of scenarios */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_Real coeff;
   int edgetail;
   int edgehead;
   int i;
   int j;
   char name[SCIP_MAXSTRLEN];
   assert(scip != NULL);
   assert(subproblems != NULL);

   /* adds the sink node constraints to the scenario */
   for( i = 0; i < nscenarios; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sinkcons_%d", i);
      SCIP_CALL( SCIPcreateConsBasicLinear(subproblems[i], &cons, name, 0, NULL, NULL, 1.0, 1.0) );

      SCIP_CALL( SCIPaddCons(subproblems[i], cons) );

      sinkconss[i] = cons;
   }

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < nsensors; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sensorarccons_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(subproblems[j], &cons, name, 0, NULL, NULL, 0.0, SCIPinfinity(subproblems[j])) );

         SCIP_CALL( SCIPaddCons(subproblems[j], cons) );

         sensorarcconss[i][j] = cons;
      }
   }

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < narcs - nsensors; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "arccons_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(subproblems[j], &cons, name, 0, NULL, NULL, 0.0, SCIPinfinity(subproblems[j])) );

         SCIP_CALL( SCIPaddCons(subproblems[j], cons) );

         arcconss[i][j] = cons;
      }
   }

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < nsensors; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "probarccons_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(subproblems[j], &cons, name, 0, NULL, NULL, 0.0, SCIPinfinity(subproblems[j])) );

         SCIP_CALL( SCIPaddCons(subproblems[j], cons) );

         probarcconss[i][j] = cons;
      }
   }

   for( i = 0; i < nsensors; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sensor_%d", i);
         SCIP_CALL( SCIPcreateVarBasic(subproblems[j], &var, name, 0, 1, 0.0, SCIP_VARTYPE_BINARY) );

         SCIP_CALL( SCIPaddVar(subproblems[j], var) );

         /* storing the variable in the facility variable list */
         subsensorvars[i][j] = var;

         /* adding the sensor vars to the probability arc constraints */
         edgehead = intdictarcids[i] % 1000;
         coeff = intdictwosensor[i]*(1 - multiplier)*shortestpaths[j][nodemapping[edgehead]];
         SCIP_CALL( SCIPaddCoefLinear(subproblems[j], probarcconss[i][j], var, coeff) );
      }
   }

   for( i = 0; i < nnodes; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         edgetail = (int) scenarioarcids[j]/1000;
         edgehead = scenarioarcids[j] % 1000;

         coeff = 0.0;
         if( i == nodemapping[edgetail] )
            coeff = scenariocost[j] * 1000;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "node_%d_%d", i, j);
         SCIP_CALL( SCIPcreateVarBasic(subproblems[j], &var, name, 0, SCIPinfinity(subproblems[j]), coeff, SCIP_VARTYPE_CONTINUOUS) );

         SCIP_CALL( SCIPaddVar(subproblems[j], var) );

         /* storing the variable in the facility variable list */
         scenarionodevars[i][j] = var;

         /* adding the node variable to the scenario sink constraints */
         if( i == nodemapping[edgehead] )
            SCIP_CALL( SCIPaddCoefLinear(subproblems[j], sinkconss[j], var, 1.0) );
      }
   }

   /* adding the node variables to the constraints */

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < nsensors; i++ )
   {
      edgetail = (int) intdictarcids[i]/1000;
      edgehead = intdictarcids[i] % 1000;

      for( j = 0; j < nscenarios; j++ )
      {
         /* adding variable to interdiction arc constraints */
         SCIP_CALL( SCIPaddCoefLinear(subproblems[j], sensorarcconss[i][j], scenarionodevars[nodemapping[edgetail]][j], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(subproblems[j], sensorarcconss[i][j], scenarionodevars[nodemapping[edgehead]][j],
               -intdictwosensor[i]*multiplier) );

         /* adding variable to probability arc constraints */
         SCIP_CALL( SCIPaddCoefLinear(subproblems[j], probarcconss[i][j], scenarionodevars[nodemapping[edgetail]][j], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(subproblems[j], probarcconss[i][j], scenarionodevars[nodemapping[edgehead]][j],-intdictwosensor[i]) );
      }
   }

   /* adds the sensor arc constraints to the scenario */
   for( i = 0; i < narcs - nsensors; i++ )
   {
      edgetail = (int) arcids[i]/1000;
      edgehead = arcids[i] % 1000;

      for( j = 0; j < nscenarios; j++ )
      {
         /* adding variable to arc constraints */
         SCIP_CALL( SCIPaddCoefLinear(subproblems[j], arcconss[i][j], scenarionodevars[nodemapping[edgetail]][j], 1.0) );
         SCIP_CALL( SCIPaddCoefLinear(subproblems[j], arcconss[i][j], scenarionodevars[nodemapping[edgehead]][j],-probwosensor[i]) );
      }
   }

   return SCIP_OKAY;
}


/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   int                   nscenarios,         /**< the number of scenarios */
   SCIP_Bool             usebenders          /**< whether Benders' decomposition is used */
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );

   /* copying the subproblem information */
   if( usebenders )
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->subproblems, subproblems, nscenarios) );

   (*probdata)->nscenarios = nscenarios;
   (*probdata)->usebenders = usebenders;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* freeing the subproblem information */
   if( (*probdata)->usebenders )
   {
      for( i = 0; i < (*probdata)->nscenarios; i++ )
         SCIP_CALL( SCIPfree(&(*probdata)->subproblems[i]) );
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->subproblems, (*probdata)->nscenarios);
   }

   /* free probdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigSnip)
{
   assert(scip != NULL);
   assert(probdata != NULL);

   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransSnip)
{
   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransSnip)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   SCIP_Real*            scenariocost,       /**< the costs for the scenarios */
   SCIP_Real*            probwosensor,       /**< the probability of detection without a sensor */
   SCIP_Real*            intdictwosensor,    /**< the probability of detection without a sensor */
   SCIP_Real**           shortestpaths,      /**< the shortest paths for each scenario */
   int*                  scenarioarcids,     /**< the scenario arc ids */
   int*                  arcids,             /**< the arc ids */
   int*                  intdictarcids,      /**< the interdiction arc ids */
   int*                  nodemapping,        /**< mapping from the node ids to the node index */
   SCIP_Real             budget,             /**< the sensor budget */
   SCIP_Real             multiplier,         /**< the probability multiplier */
   int                   narcs,              /**< the number of arcs */
   int                   nnodes,             /**< the number of nodes */
   int                   nsensors,           /**< the number of sensors */
   int                   nscenarios,         /**< the number of scenarios */
   SCIP_Bool             usebenders          /**< will Benders' decomposition be used to solve the problem */
   )
{
   SCIP** subproblems;
   SCIP_PROBDATA* probdata;
   SCIP_CONS* budgetcons;
   SCIP_CONS** sinkconss;
   SCIP_CONS*** sensorarcconss;
   SCIP_CONS*** arcconss;
   SCIP_CONS*** probarcconss;
   SCIP_VAR** sensorvars;
   SCIP_VAR*** subsensorvars;
   SCIP_VAR*** scenarionodevars;
   int i;

   assert(scip != NULL);

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigSnip) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransSnip) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransSnip) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &sinkconss, nscenarios) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sensorarcconss, nsensors) );
   for( i = 0; i < nsensors; i++ )
      SCIP_CALL( SCIPallocBufferArray(scip, &sensorarcconss[i], nscenarios) );
   SCIP_CALL( SCIPallocBufferArray(scip, &arcconss, narcs - nsensors) );
   for( i = 0; i < narcs - nsensors; i++ )
      SCIP_CALL( SCIPallocBufferArray(scip, &arcconss[i], nscenarios) );
   SCIP_CALL( SCIPallocBufferArray(scip, &probarcconss, nsensors) );
   for( i = 0; i < nsensors; i++ )
      SCIP_CALL( SCIPallocBufferArray(scip, &probarcconss[i], nscenarios) );

   SCIP_CALL( SCIPallocBufferArray(scip, &sensorvars, nsensors) );
   SCIP_CALL( SCIPallocBufferArray(scip, &scenarionodevars, nnodes) );
   for( i = 0; i < nnodes; i++ )
      SCIP_CALL( SCIPallocBufferArray(scip, &scenarionodevars[i], nscenarios) );

   subproblems = NULL;
   subsensorvars = NULL;

   if( usebenders )
   {
      char subprobname[SCIP_MAXSTRLEN];

      /* allocting the memory for the subproblem specific information */
      SCIP_CALL( SCIPallocBufferArray(scip, &subproblems, nscenarios) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subsensorvars, nsensors) );
      for( i = 0; i < nsensors; i++ )
         SCIP_CALL( SCIPallocBufferArray(scip, &subsensorvars[i], nscenarios) );

      /* creating the subproblems */
      for( i = 0; i < nscenarios; i++ )
      {
         SCIP_CALL( SCIPcreate(&subproblems[i]) );

         /* include default SCIP plugins */
         SCIP_CALL( SCIPincludeDefaultPlugins(subproblems[i]) );

         (void) SCIPsnprintf(subprobname, SCIP_MAXSTRLEN, "sub_%s_%d", probname, i);
         SCIP_CALL( SCIPcreateProbBasic(subproblems[i], subprobname) );
      }

      /* creating the master problem */
      SCIP_CALL( createMasterproblem(scip, sensorvars, &budgetcons, budget, nsensors) );
      SCIP_CALL( createSubproblems(scip, subproblems, subsensorvars, scenarionodevars, sinkconss,
            sensorarcconss, arcconss, probarcconss, scenariocost, probwosensor, intdictwosensor,
            shortestpaths, scenarioarcids, arcids, intdictarcids, nodemapping, multiplier, narcs, nnodes, nsensors, nscenarios) );

      /* including the Benders' decomposition plugin */
      SCIP_CALL( SCIPcreateBendersDefault(scip, subproblems, nscenarios) );
   }
   else
   {
      /* creating the original problem */
      SCIP_CALL( createOriginalproblem(scip, sensorvars, scenarionodevars, &budgetcons, sinkconss, sensorarcconss,
            arcconss, probarcconss, scenariocost, probwosensor, intdictwosensor, shortestpaths, scenarioarcids,
            arcids, intdictarcids, nodemapping, budget, multiplier, narcs, nnodes, nsensors, nscenarios) );
   }

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, subproblems, nscenarios, usebenders) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   SCIP_CALL( releaseVariablesAndConstraints(scip, subproblems, sensorvars, subsensorvars, scenarionodevars,
         budgetcons, sinkconss, sensorarcconss, arcconss, probarcconss, narcs, nnodes, nsensors, nscenarios,
         usebenders) );

   /* free local buffer arrays */
   if( usebenders )
   {
      SCIPfreeBufferArray(scip, &subproblems);

      for( i = nsensors - 1; i >= 0; i-- )
         SCIPfreeBufferArray(scip, &subsensorvars[i]);
      SCIPfreeBufferArray(scip, &subsensorvars);
   }

   for( i = nnodes - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &scenarionodevars[i]);
   SCIPfreeBufferArray(scip, &scenarionodevars);
   SCIPfreeBufferArray(scip, &sensorvars);


   for( i = nsensors - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &probarcconss[i]);
   SCIPfreeBufferArray(scip, &probarcconss);
   for( i = narcs - nsensors - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &arcconss[i]);
   SCIPfreeBufferArray(scip, &arcconss);
   for( i = nsensors - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &sensorarcconss[i]);
   SCIPfreeBufferArray(scip, &sensorarcconss);
   SCIPfreeBufferArray(scip, &sinkconss);



   return SCIP_OKAY;
}

/**@} */
