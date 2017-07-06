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

/**@file   probdata_cap.c
 * @brief  Problem data for cap problem
 * @author Stephen J. Maher
 *
 * This file handles the main problem data used in that project. For more details see \ref CAP_PROBLEMDATA page.
 *
 * @page CAP_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the cap problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * \code
 *  ** @brief Problem data which is accessible in all places
 *  *
 *  *   This problem data is used to store the input of the cap instance, all variables which are created, and all
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
 * for each item of the cap problem one set covering constraint and creates an initial set of variables for the
 * packings. Note that the set covering constraints have to have the <code>modifiable</code>-flag set to TRUE. This is
 * necessary to tell the solver that these constraints are not completed yet. This means, during the search new
 * variables/packings might be added.  The solver needs this information because certain reductions are not allowed.
 * See the body of the function SCIPprobdataCreate() for more details.
 *
 * A list of all interface methods can be found in probdata_cap.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "probdata_cap.h"
#include "vardata_cap.h"

#include "scip/cons_linear.h"
#include "scip/benders_default.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the cap, all variables which are created, and all
 * constrsaints.
 */
struct SCIP_ProbData
{
   SCIP**                subproblems;        /**< the Benders' decomposition subproblems */
   SCIP_VAR**            facilityvars;       /**< all variables representing facilities */
   SCIP_VAR***           subfacilityvars;    /**< duplicates of the facility variables in the subproblems */
   SCIP_VAR****          customervars;       /**< all variables representing the satisfaction of demand per scenario */
   SCIP_CONS***          capconss;           /**< capacity constraints per facility per scenario */
   SCIP_CONS***          demandconss;        /**< demand constraints per customer per scenario */
   SCIP_CONS*            sufficientcap;      /**< ensuring sufficient capacity is provided to satisfy demand (relatively complete recourse) */
   SCIP_Real**           costs;              /**< the transportation costs to a customer from a facility */
   SCIP_Real**           demands;            /**< the customer demands per scenario */
   SCIP_Real*            capacity;           /**< the capacity of each facility */
   SCIP_Real*            fixedcost;          /**< the fixed cost of openning each facility */
   int                   ncustomers;         /**< the number of customers */
   int                   nfacilities;        /**< the number of facilities */
   int                   nscenarios;         /**< the number of scenarios */
   SCIP_Bool             usebenders;         /**< whether Benders' decomposition is used */
};



/**@name Local methods
 *
 * @{
 */

/** creates the original problem */
static
SCIP_RETCODE createOriginalproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            facilityvars,       /**< all variables representing facilities */
   SCIP_VAR****          customervars,       /**< all variables representing the satisfaction of demand */
   SCIP_CONS***          capconss,           /**< capacity constraints per facility */
   SCIP_CONS***          demandconss,        /**< demand constraints per customer */
   SCIP_CONS**           sufficientcap,      /**< ensuring sufficient capacity is provided to satisfy demand */
   SCIP_Real**           costs,              /**< the transportation costs from a facility to a customer */
   SCIP_Real**           demands,            /**< the customer demands */
   SCIP_Real*            capacity,           /**< the capacity of each facility */
   SCIP_Real*            fixedcost,          /**< the fixed cost of openning a facility */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nscenarios          /**< the number of scenarios */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_VARDATA* vardata;
   SCIP_Real maxdemand;
   SCIP_Real coeff;
   int i;
   int j;
   int k;
   char name[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   /* adding the sufficient capacity constraints */
   maxdemand = 0;
   for( i = 0; i < nscenarios; i++)
   {
      SCIP_Real sumdemand = 0;
      for( j = 0; j < ncustomers; j++ )
         sumdemand += demands[j][i];

      if( sumdemand > maxdemand )
         maxdemand = sumdemand;
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sufficientcapacity");
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, sufficientcap, name, 0, NULL, NULL, maxdemand, SCIPinfinity(scip)) );

   SCIP_CALL( SCIPaddCons(scip, (*sufficientcap)) );

   /* adds the capacity constraints to the scenario */
   for( i = 0; i < nfacilities; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "capacity_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         capconss[i][j] = cons;
      }
   }

   /* adds the demand constraints to the scenario */
   for( i = 0; i < ncustomers; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "demand_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, demands[i][j], SCIPinfinity(scip)) );

         SCIP_CALL( SCIPaddCons(scip, cons) );

         demandconss[i][j] = cons;
      }
   }

   for( i = 0; i < nfacilities; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "facility_%d", i);
      SCIP_CALL( SCIPcreateVarCAP(scip, &var, name, 0, 1, fixedcost[i], SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, var) );

      /* storing the variable in the facility variable list */
      facilityvars[i] = var;

      /* creates the variable data */
      SCIP_CALL( SCIPvardataCreateCAP(scip, &vardata, MASTER, nscenarios) );

      /* add the variable data to the variable */
      SCIPvarSetData(var, vardata);

      /* adding the variable to the capacity constriants */
      for( j = 0; j < nscenarios; j++ )
         SCIP_CALL( SCIPaddCoefLinear(scip, capconss[i][j], var, -capacity[i]) );

      /* adding the variable to the sufficient capacity constriants */
      SCIP_CALL( SCIPaddCoefLinear(scip, (*sufficientcap), var, capacity[i]) );

      /* releases the variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   /* adding the customer variables to the scenario */
   for( i = 0; i < ncustomers; i++ )
   {
      for( j = 0; j < nfacilities; j++ )
      {
         for( k = 0; k < nscenarios; k++ )
         {
            coeff = costs[i][j]/(SCIP_Real)nscenarios;
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "customer(%d,%d,%d)", i, j, k);
            SCIP_CALL( SCIPcreateVarCAP(scip, &var, name, 0, SCIPinfinity(scip), coeff, SCIP_VARTYPE_CONTINUOUS) );

            SCIP_CALL( SCIPaddVar(scip, var) );

            /* storing the customer variable in the list */
            customervars[i][j][k] = var;

            /* creates the variable data */
            SCIP_CALL( SCIPvardataCreateCAP(scip, &vardata, SUBPROB, 1) );

            /* add the variable data to the variable */
            SCIPvarSetData(var, vardata);

            if( costs[i][j] > 0 )
            {
               /* adding the variable to the capacity constriants */
               SCIP_CALL( SCIPaddCoefLinear(scip, capconss[j][k], customervars[i][j][k], 1.0) );

               /* adding the variable to the demand constraints */
               SCIP_CALL( SCIPaddCoefLinear(scip, demandconss[i][k], customervars[i][j][k], 1.0) );
            }

            /* releases the variable */
            SCIP_CALL( SCIPreleaseVar(scip, &var) );
         }
      }
   }

   return SCIP_OKAY;
}

/** creates the Benders' decomposition master problem */
static
SCIP_RETCODE createMasterproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            facilityvars,       /**< all variables representing facilities */
   SCIP_CONS**           sufficientcap,      /**< ensuring sufficient capacity is provided to satisfy demand */
   SCIP_Real*            capacity,           /**< the capacity of each facility */
   SCIP_Real*            fixedcost,          /**< the fixed cost of openning a facility */
   SCIP_Real**           demands,            /**< the customer demands */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nscenarios          /**< the number of scenarios */
   )
{
   SCIP_VAR* var;
   SCIP_VARDATA* vardata;
   SCIP_Real maxdemand;
   int i;
   int j;
   char name[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   /* adding the sufficient capacity constraints */
   maxdemand = 0;
   for( i = 0; i < nscenarios; i++)
   {
      SCIP_Real sumdemand = 0;
      for( j = 0; j < ncustomers; j++ )
         sumdemand += demands[j][i];

      if( sumdemand > maxdemand )
         maxdemand = sumdemand;
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sufficientcapacity");
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, sufficientcap, name, 0, NULL, NULL, maxdemand, SCIPinfinity(scip)) );

   SCIP_CALL( SCIPaddCons(scip, (*sufficientcap)) );


   /* adding the facility variables */
   for( i = 0; i < nfacilities; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "facility_%d", i);
      SCIP_CALL( SCIPcreateVarCAP(scip, &var, name, 0, 1, fixedcost[i], SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, var) );

      /* storing the variable in the facility variable list */
      facilityvars[i] = var;

      /* creates the variable data */
      SCIP_CALL( SCIPvardataCreateCAP(scip, &vardata, MASTER, nscenarios) );

      /* add the variable data to the variable */
      SCIPvarSetData(var, vardata);

      /* adding the variable to the sufficient capacity constriants */
      SCIP_CALL( SCIPaddCoefLinear(scip, (*sufficientcap), var, capacity[i]) );

      /* releases the variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   return SCIP_OKAY;
}

/** creates the scenario subproblems */
static
SCIP_RETCODE createSubproblems(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP**                subproblems,        /**< the Benders' decomposition subproblems */
   SCIP_VAR**            facilityvars,       /**< all variables representing facilities */
   SCIP_VAR***           subfacilityvars,    /**< the copies of the facility variables in the subproblems */
   SCIP_VAR****          customervars,       /**< all variables representing the satisfaction of demand */
   SCIP_CONS***          capconss,           /**< capacity constraints per facility */
   SCIP_CONS***          demandconss,        /**< demand constraints per customer */
   SCIP_Real**           costs,              /**< the transportation costs from a facility to a customer */
   SCIP_Real**           demands,            /**< the customer demands */
   SCIP_Real*            capacity,           /**< the capacity of each facility */
   SCIP_Real*            fixedcost,          /**< the fixed cost of openning a facility */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nscenarios          /**< the number of scenarios */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_VARDATA* vardata;
   SCIP_Real coeff;
   int i;
   int j;
   int k;
   char name[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   /* adds the capacity constraints to the scenario */
   for( i = 0; i < nfacilities; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "capacity_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(subproblems[j], &cons, name, 0, NULL, NULL, -SCIPinfinity(subproblems[j]), 0.0) );

         SCIP_CALL( SCIPaddCons(subproblems[j], cons) );

         capconss[i][j] = cons;
      }
   }

   /* adds the demand constraints to the scenario */
   for( i = 0; i < ncustomers; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "demand_%d_%d", i, j);
         SCIP_CALL( SCIPcreateConsBasicLinear(subproblems[j], &cons, name, 0, NULL, NULL, demands[i][j], SCIPinfinity(subproblems[j])) );

         SCIP_CALL( SCIPaddCons(subproblems[j], cons) );

         demandconss[i][j] = cons;
      }
   }

   for( i = 0; i < nfacilities; i++ )
   {
      for( j = 0; j < nscenarios; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "facility_%d", i);
         SCIP_CALL( SCIPcreateVarCAP(subproblems[j], &var, name, 0.0, 1.0, 0.0, SCIP_VARTYPE_CONTINUOUS) );

         SCIP_CALL( SCIPaddVar(subproblems[j], var) );

         /* storing the variable in the facility variable list */
         subfacilityvars[i][j] = var;

         /* creates the variable data */
         SCIP_CALL( SCIPvardataCreateCAP(subproblems[j], &vardata, SUBPROB, 1) );

         /* adds the master variable to the variable mapping */
         SCIPvardataAddVarMapping(vardata, facilityvars[i], -1);

         /* add the variable data to the variable */
         SCIPvarSetData(var, vardata);

         /* getting the variable data for the master variable */
         vardata = SCIPvarGetData(facilityvars[i]);

         /* adds the subproblem variable to the variable mapping */
         SCIPvardataAddVarMapping(vardata, var, i);

         /* adding the variable to the capacity constriants */
         SCIP_CALL( SCIPaddCoefLinear(subproblems[j], capconss[i][j], subfacilityvars[i][j], -capacity[i]) );

         /* releases the variable */
         SCIP_CALL( SCIPreleaseVar(subproblems[j], &var) );
      }
   }

   /* adding the customer variables to the scenario */
   for( i = 0; i < ncustomers; i++ )
   {
      for( j = 0; j < nfacilities; j++ )
      {
         for( k = 0; k < nscenarios; k++ )
         {
            coeff = costs[i][j]/(SCIP_Real)nscenarios;
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "customer(%d,%d,%d)", i, j, k);
            SCIP_CALL( SCIPcreateVarCAP(subproblems[k], &var, name, 0, SCIPinfinity(subproblems[k]), coeff, SCIP_VARTYPE_CONTINUOUS) );

            SCIP_CALL( SCIPaddVar(subproblems[k], var) );

            /* storing the customer variable in the list */
            customervars[i][j][k] = var;

            /* creates the variable data */
            SCIP_CALL( SCIPvardataCreateCAP(subproblems[k], &vardata, SUBPROB, 1) );

            /* add the variable data to the variable */
            SCIPvarSetData(var, vardata);

            if( costs[i][j] > 0 )
            {
               /* adding the variable to the capacity constriants */
               SCIP_CALL( SCIPaddCoefLinear(subproblems[k], capconss[j][k], customervars[i][j][k], 1.0) );

               /* adding the variable to the demand constraints */
               SCIP_CALL( SCIPaddCoefLinear(subproblems[k], demandconss[i][k], customervars[i][j][k], 1.0) );
            }

            /* releases the variable */
            SCIP_CALL( SCIPreleaseVar(subproblems[k], &var) );
         }
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
   SCIP_VAR**            facilityvars,       /**< all variables representing facilities */
   SCIP_VAR***           subfacilityvars,    /**< the copies of the facility variables in the subproblems */
   SCIP_VAR****          customervars,       /**< all variables representing the satisfaction of demand */
   SCIP_CONS***          capconss,           /**< capacity constraints per facility per scenario */
   SCIP_CONS***          demandconss,        /**< demand constraints per customer per scenario */
   SCIP_CONS*            sufficientcap,      /**< ensuring sufficient capacity is provided to satisfy demand */
   SCIP_Real**           costs,              /**< the transportation costs to a customer from a facility */
   SCIP_Real**           demands,            /**< the customer demands per scenario */
   SCIP_Real*            capacity,           /**< the capacity of each facility */
   SCIP_Real*            fixedcost,          /**< the fixed cost of openning a facility */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nscenarios,         /**< the number of scenarios */
   SCIP_Bool             usebenders          /**< whether Benders' decomposition is used */
   )
{
   int i;
   int j;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );

   /* copying the subproblem information */
   if( usebenders )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->subproblems, subproblems, nscenarios) );

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->subfacilityvars, nfacilities) );
      for( i = 0; i < nfacilities; i++ )
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->subfacilityvars[i], subfacilityvars[i], nscenarios) );
   }

   /* copy variable arrays */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->facilityvars, facilityvars, nfacilities) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->customervars, ncustomers) );
   for( i = 0; i < ncustomers; i++ )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->customervars[i], nfacilities) );
      for( j = 0; j < nfacilities; j++ )
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->customervars[i][j], customervars[i][j],
               nscenarios) );
   }

   /* duplicate the constraint arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->capconss, nfacilities) );
   for( i = 0; i < nfacilities; i++ )
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->capconss[i], capconss[i], nscenarios) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->demandconss, ncustomers) );
   for( i = 0; i < ncustomers; i++ )
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->demandconss[i], demandconss[i], nscenarios) );

   /* duplicate the data arrays */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->demands, ncustomers) );
   for( i = 0; i < ncustomers; i++ )
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->demands[i], demands[i], nscenarios) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->costs, ncustomers) );
   for( i = 0; i < ncustomers; i++ )
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->costs[i], costs[i], nfacilities) );

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->capacity, capacity, nfacilities) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->fixedcost, fixedcost, nfacilities) );

   (*probdata)->sufficientcap = sufficientcap;
   (*probdata)->ncustomers = ncustomers;
   (*probdata)->nfacilities = nfacilities;
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
   int j;
   int k;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* release all variables */
   for( i = 0; i < (*probdata)->nfacilities; i++ )
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->facilityvars[i]) );

   for( i = 0; i < (*probdata)->nscenarios; i++ )
   {
      SCIP* varscip;
      if( (*probdata)->usebenders )
         varscip = (*probdata)->subproblems[i];
      else
         varscip = scip;

      for( j = 0; j < (*probdata)->nfacilities; j++ )
      {
         for( k = 0; k < (*probdata)->ncustomers; k++ )
           SCIP_CALL( SCIPreleaseVar(varscip, &(*probdata)->customervars[k][j][i]) );
      }
   }

   /* release all constraints */
   for( i = 0; i < (*probdata)->nscenarios; ++i )
   {
      SCIP* consscip;
      if( (*probdata)->usebenders )
         consscip = (*probdata)->subproblems[i];
      else
         consscip = scip;

      for( j = 0; j < (*probdata)->ncustomers; j++ )
         SCIP_CALL( SCIPreleaseCons(consscip, &(*probdata)->demandconss[j][i]) );
   }

   for( i = 0; i < (*probdata)->nscenarios; ++i )
   {
      SCIP* consscip;
      if( (*probdata)->usebenders )
         consscip = (*probdata)->subproblems[i];
      else
         consscip = scip;

      for( j = 0; j < (*probdata)->nfacilities; ++j )
         SCIP_CALL( SCIPreleaseCons(consscip, &(*probdata)->capconss[j][i]) );
   }

   SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->sufficientcap) );

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->fixedcost, (*probdata)->nfacilities);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->capacity, (*probdata)->nfacilities);

   for( i = (*probdata)->nfacilities - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->costs[i], (*probdata)->ncustomers);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->costs, (*probdata)->nfacilities);


   for( i = (*probdata)->ncustomers - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->demands[i], (*probdata)->nfacilities);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->demands, (*probdata)->ncustomers);

   /* freeing the constraint memory arrays */
   for( i = (*probdata)->ncustomers - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->demandconss[i], (*probdata)->nscenarios);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->demandconss, (*probdata)->ncustomers);

   for( i = (*probdata)->nfacilities - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->capconss[i], (*probdata)->nscenarios);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->capconss, (*probdata)->nfacilities);

   /* freeing the variable memory arrays */
   for( i = (*probdata)->ncustomers - 1; i >= 0; i-- )
   {
      for( j = (*probdata)->nfacilities - 1; j >= 0; j-- )
         SCIPfreeBlockMemoryArray(scip, &(*probdata)->customervars[i][j], (*probdata)->nscenarios);

      SCIPfreeBlockMemoryArray(scip, &(*probdata)->customervars[i], (*probdata)->nfacilities);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->customervars, (*probdata)->ncustomers);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->facilityvars, (*probdata)->nfacilities);

   /* freeing the subproblem information */
   if( (*probdata)->usebenders )
   {
      /* freeing the sub facility variables */
      for( i = 0; i < (*probdata)->nscenarios; i++ )
      {
         for( j = 0; j < (*probdata)->nfacilities; j++ )
            SCIP_CALL( SCIPreleaseVar((*probdata)->subproblems[i], &(*probdata)->subfacilityvars[j][i]) );
      }

      for( i = (*probdata)->nfacilities - 1; i >= 0; i-- )
         SCIPfreeBlockMemoryArray(scip, &(*probdata)->subfacilityvars[i], (*probdata)->nscenarios);

      SCIPfreeBlockMemoryArray(scip, &(*probdata)->subfacilityvars, (*probdata)->nfacilities);


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
SCIP_DECL_PROBDELORIG(probdelorigCap)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransCap)
{
   /* int i; */

   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->subproblems, sourcedata->facilityvars,
         sourcedata->subfacilityvars, sourcedata->customervars, sourcedata->capconss, sourcedata->demandconss,
         sourcedata->sufficientcap, sourcedata->costs, sourcedata->demands, sourcedata->capacity,
         sourcedata->fixedcost, sourcedata->ncustomers, sourcedata->nfacilities, sourcedata->nscenarios,
         sourcedata->usebenders) );

#if 0
   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->ndemandconss, (*targetdata)->demandconss,
         (*targetdata)->demandconss) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->ncapconss, (*targetdata)->capconss,
         (*targetdata)->capconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nfacilities, (*targetdata)->facilityvars,
         (*targetdata)->facilityvars) );

   for( i = 0; i < (*targetdata)->ncustomers; i++ )
      SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->ncustomervars[i], (*targetdata)->customervars[i],
            (*targetdata)->customervars[i]) );
#endif

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransCap)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

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
   SCIP_Real**           costs,              /**< the transportation costs from a facility to a customer */
   SCIP_Real**           demands,            /**< the customer demands */
   SCIP_Real*            capacity,           /**< the capacity of each facility */
   SCIP_Real*            fixedcost,          /**< the fixed cost of openning a facility */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nscenarios,       /**< the number of Benders' decomposition scenarios */
   SCIP_Bool             usebenders          /**< will Benders' decomposition be used to solve the problem */
   )
{
   SCIP** subproblems;
   SCIP_PROBDATA* probdata;
   SCIP_CONS*** demandconss;
   SCIP_CONS*** capconss;
   SCIP_CONS* sufficientcap;
   SCIP_VAR** facilityvars;
   SCIP_VAR*** subfacilityvars;
   SCIP_VAR**** customervars;
   int i;
   int j;

   assert(scip != NULL);

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigCap) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransCap) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransCap) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &demandconss, ncustomers) );
   for( i = 0; i < ncustomers; i++ )
      SCIP_CALL( SCIPallocBufferArray(scip, &demandconss[i], nscenarios) );
   SCIP_CALL( SCIPallocBufferArray(scip, &capconss, nfacilities) );
   for( i = 0; i < nfacilities; i++ )
      SCIP_CALL( SCIPallocBufferArray(scip, &capconss[i], nscenarios) );

   SCIP_CALL( SCIPallocBufferArray(scip, &facilityvars, nfacilities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &customervars, ncustomers) );
   for( i = 0; i < ncustomers; i++ )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &customervars[i], nfacilities) );
      for( j = 0; j < nfacilities; j++ )
         SCIP_CALL( SCIPallocBufferArray(scip, &customervars[i][j], nscenarios) );
   }

   sufficientcap = NULL;

   subproblems = NULL;
   subfacilityvars = NULL;

   if( usebenders )
   {
      char subprobname[SCIP_MAXSTRLEN];

      /* allocting the memory for the subproblem specific information */
      SCIP_CALL( SCIPallocBufferArray(scip, &subproblems, nscenarios) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subfacilityvars, nfacilities) );
      for( i = 0; i < nfacilities; i++ )
         SCIP_CALL( SCIPallocBufferArray(scip, &subfacilityvars[i], nscenarios) );

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
      SCIP_CALL( createMasterproblem(scip, facilityvars, &sufficientcap, capacity, fixedcost, demands, ncustomers,
            nfacilities, nscenarios) );
      SCIP_CALL( createSubproblems(scip, subproblems, facilityvars, subfacilityvars, customervars, capconss,
            demandconss, costs, demands, capacity, fixedcost, ncustomers, nfacilities, nscenarios) );

      /* including the Benders' decomposition plugin */
      SCIP_CALL( SCIPincludeBendersDefault(scip, subproblems, nscenarios) );
   }
   else
   {
      /* creating the original problem */
      SCIP_CALL( createOriginalproblem(scip, facilityvars, customervars, capconss, demandconss, &sufficientcap, costs,
            demands, capacity, fixedcost, ncustomers, nfacilities, nscenarios) );
   }

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, subproblems, facilityvars, subfacilityvars, customervars, capconss,
         demandconss, sufficientcap, costs, demands, capacity, fixedcost, ncustomers, nfacilities, nscenarios, usebenders) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   /* free local buffer arrays */
   if( usebenders )
   {
      SCIPfreeBufferArray(scip, &subproblems);

      for( i = nfacilities - 1; i >= 0; i-- )
         SCIPfreeBufferArray(scip, &subfacilityvars[i]);
      SCIPfreeBufferArray(scip, &subfacilityvars);
   }

   for( i = ncustomers - 1; i >= 0; i-- )
   {
      for( j = nfacilities - 1; j >= 0; j-- )
         SCIPfreeBufferArray(scip, &customervars[i][j]);
      SCIPfreeBufferArray(scip, &customervars[i]);
   }
   SCIPfreeBufferArray(scip, &customervars);
   SCIPfreeBufferArray(scip, &facilityvars);

   for( i = nfacilities - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &capconss[i]);
   SCIPfreeBufferArray(scip, &capconss);

   for( i = ncustomers - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &demandconss[i]);
   SCIPfreeBufferArray(scip, &demandconss);

   return SCIP_OKAY;
}

/** returns the number of facilities */
int SCIPprobdataGetNFacilities(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->nfacilities;
}

/** returns the number of customers  */
int SCIPprobdataGetNCustomers(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->ncustomers;
}

/** returns the facility variables */
SCIP_VAR** SCIPprobdataGetFacilityVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->facilityvars;
}

#if 0
/** returns the capacity of the facilities */
SCIP_Longint SCIPprobdataGetCapacity(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->capacity;
}

/** returns the customer costs*/
SCIP_Real** SCIPprobdataGetCosts(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->costs;
}

/** returns the customer demands */
SCIP_Real* SCIPprobdataGetDemands(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->demands;
}

/** returns the facilitty fixed cost */
SCIP_Real SCIPprobdataGetFixedCost(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert(probdata != NULL);

   return probdata->fixedcost;
}
#endif

/**@} */
