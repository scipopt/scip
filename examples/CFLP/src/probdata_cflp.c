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

/**@file   probdata_cflp.c
 * @brief  Problem data for cflp problem
 * @author Stephen J. Maher
 *
 * This file handles the main problem data used in that project. For more details see \ref CFLP_PROBLEMDATA page.
 *
 * @page CFLP_PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the cflp problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * \code
 *  ** @brief Problem data which is accessible in all places
 *  *
 *  *   This problem data is used to store the input of the cflp instance, all variables which are created, and all
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
 * for each item of the cflp problem one set covering constraint and creates an initial set of variables for the
 * packings. Note that the set covering constraints have to have the <code>modifiable</code>-flag set to TRUE. This is
 * necessary to tell the solver that these constraints are not completed yet. This means, during the search new
 * variables/packings might be added.  The solver needs this information because certain reductions are not allowed.
 * See the body of the function SCIPprobdataCreate() for more details.
 *
 * A list of all interface methods can be found in probdata_cflp.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "probdata_cflp.h"
#include "vardata_cflp.h"

#include "scip/cons_linear.h"
#include "scip/scip.h"
#include "scip/scipdefplugins.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the cflp, all variables which are created, and all
 * constrsaints.
 */
struct SCIP_ProbData
{
   SCIP_VAR**            facilityvars; /**< all variables representing facilities */
   SCIP_VAR***           customervars; /**< all variables representing the satisfaction of demand */
   SCIP_CONS**           capconss;     /**< capacity constraints per facility */
   SCIP_CONS**           demandconss;  /**< demand constraints per customer */
   SCIP_Real**           costs;        /**< the transportation costs from a facility to a customer */
   SCIP_Real*            demands;      /**< the customer demands */
   SCIP_Real             capacity;     /**< the capacity of each facility */
   SCIP_Real             fixedcost;    /**< the fixed cost of openning a facility */
   int                   ncustomers;   /**< the number of customers */
   int                   nfacilities;  /**< the number of facilities */
   int                   nfacilityvars;/**< the number of facility variables */
   int*                  ncustomervars;/**< the number of customer variables */
   int                   ncapconss;    /**< the number of capacity constraints */
   int                   ndemandconss; /**< the number of demand constraints */
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
   SCIP_VAR***           customervars,       /**< all variables representing the satisfaction of demand */
   SCIP_CONS**           capconss,           /**< capacity constraints per facility */
   SCIP_CONS**           demandconss,        /**< demand constraints per customer */
   SCIP_Real**           costs,              /**< the transportation costs from a facility to a customer */
   SCIP_Real*            demands,            /**< the customer demands */
   SCIP_Real             capacity,           /**< the capacity of each facility */
   SCIP_Real             fixedcost,          /**< the fixed cost of openning a facility */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nsubproblems,       /**< the number of subproblems */
   int*                  nfacilityvars,      /**< the number of facility variables added to the problem */
   int**                 ncustomervars,      /**< the number of customer variables */
   int*                  ncapconss,          /**< the number of capacity constraints */
   int*                  ndemandconss        /**< the number of demand constraints */
   )
{
   SCIP_CONS* cons;
   SCIP_VAR* var;
   SCIP_VARDATA* vardata;
   int i;
   int j;
   char name[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   (*nfacilityvars) = 0;

   /* adds the capacity constraints to the subproblem */
   for( i = 0; i < nfacilities; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "capacity_%d", i);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, -SCIPinfinity(scip), 0.0) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

      capconss[i] = cons;

      (*ncapconss)++;
   }

   /* adds the demand constraints to the subproblem */
   for( i = 0; i < ncustomers; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "demand_%d", i);
      SCIP_CALL( SCIPcreateConsBasicLinear(scip, &cons, name, 0, NULL, NULL, demands[i], SCIPinfinity(scip)) );

      SCIP_CALL( SCIPaddCons(scip, cons) );

      demandconss[i] = cons;

      (*ndemandconss)++;
   }

   for( i = 0; i < nfacilities; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "facility_%d", i);
      SCIP_CALL( SCIPcreateVarCFLP(scip, &var, name, 0, 1, fixedcost, SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, var) );

      /* storing the variable in the facility variable list */
      facilityvars[i] = var;

      /* creates the variable data */
      SCIP_CALL( SCIPvardataCreateCFLP(scip, &vardata, MASTER, nsubproblems) );

      /* add the variable data to the variable */
      SCIPvarSetData(var, vardata);

      /* adding the variable to the capacity constriants */
      SCIP_CALL( SCIPaddCoefLinear(scip, capconss[i], var, -capacity) );

      /* releases the variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );

      (*nfacilityvars)++;
   }

   /* adding the customer variables to the subproblem */
   for( i = 0; i < ncustomers; i++ )
   {
      for( j = 0; j < nfacilities; j++ )
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "customer(%d,%d)", i, j);
         SCIP_CALL( SCIPcreateVarCFLP(scip, &var, name, 0, SCIPinfinity(scip), costs[j][i], SCIP_VARTYPE_CONTINUOUS) );

         SCIP_CALL( SCIPaddVar(scip, var) );

         /* storing the customer variable in the list */
         customervars[i][j] = var;

         /* creates the variable data */
         SCIP_CALL( SCIPvardataCreateCFLP(scip, &vardata, SUBPROB, 1) );

         /* add the variable data to the variable */
         SCIPvarSetData(var, vardata);

         if( costs[j][i] > 0 )
         {
            /* adding the variable to the capacity constriants */
            SCIP_CALL( SCIPaddCoefLinear(scip, capconss[j], var, 1.0) );

            /* adding the variable to the demand constraints */
            SCIP_CALL( SCIPaddCoefLinear(scip, demandconss[i], var, 1.0) );
         }

         /* releases the variable */
         SCIP_CALL( SCIPreleaseVar(scip, &var) );

         (*ncustomervars)[i]++;
      }
   }

   return SCIP_OKAY;
}

/** creates the Benders' decomposition master problem */
static
SCIP_RETCODE createMasterproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            facilityvars,       /**< all variables representing facilities */
   SCIP_Real             fixedcost,          /**< the fixed cost of openning a facility */
   int                   nfacilities,        /**< the number of facilities */
   int                   nsubproblems,       /**< the number of subproblems */
   int*                  nfacilityvars       /**< the number of facility variables added to the problem */
   )
{
   SCIP_VAR* var;
   SCIP_VARDATA* vardata;
   int i;
   char name[SCIP_MAXSTRLEN];
   assert(scip != NULL);

   (*nfacilityvars) = 0;

   for( i = 0; i < nfacilities; i++ )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "facility_%d", i);
      SCIP_CALL( SCIPcreateVarCFLP(scip, &var, name, 0, 1, fixedcost, SCIP_VARTYPE_BINARY) );

      SCIP_CALL( SCIPaddVar(scip, var) );

      /* storing the variable in the facility variable list */
      facilityvars[i] = var;

      /* creates the variable data */
      SCIP_CALL( SCIPvardataCreateCFLP(scip, &vardata, MASTER, nsubproblems) );

      /* add the variable data to the variable */
      SCIPvarSetData(var, vardata);

      /* releases the variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );

      (*nfacilityvars)++;
   }

   return SCIP_OKAY;
}


/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   SCIP_VAR**            facilityvars,       /**< all variables representing facilities */
   SCIP_VAR***           customervars,       /**< all variables representing the satisfaction of demand */
   SCIP_CONS**           capconss,           /**< capacity constraints per facility */
   SCIP_CONS**           demandconss,        /**< demand constraints per customer */
   SCIP_Real**           costs,              /**< the transportation costs from a facility to a customer */
   SCIP_Real*            demands,            /**< the customer demands */
   SCIP_Real             capacity,           /**< the capacity of each facility */
   SCIP_Real             fixedcost,          /**< the fixed cost of openning a facility */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nfacilityvars,      /**< the number of facility variables */
   int*                  ncustomervars,      /**< the number of customer variables */
   int                   ncapconss,          /**< the number of capacity constraints */
   int                   ndemandconss        /**< the number of demand constraints */
   )
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );

   /* copy variable arrays */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->facilityvars, facilityvars, nfacilities) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->customervars, ncustomers) );
   for( i = 0; i < ncustomers; i++ )
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->customervars[i], customervars[i], nfacilities) );

   /* duplicate arrays */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->capconss, capconss, nfacilities) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->demandconss, demandconss, ncustomers) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->demands, demands, ncustomers) );

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->costs, nfacilities) );
   for( i = 0; i < nfacilities; i++ )
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->costs[i], costs[i], ncustomers) );

   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->ncustomervars, ncustomervars, ncustomers) );

   (*probdata)->capacity = capacity;
   (*probdata)->fixedcost = fixedcost;
   (*probdata)->ncustomers = ncustomers;
   (*probdata)->nfacilities = nfacilities;
   (*probdata)->nfacilityvars = nfacilityvars;
   (*probdata)->ncapconss = ncapconss;
   (*probdata)->ndemandconss = ndemandconss;

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

   assert(scip != NULL);
   assert(probdata != NULL);

   /* release all variables */
   for( i = 0; i < (*probdata)->nfacilityvars; i++ )
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->facilityvars[i]) );

   for( i = 0; i < (*probdata)->nfacilities; i++ )
   {
      for( j = 0; j < (*probdata)->ncustomervars[i]; j++ )
        SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->customervars[j][i]) );
   }

   /* release all constraints */
   for( i = 0; i < (*probdata)->ndemandconss; ++i )
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->demandconss[i]) );

   for( i = 0; i < (*probdata)->ncapconss; ++i )
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->capconss[i]) );

   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->ncustomervars, (*probdata)->nfacilities);

   for( i = (*probdata)->nfacilities - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->costs[i], (*probdata)->ncustomers);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->costs, (*probdata)->nfacilities);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->demands, (*probdata)->ncustomers);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->demandconss, (*probdata)->ncustomers);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->capconss, (*probdata)->nfacilities);

   for( i = (*probdata)->ncustomers - 1; i >= 0; i-- )
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->customervars[i], (*probdata)->nfacilities);

   SCIPfreeBlockMemoryArray(scip, &(*probdata)->customervars, (*probdata)->ncustomers);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->facilityvars, (*probdata)->nfacilities);

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
SCIP_DECL_PROBDELORIG(probdelorigCflp)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransCflp)
{
   int i;

   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->facilityvars, sourcedata->customervars, sourcedata->capconss,
         sourcedata->demandconss, sourcedata->costs, sourcedata->demands, sourcedata->capacity, sourcedata->fixedcost,
         sourcedata->ncustomers, sourcedata->nfacilities, sourcedata->nfacilityvars, sourcedata->ncustomervars,
         sourcedata->ncapconss, sourcedata->ndemandconss) );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->ndemandconss, (*targetdata)->demandconss,
         (*targetdata)->demandconss) );
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->ncapconss, (*targetdata)->capconss,
         (*targetdata)->capconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nfacilityvars, (*targetdata)->facilityvars,
         (*targetdata)->facilityvars) );

   for( i = 0; i < (*targetdata)->ncustomers; i++ )
      SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->ncustomervars[i], (*targetdata)->customervars[i],
            (*targetdata)->customervars[i]) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransCflp)
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
   SCIP_Real*            demands,            /**< the customer demands */
   SCIP_Real             capacity,           /**< the capacity of each facility */
   SCIP_Real             fixedcost,          /**< the fixed cost of openning a facility */
   int                   ncustomers,         /**< the number of customers */
   int                   nfacilities,        /**< the number of facilities */
   int                   nsubproblems,       /**< the number of Benders' decomposition subproblems */
   SCIP_Bool             usebenders          /**< will Benders' decomposition be used to solve the problem */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS** demandconss;
   SCIP_CONS** capconss;
   SCIP_VAR** facilityvars;
   SCIP_VAR*** customervars;
   int nfacilityvars;
   int* ncustomervars;
   int ndemandconss;
   int ncapconss;
   int i;

   assert(scip != NULL);

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigCflp) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransCflp) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransCflp) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   SCIP_CALL( SCIPallocBufferArray(scip, &demandconss, ncustomers) );
   SCIP_CALL( SCIPallocBufferArray(scip, &capconss, nfacilities) );

   SCIP_CALL( SCIPallocBufferArray(scip, &facilityvars, nfacilities) );
   SCIP_CALL( SCIPallocBufferArray(scip, &customervars, ncustomers) );
   for( i = 0; i < ncustomers; i++ )
      SCIP_CALL( SCIPallocBufferArray(scip, &customervars[i], nfacilities) );

   SCIP_CALL( SCIPallocBufferArray(scip, &ncustomervars, ncustomers) );

   nfacilityvars = 0;
   ndemandconss = 0;
   ncapconss = 0;
   for( i = 0; i < ncustomers; i++ )
      ncustomervars[i] = 0;

   if( usebenders )
   {
      /* creating the master problem */
      SCIP_CALL( createMasterproblem(scip, facilityvars, fixedcost, nfacilities, nsubproblems, &nfacilityvars) );
   }
   else
   {
      /* creating the original problem */
      SCIP_CALL( createOriginalproblem(scip, facilityvars, customervars, capconss, demandconss, costs, demands,
            capacity, fixedcost, ncustomers, nfacilities, nsubproblems, &nfacilityvars, &ncustomervars, &ncapconss,
            &ndemandconss) );
   }

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, facilityvars, customervars, capconss, demandconss, costs, demands,
         capacity, fixedcost, ncustomers, nfacilities, nfacilityvars, ncustomervars, ncapconss, ndemandconss) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   /* free local buffer arrays */
   SCIPfreeBufferArray(scip, &ncustomervars);
   for( i = ncustomers - 1; i >= 0; i-- )
      SCIPfreeBufferArray(scip, &customervars[i]);
   SCIPfreeBufferArray(scip, &customervars);
   SCIPfreeBufferArray(scip, &facilityvars);
   SCIPfreeBufferArray(scip, &capconss);
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

/**@} */
