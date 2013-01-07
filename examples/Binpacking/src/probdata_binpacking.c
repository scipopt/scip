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

/**@file   probdata_binpacking.c
 * @brief  Problem data for binpacking problem
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This file handles the main problem data used in that project. For more details see \ref PROBLEMDATA page.
 *
 * @page PROBLEMDATA Main problem data
 *
 * The problem data is accessible in all plugins. The function SCIPgetProbData() returns the pointer to that
 * structure. We use this data structure to store all the information of the binpacking problem. Since this structure is
 * not visible in the other plugins, we implemented setter and getter functions to access this data. The problem data
 * structure SCIP_ProbData is shown below.
 *
 * \code
 *  ** @brief Problem data which is accessible in all places
 *  *
 *  *   This problem data is used to store the input of the binpacking instance, all variables which are created, and all
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
 * for each item of the binpacking problem one set covering constraint and creates an initial set of variables for the
 * packings. Note that the set covering constraints have to have the <code>modifiable</code>-flag set to TRUE. This is
 * necessary to tell the solver that these constraints are not completed yet. This means, during the search new
 * variables/packings might be added.  The solver needs this information because certain reductions are not allowed.
 * See the body of the function SCIPprobdataCreate() for more details.
 *
 * A list of all interface methods can be found in probdata_binpacking.h.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "probdata_binpacking.h"
#include "vardata_binpacking.h"
#include "pricer_binpacking.h"

#include "scip/cons_setppc.h"
#include "scip/scip.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the binpacking, all variables which are created, and all
 * constrsaints.
 */
struct SCIP_ProbData
{
   SCIP_VAR**            vars;         /**< all exiting variables in the problem */
   SCIP_CONS**           conss;        /**< set partitioning constraints for each item exactly one */
   SCIP_Longint*         weights;      /**< array of item weights */
   int*                  ids;          /**< array of item ids */
   int                   nvars;        /**< number of generated variables */
   int                   varssize;     /**< size of the variable array */
   int                   nitems;       /**< number of items */
   SCIP_Longint          capacity;     /**< bin capacity */
};


/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "addedvar"
#define EVENTHDLR_DESC         "event handler for catching added variables"

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecAddedVar)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   SCIPdebugMessage("exec method of event handler for added variable to probdata\n");

   /* add new variable to probdata */
   SCIP_CALL( SCIPprobdataAddVar(scip, SCIPgetProbData(scip), SCIPeventGetVar(event)) );

   return SCIP_OKAY;
}

/**@} */


/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   SCIP_VAR**            vars,               /**< all exist variables */
   SCIP_CONS**           conss,              /**< set partitioning constraints for each job exactly one */
   SCIP_Longint*         weights,            /**< array containing the item weights */
   int*                  ids,                /**< array of item ids */
   int                   nvars,              /**< number of variables */
   int                   nitems,             /**< number of items */
   SCIP_Longint          capacity            /**< bin capacity */
   )
{
   assert(scip != NULL);
   assert(probdata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocMemory(scip, probdata) );

   if( nvars > 0 )
   {
      /* copy variable array */
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*probdata)->vars, vars, nvars) );
   }
   else
      (*probdata)->vars = NULL;

   /* duplicate arrays */
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*probdata)->conss, conss, nitems) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*probdata)->weights, weights, nitems) );
   SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*probdata)->ids, ids, nitems) );

   (*probdata)->nvars = nvars;
   (*probdata)->varssize = nvars;
   (*probdata)->nitems = nitems;
   (*probdata)->capacity = capacity;

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

   /* release all variables */
   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }

   /* release all constraints */
   for( i = 0; i < (*probdata)->nitems; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->conss[i]) );
   }

   /* free memory of arrays */
   SCIPfreeMemoryArray(scip, &(*probdata)->vars);
   SCIPfreeMemoryArray(scip, &(*probdata)->conss);
   SCIPfreeMemoryArray(scip, &(*probdata)->weights);
   SCIPfreeMemoryArray(scip, &(*probdata)->ids);

   /* free probdata */
   SCIPfreeMemory(scip, probdata);

   return SCIP_OKAY;
}

/** create initial columns */
static
SCIP_RETCODE createInitialColumns(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   SCIP_CONS** conss;
   SCIP_VARDATA* vardata;
   SCIP_VAR* var;
   char name[SCIP_MAXSTRLEN];

   int* ids;
   SCIP_Longint* weights;
   int nitems;

   int i;

   conss = probdata->conss;
   ids = probdata->ids;
   weights = probdata->weights;
   nitems = probdata->nitems;

   /* create start solution each item in exactly one bin */
   for( i = 0; i < nitems; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "item_%d", ids[i]);

      SCIPdebugMessage("create variable for item %d with weight = %"SCIP_LONGINT_FORMAT"\n", ids[i], weights[i]);

      /* create variable for the packing pattern which contains only this item */
      SCIP_CALL( SCIPcreateVarBinpacking(scip, &var, name, 1.0, TRUE, TRUE, NULL) );

      /* add variable to the problem */
      SCIP_CALL( SCIPaddVar(scip, var) );

      /* store variable in the problme data */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, var) );

      /* add variable to corresponding set covering constraint */
      SCIP_CALL( SCIPaddCoefSetppc(scip, conss[i], var) );

      /* create the variable data for the variable; the variable data contains the information in which constraints the
       * variable appears */
      SCIP_CALL( SCIPvardataCreateBinpacking(scip, &vardata, &i, 1) );

      /* add the variable data to the variable */
      SCIPvarSetData(var, vardata);

      /* change the upper bound of the binary variable to lazy since the upper bound is already enforced
       * due to the objective function the set covering constraint;
       * The reason for doing is that, is to avoid the bound of x <= 1 in the LP relaxation since this bound
       * constraint would produce a dual variable which might have a positive reduced cost
       */
      SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   return SCIP_OKAY;
}

/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigBinpacking)
{
   SCIPdebugMessage("free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransBinpacking)
{
   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->vars, sourcedata->conss, sourcedata->weights, sourcedata->ids,
         sourcedata->nvars, sourcedata->nitems, sourcedata->capacity) );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nitems, (*targetdata)->conss, (*targetdata)->conss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nvars, (*targetdata)->vars, (*targetdata)->vars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransBinpacking)
{
   SCIPdebugMessage("free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolBinpacking)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* catch variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolBinpacking)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* drop variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );


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
   int*                  ids,                /**< array of item ids */
   SCIP_Longint*         weights,            /**< array containing the item weights */
   int                   nitems,             /**< number of items */
   SCIP_Longint          capacity            /**< bin capacity */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS** conss;
   char name[SCIP_MAXSTRLEN];
   int i;

   assert(scip != NULL);

   /* create event handler if it does not exist yet */
   if( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL )
   {
      SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL) );
   }

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigBinpacking) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransBinpacking) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransBinpacking) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolBinpacking) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolBinpacking) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* tell SCIP that the objective will be always integral */
   SCIP_CALL( SCIPsetObjIntegral(scip) );

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, nitems) );

   /* create set covering constraints for each item */
   for( i = 0; i < nitems; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "item_%d", ids[i]);

      SCIP_CALL( SCIPcreateConsBasicSetcover(scip, &conss[i], name, 0, NULL) );

      /* declare constraint modifiable for adding variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, conss[i], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, conss[i]) );   
   }   
   
   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, NULL, conss, weights, ids, 0, nitems, capacity) );

   SCIP_CALL( createInitialColumns(scip, probdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   SCIP_CALL( SCIPpricerBinpackingActivate(scip, conss, weights, ids, nitems, capacity) );

   /* free local buffer arrays */
   SCIPfreeBufferArray(scip, &conss);

   return SCIP_OKAY;
}

/** returns array of item ids */
int* SCIPprobdataGetIds(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->ids;
}

/** returns array of item weights */
SCIP_Longint* SCIPprobdataGetWeights(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->weights;
}

/** returns number of items */
int SCIPprobdataGetNItems(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->nitems;
}

/** returns bin capacity */
SCIP_Longint SCIPprobdataGetCapacity(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->capacity;
}

/** returns array of all variables itemed in the way they got generated */
SCIP_VAR** SCIPprobdataGetVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->vars;
}

/** returns number of variables */
int SCIPprobdataGetNVars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->nvars;
}

/** returns array of set partitioning constrains */
SCIP_CONS** SCIPprobdataGetConss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   return probdata->conss;
}

/** adds given variable to the problem data */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var                 /**< variables to add */
   )
{
   /* check if enough memory is left */
   if( probdata->varssize == probdata->nvars )
   {
      probdata->varssize = MAX(100, probdata->varssize * 2);
      SCIP_CALL( SCIPreallocMemoryArray(scip, &probdata->vars, probdata->varssize) );
   }

   /* caputure variables */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   probdata->vars[probdata->nvars] = var;
   probdata->nvars++;

   SCIPdebugMessage("added variable to probdata; nvars = %d\n", probdata->nvars);

   return SCIP_OKAY;
}

/**@} */
