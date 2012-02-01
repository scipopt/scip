/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2011 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*#define SCIP_DEBUG*/

/**@file   presol_components.c
 * @brief  components presolver
 * @author Dieter Weninger
 * @author Gerald Gamrath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_components.h"


#define PRESOL_NAME            "components"
#define PRESOL_DESC            "components presolver"
#define PRESOL_PRIORITY        -9200000 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY               TRUE /**< should presolver be delayed, if other presolvers found reductions? */


#define DEFAULT_SEARCH             FALSE /**< should be searched for components? */
#define DEFAULT_WRITEPROBLEMS      FALSE /**< should the single components be written as an .lp-file? */
#define DEFAULT_MAXINTVARS            20 /**< maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
#define DEFAULT_NODELIMIT          10000 /**< maximum number of nodes to be solved in subproblems */
#define DEFAULT_INTFACTOR              1 /**< the weight of an integer variable compared to binary variables */

/*
 * Data structures
 */

/** control parameters */
struct SCIP_PresolData
{
   SCIP_Bool             dosearch;           /** should be searched for components? */
   SCIP_Bool             writeproblems;      /** should the single components be written as an .lp-file? */
   int                   maxintvars;         /** maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving) */
   SCIP_Longint          nodelimit;          /** maximum number of nodes to be solved in subproblems */
   SCIP_Real             intfactor;          /** the weight of an integer variable compared to binary variables */

   SCIP**                components;         /** sub-SCIPs storing the components */
   SCIP_HASHMAP**        varmaps;            /** hashmaps mapping from original variables to variables in the sub-SCIPs */
   SCIP_HASHMAP*         consmap;            /** hashmaps mapping from original constraints to constraints in the sub-SCIPs
                                              *  (needed only for performance reasons)
                                              */
   int                   componentssize;     /** size of arrays components and varmaps */
   int                   ncomponents;        /** number of components */
};


/*
 * Local methods
 */

/** initializes presolver data */
static
void initPresoldata(
   SCIP_PRESOLDATA*      presoldata          /**< presolver data */
   )
{
   assert(presoldata != NULL);

   presoldata->dosearch = 0;
   presoldata->components = NULL;
   presoldata->varmaps = NULL;
   presoldata->consmap = NULL;
   presoldata->componentssize = 0;
   presoldata->ncomponents = 0;
}


/** performs presolving by searching for components */
static
SCIP_RETCODE presolComponents(
   SCIP*                 scip,               /**< SCIP main data structure */
   SCIP_PRESOL*          presol,             /**< the presolver itself */
   int*                  nfixedvars,         /**< pointer to total number of variables fixed of all presolvers */
   int*                  ndelconss,          /**< pointer to total number of deleted constraints of all presolvers */
   SCIP_RESULT*          result              /**< pointer to store the result of the presolving call */
   )
{
   SCIP_PRESOLDATA* presoldata;

   assert(result != NULL);
   *result = SCIP_DIDNOTRUN;

   if( SCIPgetStage(scip) != SCIP_STAGE_PRESOLVING || SCIPinProbing(scip) )
      return SCIP_OKAY;

   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);
   if( !presoldata->dosearch )
   {
      /* do not search for components  */
      return SCIP_OKAY;
   }

   *result = SCIP_DIDNOTFIND;


   /* TODO: do work here ... */


   return SCIP_OKAY;
}


/*
 * Callback methods of presolver
 */

/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRESOLCOPY(presolCopyComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolCopyComponents NULL
#endif


/** destructor of presolver to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRESOLFREE(presolFreeComponents)
{  /*lint --e{715}*/
   SCIP_PRESOLDATA* presoldata;

   /* free presolver data */
   presoldata = SCIPpresolGetData(presol);
   assert(presoldata != NULL);

   SCIPfreeMemory(scip, &presoldata);
   SCIPpresolSetData(presol, NULL);

   return SCIP_OKAY;
}


/** initialization method of presolver (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRESOLINIT(presolInitComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitComponents NULL
#endif


/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitComponents NULL
#endif


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreComponents)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of components presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreComponents NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreComponents)
{  /*lint --e{715}*/

   SCIP_CALL( presolComponents(scip, presol, NULL, NULL, result) );

   *result = SCIP_FEASIBLE;

   return SCIP_OKAY;
}


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecComponents)
{  /*lint --e{715}*/

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the components presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolComponents(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;

   /* create components presolver data */
   SCIP_CALL( SCIPallocMemory(scip, &presoldata) );
   initPresoldata(presoldata);

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip,
         PRESOL_NAME,
         PRESOL_DESC,
         PRESOL_PRIORITY,
         PRESOL_MAXROUNDS,
         PRESOL_DELAY,
         presolCopyComponents,
         presolFreeComponents,
         presolInitComponents,
         presolExitComponents,
         presolInitpreComponents,
         presolExitpreComponents,
         presolExecComponents,
         presoldata) );

   /* add presolver parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/components/dosearch",
         "search for components (0: no search, 1: do search)",
         &presoldata->dosearch, FALSE, DEFAULT_SEARCH, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip,
         "presolving/components/writeproblems",
         "should the single components be written as an .lp-file?",
         &presoldata->writeproblems, FALSE, DEFAULT_WRITEPROBLEMS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip,
         "presolving/components/maxintvars",
         "maximum number of integer (or binary) variables to solve a subproblem directly (-1: no solving)",
         &presoldata->maxintvars, FALSE, DEFAULT_MAXINTVARS, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddLongintParam(scip,
         "presolving/components/nodelimit",
         "maximum number of nodes to be solved in subproblems",
         &presoldata->nodelimit, FALSE, DEFAULT_NODELIMIT, -1, SCIP_LONGINT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip,
         "presolving/components/intfactor",
         "the weight of an integer variable compared to binary variables",
         &presoldata->intfactor, FALSE, DEFAULT_INTFACTOR, 0, SCIP_REAL_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
