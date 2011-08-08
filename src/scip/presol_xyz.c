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

/**@file   presol_xyz.c
 * @ingroup PRESOLVERS
 * @brief  xyz presolver
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_xyz.h"


#define PRESOL_NAME            "xyz"
#define PRESOL_DESC            "presolver template"
#define PRESOL_PRIORITY               0 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY              FALSE /**< should presolver be delayed, if other presolvers found reductions? */




/*
 * Data structures
 */

/* TODO: fill in the necessary presolver data */

/** presolver data */
struct SCIP_PresolData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of presolver
 */

/* TODO: Implement all necessary presolver methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for constraint handler plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRESOLCOPY(presolCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define presolCopyXyz NULL
#endif


/** destructor of presolver to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PRESOLFREE(presolFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolFreeXyz NULL
#endif


/** initialization method of presolver (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRESOLINIT(presolInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitXyz NULL
#endif


/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitXyz NULL
#endif


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreXyz NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#if 0
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitpreXyz NULL
#endif


/** execution method of presolver */
static
SCIP_DECL_PRESOLEXEC(presolExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}





/*
 * presolver specific interface methods
 */

/** creates the xyz presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRESOLDATA* presoldata;

   /* create xyz presolver data */
   presoldata = NULL;
   /* TODO: (optional) create presolver specific data here */

   /* include presolver */
   SCIP_CALL( SCIPincludePresol(scip, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS, PRESOL_DELAY,
         presolCopyXyz,
         presolFreeXyz, presolInitXyz, presolExitXyz, 
         presolInitpreXyz, presolExitpreXyz, presolExecXyz,
         presoldata) );

   /* add xyz presolver parameters */
   /* TODO: (optional) add presolver specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
