/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file   presol_implint.c
 * @ingroup DEFPLUGINS_PRESOL
 * @brief  Presolver that detects implied integer variables
 * @author Rolf van der Hulst
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/presol_implint.h"
#include "scip/scip_presol.h"

#include "scip/pub_message.h"

#define PRESOL_NAME            "implint"
#define PRESOL_DESC            "detects implied integer variables"
#define PRESOL_PRIORITY               0 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers); combined with propagators */
#define PRESOL_MAXROUNDS             -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_TIMING           SCIP_PRESOLTIMING_MEDIUM /* timing of the presolver (fast, medium, or exhaustive) */


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
SCIP_DECL_PRESOLCOPY(presolCopyImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolCopyImplint NULL
#endif


/** destructor of presolver to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PRESOLFREE(presolFreeImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolFreeImplint NULL
#endif


/** initialization method of presolver (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRESOLINIT(presolInitImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitImplint NULL
#endif


/** deinitialization method of presolver (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRESOLEXIT(presolExitImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitImplint NULL
#endif


/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolInitpreImplint NULL
#endif


/** presolving deinitialization method of presolver (called after presolving has been finished) */
#if 0
static
SCIP_DECL_PRESOLEXITPRE(presolExitpreImplint)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of implint presolver not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define presolExitpreImplint NULL
#endif


/** execution method of presolver */

static
SCIP_DECL_PRESOLEXEC(presolExecImplint)
{

   return SCIP_OKAY;
}


/*
 * presolver specific interface methods
 */

/** creates the implint presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludePresolImplint(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_PRESOLDATA* presoldata;
   SCIP_PRESOL* presol;

   /* create implint presolver data */
   presoldata = NULL;
   /* TODO: (optional) create presolver specific data here */

   presol = NULL;

   /* use SCIPincludePresolBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludePresolBasic(scip, &presol, PRESOL_NAME, PRESOL_DESC, PRESOL_PRIORITY, PRESOL_MAXROUNDS,
                                     PRESOL_TIMING, presolExecImplint, presoldata) );

   assert(presol != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPresolCopy(scip, presol, presolCopyImplint) );
   SCIP_CALL( SCIPsetPresolFree(scip, presol, presolFreeImplint) );
   SCIP_CALL( SCIPsetPresolInit(scip, presol, presolInitImplint) );
   SCIP_CALL( SCIPsetPresolExit(scip, presol, presolExitImplint) );
   SCIP_CALL( SCIPsetPresolInitpre(scip, presol, presolInitpreImplint) );
   SCIP_CALL( SCIPsetPresolExitpre(scip, presol, presolExitpreImplint) );

   /* add implint presolver parameters */
   /* TODO: (optional) add presolver specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
