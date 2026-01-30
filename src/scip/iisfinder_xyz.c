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

/**@file   iisfinder_xyz.c
 * @ingroup DEFPLUGINS_IISFINDER
 * @brief  xyz iis finder
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/iisfinder.h"
#include "scip/iisfinder_xyz.h"


#define IISFINDER_NAME           "xyz"
#define IISFINDER_DESC           "iis finder template"
#define IISFINDER_PRIORITY       0


/*
 * Data structures
 */

/* TODO: fill in the necessary IIS data
 * when this plugin is written in C++, rename SCIP_IISfinderData to something unique (C++ One-Definition-Rule)
 */

/** IIS data */
struct SCIP_IISfinderData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of IIS finders
 */

/* TODO: Implement all necessary IIS finder methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for IIS finder plugin (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_IISFINDERCOPY(iisfinderCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS finder not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define iisfinderCopyXyz NULL
#endif


/** destructor of IIS finder to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_IISFINDERFREE(iisfinderFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS finder not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define iisfinderFreeXyz NULL
#endif

/** execution method of IIS */
static
SCIP_DECL_IISFINDEREXEC(iisfinderExecXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS finder not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/*
 * IIS finder specific interface methods
 */

/** creates the xyz IIS finder and includes it in SCIP */
SCIP_RETCODE SCIPincludeIISfinderXyz(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_IISFINDERDATA* iisfinderdata = NULL;
   SCIP_IISFINDER* iisfinder = NULL;

   /* create xyz IIS finder data */

   /* TODO: (optional) create IIS specific data here */

   /* include IIS */
#if 0
   /* use SCIPincludeIISfinder() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeIISfinder(scip, IISFINDER_NAME, IISFINDER_DESC, IISFINDER_PRIORITY,
         iisfinderCopyXyz, iisfinderFreeXyz, iisfinderExecXyz,
         iisdata) );
#else
   /* use SCIPincludeIISfinderBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independently of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeIISfinderBasic(scip, &iisfinder, IISFINDER_NAME, IISFINDER_DESC, IISFINDER_PRIORITY,
         iisfinderExecXyz, iisfinderdata) );

   assert(iisfinder != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetIISfinderCopy(scip, iisfinder, iisfinderCopyXyz) );
   SCIP_CALL( SCIPsetIISfinderFree(scip, iisfinder, iisfinderFreeXyz) );
#endif

   /* add xyz IIS parameters */
   /* TODO: (optional) add IIS finder specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
