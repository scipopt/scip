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

/**@file   iis_xyz.c
 * @ingroup DEFPLUGINS_IIS
 * @brief  xyz iis
 * @author Mark Turner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/iis_xyz.h"


#define IIS_NAME                 "xyz"
#define IIS_DESC                 "iis template"
#define IIS_PRIORITY             0


/*
 * Data structures
 */

/* TODO: fill in the necessary IIS data */

/** IIS data */
struct SCIP_IISData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of IIS
 */

/* TODO: Implement all necessary IIS methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for IIS plugin (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_IISCOPY(iisCopyXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define iisCopyXyz NULL
#endif


/** destructor of IIS to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_IISFREE(iisFreeXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define iisFreeXyz NULL
#endif


/** initialization method of IIS (called after problem was transformed) */
#if 0
static
SCIP_DECL_IISINIT(iisInitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define iisInitXyz NULL
#endif


/** deinitialization method of IIS (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_IISEXIT(iisExitXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define iisExitXyz NULL
#endif


/** solving process initialization method of IIS (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_IISINITSOL(iisInitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define iisInitsolXyz NULL
#endif


/** solving process deinitialization method of IIS (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_IISEXITSOL(iisExitsolXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define iisExitsolXyz NULL
#endif


/** execution method of IIS */
static
SCIP_DECL_IISGENERATE(iisGenerateXyz)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xyz IIS not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
   
   return SCIP_OKAY;
}


/*
 * IIS specific interface methods
 */

/** creates the xyz IIS and includes it in SCIP */
SCIP_RETCODE SCIPincludeIISXyz(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_IISDATA* iisdata = NULL;
   SCIP_IIS* iis = NULL;
   
   /* create xyz IIS data */
   
   /* TODO: (optional) create IIS specific data here */
   
   /* include IIS */
#if 0
   /* use SCIPincludeIIS() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeIIS(scip, IIS_NAME, IIS_DESC, IIS_PRIORITY,
         iisCopyXyz, iisFreeXyz, iisInitXyz, iisExitXyz, iisInitsolXyz, iisExitsolXyz, iisGenerateXyz,
         iisdata) );
#else
   /* use SCIPincludeIISBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independently of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeIISBasic(scip, &iis, IIS_NAME, IIS_DESC, IIS_PRIORITY, iisGenerateXyz,
                                     iisdata) );
   
   assert(iis != NULL);
   
   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetIISCopy(scip, iis, iisCopyXyz) );
   SCIP_CALL( SCIPsetIISFree(scip, iis, iisFreeXyz) );
#endif
   
   /* add xyz IIS parameters */
   /* TODO: (optional) add IIS specific parameters with SCIPaddTypeParam() here */
   
   return SCIP_OKAY;
}