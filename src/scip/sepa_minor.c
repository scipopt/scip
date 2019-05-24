/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_minor.c
 * @ingroup DEFPLUGINS_SEPA
 * @brief  minor separator
 * @author Benjamin Mueller
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_minor.h"


#define SEPA_NAME              "minor"
#define SEPA_DESC              "separator template"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */


/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

/** separator data */
struct SCIP_SepaData
{
};


/*
 * Local methods
 */

/* put your local methods here, and declare them static */


/*
 * Callback methods of separator
 */

/* TODO: Implement all necessary separator methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for separator plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_SEPACOPY(sepaCopyMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaCopyMinor NULL
#endif

/** destructor of separator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_SEPAFREE(sepaFreeMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaFreeMinor NULL
#endif


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitMinor NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitMinor NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolMinor NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitsolMinor NULL
#endif


/** LP solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECLP(sepaExeclpMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExeclpMinor NULL
#endif


/** arbitrary primal solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolMinor)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of minor separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExecsolMinor NULL
#endif


/*
 * separator specific interface methods
 */

/** creates the minor separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaMinor(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create minor separator data */
   sepadata = NULL;
   /* TODO: (optional) create separator specific data here */

   sepa = NULL;

   /* include separator */
#if 0
   /* use SCIPincludeSepa() if you want to set all callbacks explicitly and realize (by getting compiler errors) when
    * new callbacks are added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyMinor, sepaFreeMinor, sepaInitMinor, sepaExitMinor, sepaInitsolMinor, sepaExitsolMinor, sepaExeclpMinor, sepaExecsolMinor,
         sepadata) );
#else
   /* use SCIPincludeSepaBasic() plus setter functions if you want to set callbacks one-by-one and your code should
    * compile independent of new callbacks being added in future SCIP versions
    */
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaExeclpMinor, sepaExecsolMinor,
         sepadata) );

   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaCopy(scip, sepa, sepaCopyMinor) );
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeMinor) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitMinor) );
   SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitMinor) );
   SCIP_CALL( SCIPsetSepaInitsol(scip, sepa, sepaInitsolMinor) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolMinor) );
#endif

   /* add minor separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
