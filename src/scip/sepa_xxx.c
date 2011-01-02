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
#pragma ident "@(#) $Id: sepa_xxx.c,v 1.26 2011/01/02 11:10:42 bzfheinz Exp $"

/**@file   sepa_xxx.c
 * @ingroup SEPARATORS
 * @brief  xxx separator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/sepa_xxx.h"


#define SEPA_NAME              "xxx"
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
SCIP_DECL_SEPACOPY(sepaCopyXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaCopyXxx NULL
#endif

/** destructor of separator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_SEPAFREE(sepaFreeXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaFreeXxx NULL
#endif


/** initialization method of separator (called after problem was transformed) */
#if 0
static
SCIP_DECL_SEPAINIT(sepaInitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitXxx NULL
#endif


/** deinitialization method of separator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_SEPAEXIT(sepaExitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitXxx NULL
#endif


/** solving process initialization method of separator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_SEPAINITSOL(sepaInitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaInitsolXxx NULL
#endif


/** solving process deinitialization method of separator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_SEPAEXITSOL(sepaExitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExitsolXxx NULL
#endif


/** LP solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECLP(sepaExeclpXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExeclpXxx NULL
#endif


/** arbitrary primal solution separation method of separator */
#if 0
static
SCIP_DECL_SEPAEXECSOL(sepaExecsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx separator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define sepaExecsolXxx NULL
#endif




/*
 * separator specific interface methods
 */

/** creates the xxx separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaXxx(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_SEPADATA* sepadata;

   /* create xxx separator data */
   sepadata = NULL;
   /* TODO: (optional) create separator specific data here */

   /* include separator */
   SCIP_CALL( SCIPincludeSepa(scip, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
         SEPA_USESSUBSCIP, SEPA_DELAY,
         sepaCopyXxx, sepaFreeXxx, sepaInitXxx, sepaExitXxx, 
         sepaInitsolXxx, sepaExitsolXxx, sepaExeclpXxx, sepaExecsolXxx,
         sepadata) );

   /* add xxx separator parameters */
   /* TODO: (optional) add separator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
