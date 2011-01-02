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
#pragma ident "@(#) $Id: prop_xxx.c,v 1.18 2011/01/02 11:10:44 bzfheinz Exp $"

/**@file   prop_xxx.c
 * @ingroup PROPAGATORS
 * @brief  xxx propagator
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/prop_xxx.h"


#define PROP_NAME              "xxx"
#define PROP_DESC              "propagator template"
#define PROP_PRIORITY                 0
#define PROP_FREQ                    10
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */




/*
 * Data structures
 */

/* TODO: fill in the necessary propagator data */

/** propagator data */
struct SCIP_PropData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of propagator
 */

/* TODO: Implement all necessary propagator methods. The methods with an #if 0 ... #else #define ... are optional */


/** copy method for propagator plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PROPCOPY(propCopyXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define propCopyXxx NULL
#endif

/** destructor of propagator to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PROPFREE(propFreeXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propFreeXxx NULL
#endif


/** initialization method of propagator (called after problem was transformed) */
#if 0
static
SCIP_DECL_PROPINIT(propInitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propInitXxx NULL
#endif


/** deinitialization method of propagator (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PROPEXIT(propExitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propExitXxx NULL
#endif


/** solving process initialization method of propagator (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_PROPINITSOL(propInitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propInitsolXxx NULL
#endif


/** solving process deinitialization method of propagator (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_PROPEXITSOL(propExitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define propExitsolXxx NULL
#endif


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx propagator not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}




/*
 * propagator specific interface methods
 */

/** creates the xxx propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropXxx(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROPDATA* propdata;

   /* create xxx propagator data */
   propdata = NULL;
   /* TODO: (optional) create propagator specific data here */

   /* include propagator */
   SCIP_CALL( SCIPincludeProp(scip, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY,
         propCopyXxx,
         propFreeXxx, propInitXxx, propExitXxx, 
         propInitsolXxx, propExitsolXxx, propExecXxx, propRespropXxx,
         propdata) );

   /* add xxx propagator parameters */
   /* TODO: (optional) add propagator specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
