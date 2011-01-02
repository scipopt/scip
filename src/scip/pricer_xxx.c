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
#pragma ident "@(#) $Id: pricer_xxx.c,v 1.22 2011/01/02 11:10:44 bzfheinz Exp $"

/**@file   pricer_xxx.c
 * @ingroup PRICERS
 * @brief  xxx variable pricer
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "scip/pricer_xxx.h"


#define PRICER_NAME            "xxx"
#define PRICER_DESC            "variable pricer template"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */




/*
 * Data structures
 */

/* TODO: fill in the necessary variable pricer data */

/** variable pricer data */
struct SCIP_PricerData
{
};




/*
 * Local methods
 */

/* put your local methods here, and declare them static */




/*
 * Callback methods of variable pricer
 */

/* TODO: Implement all necessary variable pricer methods. The methods with an #if 0 ... #else #define ... are optional */

/** copy method for pricer plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_PRICERCOPY(pricerCopyXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/
 
   return SCIP_OKAY;
}
#else
#define pricerCopyXxx NULL
#endif

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PRICERFREE(pricerFreeXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerFreeXxx NULL
#endif


/** initialization method of variable pricer (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRICERINIT(pricerInitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerInitXxx NULL
#endif


/** deinitialization method of variable pricer (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRICEREXIT(pricerExitXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerExitXxx NULL
#endif


/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_PRICERINITSOL(pricerInitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerInitsolXxx NULL
#endif


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerExitsolXxx NULL
#endif


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}


#if 0
/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasXxx)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of xxx variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerFarkasXxx NULL
#endif




/*
 * variable pricer specific interface methods
 */

/** creates the xxx variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerXxx(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;

   /* create xxx variable pricer data */
   pricerdata = NULL;
   /* TODO: (optional) create variable pricer specific data here */

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerCopyXxx,
         pricerFreeXxx, pricerInitXxx, pricerExitXxx, 
         pricerInitsolXxx, pricerExitsolXxx, pricerRedcostXxx, pricerFarkasXxx,
         pricerdata) );

   /* add xxx variable pricer parameters */
   /* TODO: (optional) add variable pricer specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
