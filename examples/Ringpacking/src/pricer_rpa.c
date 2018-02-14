/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pricer_ringpacking.c
 * @brief  Ringpacking variable pricer
 * @author Benjamin Mueller
 *
 * This file implements the variable pricer which check if variables exist with negative reduced cost. See
 * @ref PRICER for more details.
 *
 * @page PRICER Pricing new variables
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "scip/scipdefplugins.h"
#include "pricer_rpa.h"

/**@name Pricer properties
 *
 * @{
 */

#define PRICER_NAME            "ringpacking"
#define PRICER_DESC            "pricer for ringpacking"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE          /* only call pricer if all problem variables have non-negative reduced costs */

/**@} */

#ifndef M_PI
#define M_PI           3.141592653589793238462643
#endif

/*
 * Data structures
 */

/** @brief Variable pricer data used in the \ref pricer_ringpacking.c "pricer" */
struct SCIP_PricerData
{

};


/**@name Local methods
 *
 * @{
 */

/**@} */

/**name Callback methods
 *
 * @{
 */

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeRingpacking)
{
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   SCIPfreeBlockMemoryNull(scip, &pricerdata);

   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitRingpacking)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolRingpacking)
{
   SCIP_PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostRingpacking)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   *result = SCIP_SUCCESS;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return SCIP_OKAY; /*lint !e527*/
}

/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasRingpacking)
{  /*lint --e{715}*/
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** creates the ringpacking variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerRingpacking(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;

   /* create ringpacking variable pricer data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &pricerdata) );
   BMSclearMemory(pricerdata);

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostRingpacking, pricerFarkasRingpacking, pricerdata) );

   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeRingpacking) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitRingpacking) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolRingpacking) );

   /* variable pricer parameters */

   return SCIP_OKAY;
}

/** added problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerRingpackingActivate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;
}

/**@} */
