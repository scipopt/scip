/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pricer_xxx.c,v 1.1 2003/11/27 17:48:45 bzfpfend Exp $"

/**@file   pricer_xxx.c
 * @brief  xxx variable pricer
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "pricer_xxx.h"


#define PRICER_NAME            "xxx"
#define PRICER_DESC            "variable pricer template"
#define PRICER_PRIORITY        0




/*
 * Data structures
 */

/* TODO: fill in the necessary variable pricer data */

/** variable pricer data */
struct PricerData
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

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
#if 0
static
DECL_PRICERFREE(pricerFreeXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx variable pricer not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerFreeXxx NULL
#endif


/** initialization method of variable pricer (called when problem solving starts) */
#if 0
static
DECL_PRICERINIT(pricerInitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx variable pricer not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerInitXxx NULL
#endif


/** deinitialization method of variable pricer (called when problem solving exits) */
#if 0
static
DECL_PRICEREXIT(pricerExitXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx variable pricer not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerExitXxx NULL
#endif


/** reduced cost pricing method of variable pricer for feasible LPs */
static
DECL_PRICERREDCOST(pricerRedcostXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx variable pricer not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}


#if 0
/** farkas pricing method of variable pricer for infeasible LPs */
static
DECL_PRICERREDCOST(pricerFarkasXxx)
{  /*lint --e{715}*/
   errorMessage("method of xxx variable pricer not implemented yet\n");
   abort(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerFarkasXxx NULL
#endif




/*
 * variable pricer specific interface methods
 */

/** creates the xxx variable pricer and includes it in SCIP */
RETCODE SCIPincludePricerXxx(
   SCIP*            scip                /**< SCIP data structure */
   )
{
   PRICERDATA* pricerdata;

   /* create xxx variable pricer data */
   pricerdata = NULL;
   /* TODO: (optional) create variable pricer specific data here */

   /* include variable pricer */
   CHECK_OKAY( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY,
                  pricerFreeXxx, pricerInitXxx, pricerExitXxx, pricerRedcostXxx, pricerFarkasXxx,
                  pricerdata) );

   /* add xxx variable pricer parameters */
   /* TODO: (optional) add variable pricer specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
