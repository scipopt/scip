/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2006 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2006 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pricer_healthcare.c,v 1.2 2006/01/03 12:22:39 bzfpfend Exp $"

/**@file   pricer_healthcare.c
 * @brief  healthcare variable pricer
 * @author Arne Nielsen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "pricer_healthcare.h"
#include "probdata_healthcare.h"
#include "scip/cons_setppc.h"


#define PRICER_NAME            "healthcare"
#define PRICER_DESC            "pricer for healthcare tours"
#define PRICER_PRIORITY        0




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

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
#if 0
static
SCIP_DECL_PRICERFREE(pricerFreeHealthcare)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of healthcare variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerFreeHealthcare NULL
#endif


/** initialization method of variable pricer (called after problem was transformed) */
#if 0
static
SCIP_DECL_PRICERINIT(pricerInitHealthcare)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of healthcare variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerInitHealthcare NULL
#endif


/** deinitialization method of variable pricer (called before transformed problem is freed) */
#if 0
static
SCIP_DECL_PRICEREXIT(pricerExitHealthcare)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of healthcare variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerExitHealthcare NULL
#endif


/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
#if 0
static
SCIP_DECL_PRICERINITSOL(pricerInitsolHealthcare)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of healthcare variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerInitsolHealthcare NULL
#endif


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
#if 0
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolHealthcare)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of healthcare variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerExitsolHealthcare NULL
#endif


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostHealthcare)
{  /*lint --e{715}*/
   SCIP_CONS** cons_servejobs;
   SCIP_CONS** cons_workers;
   int njobs;
   int nworkers;
   SCIP_Bool found = FALSE;

   /* get problem data */
   cons_servejobs = HCPgetConsServejobs(scip);
   cons_workers = HCPgetConsWorkers(scip);
   njobs = HCPgetNJobs(scip);
   nworkers = HCPgetNWorkers(scip);

   /* Arne: Hier pricen, um Variablen zu generieren */
   if( found )
   {
      SCIP_VAR* var;
      SCIP_Real obj = 0.0;
      int i;

      SCIP_CALL( SCIPcreateVar(scip, &var, "varname", 0.0, 1.0, obj, SCIP_VARTYPE_BINARY, FALSE, FALSE,
            NULL, NULL, NULL, NULL) );

      /* Arne: hier die Variable in die Constraints eintragen (natuerlich nur in die Constraints, wo die Variable
       * hingehoert!)
       */
      for( i = 0; i < njobs; ++i )
      {
         SCIP_Bool tourbelongstojob = TRUE;
         if( tourbelongstojob )
         {
            SCIP_CALL( SCIPaddCoefSetppc(scip, cons_servejobs[i], var) );
         }
      }
      for( i = 0; i < nworkers; ++i )
      {
         SCIP_Bool tourbelongstoworker = TRUE;
         if( tourbelongstoworker )
         {
            SCIP_CALL( SCIPaddCoefSetppc(scip, cons_workers[i], var) );
         }
      }

      SCIP_CALL( SCIPaddPricedVar(scip, var, 0.0) );
      SCIP_CALL( SCIPreleaseVar(scip, &var) );
   }

   return SCIP_OKAY;
}


#if 0
/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasHealthcare)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of healthcare variable pricer not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define pricerFarkasHealthcare NULL
#endif




/*
 * variable pricer specific interface methods
 */

/** creates the healthcare variable pricer and includes it in SCIP */
SCIP_RETCODE HCPincludePricerHealthcare(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;

   /* create healthcare variable pricer data */
   pricerdata = NULL;
   /* TODO: (optional) create variable pricer specific data here */

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricer(scip, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY,
         pricerFreeHealthcare, pricerInitHealthcare, pricerExitHealthcare, 
         pricerInitsolHealthcare, pricerExitsolHealthcare, pricerRedcostHealthcare, pricerFarkasHealthcare,
         pricerdata) );

   /* add healthcare variable pricer parameters */
   /* TODO: (optional) add variable pricer specific parameters with SCIPaddTypeParam() here */

   return SCIP_OKAY;
}
