/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License.             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: objpricer.cpp,v 1.12 2005/05/31 17:20:09 bzfpfend Exp $"

/**@file   objpricer.cpp
 * @brief  C++ wrapper for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <cassert>

#include "objpricer.h"




/*
 * Data structures
 */

/** variable pricer data */
struct PricerData
{
   scip::ObjPricer* objpricer;          /**< variable pricer object */
   Bool             deleteobject;       /**< should the pricer object be deleted when pricer is freed? */
};




/*
 * Callback methods of variable pricer
 */

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
DECL_PRICERFREE(pricerFreeObj)
{  /*lint --e{715}*/
   PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   CHECK_OKAY( pricerdata->objpricer->scip_free(scip, pricer) );

   /* free pricer object */
   if( pricerdata->deleteobject )
      delete pricerdata->objpricer;

   /* free pricer data */
   delete pricerdata;
   SCIPpricerSetData(pricer, NULL);
   
   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
DECL_PRICERINIT(pricerInitObj)
{  /*lint --e{715}*/
   PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   CHECK_OKAY( pricerdata->objpricer->scip_init(scip, pricer) );

   return SCIP_OKAY;
}


/** deinitialization method of variable pricer (called before transformed problem is freed) */
static
DECL_PRICEREXIT(pricerExitObj)
{  /*lint --e{715}*/
   PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   CHECK_OKAY( pricerdata->objpricer->scip_exit(scip, pricer) );

   return SCIP_OKAY;
}


/** solving process initialization method of variable pricer (called when branch and bound process is about to begin) */
static
DECL_PRICERINITSOL(pricerInitsolObj)
{  /*lint --e{715}*/
   PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   CHECK_OKAY( pricerdata->objpricer->scip_initsol(scip, pricer) );

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
DECL_PRICEREXITSOL(pricerExitsolObj)
{  /*lint --e{715}*/
   PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   CHECK_OKAY( pricerdata->objpricer->scip_exitsol(scip, pricer) );

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
static
DECL_PRICERREDCOST(pricerRedcostObj)
{  /*lint --e{715}*/
   PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   CHECK_OKAY( pricerdata->objpricer->scip_redcost(scip, pricer) );

   return SCIP_OKAY;
}


/** farkas pricing method of variable pricer for infeasible LPs */
static
DECL_PRICERREDCOST(pricerFarkasObj)
{  /*lint --e{715}*/
   PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   assert(pricerdata->objpricer != NULL);

   /* call virtual method of pricer object */
   CHECK_OKAY( pricerdata->objpricer->scip_farkas(scip, pricer) );

   return SCIP_OKAY;
}




/*
 * variable pricer specific interface methods
 */

/** creates the variable pricer for the given variable pricer object and includes it in SCIP */
RETCODE SCIPincludeObjPricer(
   SCIP*            scip,               /**< SCIP data structure */
   scip::ObjPricer* objpricer,          /**< variable pricer object */
   Bool             deleteobject        /**< should the pricer object be deleted when pricer is freed? */
   )
{
   PRICERDATA* pricerdata;

   /* create variable pricer data */
   pricerdata = new PRICERDATA;
   pricerdata->objpricer = objpricer;
   pricerdata->deleteobject = deleteobject;

   /* include variable pricer */
   CHECK_OKAY( SCIPincludePricer(scip, objpricer->scip_name_, objpricer->scip_desc_, objpricer->scip_priority_,
         pricerFreeObj, pricerInitObj, pricerExitObj, 
         pricerInitsolObj, pricerExitsolObj, pricerRedcostObj, pricerFarkasObj,
         pricerdata) );

   return SCIP_OKAY;
}

/** returns the variable pricer object of the given name, or NULL if not existing */
scip::ObjPricer* SCIPfindObjPricer(
   SCIP*            scip,               /**< SCIP data structure */
   const char*      name                /**< name of variable pricer */
   )
{
   PRICER* pricer;
   PRICERDATA* pricerdata;

   pricer = SCIPfindPricer(scip, name);
   if( pricer == NULL )
      return NULL;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->objpricer;
}
   
/** returns the variable pricer object for the given pricer */
scip::ObjPricer* SCIPgetObjPricer(
   SCIP*            scip,               /**< SCIP data structure */
   PRICER*          pricer              /**< pricer */
   )
{
   PRICERDATA* pricerdata;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   return pricerdata->objpricer;
}
