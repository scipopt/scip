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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pricer.h,v 1.11 2005/01/21 09:17:01 bzfpfend Exp $"

/**@file   pricer.h
 * @brief  internal methods for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRICER_H__
#define __PRICER_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_lp.h"
#include "type_prob.h"
#include "type_scip.h"
#include "type_pricer.h"
#include "pub_pricer.h"



/** creates a variable pricer */
extern
RETCODE SCIPpricerCreate(
   PRICER**         pricer,             /**< pointer to variable pricer data structure */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of variable pricer */
   const char*      desc,               /**< description of variable pricer */
   int              priority,           /**< priority of the variable pricer */
   DECL_PRICERFREE  ((*pricerfree)),    /**< destructor of variable pricer */
   DECL_PRICERINIT  ((*pricerinit)),    /**< initialize variable pricer */
   DECL_PRICEREXIT  ((*pricerexit)),    /**< deinitialize variable pricer */
   DECL_PRICERREDCOST((*pricerredcost)),/**< reduced cost pricing method of variable pricer for feasible LPs */
   DECL_PRICERFARKAS((*pricerfarkas)),  /**< farkas pricing method of variable pricer for infeasible LPs */
   PRICERDATA*      pricerdata          /**< variable pricer data */
   );

/** calls destructor and frees memory of variable pricer */
extern
RETCODE SCIPpricerFree(
   PRICER**         pricer,             /**< pointer to variable pricer data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes variable pricer */
extern
RETCODE SCIPpricerInit(
   PRICER*          pricer,             /**< variable pricer */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** calls exit method of variable pricer */
extern
RETCODE SCIPpricerExit(
   PRICER*          pricer,             /**< variable pricer */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** activates pricer such that it is called in LP solving loop */
extern
RETCODE SCIPpricerActivate(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set                 /**< global SCIP settings */
   );

/** deactivates pricer such that it is no longer called in LP solving loop */
extern
RETCODE SCIPpricerDeactivate(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set                 /**< global SCIP settings */
   );

/** calls reduced cost pricing method of variable pricer */
extern
RETCODE SCIPpricerRedcost(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob                /**< transformed problem */
   );

/** calls farkas pricing method of variable pricer */
extern
RETCODE SCIPpricerFarkas(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob                /**< transformed problem */
   );

/** depending on the LP's solution status, calls reduced cost or farkas pricing method of variable pricer */
extern
RETCODE SCIPpricerExec(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   PROB*            prob,               /**< transformed problem */
   LP*              lp                  /**< LP data */
   );

/** sets priority of variable pricer */
extern
void SCIPpricerSetPriority(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the variable pricer */
   );


#endif
