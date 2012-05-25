/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   pub_pricer.h
 * @ingroup PUBLICMETHODS
 * @brief  public methods for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_PRICER_H__
#define __SCIP_PUB_PRICER_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_pricer.h"

#ifdef __cplusplus
extern "C" {
#endif

/** compares two pricers w. r. to their priority */
extern
SCIP_DECL_SORTPTRCOMP(SCIPpricerComp);

/** comparison method for sorting pricers w.r.t. to their name */
extern
SCIP_DECL_SORTPTRCOMP(SCIPpricerCompName);

/** gets user data of variable pricer */
extern
SCIP_PRICERDATA* SCIPpricerGetData(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** sets user data of variable pricer; user has to free old data in advance! */
extern
void SCIPpricerSetData(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_PRICERDATA*      pricerdata          /**< new variable pricer user data */
   );

/** sets copy callback of pricer */
extern
void SCIPpricerSetCopy(
   SCIP_PRICER*          pricer,             /**< variable pricer */
   SCIP_DECL_PRICERCOPY  ((*pricercopy))      /**< copy callback of pricer */
   );

/** sets destructor callback of pricer */
extern
void SCIPpricerSetFree(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICERFREE  ((*pricerfree))     /**< destructor of pricer */
   );

/** sets initialization callback of pricer */
extern
void SCIPpricerSetInit(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICERINIT ((*pricerinit))     /**< initialize pricer */
   );

/** sets deinitialization callback of pricer */
extern
void SCIPpricerSetExit(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICEREXIT ((*pricerexit))     /**< deinitialize pricer */
   );

/** sets solving process initialization callback of pricer */
extern
void SCIPpricerSetInitsol(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICERINITSOL ((*pricerinitsol))/**< solving process initialization callback of pricer */
   );

/** sets solving process deinitialization callback of pricer */
extern
void SCIPpricerSetExitsol(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICEREXITSOL ((*pricerexitsol))/**< solving process deinitialization callback of pricer */
   );

/** sets Farkas pricing method of variable pricer for infeasible LPs */
extern
void SCIPpricerSetFarkas(
   SCIP_PRICER*          pricer,             /**< pricer */
   SCIP_DECL_PRICERFARKAS((*pricerfarkas))   /**< Farkas pricing method of variable pricer for infeasible LPs */
   );

/** gets name of variable pricer */
extern
const char* SCIPpricerGetName(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** gets description of variable pricer */
extern
const char* SCIPpricerGetDesc(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** gets priority of variable pricer */
extern
int SCIPpricerGetPriority(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** gets the number of times, the pricer was called and tried to find a variable with negative reduced costs */
extern
int SCIPpricerGetNCalls(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** gets the number of variables with negative reduced costs found by this pricer */
extern
int SCIPpricerGetNVarsFound(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** gets time in seconds used in this pricer for setting up for next stages */
extern
SCIP_Real SCIPpricerGetSetupTime(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** gets time in seconds used in this pricer */
extern
SCIP_Real SCIPpricerGetTime(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** returns whether the given pricer is in use in the current problem */
extern
SCIP_Bool SCIPpricerIsActive(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** returns whether the pricer should be delayed until no other pricer finds a new variable */
extern
SCIP_Bool SCIPpricerIsDelayed(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

/** is variable pricer initialized? */
extern
SCIP_Bool SCIPpricerIsInitialized(
   SCIP_PRICER*          pricer              /**< variable pricer */
   );

#ifdef __cplusplus
}
#endif

#endif
