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
#pragma ident "@(#) $Id: pricer.h,v 1.4 2003/11/26 16:09:01 bzfpfend Exp $"

/**@file   pricer.h
 * @brief  methods and datastructures for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRICER_H__
#define __PRICER_H__


typedef struct Pricer PRICER;           /**< variable pricer data */
typedef struct PricerData PRICERDATA;   /**< locally defined variable pricer data */


/** destructor of variable pricer to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define DECL_PRICERFREE(x) RETCODE x (SCIP* scip, PRICER* pricer)

/** initialization method of variable pricer (called when problem solving starts)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define DECL_PRICERINIT(x) RETCODE x (SCIP* scip, PRICER* pricer)

/** deinitialization method of variable pricer (called when problem solving exits)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define DECL_PRICEREXIT(x) RETCODE x (SCIP* scip, PRICER* pricer)

/** reduced cost pricing method of variable pricer for feasible LPs
 *
 *  Searches for variables that can contribute to improve the current LP's solution value.
 *  In standard branch-and-price, these are variables with negative feasibility, that is negative
 *  reduced costs for non-negative variables, and non-zero reduced costs for variables that can be
 *  negative.
 *
 *  The method is called in the LP solving loop after an LP was proven to be feasible.
 *
 *  Whenever the pricer finds a variable with negative feasibility, it should call SCIPcreateVar()
 *  and SCIPaddVar() to add the variable to the problem. Furthermore, it should call the appropriate
 *  methods of the constraint handlers to add the necessary variable entries to the constraints.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define DECL_PRICERREDCOST(x) RETCODE x (SCIP* scip, PRICER* pricer)

/** farkas pricing method of variable pricer for infeasible LPs
 *
 *  Searches for variables that can contribute to the feasibility of the current LP.
 *  In standard branch-and-price, these are variables with positive farkas values:
 *
 *  The LP was proven infeasible, so we have an infeasibility proof by the dual farkas values y.
 *  The valid inequality  y^T A x >= y^T b  is violated by all x, especially by the (for this
 *  inequality most feasible solution) x' defined by 
 *     x'_i = ub_i, if y^T A_i > 0
 *     x'_i = 0   , if y^T A_i = 0
 *     x'_i = lb_i, if y^T A_i < 0.
 *  Pricing in this case means to add variables i with positive farkas value, i.e. y^T A_i x'_i > 0
 *
 *  The method is called in the LP solving loop after an LP was proven to be infeasible.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define DECL_PRICERFARKAS(x) RETCODE x (SCIP* scip, PRICER* pricer)



#include "scip.h"
#include "def.h"
#include "retcode.h"
#include "set.h"
#include "memory.h"
#include "prob.h"
#include "lp.h"



/** compares two pricers w. r. to their priority */
extern
DECL_SORTPTRCOMP(SCIPpricerComp);

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

/** calls reduced cost pricing method of variable pricer */
extern
RETCODE SCIPpricerRedcost(
   PRICER*          pricer,             /**< variable pricer */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< transformed problem */
   );

/** calls farkas pricing method of variable pricer */
extern
RETCODE SCIPpricerFarkas(
   PRICER*          pricer,             /**< variable pricer */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob                /**< transformed problem */
   );

/** depending on the LP's solution status, calls reduced cost or farkas pricing method of variable pricer */
extern
RETCODE SCIPpricerExec(
   PRICER*          pricer,             /**< variable pricer */
   const SET*       set,                /**< global SCIP settings */
   PROB*            prob,               /**< transformed problem */
   LP*              lp                  /**< LP data */
   );

/** gets user data of variable pricer */
extern
PRICERDATA* SCIPpricerGetData(
   PRICER*          pricer              /**< variable pricer */
   );

/** sets user data of variable pricer; user has to free old data in advance! */
extern
void SCIPpricerSetData(
   PRICER*          pricer,             /**< variable pricer */
   PRICERDATA*      pricerdata          /**< new variable pricer user data */
   );

/** gets name of variable pricer */
extern
const char* SCIPpricerGetName(
   PRICER*          pricer              /**< variable pricer */
   );

/** gets description of variable pricer */
extern
const char* SCIPpricerGetDesc(
   PRICER*          pricer              /**< variable pricer */
   );

/** gets priority of variable pricer */
extern
int SCIPpricerGetPriority(
   PRICER*          pricer              /**< variable pricer */
   );

/** sets priority of variable pricer */
extern
void SCIPpricerSetPriority(
   PRICER*          pricer,             /**< variable pricer */
   SET*             set,                /**< global SCIP settings */
   int              priority            /**< new priority of the variable pricer */
   );

/** gets the number of times, the pricer was called and tried to find a variable with negative reduced costs */
extern
int SCIPpricerGetNCalls(
   PRICER*          pricer              /**< variable pricer */
   );

/** gets the number of variables with negative reduced costs found by this pricer */
extern
int SCIPpricerGetNVarsFound(
   PRICER*          pricer              /**< variable pricer */
   );

/** is variable pricer initialized? */
extern
Bool SCIPpricerIsInitialized(
   PRICER*            pricer                /**< variable pricer */
   );

/** gets time in seconds used in this pricer */
extern
Real SCIPpricerGetTime(
   PRICER*            pricer                /**< variable pricer */
   );


#endif
