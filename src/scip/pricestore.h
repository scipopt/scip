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
#pragma ident "@(#) $Id: pricestore.h,v 1.1 2003/11/26 16:09:02 bzfpfend Exp $"

/**@file   pricestore.h
 * @brief  methods and datastructures for storing priced variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRICESTORE_H__
#define __PRICESTORE_H__


typedef struct Pricestore PRICESTORE;   /**< storage for priced variables */


#include "def.h"
#include "retcode.h"


/** creates pricing storage */
extern
RETCODE SCIPpricestoreCreate(
   PRICESTORE**     pricestore          /**< pointer to store pricing storage */
   );

/** frees pricing storage */
extern
RETCODE SCIPpricestoreFree(
   PRICESTORE**     pricestore          /**< pointer to store pricing storage */
   );

/** adds variable to pricing storage and capture it */
extern
RETCODE SCIPpricestoreAddVar(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   VAR*             var,                /**< priced variable */
   Real             score,              /**< pricing score of variable (the larger, the better the variable) */
   Bool             root                /**< are we at the root node? */
   );

/** adds variable where zero violates the bounds to pricing storage, capture it */
extern
RETCODE SCIPpricestoreAddBdviolvar(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue,         /**< event queue */
   VAR*             var                 /**< variable, where zero violates the bounds */
   );

/** adds problem variables with negative reduced costs to pricing storage */
extern
RETCODE SCIPpricestoreAddProbVars(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory buffers */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** adds priced variables to the LP */
extern
RETCODE SCIPpricestoreApplyVars(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory buffers */
   SET*             set,                /**< global SCIP settings */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** reset variables' bounds violated by zero to its original value */
extern
RETCODE SCIPpricestoreResetBounds(
   PRICESTORE*      pricestore,         /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics */
   TREE*            tree,               /**< branch-and-bound tree */
   LP*              lp,                 /**< LP data */
   BRANCHCAND*      branchcand,         /**< branching candidate storage */
   EVENTQUEUE*      eventqueue          /**< event queue */
   );

/** gets number of variables in pricing storage */
extern
RETCODE SCIPpricestoreGetNVars(
   PRICESTORE*      pricestore          /**< pricing storage */
   );

/** gets number of variables in pricing storage whose bounds must be reset */
extern
RETCODE SCIPpricestoreGetNBoundResets(
   PRICESTORE*      pricestore          /**< pricing storage */
   );

/** gets time needed to price existing problem variables */
extern
Real SCIPpricestoreGetProbPricingTime(
   PRICESTORE*      pricestore          /**< pricing storage */
   );

/** gets total number of calls to problem variable pricing */
extern
int SCIPpricestoreGetNProbPricings(
   PRICESTORE*      pricestore          /**< pricing storage */
   );

/** gets total number of times, a problem variable was priced in */
extern
int SCIPpricestoreGetNProbvarsFound(
   PRICESTORE*      pricestore          /**< pricing storage */
   );

/** get total number of variables found so far in pricing */
extern
int SCIPpricestoreGetNVarsFound(
   PRICESTORE*      pricestore          /**< pricing storage */
   );

/** get total number of variables priced into the LP so far */
extern
int SCIPpricestoreGetNVarsApplied(
   PRICESTORE*      pricestore          /**< pricing storage */
   );


#endif
