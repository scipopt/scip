/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   price.h
 * @brief  methods and datastructures for pricing variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRICE_H__
#define __PRICE_H__


typedef struct Price PRICE;             /**< storage for priced variables */


#include "def.h"
#include "retcode.h"


extern
RETCODE SCIPpriceCreate(                /**< creates pricing storage */
   PRICE**          price               /**< pointer to store pricing storage */
   );

extern
RETCODE SCIPpriceFree(                  /**< frees pricing storage */
   PRICE**          price               /**< pointer to store pricing storage */
   );

extern
RETCODE SCIPpriceAddVar(                /**< adds variable to pricing storage and capture it */
   PRICE*           price,              /**< pricing storage */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< priced variable */
   Real             score               /**< pricing score of variable (the larger, the better the variable) */
   );

extern
RETCODE SCIPpriceAddBdviolvar(          /**< adds variable where zero violates the bounds to pricing storage, capture it */
   PRICE*           price,              /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   TREE*            tree,               /**< branch-and-bound tree */
   VAR*             var                 /**< variable, where zero violates the bounds */
   );

extern
RETCODE SCIPpriceVars(                  /**< calls all external pricers, prices problem variables, and adds some columns
                                           with negative reduced costs to the LP */
   PRICE*           price,              /**< pricing storage */
   const SET*       set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory buffers */
   STAT*            stat,               /**< dynamic problem statistics */
   PROB*            prob,               /**< transformed problem after presolve */
   LP*              lp,                 /**< LP data */
   TREE*            tree                /**< branch-and-bound tree */
   );

extern
RETCODE SCIPpriceResetBounds(           /**< reset variables' bounds violated by zero to its original value */
   PRICE*           price,              /**< pricing storage */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp,                 /**< LP data */
   TREE*            tree                /**< branch-and-bound tree */
   );



#endif
