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
#pragma ident "@(#) $Id: struct_pricer.h,v 1.1 2003/12/01 14:41:34 bzfpfend Exp $"

/**@file   struct_pricer.h
 * @brief  datastructures for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_PRICER_H__
#define __STRUCT_PRICER_H__


#include "def.h"
#include "type_clock.h"
#include "type_pricer.h"


/** variable pricers data */
struct Pricer
{
   char*            name;               /**< name of variable pricer */
   char*            desc;               /**< description of variable pricer */
   int              priority;           /**< priority of the variable pricer */
   DECL_PRICERFREE  ((*pricerfree));    /**< destructor of variable pricer */
   DECL_PRICERINIT  ((*pricerinit));    /**< initialize variable pricer */
   DECL_PRICEREXIT  ((*pricerexit));    /**< deinitialize variable pricer */
   DECL_PRICERREDCOST((*pricerredcost));/**< reduced cost pricing method of variable pricer for feasible LPs */
   DECL_PRICERFARKAS((*pricerfarkas));  /**< farkas pricing method of variable pricer for infeasible LPs */
   PRICERDATA*      pricerdata;         /**< variable pricers local data */
   CLOCK*           clock;              /**< pricer execution time */
   int              ncalls;             /**< number of times, this pricer was called */
   int              nvarsfound;         /**< number of variables priced in found so far by this pricer */
   unsigned int     initialized:1;      /**< is variable pricer initialized? */
};


#endif
