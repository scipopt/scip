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
#pragma ident "@(#) $Id: struct_pricestore.h,v 1.1 2003/12/01 14:41:34 bzfpfend Exp $"

/**@file   struct_pricestore.h
 * @brief  datastructures for storing priced variables
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_PRICESTORE_H__
#define __STRUCT_PRICESTORE_H__


#include "def.h"
#include "type_clock.h"
#include "type_var.h"
#include "type_pricestore.h"


/** storage for priced variables */
struct Pricestore
{
   CLOCK*           probpricingtime;    /**< time needed to price existing problem variables */
   VAR**            vars;               /**< array with priced variables with violated reduced costs sorted by score */
   Real*            scores;             /**< score for each priced variable (e.g. |redcost|/#nonzeros) */
   VAR**            bdviolvars;         /**< variables where zero violates the bounds */
   Real*            bdviolvarslb;       /**< lower bounds of bdviolvars */
   Real*            bdviolvarsub;       /**< upper bounds of bdbiolvars */
   int              varssize;           /**< size of vars and score arrays */
   int              nvars;              /**< number of priced variables (max. is set->maxpricevars) */
   int              bdviolvarssize;     /**< size of bdviolvars, bdviolvarslb, and bdviolvarsub arrays */
   int              nbdviolvars;        /**< number of variables, where zero violates the bounds */
   int              naddedbdviolvars;   /**< number of bound violated variables already added to the LP */
   int              nprobpricings;      /**< total number of calls to problem variable pricing */
   int              nprobvarsfound;     /**< total number of problem variables, that were added (and evtl. thrown away) */
   int              nvarsfound;         /**< total number of variables, that were added (and evtl. thrown away) */
   int              nvarsapplied;       /**< total number of variables, that were added to the LP */
};


#endif
