/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_knapsack.h,v 1.15 2004/07/06 17:04:13 bzfpfend Exp $"

/**@file   cons_knapsack.h
 * @brief  constraint handler for knapsack constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __CONS_KNAPSACK_H__
#define __CONS_KNAPSACK_H__


#include "scip.h"


/** creates the handler for knapsack constraints and includes it in SCIP */
extern
RETCODE SCIPincludeConshdlrKnapsack(
   SCIP*            scip                /**< SCIP data structure */
   );

/** creates and captures a knapsack constraint */
extern
RETCODE SCIPcreateConsKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS**           cons,               /**< pointer to hold the created constraint */
   const char*      name,               /**< name of constraint */
   int              nvars,              /**< number of items in the knapsack */
   VAR**            vars,               /**< array with item variables */
   Longint*         weights,            /**< array with item weights */
   Longint          capacity,           /**< capacity of knapsack (right hand side of inequality) */
   Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP? */
   Bool             relax,              /**< should the LP relaxation be separated during LP processing? */
   Bool             separate,           /**< should additional cutting planes be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             removeable          /**< should the constraint be removed from the LP due to aging or cleanup? */
   );

/** solves knapsack problem with dynamic programming;
 *  if needed, one can provide arrays to store all selected items and all not selected items
 */
extern
RETCODE SCIPsolveKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   int              nitems,             /**< number of available items */
   Longint*         weights,            /**< item weights */
   Real*            profits,            /**< item profits */
   Longint          capacity,           /**< capacity of knapsack */
   int*             items,              /**< item numbers, or NULL */
   int*             solitems,           /**< array to store items in solution, or NULL */
   int*             nonsolitems,        /**< array to store items not in solution, or NULL */
   int*             nsolitems,          /**< pointer to store number of items in solution, or NULL */
   int*             nnonsolitems,       /**< pointer to store number of items not in solution, or NULL */
   Real*            solval              /**< pointer to store optimal solution value, or NULL */
   );

/** lifts given cardinality inequality sum(x_i) <= c */
extern
RETCODE SCIPliftKnapsackCardinality(
   SCIP*            scip,               /**< SCIP data structure */
   int*             liftcoefs,          /**< to store lifting coefficient of non-set elements */
   Real*            solvals,            /**< LP solution values of variables */
   Longint*         weights,            /**< item weights */
   Longint          capacity,           /**< capacity of knapsack */
   int*             setvars,            /**< set elements */
   int*             nonsetvars,         /**< non-set elements */
   int              nsetvars,           /**< number of set elements */
   int              nnonsetvars,        /**< number of non-set elements */
   int              maxcardinality,     /**< maximal cardinality of selected subset in given set */ 
   Real*            liftlpval           /**< pointer to store LP solution value of lifted elements */  
   );

/** separates lifted cardinality inequalities for given knapsack problem */
extern
RETCODE SCIPseparateKnapsackCardinality(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint that originates the knapsack problem */
   int              nvars,              /**< number of variables in the knapsack constraint */
   VAR**            vars,               /**< item variables of the knapsack constraint */
   Longint*         weights,            /**< item weights */
   Longint          capacity,           /**< capacity of knapsack */
   int*             ncuts               /**< pointer to add up the number of found cuts */
   );

#endif
