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
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: cons_knapsack.h,v 1.25 2005/08/04 10:38:56 bzfwolte Exp $"

/**@file   cons_knapsack.h
 * @brief  constraint handler for knapsack constraints
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_KNAPSACK_H__
#define __SCIP_CONS_KNAPSACK_H__


#include "scip/scip.h"


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
   Bool             separate,           /**< should the constraint be separated during LP processing? */
   Bool             enforce,            /**< should the constraint be enforced during node processing? */
   Bool             check,              /**< should the constraint be checked for feasibility? */
   Bool             propagate,          /**< should the constraint be propagated during node processing? */
   Bool             local,              /**< is constraint only valid locally? */
   Bool             modifiable,         /**< is constraint modifiable (subject to column generation)? */
   Bool             dynamic,            /**< is constraint subject to aging? */
   Bool             removeable          /**< should the relaxation be removed from the LP due to aging or cleanup? */
   );

/** adds new item to knapsack constraint */
extern
RETCODE SCIPaddCoefKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint data */
   VAR*             var,                /**< item variable */
   Longint          weight              /**< item weight */
   );

/** gets the dual solution of the knapsack constraint in the current LP */
extern
Real SCIPgetDualsolKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
   );

/** gets the dual farkas value of the knapsack constraint in the current infeasible LP */
extern
Real SCIPgetDualfarkasKnapsack(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons                /**< constraint data */
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

/** lifts given cardinality inequality sum(j in C1) x_j <= |C1| to a valid inequality of the full dimensional knapsack 
 *  polytop by using uplifting for all variables not in the cover and downlifting for all variables in the cover that 
 *  are fixed to one (C2)
 */
RETCODE SCIPliftKnapsackCover(
   SCIP*            scip,               /**< SCIP data structure */
   VAR**            vars,               /**< variables in knapsack constraint */
   int              nvars,              /**< number of variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   Real*            solvals,            /**< LP values of all problem variables */
   int*             covervars,          /**< cover variables C = C2 & C1 (C2, C1 sorted by non-incr LP val then weight) */
   int*             noncovervars,       /**< noncover variables (sorted by non-incr LP val then weight) */
   int              ncovervars,         /**< number of cover variables */
   int              ncovervarsc1,       /**< number of cover variables in C1 (at the end of covervars) */
   int              ncovervarsc2,       /**< number of cover variables in C2 (at the beginning of covervars) */
   int              nnoncovervars,      /**< number of noncover variables */
   int*             liftcoefs,          /**< pointer to store lifting coefficient of variables in knapsack constraint */
   int*             liftrhs,            /**< pointer to store right hand side of the lifted cover inequality */
   Real*            liftlpval           /**< pointer to store LP solution value of lifted variables */  
   );

/** separates lifted cover inequalities for given knapsack problem */
RETCODE SCIPseparateKnapsackCover(
   SCIP*            scip,               /**< SCIP data structure */
   CONS*            cons,               /**< constraint that originates the knapsack problem */
   VAR**            vars,               /**< variables in knapsack constraint */
   int              nvars,              /**< number of variables in knapsack constraint */
   Longint*         weights,            /**< weights of variables in knapsack constraint */
   Longint          capacity,           /**< capacity of knapsack */
   int              maxnumcardlift,     /**< maximal number of cardinality inequ. lifted per sepa round (-1: unlimited) */
   int*             ncuts               /**< pointer to add up the number of found cuts */
   );

#endif
