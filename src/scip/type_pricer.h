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
#pragma ident "@(#) $Id: type_pricer.h,v 1.2 2003/12/08 11:51:05 bzfpfend Exp $"

/**@file   type_pricer.h
 * @brief  type definitions for variable pricers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_PRICER_H__
#define __TYPE_PRICER_H__


typedef struct Pricer PRICER;           /**< variable pricer data */
typedef struct PricerData PRICERDATA;   /**< locally defined variable pricer data */


/** destructor of variable pricer to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define DECL_PRICERFREE(x) RETCODE x (SCIP* scip, PRICER* pricer)

/** initialization method of variable pricer (called when problem solving starts and pricer is active)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - pricer          : the variable pricer itself
 */
#define DECL_PRICERINIT(x) RETCODE x (SCIP* scip, PRICER* pricer)

/** deinitialization method of variable pricer (called when problem solving exits and pricer is active)
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



#include "def.h"
#include "type_retcode.h"
#include "type_scip.h"


#endif
