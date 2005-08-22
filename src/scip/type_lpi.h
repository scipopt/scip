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
#pragma ident "@(#) $Id: type_lpi.h,v 1.17 2005/08/22 18:35:53 bzfpfend Exp $"

/**@file   type_lpi.h
 * @brief  type definitions for specific SCIP_LP solvers interface
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_LPI_H__
#define __SCIP_TYPE_LPI_H__


/** objective sense */
enum SCIP_ObjSen
{
   SCIP_OBJSEN_MAXIMIZE = -1,           /**< maximize objective function */
   SCIP_OBJSEN_MINIMIZE = +1            /**< minimize objective function */
};
typedef enum SCIP_ObjSen SCIP_OBJSEN;

/** SCIP_LP solver parameters */
enum SCIP_LPParam
{
   SCIP_LPPAR_FROMSCRATCH    =  0,      /**< solver should start from scratch at next call? */
   SCIP_LPPAR_FASTMIP        =  1,      /**< fast mip setting of SCIP_LP solver */
   SCIP_LPPAR_SCALING        =  2,      /**< should SCIP_LP solver use scaling? */
   SCIP_LPPAR_PRESOLVING     =  3,      /**< should SCIP_LP solver use presolving? */
   SCIP_LPPAR_PRICING        =  4,      /**< pricing strategy */
   SCIP_LPPAR_LPINFO         =  5,      /**< should SCIP_LP solver output information to the screen? */
   SCIP_LPPAR_FEASTOL        =  6,      /**< feasibility tolerance for primal variables and slacks */
   SCIP_LPPAR_DUALFEASTOL    =  7,      /**< feasibility tolerance for dual variables and reduced costs */
   SCIP_LPPAR_BARRIERCONVTOL =  8,      /**< convergence tolerance used in barrier algorithm */
   SCIP_LPPAR_LOBJLIM        =  9,      /**< lower objective limit */
   SCIP_LPPAR_UOBJLIM        = 10,      /**< upper objective limit */
   SCIP_LPPAR_LPITLIM        = 11,      /**< SCIP_LP iteration limit */
   SCIP_LPPAR_LPTILIM        = 12,      /**< SCIP_LP time limit */
   SCIP_LPPAR_MARKOWITZ      = 13       /**< Markowitz tolerance */
};
typedef enum SCIP_LPParam SCIP_LPPARAM;

/** SCIP_LP pricing strategy */
enum SCIP_Pricing
{
   SCIP_PRICING_AUTO        = 0,        /**< the SCIP_LP solver should use its preferred strategy */
   SCIP_PRICING_FULL        = 1,        /**< full pricing */
   SCIP_PRICING_PARTIAL     = 2,        /**< partial pricing */
   SCIP_PRICING_STEEP       = 3,        /**< steepest edge pricing */
   SCIP_PRICING_STEEPQSTART = 4,        /**< steepest edge pricing without initial dual norms */
   SCIP_PRICING_DEVEX       = 5         /**< devex pricing */
};
typedef enum SCIP_Pricing SCIP_PRICING;

/** basis status for columns and rows */
enum SCIP_BaseStat
{
   SCIP_BASESTAT_LOWER = 0,             /**< (slack) variable is at its lower bound */
   SCIP_BASESTAT_BASIC = 1,             /**< (slack) variable is basic */
   SCIP_BASESTAT_UPPER = 2,             /**< (slack) variable is at its upper bound */
   SCIP_BASESTAT_ZERO  = 3              /**< free variable is non-basic and set to zero */
};
typedef enum SCIP_BaseStat SCIP_BASESTAT;

typedef struct SCIP_LPi SCIP_LPI;                 /**< solver dependent SCIP_LP interface */
typedef struct SCIP_LPiState SCIP_LPISTATE;       /**< complete SCIP_LP state (i.e. basis information, dual norms) */


#endif
