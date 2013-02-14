/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2013 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   sepa_cmir.h
 * @ingroup SEPARATORS
 * @brief  complemented mixed integer rounding cuts separator (Marchand's version)
 * @author Kati Wolter
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SEPA_CMIR_H__
#define __SCIP_SEPA_CMIR_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** Performs the cut generation heuristic of the c-MIR separation algorithm, i.e., tries to generate a c-MIR cut which is
 *  valid for the mixed knapsack set corresponding to the current aggregated constraint. Cuts will only be added here if 
 *  no pointer to store best scaling factor delta is given.
 */
EXTERN
SCIP_RETCODE SCIPcutGenerationHeuristicCmir(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPA*            sepa,               /**< separator */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real*            varsolvals,         /**< LP solution value of all variables in LP */
   int                   maxtestdelta,       /**< maximal number of different deltas to try (-1: unlimited) */
   SCIP_Real*            rowweights,         /**< weight of rows in aggregated row */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Bool             trynegscaling,      /**< should negative values also be tested in scaling? */
   SCIP_Bool             cutremovable,       /**< should the cut be removed from the LP due to aging or cleanup? */
   const char*           cutclassname,       /**< name of cut class to use for row names */
   int*                  ncuts,              /**< pointer to count the number of generated cuts */
   SCIP_Real*            delta,              /**< pointer to store best delta found; NULL, if cut should be added here */
   SCIP_Bool*            deltavalid          /**< pointer to store whether best delta value is valid or NULL */
   );

/** creates the cmir separator and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeSepaCmir(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
