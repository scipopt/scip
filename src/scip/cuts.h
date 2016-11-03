/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cuts.h
 * @brief  Header file for methods used to generate and strengthen cuts
 * @author Jakob Witzig
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTS_H__
#define __SCIP_CUTS_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_cons.h"
#include "scip/type_set.h"
#include "scip/type_stat.h"
#include "scip/type_misc.h"
#include "scip/type_prob.h"
#include "scip/type_sol.h"
#include "scip/type_lp.h"

#ifdef __cplusplus
extern "C" {
#endif

/** calculates a strong CG cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0 because
 *  these rows cannot participate in an MIR cut.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
extern
SCIP_RETCODE SCIPcutsCalcStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce strong CG cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   int*                  inds,               /**< indices of non-zero entries in weights array, or NULL */
   int                   ninds,              /**< number of indices of non-zero entries in weights array, -1 if inds is
                                              *   NULL */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mircoef,            /**< array to store strong CG coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the strong CG row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid strong CG cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the returned cut is only valid locally */
   int*                  cutrank             /**< pointer to store the rank of the returned cut; or NULL */
   );

/** calculates an MIR cut out of the weighted sum of LP rows; The weights of modifiable rows are set to 0.0, because
 *  these rows cannot participate in an MIR cut.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
extern
SCIP_RETCODE SCIPcutsCalcLpMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   int                   maxmksetcoefs,      /**< maximal number of nonzeros allowed in aggregated base inequality */
   SCIP_Real             maxweightrange,     /**< maximal valid range max(|weights|)/min(|weights|) of row weights */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real*            weights,            /**< row weights in row summation; some weights might be set to zero */
   SCIP_Real             maxweight,          /**< largest magnitude of weights; set to -1.0 if sparsity information is
                                              *   unknown */
   int*                  weightinds,         /**< sparsity pattern of weights; size nrowinds; NULL if sparsity info is
                                              *   unknown */
   int                   nweightinds,        /**< number of nonzeros in weights; -1 if rowinds is NULL */
   int                   rowlensum,          /**< total number of nonzeros in used rows (row associated with nonzero weight coefficient); -1 if unknown */
   int*                  sidetypes,          /**< specify row side type (-1 = lhs, 0 = unkown, 1 = rhs) or NULL for automatic choices */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Bool*            mksetcoefsvalid,    /**< pointer to store whether mixed knapsack set coefficients are valid; or NULL */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   SCIP_Real*            cutactivity,        /**< pointer to store the activity of the resulting cut */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the returned cut is only valid locally */
   int*                  cutrank             /**< pointer to store the rank of the returned cut; or NULL */
   );

/** applies the MIR function on a constraint; the constraint is given by pairs of variables and coefficients and a rhs.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
extern
SCIP_RETCODE SCIPcutsApplyMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Bool             fixintegralrhs,     /**< should complementation tried to be adjusted such that rhs gets fractional? */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: 0 vlb_idx/vub_idx,
                                              *   -1 for global lb/ub or -2 for local lb/ub
                                              */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_Real*            mksetcoefs,         /**< array to store mixed knapsack set coefficients: size nvars; or NULL */
   SCIP_Bool*            mksetcoefsvalid,    /**< pointer to store whether mixed knapsack set coefficients are valid; or NULL */
   SCIP_Real*            mircoef,            /**< array to store MIR coefficients: must be of size SCIPgetNVars() */
   SCIP_Real*            mirrhs,             /**< pointer to store the right hand side of the MIR row */
   int*                  varinds,            /**< array of variable indices with a mircoef != 0 */
   int*                  nvarinds,           /**< number of variables indices in varinds array */
   SCIP_Real*            minact,             /**< pointer to store the minimal activity */
   SCIP_Bool*            varused,            /**< array to store whether a variable has a mircoef != 0 */
   SCIP_Bool*            success,            /**< pointer to store whether the returned coefficients are a valid MIR cut */
   SCIP_Bool*            islocal             /**< pointer to store whether the returned constraint is only valid locally */
   );

/** removes all nearly-zero coefficients from MIR row and relaxes the right hand side accordingly in order to prevent
 *  numerical rounding errors
 */
void SCIPcutsCleanupRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Real*            coefs,              /**< array to store MIR coefficients: must be of size nvars */
   SCIP_Real*            rhs,                /**< pointer to store the right hand side of the MIR row */
   SCIP_Bool*            varused,            /**< array to flag variables that appear in the MIR constraint */
   int*                  varinds,            /**< sparsity pattern of non-zero MIR coefficients */
   int*                  nvarinds,           /**< pointer to number of non-zero MIR coefficients */
   SCIP_Bool             islocal             /**< is the row only valid locally? */
   );

#ifdef __cplusplus
}
#endif

#endif
