/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2017 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   aggrrow.h
 * @ingroup PUBLICCOREAPI
 * @brief  methods for the aggregation rows
 * @author Jakob Witzig
 * @author Robert Lion Gottwald
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTS_H__
#define __SCIP_CUTS_H__

#include "scip/def.h"
#include "scip/set.h"
#include "scip/type_cuts.h"

#ifdef __cplusplus
extern "C" {
#endif

/**@addtogroup PublicCutMethods
 *
 * @{
 */

/** create an empty the aggregation row */
extern
SCIP_RETCODE SCIPaggrRowCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow             /**< pointer to return the aggregation row */
   );

/** free a the aggregation row */
extern
void SCIPaggrRowFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow             /**< pointer to the aggregation row that should be freed */
   );

/** output aggregation row to file stream */
extern
void SCIPaggrRowPrint(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< pointer to return aggregation row */
   FILE*                 file                /**< output file (or NULL for standard output) */
   );

/** copy the aggregation row */
extern
SCIP_RETCODE SCIPaggrRowCopy(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW**        aggrrow,            /**< pointer to return the aggregation row */
   SCIP_AGGRROW*         source              /**< source the aggregation row */
   );

/** adds given value to the right-hand side of the aggregation row */
extern
SCIP_RETCODE SCIPaggrRowAddData(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_VAR**            vars,               /**< variable array */
   SCIP_Real*            coefs,              /**< variable coefficients */
   int                   nvars,              /**< size of variable and coefficient array */
   SCIP_Real             rhs,                /**< right-hand side of the row */
   SCIP_Real             scale               /**< scalar to apply */
   );

/** add weighted row to the aggregation row */
extern
SCIP_RETCODE SCIPaggrRowAddRow(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_ROW*             row,                /**< row to add to the aggregation row */
   SCIP_Real             weight,             /**< scale for adding given row to the aggregation row */
   int                   sidetype            /**< specify row side type (-1 = lhs, 0 = automatic, 1 = rhs) */
   );

/** add weighted constraint to the aggregation row */
extern
SCIP_RETCODE SCIPaggrRowAddCustomCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   int*                  inds,               /**< variable problem indices in constraint to add to the aggregation row */
   SCIP_Real*            vals,               /**< values of constraint to add to the aggregation row */
   int                   len,                /**< length of constraint to add to the aggregation row */
   SCIP_Real             rhs,                /**< right hand side of constraint to add to the aggregation row */
   SCIP_Real             weight,             /**< (positive) scale for adding given constraint to the aggregation row */
   int                   rank,               /**< rank to use for given constraint */
   SCIP_Bool             local               /**< is constraint only valid locally */
   );

/** deletes variable at position @pos and updates mapping between variable indices and sparsity pattern */
extern
void SCIPaggrRowDelCoef(
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   int                   pos,                /**< position that should be removed */
   int*                  positions           /**< mapping between variable indices and sparsity pattern (or NULL) */
   );

/** change the right-hand side of the aggregation row */
extern
void SCIPaggrRowAddRhs(
   SCIP_AGGRROW*         aggrrow,            /**< aggregation row */
   SCIP_Real             value               /**< value to add to the right-hand side */
   );

/** clear all entries in the aggregation row but do not free the internal memory */
extern
void SCIPaggrRowClear(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** aggregate rows using the given weights; the current content of the aggregation
 *  row, @aggrow, gets overwritten
 */
extern
SCIP_RETCODE SCIPaggrRowSumRows(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_Real*            weights,            /**< row weights in row summation */
   int*                  rowinds,            /**< array to store indices of non-zero entries of the weights array, or NULL */
   int                   nrowinds,           /**< number of non-zero entries in weights array, -1 if rowinds is NULL */
   SCIP_Bool             sidetypebasis,      /**< choose sidetypes of row (lhs/rhs) based on basis information? */
   SCIP_Bool             allowlocal,         /**< should local rows be used? */
   int                   negslack,           /**< should negative slack variables be used? (0: no, 1: only for integral rows, 2: yes) */
   int                   maxaggrlen,         /**< maximal number of non-zeros in the aggregation row */
   SCIP_Bool*            valid               /**< is the aggregation valid */
   );

/** removes all (close enough to) zero entries in the aggregation row */
extern
void SCIPaggrRowRemoveZeros(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** safely removes variables with small coefficients from the aggregation row */
extern
void SCIPaggrRowCleanup(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** get number of aggregated rows */
extern
int SCIPaggrRowGetNRows(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** get array with lp positions of aggregated rows */
extern
int* SCIPaggrRowGetRowInds(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** checks whether a given row has been added to the aggregation row */
extern
SCIP_Bool SCIPaggrRowHasRowBeenAdded(
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_ROW*             row                 /**< row for which it is checked whether it has been added to the aggregation */
   );

/** gets the min and max absolute value of the weights used to aggregate the rows;
 *  must not be called for empty aggregation rows
 */
extern
void SCIPaggrRowGetAbsWeightRange(
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_Real*            minabsrowweight,    /**< pointer to store smallest absolute value of weights used for aggregating rows */
   SCIP_Real*            maxabsrowweight     /**< pointer to store largest absolute value of weights used for aggregating rows */
   );

/** gets the array of corresponding variable problem indices for each non-zero in the aggregation row */
extern
int* SCIPaggrRowGetInds(
    SCIP_AGGRROW*        aggrrow
   );

/** gets the array of non-zero values in the aggregation row */
extern
SCIP_Real* SCIPaggrRowGetVals(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** gets the number of non-zeros in the aggregation row */
extern
int SCIPaggrRowGetNNz(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** gets the rank of the aggregation row */
extern
int SCIPaggrRowGetRank(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** checks if the aggregation row is only valid locally */
extern
SCIP_Bool SCIPaggrRowIsLocal(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** gets the right hand side of the aggregation row */
extern
SCIP_Real SCIPaggrRowGetRhs(
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
   );

/** calculates an MIR cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in an MIR cut.
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_RETCODE SCIPcalcMIR(
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
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to the aggrrow; must be positive */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute an MIR cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store the efficacy of the cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether the returned coefficients are a valid MIR cut */
   );

/** calculates an MIR cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in an MIR cut. The function uses a cut generation heuristic which tries different scaling
 *  factors and complementations of the variables to improve the cut's efficacy.
 *  For further details we refer to:
 *
 *  Marchand, H., & Wolsey, L. A. (2001). Aggregation and mixed integer rounding to solve MIPs.
 *  Operations research, 49(3), 363-371.
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
SCIP_RETCODE SCIPcutGenerationHeuristicCMIR(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   int                   maxtestdelta,       /**< maximum number of deltas to test */
   int*                  boundsfortrans,     /**< bounds that should be used for transformed variables: vlb_idx/vub_idx,
                                              *   -1 for global lb/ub, -2 for local lb/ub, or -3 for using closest bound;
                                              *   NULL for using closest bound for all variables */
   SCIP_BOUNDTYPE*       boundtypesfortrans, /**< type of bounds that should be used for transformed variables;
                                              *   NULL for using closest bound for all variables */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce MIR cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce MIR cut for */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute MIR cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store efficacy of best cut; only cuts that are strictly better than the value of
                                              *   this efficacy on input to this function are returned */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid and efficacious cut was returned */
   );

/** calculates a lifted simple generalized flow cover cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in an MIR cut.
 *  For further details we refer to:
 *
 *  Gu, Z., Nemhauser, G. L., & Savelsbergh, M. W. (1999). Lifted flow cover inequalities for mixed 0-1 integer programs.
 *  Mathematical Programming, 85(3), 439-467.
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
SCIP_RETCODE SCIPcalcFlowCover(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute flow cover cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store the efficacy of the cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid cut was returned */
   );

/** calculates a strong CG cut out of the weighted sum of LP rows given by an aggregation row; the
 *  aggregation row must not contain non-zero weights for modifiable rows, because these rows cannot
 *  participate in a strongcg cut
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
SCIP_RETCODE SCIPcalcStrongCG(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol,                /**< the solution that should be separated, or NULL for LP solution */
   SCIP_Real             boundswitch,        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Bool             usevbds,            /**< should variable bounds be used in bound transformation? */
   SCIP_Bool             allowlocal,         /**< should local information allowed to be used, resulting in a local cut? */
   SCIP_Real             minfrac,            /**< minimal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             maxfrac,            /**< maximal fractionality of rhs to produce strong CG cut for */
   SCIP_Real             scale,              /**< additional scaling factor multiplied to all rows */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row to compute a flow cover cut for */
   SCIP_Real*            cutcoefs,           /**< array to store the non-zero coefficients in the cut */
   SCIP_Real*            cutrhs,             /**< pointer to store the right hand side of the cut */
   int*                  cutinds,            /**< array to store the problem indices of variables with a non-zero coefficient in the cut */
   int*                  cutnnz,             /**< pointer to store the number of non-zeros in the cut */
   SCIP_Real*            cutefficacy,        /**< pointer to store the efficacy of the cut, or NULL */
   int*                  cutrank,            /**< pointer to return rank of generated cut */
   SCIP_Bool*            cutislocal,         /**< pointer to store whether the generated cut is only valid locally */
   SCIP_Bool*            success             /**< pointer to store whether a valid cut was returned */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
