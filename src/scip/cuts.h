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
 * @author Robert Lion Gottwald
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CUTROW_H__
#define __SCIP_CUTROW_H__

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
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_AGGRROW**        aggrrow              /**< pointer to return the aggregation row */
   );

/** free a the aggregation row */
extern
void SCIPaggrRowFree(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_AGGRROW**        aggrrow              /**< pointer to the aggregation row that should be freed */
   );

/** copy a the aggregation row */
extern
SCIP_RETCODE SCIPaggrRowCopy(
   SCIP*                 scip,               /**< SCIP datastructure */
    SCIP_AGGRROW**         aggrrow,             /**< pointer to return the aggregation row */
    SCIP_AGGRROW*          source              /**< source the aggregation row */
   );

/** add scaled row to the aggregation row */
extern
SCIP_RETCODE SCIPaggrRowAddRow(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_AGGRROW*         aggrrow,            /**< the aggregation row */
   SCIP_ROW*             row,                /**< row to add to the aggregation row */
   SCIP_Real             scale,              /**< scale for adding given row to the aggregation row */
   int                   sidetype            /**< specify row side type (-1 = lhs, 0 = automatic, 1 = rhs) */
   );

/** removes almost zero entries and relaxes the sides of the aggregation row accordingly */
extern
void SCIPaggrRowCleanup(
   SCIP*                 scip,               /**< SCIP datastructure */
   SCIP_AGGRROW*         aggrrow             /**< the aggregation row */
  );

/** gets the array of corresponding variable problem indices for each non-zero in the aggregation row */
extern
int* SCIPaggrRowGetInds(
    SCIP_AGGRROW*          aggrrow
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

/** computes an upper bound for the number of non-zeros of cuts computed from this aggregation row; the arrays
 *  for returning a cut passed to any of the functions in this file must have at least this size. */
extern
int SCIPaggrRowGetCutsMaxNNz(
   SCIP*                 scip,               /**< SCIP datastructure */
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
   SCIP_Bool*            success             /**< pointer to store whether a valid cut was returned */
   );

/* @} */

#ifdef __cplusplus
}
#endif

#endif
