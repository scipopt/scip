/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_certificate.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for certified solving
 * @author Leon Eifler
 * @author Ambros Gleixner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_CERTIFICATE_H__
#define __SCIP_SCIP_CERTIFICATE_H__


#include "scip/def.h"
#include "scip/type_cuts.h"
#include "scip/type_cons.h"
#include "scip/type_heur.h"
#include "scip/type_retcode.h"
#include "scip/type_scip.h"
#include "scip/type_sol.h"
#include "scip/type_var.h"
#include "scip/type_certificate.h"
#include "scip/type_lpexact.h"

#ifdef __cplusplus
extern "C" {
#endif

/** returns whether certificate output is activated */
SCIP_EXPORT
SCIP_Bool SCIPisCertified(
   SCIP*                 scip                /**< certificate information */
   );

/** should the certificate track bound changes?
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_Bool SCIPshouldCertificateTrackBounds(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** prints constraint to certificate
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_INITSOLVE
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcertifyCons(
   SCIP*                 scip,               /**< certificate information */
   SCIP_Bool             isorigfile,         /**< should the original solution be printed or in transformed space */
   const char*           consname,           /**< name of the constraint */
   const char            sense,              /**< sense of the constraint, i.e., G, L, or E */
   SCIP_RATIONAL*        side,               /**< left/right-hand side */
   int                   len,                /**< number of nonzeros */
   int*                  ind,                /**< index array */
   SCIP_RATIONAL**       val                 /**< coefficient array */
   );

/** prints activity bound to proof section
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcertifyActivityVarBound(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound (upper/lower) */
   SCIP_Real             newbound,           /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   SCIP_Bool             ismaxactivity,      /**< TRUE for maxactivity, FALSE for minactivity */
   SCIP_CONS*            constraint,         /**< the constraint */
   SCIP_VAR*             variable,           /**< the variable */
   SCIP_ROWEXACT*        row,                /**< the  corresponding row, or NULL if constraint has no row representation */
   SCIP_RATIONAL**       vals,               /**< value array */
   SCIP_RATIONAL*        lhs,                /**< lhs of the constraint */
   SCIP_RATIONAL*        rhs,                /**< rhs of the constraint */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< number of values */
   );

/** prints activity bound to proof section
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcertifyActivityVarBoundExact(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           linename,           /**< name of the unsplitting line */
   SCIP_BOUNDTYPE        boundtype,          /**< type of bound (upper/lower) */
   SCIP_RATIONAL*        newbound,           /**< pointer to lower bound on the objective, NULL indicating infeasibility */
   SCIP_Bool             ismaxactivity,      /**< TRUE for maxactivity, FALSE for minactivity */
   SCIP_CONS*            constraint,         /**< the constraint */
   SCIP_VAR*             variable,           /**< the variable */
   SCIP_ROWEXACT*        row,                /**< the  corresponding row, or NULL if constraint has no row representation */
   SCIP_RATIONAL**       vals,               /**< value array */
   SCIP_RATIONAL*        lhs,                /**< lhs of the constraint */
   SCIP_RATIONAL*        rhs,                /**< rhs of the constraint */
   SCIP_VAR**            vars,               /**< variable array */
   int                   nvars               /**< number of values */
   );

/** prints activity conflict to  certificate file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcertifyActivityConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint */
   SCIP_ROWEXACT*        row,                /**< row */
   SCIP_RATIONAL*        lhs,                /**< lhs of the constraint */
   SCIP_RATIONAL*        rhs,                /**< rhs of the constraint */
   int                   nvals,              /**< number of values */
   SCIP_RATIONAL**       vals,               /**< values */
   SCIP_VAR**            vars,               /**< variables */
   SCIP_RATIONAL*        diff,               /**< difference */
   SCIP_Bool             userhs              /**< is rhs */
   );

/** adds aggregation information to certificate for one row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddCertificateAggrInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_ROW**            negslackrows,       /**< array of rows that are added implicitly with negative slack */
   SCIP_Real*            negslackweights,    /**< array of negative slack weights */
   int                   nnegslackrows       /**< length of the negative slack array */
   );

/** stores the active aggregation information in the certificate data structures for a row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPstoreCertificateActiveAggrInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that aggregation information is stored for */
   );

/** frees the active aggregation information
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeCertificateActiveAggrInfo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds mir information (split, etc) to certificate for one row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddCertificateMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** stores the active mir information in the certificate data structures for a row
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPstoreCertificateActiveMirInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that mirinfo is stored for */
   );

/** print MIR cut to certificate file
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcertifyMirCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that needs to be certified */
   );

/** frees the active mir information
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 *
 *  See \ref SCIP_Stage "SCIP_STAGE" for a complete list of all possible solving stages.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeCertificateActiveMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
