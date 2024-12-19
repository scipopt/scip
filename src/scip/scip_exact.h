/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   scip_sol.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for solutions
 * @author Tobias Achterberg
 * @author Timo Berthold
 * @author Thorsten Koch
 * @author Alexander Martin
 * @author Marc Pfetsch
 * @author Kati Wolter
 * @author Gregor Hendel
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SCIP_EXACT_H__
#define __SCIP_SCIP_EXACT_H__


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

/* In debug mode, we include the SCIP's structure in scip.c, such that no one can access
 * this structure except the interface methods in scip.c.
 * In optimized mode, the structure is included in scip.h, because some of the methods
 * are implemented as defines for performance reasons (e.g. the numerical comparisons).
 * Additionally, the internal "set.h" is included, such that the defines in set.h are
 * available in optimized mode.
 */
#ifdef NDEBUG
#include "scip/struct_scip.h"
#include "scip/struct_stat.h"
#include "scip/set.h"
#include "scip/tree.h"
#include "scip/misc.h"
#include "scip/var.h"
#include "scip/cons.h"
#include "scip/solve.h"
#include "scip/debug.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** returns whether the solution process is arithmetically exact, i.e., not subject to roundoff errors
 *
 *  @note This feature is not supported yet!
 *
 *  @return Returns TRUE if \SCIP is exact solving mode, otherwise FALSE
 */
SCIP_EXPORT
SCIP_Bool SCIPisExactSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns whether aggreagtion is allowed to use negative slack */
SCIP_EXPORT
SCIP_Bool SCIPallowNegSlack(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** returns which method is used for computing truely valid dual bounds at the nodes ('n'eumaier and shcherbina,
 *  'v'erify LP basis, 'r'epair LP basis, 'p'roject and scale, 'e'xact LP,'i'nterval neumaier and shcherbina,
 *  e'x'act neumaier and shcherbina, 'a'utomatic); only relevant for solving the problem provably correct
 */
SCIP_EXPORT
char SCIPdualBoundMethod(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** enforce integrality of the current exact rational lp solution */ 
SCIP_EXPORT
SCIP_RETCODE SCIPcheckIntegralityExact(
   SCIP*                 scip,
   SCIP_RESULT*          result
   );

/** returns whether the certificate output is activated */
SCIP_EXPORT
SCIP_Bool SCIPisCertificateActive(
   SCIP*                 scip                /**< certificate information */
   );

/** returns whether the certificate output is activated? */
SCIP_EXPORT
void SCIPcertificateExit(
   SCIP*                 scip                /**< certificate information */
   );


/** returns certificate data structure
 *
 *  @return tolerance certificate data structure
 */
SCIP_EXPORT
SCIP_CERTIFICATE* SCIPgetCertificate(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds aggregation information to certificate for one row */
SCIP_EXPORT
SCIP_RETCODE SCIPaddCertificateAggregation(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_AGGRROW*         aggrrow,            /**< agrrrow that results from the aggregation */
   SCIP_ROW**            aggrrows,           /**< array of rows used fo the aggregation */
   SCIP_Real*            weights,            /**< array of weights */
   int                   naggrrows,          /**< length of the arrays */
   SCIP_ROW**            negslackrows,       /**< array of rows that are added implicitly with negative slack */
   SCIP_Real*            negslackweights,    /**< array of negative slack weights */
   int                   nnegslackrows       /**< length of the negative slack array */
   );

/** adds mir information (split, etc) to certificate for one row */
SCIP_EXPORT
SCIP_RETCODE SCIPaddCertificateMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** print MIR cut to certificate file */
SCIP_EXPORT
SCIP_RETCODE SCIPprintCertificateMirCut(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that needs to be certified */
   );

/** free information that is possibly still stored about this row in the certifacte structure */
SCIP_RETCODE SCIPfreeRowCertInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< a SCIP row */
   );

/** stores the active aggregation information in the certificate data structures for a row */
SCIP_EXPORT
SCIP_RETCODE SCIPstoreCertificateActiveAggregationInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that aggregation information is stored for */
   );

/** stores the active mir information in the certificate data structures for a row */
SCIP_EXPORT
SCIP_RETCODE SCIPstoreCertificateActiveMirInfo(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROW*             row                 /**< row that mirinfo is stored for */
   );

/** frees the active mir information */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeCertificateActiveMirInfo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** frees the active aggregation information */
SCIP_EXPORT
SCIP_RETCODE SCIPfreeCertificateActiveAggregationInfo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** computes a safe bound that is valid in exact rational arithmetic */
SCIP_EXPORT
SCIP_RETCODE SCIPcomputeSafeBound(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             proveinfeas,        /**< should infeasibility be proven instead */
   SCIP_Real*            safebound           /**< store the safe bound */
   );

/** forces the next lp to be solved by a rational lp solver */
SCIP_EXPORT
SCIP_RETCODE SCIPforceExactSolve(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** branches on an LP solution exactly; does not call branching rules, since fractionalities are assumed to small;
 *  if no fractional variables exist, the result is SCIP_DIDNOTRUN;
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
SCIP_RETCODE SCIPbranchLPexact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< pointer to store the result of the branching (s. branch.h) */
   );

/** adds row to exact separation storage
 *
 *  @return \ref SCIP_OKAY is returned if everything worked. Otherwise a suitable error code is passed. See \ref
 *          SCIP_Retcode "SCIP_RETCODE" for a complete list of error codes.
 *
 *  @pre This method can be called if @p scip is in one of the following stages:
 *       - \ref SCIP_STAGE_SOLVING
 */
SCIP_EXPORT
SCIP_RETCODE SCIPaddRowExact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_ROWEXACT*        rowexact            /**< exact row to add */
   );

#ifdef __cplusplus
}
#endif

#endif
