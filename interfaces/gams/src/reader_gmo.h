/* Copyright (C) GAMS Development and others 2011
 * All Rights Reserved.
 * This code is published under the Eclipse Public License.
 *
 * Author: Stefan Vigerske
 */

/**@file   reader_gmo.h
 * @brief  GMO file reader
 * @author Stefan Vigerske
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_READER_GMO_H__
#define __SCIP_READER_GMO_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gmoRec gmoRec_t;

/** includes the gmo file reader into SCIP */
extern
SCIP_RETCODE SCIPincludeReaderGmo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** sets the GMO object to use in reader
 * If GMO is set in reader, then reader does not read from file when executed, but sets up problem from GMO
 */
extern
void SCIPsetGMOReaderGmo(
   SCIP*                 scip,               /**< SCIP data structure */
   gmoRec_t*             gmo                 /**< GMO object, or NULL to reset to default behaviour */
   );

/** passes GAMS options to SCIP and initiates reading of user options file, if given in GMO */
extern
SCIP_RETCODE SCIPreadParamsReaderGmo(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** creates a SCIP problem from a GMO */
extern
SCIP_RETCODE SCIPcreateProblemReaderGmo(
   SCIP*                 scip,               /**< SCIP data structure */
   gmoRec_t*             gmo,                /**< GAMS Model Object */
   const char*           indicatorfile,      /**< name of file with indicator specification, or NULL */
   int                   mipstart            /**< how to pass initial variable levels from GMO to SCIP */
);

#ifdef __cplusplus
}
#endif

#endif
