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
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   type_certificate.h
 * @brief  type definitions for certificate output
 * @author Ambros Gleixner
 * @author Daniel Steffy
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_CERTIFICATE_H__
#define __SCIP_TYPE_CERTIFICATE_H__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Certnodedata SCIP_CERTNODEDATA;               /**< branching information for nodes */
typedef struct SCIP_Certificate SCIP_CERTIFICATE;                 /**< certificate data structure */
typedef struct SCIP_CertificateBound SCIP_CERTIFICATEBOUND;       /**< data structure for hashing bounds of variables in
                                                                   *   a certificate file
                                                                   */

#ifdef __cplusplus
}
#endif

#endif
