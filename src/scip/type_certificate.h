/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2026 Zuse Institute Berlin (ZIB)                      */
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
typedef struct SCIP_AggregationInfo SCIP_AGGREGATIONINFO;         /**< aggregation information of rows */
typedef struct SCIP_MirInfo SCIP_MIRINFO;                     /**< information for split-disjunction */
typedef struct SCIP_Certificate SCIP_CERTIFICATE;                 /**< certificate data structure */
typedef struct SCIP_CertificateBound SCIP_CERTIFICATEBOUND;       /**< data structure for hashing bounds of variables in
                                                                   *   a certificate file
                                                                   */

#ifdef __cplusplus
}
#endif

#endif
