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

/**@file   type_dcmp.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for decompositions and the decomposition store
 * @author Gregor Hendel
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef SRC_SCIP_TYPE_DECOMP_H_
#define SRC_SCIP_TYPE_DECOMP_H_

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Decomp SCIP_DECOMP;
typedef struct SCIP_DecompStore SCIP_DECOMPSTORE;


#define SCIP_DECOMP_LINKVAR -1      /**< special label for linking variables */
#define SCIP_DECOMP_LINKCONS -2     /**< special label for linking constraints */

#ifdef __cplusplus
}
#endif

#endif
