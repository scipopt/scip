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

/**@file   type_bandit.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for primal heuristics
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 *  This file defines the interface for primal heuristics implemented in C.
 *
 *  - \ref BANDIT "Instructions for implementing a primal heuristic"
 *  - \ref PRIMALHEURISTICS "List of available primal heuristics"
 *  - \ref scip::ObjHeur "C++ wrapper class"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_BANDIT_H__
#define __SCIP_TYPE_BANDIT_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "scip/type_timing.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data structure for bandit algorithms */
typedef struct SCIP_Bandit SCIP_BANDIT;

/** virtual table for bandit callbacks */
typedef struct SCIP_BanditVTable SCIP_BANDITVTABLE;

/** data structure for specific bandit algorithm implementation */
typedef struct SCIP_BanditData SCIP_BANDITDATA;

/*
 * callbacks for bandit VTable
 */

/** callback to free bandit specific data structures */
#define SCIP_DECL_BANDITFREE(x) SCIP_RETCODE x (  \
   SCIP_BANDIT*          bandit,                  \
   BMS_BLKMEM*           blkmem                   \
)

/** selection callback for bandit selector */
#define SCIP_DECL_BANDITSELECT(x) SCIP_RETCODE x ( \
   SCIP_BANDIT*          bandit,                   \
   int*                  selection                 \
)

/** update callback for bandit algorithms */
#define SCIP_DECL_BANDITUPDATE(x) SCIP_RETCODE x ( \
   SCIP_BANDIT*          bandit,                   \
   int                   selection,                \
   SCIP_Real             score                     \
)

/** reset callback for bandit algorithms */
#define SCIP_DECL_BANDITRESET(x) SCIP_RETCODE x (  \
   BMS_BUFMEM*           bufmem,                   \
   SCIP_BANDIT*          bandit,                   \
   SCIP_Real*            priorities                \
)

#ifdef __cplusplus
}
#endif

#endif
