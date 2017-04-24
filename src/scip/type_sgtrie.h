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

/**@file   type_sgtrie.h
 * @brief  type definitions for signature trie
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_SGTRIE_H__
#define __SCIP_TYPE_SGTRIE_H__

#include "scip/def.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_SgTrie SCIP_SGTRIE;
typedef struct SCIP_SgTrieNode SCIP_SGTRIENODE;

typedef enum {
   SCIP_SGTRIE_EQUAL,
   SCIP_SGTRIE_SUBSET,
   SCIP_SGTRIE_SUPERSET
} SCIP_SGTRIE_QUERYTYPE;

#define SCIP_DECL_ISSETEQ(x) SCIP_Bool x (void* a, void* b, int maxdist)
#define SCIP_DECL_ISSUBSET(x) SCIP_Bool x (void* a, void* b, int maxdist)
#define SCIP_DECL_GETSIGNATURE(x) uint64_t x (void* a)

#ifdef __cplusplus
}
#endif

#endif
