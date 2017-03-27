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

/**@file   sgtrie.h
 * @brief  interface for signature trie
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_SGTRIE_H__
#define __SCIP_SGTRIE_H__

#include "scip/type_sgtrie.h"
#include "blockmemshell/memory.h"

#ifdef __cplusplus
extern "C" {
#endif

/* macro to update a signature with a given element identifier */
#define UPDATE_SIGNATURE(signature, elemid) do { \
   (signature) |= (UINT64_C(0x8000000000000000)>>((UINT32_C(0x9e3779b9) * ((uint32_t)(elemid)))>>26)); } while(0)

SCIP_RETCODE SCIPsgtrieInsert(
   SCIP_SGTRIE*          sgtrie,
   uint64_t              signature,
   void*                 data
   );

SCIP_RETCODE SCIPsgtrieRemove(
   SCIP_SGTRIE*          sgtrie,
   uint64_t              signature,
   void*                 data
   );

/** returns number of elements currently stored in the trie */
int SCIPsgtrieGetNElems(
   SCIP_SGTRIE*          sgtrie
   );

SCIP_RETCODE SCIPsgtrieFindSubsetCands(
   SCIP_SGTRIE*          sgtrie,
   uint64_t              signature,
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   );

SCIP_RETCODE SCIPsgtrieFindSubsetPlusOneCands(
   SCIP_SGTRIE*          sgtrie,
   uint64_t              signature,
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   );

void* SCIPsgtrieFindEq(
   SCIP_SGTRIE*          sgtrie,
   uint64_t              signature,
   void*                 data
   );

SCIP_RETCODE SCIPsgtrieCreate(
   SCIP_SGTRIE**         sgtrie,
   BMS_BLKMEM*           blkmem,
   BMS_BUFMEM*           bufmem
   );

void SCIPsgtrieFree(
   SCIP_SGTRIE**         sgtrie
   );

#ifdef __cplusplus
}
#endif

#endif
