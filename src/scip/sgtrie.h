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

/** creates the signature trie datastructure */
extern
SCIP_RETCODE SCIPsgtrieCreate(
   SCIP_SGTRIE**         sgtrie,             /**< pointer to return the signature trie datastructure */
   BMS_BLKMEM*           blkmem,             /**< block memory */
   BMS_BUFMEM*           bufmem,             /**< buffer memory */
   SCIP_DECL_GETSIGNATURE ((*getsignature)), /**< callback to retrieve the signature of a set */
   SCIP_DECL_SETCMP      ((*setcmp))         /**< callback to perform different comparison operations on two sets, may be NULL
                                              *   if FALSE positives are acceptable */
   );

/** frees the signature trie datastructure, does not free the elements contained in the signature trie */
extern
void SCIPsgtrieFree(
   SCIP_SGTRIE**         sgtrie              /**< pointer to the signature trie datastructure */
   );

/** inserts a set into the signature trie data structure */
extern
SCIP_RETCODE SCIPsgtrieInsert(
   SCIP_SGTRIE*          sgtrie,             /**< the signature trie data structure */
   void*                 set                 /**< the set */
   );

/** removes a set from the signature trie data structure */
extern
SCIP_RETCODE SCIPsgtrieRemove(
   SCIP_SGTRIE*          sgtrie,             /**< the signature trie data structure */
   void*                 set                 /**< the set */
   );

/** returns number of sets currently stored in the trie */
extern
int SCIPsgtrieGetNElems(
   SCIP_SGTRIE*          sgtrie              /**< the signature trie data structure */
   );

/** finds all subsets of a given set that are stored in the signature trie datastructure.
 *  If the setcmp callback given upon creation of the signature trie was NULL, the matches
 *  may contain false positives.
 */
extern
SCIP_RETCODE SCIPsgtrieFindSubsets(
   SCIP_SGTRIE*          sgtrie,             /**< the signature trie data structure */
   void*                 set,                /**< the set */
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   );

/** finds all sets stored in the signature trie datastructure that are subsets of a given set
 *  after removing at most one element.
 *  If the setcmp callback was set to NULL upon creation of the signature trie, the matches
 *  may contain false positives.
 */
extern
SCIP_RETCODE SCIPsgtrieFindSubsetsPlusOne(
   SCIP_SGTRIE*          sgtrie,             /**< the signature trie data structure */
   void*                 set,                /**< the set */
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   );

/** finds all supersets of a given set that are stored in the signature trie datastructure.
 *  If the setcmp callback given upon creation of the signature trie was NULL, the matches
 *  may contain false positives.
 */
extern
SCIP_RETCODE SCIPsgtrieFindSupersets(
   SCIP_SGTRIE*          sgtrie,             /**< the signature trie data structure */
   void*                 set,                /**< the set */
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   );

/** finds all sets stored in the signature trie datastructure that are supersets of a given set
 *  after removing at most one element.
 *  If the setcmp callback was set to NULL upon creation of the signature trie, the matches
 *  may contain false positives.
 */
extern
SCIP_RETCODE SCIPsgtrieFindSupersetsPlusOne(
   SCIP_SGTRIE*          sgtrie,             /**< the signature trie data structure */
   void*                 set,                /**< the set */
   void**                matches,            /**< buffer to store matches, must be big enough to hold all elements currently stored in the sgtrie */
   int*                  nmatches            /**< pointer to store how many matches where found */
   );

/** searches for a set stored in the signature trie datastructure that is equal to a given set.
 *  If the setcmp callback was set to NULL upon creation of the signature trie, the match
 *  may be a false positives.
 */
extern
void* SCIPsgtrieFindEq(
   SCIP_SGTRIE*          sgtrie,             /**< the signature trie data structure */
   void*                 set                 /**< the set */
   );

#ifdef __cplusplus
}
#endif

#endif
