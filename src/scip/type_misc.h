/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: type_misc.h,v 1.2 2004/02/04 17:27:49 bzfpfend Exp $"

/**@file   type_misc.h
 * @brief  type definitions for miscellaneous datastructures
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_MISC_H__
#define __TYPE_MISC_H__


typedef struct PQueue PQUEUE;           /**< priority queue */
typedef struct HashTable HASHTABLE;     /**< hash table */
typedef struct HashList HASHLIST;       /**< element list to store in a hash table */
typedef struct RealArray REALARRAY;     /**< dynamic array for storing Real values */
typedef struct IntArray INTARRAY;       /**< dynamic array for storing int values */
typedef struct BoolArray BOOLARRAY;     /**< dynamic array for storing Bool values */



/** compares two element indices
 *  result:
 *    < 0: ind1 comes before (is better than) ind2
 *    = 0: both indices have the same value
 *    > 0: ind2 comes after (is worse than) ind2
 */
#define DECL_SORTINDCOMP(x) int x (void* dataptr, int ind1, int ind2)

/** compares two data element pointers
 *  result:
 *    < 0: elem1 comes before (is better than) elem2
 *    = 0: both elements have the same value
 *    > 0: elem2 comes after (is worse than) elem2
 */
#define DECL_SORTPTRCOMP(x) int x (void* elem1, void* elem2)

/** gets the key of the given element */
#define DECL_HASHGETKEY(x) void* x (void* elem)

/** returns TRUE iff both keys are equal */
#define DECL_HASHKEYEQ(x) Bool x (void* key1, void* key2)

/** returns the hash value of the key */
#define DECL_HASHKEYVAL(x) unsigned int x (void* key)



#include "def.h"


#endif
