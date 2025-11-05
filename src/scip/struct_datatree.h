/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
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

/**@file   struct_datatree.h
 * @ingroup INTERNALAPI
 * @brief  data structures for data trees
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_DATATREE_H__
#define __SCIP_STRUCT_DATATREE_H__

#include "scip/def.h"
#include "scip/type_set.h"

#ifdef __cplusplus
extern "C" {
#endif

/** union to store any possible value of a SCIP_DATATREE */
typedef union
{
   SCIP_Bool            as_bool;            /**< value interpreted as a boolean */
   SCIP_Longint         as_long;            /**< value interpreted as a long integer */
   SCIP_Real            as_real;            /**< value interpreted as a SCIP_Real */
   char*                as_string;          /**< value interpreted as a string */
   SCIP_Bool*           as_boolarray;       /**< value interpreted as an array of boolean */
   SCIP_Longint*        as_longarray;       /**< value interpreted as an array of long integer */
   SCIP_Real*           as_realarray;       /**< value interpreted as an array of SCIP_Real */
   char**               as_stringarray;     /**< value interpreted as an array of string */
   SCIP_DATATREE*       as_dtree;           /**< value interpreted as a SCIP_DATATREE */
} SCIP_DATATREEVALUEUNION;

/** structure representing a value in a SCIP_DATATREE */
typedef struct
{
   SCIP_DATATREE_VALUETYPE type;            /**< type of the value */
   SCIP_DATATREEVALUEUNION data;            /**< storage for the value (as union) */
   int                     nvalues;         /**< length of arrays if array type */
} SCIP_DATATREEVALUE;

/** structure representing an item in a SCIP_DATATREE */
typedef struct
{
   char*                name;               /**< name of data item */
   SCIP_DATATREEVALUE   value;              /**< value of data item */
} SCIP_DATATREEITEM;

/** structure representing a node in a tree of data */
struct SCIP_Datatree
{
   SCIP_DATATREEITEM*   items;              /**< array of data items */
   int                  nitems;             /**< number of items currently stored */
   int                  itemssize;          /**< capacity (size) of the data tree (maximum number of items it can hold) */
};

#ifdef __cplusplus
}
#endif

#endif
