/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
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

/**@file    type_datatree.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for data tree
 * @author Mohammed Ghannam
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_DATATREE_H__
#define __SCIP_TYPE_DATATREE_H__

#ifdef __cplusplus
extern "C" {
#endif

/** type of values stored in a SCIP_DATATREE */
enum SCIP_Datatree_Valuetype
{
    SCIP_DATATREE_BOOL,     /**< a SCIP_Bool value */
    SCIP_DATATREE_LONG,     /**< a SCIP_Longint integer value */
    SCIP_DATATREE_REAL,     /**< a SCIP_Real floating point value */
    SCIP_DATATREE_STRING,   /**< a C string */
    SCIP_DATATREE_BOOLARRAY,/**< an array of SCIP_Bool values */
    SCIP_DATATREE_LONGARRAY,/**< an array of SCIP_Longint values */
    SCIP_DATATREE_REALARRAY,/**< an array of SCIP_Real values */
    SCIP_DATATREE_STRINGARRAY,/**< an array of C strings */
    SCIP_DATATREE_DATATREE, /**< a SCIP_DATATREE object */
};

/** type of values stored in a SCIP_DATATREE */
typedef enum SCIP_Datatree_Valuetype SCIP_DATATREE_VALUETYPE;

/** generic hierarchical data storage */
typedef struct SCIP_Datatree SCIP_DATATREE;

#ifdef __cplusplus
}
#endif

#endif
