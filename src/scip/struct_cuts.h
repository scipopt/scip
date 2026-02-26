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

/**@file   struct_cuts.h
 * @ingroup PUBLICCOREAPI
 * @brief  struct definitions for cuts
 * @author Leona Gottwald
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_STRUCT_CUTS_H__
#define __SCIP_STRUCT_CUTS_H__

#include "scip/def.h"
#include "scip/dbldblarith.h"
#include "scip/type_lp.h"
#include "scip/type_cuts.h"

struct SCIP_AggrRow
{
   SCIP_Real*            vals;               /**< non-zero coefficients of the cut row */
   int*                  inds;               /**< problem indices of variables with a non-zero coefficient in the cut row */
   int*                  rowsinds;           /**< lpposition of rows that have been added to the cutrow */
   int*                  slacksign;          /**< slacksign of rows that have been added to the cutrow */
   SCIP_Real*            rowweights;         /**< weights of rows that have been added to the cutrow */
   QUAD_MEMBER(SCIP_Real rhs);               /**< right hand side of the cut row */
   int                   nnz;                /**< number of non-zeros in the cut row */
   int                   nrows;              /**< number of rows that have been added to the cutrow */
   int                   rowssize;           /**< size of the row and slacksign array */
   int                   rank;               /**< rank of the cut row */
   SCIP_Bool             local;              /**< is the cut row only valid locally? */
   SCIP_Longint          certificateline;    /**< proof index in certificate or SCIP_LONGINT_MAX */
};

/** parameters for cut generation methods */
struct SCIP_CutGenParams
{
   SCIP_Real             boundswitch;        /**< fraction of domain up to which lower bound is used in transformation */
   SCIP_Real             minfrac;            /**< minimal fractionality of rhs to produce cut for */
   SCIP_Real             maxfrac;            /**< maximal fractionality of rhs to produce cut for */
   int                   vartypeusevbds;     /**< variable types for which variable bound substitution is allowed */
   int                   maxtestdelta;       /**< maximum number of deltas to test (CMIR) */
   int*                  boundsfortrans;     /**< bounds that should be used for transformed variables (CMIR) */
   SCIP_BOUNDTYPE*       boundtypesfortrans; /**< type of bounds for transformed variables (CMIR) */
   SCIP_Bool             postprocess;        /**< apply post-processing step? */
   SCIP_Bool             allowlocal;         /**< should local information be allowed, resulting in a local cut? */
};

/** result of cut generation attempt */
struct SCIP_CutGenResult
{
   SCIP_Real*            cutcoefs;           /**< array of non-zero coefficients in the cut (pre-allocated) */
   int*                  cutinds;            /**< array of variable indices of non-zero coefficients (pre-allocated) */
   SCIP_CUTGENMETHOD     winningmethod;      /**< cut generation method which produced the best cut */
   SCIP_Real             cutefficacy;        /**< efficacy of the best cut */
   SCIP_Real             cutrhs;             /**< right hand side of the best cut */
   int                   cutnnz;             /**< number of non-zeros in the best cut */
   int                   cutrank;            /**< rank of the best cut */
   SCIP_Bool             cutislocal;         /**< is the best cut only valid locally? */
   SCIP_Bool             success;            /**< is a valid cut found? */
};

#endif
