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

/**@file   compute_symmetry_none.cpp
 * @brief  interface for no symmetry computations
 * @author Marc Pfetsch
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "compute_symmetry.h"

/** return whether symmetry can be computed */
SCIP_Bool SYMcanComputeSymmetry(void)
{
   return FALSE;
}

/** return name of external program used to compute generators */
const char* SYMsymmetryGetName(void)
{
   return "none";
}

/** return description of external program used to compute generators */
const char* SYMsymmetryGetDesc(void)
{
   return "";
}

/** return name of additional external program used for computing symmetries */
const char* SYMsymmetryGetAddName(void)
{
   return "";
}

/** return description of additional external program used to compute symmetries */
const char* SYMsymmetryGetAddDesc(void)
{
   return "";
}

/** compute generators of symmetry group */ /*lint -e{715}*/
SCIP_RETCODE SYMcomputeSymmetryGenerators(
   SCIP*                 scip,               /**< SCIP pointer */
   int                   maxgenerators,      /**< maximal number of generators constructed (= 0 if unlimited) */
   SYM_GRAPH*            graph,              /**< symmetry detection graph */
   int*                  nperms,             /**< pointer to store number of permutations */
   int*                  nmaxperms,          /**< pointer to store maximal number of permutations (needed for freeing storage) */
   int***                perms,              /**< pointer to store permutation generators as (nperms x npermvars) matrix */
   SCIP_Real*            log10groupsize,     /**< pointer to store size of group */
   SCIP_Real*            symcodetime         /**< pointer to store the time for symmetry code */
   )
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( graph != NULL );
   assert( nperms != NULL );
   assert( nmaxperms != NULL );
   assert( perms != NULL );
   assert( log10groupsize != NULL );
   assert( symcodetime != NULL );

   /* init */
   *nperms = 0;
   *nmaxperms = 0;
   *perms = NULL;
   *log10groupsize = 0;
   *symcodetime = 0.0;

   return SCIP_OKAY;
}

/** returns whether two given graphs are identical */
SCIP_Bool SYMcheckGraphsAreIdentical(
   SCIP*                 scip,               /**< SCIP pointer */
   SYM_SYMTYPE           symtype,            /**< type of symmetries to be checked */
   SYM_GRAPH*            G1,                 /**< first graph */
   SYM_GRAPH*            G2                  /**< second graph */
   )
{
   return FALSE;
}
