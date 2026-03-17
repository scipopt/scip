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

/**@file   pub_sym.h
 * @ingroup PUBLICCOREAPI
 * @brief  public methods for symmetry handlers
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_SYM_H__
#define __SCIP_PUB_SYM_H__


#include "scip/def.h"
#include "scip/type_misc.h"
#include "scip/type_sym.h"

#ifdef __cplusplus
extern "C" {
#endif

/** gets data of symmetry component */
SCIP_EXPORT
SCIP_SYMCOMPDATA* SCIPsymcompGetData(
   SCIP_SYMCOMP*         symcomp             /**< symmetry component */
   );

/** gets symmetry handler of symmetry component */
SCIP_EXPORT
SCIP_SYMHDLR* SCIPsymcompGetHdlr(
   SCIP_SYMCOMP*         symcomp             /**< symmetry component */
   );

/** gets name of symmetry component */
SCIP_EXPORT
const char* SCIPsymcompGetName(
   SCIP_SYMCOMP*         symcomp             /**< symmetry component */
   );

/** gets name of symmetry handler */
SCIP_EXPORT
const char* SCIPsymhdlrGetName(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets user data of symmetry handler */
SCIP_EXPORT
SCIP_SYMHDLRDATA* SCIPsymhdlrGetData(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets priority of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets presolving priority of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetPresolPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets description of symmetry handler */
SCIP_EXPORT
const char* SCIPsymhdlrGetDesc(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** does the symmetry handler perform presolving? */
SCIP_EXPORT
SCIP_Bool SCIPsymhdlrDoesPresolve(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets time in seconds used for setting up this symmetry handler for new stages */
SCIP_EXPORT
SCIP_Real SCIPsymhdlrGetSetupTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets time in seconds used in this symmetry handler for presolving */
SCIP_EXPORT
SCIP_Real SCIPsymhdlrGetPresolTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of times the symmetry handler was called in presolving and tried to find reductions */
SCIP_EXPORT
SCIP_Longint SCIPsymhdlrGetNPresolCalls(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of variables fixed during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNFixedVars(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of variables aggregated during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNAggrVars(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of variable types changed during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNChgVarTypes(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of bounds changed during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNChgBds(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of holes added to domains of variables during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNAddHoles(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of constraints deleted during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNDelConss(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of constraints added during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNAddConss(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of coefficients changed during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNChgCoefs(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets priority of separation method of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetSepaPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets priority of propagation method of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetPropPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets number of constraint sides changed during presolving of symmetry handler */
SCIP_EXPORT
int SCIPsymhdlrGetNChgSides(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets the total number of times, the propagator of the symmetry handler was called */
SCIP_EXPORT
SCIP_Longint SCIPsymhdlrGetNPropCalls(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets the total number of times, the propagator of the symmetry handler was called for resolving a propagation */
SCIP_EXPORT
SCIP_Longint SCIPsymhdlrGetNRespropCalls(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets total number of times, the symmetry handler's propagator detected a cutoff */
SCIP_EXPORT
SCIP_Longint SCIPsymhdlrGetNCutoffs(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets total number of domain reductions found by this symmetry handler's propagator */
SCIP_EXPORT
SCIP_Longint SCIPsymhdlrGetNDomredsFound(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** returns the timing mask of the symmetry handler's propagator */
SCIP_EXPORT
SCIP_PROPTIMING SCIPsymhdlrPropGetTimingmask(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets priority of symmetry handler's propagator */
SCIP_EXPORT
int SCIPsymhdlrPropGetPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** was symmetry handler's propagator delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPsymhdlrPropWasDelayed(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets priority of symmetry handlers' separator */
SCIP_EXPORT
int SCIPsymhdlrSepaGetPriority(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** was symmetry handler's separator delayed at the last call? */
SCIP_EXPORT
SCIP_Bool SCIPsymhdlrSepaWasLPDelayed(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets the total number of times, the separator of the symmetry handler was called */
SCIP_EXPORT
SCIP_Longint SCIPsymhdlrGetNSepaCalls(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets the total number of times the symmetry handler's separator was called */
SCIP_EXPORT
SCIP_Longint SCIPsymhdlrGetNCutsFound(
      SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets the total number of cutting planes added by the symmetry handler's separator */
SCIP_EXPORT
SCIP_Longint SCIPsymhdlrGetNCutsAdded(
      SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets time in seconds used in the symmetry handler's separator */
SCIP_EXPORT
SCIP_Real SCIPsymhdlrGetSepaTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets time in seconds used in the symmetry handler's propagator */
SCIP_EXPORT
SCIP_Real SCIPsymhdlrGetPropTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets time in seconds used in the symmetry handler's propagator during strong branching */
SCIP_EXPORT
SCIP_Real SCIPsymhdlrGetStrongBranchPropTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** gets time in seconds used in th symmetry handler's propagator for resolve propagation */
SCIP_EXPORT
SCIP_Real SCIPsymhdlrGetRespropTime(
   SCIP_SYMHDLR*         symhdlr             /**< symmetry handler */
   );

/** compares two symmetry handlers with respect to their try-add priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompTryadd);

/** compares two symmetry handlers with respect to their separation priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompSepa);

/** compares two symmetry handlers with respect to their propagation priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompProp);

/** compares two symmetry handlers w.r.t. their presolving priority */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompPresol);

/** comparison method for sorting symmetry handlers w.r.t. to their name */
SCIP_EXPORT
SCIP_DECL_SORTPTRCOMP(SCIPsymhdlrCompName);

/** creates new operator node type (used for symmetry detection) and returns its representation
 *
 *  If the operator node already exists, the function terminates with SCIP_INVALIDDATA.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPcreateSymOpNodeType(
   SCIP*                 scip,               /**< SCIP pointer */
   const char*           opnodename,         /**< name of new operator node type */
   int*                  nodetype            /**< pointer to store the new node type */
   );

/** returns representation of an operator node type.
 *
 *  If the node type does not already exist, a new node type will be created.
 */
SCIP_EXPORT
SCIP_RETCODE SCIPgetSymOpNodeType(
   SCIP*                 scip,               /**< SCIP pointer */
   const char*           opnodename,         /**< name of new operator node type */
   int*                  nodetype            /**< pointer to store the node type */
   );

#ifdef __cplusplus
}
#endif

#endif
