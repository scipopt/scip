/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2004 Tobias Achterberg                              */
/*                  2002-2004 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_branch.h,v 1.3 2004/02/04 17:27:36 bzfpfend Exp $"

/**@file   pub_branch.h
 * @brief  public methods for branching rules
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_BRANCH_H__
#define __PUB_BRANCH_H__


#include "def.h"
#include "type_misc.h"
#include "type_branch.h"



/** compares two branching rules w. r. to their priority */
extern
DECL_SORTPTRCOMP(SCIPbranchruleComp);

/** gets user data of branching rule */
extern
BRANCHRULEDATA* SCIPbranchruleGetData(
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** sets user data of branching rule; user has to free old data in advance! */
extern
void SCIPbranchruleSetData(
   BRANCHRULE*      branchrule,         /**< branching rule */
   BRANCHRULEDATA*  branchruledata      /**< new branching rule user data */
   );

/** gets name of branching rule */
extern
const char* SCIPbranchruleGetName(
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets description of branching rule */
extern
const char* SCIPbranchruleGetDesc(
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** gets priority of branching rule */
extern
int SCIPbranchruleGetPriority(
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** is branching rule initialized? */
extern
Bool SCIPbranchruleIsInitialized(
   BRANCHRULE*      branchrule          /**< branching rule */
   );


#endif
