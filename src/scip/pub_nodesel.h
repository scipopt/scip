/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2005 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2005 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_nodesel.h,v 1.4 2005/01/21 09:17:03 bzfpfend Exp $"

/**@file   pub_nodesel.h
 * @brief  public methods for node selectors
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_NODESEL_H__
#define __PUB_NODESEL_H__


#include "def.h"
#include "type_nodesel.h"



/** gets name of node selector */
extern
const char* SCIPnodeselGetName(
   NODESEL*         nodesel             /**< node selector */
   );

/** gets description of node selector */
extern
const char* SCIPnodeselGetDesc(
   NODESEL*         nodesel             /**< node selector */
   );

/** gets priority of node selector in standard mode */
extern
int SCIPnodeselGetStdPriority(
   NODESEL*         nodesel             /**< node selector */
   );

/** gets priority of node selector in memory saving mode */
extern
int SCIPnodeselGetMemsavePriority(
   NODESEL*         nodesel             /**< node selector */
   );

/** gets user data of node selector */
extern
NODESELDATA* SCIPnodeselGetData(
   NODESEL*         nodesel             /**< node selector */
   );

/** sets user data of node selector; user has to free old data in advance! */
extern
void SCIPnodeselSetData(
   NODESEL*         nodesel,            /**< node selector */
   NODESELDATA*     nodeseldata         /**< new node selector user data */
   );

/** is node selector initialized? */
extern
Bool SCIPnodeselIsInitialized(
   NODESEL*         nodesel             /**< node selector */
   );


#endif
