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
#pragma ident "@(#) $Id: pub_relax.h,v 1.2 2004/11/19 14:45:12 bzfpfend Exp $"

/**@file   pub_relax.h
 * @brief  public methods for relaxators
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_RELAX_H__
#define __PUB_RELAX_H__


#include "def.h"
#include "type_misc.h"
#include "type_relax.h"



/** compares two relaxators w. r. to their priority */
extern
DECL_SORTPTRCOMP(SCIPrelaxComp);

/** gets user data of relaxator */
extern
RELAXDATA* SCIPrelaxGetData(
   RELAX*           relax               /**< relaxator */
   );

/** sets user data of relaxator; user has to free old data in advance! */
extern
void SCIPrelaxSetData(
   RELAX*           relax,              /**< relaxator */
   RELAXDATA*       relaxdata           /**< new relaxator user data */
   );

/** gets name of relaxator */
extern
const char* SCIPrelaxGetName(
   RELAX*           relax               /**< relaxator */
   );

/** gets description of relaxator */
extern
const char* SCIPrelaxGetDesc(
   RELAX*           relax               /**< relaxator */
   );

/** gets priority of relaxator */
extern
int SCIPrelaxGetPriority(
   RELAX*           relax               /**< relaxator */
   );

/** gets frequency of relaxator */
extern
int SCIPrelaxGetFreq(
   RELAX*           relax               /**< relaxator */
   );

/** gets time in seconds used in this relaxator */
extern
Real SCIPrelaxGetTime(
   RELAX*           relax               /**< relaxator */
   );

/** gets the total number of times, the relaxator was called */
extern
Longint SCIPrelaxGetNCalls(
   RELAX*           relax               /**< relaxator */
   );

/** is relaxator initialized? */
extern
Bool SCIPrelaxIsInitialized(
   RELAX*           relax               /**< relaxator */
   );

/** marks the current relaxation unsolved, s.t. the relaxator is called again in the next solving round */
extern
void SCIPrelaxMarkUnsolved(
   RELAX*           relax               /**< relaxator */
   );


#endif
