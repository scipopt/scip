/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_presol.h,v 1.1 2003/12/01 14:41:30 bzfpfend Exp $"

/**@file   pub_presol.h
 * @brief  public methods for presolvers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PUB_PRESOL_H__
#define __PUB_PRESOL_H__


#include "def.h"
#include "type_presol.h"



/** gets user data of presolver */
extern
PRESOLDATA* SCIPpresolGetData(
   PRESOL*          presol              /**< presolver */
   );

/** sets user data of presolver; user has to free old data in advance! */
extern
void SCIPpresolSetData(
   PRESOL*          presol,             /**< presolver */
   PRESOLDATA*      presoldata          /**< new presolver user data */
   );

/** gets name of presolver */
extern
const char* SCIPpresolGetName(
   PRESOL*          presol              /**< presolver */
   );

/** gets description of presolver */
extern
const char* SCIPpresolGetDesc(
   PRESOL*          presol              /**< presolver */
   );

/** gets priority of presolver */
extern
int SCIPpresolGetPriority(
   PRESOL*          presol              /**< presolver */
   );

/** is presolver initialized? */
extern
Bool SCIPpresolIsInitialized(
   PRESOL*          presol              /**< presolver */
   );

/** gets time in seconds used in this presolver */
extern
Real SCIPpresolGetTime(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of variables fixed in presolver */
extern
int SCIPpresolGetNFixedVars(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of variables aggregated in presolver */
extern
int SCIPpresolGetNAggrVars(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of variable types changed in presolver */
extern
int SCIPpresolGetNVarTypes(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of bounds changed in presolver */
extern
int SCIPpresolGetNChgBds(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of holes added to domains of variables in presolver */
extern
int SCIPpresolGetNAddHoles(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of constraints deleted in presolver */
extern
int SCIPpresolGetNDelConss(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of constraints upgraded in presolver */
extern
int SCIPpresolGetNUpgdConss(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of coefficients changed in presolver */
extern
int SCIPpresolGetNChgCoefs(
   PRESOL*          presol              /**< presolver */
   );

/** gets number of constraint sides changed in presolver */
extern
int SCIPpresolGetNChgSides(
   PRESOL*          presol              /**< presolver */
   );


#endif
