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
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: pub_disp.h,v 1.8 2005/07/15 17:20:14 bzfpfend Exp $"

/**@file   pub_disp.h
 * @brief  public methods for displaying runtime statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_PUB_DISP_H__
#define __SCIP_PUB_DISP_H__


#include <stdio.h>

#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_disp.h"



/** gets user data of display column */
extern
DISPDATA* SCIPdispGetData(
   DISP*            disp                /**< display column */
   );

/** sets user data of display column; user has to free old data in advance! */
extern
void SCIPdispSetData(
   DISP*            disp,               /**< display column */
   DISPDATA*        dispdata            /**< new display column user data */
   );

/** gets name of display column */
extern
const char* SCIPdispGetName(
   DISP*            disp                /**< display column */
   );

/** gets description of display column */
extern
const char* SCIPdispGetDesc(
   DISP*            disp                /**< display column */
   );

/** gets head line of display column */
extern
const char* SCIPdispGetHeader(
   DISP*            disp                /**< display column */
   );

/** gets width of display column */
extern
int SCIPdispGetWidth(
   DISP*            disp                /**< display column */
   );

/** gets priority of display column */
extern
int SCIPdispGetPriority(
   DISP*            disp                /**< display column */
   );

/** gets position of display column */
extern
int SCIPdispGetPosition(
   DISP*            disp                /**< display column */
   );

/** gets status of display column */
extern
DISPSTATUS SCIPdispGetStatus(
   DISP*            disp                /**< display column */
   );

/** is display column initialized? */
extern
Bool SCIPdispIsInitialized(
   DISP*            disp                /**< display column */
   );

/** displays a long integer in decimal form fitting in a given width */
extern
void SCIPdispLongint(
   FILE*            file,               /**< output stream */
   Longint          val,                /**< value to display */
   int              width               /**< width to fit into */
   );

/** displays an integer in decimal form fitting in a given width */
extern
void SCIPdispInt(
   FILE*            file,               /**< output stream */
   int              val,                /**< value to display */
   int              width               /**< width to fit into */
   );

/** displays a time value fitting in a given width */
extern
void SCIPdispTime(
   FILE*            file,               /**< output stream */
   Real             val,                /**< value in seconds to display */
   int              width               /**< width to fit into */
   );


#endif
