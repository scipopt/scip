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
#pragma ident "@(#) $Id: stat.h,v 1.28 2004/01/07 13:14:15 bzfpfend Exp $"

/**@file   stat.h
 * @brief  internal methods for problem statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STAT_H__
#define __STAT_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_mem.h"

#include "struct_stat.h"



/** creates problem statistics data */
extern
RETCODE SCIPstatCreate(
   STAT**           stat,               /**< pointer to problem statistics data */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set                 /**< global SCIP settings */
   );

/** frees problem statistics data */
extern
RETCODE SCIPstatFree(
   STAT**           stat,               /**< pointer to problem statistics data */
   MEMHDR*          memhdr              /**< block memory */
   );

/** marks statistics to be able to reset them when solving process is freed */
extern
void SCIPstatMark(
   STAT*            stat                /**< problem statistics data */
   );

/** reset statistics to the data before solving started */
extern
void SCIPstatReset(
   STAT*            stat                /**< problem statistics data */
   );

/** resets display statistics, such that a new header line is displayed before the next display line */
extern
void SCIPstatResetDisplay(
   STAT*            stat                /**< problem statistics data */
   );

/** depending on the current memory usage, switches mode flag to standard or memory saving mode */
extern
void SCIPstatUpdateMemsaveMode(
   STAT*            stat,               /**< problem statistics data */
   SET*             set,                /**< global SCIP settings */
   MEM*             mem                 /**< block memory pools */
   );

#endif
