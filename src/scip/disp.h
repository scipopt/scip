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
#pragma ident "@(#) $Id: disp.h,v 1.19 2004/03/15 15:40:18 bzfpfend Exp $"

/**@file   disp.h
 * @brief  internal methods for displaying runtime statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DISP_H__
#define __DISP_H__


#include <stdio.h>

#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_set.h"
#include "type_stat.h"
#include "type_disp.h"
#include "pub_disp.h"



/** parameter change information method to autoselect display columns again */
extern
DECL_PARAMCHGD(SCIPparamChgdDispActive);

/** creates a display column */
extern
RETCODE SCIPdispCreate(
   DISP**           disp,               /**< pointer to store display column */
   SET*             set,                /**< global SCIP settings */
   MEMHDR*          memhdr,             /**< block memory for parameter settings */
   const char*      name,               /**< name of display column */
   const char*      desc,               /**< description of display column */
   const char*      header,             /**< head line of display column */
   DISPSTATUS       dispstatus,         /**< display activation status of display column */
   DECL_DISPFREE    ((*dispfree)),      /**< destructor of display column */
   DECL_DISPINIT    ((*dispinit)),      /**< initialize display column */
   DECL_DISPEXIT    ((*dispexit)),      /**< deinitialize display column */
   DECL_DISPOUTPUT  ((*dispoutput)),    /**< output method */
   DISPDATA*        dispdata,           /**< display column data */
   int              width,              /**< width of display column (no. of chars used) */
   int              priority,           /**< priority of display column */
   int              position,           /**< relative position of display column */
   Bool             stripline           /**< should the column be separated with a line from its right neighbour? */
   );

/** frees memory of display column */
extern
RETCODE SCIPdispFree(
   DISP**           disp,               /**< pointer to display column data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** initializes display column */
extern
RETCODE SCIPdispInit(
   DISP*            disp,               /**< display column */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** deinitializes display column */
extern
RETCODE SCIPdispExit(
   DISP*            disp,               /**< display column */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** output display column to screen */
extern
RETCODE SCIPdispOutput(
   DISP*            disp,               /**< display column */
   const SET*       set                 /**< global SCIP settings */
   );

/** prints one line of output with the active display columns */
extern
RETCODE SCIPdispPrintLine(
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   Bool             forcedisplay        /**< should the line be printed without regarding frequency? */
   );

/** activates all display lines fitting in the display w.r. to priority */
extern
RETCODE SCIPdispAutoActivate(
   const SET*       set                 /**< global SCIP settings */
   );


#endif
