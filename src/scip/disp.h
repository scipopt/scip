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
#pragma ident "@(#) $Id: disp.h,v 1.15 2003/11/21 10:35:35 bzfpfend Exp $"

/**@file   disp.h
 * @brief  methods and datastructures for displaying runtime statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DISP_H__
#define __DISP_H__

/** display activation status of display column */
enum DispStatus
{
   SCIP_DISPSTATUS_OFF  = 0,            /**< display column is not displayed */
   SCIP_DISPSTATUS_AUTO = 1,            /**< display column is switched on and off automatically */
   SCIP_DISPSTATUS_ON   = 2             /**< display column is displayed */
};
typedef enum DispStatus DISPSTATUS;

typedef struct Disp DISP;               /**< display column data structure */
typedef struct DispData DISPDATA;       /**< display column specific data */


/** destructor of display column to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define DECL_DISPFREE(x) RETCODE x (SCIP* scip, DISP* disp)

/** initialization method of display column (called when problem solving starts)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define DECL_DISPINIT(x) RETCODE x (SCIP* scip, DISP* disp)

/** deinitialization method of display column (called when problem solving exits)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define DECL_DISPEXIT(x) RETCODE x (SCIP* scip, DISP* disp)

/** output method of display column to output file stream 'file'
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 *  - file            : file stream for output
 */
#define DECL_DISPOUTPUT(x) RETCODE x (SCIP* scip, DISP* disp, FILE* file)



#include "scip.h"
#include "def.h"
#include "retcode.h"
#include "set.h"
#include "stat.h"
#include "tree.h"
#include "lp.h"


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

/** gets position of display column */
extern
int SCIPdispGetPosition(
   DISP*            disp                /**< display column */
   );

/** is display column initialized? */
extern
Bool SCIPdispIsInitialized(
   DISP*            disp                /**< display column */
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

/** displays an integer in decimal form fitting in a given width */
extern
void SCIPdispDecimal(
   FILE*            file,               /**< output stream */
   Longint          val,                /**< value to display */
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
