/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2002 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                            Alexander Martin                               */
/*                  2002-2002 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   disp.h
 * @brief  datastructures and methods for displaying runtime statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __DISP_H__
#define __DISP_H__


typedef struct Disp DISP;               /**< display column data structure */
typedef struct DispData DISPDATA;       /**< display column specific data */


/** destructor of display column to free user data (called when SCIP is exiting)
 *  possible return values:
 *    SCIP_OKAY       : normal termination
 *    neg. values     : error codes
 */
#define DECL_DISPFREE(x) RETCODE x (DISP* disp, SCIP* scip)

/** initialization method of display column (called at problem creation)
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_DISPINIT(x) RETCODE x (DISP* disp, SCIP* scip)

/** deinitialization method of display column (called at problem destruction)
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_DISPEXIT(x) RETCODE x (DISP* disp, SCIP* scip)

/** output method of display column to output file stream 'file'
 *  possible return values:
 *    SCIP_OKAY   : normal termination
 *    neg. values : error codes
 */
#define DECL_DISPOUTP(x) RETCODE x (DISP* disp, SCIP* scip, FILE* file)



#include "scip.h"
#include "retcode.h"
#include "set.h"
#include "stat.h"
#include "tree.h"
#include "lp.h"


extern
RETCODE SCIPdispCreate(                 /**< creates a display column */
   DISP**           disp,               /**< pointer to store display column */
   const char*      name,               /**< name of display column */
   const char*      desc,               /**< description of display column */
   const char*      header,             /**< head line of display column */
   DECL_DISPFREE((*dispfree)),          /**< destructor of display column */
   DECL_DISPINIT((*dispinit)),          /**< initialise display column */
   DECL_DISPEXIT((*dispexit)),          /**< deinitialise display column */
   DECL_DISPOUTP((*dispoutp)),          /**< output method */
   DISPDATA*        dispdata,           /**< display column data */
   int              width,              /**< width of display column (no. of chars used) */
   int              priority,           /**< priority of display column */
   int              position,           /**< relative position of display column */
   Bool             stripline           /**< should the column be separated with a line from its right neighbour? */
   );

extern
RETCODE SCIPdispFree(                   /**< frees memory of display column */
   DISP**           disp,               /**< pointer to display column data structure */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPdispInit(                   /**< initializes display column */
   DISP*            disp,               /**< display column */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPdispExit(                   /**< deinitializes display column */
   DISP*            disp,               /**< display column */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPdispOutput(                 /**< output display column to screen */
   DISP*            disp,               /**< display column */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
const char* SCIPdispGetName(            /**< gets name of display column */
   DISP*            disp                /**< display column */
   );

extern
DISPDATA* SCIPdispGetData(              /**< gets user data of display column */
   DISP*            disp                /**< display column */
   );

extern
void SCIPdispSetData(                   /**< sets user data of display column; user has to free old data in advance! */
   DISP*            disp,               /**< display column */
   DISPDATA*        dispdata            /**< new display column user data */
   );

extern
int SCIPdispGetPosition(                /**< gets position of display column */
   DISP*            disp                /**< display column */
   );

extern
Bool SCIPdispIsInitialized(             /**< is display column initialized? */
   DISP*            disp                /**< display column */
   );

extern
RETCODE SCIPdispPrintLine(              /**< prints one line of output with the active display columns */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   Bool             forcedisplay        /**< should the line be printed without regarding frequency? */
   );

extern
RETCODE SCIPdispAutoActivate(           /**< activates all display lines fitting in the display w.r. to priority */
   const SET*       set                 /**< global SCIP settings */
   );


#endif
