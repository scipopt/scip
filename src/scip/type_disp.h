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
#pragma ident "@(#) $Id: type_disp.h,v 1.3 2004/04/27 15:50:06 bzfpfend Exp $"

/**@file   type_disp.h
 * @brief  type definitions for displaying runtime statistics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __TYPE_DISP_H__
#define __TYPE_DISP_H__

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

/** initialization method of display column (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - disp            : the display column itself
 */
#define DECL_DISPINIT(x) RETCODE x (SCIP* scip, DISP* disp)

/** deinitialization method of display column (called before transformed problem is freed)
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


#include <stdio.h>

#include "def.h"
#include "type_retcode.h"
#include "type_scip.h"


#endif
