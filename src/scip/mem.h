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
#pragma ident "@(#) $Id: mem.h,v 1.13 2004/02/04 17:27:27 bzfpfend Exp $"

/**@file   mem.h
 * @brief  methods and datastructures for block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MEM_H__
#define __MEM_H__


#include "def.h"
#include "memory.h"
#include "type_retcode.h"
#include "type_mem.h"


/** various block memory buffers */
struct Mem
{
   MEMHDR*          setmem;             /**< memory blocks for parameter settings */
   MEMHDR*          probmem;            /**< memory blocks for original problem */
   MEMHDR*          solvemem;           /**< memory blocks for solution process: preprocessing, bab-tree, ... */
};



/** creates block memory structures */
extern
RETCODE SCIPmemCreate(
   MEM**            mem                 /**< pointer to block memory structure */
   );

/** frees block memory structures */
extern
RETCODE SCIPmemFree(
   MEM**            mem                 /**< pointer to block memory structure */
   );

/** returns the total number of bytes used in block memory */
extern
Longint SCIPmemGetUsed(
   MEM*             mem                 /**< pointer to block memory structure */
   );


#endif
