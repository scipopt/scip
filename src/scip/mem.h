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
#pragma ident "@(#) $Id: mem.h,v 1.18 2005/05/31 17:20:16 bzfpfend Exp $"

/**@file   mem.h
 * @brief  methods for block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MEM_H__
#define __MEM_H__


#include "scip/def.h"
#include "scip/type_retcode.h"
#include "scip/type_mem.h"
#include "scip/struct_mem.h"


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
