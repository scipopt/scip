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
/*  SCIP is distributed under the terms of the SCIP Academic License.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: struct_mem.h,v 1.2 2005/02/14 13:35:52 bzfpfend Exp $"

/**@file   struct_mem.h
 * @brief  datastructures for block memory pools and memory buffers
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_MEM_H__
#define __STRUCT_MEM_H__


#include "scip/def.h"
#include "scip/memory.h"
#include "scip/type_mem.h"


/** various block memory buffers */
struct Mem
{
   BLKMEM*          setmem;             /**< memory blocks for parameter settings */
   BLKMEM*          probmem;            /**< memory blocks for original problem */
   BLKMEM*          solvemem;           /**< memory blocks for solution process: preprocessing, bab-tree, ... */
};


#endif
