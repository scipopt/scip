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

/**@file   mem.h
 * @brief  block memory pools
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __MEM_H__
#define __MEM_H__

#include "memory.h"

struct Mem                              /**< various block memory buffers */
{
   MEMHDR*          treemem;            /**< ptr to memory blocks for the tree */
   MEMHDR*          statemem;           /**< ptr to memory blocks for LP states */
   MEMHDR*          lpmem;              /**< ptr to memory blocks for LP data */
   MEMHDR*          primalmem;          /**< ptr to memory blocks for primal solutions */
   MEMHDR*          tempmem;            /**< ptr to memory blocks for short living objects */
};
typedef struct Mem MEM;


#endif
