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

/**@file   sol.h
 * @brief  datastructures and methods for storing primal IP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SOL_H__
#define __SOL_H__


typedef struct Sol SOL;                 /**< primal IP solution */


#include "def.h"
#include "retcode.h"
#include "memory.h"
#include "var.h"


/** primal IP solution (variables with index < firstindex or index >= firstindex+nvals have solution value 0.0) */
struct Sol
{
   VAR**            vars;               /**< variables in the index range */
   Real*            vals;               /**< solution values for variables in the index range */
   Real             obj;                /**< objective value of solution */
   int              nvals;              /**< number of values in the index range of solution */
   int              valssize;           /**< size of vars and vals array */
   int              firstindex;         /**< first index of the index range */
};


extern
RETCODE SCIPsolCreate(                  /**< creates primal IP solution */
   SOL**            sol,                /**< pointer to primal IP solution */
   MEMHDR*          memhdr              /**< block memory */
   );

extern
void SCIPsolFree(                       /**< frees primal IP solution */
   SOL**            sol,                /**< pointer to primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's an original variable) */
   );

extern
void SCIPsolClear(                      /**< clears primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data (or NULL, if it's an original variable) */
   );

extern
RETCODE SCIPsolSetVal(                  /**< sets value of variable in primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   );

extern
RETCODE SCIPsolIncVal(                  /**< increases value of variable in primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to solution */
   Real             incval              /**< increment for solution value of variable */
   );

extern
RETCODE SCIPsolCopyLPSol(               /**< copys LP solution to primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
Real SCIPsolGetVal(                     /**< returns value of variable in primal IP solution */
   SOL*             sol,                /**< primal IP solution */
   VAR*             var                 /**< variable to get value for */
   );

extern
void SCIPsolPrint(                      /**< outputs non-zero elements of solution to file stream */
   SOL*             sol,                /**< primal IP solution */
   const SET*       set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

#endif
