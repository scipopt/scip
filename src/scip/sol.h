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
 * @brief  datastructures and methods for storing primal CIP solutions
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SOL_H__
#define __SOL_H__


typedef struct Sol SOL;                 /**< primal CIP solution */


#include "def.h"
#include "retcode.h"
#include "memory.h"
#include "stat.h"
#include "var.h"
#include "heur.h"
#include "prob.h"


/** primal CIP solution (variables with index < firstindex or index >= firstindex+nvals have solution value 0.0) */
struct Sol
{
   HEUR*            heur;               /**< heuristic that found the solution (or NULL if it's an LP solution) */
   VAR**            vars;               /**< variables in the index range */
   Real*            vals;               /**< solution values for variables in the index range */
   Real             obj;                /**< objective value of solution */
   int              nvals;              /**< number of values in the index range of solution */
   int              valssize;           /**< size of vars and vals array */
   int              firstindex;         /**< first index of the index range */
   int              nuses;              /**< number of times, this solution is referenced */
   int              nodenum;            /**< node number, where this solution was found */
};


extern
RETCODE SCIPsolCreate(                  /**< creates and captures primal CIP solution */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   STAT*            stat,               /**< problem statistics data */
   HEUR*            heur                /**< heuristic that found the solution (or NULL if it's an LP solution) */
   );

extern
RETCODE SCIPsolCreateLPSol(             /**< copys LP solution to primal CIP solution, and captures solution */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPsolCreatePseudoSol(         /**< copys pseudo solution to primal CIP solution, and captures solution */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   STAT*            stat,               /**< problem statistics data */
   PROB*            prob                /**< problem data */
   );

extern
RETCODE SCIPsolFree(                    /**< frees primal CIP solution */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
void SCIPsolCapture(                    /**< increases usage counter of primal CIP solution */
   SOL*             sol                 /**< primal CIP solution */
   );

extern
RETCODE SCIPsolRelease(                 /**< decreases usage counter of primal CIP solution, frees memory if necessary */
   SOL**            sol,                /**< pointer to primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPsolClear(                   /**< clears primal CIP solution */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   LP*              lp                  /**< actual LP data */
   );

extern
RETCODE SCIPsolSetVal(                  /**< sets value of variable in primal CIP solution */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to solution */
   Real             val                 /**< solution value of variable */
   );

extern
RETCODE SCIPsolIncVal(                  /**< increases value of variable in primal CIP solution */
   SOL*             sol,                /**< primal CIP solution */
   MEMHDR*          memhdr,             /**< block memory */
   const SET*       set,                /**< global SCIP settings */
   VAR*             var,                /**< variable to add to solution */
   Real             incval              /**< increment for solution value of variable */
   );

extern
Real SCIPsolGetVal(                     /**< returns value of variable in primal CIP solution */
   SOL*             sol,                /**< primal CIP solution */
   VAR*             var                 /**< variable to get value for */
   );

extern
void SCIPsolPrint(                      /**< outputs non-zero elements of solution to file stream */
   SOL*             sol,                /**< primal CIP solution */
   const SET*       set,                /**< global SCIP settings */
   FILE*            file                /**< output file (or NULL for standard output) */
   );

extern
int SCIPsolGetNodenum(                  /**< gets node number, where this solution was found */
   SOL*             sol                 /**< primal CIP solution */
   );

extern
HEUR* SCIPsolGetHeur(                   /**< gets heuristic, that found this solution (or NULL if it's an LP solution) */
   SOL*             sol                 /**< primal CIP solution */
   );


#endif
