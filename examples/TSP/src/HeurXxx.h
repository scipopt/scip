/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2007 Tobias Achterberg                              */
/*                                                                           */
/*                  2002-2007 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur_xxx.h
 * @brief  xxx primal heuristic
 * @author
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __HEURXXX_H__
#define __HEURXXX_H__


#include "objscip/objheur.h"
#include "GomoryHuTree.h"

namespace tsp
{

   /** C++ primal heuristic */
   class HeurXxx : public scip::ObjHeur
   {
      
   public:
      /** default constructor */
      HeurXxx(
      )
         : ObjHeur("xxx", "", '?',0, 1, 0, -1, SCIP_HEURTIMING_AFTERNODE)
      {
      }

      /** destructor */
      virtual ~HeurXxx()
      {
      }


      /** destructor of primal heuristic to free user data (called when SCIP is exiting) */
      virtual SCIP_RETCODE scip_free(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** initialization method of primal heuristic (called after problem was transformed) */
      virtual SCIP_RETCODE scip_init(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** deinitialization method of primal heuristic (called before transformed problem is freed) */
      virtual SCIP_RETCODE scip_exit(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
       *
       *  This method is called when the presolving was finished and the branch and bound process is about to begin.
       *  The primal heuristic may use this call to initialize its branch and bound specific data.
       *
       */
      virtual SCIP_RETCODE scip_initsol(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
       *
       *  This method is called before the branch and bound process is freed.
       *  The primal heuristic should use this call to clean up its branch and bound data.
       */
      virtual SCIP_RETCODE scip_exitsol(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur                /**< the primal heuristic itself */
         );
   
      /** execution method of primal heuristic
       *
       *  Searches for feasible primal solutions. The method is called in the node processing loop.
       *
       *  possible return values for *result:
       *  - SCIP_FOUNDSOL   : at least one feasible primal solution was found
       *  - SCIP_DIDNOTFIND : the heuristic searched, but did not find a feasible solution
       *  - SCIP_DIDNOTRUN  : the heuristic was skipped
       *  - SCIP_DELAYED    : the heuristic was skipped, but should be called again as soon as possible, disregarding
       *                      its frequency
       */
      virtual SCIP_RETCODE scip_exec(
         SCIP*              scip,               /**< SCIP data structure */
         SCIP_HEUR*         heur,               /**< the primal heuristic itself */
         SCIP_RESULT*       result              /**< pointer to store the result of the heuristic call */
         );
   };

}
#endif
