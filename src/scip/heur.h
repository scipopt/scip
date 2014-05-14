/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   heur.h
 * @brief  internal methods for primal heuristics
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_HEUR_H__
#define __SCIP_HEUR_H__


#include "scip/def.h"
#include "blockmemshell/memory.h"
#include "scip/type_retcode.h"
#include "scip/type_result.h"
#include "scip/type_set.h"
#include "scip/type_primal.h"
#include "scip/type_heur.h"
#include "scip/pub_heur.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create a set of diving heuristic settings */
extern
SCIP_RETCODE SCIPdivesetCreate(
      SCIP_DIVESET**        diveset,            /**< pointer to the freshly created diveset */
      SCIP_HEUR*            heur,               /**< the heuristic to which this dive setting belongs */
      SCIP_SET*             set,                /**< global SCIP settings */
      SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
      BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
      SCIP_Real             minreldepth,        /**< minimal relative depth to start diving */
      SCIP_Real             maxreldepth,        /**< maximal relative depth to start diving */
      SCIP_Real             maxlpiterquot,      /**< maximal fraction of diving LP iterations compared to node LP iterations */
      SCIP_Real             maxdiveubquot,      /**< maximal quotient (curlowerbound - lowerbound)/(cutoffbound - lowerbound)
                                                 *   where diving is performed (0.0: no limit) */
      SCIP_Real             maxdiveavgquot,     /**< maximal quotient (curlowerbound - lowerbound)/(avglowerbound - lowerbound)
                                                 *   where diving is performed (0.0: no limit) */
      SCIP_Real             maxdiveubquotnosol, /**< maximal UBQUOT when no solution was found yet (0.0: no limit) */
      SCIP_Real             maxdiveavgquotnosol,/**< maximal AVGQUOT when no solution was found yet (0.0: no limit) */
      int                   maxlpiterofs,       /**< additional number of allowed LP iterations */
      SCIP_Bool             backtrack,          /**< use one level of backtracking if infeasibility is encountered? */
      SCIP_DECL_DIVESETGETSCORE((*divesetgetscore)), /**< get candidate score */
      SCIP_DECL_DIVESETCANDBRANCHDIR ((*divesetcandbranchdir)), /**< get preferred branching direction for a candidate */
      SCIP_DECL_DIVESETGETCANDS ((*divesetgetcands)), /**< allocate and get candidate variables for diving */
      SCIP_DECL_DIVESETFREECANDS ((*divesetfreecands))  /**< free previously allocated variables for diving */
   );

/** get the target depth fraction of the diving settings  */
extern
SCIP_Real SCIPdivesetGetTargetdepthfrac(
   SCIP_DIVESET*         diveset             /**< diving settings */
   );

/** set the target depth fraction of the diving settings  */
extern
void SCIPdivesetSetTargetdepthfrac(
   SCIP_DIVESET*         diveset,            /**< diving settings */
   SCIP_Real             newval              /**< new value for target depth frac */
   );

/** get branching candidates defined by the diveset. This will call the divesetGetCandsXyz callback
 *  of the diving settings. If the diving setting has no such callback, this is indicated by
 *  assigning the value -1. The calling method can then use SCIPgetLPBranchCands() instead
 */
SCIP_RETCODE SCIPdivesetGetCands(
   SCIP_DIVESET*         diveset,            /**< general diving settings */
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_VAR***           divecands,          /**< pointer to store the diving candidates */
   SCIP_Real**           divecandssol,       /**< pointer to store the LP solution values of the candidates */
   SCIP_Real**           divecandsfrac,      /**< pointer to store the fractionalities of the candidates */
   int*                  ndivecands          /**< pointer to store the number of candidates, or -1 */
   );

/** free memory allocated during SCIPdivesetGetCands() */
extern
SCIP_RETCODE SCIPdivesetFreeCands(
   SCIP_DIVESET*        diveset,             /**< general diving settings */
   SCIP_SET*            set,                 /**< SCIP settings */
   SCIP_VAR***          divecands,           /**< pointer to the diving candidates */
   SCIP_Real**          divecandssol,        /**< pointer to LP solution values of the candidates */
   SCIP_Real**          divecandsfrac,       /**< pointer to the fractionalities of the candidates */
   int                  ndivecands           /**< the number of candidates that were allocated */
   );

/** request the (preferred) branching direction of a candidate variable */
extern
SCIP_BRANCHDIR SCIPdivesetCandBranchdir(
   SCIP_DIVESET*         diveset,            /**< general diving settings */
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_VAR*             divecand,           /**< the candidate for which the branching direction is requested */
   SCIP_Real             divecandsol,        /**< LP solution value of the candidate */
   SCIP_Real             divecandfrac        /**< fractionality of the candidate */
   );

/** returns the candidate score for diving */
extern
SCIP_Real SCIPdivesetGetScore(
   SCIP_DIVESET*         diveset,            /**< general diving settings */
   SCIP_SET*             set,                /**< SCIP settings */
   SCIP_VAR*             divecand,           /**< the candidate for which the branching direction is requested */
   SCIP_Real             divecandsol,        /**< LP solution value of the candidate */
   SCIP_Real             divecandfrac        /**< fractionality of the candidate */
   );
/** copies the given primal heuristic to a new scip */
extern
SCIP_RETCODE SCIPheurCopyInclude(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< SCIP_SET of SCIP to copy to */
   );

/** creates a primal heuristic */
extern
SCIP_RETCODE SCIPheurCreate(
   SCIP_HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_MESSAGEHDLR*     messagehdlr,        /**< message handler */
   BMS_BLKMEM*           blkmem,             /**< block memory for parameter settings */
   const char*           name,               /**< name of primal heuristic */
   const char*           desc,               /**< description of primal heuristic */
   char                  dispchar,           /**< display character of primal heuristic */
   int                   priority,           /**< priority of the primal heuristic */
   int                   freq,               /**< frequency for calling primal heuristic */
   int                   freqofs,            /**< frequency offset for calling primal heuristic */
   int                   maxdepth,           /**< maximal depth level to call heuristic at (-1: no limit) */
   unsigned int          timingmask,         /**< positions in the node solving loop where heuristic should be executed */
   SCIP_Bool             usessubscip,        /**< does the separator use a secondary SCIP instance? */
   SCIP_DECL_HEURCOPY    ((*heurcopy)),      /**< copy method of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   SCIP_DECL_HEURFREE    ((*heurfree)),      /**< destructor of primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit)),      /**< initialize primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit)),      /**< deinitialize primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol)),   /**< solving process initialization method of primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol)),   /**< solving process deinitialization method of primal heuristic */
   SCIP_DECL_HEUREXEC    ((*heurexec)),      /**< execution method of primal heuristic */
   SCIP_HEURDATA*        heurdata            /**< primal heuristic data */
   );

/** calls destructor and frees memory of primal heuristic */
extern
SCIP_RETCODE SCIPheurFree(
   SCIP_HEUR**           heur,               /**< pointer to primal heuristic data structure */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** initializes primal heuristic */
extern
SCIP_RETCODE SCIPheurInit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** calls exit method of primal heuristic */
extern
SCIP_RETCODE SCIPheurExit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs primal heuristic that the branch and bound process is being started */
extern
SCIP_RETCODE SCIPheurInitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** informs primal heuristic that the branch and bound process data is being freed */
extern
SCIP_RETCODE SCIPheurExitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set                 /**< global SCIP settings */
   );

/** should the heuristic be executed at the given depth, frequency, timing, ... */
extern
SCIP_Bool SCIPheurShouldBeExecuted(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   int                   depth,              /**< depth of current node */
   int                   lpstateforkdepth,   /**< depth of the last node with solved LP */
   SCIP_HEURTIMING       heurtiming,         /**< current point in the node solving process */
   SCIP_Bool*            delayed             /**< pointer to store whether the heuristic should be delayed */
   );

/** calls execution method of primal heuristic */
extern
SCIP_RETCODE SCIPheurExec(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set,                /**< global SCIP settings */
   SCIP_PRIMAL*          primal,             /**< primal data */
   int                   depth,              /**< depth of current node */
   int                   lpstateforkdepth,   /**< depth of the last node with solved LP */
   SCIP_HEURTIMING       heurtiming,         /**< current point in the node solving process */
   SCIP_Bool             nodeinfeasible,     /**< was the current node already detected to be infeasible? */
   int*                  ndelayedheurs,      /**< pointer to count the number of delayed heuristics */
   SCIP_RESULT*          result              /**< pointer to store the result of the callback method */
   );

/** sets priority of primal heuristic */
extern
void SCIPheurSetPriority(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_SET*             set,                /**< global SCIP settings */
   int                   priority            /**< new priority of the primal heuristic */
   );

/** sets copy callback of primal heuristic */
extern
void SCIPheurSetCopy(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURCOPY    ((*heurcopy))       /**< copy callback of primal heuristic or NULL if you don't want to copy your plugin into sub-SCIPs */
   );

/** sets destructor callback of primal heuristic */
extern
void SCIPheurSetFree(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURFREE    ((*heurfree))       /**< destructor of primal heuristic */
   );

/** sets initialization callback of primal heuristic */
extern
void SCIPheurSetInit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINIT    ((*heurinit))       /**< initialize primal heuristic */
   );

/** sets deinitialization callback of primal heuristic */
extern
void SCIPheurSetExit(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXIT    ((*heurexit))       /**< deinitialize primal heuristic */
   );

/** sets solving process initialization callback of primal heuristic */
extern
void SCIPheurSetInitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEURINITSOL ((*heurinitsol))    /**< solving process initialization callback of primal heuristic */
   );

/** sets solving process deinitialization callback of primal heuristic */
extern
void SCIPheurSetExitsol(
   SCIP_HEUR*            heur,               /**< primal heuristic */
   SCIP_DECL_HEUREXITSOL ((*heurexitsol))    /**< solving process deinitialization callback of primal heuristic */
   );

#ifdef __cplusplus
}
#endif

#endif
