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

/**@file   type_heur.h
 * @ingroup TYPEDEFINITIONS
 * @brief  type definitions for primal heuristics
 * @author Tobias Achterberg
 * @author Timo Berthold
 *
 *  This file defines the interface for primal heuristics implemented in C.
 *
 *  - \ref HEUR "Instructions for implementing a primal heuristic"
 *  - \ref PRIMALHEURISTICS "List of available primal heuristics"
 *  - \ref scip::ObjHeur "C++ wrapper class"
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_TYPE_HEUR_H__
#define __SCIP_TYPE_HEUR_H__

#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_result.h"
#include "scip/type_timing.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct SCIP_Heur SCIP_HEUR;               /**< primal heuristic */
typedef struct SCIP_HeurData SCIP_HEURDATA;       /**< locally defined primal heuristic data */
typedef struct SCIP_Diveset SCIP_DIVESET;         /**< common parameters for all diving heuristics */

/** copy method for heuristic plugins (called when SCIP copies plugins)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEURCOPY(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** destructor of primal heuristic to free user data (called when SCIP is exiting)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEURFREE(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** initialization method of primal heuristic (called after problem was transformed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEURINIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** deinitialization method of primal heuristic (called before transformed problem is freed)
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEUREXIT(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** solving process initialization method of primal heuristic (called when branch and bound process is about to begin)
 *
 *  This method is called when the presolving was finished and the branch and bound process is about to begin.
 *  The primal heuristic may use this call to initialize its branch and bound specific data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEURINITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** solving process deinitialization method of primal heuristic (called before branch and bound process data is freed)
 *
 *  This method is called before the branch and bound process is freed.
 *  The primal heuristic should use this call to clean up its branch and bound data.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 */
#define SCIP_DECL_HEUREXITSOL(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur)

/** execution method of primal heuristic
 *
 *  Searches for feasible primal solutions. The method is called in the node processing loop.
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - heur            : the primal heuristic itself
 *  - heurtiming      : current point in the node solving loop
 *  - nodeinfeasible  : was the current node already detected to be infeasible?
 *  - result          : pointer to store the result of the heuristic call
 *
 *  possible return values for *result:
 *  - SCIP_FOUNDSOL   : at least one feasible primal solution was found
 *  - SCIP_DIDNOTFIND : the heuristic searched, but did not find a feasible solution
 *  - SCIP_DIDNOTRUN  : the heuristic was skipped
 *  - SCIP_DELAYED    : the heuristic was skipped, but should be called again as soon as possible, disregarding
 *                      its frequency
 */
#define SCIP_DECL_HEUREXEC(x) SCIP_RETCODE x (SCIP* scip, SCIP_HEUR* heur, SCIP_HEURTIMING heurtiming, \
      SCIP_Bool nodeinfeasible, SCIP_RESULT* result)


/* callbacks for diving heuristic specific settings */

/** returns a score for the given candidate -- the best candidate minimizes the diving score
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cand            : Branch candidate for which the score should be determined
 *  - candsol         : solution value of candidate in LP relaxation solution
 *  - candsfrac       : fractional part of solution value of candidate
 */
#define SCIP_DECL_DIVESETGETSCORE(x) SCIP_Real x (SCIP* scip, SCIP_VAR* cand, SCIP_Real candsol, SCIP_Real candsfrac)

/** returns the preferred branching direction of \p cand
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - cand            : Branch candidate for which the score should be determined
 *  - candsol         : solution value of candidate in LP relaxation solution
 *  - candsfrac       : fractional part of solution value of candidate
 */
#define SCIP_DECL_DIVESETCANDBRANCHDIR(x) SCIP_BRANCHDIR x (SCIP* scip, SCIP_VAR* cand, SCIP_Real candsol, SCIP_Real candsfrac)

/** return arrays of all diving candidates
 *
 * This callback is optional. If not declared, the diving will be performed
 *  on the set of LP branch candidates. This callback enables to consider other candidates than only LP branching candidates
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - diveset         : the dive set
 *  - divecands       : reference to array of possible diving candidates
 *  - divecandssol    : reference to array of diving candidate solution values
 *  - divecandsfrac   : reference to array of diving candidate fractionalities
 *  - ndivecands      : pointer to store the number of diving candidates
 */
#define SCIP_DECL_DIVESETGETCANDS(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIVESET* diveset, SCIP_VAR*** divecands, \
      SCIP_Real** divecandssol, SCIP_Real** divecandsfrac, int* ndivecands)

/** free candidates storage allocated before
 *
 *  this method should free all memory that was allocated during the divesetGetCandsXYZ callback. It is optional and only
 *  required if memory was allocated within the (optinal )divesetGetCandsXYZ callback
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - diveset         : the dive set
 *  - divecands       : reference to array of possible diving candidates
 *  - divecandssol    : reference to array of diving candidate solution values
 *  - divecandsfrac   : reference to array of diving candidate fractionalities
 *  - ndivecands      : the number of diving candidates that was previously allocated
 */
#define SCIP_DECL_DIVESETFREECANDS(x) SCIP_RETCODE x (SCIP* scip, SCIP_DIVESET* diveset, SCIP_VAR*** divecands, \
      SCIP_Real** divecandssol, SCIP_Real** divecandsfrac, int ndivecands)



#ifdef __cplusplus
}
#endif

#endif
