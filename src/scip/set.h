/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2003 Tobias Achterberg                              */
/*                            Thorsten Koch                                  */
/*                  2002-2003 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the SCIP Academic Licence.        */
/*                                                                           */
/*  You should have received a copy of the SCIP Academic License             */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   set.h
 * @brief  global SCIP settings
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SET_H__
#define __SET_H__


/** possible settings for enabling/disabling algorithms and other features */
enum Setting
{
   SCIP_UNDEFINED = 0,                  /**< undefined setting */
   SCIP_DISABLED  = 1,                  /**< feature is disabled */
   SCIP_AUTO      = 2,                  /**< feature is set to automatic mode */
   SCIP_ENABLED   = 3                   /**< feature is enabled */
};
typedef enum Setting SETTING;

typedef struct Set SET;                 /**< global SCIP settings */

#include <math.h>

#include "def.h"
#include "sort.h"
#include "scip.h"
#include "reader.h"
#include "cons.h"
#include "event.h"
#include "nodesel.h"
#include "disp.h"
#include "branch.h"
#include "lp.h"
#include "message.h"
#include "buffer.h"


/** global SCIP settings */
struct Set
{
   SCIP*            scip;               /**< very ugly: pointer to scip main data structure for callback methods */
   VERBLEVEL        verblevel;          /**< verbosity level of output */
   Real             epsilon;            /**< absolute values smaller than this are considered zero */
   Real             sumepsilon;         /**< absolute values of sums smaller than this are considered zero */
   Real             infinity;           /**< values larger than this are considered infinity */
   Real             feastol;            /**< LP feasibility tolerance */
   Real             memGrowFac;         /**< memory growing factor for dynamically allocated arrays */
   int              memGrowInit;        /**< initial size of dynamically allocated arrays */
   Real             treeGrowFac;        /**< memory growing factor for tree array */
   int              treeGrowInit;       /**< initial size of tree array */
   Real             pathGrowFac;        /**< memory growing factor for path array */
   int              pathGrowInit;       /**< initial size of path array */
   BUFFER*          buffer;             /**< memory buffers for short living temporary objects */
   READER**         readers;            /**< file readers */
   int              nreaders;           /**< number of file readers */
   int              readerssize;        /**< size of readers array */
   CONSHDLR**       conshdlrs;          /**< constraint handlers */
   int              nconshdlrs;         /**< number of constraint handlers */
   int              conshdlrssize;      /**< size of conshdlrs array */
   HEUR**           heurs;              /**< primal heuristics */
   int              nheurs;             /**< number of primal heuristics */
   int              heurssize;          /**< size of heurs array */
   EVENTHDLR**      eventhdlrs;         /**< event handlers */
   int              neventhdlrs;        /**< number of event handlers */
   int              eventhdlrssize;     /**< size of eventhdlrs array */
   NODESEL**        nodesels;           /**< node selectors */
   int              nnodesels;          /**< number of node selectors */
   int              nodeselssize;       /**< size of nodesels array */
   NODESEL*         nodesel;            /**< active node selector */
   BRANCHRULE**     branchrules;        /**< branching rules */
   int              nbranchrules;       /**< number of branching rules */
   int              branchrulessize;    /**< size of branchrules array */
   DISP**           disps;              /**< display columns */
   int              ndisps;             /**< number of display columns */
   int              dispssize;          /**< size of disps array */
   int              dispwidth;          /**< maximal number of characters in a node information line */
   int              dispfreq;           /**< frequency for displaying node information lines */
   int              dispheaderfreq;     /**< frequency for displaying header lines (every n'th node information line) */
   int              maxpricevars;       /**< maximal number of variables priced in per pricing round */
   int              maxpricevarsroot;   /**< maximal number of priced variables at the root node */
   Real             abortpricevarsfac;  /**< pricing is aborted, if fac * maxpricevars pricing candidates were found */
   int              maxsepacuts;        /**< maximal number of cuts separated per separation round */
   int              maxsepacutsroot;    /**< maximal number of separated cuts at the root node */
   int              agelimit;           /**< maximum age a cut can reach before it is deleted from the global cut pool */
   int              maxsol;             /**< maximal number of solutions to store in the solution storage */
   Longint          nodelimit;          /**< maximal number of nodes to process */
   int              lpsolvefreq;        /**< frequency for solving LP at the nodes */
   int              lpsolvedepth;       /**< maximal depth for solving LP at the nodes */
   unsigned int     usepricing:1;       /**< use pricing of variables */
};


/** creates global SCIP settings */
extern
RETCODE SCIPsetCreate(
   SET**            set,                /**< pointer to SCIP settings */
   SCIP*            scip                /**< SCIP data structure */   
   );

/** frees global SCIP settings */
extern
RETCODE SCIPsetFree(
   SET**            set                 /**< pointer to SCIP settings */
   );

/** inserts file reader in file reader list */
extern
RETCODE SCIPsetIncludeReader(
   SET*             set,                /**< global SCIP settings */
   READER*          reader              /**< file reader */
   );

/** finds the file reader of the given name */
extern
RETCODE SCIPsetFindReader(
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of file reader */
   READER**         reader              /**< pointer for storing the file reader (returns NULL, if not found) */
   );

/** inserts constraint handler in constraint handler list */
extern
RETCODE SCIPsetIncludeConsHdlr(
   SET*             set,                /**< global SCIP settings */
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

/** finds the constraint handler of the given name */
extern
RETCODE SCIPsetFindConsHdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer for storing the constraint handler (returns NULL, if not found) */
   );

/** inserts primal heuristic in primal heuristic list */
extern
RETCODE SCIPsetIncludeHeur(
   SET*             set,                /**< global SCIP settings */
   HEUR*            heur                /**< primal heuristic */
   );

/** finds the primal heuristic of the given name */
extern
RETCODE SCIPsetFindHeur(
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of primal heuristic */
   HEUR**           heur                /**< pointer for storing the primal heuristic (returns NULL, if not found) */
   );

/** inserts event handler in event handler list */
extern
RETCODE SCIPsetIncludeEventHdlr(
   SET*             set,                /**< global SCIP settings */
   EVENTHDLR*       eventhdlr           /**< event handler */
   );

/** finds the event handler of the given name */
extern
RETCODE SCIPsetFindEventHdlr(
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of event handler */
   EVENTHDLR**      eventhdlr           /**< pointer for storing the event handler (returns NULL, if not found) */
   );

/** inserts node selector in node selector list */
extern
RETCODE SCIPsetIncludeNodesel(
   SET*             set,                /**< global SCIP settings */
   NODESEL*         nodesel             /**< node selector */
   );

/** inserts branching rule in branching rule list */
extern
RETCODE SCIPsetIncludeBranchrule(
   SET*             set,                /**< global SCIP settings */
   BRANCHRULE*      branchrule          /**< branching rule */
   );

/** inserts display column in display column list */
extern
RETCODE SCIPsetIncludeDisp(
   SET*             set,                /**< global SCIP settings */
   DISP*            disp                /**< display column */
   );

/** initializes all user callback functions */
extern
RETCODE SCIPsetInitCallbacks(
   const SET*       set                 /**< global SCIP settings */
   );

/** calls exit methods of all user callback functions */
extern
RETCODE SCIPsetExitCallbacks(
   const SET*       set                 /**< global SCIP settings */
   );

/** calculate memory size for dynamically allocated arrays */
extern
int SCIPsetCalcMemGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** calculate memory size for tree array */
extern
int SCIPsetCalcTreeGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** calculate memory size for path array */
extern
int SCIPsetCalcPathGrowSize(
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

/** sets verbosity level for message output */
extern
RETCODE SCIPsetSetVerbLevel(
   SET*             set,                /**< global SCIP settings */
   VERBLEVEL        verblevel           /**< verbosity level for message output */
   );

/** sets LP feasibility tolerance */
extern
RETCODE SCIPsetSetFeastol(
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data (or NULL) */
   Real             feastol             /**< new feasibility tolerance */
   );


/** returns the relative difference: (val1-val2)/max(|val1|,|val2|,1.0) */
extern
Real SCIPsetRelDiff(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

/** checks, if values are in range of epsilon */
extern
Bool SCIPsetIsEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) lower than val2 */
extern
Bool SCIPsetIsLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) greater than val2 */
extern
Bool SCIPsetIsLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than epsilon) greater than val2 */
extern
Bool SCIPsetIsGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than epsilon) lower than val2 */
extern
Bool SCIPsetIsGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range epsilon of 0.0 */
extern
Bool SCIPsetIsZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than epsilon */
extern
Bool SCIPsetIsPos(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -epsilon */
extern
Bool SCIPsetIsNeg(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if values are in range of sumepsilon */
extern
Bool SCIPsetIsSumEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) lower than val2 */
extern
Bool SCIPsetIsSumLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) greater than val2 */
extern
Bool SCIPsetIsSumLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than sumepsilon) greater than val2 */
extern
Bool SCIPsetIsSumGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than sumepsilon) lower than val2 */
extern
Bool SCIPsetIsSumGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range sumepsilon of 0.0 */
extern
Bool SCIPsetIsSumZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than sumepsilon */
extern
Bool SCIPsetIsSumPos(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -sumepsilon */
extern
Bool SCIPsetIsSumNeg(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if values are in range of feasibility tolerance */
extern
Bool SCIPsetIsFeasEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than feasibility tolerance) lower than val2 */
extern
Bool SCIPsetIsFeasLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than feasibility tolerance) greater than val2 */
extern
Bool SCIPsetIsFeasLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is (more than feasibility tolerance) greater than val2 */
extern
Bool SCIPsetIsFeasGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if val1 is not (more than feasibility tolerance) lower than val2 */
extern
Bool SCIPsetIsFeasGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is in range feasibility tolerance of 0.0 */
extern
Bool SCIPsetIsFeasZero(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is greater than feasibility tolerance */
extern
Bool SCIPsetIsFeasPos(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if value is lower than -feasibility tolerance */
extern
Bool SCIPsetIsFeasNeg(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if relative difference of values is in range of epsilon */
extern
Bool SCIPsetIsRelEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than epsilon */
extern
Bool SCIPsetIsRelLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than epsilon */
extern
Bool SCIPsetIsRelLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than epsilon */
extern
Bool SCIPsetIsRelGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -epsilon */
extern
Bool SCIPsetIsRelGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of values is in range of sumepsilon */
extern
Bool SCIPsetIsSumRelEQ(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is lower than sumepsilon */
extern
Bool SCIPsetIsSumRelLT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not greater than sumepsilon */
extern
Bool SCIPsetIsSumRelLE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is greater than sumepsilon */
extern
Bool SCIPsetIsSumRelGT(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if relative difference of val1 and val2 is not lower than -sumepsilon */
extern
Bool SCIPsetIsSumRelGE(
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

/** checks, if value is (positive) infinite */
extern
Bool SCIPsetIsInfinity(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against infinity */
   );

/** checks, if value is non-negative within the LP feasibility bounds */
extern
Bool SCIPsetIsFeasible(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value down to the next integer */
extern
Real SCIPsetFloor(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** rounds value up to the next integer */
extern
Real SCIPsetCeil(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** returns fractional part of value, i.e. ceil(x) - x */
extern
Real SCIPsetFrac(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to return fractional part for */
   );

/** checks, if value is integral within the LP feasibility bounds */
extern
Bool SCIPsetIsIntegral(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if given fractional part is smaller than feastol */
extern
Bool SCIPsetIsFracIntegral(
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

/** checks, if the given integer bounds correspond to a fixed interval */
extern
Bool SCIPsetIsFixed(
   const SET*       set,                /**< global SCIP settings */
   Real             lb,                 /**< lower integer bound */
   Real             ub                  /**< upper integer bound */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsetIsEQ(set, val1, val2)       ( EPSEQ(val1, val2, (set)->epsilon) )
#define SCIPsetIsLT(set, val1, val2)       ( EPSLT(val1, val2, (set)->epsilon) )
#define SCIPsetIsLE(set, val1, val2)       ( EPSLE(val1, val2, (set)->epsilon) )
#define SCIPsetIsGT(set, val1, val2)       ( EPSGT(val1, val2, (set)->epsilon) )
#define SCIPsetIsGE(set, val1, val2)       ( EPSGE(val1, val2, (set)->epsilon) )
#define SCIPsetIsZero(set, val)            ( EPSZ(val, (set)->epsilon) )
#define SCIPsetIsPos(set, val)             ( EPSP(val, (set)->epsilon) )
#define SCIPsetIsNeg(set, val)             ( EPSN(val, (set)->epsilon) )

#define SCIPsetIsSumEQ(set, val1, val2)    ( EPSEQ(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumLT(set, val1, val2)    ( EPSLT(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumLE(set, val1, val2)    ( EPSLE(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumGT(set, val1, val2)    ( EPSGT(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumGE(set, val1, val2)    ( EPSGE(val1, val2, (set)->sumepsilon) )
#define SCIPsetIsSumZero(set, val)         ( EPSZ(val, (set)->sumepsilon) )
#define SCIPsetIsSumPos(set, val)          ( EPSP(val, (set)->sumepsilon) )
#define SCIPsetIsSumNeg(set, val)          ( EPSN(val, (set)->sumepsilon) )

#define SCIPsetIsFeasEQ(set, val1, val2)   ( EPSEQ(val1, val2, (set)->feastol) )
#define SCIPsetIsFeasLT(set, val1, val2)   ( EPSLT(val1, val2, (set)->feastol) )
#define SCIPsetIsFeasLE(set, val1, val2)   ( EPSLE(val1, val2, (set)->feastol) )
#define SCIPsetIsFeasGT(set, val1, val2)   ( EPSGT(val1, val2, (set)->feastol) )
#define SCIPsetIsFeasGE(set, val1, val2)   ( EPSGE(val1, val2, (set)->feastol) )
#define SCIPsetIsFeasZero(set, val)        ( EPSZ(val, (set)->feastol) )
#define SCIPsetIsFeasPos(set, val)         ( EPSP(val, (set)->feastol) )
#define SCIPsetIsFeasNeg(set, val)         ( EPSN(val, (set)->feastol) )

#define SCIPsetIsRelEQ(set, val1, val2)    ( EPSZ(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )
#define SCIPsetIsRelLT(set, val1, val2)    ( EPSN(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )
#define SCIPsetIsRelLE(set, val1, val2)    ( !EPSP(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )
#define SCIPsetIsRelGT(set, val1, val2)    ( EPSP(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )
#define SCIPsetIsRelGE(set, val1, val2)    ( !EPSN(SCIPsetRelDiff(set, val1, val2), (set)->epsilon) )

#define SCIPsetIsSumRelEQ(set, val1, val2) ( EPSZ(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )
#define SCIPsetIsSumRelLT(set, val1, val2) ( EPSN(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )
#define SCIPsetIsSumRelLE(set, val1, val2) ( !EPSP(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )
#define SCIPsetIsSumRelGT(set, val1, val2) ( EPSP(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )
#define SCIPsetIsSumRelGE(set, val1, val2) ( !EPSN(SCIPsetRelDiff(set, val1, val2), (set)->sumepsilon) )

#define SCIPsetIsInfinity(set, val)        ( (val) >= (set)->infinity )
#define SCIPsetIsFeasible(set, val)        ( (val) >= -(set)->feastol )
#define SCIPsetFloor(set, val)             ( floor((val) + (set)->feastol) )
#define SCIPsetCeil(set, val)              ( ceil((val) - (set)->feastol) )
#define SCIPsetFrac(set, val)              ( (val) - SCIPsetFloor(set, val) )
#define SCIPsetIsIntegral(set, val)        ( EPSLE(SCIPsetCeil(set, val), val, (set)->feastol) )
#define SCIPsetIsFracIntegral(set, val)    ( !EPSP(val, (set)->feastol) )
#define SCIPsetIsFixed(set, lb, ub)        ( SCIPsetIsEQ(set, lb, ub) )

#endif


#define SCIPsetCaptureBufferArray(set,ptr,num)   ( SCIPbufferCapture((set)->buffer, set, (void**)(ptr), \
                                                   (num)*sizeof(**(ptr))) )
#define SCIPsetReleaseBufferArray(set,ptr)       ( SCIPbufferRelease((set)->buffer, (void**)(ptr), 0*sizeof(**(ptr))) )
#define SCIPsetCaptureBufferSize(set,ptr,size)   ( SCIPbufferCapture((set)->buffer, set, (void**)(ptr), size) )
#define SCIPsetReleaseBufferSize(set,ptr)        ( SCIPbufferRelease((set)->buffer, (void**)(ptr), 0) )


#endif
