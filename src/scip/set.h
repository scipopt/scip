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
#include "nodesel.h"
#include "disp.h"
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
   NODESEL**        nodesels;           /**< node selectors */
   int              nnodesels;          /**< number of node selectors */
   int              nodeselssize;       /**< size of nodesels array */
   NODESEL*         nodesel;            /**< active node selector */
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
   int              nodelimit;          /**< maximal number of nodes to process */
   int              lpsolvefreq;        /**< frequency for solving LP at the nodes */
   unsigned int     usepricing:1;       /**< use pricing of variables */
};


extern
RETCODE SCIPsetCreate(                  /**< creates global SCIP settings */
   SET**            set,                /**< pointer to SCIP settings */
   SCIP*            scip                /**< SCIP data structure */   
   );

extern
RETCODE SCIPsetFree(                    /**< frees global SCIP settings */
   SET**            set                 /**< pointer to SCIP settings */
   );

extern
RETCODE SCIPsetIncludeReader(           /**< inserts file reader in file reader list */
   SET*             set,                /**< global SCIP settings */
   READER*          reader              /**< file reader */
   );

extern
RETCODE SCIPsetFindReader(              /**< finds the file reader of the given name */
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of file reader */
   READER**         reader              /**< pointer for storing the file reader (returns NULL, if not found) */
   );

extern
RETCODE SCIPsetIncludeConsHdlr(         /**< inserts constraint handler in constraint handler list */
   SET*             set,                /**< global SCIP settings */
   CONSHDLR*        conshdlr            /**< constraint handler */
   );

extern
RETCODE SCIPsetFindConsHdlr(            /**< finds the constraint handler of the given name */
   const SET*       set,                /**< global SCIP settings */
   const char*      name,               /**< name of constraint handler */
   CONSHDLR**       conshdlr            /**< pointer for storing the constraint handler (returns NULL, if not found) */
   );

extern
RETCODE SCIPsetIncludeNodesel(          /**< inserts node selector in node selector list */
   SET*             set,                /**< global SCIP settings */
   NODESEL*         nodesel             /**< node selector */
   );

extern
RETCODE SCIPsetIncludeDisp(             /**< inserts display column in display column list */
   SET*             set,                /**< global SCIP settings */
   DISP*            disp                /**< display column */
   );

extern
RETCODE SCIPsetInitCallbacks(           /**< initializes all user callback functions */
   const SET*       set                 /**< global SCIP settings */
   );

extern
RETCODE SCIPsetExitCallbacks(           /**< calls exit methods of all user callback functions */
   const SET*       set                 /**< global SCIP settings */
   );

extern
int SCIPsetCalcMemGrowSize(             /**< calculate memory size for dynamically allocated arrays */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

extern
int SCIPsetCalcTreeGrowSize(            /**< calculate memory size for tree array */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

extern
int SCIPsetCalcPathGrowSize(            /**< calculate memory size for path array */
   const SET*       set,                /**< global SCIP settings */
   int              num                 /**< minimum number of entries to store */
   );

extern
RETCODE SCIPsetSetVerbLevel(            /**< sets verbosity level for message output */
   SET*             set,                /**< global SCIP settings */
   VERBLEVEL        verblevel           /**< verbosity level for message output */
   );

extern
RETCODE SCIPsetSetFeastol(              /**< sets LP feasibility tolerance */
   SET*             set,                /**< global SCIP settings */
   LP*              lp,                 /**< actual LP data (or NULL) */
   Real             feastol             /**< new feasibility tolerance */
   );


#ifndef NDEBUG

/* In debug mode, the following methods are implemented as function calls to ensure
 * type validity.
 */

extern
Bool SCIPsetIsEQ(                       /**< checks, if values are in range of epsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsL(                        /**< checks, if val1 is (more than epsilon) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsLE(                       /**< checks, if val1 is not (more than epsilon) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsG(                        /**< checks, if val1 is (more than epsilon) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsGE(                       /**< checks, if val1 is not (more than epsilon) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsZero(                     /**< checks, if value is in range epsilon of 0.0 */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPsetIsPos(                      /**< checks, if value is greater than epsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPsetIsNeg(                      /**< checks, if value is lower than -epsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPsetIsSumEQ(                    /**< checks, if values are in range of sumepsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsSumL(                     /**< checks, if val1 is (more than sumepsilon) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsSumLE(                    /**< checks, if val1 is not (more than sumepsilon) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsSumG(                     /**< checks, if val1 is (more than sumepsilon) greater than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsSumGE(                    /**< checks, if val1 is not (more than sumepsilon) lower than val2 */
   const SET*       set,                /**< global SCIP settings */
   Real             val1,               /**< first value to be compared */
   Real             val2                /**< second value to be compared */
   );

extern
Bool SCIPsetIsSumZero(                  /**< checks, if value is in range sumepsilon of 0.0 */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPsetIsSumPos(                   /**< checks, if value is greater than sumepsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPsetIsSumNeg(                   /**< checks, if value is lower than -sumepsilon */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPsetIsInfinity(                 /**< checks, if value is (positive) infinite */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against infinity */
   );

extern
Bool SCIPsetIsFeasible(                 /**< checks, if value is non-negative within the LP feasibility bounds */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Real SCIPsetFloor(                      /**< rounds value down to the next integer */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Real SCIPsetCeil(                       /**< rounds value up to the next integer */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

extern
Bool SCIPsetIsIntegral(                 /**< checks, if value is integral within the LP feasibility bounds */
   const SET*       set,                /**< global SCIP settings */
   Real             val                 /**< value to be compared against zero */
   );

#else

/* In optimized mode, the methods are implemented as defines to reduce the number of function calls and
 * speed up the algorithms.
 */

#define SCIPsetIsEQ(set, val1, val2)    ( ABS((val1)-(val2)) <= (set)->epsilon )
#define SCIPsetIsL(set, val1, val2)     ( (val1) < (val2) - (set)->epsilon )
#define SCIPsetIsLE(set, val1, val2)    ( (val1) <= (val2) + (set)->epsilon )
#define SCIPsetIsG(set, val1, val2)     ( (val1) > (val2) + (set)->epsilon )
#define SCIPsetIsGE(set, val1, val2)    ( (val1) >= (val2) - (set)->epsilon )
#define SCIPsetIsZero(set, val)         ( ABS(val) <= (set)->epsilon )
#define SCIPsetIsPos(set, val)          ( (val) > (set)->epsilon )
#define SCIPsetIsNeg(set, val)          ( (val) < -(set)->epsilon )

#define SCIPsetIsSumEQ(set, val1, val2) ( ABS((val1)-(val2)) <= (set)->sumepsilon )
#define SCIPsetIsSumL(set, val1, val2)  ( (val1) < (val2) - (set)->sumepsilon )
#define SCIPsetIsSumLE(set, val1, val2) ( (val1) <= (val2) + (set)->sumepsilon )
#define SCIPsetIsSumG(set, val1, val2)  ( (val1) > (val2) + (set)->sumepsilon )
#define SCIPsetIsSumGE(set, val1, val2) ( (val1) >= (val2) - (set)->sumepsilon )
#define SCIPsetIsSumZero(set, val)      ( ABS(val) <= (set)->sumepsilon )
#define SCIPsetIsSumPos(set, val)       ( (val) > (set)->sumepsilon )
#define SCIPsetIsSumNeg(set, val)       ( (val) < -(set)->sumepsilon )

#define SCIPsetIsInfinity(set, val)     ( (val) >= (set)->infinity )
#define SCIPsetIsFeasible(set, val)     ( (val) >= -(set)->feastol )
#define SCIPsetFloor(set, val)          ( floor((val) + (set)->feastol) )
#define SCIPsetCeil(set, val)           ( ceil((val) - (set)->feastol) )
#define SCIPsetIsIntegral(set, val)     ( SCIPsetCeil(set, val) - (val) <= (set)->feastol )

#endif


#define SCIPsetCaptureBufferArray(set,ptr,num)   ( SCIPbufferCapture((set)->buffer, set, (void**)(&(ptr)), \
                                                   (num)*sizeof(*(ptr))) )
#define SCIPsetReleaseBufferArray(set,ptr)       ( SCIPbufferRelease((set)->buffer, (void**)(&ptr)) )
#define SCIPsetCaptureBufferSize(set,ptr,size)   ( SCIPbufferCapture((set)->buffer, set, (void**)(&(ptr)), size) )
#define SCIPsetReleaseBufferSize(set,ptr)        ( SCIPbufferRelease((set)->buffer, (void**)(&ptr)) )


#endif
