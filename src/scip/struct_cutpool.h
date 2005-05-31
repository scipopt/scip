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
#pragma ident "@(#) $Id: struct_cutpool.h,v 1.7 2005/05/31 17:20:22 bzfpfend Exp $"

/**@file   struct_cutpool.h
 * @brief  datastructures for storing cuts in a cut pool
 * @author Tobias Achterberg
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_CUTPOOL_H__
#define __STRUCT_CUTPOOL_H__


#include "scip/def.h"
#include "scip/type_clock.h"
#include "scip/type_misc.h"
#include "scip/type_lp.h"
#include "scip/type_cutpool.h"



/** datastructure for cuts in a cut pool */
struct Cut
{
   ROW*             row;                /**< LP row of this cut */
   int              age;                /**< age of the cut: number of successive times, the cut was not violated */
   int              processedlp;        /**< last LP, where this cut was processed */
   int              pos;                /**< position of cut in the cuts array of the cut pool */
};

/** storage for pooled cuts */
struct Cutpool
{
   Longint          ncalls;             /**< number of times, the cutpool was separated */
   Longint          ncutsfound;         /**< total number of cuts that were separated from the pool */
   CLOCK*           clock;              /**< separation time */
   HASHTABLE*       hashtable;          /**< hash table to identify already stored cuts */
   CUT**            cuts;               /**< stored cuts of the pool */
   int              cutssize;           /**< size of cuts array */
   int              ncuts;              /**< number of cuts stored in the pool */
   int              agelimit;           /**< maximum age a cut can reach before it is deleted from the pool */
   int              processedlp;        /**< last LP that has been processed */
   int              firstunprocessed;   /**< first cut that has not been processed in the last LP */
   int              maxncuts;           /**< maximal number of cuts stored in the pool at the same time */
};


#endif
