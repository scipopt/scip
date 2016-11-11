/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   struct_parallel.h
 * @ingroup PARAINTERFACE
 * @brief  the struct definitions for the SCIP parallel interface
 * @author Stephen J. Maher
 * @author Robert Lion Gottwald
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __STRUCT_SPI_H__
#define __STRUCT_SPI_H__

#include "scip/type_syncstore.h"
#include "tpi/type_tpi.h"
#include "scip/def.h"
#include "scip/type_scip.h"
#include "scip/type_stat.h"
#include "scip/type_lp.h"

#ifdef __cplusplus
extern "C" {
#endif

struct SCIP_SyncStore
{
   int                     references;         /**< reference count of spi */
   SCIP_PARALLELMODE       mode;               /**< the mode for the parallel solving */
   SCIP_Bool               initialized;
   SCIP_SYNCDATA*          syncdata;
   SCIP_SYNCDATA*          lastsync;

   SCIP*                   mainscip;
   SCIP_Bool               stopped;
   SCIP_LOCK               lock;

   /* SPI settings */
   int                     nsyncdata;
   int                     maxnsyncdelay;
   SCIP_Real               minsyncdelay;
   SCIP_Real               syncfreqinit;
   SCIP_Real               syncfreqmax;
   SCIP_Real               targetprogress;
   SCIP_Real               syncfreqfactor;
   int                     maxnsols;
   int                     nsolvers;
};


struct SCIP_SyncData
{
   SCIP_Real*              solobj;
   SCIP_Real**             sols;
   int*                    solsource;
   int                     nsols;
   SCIP_Real               bestlowerbound;
   SCIP_Real               bestupperbound;
   SCIP_Longint            syncnum;
   int                     winner;
   SCIP_STATUS             status;
   SCIP_LOCK               lock;
   int                     syncedcount;
   SCIP_CONDITION          allsynced;
   SCIP_BOUNDSTORE*        boundstore;
   SCIP_Real               syncfreq;
   SCIP_Longint            memtotal;
};


typedef struct {
   int pos[2];
} BoundPos;

typedef struct {
   int                 varidx;      /**< the variables position in the variable array of the main scip */
   SCIP_Real           newbound;    /**< the variables new bound */
   SCIP_BOUNDTYPE      boundtype;   /**< the type of the variables new bound */
} BoundChg;

struct SCIP_BoundStore {
   int                 nvars;
   BoundPos*           bndpos;
   BoundChg*           bndchg;
   int                 nbndchg;
   int                 bndchgsize;
};

#endif
