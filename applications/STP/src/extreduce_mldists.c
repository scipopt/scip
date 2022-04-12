/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2022 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   extreduce_mldists.c
 * @brief  multi-level distance storage methods for Steiner tree extended reductions
 * @author Daniel Rehfeldt
 *
 * This file implements utility methods for Steiner tree problem extended reduction techniques
 * that allow use to store and retrieve (special) distances between vertices of the extension tree.
 *
 * A list of all interface methods can be found in extreduce.h.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

// #define SCIP_DEBUG
#include "extreduce.h"
#include "misc_stp.h"
#include "portab.h"



#define MLDISTS_EMPTYSLOT_NONE -1


/** Structure for storing distances in the extension tree.
 *  Organized in slots that can be filled by the user.
 *  On each level there are a number of slots available (specified by the user).
 *  Each slots consists of a base (id) and a number of targets. Each target has a distance and an ID.
 *  Each slot on a level has the same number of targets, namely level_ntargets[level].  */
struct multi_level_distances_storage
{
   int*                  target_ids;        /**< target ids only in DEBUG mode! */
   SCIP_Real*            target_dists;      /**< target ids */
   int*                  base_ids;          /**< bases ids */
   int*                  level_basestart;   /**< start of bases for given level */
   int*                  level_targetstart; /**< start of targets for given level */
   int*                  level_ntargets;    /**< number of targets per base on given level */
   int                   level_maxntargets; /**< maximum number of targets per level */
   int                   level_maxnslots;   /**< maximum number of bases per level */
   int                   nlevels;           /**< number of levels */
   int                   maxnlevels;        /**< maximum number of levels */
   int                   maxntargets;       /**< total maximum number of targets */
   int                   maxnslots;         /**< total maximum number of bases */
   int                   emptyslot_number;  /**< number (0,...) of current empty slot, or EMPTYSLOT_NONE if none exists */
   SCIP_Bool             target_withids;    /**< use ids? */
};



/**@name Local methods
 *
 * @{
 */


/** gets current level */
static inline
int mldistsGetTopLevel(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(mldists->nlevels > 0 && mldists->nlevels <= mldists->maxnlevels);

   return (mldists->nlevels - 1);
}


/** gets start of bases */
static inline
int mldistsGetPosBasesStart(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level               /**< the level */
)
{
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(mldists->level_basestart[level] >= 0 && mldists->level_basestart[level] < mldists->maxnslots);

   return mldists->level_basestart[level];
}


/** gets start of bases */
static inline
int mldistsGetPosBasesEnd(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level               /**< the level */
)
{
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(mldists->level_basestart[level + 1] >= 0 && mldists->level_basestart[level] <= mldists->maxnslots);

   return mldists->level_basestart[level + 1];
}


/** Gets (internal) position of given base id,
 *  or -1 if it could not be found. */
static
int mldistsGetPosBase(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base id */
)
{
   const int start = mldistsGetPosBasesStart(mldists, level);
   const int end = mldistsGetPosBasesEnd(mldists, level);
   const int* const base_ids = mldists->base_ids;

   assert(baseid >= 0);

   /* top level and with empty slots? */
   if( level == mldistsGetTopLevel(mldists) && extreduce_mldistsEmptySlotExists(mldists) )
   {
      const int emptyslotnumber = mldists->emptyslot_number;

      for( int i = start, j = 0; i < end && j < emptyslotnumber; ++i, ++j )
      {
         const int id = base_ids[i];

         assert(id >= 0);

         if( id == baseid )
            return i;
      }
   }
   else
   {
      for( int i = start; i < end; ++i )
      {
         const int id = base_ids[i];

         assert(id >= 0 || id == STP_MLDISTS_ID_EMPTY);

         if( id == baseid )
            return i;
      }
   }

   return -1;
}


/** gets (internal) position of targets for given base id */
static inline
int mldistsGetPosTargetsStart(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base id */
)
{
   int targetpos;
   int offset = 0;
   const int ntargets = mldists->level_ntargets[level];

   assert(extreduce_mldistsLevelContainsBase(mldists, level, baseid));
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(baseid >= 0);

   offset = mldistsGetPosBase(mldists, level, baseid) - mldists->level_basestart[level];

   assert(offset >= 0);
   assert(ntargets >= 0);

   targetpos = mldists->level_targetstart[level] + offset * ntargets;

   assert(targetpos >= 0 && targetpos < mldists->maxntargets);

   return targetpos;
}


/** gets (internal) position of targets for given base id and target id */
static inline
int mldistsGetPosTargets(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid,             /**< the base id */
   int                   targetid            /**< the target id */
)
{
   const int start = mldistsGetPosTargetsStart(mldists, level, baseid);
   const int end = start + mldists->level_ntargets[level];
   const int* const target_ids = mldists->target_ids;
   int i;

   assert(mldists->target_withids);

   for( i = start; i != end; i++ )
   {
      if( target_ids[i] == targetid )
      {
         break;
      }
   }

   assert(i != end);

   return i;
}


/** gets targets start of current empty slot */
static inline
int mldistsGetPosEmptyTargetsStart(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);
   const int ntargets = mldists->level_ntargets[level];
   const int start = mldists->level_targetstart[level] + mldists->emptyslot_number * ntargets;

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(start >= 0 && start < mldists->maxntargets);
   assert(mldists->emptyslot_number >= 0);
   assert(ntargets >= 0 && ntargets <= mldists->level_maxntargets);

   return start;
}


/** only for debugging */
static inline
void mldistsTopLevelUnset(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
#ifndef NDEBUG
   const int nlevels = mldists->nlevels;
   const int start_base = mldists->level_basestart[nlevels - 1];
   const int end_base = mldists->level_basestart[nlevels];
   const int start_target = mldists->level_targetstart[nlevels - 1];
   const int end_target = mldists->level_targetstart[nlevels];

   assert(start_base >= 0);
   assert(start_base <= end_base);
   assert(end_base <= mldists->maxnslots);

   for( int i = start_base; i < end_base; ++i )
   {
      mldists->base_ids[i] = STP_MLDISTS_ID_UNSET;
   }

   assert(start_target >= 0);
   assert(start_target <= end_target);
   assert(end_target <= mldists->maxntargets);

   for( int i = start_target; i < end_target; ++i )
   {
      mldists->target_dists[i] = STP_MLDISTS_DIST_UNSET;
   }

   if( mldists->target_withids )
   {
      assert(mldists->target_ids);

      for( int i = start_target; i < end_target; ++i )
      {
         mldists->target_ids[i] = STP_MLDISTS_ID_UNSET;
      }
   }
#endif
}



/**@} */

/**@name Interface methods
 *
 * @{
 */



/** initializes multi-level distances structure */
SCIP_RETCODE extreduce_mldistsInit(
   SCIP*                 scip,               /**< SCIP */
   int                   maxnlevels,         /**< maximum number of levels that can be handled */
   int                   maxnslots,          /**< maximum number of of slots (per level) that can be handled */
   int                   maxntargets,        /**< maximum number of of targets (per slot) that can be handled */
   int                   emptyslot_nbuffers, /**< buffer entries for empty slot (dists and IDs array is that much longer) */
   SCIP_Bool             use_targetids,      /**< use target IDs? */
   MLDISTS**             mldistances         /**< to be initialized */
)
{
   MLDISTS* mldists;
   /* NOTE: was also need to consider dummy/empty entry! */
   const int maxnslots_internal = maxnslots + 1;

   assert(scip && mldistances);
   assert(maxnlevels >= 1 && maxnslots >= 1 && maxntargets >= 1);
   assert(emptyslot_nbuffers >= 0);

   SCIP_CALL( SCIPallocMemory(scip, mldistances) );

   mldists = *mldistances;
   mldists->nlevels = 0;
   mldists->maxnlevels = maxnlevels;
   mldists->level_maxntargets = maxntargets;
   mldists->level_maxnslots = maxnslots_internal;
   mldists->target_withids = use_targetids;
   mldists->maxnslots = maxnlevels * maxnslots_internal;
   mldists->maxntargets = maxnlevels * maxnslots_internal * maxntargets + emptyslot_nbuffers;
   mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;

   if( use_targetids )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->target_ids), mldists->maxntargets) );
   }
   else
   {
      mldists->target_ids = NULL;
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->target_dists), mldists->maxntargets) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->base_ids), mldists->maxnslots) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->level_basestart), maxnlevels + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->level_targetstart), maxnlevels + 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &(mldists->level_ntargets), maxnlevels) );

   mldists->level_basestart[0] = 0;
   mldists->level_targetstart[0] = 0;

   return SCIP_OKAY;
}


/** frees  multi-level distances structure */
void extreduce_mldistsFree(
   SCIP*                 scip,               /**< SCIP */
   MLDISTS**             mldistances         /**< to be freed */
)
{
   MLDISTS* mldists;

   assert(scip && mldistances);

   mldists = *mldistances;

   SCIPfreeMemoryArray(scip, &(mldists->level_ntargets));
   SCIPfreeMemoryArray(scip, &(mldists->level_targetstart));
   SCIPfreeMemoryArray(scip, &(mldists->level_basestart));
   SCIPfreeMemoryArray(scip, &(mldists->base_ids));
   SCIPfreeMemoryArray(scip, &(mldists->target_dists));

   if( mldists->target_withids )
   {
      assert(mldists->target_ids);

      SCIPfreeMemoryArray(scip, &(mldists->target_ids));
   }

   SCIPfreeMemory(scip, mldistances);
}

/** empty? */
SCIP_Bool extreduce_mldistsIsEmpty(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(mldists);

   return (mldists->nlevels == 0);
}


/** does an empty slot exits? (on current level) */
SCIP_Bool extreduce_mldistsEmptySlotExists(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(MLDISTS_EMPTYSLOT_NONE == mldists->emptyslot_number || mldists->emptyslot_number >= 0);

   return (mldists->emptyslot_number != MLDISTS_EMPTYSLOT_NONE);
}


/** gets targets IDs memory from clean slot (to be filled in) */
int* extreduce_mldistsEmptySlotTargetIds(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int start = mldistsGetPosEmptyTargetsStart(mldists);

#ifndef NDEBUG
   const int level = mldistsGetTopLevel(mldists);
   const int end = start + mldists->level_ntargets[level];

   assert(mldists->target_withids);
   assert(mldists->target_ids);

   for( int i = start; i < end; ++i )
   {
      assert(mldists->target_ids[i] == STP_MLDISTS_ID_UNSET);
   }
#endif

   return &(mldists->target_ids[start]);
}

/** gets targets IDs memory from clean slot (to be filled in) */
int* extreduce_mldistsEmptySlotTargetIdsDirty(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int start = mldistsGetPosEmptyTargetsStart(mldists);

   assert(mldists->target_ids);

   return &(mldists->target_ids[start]);
}


/** Gets targets distances memory from clean slot (to be filled in).
 *  NOTE: Can only be called as long as this empty slots' distances have not not modified! */
SCIP_Real* extreduce_mldistsEmptySlotTargetDists(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int start = mldistsGetPosEmptyTargetsStart(mldists);

#ifndef NDEBUG
   const int level = mldistsGetTopLevel(mldists);
   const int end = start + mldists->level_ntargets[level];

   for( int i = start; i < end; ++i )
   {
      assert(EQ(mldists->target_dists[i], STP_MLDISTS_DIST_UNSET));
   }
#endif

   return &(mldists->target_dists[start]);
}


/** Gets targets distances memory from empty slot.
 *  NOTE: This method does not make sure that the distances are clean! (i.e. not already set) */
SCIP_Real* extreduce_mldistsEmptySlotTargetDistsDirty(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   const int start = mldistsGetPosEmptyTargetsStart(mldists);

   return &(mldists->target_dists[start]);
}



/** gets level of current empty slot */
int extreduce_mldistsEmptySlotLevel(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(extreduce_mldistsEmptySlotExists(mldists));

   return mldistsGetTopLevel(mldists);
}


/** sets base of empty slot */
void extreduce_mldistsEmptySlotSetBase(
   int                   baseid,             /**< base */
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);
   const int position = mldists->level_basestart[level] + mldists->emptyslot_number;

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(position >= 0 && position < mldists->maxnslots);
   assert(baseid >= 0 || baseid == STP_MLDISTS_ID_EMPTY);
   assert(mldists->base_ids[position] == STP_MLDISTS_ID_UNSET);

   mldists->base_ids[position] = baseid;
}


/** Resets all changes (especially bases) of empty slot.
 *  NOTE: Assumes that at least the basis of the slot has been set */
void extreduce_mldistsEmptySlotReset(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);
   const int position = mldists->level_basestart[level] + mldists->emptyslot_number;

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(position >= 0 && position < mldists->maxnslots);
   assert(mldists->base_ids[position] != STP_MLDISTS_ID_UNSET);

   mldists->base_ids[position] = STP_MLDISTS_ID_UNSET;

#ifndef NDEBUG
   {
      const int target_start = mldistsGetPosEmptyTargetsStart(mldists);
      const int target_end = target_start + mldists->level_ntargets[level];

      for( int i = target_start; i != target_end; i++ )
      {
         mldists->target_dists[i] = STP_MLDISTS_DIST_UNSET;
      }

      if( mldists->target_withids )
      {
         assert(mldists->target_ids);

         for( int i = target_start; i != target_end; i++ )
         {
            mldists->target_ids[i] = STP_MLDISTS_ID_UNSET;
         }
      }
   }
#endif
}


/** marks current empty slot as filled */
void extreduce_mldistsEmptySlotSetFilled(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int level = mldistsGetTopLevel(mldists);

   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(mldists->emptyslot_number < extreduce_mldistsLevelNSlots(mldists, level));

#ifndef NDEBUG
   {
      const int position = mldists->level_basestart[level] + mldists->emptyslot_number;

      assert(mldists->base_ids[position] != STP_MLDISTS_ID_UNSET);
   }
#endif

   mldists->emptyslot_number++;

   /* all slots of current level used? */
   if( mldists->emptyslot_number >= extreduce_mldistsLevelNSlots(mldists, level) )
   {
      mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;
   }
}


/** adds another level of slots at top */
void extreduce_mldistsLevelAddTop(
   int                   nslots,             /**< number of slots per this level */
   int                   nslottargets,       /**< number of targets per slot */
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   int nlevels;

   assert(!extreduce_mldistsEmptySlotExists(mldists));
   assert(nslots > 0 && nslots <= mldists->level_maxnslots);
   assert(nslottargets >= 0 && nslottargets <= mldists->level_maxntargets);
   assert(mldists->nlevels < mldists->maxnlevels);

   mldists->emptyslot_number = 0;
   mldists->nlevels++;

   nlevels = mldists->nlevels;

   assert(nlevels > 0);

   mldists->level_basestart[nlevels] = mldists->level_basestart[nlevels - 1] + nslots;
   mldists->level_targetstart[nlevels] = mldists->level_targetstart[nlevels - 1] + nslots * nslottargets;

   mldists->level_ntargets[nlevels - 1] = nslottargets;

   assert(extreduce_mldistsEmptySlotExists(mldists));

   mldistsTopLevelUnset(mldists);
}


/** adds dummy level */
void extreduce_mldistsLevelAddAndCloseEmpty(
   int                   nslottargets,       /**< number of targets per slot */
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   assert(mldists);

   extreduce_mldistsLevelAddTop(1, nslottargets, mldists);
   extreduce_mldistsEmptySlotSetBase(STP_MLDISTS_ID_EMPTY, mldists);
   extreduce_mldistsEmptySlotSetFilled(mldists);
   extreduce_mldistsLevelCloseTop(mldists);
}


/** adds root level of slots */
void extreduce_mldistsLevelAddAndCloseRoot(
   int                   base,               /**< the base */
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(base >= 0);

   extreduce_mldistsLevelAddTop(1, 0, mldists);
   extreduce_mldistsEmptySlotSetBase(base, mldists);
   extreduce_mldistsEmptySlotSetFilled(mldists);
   extreduce_mldistsLevelCloseTop(mldists);
}


/** closes the top level for further extensions */
void extreduce_mldistsLevelCloseTop(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   assert(!extreduce_mldistsIsEmpty(mldists));

   if( extreduce_mldistsEmptySlotExists(mldists) )
   {
      const int nlevels = mldists->nlevels;
      const int emptyslot = mldists->emptyslot_number;
      int nslottargets;

      assert(nlevels >= 1);
      assert(emptyslot >= 0);

      nslottargets = mldists->level_ntargets[nlevels - 1];

      assert(nslottargets >= 1);

      mldists->level_basestart[nlevels] = mldists->level_basestart[nlevels - 1] + emptyslot;
      mldists->level_targetstart[nlevels] = mldists->level_targetstart[nlevels - 1] + emptyslot * nslottargets;

      mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;
   }

   assert(!extreduce_mldistsEmptySlotExists(mldists));
}



/** reopens top level for one further extensions */
void extreduce_mldistsLevelReopenTop(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   const int nlevels = mldists->nlevels;
   const int nslottargets = mldists->level_ntargets[nlevels - 1];

   assert(nlevels > 0);
   assert(!extreduce_mldistsIsEmpty(mldists));
   assert(!extreduce_mldistsEmptySlotExists(mldists));
   assert(mldists->nlevels <= mldists->maxnlevels);

   mldists->emptyslot_number = mldists->level_basestart[nlevels] - mldists->level_basestart[nlevels - 1];
   mldists->level_basestart[nlevels]++;
   mldists->level_targetstart[nlevels] += nslottargets;

   assert(extreduce_mldistsEmptySlotExists(mldists));

#ifndef NDEBUG
   {
      const int targetend = mldists->level_targetstart[nlevels];
      mldists->base_ids[mldists->level_basestart[nlevels] - 1] = STP_MLDISTS_ID_UNSET;

      for( int i = targetend - nslottargets; i < targetend; ++i )
         mldists->target_dists[i] = STP_MLDISTS_DIST_UNSET;

      if( mldists->target_withids )
      {
         assert(mldists->target_ids);

         for( int i = targetend - nslottargets; i < targetend; ++i )
            mldists->target_ids[i] = STP_MLDISTS_ID_UNSET;
      }
   }
#endif
}


/** removes top level of slots */
void extreduce_mldistsLevelRemoveTop(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(!extreduce_mldistsEmptySlotExists(mldists));
   assert(mldists->nlevels > 0);

   mldistsTopLevelUnset(mldists);

   mldists->nlevels--;

   assert(!extreduce_mldistsEmptySlotExists(mldists));
}


/** removes top level of slots, which is not yet closed */
void extreduce_mldistsLevelRemoveTopNonClosed(
   MLDISTS*              mldists             /**< multi-level distances */
)
{
   assert(mldists);
   assert(extreduce_mldistsEmptySlotExists(mldists));
   assert(mldists->nlevels > 0);

   mldistsTopLevelUnset(mldists);

   mldists->nlevels--;
   mldists->emptyslot_number = MLDISTS_EMPTYSLOT_NONE;

   assert(!extreduce_mldistsEmptySlotExists(mldists));
}


/** gets number of targets (per slots) for given level */
int extreduce_mldistsLevelNTargets(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level               /**< level */
)
{
   assert(mldists);
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(mldists->level_ntargets[level] >= 0 && mldists->level_ntargets[level] <= mldists->level_maxntargets);

   return mldists->level_ntargets[level];
}

/** gets number of targets (per slots) for top level */
int extreduce_mldistsLevelNTopTargets(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   assert(mldists);

   return extreduce_mldistsLevelNTargets(mldists, mldistsGetTopLevel(mldists));
}


/** gets number of slots for given level */
int extreduce_mldistsLevelNSlots(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level               /**< level */
)
{
   assert(mldists);
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(mldists->level_basestart[level + 1] - mldists->level_basestart[level] >= 0);
   assert(mldists->level_basestart[level + 1] - mldists->level_basestart[level] <= mldists->level_maxnslots);

   return (mldists->level_basestart[level + 1] - mldists->level_basestart[level]);
}


/** gets number of slots for top level */
int extreduce_mldistsTopLevelNSlots(
   const MLDISTS*        mldists             /**< multi-level distances */
)
{
   return (extreduce_mldistsLevelNSlots(mldists, extreduce_mldistsTopLevel(mldists)));
}


/** gets number of levels */
int extreduce_mldistsNlevels(
   const MLDISTS*        mldists            /**< multi-level distances */
)
{
   assert(mldists);
   assert(mldists->nlevels >= 0);

   return (mldists->nlevels);
}


/** gets top level */
int extreduce_mldistsTopLevel(
   const MLDISTS*        mldists            /**< multi-level distances */
)
{
   assert(mldists);
   assert(mldists->nlevels >= 1);

   return (mldists->nlevels - 1);
}


/** gets top level bases*/
const int* extreduce_mldistsTopLevelBases(
   const MLDISTS*        mldists            /**< multi-level distances */
)
{
   const int toplevel = mldistsGetTopLevel(mldists);
   const int basestart_pos = mldistsGetPosBasesStart(mldists, toplevel);

   return &(mldists->base_ids[basestart_pos]);
}


/** is the base contained in a slot of the given level? */
SCIP_Bool extreduce_mldistsLevelContainsBase(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base id */
)
{
   assert(mldists);
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));
   assert(baseid >= 0);

   return (mldistsGetPosBase(mldists, level, baseid) != -1);
}


/** gets targets ids */
const int* extreduce_mldistsTargetIds(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base */
)
{
   const int targetpos = mldistsGetPosTargetsStart(mldists, level, baseid);

   assert(mldists->target_withids);
   assert(mldists->target_ids);

   return &(mldists->target_ids[targetpos]);
}


/** gets targets distances */
const SCIP_Real* extreduce_mldistsTargetDists(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid              /**< the base */
)
{
   int targetpos;

   assert(mldists);
   assert(level >= 0 && baseid >= 0);

   targetpos = mldistsGetPosTargetsStart(mldists, level, baseid);

   return &(mldists->target_dists[targetpos]);
}


/** Gets targets ids */
const int* extreduce_mldistsTopTargetIds(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   baseid              /**< the base */
)
{
   return extreduce_mldistsTargetIds(mldists, mldistsGetTopLevel(mldists), baseid);
}


/** gets targets distances */
const SCIP_Real* extreduce_mldistsTopTargetDists(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   baseid              /**< the base */
)
{
   const int targetpos = mldistsGetPosTargetsStart(mldists, mldistsGetTopLevel(mldists), baseid);

   return &(mldists->target_dists[targetpos]);
}


/** gets (one!) target distance for given target ID and base ID */
SCIP_Real extreduce_mldistsTopTargetDist(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   baseid,             /**< the base */
   int                   targetid            /**< the identifier */
)
{
   int targetpos;

   assert(mldists);
   assert(baseid >= 0 && targetid >= 0);

   targetpos = mldistsGetPosTargets(mldists, mldistsGetTopLevel(mldists), baseid, targetid);

   assert(GE(mldists->target_dists[targetpos], 0.0));

   return (mldists->target_dists[targetpos]);
}

/** gets (one!) target distance for given level, target ID and base ID */
SCIP_Real extreduce_mldistsTargetDist(
   const MLDISTS*        mldists,            /**< multi-level distances */
   int                   level,              /**< level */
   int                   baseid,             /**< the base */
   int                   targetid            /**< the identifier */
)
{
   int targetpos;

   assert(mldists);
   assert(baseid >= 0 && targetid >= 0);
   assert(level >= 0 && level <= mldistsGetTopLevel(mldists));

   targetpos = mldistsGetPosTargets(mldists, level, baseid, targetid);

   assert(GE(mldists->target_dists[targetpos], 0.0));

   return (mldists->target_dists[targetpos]);
}

/**@} */
