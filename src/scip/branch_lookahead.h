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

/**@file   branch_lookahead.h
 * @ingroup BRANCHINGRULES
 * @brief  lookahead branching rule TODO CS: expand the description
 * @author Christoph Schubert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_BRANCH_LOOKAHEAD_H__
#define __SCIP_BRANCH_LOOKAHEAD_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * This enum is used to represent whether an upper bound, lower bound or both are set for a variable.
 */
typedef enum
{
   BOUNDSTATUS_NONE = 0,
   BOUNDSTATUS_UPPERBOUND,
   BOUNDSTATUS_LOWERBOUND,
   BOUNDSTATUS_BOTH
} BOUNDSTATUS;

/**
 * This struct collect the bounds, which can be used in the root problem. Therefore it contains two types of bounds:
 * - The bounds that occur when both up and down branches of a variable after a first level branch are cutoff. In this case
 *   the whole first level branch can be added as a restriction.
 * - The bounds that occur implicitly while branching on the second level. See SupposedBoundData for more information.
 */
typedef struct
{
   SCIP_Real*            upperbounds;        /**< The current upper bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_UPPERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   SCIP_Real*            lowerbounds;        /**< The current lower bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_LOWERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   BOUNDSTATUS*          boundstatus;        /**< The current boundstatus for each active variable. Depending on this value
                                              *   the corresponding upperbound and lowerbound values are meaningful.*/
   int*                  boundedvars;        /**< Contains the var indices that have entries in the other arrays. This array
                                              *   may be used to only iterate over the relevant variables. */
   int                   nboundedvars;       /**< The length of the boundedvars array. */
} ValidBoundData;

/**
 * This struct collects the bounds, that are given implicitly on the second branching level.
 * Concrete: If a variable is regarded on both sides of the second level and is infeasible (in the same bound direction) on
 * both sides, the weaker bound can be applied.
 * Even more concrete: First level branching on variable x, second level branching on variable y (and may others). If the
 * constraint y <= 3 on the up branch of x and y <= 6 on the down branch of x are both infeasible, the y <= 3 bound can be
 * applied on the first level.
 */
typedef struct
{
   SCIP_Real*            upperbounds;        /**< The current upper bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_UPPERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   int*                  nupperboundupdates; /**< The number of times the corresponding upper bound was updated. */
   SCIP_Real*            lowerbounds;        /**< The current lower bound for each active variable. Only contains meaningful
                                              *   data, if the corresponding boundstatus entry is BOUNDSTATUS_LOWERBOUND or
                                              *   BOUNDSTATUS_BOTH. */
   int*                  nlowerboundupdates; /**< The number of times the corresponding lower bound was updated. */
   int*                  boundedvars;        /**< Contains the var indices that have entries in the other arrays. This array
                                              *   may be used to only iterate over the relevant variables. */
   int                   nboundedvars;       /**< The length of the boundedvars array. */
} SupposedBoundData;

typedef struct
{
   SCIP_VAR**            eithervars;
   SCIP_VAR**            othervars;
   int                   nentries;
   int                   memsize;
} BinaryBoundData;

static
SCIP_RETCODE allocValidBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   ValidBoundData**      validbounddata
)
{
   int ntotalvars;

   ntotalvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBuffer(scip, validbounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->upperbounds, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->lowerbounds, ntotalvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(*validbounddata)->boundstatus, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*validbounddata)->boundedvars, ntotalvars) );
   return SCIP_OKAY;
}

/**
 * Clears the given struct.
 * Assumptions:
 * - The boundstatus array was cleared, when the bounds were transferred to the valid bounds data structure.
 * - The upper-/lowerbounds arrays don't have to be reset, as these are only read in connection with the boundstatus array.
 * - The boundedvars array is only read in connection with the nboundedvars value, which will be set to 0.
 */
static
void initValidBoundData(
   ValidBoundData*       validbounddata      /*< The struct that should get cleared.*/
)
{
   validbounddata->nboundedvars = 0;
}

static
void freeValidBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   ValidBoundData**      validbounddata
)
{
   SCIPfreeBufferArray(scip, &(*validbounddata)->boundedvars);
   SCIPfreeCleanBufferArray(scip, &(*validbounddata)->boundstatus);
   SCIPfreeBufferArray(scip, &(*validbounddata)->lowerbounds);
   SCIPfreeBufferArray(scip, &(*validbounddata)->upperbounds);
   SCIPfreeBuffer(scip, validbounddata);
}

static
SCIP_RETCODE allocSupposedBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   SupposedBoundData**   supposedbounddata
)
{
   int ntotalvars;

   ntotalvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBuffer(scip, supposedbounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*supposedbounddata)->upperbounds, ntotalvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(*supposedbounddata)->nupperboundupdates, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*supposedbounddata)->lowerbounds, ntotalvars) );
   SCIP_CALL( SCIPallocCleanBufferArray(scip, &(*supposedbounddata)->nlowerboundupdates, ntotalvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*supposedbounddata)->boundedvars, ntotalvars) );
   return SCIP_OKAY;
}

/**
 * Clears the given struct.
 * Assumptions:
 * - The boundstatus array was cleared, when the bounds were transferred to the valid bounds data structure.
 * - The upper-/lowerbounds arrays don't have to be reset, as these are only read in connection with the boundstatus array.
 * - The boundedvars array is only read in connection with the nboundedvars value, which will be set to 0.
 */
static
void initSupposedBoundData(
   SupposedBoundData*    supposedbounddata   /*< The struct that should get cleared.*/
)
{
   supposedbounddata->nboundedvars = 0;
}

static
void freeSupposedBoundData(
   SCIP*                 scip,               /**< SCIP data structure */
   SupposedBoundData**   supposedbounddata
)
{
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->boundedvars);
   SCIPfreeCleanBufferArray(scip, &(*supposedbounddata)->nlowerboundupdates);
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->lowerbounds);
   SCIPfreeCleanBufferArray(scip, &(*supposedbounddata)->nupperboundupdates);
   SCIPfreeBufferArray(scip, &(*supposedbounddata)->upperbounds);
   SCIPfreeBuffer(scip, supposedbounddata);
}

static
SCIP_RETCODE allocBinaryBoundData(
   SCIP*                 scip,
   BinaryBoundData**     binarybounddata,
   int                   nentries
)
{
   SCIP_CALL( SCIPallocBuffer(scip, binarybounddata) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*binarybounddata)->eithervars, nentries) );
   SCIP_CALL( SCIPallocBufferArray(scip, &(*binarybounddata)->othervars, nentries) );
   (*binarybounddata)->memsize = nentries;
   return SCIP_OKAY;
}

static
void initBinaryBoundData(
   BinaryBoundData*      binarybounddata
)
{
   binarybounddata->nentries = 0;
}

static
SCIP_RETCODE addBinaryBoundEntry(
   SCIP*                 scip,
   BinaryBoundData*      container,
   SCIP_VAR*             eithervar,
   SCIP_VAR*             othervar
)
{
   int emptyindex = container->nentries;

   if( emptyindex == container->memsize )
   {
      /* calculate new size, that can at least hold the old number of entries + 1 for the new entry */
      int newmemsize = SCIPcalcMemGrowSize(scip, emptyindex + 1);
      SCIP_CALL( SCIPreallocBufferArray(scip, &container->eithervars, newmemsize) );
      SCIP_CALL( SCIPreallocBufferArray(scip, &container->othervars, newmemsize) );
      container->memsize = newmemsize;
   }

   container->eithervars[emptyindex] = eithervar;
   container->othervars[emptyindex] = othervar;
   container->nentries = emptyindex + 1;

   return SCIP_OKAY;
}

static
SCIP_Bool isBinaryBoundDataEmpty(
   BinaryBoundData*      container
)
{
   return container->nentries == 0;
}

static
void freeBinaryBoundData(
   SCIP*                 scip,
   BinaryBoundData**     binarybounddata
)
{
   SCIPfreeBufferArray(scip, &(*binarybounddata)->othervars);
   SCIPfreeBufferArray(scip, &(*binarybounddata)->eithervars);
   SCIPfreeBuffer(scip, binarybounddata);
}

/** creates the Lookahead branching rule and includes it in SCIP */
EXTERN
SCIP_RETCODE SCIPincludeBranchruleLookahead(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
