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

/**@file   random.c
 * @brief  random number generator
 * @author Jakob Witzig
 */

#include "scip/random.h"


#define DEFAULT_SEED 123456789;
#define DEFAULT_XOR  362436000;
#define DEFAULT_MWC  521288629;
#define DEFAULT_CST  7654321;


/** initialize the random number generator with a given start seed */
static
void randomInitialze(
   SCIP_RANDGEN*         randgen,
   unsigned int          initseed
   )
{
   assert(randgen != NULL);
   assert(initseed > 0);

   randgen->seed = (uint32_t)DEFAULT_SEED;
   randgen->seed += initseed;
   /* this is a special case to avoid zero after over flowing */
   if( randgen->seed == 0 )
      randgen->seed = initseed;
   assert(randgen->seed > 0);

   randgen->xor_seed = (uint32_t)DEFAULT_XOR;
   randgen->xor_seed += initseed;
   /* this is a special case to avoid zero after over flowing */
   if( randgen->xor_seed == 0 )
      randgen->xor_seed = initseed;
   assert(randgen->xor_seed > 0);

   randgen->mwc_seed = (uint32_t)DEFAULT_MWC;
   randgen->mwc_seed += initseed;

   randgen->cst_seed = (uint32_t)DEFAULT_CST;
   randgen->cst_seed += initseed;
   /* this is a special case to avoid zero after over flowing */
   if( randgen->mwc_seed == 0 && randgen->cst_seed == 0 )
   {
      randgen->mwc_seed = (uint32_t)DEFAULT_MWC;
      randgen->mwc_seed -= initseed;

      randgen->cst_seed = (uint32_t)DEFAULT_CST;
      randgen->cst_seed -= initseed;
   }
   assert(randgen->mwc_seed > 0 || randgen->cst_seed > 0);

   return;
}

/** creates and initialzes a random number generator */
SCIP_RETCODE SCIPrandomCreate(
   SCIP_RANDGEN**        randnumgen,
   BMS_BLKMEM*           blkmem,
   unsigned int          initialseed
   )
{
   assert(randnumgen != NULL);

   SCIP_ALLOC( BMSallocBlockMemory(blkmem, randnumgen) );

   randomInitialze((*randnumgen), initialseed);

   return SCIP_OKAY;
}

/** creates and initialzes a random number generator */
void SCIPrandomFree(
   SCIP_RANDGEN**        randnumgen,
   BMS_BLKMEM*           blkmem
   )
{
   assert(randnumgen != NULL);
   assert((*randnumgen) != NULL);

   BMSfreeBlockMemory(blkmem, randnumgen);

   return;
}

/** returns a random number between 0 and INT_MAX
 *
 *  implementation of KISS random number generator developed by George Marsaglia.
 *  KISS is combination of three different random number generators:
 *   - Linear congruential generator
 *   - Xorshift
 *   - Lag-1 Multiply-with-carry
 *
 *  KISS has a period of 2^123 and passes all statistical test part of BigCrush-Test of TestU01 [1].
 *
 *  [1] http://dl.acm.org/citation.cfm?doid=1268776.1268777
 */
static
int getRand(
   SCIP_RANDGEN*         randgen             /**< random number generator */
   )
{
   unsigned long t;

   /* linear congruential */
   randgen->seed = randgen->seed * (SCIP_Longint)1103515245 + 12345;

   /* Xorshift */
   randgen->xor_seed ^= (randgen->xor_seed << 13);
   randgen->xor_seed ^= (randgen->xor_seed >> 17);
   randgen->xor_seed ^= (randgen->xor_seed << 5);

   /* Multiple-with-carry */
   t = 698769069ULL * randgen->mwc_seed + randgen->cst_seed;
   randgen->cst_seed = t >> 32;
   randgen->mwc_seed = (unsigned int) t;

   return (int)((randgen->seed + randgen->xor_seed + randgen->mwc_seed) % INT_MAX);
}

/** returns a random integer between minrandval and maxrandval */
int SCIPrandomGetInt(
   SCIP_RANDGEN*         randgen,            /**< random number generator */
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval          /**< maximal value to return */
   )
{
   SCIP_Real randnumber;

   randnumber = (SCIP_Real)getRand(randgen)/(INT_MAX+1.0);
   assert(randnumber >= 0.0);
   assert(randnumber < 1.0);

   /* we multiply minrandval and maxrandval separately by randnumber in order to avoid overflow if they are more than INT_MAX
    * apart
    */
   return (int) (minrandval*(1.0 - randnumber) + maxrandval*randnumber + randnumber);
}

/** returns a random real between minrandval and maxrandval */
SCIP_Real SCIPrandomGetReal(
   SCIP_RANDGEN*         randgen,            /**< random number generator */
   SCIP_Real             minrandval,         /**< minimal value to return */
   SCIP_Real             maxrandval          /**< maximal value to return */
   )
{
   SCIP_Real randnumber;

   randnumber = (SCIP_Real)getRand(randgen)/(SCIP_Real)INT_MAX;
   assert(randnumber >= 0.0);
   assert(randnumber <= 1.0);

   /* we multiply minrandval and maxrandval separately by randnumber in order to avoid overflow if they are more than
    * SCIP_REAL_MAX apart
    */
   return minrandval*(1.0 - randnumber) + maxrandval*randnumber;
}

/** randomly shuffles parts of an integer array using the Fisher-Yates algorithm */
void SCIPrandomPermuteIntArray(
   SCIP_RANDGEN*         randgen,            /**< random number generator data */
   int*                  array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end                 /**< last index that should be subject to shuffling (array size for whole
                                               *   array)
                                               */
   )
{
   int tmp;
   int i;

   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 )
   {
      --end;

      /* get a random position into which the last entry should be shuffled */
      i = SCIPrandomGetInt(randgen, begin, end);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}

/** randomly shuffles parts of an array using the Fisher-Yates algorithm */
void SCIPrandomPermuteArray(
   SCIP_RANDGEN*         randgen,            /**< random number generator data */
   void**                array,              /**< array to be shuffled */
   int                   begin,              /**< first index that should be subject to shuffling (0 for whole array) */
   int                   end                 /**< last index that should be subject to shuffling (array size for whole
                                              *   array)
                                              */
   )
{
   void* tmp;
   int i;

   /* loop backwards through all elements and always swap the current last element to a random position */
   while( end > begin+1 )
   {
      end--;

      /* get a random position into which the last entry should be shuffled */
      i = SCIPrandomGetInt(randgen, begin, end);

      /* swap the last element and the random element */
      tmp = array[i];
      array[i] = array[end];
      array[end] = tmp;
   }
}

/** draws a random subset of disjoint elements from a given set of disjoint elements;
 *  this implementation is suited for the case that nsubelems is considerably smaller then nelems
 */
SCIP_RETCODE SCIPrandomGetSubset(
   SCIP_RANDGEN*         randgen,            /**< random number generator data */
   void**                set,                /**< original set, from which elements should be drawn */
   int                   nelems,             /**< number of elements in original set */
   void**                subset,             /**< subset in which drawn elements should be stored */
   int                   nsubelems           /**< number of elements that should be drawn and stored */
   )
{
   int i;
   int j;

   /* if both sets are of equal size, we just copy the array */
   if( nelems == nsubelems)
   {
      BMScopyMemoryArray(subset,set,nelems);
      return SCIP_OKAY;
   }

   /* abort, if size of subset is too big */
   if( nsubelems > nelems )
   {
      SCIPerrorMessage("Cannot create %d-elementary subset of %d-elementary set.\n", nsubelems, nelems);
      return SCIP_INVALIDDATA;
   }
#ifndef NDEBUG
   for( i = 0; i < nsubelems; i++ )
      for( j = 0; j < i; j++ )
         assert(set[i] != set[j]);
#endif

   /* draw each element individually */
   i = 0;
   while( i < nsubelems )
   {
      int r;

      r = SCIPrandomGetInt(randgen, 0, nelems-1);
      subset[i] = set[r];

      /* if we get an element that we already had, we will draw again */
      for( j = 0; j < i; j++ )
      {
         if( subset[i] == subset[j] )
         {
            --i;
            break;
         }
      }
      ++i;
   }
   return SCIP_OKAY;
}
