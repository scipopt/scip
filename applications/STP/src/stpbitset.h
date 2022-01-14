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

/**@file   stpbitset.h
 * @brief  header only, simple implementation of a bitset
 * @author Daniel Rehfeldt
 *
 * Implements a simple bitset. Should be set to NULL if undefined.
 * NOTE: for efficiency reasons the bitset type is based on an STP_Vector and thus
 * uses a non-SCIP-standard allocation method. In this way we avoid indirections, because
 * we directly access the raw array.
 * todo if too slow use extra fixed-size vector that uses standard memory allocs, or even
 * cache-aligned mallocs
 *
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_STPBITSET_H_
#define APPLICATIONS_STP_SRC_STPBITSET_H_

#include "scip/scip.h"
#include "stpvector.h"

#ifdef __cplusplus
extern "C" {
#endif


#define STP_Bitset STP_Vectype(uint64_t)


/** todo: do that more efficiently; currently just a text-book implementation.
 *  probably want to have case distinction for compilers and use intrinsics
 *  at least for gcc, clang, and intel */
static
inline int bitsetinternalPopcount(
   uint64_t              number              /**< to get popcount for */
   )
{
   const uint64_t stpbit_m1 = 0x5555555555555555;
   const uint64_t stpbit_m2 = 0x3333333333333333;
   const uint64_t stpbit_m4 = 0x0f0f0f0f0f0f0f0f;
   const uint64_t stpbit_powseries = 0x0101010101010101;
   uint64_t n = number;

   n -= (n >> 1) & stpbit_m1;
   n = (n & stpbit_m2) + ((n >> 2) & stpbit_m2);
   n = (n + (n >> 4)) & stpbit_m4;

   return (int) ((n * stpbit_powseries) >> 56);
}


/*
 * Interface methods
 */

/** initializes clean (all-0) bitset and returns it */
static
inline STP_Bitset stpbitset_new(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxnbits            /**< size of bitset */
   )
{
   STP_Bitset bitset = NULL;
   const int size = (maxnbits + 63) / 64;

   assert(maxnbits > 0);
   assert(size > 0);

   StpVecReserve(scip, bitset, size);

   // todo do that more efficiently with new Stp_Vector reserve method
   for( int i = 0; i < size; i++ )
   {
      StpVecPushBack(scip, bitset, (uint64_t) 0);
   }

   return bitset;
}


/** frees */
static
inline void stpbitset_free(
   SCIP*                 scip,               /**< SCIP data structure */
   STP_Bitset*           bitset              /**< bitset pointer */
   )
{
   assert(scip && bitset);
   assert(*bitset);

   StpVecFree(scip, *bitset);

   assert(NULL == *bitset);
}


/** gets number of bits that can be stored  */
static
inline int stpbitset_getCapacity(
   STP_Bitset           bitset              /**< bitset */
   )
{
   assert(bitset);
   assert(StpVecGetcapacity(bitset) == StpVecGetSize(bitset));

   return (StpVecGetcapacity(bitset) * 64);
}


/** sets bit to TRUE (1) */
static
inline void stpbitset_setBitTrue(
   STP_Bitset           bitset,              /**< bitset */
   int                  index                /**< bit index */
   )
{
   assert(bitset);
   assert(0 <= index && index < stpbitset_getCapacity(bitset));

   bitset[index / 64] |= (uint64_t) 1 << (((uint64_t) index) & 63);
}


/** sets bit to FALSE (0) */
static
inline void stpbitset_setBitFalse(
   STP_Bitset           bitset,              /**< bitset */
   int                  index                /**< bit index */
   )
{
   assert(bitset);
   assert(0 <= index && index < stpbitset_getCapacity(bitset));

   bitset[index / 64] &= ~((uint64_t) 1 << (((uint64_t) index) & 63));
}


/** are given bitsets compatible? */
static
inline SCIP_Bool stpbitset_setsAreCompatible(
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2              /**< bitset */
   )
{
   assert(bitset1 && bitset2);

   if( stpbitset_getCapacity(bitset1) != stpbitset_getCapacity(bitset2) )
   {
      return FALSE;
   }

   assert(StpVecGetSize(bitset1) == StpVecGetSize(bitset2));

   return TRUE;
}


/** do bitsets (of same size) intersect? */
static
inline SCIP_Bool stpbitset_haveIntersection(
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2              /**< bitset */
   )
{
   const int vecsize = StpVecGetSize(bitset1);

   assert(bitset1 && bitset2);
   assert(stpbitset_setsAreCompatible(bitset1, bitset2));
   assert(vecsize > 0);

   for( int i = 0; i < vecsize; i++ )
   {
      if( (bitset1[i] & bitset2[i]) != 0 )
         return TRUE;
   }

   return FALSE;
}


/** are bitsets the same? */
static
inline SCIP_Bool stpbitset_areEqual(
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2              /**< bitset */
   )
{
   const int vecsize = StpVecGetSize(bitset1);

   assert(bitset1 && bitset2);
   assert(stpbitset_setsAreCompatible(bitset1, bitset2));
   assert(vecsize > 0);

   if( memcmp(bitset1, bitset2, sizeof(bitset1[0]) * vecsize) != 0 )
   {
      return FALSE;

   }

   return TRUE;
}


/** is bit at given index set? */
static
inline SCIP_Bool stpbitset_bitIsTrue(
   STP_Bitset           bitset,              /**< bitset */
   int                  index                /**< bit index */
   )
{
   const uint64_t i = (uint64_t) index;

   assert(bitset);
   assert(0 <= index && index < stpbitset_getCapacity(bitset));

   if( ((bitset[i / 64] >> ((i & 63)) & (uint64_t) 1)) == 0 )
   {
      return FALSE;
   }
   else
   {
      return TRUE;
   }
}



/** gets number of 1-bits */
static
inline int stpbitset_getPopcount(
   STP_Bitset           bitset               /**< bitset */
   )
{
   int popcount = 0;
   const int vecsize = StpVecGetSize(bitset);

   assert(vecsize > 0);

   for( int i = 0; i < vecsize; i++ )
   {
      popcount += bitsetinternalPopcount(bitset[i]);
   }

   return popcount;
}


/** gets copy */
static
inline STP_Bitset stpbitset_newCopy(
   SCIP*                scip,                /**< SCIP data structure */
   STP_Bitset           bitsetOrg            /**< bitset */
   )
{
   const int vecsize = StpVecGetSize(bitsetOrg);
   STP_Bitset bitset = stpbitset_new(scip, stpbitset_getCapacity(bitsetOrg));

   assert(stpbitset_setsAreCompatible(bitsetOrg, bitset));
   assert(vecsize > 0);

   BMScopyMemoryArray(bitset, bitsetOrg, vecsize);

   return bitset;
}


/** sets  AND bitset */
static
inline void stpbitset_and(
   SCIP*                scip,                /**< SCIP data structure */
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2,             /**< bitset */
   STP_Bitset           bitsetAnd            /**< bitset to set */
   )
{
   const int vecsize = StpVecGetSize(bitsetAnd);

   assert(stpbitset_setsAreCompatible(bitset1, bitset2));
   assert(stpbitset_setsAreCompatible(bitsetAnd, bitset2));
   assert(vecsize > 0);

   for( int i = 0; i < vecsize; i++ )
   {
      bitsetAnd[i] = (bitset1[i] & bitset2[i]);
   }
}


/** sets OR bitset */
static
inline void stpbitset_or(
   SCIP*                scip,                /**< SCIP data structure */
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2,             /**< bitset */
   STP_Bitset           bitsetOr             /**< bitset to set */
   )
{
   const int vecsize = StpVecGetSize(bitsetOr);

   assert(stpbitset_setsAreCompatible(bitset1, bitset2));
   assert(stpbitset_setsAreCompatible(bitsetOr, bitset2));
   assert(vecsize > 0);

   for( int i = 0; i < vecsize; i++ )
   {
      bitsetOr[i] = (bitset1[i] | bitset2[i]);
   }
}


/** gets new AND bitset */
static
inline STP_Bitset stpbitset_newAnd(
   SCIP*                scip,                /**< SCIP data structure */
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2              /**< bitset */
   )
{
   STP_Bitset bitset = stpbitset_new(scip, stpbitset_getCapacity(bitset1));

   assert(stpbitset_setsAreCompatible(bitset1, bitset2));

   stpbitset_and(scip, bitset1, bitset2, bitset);

   return bitset;
}


/** gets new OR bitset */
static
inline STP_Bitset stpbitset_newOr(
   SCIP*                scip,                /**< SCIP data structure */
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2              /**< bitset */
   )
{
   STP_Bitset bitset = stpbitset_new(scip, stpbitset_getCapacity(bitset1));

   assert(stpbitset_setsAreCompatible(bitset1, bitset2));

   stpbitset_or(scip, bitset1, bitset2, bitset);

   return bitset;
}


/** gets new XOR bitset */
static
inline STP_Bitset stpbitset_newXor(
   SCIP*                scip,                /**< SCIP data structure */
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2              /**< bitset */
   )
{
   const int vecsize = StpVecGetSize(bitset1);
   STP_Bitset bitset = stpbitset_new(scip, stpbitset_getCapacity(bitset1));

   assert(stpbitset_setsAreCompatible(bitset1, bitset2));
   assert(StpVecGetSize(bitset) == StpVecGetSize(bitset1));
   assert(vecsize > 0);

   for( int i = 0; i < vecsize; i++ )
   {
      bitset[i] = (bitset1[i] ^ bitset2[i]);
   }

   return bitset;
}


/** gets new Not bitset */
static
inline STP_Bitset stpbitset_newNot(
   SCIP*                scip,                /**< SCIP data structure */
   STP_Bitset           bitset1,             /**< bitset for not */
   int                  size                 /**< size to apply not operation to */
   )
{
   const int vecsize_cut = size / 64;
   STP_Bitset bitset = stpbitset_new(scip, stpbitset_getCapacity(bitset1));

   assert(StpVecGetSize(bitset) == StpVecGetSize(bitset1));
   assert(0 <= size && size <= stpbitset_getCapacity(bitset));

   for( int i = 0; i < vecsize_cut; i++ )
   {
      bitset[i] = ~(bitset1[i]);
   }

   if( (vecsize_cut * 64) != size )
   {
      const int nremaining = size - vecsize_cut * 64;
      assert(nremaining > 0);

      bitset[vecsize_cut] = ~bitset1[vecsize_cut];

      for( int i = nremaining; i < 64; i++ )
      {
         bitset[vecsize_cut] &= ~((uint64_t) 1 << ((uint64_t) i));
      }
   }

   assert((size == stpbitset_getCapacity(bitset))
       || !stpbitset_bitIsTrue(bitset, size));

   return bitset;
}


/** initializes all TRUE bitset and returns it
 *  todo more efficiently */
static
inline STP_Bitset stpbitset_newAllTrue(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   maxnbits            /**< size of bitset */
   )
{
   STP_Bitset bitset = stpbitset_new(scip, maxnbits);
   const int cap = stpbitset_getCapacity(bitset);

   assert(cap > 0);

   for( int i = 0; i < cap; i++ )
   {
      assert(!stpbitset_bitIsTrue(bitset, i));
      stpbitset_setBitTrue(bitset, i);
   }

   return bitset;
}


/** prints bitset (0,1) ... for debugging */
static
inline void stpbitset_print(
   STP_Bitset           bitset               /**< bitset to print*/
   )
{
   const int cap = stpbitset_getCapacity(bitset);

   assert(cap > 0);

   printf("bit-set (cap=%d): ", cap);

   for( int i = 0; i < cap; i++ )
   {
      printf("%d", stpbitset_bitIsTrue(bitset, i));
   }

   printf("\n");
}


/** bitset1 > bitset2? */
static
inline SCIP_Bool stpbitset_GT(
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2              /**< bitset */
   )
{
   const int vecsize = StpVecGetSize(bitset1);
   int ret;

   assert(stpbitset_setsAreCompatible(bitset1, bitset2));
   assert(vecsize > 0);

   ret = memcmp(bitset1, bitset2, sizeof(bitset1[0]) * vecsize);

   return (ret > 0);
}


/** bitset1 < bitset2? */
static
inline SCIP_Bool stpbitset_LT(
   STP_Bitset           bitset1,             /**< bitset */
   STP_Bitset           bitset2              /**< bitset */
   )
{
   const int vecsize = StpVecGetSize(bitset1);
   int ret;

   assert(stpbitset_setsAreCompatible(bitset1, bitset2));
   assert(vecsize > 0);

   ret = memcmp(bitset1, bitset2, sizeof(bitset1[0]) * vecsize);

   return (ret < 0);
}


#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_STPBITSET_H_ */
