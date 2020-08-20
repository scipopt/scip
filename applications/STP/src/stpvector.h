/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2019 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scip.zib.de.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   stpvector.h
 * @brief  header only, simple implementation of an STL like vector
 * @author Daniel Rehfeldt
 *
 * STP vector (basically a much simplified replication of the C++ std::vector)
 * Usage: define with
 * STP_Vectype(template) myvec = NULL;
 *
 * for example:
 * STP_Vectype(int) myvec = NULL;
 *
 * Vector needs to be freed!
 *
 */


/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef APPLICATIONS_STP_SRC_STPVECTOR_H_
#define APPLICATIONS_STP_SRC_STPVECTOR_H_

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

#define STP_Vectype(type) type*


/** internal method */
#define vecinternalSetCapacity(vec, size)    \
   do                                        \
   {                                         \
      assert(size >= 0);                     \
      assert(vec);                           \
      ((int*) (vec))[-2] = (size);           \
   } while( 0 )


/** internal method */
#define vecinternalSetSize(vec, size)        \
   do                                        \
   {                                         \
      assert(size >= 0);                     \
      assert(vec);                           \
      ((int*) (vec))[-1] = (size);           \
   } while( 0 )


/** internal method */
#define vecinternalIncrementSize(vec)        \
   do                                        \
   {                                         \
      assert(vec);                           \
      ((int*) (vec))[-1]++;                  \
   } while( 0 )


/** internal method */
#define vecinternalDecrementSize(vec)        \
   do                                        \
   {                                         \
      assert(vec);                           \
      assert(((int*) (vec))[-1] >= 1);       \
      ((int*) (vec))[-1]--;                  \
   } while( 0 )


/** internal method */
#define vecinternalGetCapacity(vec)    \
   ((vec) ? ((int*) (vec))[-2] : (int) 0)


/** internal method */
#define vecinternalGetSize(vec)      \
   ((vec) ? ((int*) (vec))[-1] : (int) 0)


/** internal method */
#define vecinternalIncreaseCapacity(scip, vec, cap)                                 \
   do                                                                               \
   {                                                                                \
      const int nbytes_new = (cap) * sizeof(*(vec)) + (sizeof(int) * 2);            \
      assert(nbytes_new >= ((int) sizeof(int) * 2));                                \
      if( !(vec) )                                                                  \
      {                                                                             \
         char* p;                                                                   \
         SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &p, nbytes_new) );             \
         (vec) = (void*) (&p[sizeof(int) * 2]);                                     \
         vecinternalSetCapacity((vec), (cap));                                      \
         vecinternalSetSize((vec), 0);                                              \
      }                                                                             \
      else                                                                          \
      {                                                                             \
         const int nbytes_old = (vecinternalGetCapacity((vec))) * sizeof(*(vec)) + (sizeof(int) * 2); \
         char* p = &((char*) (vec))[(int) sizeof(int) * (-2)];                                        \
         assert(nbytes_old < nbytes_new);                                                             \
         SCIP_CALL_ABORT( SCIPreallocBlockMemoryArray(scip, &p, nbytes_old, nbytes_new) );            \
         (vec) = (void*)(&p[sizeof(int) * 2]);                                                        \
         vecinternalSetCapacity((vec), (cap));                                                        \
      }                                                                                               \
   } while( 0 )


/** gets capacity of the vector */
#define StpVecGetCapacity(vec)      \
    (vecinternalGetCapacity(vec))


/** removes all elements from vector */
#define StpVecClear(vec)   \
   vecinternalSetSize(vec, (int)0)


/** gets size of the vector (number of elements) */
#define StpVecGetSize(vec)      \
   vecinternalGetSize(vec)


/** is the vector empty? */
#define StpVecIsEmpty(vec)    \
   (vecinternalGetSize(vec) == 0)


/** frees vector */
#define StpVecFree(scip, vec)                                                                          \
   do                                                                                                  \
   {                                                                                                   \
      if( vec )                                                                                        \
      {                                                                                                \
         char* p = &((char*) (vec))[(int) sizeof(int) * (-2)];                                         \
         const int nbytes = (vecinternalGetCapacity((vec))) * sizeof(*(vec)) + (sizeof(int) * 2);      \
         SCIPfreeBlockMemoryArray(scip, &p, nbytes);                                                   \
      }                                                                                                \
   } while( 0 )


/** adds (appends) element */
#define StpVecPushBack(scip, vec, value)                                 \
   do                                                                    \
   {                                                                     \
      const int cap = vecinternalGetCapacity(vec);                       \
      const int size = vecinternalGetSize(vec);                          \
      if( cap <= size )                                                  \
      {                                                                  \
         vecinternalIncreaseCapacity(scip, (vec), (cap == 0) ? 2 : cap * 2); \
      }                                                                  \
      vec[size] = (value);                                               \
      vecinternalIncrementSize((vec));                                   \
   } while( 0 )


/** remove last element */
#define StpVecPopBack(vec)                           \
   vecinternalDecrementSize((vec))

/** reserves space */
#define StpVecReserve(scip, vec, size)                                   \
   do                                                                    \
   {                                                                     \
      const int cap = vecinternalGetCapacity(vec);                       \
      if( cap < size )                                                   \
      {                                                                  \
         vecinternalIncreaseCapacity(scip, (vec), size);                 \
      }                                                                  \
   } while( 0 )



#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_STPVECTOR_H_ */
