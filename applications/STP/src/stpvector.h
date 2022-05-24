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
#define vecinternalSetCapacity(vec, _size_)    \
   do                                        \
   {                                         \
      assert(_size_ >= 0);                     \
      assert(vec);                           \
      ((int*) (vec))[-2] = (_size_);           \
   } while( 0 )


/** internal method */
#define vecinternalSetSize(vec, _size_)        \
   do                                        \
   {                                         \
      assert(_size_ >= 0);                     \
      assert(vec);                           \
      ((int*) (vec))[-1] = (_size_);           \
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
#define vecinternalComputeNBytes(size, vec)      \
   (size * sizeof(*(vec)) + (sizeof(int) * 2))


/** internal method */
#define vecinternalIncreaseCapacity(scip, vec, cap)                                 \
   do                                                                               \
   {                                                                                \
      const int _nbytes_new_ = vecinternalComputeNBytes((cap), (vec));                \
      assert(_nbytes_new_ >= ((int) sizeof(int) * 2));                                \
      if( !(vec) )                                                                  \
      {                                                                             \
         char* _p_;                                                                   \
         SCIP_CALL_ABORT( SCIPallocBlockMemoryArray((scip), &_p_, _nbytes_new_) );             \
         (vec) = (void*) (&_p_[sizeof(int) * 2]);                                     \
         vecinternalSetCapacity((vec), (cap));                                      \
         vecinternalSetSize((vec), 0);                                              \
      }                                                                             \
      else                                                                          \
      {                                                                             \
         const int _cap_old_ = vecinternalGetCapacity((vec));                      \
         const int _nbytes_old_ = vecinternalComputeNBytes(_cap_old_, (vec));      \
         char* _p_ = &((char*) (vec))[(int) sizeof(int) * (-2)];                                        \
         assert(_nbytes_old_ < _nbytes_new_);                                                             \
         SCIP_CALL_ABORT( SCIPreallocBlockMemoryArray((scip), &_p_, _nbytes_old_, _nbytes_new_) );            \
         (vec) = (void*)(&_p_[sizeof(int) * 2]);                                                        \
         vecinternalSetCapacity((vec), (cap));                                                        \
      }                                                                                               \
   } while( 0 )


/** gets _cap_acity of the vector */
#define StpVecGetcapacity(vec)      \
    (vecinternalGetCapacity(vec))


/** removes all elements from vector */
#define StpVecClear(vec)   \
   vecinternalSetSize(vec, (int)0)


/** gets _size_ of the vector (number of elements) */
#define StpVecGetSize(vec)      \
   vecinternalGetSize(vec)

/** gets top element of the vector */
#define StpVecTop(vec)      \
    vec[StpVecGetSize(vec) - 1]


/** is the vector empty? */
#define StpVecIsEmpty(vec)    \
   (vecinternalGetSize(vec) == 0)


/** frees vector */
#define StpVecFree(scip, vec)                                                                          \
   do                                                                                                  \
   {                                                                                                   \
      if( (vec) )                                                                                        \
      {                                                                                                \
         char* _p_ = &((char*) (vec))[(int) sizeof(int) * (-2)];                                         \
         const int _nbytes_ = vecinternalComputeNBytes(vecinternalGetCapacity((vec)), (vec));         \
         SCIPfreeBlockMemoryArray((scip), &_p_, _nbytes_);                                                   \
         (vec) = NULL; \
      }                                                                                                \
   } while( 0 )


/** adds (appends) element */
#define StpVecPushBack(scip, vec, value)                                 \
   do                                                                    \
   {                                                                     \
      const int _cap_ = vecinternalGetCapacity((vec));                       \
      const int _size_ = vecinternalGetSize((vec));                          \
      if( _cap_ <= _size_ )                                                  \
      {                                                                  \
         vecinternalIncreaseCapacity((scip), (vec), ((_cap_ == 0) ? 2 : _cap_ * 2)); \
      }                                                                  \
      vec[_size_] = (value);                                               \
      vecinternalIncrementSize((vec));                                   \
   } while( 0 )


/** remove last element */
#define StpVecPopBack(vec)                           \
   vecinternalDecrementSize((vec))

/** reserves space */
#define StpVecReserve(scip, vec, _size_)                                   \
   do                                                                    \
   {                                                                     \
      const int _cap_ = vecinternalGetCapacity(vec);                       \
      if( _cap_ < _size_ )                                                   \
      {                                                                  \
         vecinternalIncreaseCapacity(scip, (vec), _size_);                 \
      }                                                                  \
   } while( 0 )



#ifdef __cplusplus
}
#endif


#endif /* APPLICATIONS_STP_SRC_STPVECTOR_H_ */
