/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2024 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   rectlu_num.h
 * @brief  rectlu number type macros
 * @author David Applegate
 * @author Bill Cook
 * @author Sanjeeb Dash
 * @author Daniel Espinoza
 * @author Dan Steffy
 * @author Kati Wolter
 */
/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifdef SCIP_WITH_GMP
#include <gmp.h>
#endif

#ifndef __RECTLU_NUM_H
#define __RECTLU_NUM_H

#ifdef SCIP_WITH_GMP
#define QSnum_type mpq_t

/*lint --e(160)*/

/* allocates array with size elements of QSnum_type */
QSnum_type* QSnum_AllocArray(int size);

/* frees array of QSnum_type */
void QSnum_FreeArray(QSnum_type* ea, int size);


/* initializes an individual element of QSnum_type */
#define QSnum_Init(a) mpq_init(a)

/* clears an individual element of QSnum_type */
#define QSnum_Clear(a) mpq_clear(a)

/*
 * basic arithmetic operation definitions for QSnum_type
 */
#define QSnum_Copy(a,b) mpq_set(a,b)              /* a = b   */
#define QSnum_CopyDbl(a,b) mpq_set_d(a,b)         /* a = b (for double b) */
#define QSnum_CopyAbs(a,b) mpq_abs(a,b)           /* a = |b| */
#define QSnum_CopyNeg(a,b) mpq_neg(a,b)           /* a = -b  */
#define QSnum_Sign(a) mpq_neg(a,a)                /* a = -a  */
#define QSnum_CopyFrac(a,b,c) mpq_div(a,b,c)      /* a = b/c */
#define QSnum_CopyDiv(a,b) mpq_div(a,a,b)         /* a = a/b */
#define QSnum_CopyMult(a,b) mpq_mul(a,a,b)        /* a = a*b */

#define QSnum_Add(a,b,c) mpq_add(a,b,c)           /* a = b+c */
#define QSnum_Sub(a,b,c) mpq_add(a,b,c)           /* a = a+b */
#define QSnum_Mult(a,b,c) mpq_mul(a,b,c)          /* a = b*c */

/* a = a*b (for b an unsigned int) */
#define QSnum_CopyMultUi(a,b) do {                                       \
        mpz_mul_ui(mpq_numref(a),mpq_numref(a),(unsigned long int)b);  \
        mpq_canonicalize(a);} while(0)

/* a = a - b*c */
#define QSnum_CopySubProd(a, b, c) do {                                  \
        mpq_mul (__t__, b, c);  mpq_sub (a, a, __t__);  } while(0)

/* a = max(a,|b|) */
#define QSnum_CopyMaxAbs(a, b) do {                                      \
        mpq_abs (__t__, b);                                            \
        if (mpq_cmp (a, __t__) < 0) mpq_set (a, __t__); } while(0)

/* a = max(a,abs(b)), execute extra code if update is needed */
#define QSnum_CopyMaxAbsAndDo(a,b,c)                                   \
    if(mpq_sgn(b) > 0) {                                               \
        if(QSnum_Less(a,b)){ QSnum_Copy(a,b); c; }                     \
    } else {                                                           \
        QSnum_Sign(a);                                                 \
        if(QSnum_Less(b,a)){ QSnum_Copy(a,b); c; }                     \
        QSnum_Sign(a); }

/*
 * basic comparison operations for QSnum_type
 */
#define QSnum_Equal(a,b) (mpq_equal(a,b))         /* a == b ? */
#define QSnum_EqualTol(a,b,tol) (mpq_equal(a,b))  /* |a-b| <= tol ? */
#define QSnum_Leq(a,b) (mpq_cmp(a,b) <= 0)        /* a <= b ? */
#define QSnum_Less(a,b) (mpq_cmp(a,b) < 0)        /* a < b ? */
#define QSnum_LessDbl(a,b) (mpq_get_d(a) < b)     /* a < b (for double b) ? */
#define QSnum_NeqZero(a) (mpq_sgn(a))             /* nonzero if a != 0 */
#define QSnum_NeqZeroTol(a,tol) (mpq_sgn(a))      /* |a >= tol? */

/*
 * basic assignments for QSnum_type
 */
#define QSnum_SetOne(a) mpq_set_ui(a,(unsigned long)1,(unsigned long)1)
                                                  /* a = 1 */
#define QSnum_SetZero(a) mpq_set_ui(a,(unsigned long)0,(unsigned long)1)
                                                  /* a = 0 */
#define QSnum_SetEpsilon(a) mpq_set_ui(a,(unsigned long)0,(unsigned long)1)
                                                  /* a = 2^-50 */
#define QSnum_SetDelta(a) mpq_set_ui(a,(unsigned long)0,(unsigned long)1)
                                                  /* a = 2^-7 */
/* swap function for QSnum_type */
#define QSnum_SWAP(a,b,t) ((QSnum_Copy(t,a)),(QSnum_Copy(a,b)),(QSnum_Copy(b,t)))

/* print a to FILE f */
#define QSnum_Print(f,a) mpq_out_str (f, 10, a)

/* read a from FILE f */
#define QSnum_Read(f,a) { mpq_inp_str (a,f,10);  mpq_canonicalize(a); }

#endif
#endif /* __QS_NUM_H */
