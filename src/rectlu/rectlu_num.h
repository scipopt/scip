/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                        This file is part of the program                   */
/*            RECTLU --- Algorithm for Exact Rectangular LU Factorization    */
/*                                                                           */
/*    Copyright (C) 2009-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*   RECTLU is distributed under the terms of the ZIB Academic License.      */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with RECTLU; see the file COPYING.                                 */
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



#ifndef __RECTLU_NUM_H
#define __RECTLU_NUM_H

#define QSnum_type mpq_t

#define QSnum_AllocArray(size) ({                                \
        size_t __i__ = (size);                                   \
        QSnum_type *__res = (QSnum_type *) malloc(size*sizeof(QSnum_type));  \
        if (__res) while(__i__--) mpq_init(__res[__i__]);        \
        __res;})

#define QSnum_FreeArray(ea,size) ({\
        size_t __sz = size;\
        QSnum_type* __ptr__ = (ea);\
        if (ea) while(__sz--) mpq_clear(__ptr__[__sz]);\
        if (ea) { free(ea);  ea = NULL;}})

#define QSnum_Init(a) mpq_init(a)
#define QSnum_Clear(a) mpq_clear(a)

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
#define QSnum_CopyMultUi(a,b) ({                                       \
        mpz_mul_ui(mpq_numref(a),mpq_numref(a),(unsigned long int)b);  \
        mpq_canonicalize(a);(unsigned long)0;})

/* a = a - b*c */
#define QSnum_CopySubProd(a, b, c) ({                                  \
        mpq_mul (__t__, b, c);  mpq_sub (a, a, __t__);  })

/* a = max(a,|b|) */
#define QSnum_CopyMaxAbs(a, b) ({                                      \
        mpq_abs (__t__, b);                                            \
        if (mpq_cmp (a, __t__) < 0) mpq_set (a, __t__); })

/* a = max(a,abs(b)), execute extra code if update is needed */
#define QSnum_CopyMaxAbsAndDo(a,b,c)                                   \
    if(mpq_sgn(b) > 0) {                                               \
        if(QSnum_Less(a,b)){ QSnum_Copy(a,b); c; }                     \
    } else {                                                           \
        QSnum_Sign(a);                                                 \
        if(QSnum_Less(b,a)){ QSnum_Copy(a,b); c; }                     \
        QSnum_Sign(a); }

#define QSnum_Equal(a,b) (mpq_equal(a,b))         /* a == b ? */
#define QSnum_EqualTol(a,b,tol) (mpq_equal(a,b))  /* |a-b| <= tol ? */
#define QSnum_Leq(a,b) (mpq_cmp(a,b) <= 0)        /* a <= b ? */
#define QSnum_Less(a,b) (mpq_cmp(a,b) < 0)        /* a < b ? */
#define QSnum_LessDbl(a,b) (mpq_get_d(a) < b)     /* a < b (for double b) ? */
#define QSnum_NeqZero(a) (mpq_sgn(a))             /* nonzero if a != 0 */
#define QSnum_NeqZeroTol(a,tol) (mpq_sgn(a))      /* |a >= tol? */

#define QSnum_SetOne(a) mpq_set_ui(a,(unsigned long)1,(unsigned long)1)
                                                  /* a = 1 */
#define QSnum_SetZero(a) mpq_set_ui(a,(unsigned long)0,(unsigned long)1)
                                                  /* a = 0 */
#define QSnum_SetEpsilon(a) mpq_set_ui(a,(unsigned long)0,(unsigned long)1)
                                                  /* a = 2^-50 */
#define QSnum_SetDelta(a) mpq_set_ui(a,(unsigned long)0,(unsigned long)1)
                                                  /* a = 2^-7 */

#define QSnum_SWAP(a,b,t) ((QSnum_Copy(t,a)),(QSnum_Copy(a,b)),(QSnum_Copy(b,t)))

#define QSnum_Print(f,a) mpq_out_str (f, 10, a)  /* print a to FILE f */
#define QSnum_Read(f,a) { mpq_inp_str (a,f,10);  mpq_canonicalize(a); }
                                                 /* read a from FILE f */

#endif /* __QS_NUM_H */
