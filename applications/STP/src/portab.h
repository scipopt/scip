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

/**@file   portab.h
 * @brief  Portable definitions
 * @author Thorsten Koch
 * @author Daniel Rehfeldt
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef PORTAB_H
#define PORTAB_H

#include <limits.h>
#include <math.h>
#include <float.h>

#define EPSILON 1e-06

typedef unsigned char STP_Bool;


#ifndef TRUE
#define TRUE            1
#endif
#ifndef FALSE
#define FALSE           0
#endif

#ifndef SUCCESS
#define SUCCESS         0
#endif
#ifndef FAILURE
#define FAILURE       (-1)
#endif

#ifndef Max
#define Max(a, b)    (((a) >= (b)) ? (a) : (b))
#endif
#ifndef Min
#define Min(a, b)    (((a) <= (b)) ? (a) : (b))
#endif

#define EPS_ZERO  (DBL_EPSILON * 1e4)

#ifndef Fsgn
#define Fsgn(x)   ((((x) > -EPS_ZERO) && ((x) < EPS_ZERO)) ? 0 : (((x) < 0.0) ? -1 : 1))
#endif /* fsgn */


#define RELDIFF(a, b) (((a)-(b))/ MAX(MAX(fabs(a), fabs(b)), 1.0))
#define EQ_FEAS(a, b)  (fabs(RELDIFF(a, b)) <= EPS_ZERO)
#define LT_FEAS(a, b)  (RELDIFF(a, b) < -EPS_ZERO)
#define LE_FEAS(a, b)  (RELDIFF(a, b) < EPS_ZERO)
#define GT_FEAS(a, b)  (RELDIFF(a, b) > EPS_ZERO)
#define GE_FEAS(a, b)  (RELDIFF(a, b) > -EPS_ZERO)
#define EQ_FEAS_EPS(a, b, eps)  (fabs(RELDIFF(a, b)) <= eps)
#define LT_FEAS_EPS(a, b, eps)  (RELDIFF(a, b) < -eps)
#define LE_FEAS_EPS(a, b, eps)  (RELDIFF(a, b) < eps)
#define GT_FEAS_EPS(a, b, eps)  (RELDIFF(a, b) > eps)
#define GE_FEAS_EPS(a, b, eps)  (RELDIFF(a, b) > -eps)

#define EQ(a, b)   (fabs((a) - (b)) <= EPS_ZERO)
#define NE(a, b)   (fabs((a) - (b)) >  EPS_ZERO)
#define LT(a, b)   (((a) - (b))     < -EPS_ZERO)
#define LE(a, b)   (((a) - (b))     <  EPS_ZERO)
#define GT(a, b)   (((a) - (b))     >  EPS_ZERO)
#define GE(a, b)   (((a) - (b))     > -EPS_ZERO)

#ifdef __GNUC__
#define does_not_return   __attribute__((noreturn))
#else
#define does_not_return   /**/
#define inline            /**/
#endif

#if defined(MSDOS) || defined(WIN32)
#define  DIRSEP "\\"
#else
#define  DIRSEP "/"
#endif

#endif   /* PORTAB_H */
