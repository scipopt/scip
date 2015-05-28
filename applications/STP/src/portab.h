/* $Id: portab.h,v 1.3 1999/07/27 14:24:10 bzfkocht Exp $ */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Header                                                        */
/*   File....: portab.h                                                      */
/*   Name....: Portable Defintions                                           */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef PORTAB_H
#define PORTAB_H

#include <limits.h>
#include <math.h>
#include <float.h>

#define EPSILON 1e-06

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

/* Hi Fortran !
 */
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
