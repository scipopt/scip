
/*  stddef.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__STDDEF_H)
#define __STDDEF_H

#ifndef __PTRDIFF_T
#define __PTRDIFF_T
typedef int	ptrdiff_t;
#endif

#ifndef __SIZE_T
#define __SIZE_T
typedef unsigned size_t;
#endif

#ifndef __WCHAR_T
#define __WCHAR_T
typedef char wchar_t;
#endif

#ifndef NULL
#define NULL ((void *) 0)
#endif

#define offsetof(T,memb) ((size_t)&(((T *)0)->memb)-(size_t)((T *)0))

#endif
