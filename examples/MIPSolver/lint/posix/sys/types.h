
/*  sys/types.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__TYPES_H)
#define __TYPES_H

typedef short           ino_t;
typedef short           dev_t;
typedef long            off_t;

#ifndef __SIZE_T
#define __SIZE_T
typedef unsigned    size_t;
#endif

typedef  long           gid_t;
typedef  long           mode_t;
typedef  long           nlink_t;
typedef  long           pid_t;
typedef  long           uid_t;

#ifndef __SSIZE_T
#define __SSIZE_T
typedef int ssize_t;
#endif


#endif
