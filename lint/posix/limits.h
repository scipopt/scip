
/*  limits.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__LIMITS_H)
#define __LIMITS_H

#define CHAR_BIT	8
#define CHAR_MIN        (-128)
#define CHAR_MAX        127
#define SCHAR_MAX	127
#define SCHAR_MIN	(-128)
#define UCHAR_MAX	255U
#define SHRT_MAX	32767
#define SHRT_MIN	(-32767-1)
#define USHRT_MAX	65535U
#define INT_MAX         (int) 2147483647
#define INT_MIN         (int) (-2147483647-1)
#define UINT_MAX        4294967295U
#define LONG_MAX	2147483647
#define LONG_MIN	(-2147483647-1)
#define ULONG_MAX       0xFFFFFFFF
#define MB_LEN_MAX	1

#define NGROUPS_MAX        100
#define SSIZE_MAX        32767
#define _POSIX_ARG_MAX    4096
#define _POSIX_CHILD_MAX     6
#define _POSIX_LINK_MAX      8
#define _POSIX_MAX_CANON   255
#define _POSIX_MAX_INPUT   255
#define _POSIX_NAME_MAX     14
#define _POSIX_NGROUPS_MAX   0
#define _POSIX_OPEN_MAX     16
#define _POSIX_PATH_MAX    255
#define _POSIX_PIPE_BUF    512
#define _POSIX_SSIZE_MAX 32767
#define _POSIX_STREAM_MAX   20
#define _POSIX_TZNAME_MAX    3

#ifndef _USE_SYSCONF
#define ARG_MAX         5120
#define CHILD_MAX       25
#define OPEN_MAX        60
#ifndef __STREAM_MAX
#define __STREAM_MAX
#define STREAM_MAX   20
#endif
#define TZNAME_MAX      100
#endif

#ifndef _USE_PATHCONF
#define LINK_MAX        1000
#define MAX_CANON       255
#define MAX_INPUT       255
#define NAME_MAX        100
#define PATH_MAX        256
#define PIPE_BUF        5120
#endif

#endif
