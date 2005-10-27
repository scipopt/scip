
/*  errno.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__ERRNO_H)
#define __ERRNO_H

extern int volatile errno;

#define EDOM            1
#define ERANGE          2

#define E2BIG           3
#define EACCES          4
#define EAGAIN          5
#define EBADF           6
#define EBUSY           7
#define ECHILD          8
#define EDEADLK         9
#define EEXIST          10
#define EFAULT          11
#define EFBIG           12
#define EINTR           13
#define EINVAL          14
#define EIO             15
#define EISDIR          16
#define EMFILE          17
#define EMLINK          18
#define ENAMETOOLONG    19
#define ENFILE          20
#define ENODEV          21
#define ENOENT          22
#define ENOEXEC         23
#define ENOLCK          24
#define ENOMEM          25
#define ENOSPC          26
#define ENOSYS          27
#define ENOTDIR         28
#define ENOTEMPTY       29
#define ENOTTY          30
#define ENXIO           31
#define EPERM           32
#define EPIPE           33
#define EROFS           34
#define ESPIPE          35
#define ESRCH           36
#define EXDEV           37

#endif
