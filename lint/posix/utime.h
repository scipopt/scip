
/*  utime.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__UTIME_H)
#define __UTIME_H

struct utimbuf {
	time_t  actime;
	time_t  modtime;
};

int utime(const char *, const struct utimbuf * );

#endif
