
/*  sys/times.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__TIMES_H)
#define __TIMES_H

#ifndef __CLOCK_T
typedef long clock_t;
#define __CLOCK_T
#endif

struct tms {
	clock_t  tms_utime;
	clock_t  tms_stime;
	clock_t  tms_cutime;
	clock_t  tms_cstime;
};

clock_t times(struct tms *);


#endif
