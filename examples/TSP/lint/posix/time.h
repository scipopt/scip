
/*  time.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__TIME_H)
#define __TIME_H

#ifndef NULL
#define NULL ((void *) 0)
#endif

#ifndef __SIZE_T
#define __SIZE_T
typedef unsigned  size_t;
#endif

#ifndef __TIME_T
typedef long      time_t;
#define __TIME_T
#endif

#ifndef __CLOCK_T
typedef long      clock_t;
#define __CLOCK_T
#endif

struct tm {
    int tm_sec;         /* seconds after the minute - [0,59] */
    int tm_min;         /* minutes after the hour - [0,59] */
    int tm_hour;        /* hours since midnight - [0,23] */
    int tm_mday;        /* day of the month - [1,31] */
    int tm_mon;         /* months since January - [0,11] */
    int tm_year;        /* years since 1900 */
    int tm_wday;        /* days since Sunday - [0,6] */
    int tm_yday;        /* days since January 1 - [0,365] */
    int tm_isdst;       /* daylight savings time flag */
};

#define CLOCKS_PER_SEC 1000
#define CLK_TCK  60

char *asctime(const struct tm *);
char *ctime(const time_t *);
clock_t  clock(void);
double  difftime(time_t, time_t);
struct tm *gmtime(const time_t *);
struct tm *localtime(const time_t *);
time_t  mktime(struct tm *);
time_t  time(time_t *);
size_t strftime(char *, size_t, const char *, const struct tm * );
void tzset(void);

extern char *tzname[2];

#endif
