
/*  sys/utsname.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__UTSNAME_H)
#define __UTSNAME_H

struct utsname {
	char    sysname[9];
	char    nodename[9];
	char    release[9];
	char    version[9];
	char    machine[9];
};

int uname(struct utsname *);

#endif
