
/*  grp.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__GRP_H)
#define __GRP_H

struct  group {
	char *    gr_name;
	gid_t     gr_gid;
	char **   gr_mem;
};

struct group *getgrgid(gid_t );
struct group *getgrnam(const char * );

#endif
