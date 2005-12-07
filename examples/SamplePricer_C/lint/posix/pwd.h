
/*  pwd.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__PWD_H)
#define __PWD_H

struct passwd {
	char *   pw_name;
	uid_t    pw_uid;
	gid_t    pw_gid;
	char *   pw_dir;
	char *   pw_shell;
};

struct passwd *getpwuid(uid_t);
struct passwd *getpwnam(const char *);

#endif
