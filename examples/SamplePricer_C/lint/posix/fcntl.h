
/*  fcntl.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__FCNTL_H)
#define __FCNTL_H

#define O_RDONLY 0
#define O_WRONLY 1
#define O_RDWR   2
#define O_APPEND 010
#define O_ACCMODE 020
#define O_NOCTTY  030
#define O_NONBLOCK 040


#define O_CREAT 00400
#define O_TRUNC 01000
#define O_EXCL  02000


#define F_DUPFD     0
#define F_GETFD     1
#define F_SETFD     2
#define F_GETFL     3
#define F_SETFL     4
#define F_GETLK     5
#define F_SETLK     6
#define F_SETLKW    7
#define FD_CLOEXEC  8


struct flock {
	short   l_type;
	short   l_whence;
	off_t   l_start;
	off_t   l_len;
	pid_t   l_pid;
};


#define F_RDLCK 01
#define F_WRLCK 02
#define F_UNLCK 03

int creat(const char *, mode_t );
int fcntl(int, int, ...);
int open(const char *, int, ...);

#endif
