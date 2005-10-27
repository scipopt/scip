
/*  unistd.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__UNISTD_H)
#define __UNISTD_H

#define R_OK   1
#define W_OK   2
#define X_OK   3
#define F_OK   4

#ifndef NULL
#define NULL ((void *) 0)
#endif

#ifndef __SIZE_T
#define __SIZE_T
typedef unsigned size_t;
#endif

#ifndef __SSIZE_T
#define __SSIZE_T
typedef int ssize_t;
#endif

#ifndef SEEK_CUR
#define SEEK_SET                0
#define SEEK_CUR                1
#define SEEK_END                2
#endif
#define STDERR_FILENO           2
#define STDIN_FILENO            0
#define STDOUT_FILENO           1
#define _PC_CHOWN_RESTRICTED    1
#define _PC_LINK_MAX            100
#define _PC_MAX_CANON           200
#define _PC_MAX_INPUT           200
#define _PC_NAME_MAX            200
#define _PC_NO_TRUNC            2
#define _PC_PATH_MAX            200
#define _PC_PIPE_BUF            200
#define _PC_VDISABLE            3
#define _POSIX_CHOWN_RESTRICTED 1
#define _POSIX_JOB_CONTROL      1
#define _POSIX_NO_TRUNC         3
#define _POSIX_SAVED_IDS        2
#define _POSIX_VDISABLE         3
#define _POSIX_VERSION          199009L
#define _SC_ARG_MAX             1
#define _SC_CHILD_MAX           2
#define _SC_CLK_TCK             3
#define _SC_JOB_CONTROL         4
#define _SC_NGROUPS_MAX         5
#define _SC_OPEN_MAX            6
#define _SC_SAVED_IDS           7
#define _SC_STREAM_MAX          8
#define _SC_TZNAME_MAX          9
#define _SC_VERSION             10

void _exit(int );
int access(const char *, int );
unsigned int alarm(unsigned int );
int chdir(const char *);
int chown(const char *, uid_t, gid_t );
int close(int );
char *ctermid(char *);
int dup2(int, int );
int dup(int );
int execl(const char *, const char *, ...);
int execle(const char *, const char *, ...);
int execlp(const char *, const char *, ...);
int execv(const char *, char *const []);
int execve(const char *, char *const [], char *const [] );
int execvp(const char *, char *const []);
pid_t fork(void);
long fpathconf(int, int );
char *getcwd(char *, size_t );
gid_t getegid(void);
uid_t geteuid(void);
gid_t getgid(void);
int getgroups(int, gid_t []);
char *getlogin(void);
pid_t getpgrp(void);
pid_t getpid(void);
pid_t getppid(void);
uid_t getuid(void);
int isatty(int );
int link(const char *, const char *);
off_t lseek(int, off_t, int );
long pathconf(const char *, int );
int pause(void);
int pipe(int  [2]);
ssize_t read(int, void *, size_t );
int rmdir(const char *);
int setgid(gid_t );
int setpgid(pid_t, pid_t );
pid_t setsid(void);
int setuid(uid_t );
unsigned int sleep(unsigned int );
long sysconf(int );
pid_t tcgetpgrp(int );
int tcsetpgrp(int, pid_t );
char *ttyname(int );
int unlink(const char *);
ssize_t write(int , const void *, size_t );

#endif
