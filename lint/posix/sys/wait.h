
/*  sys/wait.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__WAIT_H)
#define __WAIT_H

#define WEXITSTATUS  1
#define WIFEXITED    2
#define WIFSIGNALED  3
#define WIFSTOPPED   4
#define WNOHANG      5
#define WSTOPSIG     6
#define WTERMSIG     7
#define WUNTRACED    8

pid_t wait(int * );
pid_t waitpid(pid_t, int *, int );

#endif
