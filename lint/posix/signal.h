
/*  signal.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__SIGNAL_H)
#define __SIGNAL_H

typedef void (*_sigfunc_t)(int);   /* the type returned
	    by and passed to (as second argument) the signal function */

#define SIG_DFL (_sigfunc_t)0   /* default action */
#define SIG_IGN (_sigfunc_t)1   /* ignore the signal */
#define SIG_ERR (_sigfunc_t)(-1)  /* error return */

extern _sigfunc_t signal(int, _sigfunc_t);
extern int raise(int);

typedef int sig_atomic_t;
typedef int sigset_t;

#define SIGHUP         1
#define SIGINT         2
#define SIGQUIT        3
#define SIGILL         4
#define SIGABRT        6
#define SIGFPE         8
#define SIGKILL        9
#define SIGSEGV       11
#define SIGPIPE       13
#define SIGALRM       14
#define SIGTERM       15
#define SIGUSR1       16
#define SIGUSR2       17
#define SA_NOCLDSTOP  18
#define SIGCHLD       19
#define SIGCONT       20
#define SIGSTOP       21
#define SIGTSTP       22
#define SIGTTIN       23
#define SIGTTOU       24
#define SIG_BLOCK     25
#define SIG_SETMASK   26
#define SIG_UNBLOCK   27

struct sigaction {
    int sa_flags;
    void(* sa_handler)();
    sigset_t sa_mask;
};

int kill(pid_t, int);
int sigaction(int, const struct sigaction *, struct sigaction *);
int sigaddset(sigset_t *, int);
int sigdelset(sigset_t *, int);
int sigemptyset(sigset_t *);
int sigfillset(sigset_t *);
int sigismember(const sigset_t *, int);
int sigpending(sigset_t *);
int sigprocmask(int how, const sigset_t *, sigset_t *);
int sigsuspend(const sigset_t *);

#endif
