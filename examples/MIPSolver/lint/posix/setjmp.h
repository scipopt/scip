
/*  setjmp.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__SETJMP_H)
#define __SETJMP_H

typedef int jmp_buf[9];

void longjmp(jmp_buf , int );
int setjmp(jmp_buf );

typedef int sigjmp_buf[9];

void siglongjmp(sigjmp_buf, int );
int sigsetjmp(sigjmp_buf, int );

#endif
