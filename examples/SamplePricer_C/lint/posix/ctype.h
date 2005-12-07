
/*  ctype.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__CTYPE_H)
#define __CTYPE_H

extern  int     isalnum(int);
extern  int     isalpha(int);
extern  int     iscntrl(int);
extern  int     isdigit(int);
extern  int     isgraph(int);
extern  int     islower(int);
extern  int     isprint(int);
extern  int     ispunct(int);
extern  int     isspace(int);
extern  int     isupper(int);
extern  int     isxdigit(int);
extern  int     tolower(int);
extern  int     toupper(int);

#endif
