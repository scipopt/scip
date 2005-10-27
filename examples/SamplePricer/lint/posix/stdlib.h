
/*  stdlib.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__STDLIB_H)
#define __STDLIB_H

#ifndef __SIZE_T
#define __SIZE_T
typedef unsigned size_t;
#endif

#ifndef __DIV_T
#define __DIV_T
typedef struct {
 int quot;
 int rem;
} div_t;
#endif

#ifndef __LDIV_T
#define __LDIV_T
typedef struct {
 long quot;
 long rem;
} ldiv_t;
#endif

#define MB_CUR_MAX  _Mbcurmax
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1

#define RAND_MAX 0x7FFF

typedef void (* __atexit_t)(void);

#ifndef __WCHAR_T
#define __WCHAR_T
typedef char wchar_t;
#endif

void abort (void);
int abs (int);
int atexit (__atexit_t);
double atof (const char *);
int atoi (const char *);
long atol (const char *);
void * bsearch(const void *, const void *, size_t, size_t,
	 int (*) (const void *, const void *));
void * calloc (size_t, size_t);
div_t div (int , int);
void exit (int);
void free (void *);
char * getenv (const char *);
long labs (long);
ldiv_t ldiv (long, long);
void * malloc (size_t);
void qsort (void *, size_t, size_t, int (*) (const void *, const void *));
int rand (void);
void * realloc(void *, size_t);
void srand (unsigned);
double strtod (const char *, char **);
long  strtol (const char *, char **, int);
unsigned long strtoul (const char *, char **, int);
int system (const char *);

#ifndef NULL
#define NULL ((void *) 0)
#endif

extern char _Mbcurmax;

extern int mblen(const char *, size_t);
extern int mbtowc(wchar_t *, const char *, size_t);
extern int wctomb(char *, wchar_t);
extern size_t mbstowcs(wchar_t *, const char *, size_t);
extern size_t wcstombs(char *, const wchar_t *, size_t);

#endif
