
/*  string.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint) && !defined(__STRING_H)
#define __STRING_H

#ifndef NULL
#define NULL ((void *) 0)
#endif

#ifndef __SIZE_T
#define __SIZE_T
typedef unsigned size_t;
#endif

void * memchr(const void *, int, size_t);
int memcmp(const void *, const void *, size_t);
void * memcpy(void *, const void *, size_t);
void * memmove(void *, const void *, size_t);
void * memset(void *, int, size_t);
char * strcat(char *, const char *);
char * strchr(const char *, int);
int strcmp(const char *, const char *);
int strcoll(const char *, const char *);
char * strcpy(char *, const char *);
size_t strcspn(const char *, const char *);
char * strerror(int);
size_t strlen(const char *);
char * strncat(char *, const char *, size_t);
int strncmp(const char *, const char *, size_t );
char * strncpy(char *, const char *, size_t);
char * strpbrk(const char *, const char *);
char * strrchr(const char *, int);
size_t strspn(const char *, const char *);
char * strstr(const char *, const char *);
char * strtok(char *, const char *);
size_t strxfrm(char *, const char *, size_t);

#endif
