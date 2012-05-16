
/*  stdio.h --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */
#ifndef __STDIO_H
#define __STDIO_H
#if defined(_lint)


#ifndef __SIZE_T
#define __SIZE_T
typedef unsigned size_t;
#endif

#ifndef NULL
#define NULL ((void *) 0)
#endif

typedef long    fpos_t;

typedef struct  {
	long _inside_FILE;
}       FILE;

#define _IOFBF  0
#define _IOLBF  1
#define _IONBF  2

#define EOF (-1)
#define FOPEN_MAX 20
#define FILENAME_MAX 100
#define BUFSIZ  512
#define L_tmpnam    12
#ifndef SEEK_CUR
#define SEEK_CUR    1
#define SEEK_END    2
#define SEEK_SET    0
#endif
#define TMP_MAX     30
extern  FILE    _streams[];

#define stdin   (&_streams[0])
#define stdout  (&_streams[1])
#define stderr  (&_streams[2])


void clearerr(FILE *);
int  fclose(FILE *);
int  feof(FILE *);
int  ferror(FILE *);
int  fflush(FILE *);
int  fgetc(FILE *);
int  fgetpos(FILE *, fpos_t *);
char *fgets(char *, int , FILE *);
FILE *fopen(const char *, const char *);
int  fprintf(FILE *, const char *, ...);
int  fputc(int , FILE *);
int  fputs(const char *, FILE *);
size_t fread(void *, size_t , size_t , FILE *);
FILE *freopen(const char *, const char *, FILE *);
int  fscanf(FILE *, const char *, ...);
int  fseek(FILE *, long , int );
int  fsetpos(FILE *, const fpos_t *);
long ftell(FILE *);
size_t fwrite(const void *, size_t , size_t , FILE *);
char *gets(char *);
void perror(const char *);
int  printf(const char *, ...);
int  puts(const char *);
int  remove(const char *);
int  rename(const char *,const char *);
void rewind(FILE *);
int  scanf(const char *, ...);
void setbuf(FILE *, char *);
int  setvbuf(FILE *, char *, int, size_t );
int  sprintf(char *, const char *, ...);
int  snprintf(char *, size_t, const char *, ...);
int  sscanf(const char *, const char *, ...);
FILE *tmpfile(void);
char *tmpnam(char *);
int  ungetc(int , FILE *);

int getc(FILE *);
int putc(int, FILE *);
int getchar(void);
int putchar(int);

#ifndef __VA_LIST
#define __VA_LIST
typedef char *va_list;
#endif

int  vfprintf(FILE *, const char *, va_list );
int  vprintf( const char *, va_list );
int  vsprintf(char *,  const char *, va_list );

#define L_ctermid   100

#ifndef __STREAM_MAX
#define __STREAM_MAX
#define STREAM_MAX   20
#endif

FILE *fdopen(int fildes, const char *type);
int fileno(FILE *stream);

#endif
#endif
