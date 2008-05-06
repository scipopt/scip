
/*  strings.h
    added by Tobias Achterberg
 */

#if defined(_lint) && !defined(__STRINGS_H)
#define __STRINGS_H

int bcmp (const void *s1, const void *s2, size_t n);
void bcopy (const void *src, void *dest, size_t n);
void bzero (void *s, size_t n);
int ffs (int i);
char *index (const char *s, int c);
char *rindex (const char *s, int c);

#endif
