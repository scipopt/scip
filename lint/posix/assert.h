
/*  assert.h  --  POSIX header
    Copyright 1992, Gimpel Software, All Rights Reserved
 */

#if defined(_lint)

#ifdef NDEBUG
#define assert(p) ((void)0)
#else
void __assertfail( const char *, ... );
#define assert(p)  ((p) ? (void) 0 : __assertfail(   \
		"Assertion failed: %s, file %s, line %d\n", #p,\
		__FILE__, __LINE__) )
#endif
#endif
