/* This file has all the included material used by cliquer-1.21
   and prototypes for the interface procedures in nautycliquer.c */

#ifndef NAUTYCLIQUER_H
#define NAUTYCLIQUER_H

#include "nauty.h"
#include "gtools.h"
#include <limits.h>

/**********************************************************************
#include "cliquerconf.h"
*/

/*
 * setelement is the basic memory type used in sets.  It is often fastest
 * to be as large as can fit into the CPU registers.
 *
 * ELEMENTSIZE is the size of one setelement, measured in bits.  It must
 * be either 16, 32 or 64  (otherwise additional changes must be made to
 * the source).
 *
 * The default is to use "unsigned long int" and attempt to guess the
 * size using <limits.h>, which should work pretty well.  Check functioning
 * with "make test".
 */

/* typedef unsigned long int setelement; */
/* #define ELEMENTSIZE 64 */

/*
 * INLINE is a command prepended to function declarations to instruct the
 * compiler to inline the function.  If inlining is not desired, define blank.
 *
 * The default is to use "inline", which is recognized by most compilers.
 */

/* #define INLINE */
/* #define INLINE __inline__ */


/*
 * Set handling functions are defined as static functions in set.h for
 * performance reasons.  This may cause unnecessary warnings from the
 * compiler.  Some compilers (such as GCC) have the possibility to turn
 * off the warnings on a per-function basis using a flag prepended to
 * the function declaration.
 *
 * The default is to use the correct attribute when compiling with GCC,
 * or no flag otherwise.
 */

/* #define UNUSED_FUNCTION __attribute__((unused)) */
/* #define UNUSED_FUNCTION */

/*
 * Uncommenting the following will disable all assertions  (checks that
 * function arguments and other variables are correct).  This is highly
 * discouraged, as it allows bugs to go unnoticed easier.  The assertions
 * are set so that they do not slow down programs notably.
 */

/* #define ASSERT(x) */

/**********************************************************************
#include "misc.h"
*/

/*
 * We #define boolean instead of using a typedef because nauty.h uses it
 * also.  AFAIK, there is no way to check for an existing typedef, and
 * re-typedefing is illegal (even when using exactly the same datatype!).
#ifndef boolean
#define boolean int
#endif

BDM: In nauty's version we will use nauty's boolean (which is int anyway).
 */

/*
 * Default value for UNUSED_FUNCTION:  use "__attribute__((unused))" for
 * GCC versions that support it, otherwise leave blank.
 */
#ifndef UNUSED_FUNCTION
# if     __GNUC__ > 2 || (__GNUC__ == 2 && __GNUC_MINOR__ > 4)
#  define UNUSED_FUNCTION __attribute__((unused))
# else
#  define UNUSED_FUNCTION
# endif
#endif  /* !UNUSED_FUNCTION */

/*
 * Default inlining directive:  "inline"
 */
#ifndef INLINE
#define INLINE inline
#endif

#ifndef ASSERT
#define ASSERT(expr) \
        if (!(expr)) { \
		fprintf(stderr,"cliquer file %s: line %d: assertion failed: " \
			"(%s)\n",__FILE__,__LINE__,#expr); \
		abort(); \
	}
#endif /* !ASSERT */


#ifndef FALSE
#define FALSE (0)
#endif
#ifndef TRUE
#define TRUE (!FALSE)
#endif


#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif
#ifndef ABS
#define ABS(v)  (((v)<0)?(-(v)):(v))
#endif

/**********************************************************************
#include "set.h"
*/

/*
 * This file contains the set handling routines.
 *
 * Copyright (C) 2002 Sampo Niskanen, Patric Östergård.
 * Licensed under the GNU GPL, read the file LICENSE for details.
 */

/*
 * Sets are arrays of setelement's (typically unsigned long int's) with
 * representative bits for each value they can contain.  The values
 * are numbered 0,...,n-1.
 */


/*** Variable types and constants. ***/


/*
 * If setelement hasn't been declared:
 *   - use "unsigned long int" as setelement
 *   - try to deduce size from ULONG_MAX
 */

#ifndef ELEMENTSIZE
typedef unsigned long int setelement;
# if (ULONG_MAX == 65535)
#  define ELEMENTSIZE 16
# elif (ULONG_MAX == 4294967295)
#  define ELEMENTSIZE 32
# else
#  define ELEMENTSIZE 64
# endif
#endif  /* !ELEMENTSIZE */

typedef setelement * set_t;


/*** Counting amount of 1 bits in a setelement ***/

/* Array for amount of 1 bits in a byte. */
static int set_bit_count[256] = {
	0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
	4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8 };

/* The following macros assume that all higher bits are 0.
 * They may in some cases be useful also on with other ELEMENTSIZE's,
 * so we define them all.  */
#define SET_ELEMENT_BIT_COUNT_8(a)  (set_bit_count[(a)])
#define SET_ELEMENT_BIT_COUNT_16(a) (set_bit_count[(a)>>8] + \
				     set_bit_count[(a)&0xFF])
#define SET_ELEMENT_BIT_COUNT_32(a) (set_bit_count[(a)>>24] + \
				     set_bit_count[((a)>>16)&0xFF] + \
				     set_bit_count[((a)>>8)&0xFF] + \
				     set_bit_count[(a)&0xFF])
#define SET_ELEMENT_BIT_COUNT_64(a) (set_bit_count[(a)>>56] + \
				     set_bit_count[((a)>>48)&0xFF] + \
				     set_bit_count[((a)>>40)&0xFF] + \
				     set_bit_count[((a)>>32)&0xFF] + \
				     set_bit_count[((a)>>24)&0xFF] + \
				     set_bit_count[((a)>>16)&0xFF] + \
				     set_bit_count[((a)>>8)&0xFF] + \
				     set_bit_count[(a)&0xFF])
#if (ELEMENTSIZE==64)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_64(a)
# define FULL_ELEMENT ((setelement)0xFFFFFFFFFFFFFFFF)
#elif (ELEMENTSIZE==32)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_32(a)
# define FULL_ELEMENT ((setelement)0xFFFFFFFF)
#elif (ELEMENTSIZE==16)
# define SET_ELEMENT_BIT_COUNT(a) SET_ELEMENT_BIT_COUNT_16(a)
# define FULL_ELEMENT ((setelement)0xFFFF)
#else
# error "SET_ELEMENT_BIT_COUNT(a) not defined for current ELEMENTSIZE"
#endif

/*** Macros and functions ***/

/*
 * Gives a value with bit x (counting from lsb up) set.
 *
 * Making this as a table might speed up things on some machines
 * (though on most modern machines it's faster to shift instead of
 * using memory).  Making it a macro makes it easy to change.
 */
#define SET_BIT_MASK(x) ((setelement)1<<(x))

/* Set element handling macros */

#define SET_ELEMENT_INTERSECT(a,b)  ((a)&(b))
#define SET_ELEMENT_UNION(a,b)      ((a)|(b))
#define SET_ELEMENT_DIFFERENCE(a,b) ((a)&(~(b)))
#define SET_ELEMENT_CONTAINS(e,v)   ((e)&SET_BIT_MASK(v))

/* Set handling macros */

#define SET_ADD_ELEMENT(s,a) \
                       ((s)[(a)/ELEMENTSIZE] |= SET_BIT_MASK((a)%ELEMENTSIZE))
#define SET_DEL_ELEMENT(s,a) \
                       ((s)[(a)/ELEMENTSIZE] &= ~SET_BIT_MASK((a)%ELEMENTSIZE))
#define SET_CONTAINS_FAST(s,a) (SET_ELEMENT_CONTAINS((s)[(a)/ELEMENTSIZE], \
						      (a)%ELEMENTSIZE))
#define SET_CONTAINS(s,a) (((a)<SET_MAX_SIZE(s))?SET_CONTAINS_FAST(s,a):FALSE)

/* Sets can hold values between 0,...,SET_MAX_SIZE(s)-1 */
#define SET_MAX_SIZE(s) ((s)[-1])
/* Sets consist of an array of SET_ARRAY_LENGTH(s) setelements */
#define SET_ARRAY_LENGTH(s) (((s)[-1]+ELEMENTSIZE-1)/ELEMENTSIZE)

/*
 * set_new()
 *
 * Create a new set that can hold values in the range 0,...,size-1.
 */
UNUSED_FUNCTION
static set_t set_new(int size) {
	int n;
	set_t s;

	ASSERT(size>0);

	n=(size/ELEMENTSIZE+1)+1;
	s=calloc(n,sizeof(setelement));
	s[0]=size;

	return &(s[1]);
}

/*
 * set_free()
 *
 * Free the memory associated with set s.
 */
UNUSED_FUNCTION INLINE
static void set_free(set_t s) {
	ASSERT(s!=NULL);
	free(&(s[-1]));
}

/*
 * set_resize()
 *
 * Resizes set s to given size.  If the size is less than SET_MAX_SIZE(s),
 * the last elements are dropped.
 *
 * Returns a pointer to the new set.
 */
UNUSED_FUNCTION INLINE
static set_t set_resize(set_t s, int size) {
	int n;

	ASSERT(size>0);

	n=(size/ELEMENTSIZE+1);
	s=((setelement *)realloc(s-1,(n+1)*sizeof(setelement)))+1;

	if (n>SET_ARRAY_LENGTH(s))
		memset(s+SET_ARRAY_LENGTH(s),0,
		       (n-SET_ARRAY_LENGTH(s))*sizeof(setelement));
	if (size < SET_MAX_SIZE(s))
		s[(size-1)/ELEMENTSIZE] &= (FULL_ELEMENT >>
					    (ELEMENTSIZE-size%ELEMENTSIZE));
	s[-1]=size;

	return s;
}

/*
 * set_size()
 *
 * Returns the number of elements in set s.
 */
UNUSED_FUNCTION INLINE
static int set_size(set_t s) {
	int count=0;
	setelement *c;

	for (c=s; c < s+SET_ARRAY_LENGTH(s); c++)
		count+=SET_ELEMENT_BIT_COUNT(*c);
	return count;
}

/*
 * set_duplicate()
 *
 * Returns a newly allocated duplicate of set s.
 */
UNUSED_FUNCTION INLINE
static set_t set_duplicate(set_t s) {
	set_t new;

	new=set_new(SET_MAX_SIZE(s));
	memcpy(new,s,SET_ARRAY_LENGTH(s)*sizeof(setelement));
	return new;
}

/*
 * set_copy()
 *
 * Copies set src to dest.  If dest is NULL, is equal to set_duplicate.
 * If dest smaller than src, it is freed and a new set of the same size as
 * src is returned.
 */
UNUSED_FUNCTION INLINE
static set_t set_copy(set_t dest,set_t src) {
	if (dest==NULL)
		return set_duplicate(src);
	if (SET_MAX_SIZE(dest)<SET_MAX_SIZE(src)) {
		set_free(dest);
		return set_duplicate(src);
	}
	memcpy(dest,src,SET_ARRAY_LENGTH(src)*sizeof(setelement));
	memset(dest+SET_ARRAY_LENGTH(src),0,((SET_ARRAY_LENGTH(dest) -
					      SET_ARRAY_LENGTH(src)) *
					     sizeof(setelement)));
	return dest;
}

/*
 * set_empty()
 *
 * Removes all elements from the set s.
 */
UNUSED_FUNCTION INLINE
static void set_empty(set_t s) {
	memset(s,0,SET_ARRAY_LENGTH(s)*sizeof(setelement));
	return;
}

/*
 * set_intersection()
 *
 * Store the intersection of sets a and b into res.  If res is NULL,
 * a new set is created and the result is written to it.  If res is
 * smaller than the larger one of a and b, it is freed and a new set
 * is created and the result is returned.
 *
 * Returns either res or a new set that has been allocated in its stead.
 *
 * Note:  res may not be a or b.
 */
UNUSED_FUNCTION INLINE
static set_t set_intersection(set_t res,set_t a,set_t b) {
	int i,max;

	if (res==NULL) {
		res = set_new(MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b)));
	} else if (SET_MAX_SIZE(res) < MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b))) {
		set_free(res);
		res = set_new(MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b)));
	} else {
		set_empty(res);
	}

	max=MIN(SET_ARRAY_LENGTH(a),SET_ARRAY_LENGTH(b));
	for (i=0; i<max; i++) {
		res[i]=SET_ELEMENT_INTERSECT(a[i],b[i]);
	}

	return res;
}

/*
 * set_union()
 *
 * Store the union of sets a and b into res.  If res is NULL, a new set
 * is created and the result is written to it.  If res is smaller than
 * the larger one of a and b, it is freed and a new set is created and
 * the result is returned.
 *
 * Returns either res or a new set that has been allocated in its stead.
 *
 * Note:  res may not be a or b.
 */
UNUSED_FUNCTION INLINE
static set_t set_union(set_t res,set_t a,set_t b) {
	int i,max;

	if (res==NULL) {
		res = set_new(MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b)));
	} else if (SET_MAX_SIZE(res) < MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b))) {
		set_free(res);
		res = set_new(MAX(SET_MAX_SIZE(a),SET_MAX_SIZE(b)));
	} else {
		set_empty(res);
	}

	max=MAX(SET_ARRAY_LENGTH(a),SET_ARRAY_LENGTH(b));
	for (i=0; i<max; i++) {
		res[i]=SET_ELEMENT_UNION(a[i],b[i]);
	}

	return res;
}


/*
 * set_return_next()
 *
 * Returns the smallest value in set s which is greater than n, or -1 if
 * such a value does not exist.
 *
 * Can be used to iterate through all values of s:
 *
 * int i=-1;
 * while ((i=set_return_next(s,i))>=0) {
 *         // i is in set s
 * }
 */
UNUSED_FUNCTION INLINE
static int set_return_next(set_t s, int n) {
	if (n<0)
		n=0;
	else
		n++;
	if (n >= SET_MAX_SIZE(s))
		return -1;

	while (n%ELEMENTSIZE) {
		if (SET_CONTAINS(s,n))
			return n;
		n++;
		if (n >= SET_MAX_SIZE(s))
			return -1;
	}

	while (s[n/ELEMENTSIZE]==0) {
		n+=ELEMENTSIZE;
		if (n >= SET_MAX_SIZE(s))
			return -1;
	}
	while (!SET_CONTAINS(s,n)) {
		n++;
		if (n >= SET_MAX_SIZE(s))
			return -1;
	}
	return n;
}


/*
 * set_print()
 *
 * Prints the size and contents of set s to stdout.
 * Mainly useful for debugging purposes and trivial output.
 */
UNUSED_FUNCTION
static void set_print(set_t s) {
	int i;
	printf("size=%d(max %d)",set_size(s),(int)SET_MAX_SIZE(s));
	for (i=0; i<SET_MAX_SIZE(s); i++)
		if (SET_CONTAINS(s,i))
			printf(" %d",i);
	printf("\n");
	return;
}

/********************************************************************
#include "graph.h"
*/

typedef struct _graph_t graph_t;
struct _graph_t {
	int n;             /* Vertices numbered 0...n-1 */
	set_t *edges;      /* A list of n sets (the edges). */
	int *weights;      /* A list of n vertex weights. */
};


#define GRAPH_IS_EDGE_FAST(g,i,j)  (SET_CONTAINS_FAST((g)->edges[(i)],(j)))
#define GRAPH_IS_EDGE(g,i,j) (((i)<((g)->n))?SET_CONTAINS((g)->edges[(i)], \
							  (j)):FALSE)
#define GRAPH_ADD_EDGE(g,i,j) do {            \
	SET_ADD_ELEMENT((g)->edges[(i)],(j)); \
	SET_ADD_ELEMENT((g)->edges[(j)],(i)); \
} while (FALSE)
#define GRAPH_DEL_EDGE(g,i,j) do {            \
	SET_DEL_ELEMENT((g)->edges[(i)],(j)); \
	SET_DEL_ELEMENT((g)->edges[(j)],(i)); \
} while (FALSE)


extern graph_t *graph_new(int n);
extern void graph_free(graph_t *g);
extern void graph_resize(graph_t *g, int size);
extern void graph_crop(graph_t *g);

extern boolean graph_weighted(graph_t *g);
extern int graph_edge_count(graph_t *g);

extern graph_t *graph_read_dimacs(FILE *fp);
extern graph_t *graph_read_dimacs_file(char *file);
extern boolean graph_write_dimacs_ascii(graph_t *g, char *comment,FILE *fp);
extern boolean graph_write_dimacs_ascii_file(graph_t *g,char *comment,
					     char *file);
extern boolean graph_write_dimacs_binary(graph_t *g, char *comment,FILE *fp);
extern boolean graph_write_dimacs_binary_file(graph_t *g, char *comment,
					      char *file);

extern void graph_print(graph_t *g);
extern boolean graph_test(graph_t *g, FILE *output);
extern int graph_test_regular(graph_t *g);

UNUSED_FUNCTION INLINE
static int graph_subgraph_weight(graph_t *g,set_t s) {
	int i,j;
	int count=0;
	setelement e;

	for (i=0; i<SET_ARRAY_LENGTH(s); i++) {
		if (s[i]) {
			e=s[i];
			for (j=0; j<ELEMENTSIZE; j++) {
				if (e&1)
					count+=g->weights[i*ELEMENTSIZE+j];
				e = e>>1;
			}
		}
	}
	return count;
}

UNUSED_FUNCTION INLINE
static int graph_vertex_degree(graph_t *g, int v) {
	return set_size(g->edges[v]);
}

/********************************************************************
#include "reorder.h"
*/

extern void reorder_set(set_t s,int *order);
extern void reorder_graph(graph_t *g, int *order);
extern int *reorder_duplicate(int *order,int n);
extern void reorder_invert(int *order,int n);
extern void reorder_reverse(int *order,int n);
extern int *reorder_ident(int n);
extern boolean reorder_is_bijection(int *order,int n);


#define reorder_by_default reorder_by_greedy_coloring
extern int *reorder_by_greedy_coloring(graph_t *g, boolean weighted);
extern int *reorder_by_weighted_greedy_coloring(graph_t *g, boolean weighted);
extern int *reorder_by_unweighted_greedy_coloring(graph_t *g,boolean weighted);
extern int *reorder_by_degree(graph_t *g, boolean weighted);
extern int *reorder_by_random(graph_t *g, boolean weighted);
extern int *reorder_by_ident(graph_t *g, boolean weighted);
extern int *reorder_by_reverse(graph_t *g, boolean weighted);


typedef struct _clique_options clique_options;
struct _clique_options {
	int *(*reorder_function)(graph_t *, boolean);
	int *reorder_map;

	/* arguments:  level, n, max, user_time, system_time, opts */
	boolean (*time_function)(int,int,int,int,double,double,
				 clique_options *);
	FILE *output;

	boolean (*user_function)(set_t,graph_t *,clique_options *);
	void *user_data;
	set_t *clique_list;
	int clique_list_length;
};

extern clique_options *clique_default_options;

/* Weighted clique functions */
extern int clique_max_weight(graph_t *g,clique_options *opts);
extern set_t clique_find_single(graph_t *g,int min_weight,int max_weight,
				boolean maximal, clique_options *opts);
extern int clique_find_all(graph_t *g, int req_weight, boolean exact,
			   boolean maximal, clique_options *opts);

/* Unweighted clique functions */
#define clique_unweighted_max_size clique_unweighted_max_weight
extern int clique_unweighted_max_weight(graph_t *g, clique_options *opts);
extern set_t clique_unweighted_find_single(graph_t *g,int min_size,
					   int max_size,boolean maximal,
					   clique_options *opts);
extern int clique_unweighted_find_all(graph_t *g, int min_size, int max_size,
				      boolean maximal, clique_options *opts);

/* Time printing functions */
extern boolean clique_print_time(int level, int i, int n, int max,
				 double cputime, double realtime,
				 clique_options *opts);
extern boolean clique_print_time_always(int level, int i, int n, int max,
					double cputime, double realtime,
					clique_options *opts);


/* Alternate spelling (let's be a little forgiving): */
#define cliquer_options clique_options
#define cliquer_default_options clique_default_options

/* Procedures defined in nautycliquer.c */
extern int find_clique(graph *g, int m, int n,
				     int min, int max, boolean maximal);
extern int find_indset(graph *g, int m, int n,
				     int min, int max, boolean maximal);

#endif /* !NAUTYCLIQUER_H */
