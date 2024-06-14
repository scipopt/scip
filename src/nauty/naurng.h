/* naurng.h : header for George Marsaglia's  64-bit random number generator.

   To use it:
     1.  Call ran_init(seed) with any long long seed.  (Optional,
	   but you will always get the same sequence otherwise.)
         OR
         Call ran_init_time(extra) with long long extra seed.
         In this case the real time clock is used in addition to
         the extra seed you provide.
     
     2.  Use NEXTRAN to get the next number (0..2^64-1).
         All 64-bit unsigned values are possible.
         Alternatively, use KRAN(k) to get a random number 0..k-1.
	 For very large k, KRAN(k) is not quite uniform.  In that case
         use GETKRAN(k,var) to set the variable var to a better
         random number 0..k-1.  Since NEXTRAN gibes a full 64-bit
         integer, it is not worth using this slower version unless
         k is bigger than a billion or so.
*/

#ifndef NAURNG_H
#define NAURNG_H
#include "naututil.h"

#ifdef __cplusplus
extern "C" {
#endif

extern void ran_init(unsigned long long);
extern void ran_init_2(unsigned long long, unsigned long long);
extern unsigned long long ran_init_time(unsigned long long);
extern unsigned long long ran_nextran(void);

#ifdef __cplusplus
}
#endif

#define MAXRAN (~0ULL)    /* Values are 0..MAXRAN */
#define NEXTRAN (ran_nextran())
#define KRAN(k) (NEXTRAN%(k))
#define RANREAL ((NEXTRAN+0.5)/(MAXRAN+1.0))  /* Uniform (0,1) */

#define MAXSAFE(k) (MAXRAN - MAXRAN%(k))
#define GETKRAN(k,var) {unsigned long long __getkran; \
    do {__getkran = NEXTRAN;} while (__getkran >= MAXSAFE(k)); \
    var = __getkran % (k);}
#define INITRANBYTIME ran_init_time(0)

#endif
