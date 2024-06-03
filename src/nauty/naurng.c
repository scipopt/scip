/*  This file contains a 64-bit random number generator written
by the late George Marsaglia, published online at
https://www.thecodingforums.com/threads/64-bit-kiss-rngs.673657/

Since George is no longer with us, I can't get his permission to
distribute this with nauty, but since he published it in an open
forum without mentioning any restrictions I'm confident that he
would be pleased to see his effort put to good use.

See naurng.h for usage instructions.

The rest of these comments are from George.

Use of KISS or KISS() as a general 64-bit RNG requires specifying
3*64+58=250 bits for seeds, 64 bits each for x,y,z and 58 for c,
resulting in a composite sequence with period around 2^250.
The actual period is
(2^250+2^192+2^64-2^186-2^129)/6 ~= 2^(247.42) or 10^(74.48).
We "lose" 1+1.58=2.58 bits from maximum possible period, one bit
because b=2^64, a square, cannot be a primitive root of p=ab-1,
so the best possible order for b is (p-1)/2.
The periods of MWC and XSH have gcd 3=2^1.58, so another 1.58
bits are "lost" from the best possible period we could expect
from 250 seed bits.

Some users may think 250 seed bits are an unreasonable requirement.
A good seeding procedure might be to assume the default seed
values then let the user choose none, one, two,..., or all
of x,y,z, and c to be reseeded.
*/

#include "naurng.h"

#define KISSX 1234567890987654321ULL
#define KISSY 362436362436362436ULL
#define KISSC 123456123456123456ULL
#define KISSZ 1066149217761810ULL
static TLS_ATTR unsigned long long kissx=KISSX, kissc=KISSC,
                                   kissy=KISSY, kissz=KISSZ, kisst;

#define MWC (kisst=(kissx<<58)+kissc, kissc=(kissx>>6), kissx+=kisst, \
             kissc+=(kissx<kisst), kissx)
#define XSH ( kissy^=(kissy<<13), kissy^=(kissy>>17), kissy^=(kissy<<43) )
#define CNG ( kissz=6906969069LL*kissz+1234567 )
#define KISS (MWC+XSH+CNG)

void
ran_init(unsigned long long seed)
/* Initialize random number generator */
{
   ran_init_2(seed,0);
}

void
ran_init_2(unsigned long long seed1, unsigned long long seed2)
/* Use the two seeds to initialize */
{
    int i;
    unsigned long long ul;

    kissx = KISSX + seed1;
    kissy = KISSY + (seed2 * 997);
    kissc = KISSC;
    kissz = KISSZ;

    for (i = 0; i < 1000; ++i) ul = KISS;
}

unsigned long long
ran_init_time(unsigned long long extra)
/* Use the real time and an extra seed to initialize. If val is
   the value returned, then ran_init_2(val,extra) will perform
   the same initialization. */
{
    double t;
    unsigned long long ul,seed;
    REALTIMEDEFS

    t = NAUTYREALTIME;
    if (t > 1660000000.0) seed = (unsigned long long)(t*2100001.0);
    else                  seed = (unsigned long long)(t+212300021.0);

    ran_init_2(seed,extra);

    return seed;
}

unsigned long long
ran_nextran()
/* Make a 64-bit random number */
{
    return KISS;
}
