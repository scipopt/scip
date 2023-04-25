/*  This file contains a 64-bit random number generator written
by the late George Marsaglia, published online at
https://www.thecodingforums.com/threads/64-bit-kiss-rngs.673657/

Since George is no longer with us, I can't get his permission to
distribute this with nauty, but since he published it in an open
forum without mentioning any restrictions I'm confident that he
would be pleased to see his effort put to good use.

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

#include <stdio.h>
#include <unistd.h>
#include <sys/random.h>

#define KISSX (1234567890987654321ULL)
#define KISSY (362436362436362436ULL)
static unsigned long long
kissx=KISSX, kissc=123456123456123456ULL,
kissy=KISSY,kissz=1066149217761810ULL,kisst;

#define MWC (kisst=(kissx<<58)+kissc, kissc=(kissx>>6), kissx+=kisst, \
             kissc+=(kissx<kisst), kissx)
#define XSH ( kissy^=(kissy<<13), kissy^=(kissy>>17), kissy^=(kissy<<43) )
#define CNG ( kissz=6906969069LL*kissz+1234567 )
#define KISS (MWC+XSH+CNG)

static long
ran_init_time_kiss(void)
/* Added by BDM: use the real time and getrandom() to initialize. */
{
    double t;
    unsigned long long ul;
    long seed;
    int i;
    REALTIMEDEFS

    t = NAUTYREALTIME;
    if (t > 1660000000.0) ul = (nauty_counter)(t*2100001.0);
    else                  ul = (nauty_counter)(t+212300021.0);
    kissx = KISSX + ul;
    i = getrandom(&ul,sizeof(ul),0);
    kissy = KISSY + (ul ^ getpid() * 997);

    for (i = 0; i < 1000; ++i) ul = KISS;

    return ul;
}

static long
ran_init_kiss(int seed)
/* Added by BDM: deterministically initialise using seed */
{
    double t;
    unsigned long long ul;
    int i;

    ul = 5003ULL * (unsigned long)seed;
 
    kissx = KISSX + ul;
    kissy = KISSY ^ ul;
    kissc = 123456123456123456ULL;
    kissz = 1066149217761810ULL;

    for (i = 0; i < 1000; ++i) ul = KISS;

    return (long)ul;
}

#undef KRAN
#define KRAN(mod) (KISS % (mod))
#undef INITRANBYTIME
#define INITRANBYTIME ran_init_time_kiss()
#define RAN_INIT(seed) ran_init_kiss(seed)
