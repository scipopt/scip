/* gutil3.c: Some more graph utilities, using sparse representation. */

#include "gtools.h"
#include "gutils.h"

#if MAXN
static TLS_ATTR int work1[MAXN];
static TLS_ATTR int work2[MAXN];
#else
DYNALLSTAT(int,work1,work1_sz);
DYNALLSTAT(int,work2,work2_sz);
#endif

/*****************************************************************************
*                                                                            *
*  diameter_sg(*sg, *v1, *v2)                                                *
*  For a connected graph, return the diamater, and if v1,v2!=NULL give an    *
*  example of two vertices at that distance.  For a disconnected graph,      *
*  give that information for a component of maximum diameter.                *
*                                                                            *
*****************************************************************************/

int
diameter_sg(sparsegraph *sg, int *v1, int *v2)
{
    int n,*d,*e;
    int v0,i,head,tail;
    int di,k,eg1,eg2,thisv,maxd;
    size_t *v,vi,j;

#define DIST work1
#define QUEUE work2

    n = sg->nv;
#if !MAXN
    DYNALLOC1(int,work1,work1_sz,n,"diameter_sg");
    DYNALLOC1(int,work2,work2_sz,n,"diameter_sg");
#endif

    maxd = eg1 = eg2 = 0;
    SG_VDE(sg,v,d,e);

    for (v0 = 0; v0 < n; ++v0)
    {
        if (d[v0] == 0) continue;
        for (i = 0; i < n; ++i) DIST[i] = n;
        DIST[v0] = 0;
        QUEUE[0] = v0;

        head = 0;
        tail = 1;
        thisv = v0;
        while (tail < n && head < tail)
        {
            i = QUEUE[head++];
            vi = v[i];
            di = d[i];
            for (j = 0; j < di; ++j)
            {
                k = e[vi+j];
                if (DIST[k] == n)
                {
                    thisv = k;
                    DIST[k] = DIST[i] + 1;
                    QUEUE[tail++] = k;
                }
            }
        }

        if (DIST[thisv] > maxd)
        {
            maxd = DIST[thisv];
            eg1 = v0;
            eg2 = thisv;
        }
    }

    if ((v1)) *v1 = eg1;
    if ((v2)) *v2 = eg2;
    return maxd;
}
