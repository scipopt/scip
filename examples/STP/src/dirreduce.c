/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*   Type....: Function                                                      */
/*   File....: dirreduce.c                                                   */
/*   Name....: Steiner Tree Reduction                                        */
/*   Author..: Thorsten Koch                                                 */
/*   Copyright by Author, All rights reserved                                */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*lint -esym(750,REDUCE_C) -esym(766,stdlib.h) -esym(766,string.h)           */

#define DIRREDUCE_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "grph.h"
#include "portab.h"
#include "scip/scip.h"

/* Moeglichkeiten:
 *    a 1 b 2 c
 *    *---*       : Kontraction entlang 1, a -> b
 *    *---o---*   : Kontraktion entlang 2, b -> c
 *    t---t---t   : Kontraktion entlang min(1, 2), b -> min(d(a),d(c))
 *    o---t---o   : Nichts
 *    t---t---o   : Kontraktion entlang 1, b -> a, wenn c(1) <= c(2)
 *    o---t---t   : Kontraktion entlang 2, b -> c, wenn c(2) <= c(1)
 */
SCIP_RETCODE degree_test_dir(
   SCIP* scip,
   GRAPH*  g,
   double* fixed,
   int*    count
   )
{
   int i;
   int i1;
   int i2;
   int e1;
   int e2;
   int rerun = TRUE;

   assert(g      != NULL);
   assert(fixed  != NULL);
   assert(count != NULL);

   *count = 0;
   SCIPdebugMessage("Degree Test: ");
   fflush(stdout);

   while(rerun)
   {
      rerun = FALSE;

      SCIPdebug(fputc('.', stdout));
      SCIPdebug(fflush(stdout));

      for(i = 0; i < g->knots; i++)
      {
         assert(g->grad[i] >= 0);
         if( g->grad[i] == 1 && g->cost[g->outbeg[i]] < FARAWAY && g->cost[g->inpbeg[i]] < FARAWAY )
         {
            e1  = g->inpbeg[i];
            i1  = g->tail[e1];
            assert(e1 >= 0);
            assert(e1 == Edge_anti(g->outbeg[i]));
            assert(g->ieat[e1] == EAT_LAST);
            assert(g->oeat[g->outbeg[i]] == EAT_LAST);

            if (Is_term(g->term[i]))
	    {
	       SCIP_CALL( SCIPindexListNodeAppendCopy(scip, &(g->fixedges), g->ancestors[e1]) );
               *fixed += g->cost[e1];
	    }

            SCIP_CALL( graph_knot_contract(scip, g, i1, i) );

            assert(g->grad[i] == 0);

            /* the last node? */
            if (g->grad[i1] == 0)
            {
               rerun = FALSE;
               break;
            }
            if ((i1 < i) && (g->grad[i1] < 3))
               rerun = TRUE;

            (*count)++;

            continue;
         }

         /* Note on costs in the directed graph
          * g->outbeg[i] is the outgoing directed edge
          * g->inpbeg[i] is the incoming directed edge */
         if (g->grad[i] == 2)
         {
            e1 = g->outbeg[i];
            e2 = g->oeat[e1];
            i1 = g->head[e1];
            i2 = g->head[e2];

            assert(e1 >= 0);
            assert(e2 >= 0);


            if (!Is_term(g->term[i]) )
            {
               /* both the edges are outgoing from node i
                * need to ensure that the flow of the edge costs is correct
                * Edge_anti(e2) -> e1 and Edge_anti(e1) -> e2  */
               if( !Is_term(g->term[i2]) && !Is_term(g->term[i1]) )
               {
                  g->cost[e1] += g->cost[Edge_anti(e2)];
                  g->cost[Edge_anti(e1)] += g->cost[e2];
                  SCIP_CALL( graph_knot_contract(scip, g, i2, i) );
                  (*count)++;

                  if (!Is_term(g->term[i2]) && (((i1 < i) && (g->grad[i1] < 3))
                        || ((i2 < i) && (g->grad[i2] < 3))))
                     rerun = TRUE;

               }

            }

            /* At present the degree two test for terminals does not work.
             * Needs to be researched and updated. */
#if 0
            if (Is_term(g->term[i1]) && Is_term(g->term[i2]))
            {
               if (LT(g->cost[e1], g->cost[e2]))
               {
                  SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e1]);
                  *fixed += g->cost[e1];
                  graph_knot_contract(g, i1, i);
               }
               else
               {
                  SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e2]);
                  *fixed += g->cost[e2];
                  graph_knot_contract(g, i2, i);
               }
               count++;

               break;
            }
            if (Is_term(g->term[i1]) && !Is_term(g->term[i2]) && LE(g->cost[e1], g->cost[e2]))
            {
               SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e1]);
               *fixed += g->cost[e1];
               graph_knot_contract(g, i1, i);

               count++;

               break;
            }
            if (Is_term(g->term[i2]) && !Is_term(g->term[i1]) && LE(g->cost[e2], g->cost[e1]))
            {
               SCIPindexListNodeAppendCopy(&(g->fixedges), g->ancestors[e2]);
               *fixed += g->cost[e2];
               graph_knot_contract(g, i2, i);

               count++;

               break;
            }
#endif

            /* CONSTCOND */
            /*lint -save -e717 */

            /*lint -restore */
         }
      }
   }
   SCIPdebugMessage(" %d Knots deleted\n", *count);
   printf("dirdeg %d Knots deleted\n", *count);
   assert(graph_valid(g));

   return SCIP_OKAY;
}
