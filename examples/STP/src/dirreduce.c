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
int degree_test_dir(
   GRAPH*  g,
   double* fixed)
{
   int i;
   int i1;
   int i2;
   int e1;
   int e2;
   int rerun = TRUE;
   int done  = TRUE;
   int count = 0;

   assert(g      != NULL);
   assert(fixed  != NULL);

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


         if (g->grad[i] == 1 && g->cost[g->inpbeg[i]] < FARAWAY)
         {
            e1  = g->inpbeg[i];
            i1  = g->tail[e1];

            assert(e1 >= 0);
            assert(e1 == Edge_anti(g->outbeg[i]));
            assert(g->ieat[e1] == EAT_LAST);
            assert(g->oeat[g->outbeg[i]] == EAT_LAST);

            if (Is_term(g->term[i]))
               *fixed += g->cost[e1];

            graph_knot_contract(g, i1, i);

            assert(g->grad[i] == 0);

            /* Ist es etwa der Letzte gewesen ?
             */
            if (g->grad[i1] == 0)
            {
               rerun = FALSE;
               break;
            }
            if ((i1 < i) && (g->grad[i1] < 3))
               rerun = TRUE;

            count++;

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

            do
            {
               done = TRUE;

               if (!Is_term(g->term[i]))
               {
                  /* both the edges are outgoing from node i
                   * need to ensure that the flow of the edge costs is correct
                   * Edge_anti(e2) -> e1 and Edge_anti(e1) -> e2  */
                  g->cost[e1]            += g->cost[Edge_anti(e2)];
                  g->cost[Edge_anti(e1)] += g->cost[e2];

                  graph_knot_contract(g, i2, i);

                  count++;

                  break;
               }
               assert(Is_term(g->term[i]));

               if (g->cost[g->outbeg[i]] != g->cost[g->inpbeg[i]])
                  continue;

               /* At present the degree two test for terminals does not work.
                * Needs to be researched and updated. */
               if (Is_term(g->term[i1]) && Is_term(g->term[i2]))
               {
                  if (LT(g->cost[e1], g->cost[e2]))
                  {
                     *fixed += g->cost[e1];
                     graph_knot_contract(g, i1, i);
                  }
                  else
                  {
                     *fixed += g->cost[e2];
                     graph_knot_contract(g, i2, i);
                  }
                  count++;

                  break;
               }
               if (Is_term(g->term[i1]) && !Is_term(g->term[i2]) && LE(g->cost[e1], g->cost[e2]))
               {
                  *fixed += g->cost[e1];
                  graph_knot_contract(g, i1, i);

                  count++;

                  break;
               }
               if (Is_term(g->term[i2]) && !Is_term(g->term[i1]) && LE(g->cost[e2], g->cost[e1]))
               {
                  *fixed += g->cost[e2];
                  graph_knot_contract(g, i2, i);

                  count++;

                  break;
               }
               done = FALSE;
            }
            /* CONSTCOND */
            /*lint -save -e717 */
            while(FALSE);
            /*lint -restore */

            if (done
                  && (((i1 < i) && (g->grad[i1] < 3))
                     || ((i2 < i) && (g->grad[i2] < 3))))
               rerun = TRUE;
         }
      }
   }
   SCIPdebugMessage(" %d Knots deleted\n", count);

   assert(graph_valid(g));

   return(count);
}
