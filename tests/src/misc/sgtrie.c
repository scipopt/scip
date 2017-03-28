/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2016 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   select.c
 * @brief  unit tests for selection of unweighted and weighted median
 * @author Gregor Hendel
 */

#include "scip/pub_misc.h"
#include "scip/scip.h"
#include "scip/sgtrie.h"

#include "include/scip_test.h"

/** GLOBAL VARIABLES **/
static SCIP* scip;

/* TEST SUITE */
static
void setup(void)
{
   SCIPcreate(&scip);
}

static
void teardown(void)
{
   SCIPfree(&scip);
}

TestSuite(sgtrie, .init = setup, .fini = teardown);

/* TESTS  */
Test(sgtrie, create_and_free)
{
   /* calls setup and teardown */
}

typedef struct IntSet
{
   int elems[64];
   int nelems;
} INTSET;

static
void printSet(INTSET* set)
{
   int i;

   printf("{");
   for( i = 0; i < set->nelems; ++i )
   {
      if( i+1 == set->nelems)
         printf("%i", set->elems[i]);
      else
         printf("%i, ", set->elems[i]);
   }
   printf("}\n");
}

static
SCIP_DECL_SETCMP(cmpSets)
{
   INTSET* x;
   INTSET* y;
   int i,j;

   x = (INTSET*) a;
   y = (INTSET*) b;

   i = 0;
   j = 0;
   switch(op)
   {
      case SCIP_SGTRIE_SETEQ:
         if( x->nelems != y->nelems )
            return FALSE;

         while( i < x->nelems && j < y->nelems )
         {
            if( x->elems[i++] != y->elems[j++] )
               return FALSE;
         }

         return TRUE;
      case SCIP_SGTRIE_SUBSET:
         if( x->nelems > y->nelems )
            return FALSE;

         while( i < x->nelems && j < y->nelems )
         {
            if( x->elems[i] < y->elems[j] )
               return FALSE;
            else if( x->elems[i] > y->elems[j] )
               ++j;
            else
            {
               ++i;
               ++j;
            }
         }

         if( i < x->nelems )
            return FALSE;

         return TRUE;
      case SCIP_SGTRIE_SUBSETPLUSONE:
      {
         int distance;

         if( x->nelems > y->nelems + 1 )
            return FALSE;

         distance = 0;

         while( i < x->nelems && j < y->nelems )
         {
            if( x->elems[i] < y->elems[j] )
            {
               ++distance;
               if( distance > 1 )
                  return FALSE;
               ++i;
            }
            else if( x->elems[i] > y->elems[j] )
               ++j;
            else
            {
               ++i;
               ++j;
            }
         }

         while( i < x->nelems )
         {
            ++distance;
            if( distance > 1 )
               return FALSE;
            ++i;
         }

         return TRUE;
      }
   }

   assert(FALSE);
   return FALSE;
}

static
SCIP_DECL_GETSIGNATURE(getSignature)
{
   INTSET* set;
   int k;
   uint64_t signature = 0;

   set = (INTSET*) a;

   for( k = 0; k < set->nelems; ++k )
   {
      UPDATE_SIGNATURE(signature, set->elems[k]);
   }

   return signature;
}

Test(sgtrie, first_test, .description = "tests basic functionality of sgtrie")
{
   SCIP_SGTRIE* sgtrie;
   int nmatches;
   INTSET* matches[7];
   INTSET a,b,c,d,e,f,g;
   INTSET* sets[] = {&a, &b, &c, &d, &e, &f, &g};
   int i;

   a.nelems = 4;
   a.elems[0] = 1;
   a.elems[1] = 2;
   a.elems[2] = 3;
   a.elems[3] = 4;

   b.nelems = 3;
   b.elems[0] = 1;
   b.elems[1] = 2;
   b.elems[2] = 3;

   c.nelems = 3;
   c.elems[0] = 2;
   c.elems[1] = 3;
   c.elems[2] = 4;

   d.nelems = 3;
   d.elems[0] = 1;
   d.elems[1] = 2;
   d.elems[2] = 4;

   e.nelems = 4;
   e.elems[0] = 2;
   e.elems[1] = 3;
   e.elems[2] = 4;
   e.elems[3] = 5;

   f.nelems = 3;
   f.elems[0] = 1;
   f.elems[1] = 4;
   f.elems[2] = 5;

   g.nelems = 3;
   g.elems[0] = 2;
   g.elems[1] = 4;
   g.elems[2] = 5;

   SCIP_CALL( SCIPsgtrieCreate(&sgtrie, SCIPblkmem(scip), SCIPbuffer(scip), getSignature, cmpSets) );

   for( i = 0; i < 7; ++i )
   {
      SCIP_CALL( SCIPsgtrieInsert(sgtrie, sets[i]) );
   }

   for( i = 0; i < 7; ++i )
   {
      int k;

      SCIP_CALL( SCIPsgtrieFindSubsets(sgtrie, sets[i], (void**) matches, &nmatches) );

      printf("queryset: ");
      printSet(sets[i]);
      printf("subset candidates:\n");
      for( k = 0; k < nmatches; ++k )
      {
         printf("  ");
         printSet(matches[k]);
      }

      SCIP_CALL( SCIPsgtrieFindSubsetsPlusOne(sgtrie, sets[i], (void**) matches, &nmatches) );

      printf("subset plus one candidates:\n");
      for( k = 0; k < nmatches; ++k )
      {
         printf("  ");
         printSet(matches[k]);
      }

      printf("\n");
   }

   SCIPsgtrieFree(&sgtrie);
}
