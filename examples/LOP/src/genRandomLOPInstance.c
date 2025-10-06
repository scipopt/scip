/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright (c) 2002-2025 Zuse Institute Berlin (ZIB)                      */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   genRandomLOPInstance.c
 * @brief  generate a random linear ordering problem instance
 * @author Marc Pfetsch
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


/* the following is copied from an old version of scip/misc.c */

#ifdef _WIN32
#define NO_RAND_R
#endif

/* define own random numbers or take library version depending on the following define */
#ifdef NO_RAND_R

#define SCIP_RAND_MAX 32767

/** returns a random number between 0 and SCIP_RAND_MAX */
static
int getRand(
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   SCIP_Longint nextseed;

   assert(seedp != NULL);

   nextseed = (*seedp) * 1103515245 + 12345;
   *seedp = (unsigned int)nextseed;

   return (int)((unsigned int)(nextseed/(2*(SCIP_RAND_MAX+1))) % (SCIP_RAND_MAX+1));
}

#else

#define SCIP_RAND_MAX RAND_MAX

/** returns a random number between 0 and SCIP_RAND_MAX */
static
int getRand(
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return rand_r(seedp);
}

#endif

/** returns a random integer between minrandval and maxrandval */
static
int getRandomInt(
   int                   minrandval,         /**< minimal value to return */
   int                   maxrandval,         /**< maximal value to return */
   unsigned int*         seedp               /**< pointer to seed value */
   )
{
   return minrandval + (int) ((maxrandval - minrandval + 1)*(double)getRand(seedp)/(SCIP_RAND_MAX+1.0));
}

int main(int argc, char** argv)
{
   int n;
   int d;
   int i;
   int j;
   unsigned int seed;
   FILE *file;

   if ( argc != 4 )
   {
      printf("usage: %s <filename> <n> <d>.\n", argv[0]);
      return 1;
   }

   n = atoi(argv[2]);
   d = atoi(argv[3]);
   seed = 0;
   assert( n > 0 );
   assert( d > 0 );

   /* open file */
   file = fopen(argv[1], "w");
   if ( file == NULL )
   {
      printf("Could not open file %s.\n", argv[1]);
      return 1;
   }

   /* write comment line and size*/
   fprintf(file, "Randomly generated LOP instance.\n");
   fprintf(file, "%d\n", n);
   for (i = 0; i < n; ++i)
   {
      for (j = 0; j < n; ++j)
	 fprintf(file, "%d ", getRandomInt(0, d, &seed));
      fprintf(file, "\n");
   }
   fclose(file);

   printf("Wrote random LOP instance to %s\n", argv[1]);
   printf("Size: %d\n", n);
   printf("Entries: {0, ..., %d}\n", d);

   return 0;
}
