/**@file   genRandomLOPInstance.c
 * @brief  generate a random linear ordering problem instance
 * @author Marc Pfetsch
 */

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

int main(int argc, char** argv)
{
   int n, d, i, j;
   FILE *file;

   if ( argc != 4 )
   {
      printf("usage: %s <filename> <n> <d>.\n", argv[0]);
      return 1;
   }

   n = atoi(argv[2]);
   d = atoi(argv[3]);
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
	 fprintf(file, "%d ", (int) (((double)d) * drand48()));
      fprintf(file, "\n");
   }

   printf("Wrote random LOP instance to %s\n", argv[1]);
   printf("Size: %d\n", n);
   printf("Entries: {0, ..., %d}\n", d);
}
