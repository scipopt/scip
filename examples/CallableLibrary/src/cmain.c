/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   cmain.c
 * @brief  main file for C compilation
 * @author Stefan Vigerske
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "scip/scip.h"

extern SCIP_RETCODE runString(void);

extern SCIP_RETCODE runGastrans(void);

extern SCIP_RETCODE runCircle(void);

/** main method starting SCIP */
int main(
   int                        argc,          /**< number of arguments from the shell */
   char**                     argv           /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode;
   int i;

   assert(argc >= 1);
   if( argc == 1 )
   {
      printf("Usage: %s { string | gastrans | circle }\n", argv[0]);
      return 0;
   }

   for( i = 1; i < argc; ++i )
   {
      if( strcmp(argv[i], "string") == 0 )
      {
         retcode = runString();
      }
      else if( strcmp(argv[i], "gastrans") == 0 )
      {
         retcode = runGastrans();
      }
      else if( strcmp(argv[i], "circle") == 0 )
      {
         retcode = runCircle();
      }
      else
      {
         printf("Example '%s' unknown.\n", argv[i]);
         return 1;
      }

      /* evaluate return code of the SCIP process */
      if( retcode != SCIP_OKAY )
      {
         /* write error back trace */
         SCIPprintError(retcode);
         return -1;
      }
   }

   return 0;
}
