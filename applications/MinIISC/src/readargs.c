/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2015 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   readargs.c
 * @brief  read comand line arguments
 * @author Marc Pfetsch
 */

#include "readargs.h"


/** get problem name
 *
 *  Returns 0 on error
 */
int getProblemName(
   const char*           filename,           /**< file name */
   char*                 probname,           /**< name of problem (output) */
   int                   maxsize             /**< max length of problem name */
   )
{
   int i = 0;
   int result = 1;
   int l;
   int j;

   /* first find end of string */
   while ( filename[i] != 0) /* find end of string */
      ++i;
   l = i;                /* end of string */

   /* if we found ".gz" */
   if ( l > 3 && filename[l-3] == '.' && filename[l-2] == 'g' && filename[l-1] == 'z' )
   {
      l = l - 4;
      i = l;
   }

   /* go back until '.' or '/' appears */
   while ((i > 0) && (filename[i] != '.') && (filename[i] != '/'))
      --i;
   assert(i > 0);

   /* if we found '.', search for '/' */
   if (filename[i] == '.')
   {
      l = i;
      while ((i > 0) && (filename[i] != '/'))
	 --i;
   }

   /* correct counter */
   if (filename[i] == '/')
      ++i;

   /* copy name */
   j = 0;
   while ( (i < l) && (filename[i] != 0) )
   {
      probname[j++] = filename[i++];
      if ( j >= maxsize-1 )
      {
	 result = 0;
	 break;
      }
   }
   probname[j] = 0;

   return result;
}


/** read comand line arguments */
SCIP_RETCODE readArguments(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char**          filename,           /**< file name from arguments */
   const char**          settingsname,       /**< name of settings file */
   SCIP_Real*            timelimit,          /**< time limit read from arguments */
   SCIP_Real*            memlimit,           /**< memory limit read from arguments */
   SCIP_Longint*         nodelimit,          /**< node limit read from arguments */
   int*                  dispfreq            /**< display frequency */
   )
{
   int i;
   char usage[256];
#ifndef NDEBUG
   int status;
#endif

   assert( argc > 0 );
   assert( argv != NULL );
   assert( filename != NULL );
   assert( settingsname != NULL );
   assert( timelimit != NULL );
   assert( memlimit != NULL );
   assert( nodelimit != NULL );
   assert( dispfreq != NULL );

   /* init usage text */
#ifndef NDEBUG
   status = snprintf(usage, 256, "usage: %s <file> [-s <setting file>] [-t <time limit>] [-m <mem limit>] [-n <node limit>] [-d <display frequency>]", argv[0]);
   assert( 0 <= status && status < 256 );
#else
   (void) snprintf(usage, 256, "usage: %s <file> [-s <setting file>] [-t <time limit>] [-m <mem limit>] [-n <node limit>] [-d <display frequency>]", argv[0]);
#endif

   /* init arguments */
   *timelimit = 1e20;
   *memlimit = 1e20;
   *nodelimit = SCIP_LONGINT_MAX;
   *filename = NULL;
   *settingsname = NULL;
   *dispfreq = -1;

   /* check all arguments */
   for (i = 1; i < argc; ++i)
   {
      /* check for settings */
      if ( ! strcmp(argv[i], "-s") )
      {
         if ( *settingsname != NULL )
         {
	    printf("Error: %s\n", usage);
            return SCIP_ERROR;
         }
         if ( i == argc-1 )
         {
            printf("No setting file name supplied.\n");
            printf("%s\n", usage);
            return SCIP_ERROR;
         }
         ++i;
         *settingsname = argv[i];
         assert( i < argc );
      }
      else
      {
         /* check for time limit */
         if ( ! strcmp(argv[i], "-t") )
         {
            if ( i == argc-1 )
            {
               printf("No time limit supplied.\n");
               printf("%s\n", usage);
               return SCIP_ERROR;
            }
            ++i;
            *timelimit = atof(argv[i]);
            assert( i < argc );
         }
         else
         {
            /* check for memory limit */
            if ( ! strcmp(argv[i], "-m") )
            {
               if ( i == argc-1 )
               {
                  printf("No memory limit supplied.\n");
                  printf("%s\n", usage);
                  return SCIP_ERROR;
               }
               ++i;
               *memlimit = atof(argv[i]);
               assert( i < argc );
            }
            else
            {
               /* check for memory limit */
               if ( ! strcmp(argv[i], "-n") )
               {
                  if ( i == argc-1 )
                  {
                     printf("No node limit supplied.\n");
                     printf("%s\n", usage);
                     return SCIP_ERROR;
                  }
                  ++i;
                  *nodelimit = atol(argv[i]);
                  assert( i < argc );
               }
               else
               {
                  /* check for display frequency */
                  if ( ! strcmp(argv[i], "-d") )
                  {
                     if ( i == argc-1 )
                     {
                        printf("No display frequency supplied.\n");
                        printf("%s\n", usage);
                        return SCIP_ERROR;
                     }
                     ++i;
                     *dispfreq = atoi(argv[i]);
                     assert( i < argc );
                  }
                  else
                  {
                     /* if filename is already specified */
                     if ( *filename != NULL )
                     {
                        printf("%s\n", usage);
                        return SCIP_ERROR;
                     }
                     *filename = argv[i];
                  }
               }
	    }
	 }
      }
   }

   if ( *filename == NULL )
   {
      printf("No filename supplied.\n");
      printf("%s\n", usage);
      return SCIP_ERROR;
   }

   return SCIP_OKAY;
}
