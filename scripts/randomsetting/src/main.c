/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2021 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   main.c
 * @brief  random settings generator
 * @author Ambros Gleixner
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "scip/scip.h"
#include "scip/scipgithash.h"
#include "scip/scipdefplugins.h"

#define MAX_PARAMNAME_LEN                255 /**< maximum length of a parameter name */
#define MAX_CHANGE                     10000 /**< maximum number of parameters changed */
#define CHANGE_ADVANCED                 TRUE /**< should advanced parameters also be changed? */


/** checks whether parameter is in the list of excluded parameters
 *
 *  @note this should be parameters for which SCIP most of the time aborts with a proper error in optimized mode;
 *        ideally there should be no random parameter settings which cause an assert in debug mode and undefined
 *        behaviour in optimized mode
 */
static
SCIP_Bool paramIsExcluded(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAM*           param               /**< pointer to parameter */
   )
{
   const char* paramname;
   SCIP_Bool excluded;

   assert(scip != NULL);
   assert(param != NULL);

   paramname = SCIPparamGetName(param);
   excluded = FALSE;

   /* otherwise the readers do not work */
   excluded = excluded
      || strncmp(paramname, "misc/usevartable", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "misc/useconstable", MAX_PARAMNAME_LEN) == 0;

   /* too large initial values easily yield a memory out */
   excluded = excluded
      || strncmp(paramname, "memory/arraygrowinit", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "memory/treegrowinit", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "memory/pathgrowinit", MAX_PARAMNAME_LEN) == 0;

   /* already larger values inside the feasible interval [1,10] yield many memory outs */
   excluded = excluded
      || strncmp(paramname, "memory/arraygrowfac", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "memory/treegrowfac", MAX_PARAMNAME_LEN) == 0
      || strncmp(paramname, "memory/pathgrowfac", MAX_PARAMNAME_LEN) == 0;

   return excluded;
}

/** sets each parameter to a random value within its domain */
static
SCIP_RETCODE generateRandomSettings(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RANDNUMGEN*      randgen,            /**< random number generator */
   int                   maxchange,          /**< maximum number of settings to change */
   SCIP_Bool             advanced            /**< should advanced parameters also be changed? */
   )
{
   SCIP_PARAM** params;
   int nparams;
   int nchangeparams;
   int nchanged;
   int start;
   int step;
   int i;

   assert(scip != NULL);
   assert(randgen != NULL);

   params = SCIPgetParams(scip);
   nparams = SCIPgetNParams(scip);
   if( nparams <= 0 )
      return SCIP_OKAY;

   start = SCIPrandomGetInt(randgen, 0, nparams - 1);
   step = SCIPrandomGetInt(randgen, 0, nparams - 1);
   nchanged = 0;

   assert(params != NULL);
   assert(nparams > 0);
   assert(start >= 0);
   assert(start < nparams);
   assert(step >= 0);
   assert(step < nparams);

   /* make sure all parameters are at their default values */
   SCIP_CALL( SCIPresetParams(scip) );

   /* count number of changable parameters */
   nchangeparams = 0;

   for( i = 0; i < nparams; ++i )
   {
      SCIP_PARAM* param = params[i];
      assert(param != NULL);

      if( !advanced && SCIPparamIsAdvanced(param) )
      {
         SCIPdebugMessage("excluding advanced parameter <%s>\n", SCIPparamGetName(param));
      }
      else if( SCIPparamGetType(param) == SCIP_PARAMTYPE_STRING )
      {
         SCIPdebugMessage("excluding string parameter <%s>\n", SCIPparamGetName(param));
      }
      else if( paramIsExcluded(scip, param) )
      {
         SCIPdebugMessage("excluding parameter <%s>\n", SCIPparamGetName(param));
      }
      else
         nchangeparams++;
   }

   SCIPdebugMessage("counted %d changable parameters (of %d)\n", nchangeparams, nparams);

   /* reduce limit if there are not enough changable parameters */
   maxchange = MIN(maxchange, nchangeparams);

   /* loop through parameters; if a parameter has been changed, we make a big step, otherwise we just increment */
   for( i = start; nchanged < maxchange; i = (i + 1) % nparams )
   {
      SCIP_PARAM* param;

      assert(i >= 0);
      assert(i < nparams);

      param = params[i];
      assert(param != NULL);

      if( advanced == FALSE && SCIPparamIsAdvanced(param) )
         continue;

      if( SCIPparamGetType(param) == SCIP_PARAMTYPE_STRING )
         continue;

      if( paramIsExcluded(scip, param) )
         continue;

      switch( SCIPparamGetType(param) )
      {
      case SCIP_PARAMTYPE_BOOL:
         {
            SCIP_RETCODE retcode;
            int value;

            /* if the value is already changed, continue looking for an unchanged parameter */
            if( SCIPparamGetBool(param) != SCIPparamGetBoolDefault(param) )
               continue;

            value = 1 - SCIPparamGetBool(param);
            assert(value >= 0);
            assert(value <= 1);

            retcode = SCIPchgBoolParam(scip, param, (SCIP_Bool)value);
            if( retcode != SCIP_OKAY )
            {
               SCIP_CALL( retcode );
            }
            else
            {
               SCIPdebugMessage("changing bool parameter <%s> to %s\n", SCIPparamGetName(param),
                  SCIPparamGetBool(param) ? "TRUE" : "FALSE");
               nchanged++;
               i += step;
            }

            break;
         }

      case SCIP_PARAMTYPE_INT:
         {
            SCIP_RETCODE retcode;
            int oldvalue;

            /* if the value is already changed, continue looking for an unchanged parameter */
            if( SCIPparamGetInt(param) != SCIPparamGetIntDefault(param) )
               continue;

            /* if there is only one feasible value, we pretend it is changed */
            if( SCIPparamGetIntMin(param) == SCIPparamGetIntMax(param) )
            {
               SCIPdebugMessage("leaving int parameter <%s> at %d\n", SCIPparamGetName(param), SCIPparamGetInt(param));
               nchanged++;
               i += step;
               continue;
            }

            retcode = SCIP_OKAY;
            oldvalue = SCIPparamGetInt(param);

            /* try to change parameter until we reach a different value; give up if there is an error (other than
             * SCIP_PARAMETERWRONGVAL) or there is only one feasible value
             */
            while( ((retcode == SCIP_OKAY && SCIPparamGetInt(param) == oldvalue) || retcode == SCIP_PARAMETERWRONGVAL) )
            {
               int newvalue;

               newvalue = SCIPrandomGetInt(randgen, SCIPparamGetIntMin(param), SCIPparamGetIntMax(param));
               assert(newvalue >= SCIPparamGetIntMin(param));
               assert(newvalue <= SCIPparamGetIntMax(param));

               retcode = SCIPchgIntParam(scip, param, newvalue);
               if( retcode == SCIP_PARAMETERWRONGVAL )
               {
                  SCIPinfoMessage(scip, NULL, "could not set parameter to random value within range - trying again");
               }
            }

            if( retcode != SCIP_OKAY )
            {
               SCIP_CALL( retcode );
            }
            else
            {
               assert(SCIPparamGetInt(param) != oldvalue);

               SCIPdebugMessage("changing int parameter <%s> from %d to %d\n", SCIPparamGetName(param), oldvalue,
                  SCIPparamGetInt(param));
               nchanged++;
               i += step;
            }

            break;
         }

      case SCIP_PARAMTYPE_LONGINT:
         {
            SCIP_RETCODE retcode;
            SCIP_Longint oldvalue;

            /* if the value is already changed, continue looking for an unchanged parameter */
            if( SCIPparamGetLongint(param) != SCIPparamGetLongintDefault(param) )
               continue;

            /* if there is only one feasible value, we pretend it is changed */
            if( SCIPparamGetLongintMin(param) == SCIPparamGetLongintMax(param) )
            {
               SCIPdebugMessage("leaving longint parameter <%s> at %"SCIP_LONGINT_FORMAT"\n", SCIPparamGetName(param), SCIPparamGetLongint(param));
               nchanged++;
               i += step;
               continue;
            }

            retcode = SCIP_OKAY;
            oldvalue = SCIPparamGetLongint(param);

            /* try to change parameter until we reach a different value; give up if there is an error (other than
             * SCIP_PARAMETERWRONGVAL) or there is only one feasible value
             */
            while( ((retcode == SCIP_OKAY && SCIPparamGetLongint(param) == oldvalue) || retcode == SCIP_PARAMETERWRONGVAL) )
            {
               SCIP_Longint newvalue;
               int minvalue;
               int maxvalue;

               minvalue = (SCIPparamGetLongintMin(param) >= (SCIP_Longint)INT_MIN)
                  ? (int)SCIPparamGetLongintMin(param)
                  : INT_MIN;
               maxvalue = (SCIPparamGetLongintMax(param) <= (SCIP_Longint)INT_MAX)
                  ? (int)SCIPparamGetLongintMax(param)
                  : INT_MAX;
               newvalue = (SCIP_Longint)SCIPrandomGetInt(randgen, minvalue, maxvalue);
               assert(newvalue >= SCIPparamGetLongintMin(param));
               assert(newvalue <= SCIPparamGetLongintMax(param));

               retcode = SCIPchgLongintParam(scip, param, newvalue);
               if( retcode == SCIP_PARAMETERWRONGVAL )
               {
                  SCIPinfoMessage(scip, NULL, "could not set parameter to random value within range - trying again");
               }
            }

            if( retcode != SCIP_OKAY )
            {
               SCIP_CALL( retcode );
            }
            else
            {
               assert(SCIPparamGetLongint(param) != oldvalue);

               SCIPdebugMessage("changing longint parameter <%s> from %"SCIP_LONGINT_FORMAT" to %"SCIP_LONGINT_FORMAT"\n",
                  SCIPparamGetName(param), oldvalue, SCIPparamGetLongint(param));

               nchanged++;
               i += step;
            }

            break;
         }

      case SCIP_PARAMTYPE_REAL:
         {
            SCIP_RETCODE retcode;
            SCIP_Real oldvalue;

            /* if the value is already changed, continue looking for an unchanged parameter */
            if( SCIPparamGetReal(param) != SCIPparamGetRealDefault(param) )
               continue;

            /* if there is only one feasible value, we pretend it is changed */
            if( SCIPparamGetRealMin(param) == SCIPparamGetRealMax(param) )
            {
               SCIPdebugMessage("leaving real parameter <%s> at %g\n", SCIPparamGetName(param), SCIPparamGetReal(param));
               nchanged++;
               i += step;
               continue;
            }

            retcode = SCIP_OKAY;
            oldvalue = SCIPparamGetReal(param);

            /* try to change parameter until we reach a different value; give up if there is an error (other than
             * SCIP_PARAMETERWRONGVAL) or there is only one feasible value
             */
            while( ((retcode == SCIP_OKAY && SCIPparamGetReal(param) == oldvalue) || retcode == SCIP_PARAMETERWRONGVAL) )
            {
               SCIP_Real newvalue;

               newvalue = SCIPrandomGetReal(randgen, SCIPparamGetRealMin(param), SCIPparamGetRealMax(param));
               assert(newvalue >= SCIPparamGetRealMin(param));
               assert(newvalue <= SCIPparamGetRealMax(param));

               retcode = SCIPchgRealParam(scip, param, newvalue);
               if( retcode == SCIP_PARAMETERWRONGVAL )
               {
                  SCIPinfoMessage(scip, NULL, "could not set parameter to random value within range - trying again\n");
               }
            }

            if( retcode != SCIP_OKAY )
            {
               SCIP_CALL( retcode );
            }
            else
            {
               assert(SCIPparamGetReal(param) != oldvalue);

               SCIPdebugMessage("changing real parameter <%s> from %g to %g\n", SCIPparamGetName(param), oldvalue, SCIPparamGetReal(param));

               nchanged++;
               i += step;
            }
            break;
         }

      case SCIP_PARAMTYPE_CHAR:
         {
            SCIP_RETCODE retcode;
            char* allowedvalues;
            char oldvalue;

            /* if the value is already changed, continue looking for an unchanged parameter */
            if( SCIPparamGetChar(param) != SCIPparamGetCharDefault(param) )
               continue;

            retcode = SCIP_OKAY;
            allowedvalues = SCIPparamGetCharAllowedValues(param);
            oldvalue = SCIPparamGetChar(param);

            /* try to change parameter until we reach a different value; give up if there is an error (other than
             * SCIP_PARAMETERWRONGVAL) or there is only one feasible value
             */
            while( ((retcode == SCIP_OKAY && SCIPparamGetChar(param) == oldvalue) || retcode == SCIP_PARAMETERWRONGVAL)
               && (allowedvalues == NULL || strlen(allowedvalues) > 1) )
            {
               char newvalue;

               if( allowedvalues == NULL )
               {
                  newvalue = (char)SCIPrandomGetInt(randgen, (int)CHAR_MIN, (int)CHAR_MAX);
               }
               else
               {
                  int length;
                  int pos;

                  length = strlen(allowedvalues);
                  assert(length >= 1);

                  pos = SCIPrandomGetInt(randgen, 0, length - 1);
                  newvalue = allowedvalues[pos];
               }

               retcode = SCIPchgCharParam(scip, param, newvalue);
               if( retcode == SCIP_PARAMETERWRONGVAL )
               {
                  SCIPinfoMessage(scip, NULL, "could not set parameter to random value within range - trying again");
               }
            }

            if( retcode != SCIP_OKAY )
            {
               SCIP_CALL( retcode );
            }
            else if( SCIPparamGetChar(param) != oldvalue )
            {
               SCIPdebugMessage("changing char parameter <%s> from %c to %c\n", SCIPparamGetName(param), oldvalue,
                  SCIPparamGetChar(param));
               nchanged++;
               i += step;
            }
            else
            {
               SCIPdebugMessage("leaving char parameter <%s> at %c\n", SCIPparamGetName(param), SCIPparamGetChar(param));
               nchanged++;
               i += step;
            }
            break;
         }

      case SCIP_PARAMTYPE_STRING:
         /* we currently exclude string parameters; this would need to be done manually */
         SCIPABORT();
         break;

      default:
         SCIPABORT();
         break;
      }
   }

   return SCIP_OKAY;
}

/** creates a settings file with random parameters */
static
SCIP_RETCODE createSettingsFile(
   unsigned int          seed                /**< random seed */
   )
{
   SCIP_RANDNUMGEN* randgen;
   SCIP* scip = NULL;
   char filename[256];
   int status;

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
   SCIPdebug( SCIPprintVersion(scip, NULL) );

   /* create file name */
   status = snprintf(filename, 255, "%s-%u.set", SCIPgetGitHash(), seed);
   if( status < 0 || status >= 256 )
   {
      SCIPerrorMessage("error creating filename: snprintf() returned %d\n", status);
      return SCIP_ERROR;
   }

   /* generate random number generator */
   SCIP_CALL( SCIPcreateRandom(scip, &randgen, seed, TRUE) );

   /* generate random setting and write it to file */
   SCIP_CALL( generateRandomSettings(scip, randgen, MAX_CHANGE, CHANGE_ADVANCED) );
   SCIP_CALL( SCIPwriteParams(scip, filename, TRUE, TRUE) );
   SCIPinfoMessage(scip, NULL, "-> saved random settings to file <%s>\n", filename);

   /* close SCIP */
   SCIPfreeRandom(scip, &randgen);
   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main method  */
int main(
   int                   argc,               /**< number of arguments from the shell */
   char**                argv                /**< array of shell arguments */
   )
{
   SCIP_RETCODE retcode = SCIP_ERROR;

   if( argc >= 2 )
   {
      /* because of the filename we require a non-negative seed */
      if( strtol(argv[1], NULL, 10) < 0 )
      {
         SCIPerrorMessage("error creating filename: seed value is negative\n");
         return -1;
      }

      retcode = createSettingsFile((unsigned int)strtol(argv[1], NULL, 10));
   }
   else
      retcode = createSettingsFile((unsigned int)time(NULL));

   if( retcode != SCIP_OKAY )
   {
      printf("usage: %s [<seed>]\n\n", argv[0]);
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}
