// Copyright (C) GAMS Development and others 2009-2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

#include "GamsSolver.hpp"
#ifdef HAVE_CONFIG_H
#include "GAMSlinksConfig.h"
#endif

#include <cstdio>
#include <cassert>

#include "gmomcc.h"
#include "gevmcc.h"

#ifdef GAMS_BUILD
#include "palmcc.h"
#endif

#include "GamsCompatibility.h"

#ifdef HAVE_GOTO_SETNUMTHREADS
extern "C" void goto_set_num_threads(int);
#endif

#ifdef HAVE_MKL_SETNUMTHREADS
extern "C" void MKL_Set_Num_Threads(int);
#endif

void GamsSolver::setNumThreadsLinearAlgebra(
   struct gevRec*      gev,                /**< GAMS environment */
   int                 nthreads            /**< number of threads for BLAS routines */
)
{
#ifdef HAVE_GOTO_SETNUMTHREADS
   if( gev != NULL && nthreads > 1 )
   {
      char msg[100];
      sprintf(msg, "Limit number of threads in GotoBLAS to %d.\n", nthreads);
      gevLogPChar(gev, msg);
   }
   goto_set_num_threads(nthreads);
#endif

#ifdef HAVE_MKL_SETNUMTHREADS
   if( gev != NULL && nthreads > 1 )
   {
      char msg[100];
      sprintf(msg, "Limit number of threads in MKL BLAS to %d.\n", nthreads);
      gevLogPChar(gev, msg);
   }
   MKL_Set_Num_Threads(nthreads);
#endif
}

int GamsSolver::getGmoReady()
{
   if( !gmoLibraryLoaded() )
   {
      char msg[256];

      if( !gmoGetReady(msg, sizeof(msg)) )
      {
         fprintf(stderr, "Error loading GMO library: %s\n",msg);
         return 1;
      }
   }

   return 0;
}

int GamsSolver::getGevReady()
{
   if( !gevLibraryLoaded() )
   {
      char msg[256];

      if( !gevGetReady(msg, sizeof(msg)) )
      {
         fprintf(stderr, "Error loading GEV library: %s\n",msg);
         return 1;
      }
   }

   return 0;
}

void GamsSolver::initLicensing(
   struct gmoRec*     gmo,                /**< GAMS modeling object */
   struct palRec*     pal                 /**< GAMS audit and license object */
   )
{
#ifdef GAMS_BUILD
   char ll[80];
   gevRec* gev;

   assert(pal != NULL);

   gev = (gevRec*)gmoEnvironment(gmo);

   palLicenseRegisterGAMS(pal, 1, gevGetStrOpt(gev, "License1", ll));
   palLicenseRegisterGAMS(pal, 2, gevGetStrOpt(gev, "License2", ll));
   palLicenseRegisterGAMS(pal, 3, gevGetStrOpt(gev, "License3", ll));
   palLicenseRegisterGAMS(pal, 4, gevGetStrOpt(gev, "License4", ll));
   palLicenseRegisterGAMS(pal, 5, gevGetStrOpt(gev, "License5", ll));
   palLicenseRegisterGAMSDone(pal);

   palLicenseCheck(pal,gmoM(gmo),gmoN(gmo),gmoNZ(gmo),gmoNLNZ(gmo),gmoNDisc(gmo));
#endif
}

#if 0
bool GamsSolver::checkGamsLicense(
   struct gmoRec*     gmo,                /**< GAMS modeling object */
   struct palRec*     pal                 /**< GAMS audit and license object */
)
{
#ifdef GAMS_BUILD
   gevRec* gev = (gevRec*)gmoEnvironment(gmo);

   if( palLicenseCheck(pal,gmoM(gmo),gmoN(gmo),gmoNZ(gmo),gmoNLNZ(gmo),gmoNDisc(gmo)) )
   {
      char msg[256];
      gevLogStat(gev, "The license check failed:");
      while( palLicenseGetMessage(pal, msg, sizeof(msg)) )
         gevLogStat(gev, msg);
      return false;
   }
#endif
   return true;
}
#endif

bool GamsSolver::checkCplexLicense(
   struct gmoRec*     gmo,                /**< GAMS modeling object */
   struct palRec*     pal                 /**< GAMS audit and license object */
)
{
   assert(gmo != NULL);
#ifdef GAMS_BUILD
   gevRec* gev = (gevRec*)gmoEnvironment(gmo);

   if( !palLicenseIsDemoCheckout(pal) )
   {
      if( palLicenseCheckSubSys(pal, const_cast<char*>("OCCPCL")) )
      {
         gevLogStat(gev,"***");
         gevLogStat(gev,"*** LICENSE ERROR:");
         gevLogStat(gev,"*** See http://www.gams.com/osicplex/ for OsiCplex licensing information.");
         gevLogStat(gev,"***");
         return false;
      }
   }
#endif
   return true;
}

bool GamsSolver::checkIpoptLicense(
   struct gmoRec*     gmo,                /**< GAMS modeling object */
   struct palRec*     pal                /**< GAMS audit and license object */
)
{
#ifdef GAMS_BUILD
   assert(gmo != NULL);
   if( !palLicenseIsDemoCheckout(pal) )
   {
      if( palLicenseCheckSubSys(pal, const_cast<char*>("IP")) )
         return false;
      else
         return true;
   }
#endif

   return true;
}

bool GamsSolver::checkScipLicense(
   struct gmoRec*     gmo,                /**< GAMS modeling object */
   struct palRec*     pal                /**< GAMS audit and license object */
)
{
#ifdef GAMS_BUILD
   assert(gmo != NULL);
   if( !palLicenseIsDemoCheckout(pal) && !palLicenseIsAcademic(pal) )
   {
      if( palLicenseCheckSubSys(pal, const_cast<char*>("SC")) )
         return false;
      else
         return true;
   }
#endif

   return true;
}
