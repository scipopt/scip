// Copyright (C) GAMS Development and others 2009-2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

#include "GamsOptions.hpp"
#include "GAMSlinksConfig.h"

#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>

#include "gclgms.h"
#include "gevmcc.h"
#include "optcc.h"

GamsOptions::~GamsOptions()
{
   if( opt_is_own )
   {
      if( optionshandle != NULL )
         optFree(&optionshandle);

      /* try to unload library, since it seems to have an instance counter, nothing should happen if it is still in use somewhere else */
      optLibraryUnload();
   }
}

bool GamsOptions::initOpt(
   const char*        solvername          /**< name of solver which options definitions file should be read */
)
{
   assert(gev != NULL);

   char buffer[512];

   /* get the Option File Handling set up */
   if( !optCreate(&optionshandle, buffer, 512) )
   {
      gevLogStatPChar(gev, "\n*** Could not create optionfile handle: ");
      gevLogStat(gev, buffer);
      return false;
   }
   opt_is_own = true;

   assert(sizeof(buffer) >= GMS_SSSIZE);
   gevGetStrOpt(gev, gevNameSysDir, buffer);
   buffer[sizeof(buffer)-1] = '\0';
   size_t len = strlen(buffer);
   if( snprintf(buffer+len, 512-len, "opt%s.def", solvername) >= 512 )
   {
      gevLogStat(gev, "\n*** Path to GAMS system directory too long.");
      return false;
   }

   if( optReadDefinition(optionshandle, buffer) )
   {
      printOptMessages();
      return false;
   }
   optEOLOnlySet(optionshandle, 1);

   return true;
}

void GamsOptions::printOptMessages()
{
   assert(optionshandle != NULL);
   assert(gev != NULL);

   if( optMessageCount(optionshandle) == 0 )
      return;

   char buffer[GMS_SSSIZE];
   int itype;
   for( int i = 1; i <= optMessageCount(optionshandle); ++i )
   {
      optGetMessage(optionshandle, i, buffer, &itype);
      if( itype <= optMsgFileLeave || itype == optMsgUserError )
         gevLogStat(gev, buffer);
   }
   optClearMessages(optionshandle);
}

void GamsOptions::setOpt(
   struct optRec*     opt_                /**< GAMS options handler */
)
{
   assert(optionshandle == NULL);

   char buffer[512];
   if( !optGetReady(buffer, 512) )
   {
      gevLogStatPChar(gev, "\n*** Could not load optionfile library: ");
      gevLogStat(gev, buffer);
      return;
   }

   opt_is_own = false;
   optionshandle = opt_;
}

bool GamsOptions::readOptionsFile(
   const char*        solvername,         /**< name of solver */
   const char*        optfilename         /**< name of option file, or NULL to read nothing */
)
{
   if( optionshandle == NULL && !initOpt(solvername) )
      return false;

   if( optfilename == NULL )
      return true;

   /* read option file */
   optEchoSet(optionshandle, 1);
   optReadParameterFile(optionshandle, optfilename);
   if( optMessageCount(optionshandle) )
   {
      printOptMessages();
      optEchoSet(optionshandle, 0);
      return false;
   }
   else
   {
      optEchoSet(optionshandle, 0);
      return true;
   }
}

bool GamsOptions::isKnown(
   const char*        optname             /**< name of option to check */
)
{
   assert(optionshandle != NULL);

   int i;
   int refNum;

   return optFindStr(optionshandle, optname, &i, &refNum) > 0;
}

bool GamsOptions::isDefined(
   const char*        optname             /**< name of option to check */
)
{
   assert(optionshandle != NULL);

   int i;
   int refNum;

   if( optFindStr(optionshandle, optname, &i, &refNum) ) {
      int isdefined;
      int isdefinedrecent;
      int dataType;
      int optType;
      int subType;

      optGetInfoNr(optionshandle, i, &isdefined, &isdefinedrecent, &refNum, &dataType, &optType, &subType);

      return isdefined;
   }

   gevLogStatPChar(gev, "*** Internal Error. Unknown option ");
   gevLogStat(gev, optname);

   return false;
}

int GamsOptions::getInteger(
   const char*        optname             /**< name of option to check */
)
{
#ifndef NDEBUG
   assert(optionshandle != NULL);

   int i, refNum, isdefined, isDefinedRecent, dataType, optType, subType;

   assert(optFindStr(optionshandle, optname, &i, &refNum) > 0);
   optGetInfoNr(optionshandle, i, &isdefined, &isDefinedRecent, &refNum, &dataType, &optType, &subType);
   assert(dataType == optDataInteger);
#endif

   return optGetIntStr(optionshandle, optname);
}

double GamsOptions::getDouble(
   const char*        optname             /**< name of option to check */
)
{
#ifndef NDEBUG
   assert(optionshandle != NULL);

   int i, refNum, isdefined, isDefinedRecent, dataType, optType, subType;

   assert(optFindStr(optionshandle, optname, &i, &refNum) > 0);
   optGetInfoNr(optionshandle, i, &isdefined, &isDefinedRecent, &refNum, &dataType, &optType, &subType);
   assert(dataType == optDataDouble);
#endif

   return optGetDblStr(optionshandle, optname);
}

char* GamsOptions::getString(
   const char*        optname,            /**< name of option to check */
   char*              buffer              /**< buffer where value can be stored, should have size at least GMS_SSSIZE (256) */
)
{
#ifndef NDEBUG
   assert(optionshandle != NULL);

   int i, refNum, isdefined, isDefinedRecent, dataType, optType, subType;

   assert(optFindStr(optionshandle, optname, &i, &refNum) > 0);
   optGetInfoNr(optionshandle, i, &isdefined, &isDefinedRecent, &refNum, &dataType, &optType, &subType);
   assert(dataType == optDataString);
#endif

   return optGetStrStr(optionshandle, optname, buffer);
}

void GamsOptions::setInteger(
   const char*        optname,            /**< name of option to set */
   int                ival                /**< new value of option */
)
{
#ifndef NDEBUG
   assert(optionshandle != NULL);

   int i, refNum, isdefined, isDefinedRecent, dataType, optType, subType;

   assert(optFindStr(optionshandle, optname, &i, &refNum) > 0);
   optGetInfoNr(optionshandle, i, &isdefined, &isDefinedRecent, &refNum, &dataType, &optType, &subType);
   assert(dataType == optDataInteger);
#endif

   optSetIntStr(optionshandle, optname, ival);
   printOptMessages();
}

void GamsOptions::setDouble(
   const char*        optname,            /**< name of option to set */
   double             dval                /**< new value of option */
)
{
#ifndef NDEBUG
   assert(optionshandle != NULL);

   int i, refNum, isdefined, isDefinedRecent, dataType, optType, subType;

   assert(optFindStr(optionshandle, optname, &i, &refNum) > 0);
   optGetInfoNr(optionshandle, i, &isdefined, &isDefinedRecent, &refNum, &dataType, &optType, &subType);
   assert(dataType == optDataDouble);
#endif

   optSetDblStr(optionshandle, optname, dval);
   printOptMessages();
}

void GamsOptions::setString(
   const char*        optname,            /**< name of option to set */
   const char*        sval                /**< new value of option */
)
{
#ifndef NDEBUG
   assert(optionshandle != NULL);

   int i, refNum, isdefined, isDefinedRecent, dataType, optType, subType;

   assert(optFindStr(optionshandle, optname, &i, &refNum) > 0);
   optGetInfoNr(optionshandle, i, &isdefined, &isDefinedRecent, &refNum, &dataType, &optType, &subType);
   assert(dataType == optDataString);
#endif

   optSetStrStr(optionshandle, optname, sval);
   printOptMessages();
}

