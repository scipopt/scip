// Copyright (C) GAMS Development and others 2009-2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

#ifndef GAMSOPTIONS_HPP_
#define GAMSOPTIONS_HPP_

#include <cstdlib>

struct gevRec;
struct optRec;

/** Wrapper class for GAMS option file handler.
 * Provides methods to simplify reading and setting of options.
 * Reads the file "<systemdir>/opt<solvername>.def" to retrieve list of supported options.
 */
class GamsOptions
{
private:
   struct gevRec*        gev;                /**< GAMS environment */
   struct optRec*        optionshandle;      /**< GAMS options handler object */

   bool                  opt_is_own;         /**< have we created the optionshandle? */

   /** initializes options object (create options handler, read options definitions file) */
   bool initOpt(
      const char*        solvername          /**< name of solver which options definitions file should be read */
   );

   /** print message queue from options handler */
   void printOptMessages();

public:
   GamsOptions(
      struct gevRec*     gev_ = NULL,        /**< GAMS environment to get access to the system directory name and other stuff, can be NULL and set later by setGEV */
      struct optRec*     opt_ = NULL         /**< GAMS options handler, if NULL and setOpt is not called, an own one might be created */
   )
   : gev(gev_),
     optionshandle(NULL),
     opt_is_own(true)
   {
      if( opt_ != NULL )
         setOpt(opt_);
   }

   ~GamsOptions();

   /** set GAMS environment */
   void setGEV(
      struct gevRec*     gev_                /**< GAMS environment */
   )
   {
      gev = gev_;
   }

   /** set GAMS options handler */
   void setOpt(
      struct optRec*     opt_                /**< GAMS options handler */
   );

   /** read an options file */
   bool readOptionsFile(
      const char*        solvername,         /**< name of solver */
      const char*        optfilename         /**< name of option file, or NULL to read nothing */
   );

   /** checks whether an option exists */
   bool isKnown(
      const char*        optname             /**< name of option to check */
   );

   /** checks whether user specified an option */
   bool isDefined(
      const char*        optname             /**< name of option to check */
   );

   /** get value of a boolean option */
   inline bool getBool(
      const char*        optname             /**< name of option to check */
   )
   {
      return getInteger(optname) != 0;
   }

   /** get value of an integer option */
   int getInteger(
      const char*        optname             /**< name of option to check */
   );

   /** gets value of a real (double) option */
   double getDouble(
      const char*        optname             /**< name of option to check */
   );

   /** gets value of a string option */
   char* getString(
      const char*        optname,            /**< name of option to check */
      char*              buffer              /**< buffer where value can be stored, should have size at least GMS_SSSIZE (256) */
   );

   /** sets value of a boolean option */
   inline void setBool(
      const char*        optname,            /**< name of option to set */
      bool               bval                /**< new value of option */
   )
   {
      setInteger(optname, bval != 0 ? 1 : 0);
   }

   /** sets the value of an integer option */
   void setInteger(
      const char*        optname,            /**< name of option to set */
      int                ival                /**< new value of option */
   );

   /** sets the value of a double option */
   void setDouble(
      const char*        optname,            /**< name of option to set */
      double             dval                /**< new value of option */
   );

   /** sets the value of a string option */
   void setString(
      const char*        optname,            /**< name of option to set */
      const char*        sval                /**< new value of option */
   );
};

#endif /*GAMSOPTIONS_HPP_*/
