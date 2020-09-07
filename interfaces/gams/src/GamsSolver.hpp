// Copyright (C) GAMS Development and others 2009-2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

#ifndef GAMSSOLVER_HPP_
#define GAMSSOLVER_HPP_

typedef struct gmoRec* gmoHandle_t;
typedef struct gevRec* gevHandle_t;
typedef struct optRec* optHandle_t;
typedef struct palRec* palHandle_t;

/** abstract interface to a solver that takes Gams Modeling Objects (GMO) as input */
class GamsSolver
{
protected:

   void initLicensing(
      struct gmoRec*     gmo,                /**< GAMS modeling object */
      struct palRec*     pal                 /**< GAMS audit and license object */
      );

   /** initializes GAMS licensed HSL routines if a commercial Ipopt license is available
    *
    * Returns true if a commercial Ipopt license is available, and false otherwise.
    */
#ifdef GAMS_BUILD
   bool HSLInit(
      struct gmoRec*     gmo,                /**< GAMS modeling object */
      struct palRec*     pal                 /**< GAMS audit and license object */
   );
#else
   bool HSLInit(
      struct gmoRec*     gmo,                /**< GAMS modeling object */
      struct palRec*     pal                 /**< GAMS audit and license object */
   ) const
   { /*lint --e{715}*/
      return false;
   }
#endif

public:
   virtual ~GamsSolver() {};

   /** sets number of OpenMP/GotoBlas threads */
   static void setNumThreads(
      struct gevRec*      gev,                /**< GAMS environment */
      int                 nthreads            /**< number of threads for OpenMP/GotoBlas */
   );

   /** ensure that GMO library has been loaded
    * @return nonzero on failure, 0 on success
    */
   static int getGmoReady();

   /** ensure that GEV library has been loaded
    * @return nonzero on failure, 0 on success
    */
   static int getGevReady();

#if 0
   /** calls GAMS license check, if build by GAMS
    * @return True if license check was skipped or successful or model fits into demo size restrictions
    */
   bool checkGamsLicense(
      struct gmoRec*     gmo,                /**< GAMS modeling object */
      struct palRec*     pal                 /**< GAMS audit and license object */
   );
#endif

   /** registers a GAMS/CPLEX license, if build by GAMS
    * @return True if license was registered or no CPLEX available
    */
   bool checkCplexLicense(
      struct gmoRec*     gmo,                /**< GAMS modeling object */
      struct palRec*     pal                 /**< GAMS audit and license object */
   );

   /** checks for GAMS/Ipopt commercial license
    * @return True if Ipopt commercial license was found, false otherwise (even for demo models).
    */
   bool checkIpoptLicense(
      struct gmoRec*     gmo,                /**< GAMS modeling object */
      struct palRec*     pal                 /**< GAMS audit and license object */
   );

   /** checks for GAMS/Scip commercial license
    * @return True if Scip commercial license was found, false otherwise (even for demo models).
    */
   bool checkScipLicense(
      struct gmoRec*     gmo,                /**< GAMS modeling object */
      struct palRec*     pal                 /**< GAMS audit and license object */
   );

   /** initialization of solver interface and solver
    * The method should do initializes necessary for solving a given model.
    * In case of an error, the method should return with a nonzero return code.
    */
   virtual int readyAPI(
      struct gmoRec*     gmo,                /**< GAMS modeling object */
      struct optRec*     opt                 /**< GAMS options object */
   ) = 0;

   /** notifies solver that the GMO object has been modified and changes should be passed forward to the solver */
   virtual int modifyProblem()
   {
      return -1;
   }

   /** initiate solving and storing solution information in GMO object */
   virtual int callSolver() = 0;
};

#endif /*GAMSSOLVER_HPP_*/
