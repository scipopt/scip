// Copyright (C) GAMS Development and others 2009-2011
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

#ifndef GAMSSCIP_HPP_
#define GAMSSCIP_HPP_

#include "GamsSolver.hpp"

#include <cstdlib>

#include "scip/type_retcode.h"

typedef struct Scip             SCIP;
typedef struct SCIP_Messagehdlr SCIP_MESSAGEHDLR;

/** GAMS interface to SCIP */
class GamsScip : public GamsSolver
{
   friend void printSCIPOptions();
private:
   struct gmoRec*        gmo;                /**< GAMS modeling object */
   struct gevRec*        gev;                /**< GAMS environment */
   struct palRec*        pal;                /**< GAMS audit and license object */

   SCIP*                 scip;               /**< SCIP structure */
   bool                  ipoptlicensed;      /**< whether a commercial Ipopt license is available */

   SCIP_RETCODE setupSCIP();
   SCIP_RETCODE freeSCIP();

public:
   GamsScip()
   : gmo(NULL),
     gev(NULL),
     pal(NULL),
     scip(NULL),
     ipoptlicensed(false)
   { }

   ~GamsScip();

   int readyAPI(
      struct gmoRec*     gmo_,               /**< GAMS modeling object */
      struct optRec*     opt_                /**< GAMS options object */
   );

   int callSolver();

   /** indicates whether the solver interface and solver supports the modifyProblem call */
   static int haveModifyProblem()
   {
      return -1;
   }
};

#endif /*GAMSSCIP_HPP_*/
