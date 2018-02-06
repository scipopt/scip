// Copyright (C) GAMS Development and others 2009-2012
// All Rights Reserved.
// This code is published under the Eclipse Public License.
//
// Author: Stefan Vigerske

#include "GamsScip.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>

// GAMS
#include "gmomcc.h"
#include "gevmcc.h"
#ifdef GAMS_BUILD
#include "gevlice.h" // for xpress license setup
#include "palmcc.h"
#endif

#include "GAMSlinksConfig.h"
#include "GamsCompatibility.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/githash.c"
#include "nlpi/nlpi_ipopt.h"
#include "reader_gmo.h"
#include "event_solvetrace.h"
/* #include "prop_defaultbounds.h" */

/* #include "lpiswitch.h" */

static
SCIP_DECL_ERRORPRINTING(printErrorGev)
{
   assert(data != NULL);
   gevLogStatPChar((gevHandle_t)data, msg);
}

static
SCIP_DECL_MESSAGEWARNING(GamsScipPrintLogStat)
{
   assert(SCIPmessagehdlrGetData(messagehdlr) != NULL);
   assert(file != NULL);

   if( file != stderr )
      fputs(msg, file);
   else
   {
      void* gev = (void*)SCIPmessagehdlrGetData(messagehdlr);
      gevLogStatPChar((gevHandle_t)gev, msg);
   }
}

static
SCIP_DECL_MESSAGEINFO(GamsScipPrintLog)
{
   assert(SCIPmessagehdlrGetData(messagehdlr) != NULL);
   assert(file != NULL);

   if( file != stdout )
      fputs(msg, file);
   else
   {
      void* gev = (void*)SCIPmessagehdlrGetData(messagehdlr);
      gevLogPChar((gevHandle_t)gev, msg);
   }
}

/*
static
SCIP_DECL_PARAMCHGD(GamsScipParamChgdLpSolver)
{
   SCIP_LPISW_LPSOLVER lpsolver;

   assert(scip != NULL);

   if( SCIPgetStage(scip) >= SCIP_STAGE_INITPRESOLVE && SCIPgetStage(scip) <= SCIP_STAGE_EXITSOLVE )
   {
      SCIPerrorMessage("Cannot change SCIP LP Solver during presolve or solve.\n");
      return SCIP_INVALIDCALL;
   }

   for( lpsolver = SCIP_LPISW_LPSOLVER(0); lpsolver < SCIP_LPISW_NSOLVERS; lpsolver = SCIP_LPISW_LPSOLVER((int)lpsolver + 1) )
      if( strcmp(SCIPparamGetString(param), SCIP_LPISW_SOLVERNAMES[lpsolver]) == 0 )
         break;
   if( lpsolver == SCIP_LPISW_NSOLVERS )
   {
      SCIPerrorMessage("Value '%s' for parameter lp/solver not understood.\n", SCIPparamGetString(param));
      return SCIP_PARAMETERWRONGVAL;
   }

   if( SCIPlpiSwitchGetCurrentSolver() != lpsolver )
   {
      SCIP_CALL( SCIPlpiSwitchSetSolver(lpsolver) );

      SCIPinfoMessage(scip, NULL, "Changed LP solver to\n  %-20s %s\n", SCIPlpiGetSolverName(), SCIPlpiGetSolverDesc());
   }

   return SCIP_OKAY;
}
*/

GamsScip::~GamsScip()
{
   SCIP_CALL_ABORT( freeSCIP() );

#ifdef GAMS_BUILD
   if( pal != NULL )
      palFree(&pal);

   if( calledxprslicense )
      gevxpressliceFreeTS();
#endif

}

int GamsScip::readyAPI(
   struct gmoRec*     gmo_,               /**< GAMS modeling object */
   struct optRec*     opt_                /**< GAMS options object */
)
{
   char buffer[512];

   gmo = gmo_;
   assert(gmo != NULL);

   if( getGmoReady() || getGevReady() )
      return 1;

   gev = (gevRec*)gmoEnvironment(gmo);
   assert(gev != NULL);

   ipoptlicensed = false;
#ifdef GAMS_BUILD
   if( pal == NULL && !palCreate(&pal, buffer, sizeof(buffer)) )
      return 1;

#define PALPTR pal
#include "coinlibdCL5svn.h" 
   palGetAuditLine(pal, buffer);
   gevLogStat(gev, "");
   gevLogStat(gev, buffer);
   gevStatAudit(gev, buffer);

   initLicensing(gmo, pal);

#ifdef COIN_HAS_OSIXPR
   /* Xpress license setup - don't say anything if failing, since Xpress is not used by default */
   if( !calledxprslicense )
   {
      XPlicenseInit_t initType;
      int initRC;
      gevxpressliceInitTS(gev, pal, gmoM(gmo), gmoN(gmo), gmoNZ(gmo), gmoNLNZ(gmo), gmoNDisc(gmo), 0, &initType, &initRC, buffer, sizeof(buffer));
      calledxprslicense = true;
   }
#endif
#endif

   // check for academic license, or if we run in demo mode
   if( !checkScipLicense(gmo, pal) )
   {
      gevLogStat(gev, "*** No SCIP license available.");
      gevLogStat(gev, "*** Please contact sales@gams.com to arrange for a license.");
      gmoSolveStatSet(gmo, gmoSolveStat_License);
      gmoModelStatSet(gmo, gmoModelStat_LicenseError);
      return 1;
   }

   ipoptlicensed = HSLInit(gmo, pal);

   // print version info and copyright
   if( SCIPsubversion() > 0 )
      sprintf(buffer, "SCIP version %d.%d.%d.%d (" SCIP_GITHASH ")\n", SCIPmajorVersion(), SCIPminorVersion(), SCIPtechVersion(), SCIPsubversion());
   else
      sprintf(buffer, "SCIP version %d.%d.%d (" SCIP_GITHASH ")\n", SCIPmajorVersion(), SCIPminorVersion(), SCIPtechVersion());
   gevLogStatPChar(gev, buffer);
   gevLogStatPChar(gev, SCIP_COPYRIGHT"\n\n");

   // install or update error printing callback in SCIP to use current gev
   SCIPmessageSetErrorPrinting(printErrorGev, (void*)gev);

   // setup (or reset) SCIP instance
   SCIP_RETCODE scipret;
   scipret = setupSCIP();

   if( scipret != SCIP_OKAY )
   {
      snprintf(buffer, sizeof(buffer), "Error %d in call of SCIP function\n", scipret);
      gevLogStatPChar(gev, buffer);
      gmoSolveStatSet(gmo, gmoSolveStat_SystemErr);
      gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
      return 1;
   }

   assert(scip != NULL);

   // print info on used external codes
   SCIPprintExternalCodes(scip, NULL);

   return 0;
}

int GamsScip::callSolver()
{
   assert(gmo  != NULL);
   assert(gev  != NULL);
   assert(scip != NULL);

   /* set interface type so we see =B= and =X= equations */
   gmoInterfaceSet(gmo, gmoIFace_Raw);

   if( gmoGetEquTypeCnt(gmo, gmoequ_C) || gmoGetEquTypeCnt(gmo, gmoequ_B) || gmoGetEquTypeCnt(gmo, gmoequ_X) )
   {
      gevLogStat(gev, "ERROR: Conic and logic constraints and external functions not supported by SCIP interface.\n");
      gmoSolveStatSet(gmo, gmoSolveStat_Capability);
      gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
      return 1;
   }

   // set number of threads for linear algebra routines used in Ipopt
   setNumThreads(gev, gevThreads(gev));

   // update error printing callback in SCIP to use current gev
   SCIPmessageSetErrorPrinting(printErrorGev, (void*)gev);

   SCIP_RETCODE scipret;

   // let GMO reader setup SCIP parameters and read options file
   // do this here already so we know how to assemble dialog
   scipret = SCIPreadParamsReaderGmo(scip);
   if( scipret != SCIP_OKAY )
   {
      char buffer[256];
      sprintf(buffer, "Error %d in call of SCIP function\n", scipret);
      gevLogStatPChar(gev, buffer);
      gmoSolveStatSet(gmo, gmoSolveStat_SystemErr);
      gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
      return 1;
   }
   SCIPinfoMessage(scip, NULL, "non-default parameter settings:\n");
   SCIPwriteParams(scip, NULL, FALSE, TRUE);

   char* interactive = NULL;
   SCIP_CALL_ABORT( SCIPgetStringParam(scip, "gams/interactive", &interactive) );
   assert(interactive != NULL);
#ifdef GAMS_BUILD
   if( interactive[0] != '\0' && !palLicenseIsAcademic(pal) && palLicenseCheckSubSys(pal, const_cast<char*>("SC")) )
   {
      gevLogStat(gev, "SCIP interactive shell not available in demo mode.\n");
      interactive[0] = '\0';
   }
#endif

   SCIP_Bool printstat;
   SCIP_CALL_ABORT( SCIPgetBoolParam(scip, "display/statistics", &printstat) );

   char* attrfile = NULL;
#if 0
   SCIP_CALL( SCIPgetStringParam(scip, "constraints/attrfile", &attrfile) );
#endif

   // setup commands to be executed by SCIP
   SCIP_CALL_ABORT( SCIPaddDialogInputLine(scip, "readgams") );              // setup model

   if( attrfile != NULL && *attrfile != '\0' )
   {
      char buffer[SCIP_MAXSTRLEN + 10];
      size_t len;

      len = strlen(attrfile);
      if( len >= 3 && strcmp(&attrfile[len-3], ".ca") == 0 )
         (void) SCIPsnprintf(buffer, sizeof(buffer), "read %g", attrfile);
      else
         (void) SCIPsnprintf(buffer, sizeof(buffer), "read %g ca", attrfile);
      SCIP_CALL_ABORT( SCIPaddDialogInputLine(scip, buffer) );               // process constraints attribute file
   }

   if( interactive[0] == '\0' )
   {
      SCIP_CALL_ABORT( SCIPaddDialogInputLine(scip, "optimize") );           // solve model

      if( printstat )
      {
         SCIP_CALL_ABORT( SCIPaddDialogInputLine(scip, "disp statistics") ); // display solution statistics
      }
      SCIP_CALL_ABORT( SCIPaddDialogInputLine(scip, "write gamssol") );      // pass solution to GMO

      SCIP_CALL_ABORT( SCIPaddDialogInputLine(scip, "quit") );               // quit shell
   }
   else
   {
      // pass user commands to shell
      SCIP_CALL_ABORT( SCIPaddDialogInputLine(scip, interactive) );
   }


   // run SCIP
   scipret = SCIPstartInteraction(scip);

   // evaluate SCIP return code
   switch( scipret )
   {
      case SCIP_OKAY:
         break;

      case SCIP_READERROR:
         /* if it's readerror, then we guess that it comes from encountering an unsupported gams instruction in the gmo readers makeExprtree method
          * we still return with zero then
          */
         gmoModelStatSet(gmo, gmoModelStat_NoSolutionReturned);
         gmoSolveStatSet(gmo, gmoSolveStat_Capability);
         break;

      case SCIP_LPERROR:
      case SCIP_MAXDEPTHLEVEL:
         /* if SCIP failed due to internal error (forced LP solve failed, max depth level reached), also return zero */
         gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
         gmoSolveStatSet(gmo, gmoSolveStat_SolverErr);
         break;

      case SCIP_NOMEMORY:
         /* there is no extra solver status for running out of memory, but memory is a resource, so return this */
         gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
         gmoSolveStatSet(gmo, gmoSolveStat_Resource);
         break;

      default:
      {
         char buffer[256];
         sprintf(buffer, "Error %d in call of SCIP function\n", scipret);
         gevLogStatPChar(gev, buffer);

         gmoModelStatSet(gmo, gmoModelStat_ErrorNoSolution);
         gmoSolveStatSet(gmo, gmoSolveStat_SystemErr);
         return 1;
      }
   }

   return 0;
}

SCIP_RETCODE GamsScip::setupSCIP()
{
/*
#ifdef COIN_HAS_OSICPX
   // change default LP solver to CPLEX, if license available
   if( gmo != NULL && checkCplexLicense(gmo, pal) )
   {
      SCIP_CALL( SCIPlpiSwitchSetSolver(SCIP_LPISW_CPLEX) );
   }
#endif
*/

   if( scip == NULL )
   {
      // if called first time, create a new SCIP instance and include all plugins that we need and setup interface parameters
      SCIP_MESSAGEHDLR* messagehdlr;
      SCIP_NLPI* nlpiipopt;

      SCIP_CALL( SCIPcreate(&scip) );

      // create and install our message handler
      SCIP_CALL( SCIPmessagehdlrCreate(&messagehdlr, FALSE, NULL, FALSE,
         GamsScipPrintLogStat, GamsScipPrintLog, GamsScipPrintLog, NULL,
         (SCIP_MESSAGEHDLRDATA*)gev) );
      SCIP_CALL( SCIPsetMessagehdlr(scip, messagehdlr) );
      SCIP_CALL( SCIPmessagehdlrRelease(&messagehdlr) );

      SCIP_CALL( SCIPincludeDefaultPlugins(scip) );
      SCIP_CALL( SCIPincludeReaderGmo(scip) );
      SCIP_CALL( SCIPincludeEventHdlrSolveTrace(scip, gmo) );
      /* SCIP_CALL( SCIPincludePropDefaultBounds(scip) ); */

      if( ipoptlicensed )
      {
         nlpiipopt = SCIPfindNlpi(scip, "ipopt");
         if( nlpiipopt != NULL )
         {
            SCIPsetModifiedDefaultSettingsIpopt(nlpiipopt, "linear_solver ma27\nlinear_system_scaling mc19\n");
            SCIP_CALL( SCIPincludeExternalCodeInformation(scip, "HSL MA27 and MC19", "Harwell Subroutine Libraries (www.hsl.rl.ac.uk) from commercially supported Ipopt") );
         }
      }
/*
      else
      {
         nlpiipopt = SCIPfindNlpi(scip, "ipopt");
         if( nlpiipopt != NULL )
         {
            SCIPsetModifiedDefaultSettingsIpopt(nlpiipopt, "linear_solver mumps\n");
         }
      }
*/

      /* SCIP_CALL( SCIPaddBoolParam(scip, "gams/solvefinal",
       * "whether the problem should be solved with fixed discrete variables to get dual values",
       * NULL, FALSE, TRUE,  NULL, NULL) );
       */
      SCIP_CALL( SCIPaddBoolParam(scip, "display/statistics",
         "whether to print statistics on a solve",
         NULL, FALSE, FALSE, NULL, NULL) );
      SCIP_CALL( SCIPaddStringParam(scip, "gams/interactive",
         "command to be issued to the SCIP shell instead of issuing a solve command",
         NULL, FALSE, "", NULL, NULL) );
#if 0
      SCIP_CALL( SCIPaddStringParam(scip, "constraints/attrfile",
         "name of file that specifies constraint attributes",
         NULL, FALSE, "", NULL, NULL) );
#endif
/*
      SCIP_CALL( SCIPaddStringParam(scip, "lp/solver",
         "LP solver to use (clp, cplex, mosek, soplex, gurobi, xpress)",
         NULL, FALSE, SCIP_LPISW_SOLVERNAMES[SCIPlpiSwitchGetCurrentSolver()], GamsScipParamChgdLpSolver, NULL) );
*/
   }
   else
   {
      // if called before, only clear up problem and reset parameters
      SCIP_CALL( SCIPfreeProb(scip) );
      SCIP_CALL( SCIPresetParams(scip) );
   }

   /** pass current GMO into GMO reader, so it does not read instance from file */
   SCIPsetGMOReaderGmo(scip, gmo);

   return SCIP_OKAY;
}

static
void interruptDuringFree(void)
{
   printf("SCIP is freeing its memory. Please be patient...\n");
}

SCIP_RETCODE GamsScip::freeSCIP()
{
   if( scip != NULL )
   {
      /* catch Ctrl+C while freeing SCIP
       * This can take several seconds after a longer run and the SCIP
       * signal handler is not in place at this point anymore.
       * Without it, the run will be aborted here, thereby not returning
       * results back to GAMS.
       */
      if( gev != NULL )
         gevTerminateSet(gev, NULL, (void*)&interruptDuringFree);

      SCIP_CALL( SCIPfree(&scip) );

      if( gev != NULL )
         gevTerminateUninstall(gev);
   }

   return SCIP_OKAY;
}

#define GAMSSOLVERC_ID         scp
#define GAMSSOLVERC_CLASS      GamsScip
#include "GamsSolverC_tpl.cpp"
